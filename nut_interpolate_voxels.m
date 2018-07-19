function beamout=nut_interpolate_voxels(beam)
% NUT_INTERPOLATE_VOXELS  estimates activation values of missing voxels by
% using a linear interpolation of all existing neighbouring voxels.
%
% Usages:
%   nut_interpolate_voxels s_beamfile_spatnorm
%   beam = nut_interpolate_voxels(beam)
%
% beam   structure containing SPATIALLY NORMALIZED (!!!) beamforming data
%
% (c) Adrian G. Guggisberg


if ischar(beam)
    beamfile=beam;
    beam=load(beam);
    beam=nut_beam_legacy_compatibility(beam);
end

% ispop=~iscell(beam.s);

xx = unique(beam.voxels(:,1))';     % Make sure the template grid matches the data voxel grid
yy = unique(beam.voxels(:,2))';
zz = unique(beam.voxels(:,3))';
xx = [ fliplr(xx(1)-beam.voxelsize(1):-beam.voxelsize(1): -90) xx xx(end)+beam.voxelsize(1):beam.voxelsize(1): 90 ];
yy = [ fliplr(yy(1)-beam.voxelsize(2):-beam.voxelsize(2):-130) yy yy(end)+beam.voxelsize(2):beam.voxelsize(2): 90 ];
zz = [ fliplr(zz(1)-beam.voxelsize(3):-beam.voxelsize(3): -80) zz zz(end)+beam.voxelsize(3):beam.voxelsize(3):110 ];
voxels3d = nut_coordgrid(xx,yy,zz);

load MNIvoxels
MNIvoxels=intersect(voxels3d,MNIvoxels,'rows'); 

%MNItfm = [beam.voxelsize(1) 0 0 min(MNIvoxels(:,1)); 0 beam.voxelsize(2) 0 min(MNIvoxels(:,2)); 0 0 beam.voxelsize(3) min(MNIvoxels(:,3)); 0 0 0 1];
%voxelgrid = nut_coordtfm(beam.voxels,inv(MNItfm));

missing=setdiff(MNIvoxels,beam.voxels,'rows');
num = size(missing,1);
nuv = size(beam.voxels,1);
ns  = length(beam.s);
nf  = size(beam.s{1},3);

% if ispop
% end

lim = beam.voxelsize*1.1;
sinterp = repmat({nan(num,length(beam.timepts),size(beam.bands,1))},[1 ns]);
for k=1:num
    vot = abs(beam.voxels-repmat(missing(k,:),[nuv 1]));
    v = find(vot(:,1)<lim(1) & vot(:,2)<lim(2) & vot(:,3)<lim(3));      % find neighbouring voxels
    if ~isempty(v)
        dist = (1./sqrt(sum(vot(v,:).^2,2)))';                           % calculate Euclidian distances to missing voxel
        totdist = sum(dist);                                            % shorter distances weigh more
%         if ispop
%         else
            
        for s=1:ns
            for f=1:nf
                sinterp{s}(k,:,f) = (dist./totdist) * beam.s{s}(v,:,f);     % linear interpolation
            end
        end
%         end
    end
end

v = find(~isnan(sinterp{1}(:,1,1)));
for s=1:ns
    beam.s{s}=cat(1,beam.s{s},sinterp{s}(v,:,:));
end
beam.voxels=cat(1,beam.voxels,missing(v,:));

if nargout<1
    save(beamfile,'-struct','beam')
else
    beamout=beam;
end
