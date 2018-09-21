function beam=fcm_comcohA2Lbeam(comcohfile,radius)
% FCM_COMCOHA2LBEAM  calculates mean functional connectivity of each voxel and
%       compares to homologous contralateral voxel (i.e., creates L-images).
%       A --> All voxels are seed voxels.
%
% beam = fcm_comcohA2Lbeam(comcohfile,�radius�)
%   ( parameters in �� are optional )
%
%  COMCOHFILE  name of file containing connectivity results.
%  RADIUS      Radius of contralateral region in mm to be used for
%              t-test of each voxel. Default is 20 mm.

global nuts

%if (nargin<3 || isempty(remspur)), remspur = false;  end
if nargin<2, radius=20; end

% nuts = load(sessionfile,'coreg','voxels','voxelsize');
load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);
nv=size(nuts.voxels,1);

if nc==size(CC.frq,1)     % Legacy compatibility
    CC.coh=permute(CC.coh,[2 3 1]); 
    [nc,nt,nf]=size(CC.coh);
end    

% Do Z-transform 
switch CC.method
    case {'ccohere' 'nccohere'}
        CC.coh = CC.coh./abs(CC.coh) .* atanh(abs(CC.coh));     % Do Z-transform 
        CC.coh = abs(imag(CC.coh));                             % get abs of imaginary component
    case 'ampcorr'
        CC.coh = abs(atanh(CC.coh));
    case 'pli'
        CC.coh = abs(CC.coh);
end

MRIvoxels = nut_meg2mri(nuts.voxels,nuts.coreg);
T=zeros(nv,nt,nf);
p=ones(nv,nt,nf);
numconn=nv-1;
df=numconn-1;

% Calculate distances
voxma=[abs(MRIvoxels(:,1)) MRIvoxels(:,2:3)];
voxma=repmat(reshape(voxma,[nv 1 3]),[1 nv 1]);
dist=zeros(nv,nv,3);
for k=1:3
    dist(:,:,k)=voxma(:,:,k)-voxma(:,:,k)';
end
clear voxma
dist=sqrt(sum(dist.^2,3));
dist=dist+diag(Inf(nv,1));

% Bring to matrix form for faster performance
CM=NaN(nv,nv,nt,nf);
for k=1:size(CC.comps,1)
    CM(CC.comps(k,1),CC.comps(k,2),:,:)=CC.coh(k,:,:);
    CM(CC.comps(k,2),CC.comps(k,1),:,:)=CC.coh(k,:,:);
end

voxsign=sign(MRIvoxels(:,1));
for vv=1:nv
    e = find( voxsign~=sign(MRIvoxels(vv,1)) );    
    f = find( dist(vv,e) < radius  );
    if ~isempty(f)
        I = shiftdim(CM(vv,:,:,:),1);
        I(vv,:,:)=[];
        C = shiftdim(nanmean(CM(e(f),:,:,:),1),1);
        C(vv,:,:)=[];
        test = I - C;
        T(vv,:,:)=mean(test,1) ./ (std(test,[],1) ./ sqrt(numconn));
        p(vv,:,:)=2*tcdf(-abs(T(vv,:,:)),df);
    end                
end

% Output structure
isspec = ( nt==1 & nf>1 );
if isspec   % if we have a spectrogram, put freq data in time dimension.
    T = permute(T,[1 3 2]);
    p = permute(p,[1 3 2]);
    timepts=mean(CC.frq,2);     
    timewin=CC.frq;    
    nt=nf; nf=1;
    bands=[CC.frq(1,1) CC.frq(end,2)];
else
    timepts=mean(CC.time,2);
    timewin=CC.time;
    bands=CC.frq;
end
if length(timepts)>1
    srate = 1/(timepts(2)-timepts(1));
else
    srate = 1;  % arbitrary
end
beam=struct('s',{{T}},'timewindow',timewin, ...
    'timepts',timepts,'bands',bands,'voxels',nuts.voxels, ...
    'voxelsize',nuts.voxelsize,'srate',srate,'coreg',nuts.coreg);    

beam.sinfo={'T'};
beam.ttest.tail='both';
beam.ttest.T=T;
beam.ttest.p_uncorr = p;
beam.ttest.FDR = 0.01;    

if isspec
    beam.labels.xaxis = 'Frequency (Hz)';
    beam.labels.yaxis = 'T';
else
    beam.labels.xaxis = 'Time (ms)';
    if nf>1
        beam.labels.yaxis = 'Frequency (Hz)';
        beam.labels.colorbar = 'T';
    else
        beam.labels.yaxis = 'T';
    end
end
