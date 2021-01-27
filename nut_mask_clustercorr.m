function nut_mask_clustercorr(sbeamfile, clusternum, hardthresh, masktype)
% usage:  >nut_mask_clustercorr('s_beamtf12_ttest2.mat', 25, .005, 'GWM')
%
% User restricts the viewable stat result to labeled voxels (AAL/SPM), as 
% inflated to the resolution of the lead field, and defines the cluster
% size (in # of contiguous voxels) used to threshold spatial extent of 
% activity using the uncorrected p-value array and a user-set value (e.g. p<.001).
%
%  masktype options:
%  'GWM' = Grey and White Matter    'LRC' = Left and Right Cerebrum/Cerebellum
%  'LCC' = Left Cerebrum/Cerebellum 'RCC' = Right Cerebrum/Cerebellum
%  'COR' = 4 cortical lobes         'BGS' = Basal Ganglia structures
%  'ACC' = anterior cingulate       'GMO' = Grey Matter Only
%  'Lim' = Limbic/Cingulate         'MNI' = enter limits via dialog/default
%
% Appropriate cluster size depends on masktype & leadfield size.
% Generally, one voxel directly touches 26 other voxels via sides and edges.
% Saves a voi-masked ver4 nutmeg sbeam in p_uncorr and p_clustercorr
% options of nutmeg viewer. Created for group stat sbeams.
%
% dependencies: 
% Nutmeg 4.1 or greater
% nut_clusterstats_hard.m

beam=load(sbeamfile);
if isfield(beam,'s')
    disp('Looks like a Nutmeg 4 structure!')
else
    error('This structure is not from Nutmeg 4.')
end
voxrez=round(beam.voxelsize(1)/3); % accepted distance between closest label coord and beam.voxel coord
%% based on nut_voi_corticalmask_tal 
nm;nut_results_viewer(sbeamfile);
switch masktype
    case 'GWM'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Gray Matter');
        [~, coords_tal_2]=subf_find_talDB_anatindex(sbeamfile, 'White Matter');
        coords_tal = unique([coords_tal_1;coords_tal_2],'rows');
        maskname = 'AALmaskGWMatter';
    case 'GMO'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Gray Matter');
        coords_tal = unique([coords_tal_1],'rows');
        maskname = 'AALmaskGMatterOnly';
    case 'LCC'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Left Cerebrum');
        [~, coords_tal_2]=subf_find_talDB_anatindex(sbeamfile, 'Left Cerebellum');
        coords_tal = unique([coords_tal_1;coords_tal_2],'rows');
        maskname = 'AALmaskLeft';
    case 'RCC'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Right Cerebrum');
        [~, coords_tal_2]=subf_find_talDB_anatindex(sbeamfile, 'Right Cerebellum');
        coords_tal = unique([coords_tal_1;coords_tal_2],'rows');
        maskname = 'AALmaskRight';
    case 'LRC'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Left Cerebrum');
        [~, coords_tal_2]=subf_find_talDB_anatindex(sbeamfile, 'Right Cerebrum');
        [~, coords_tal_3]=subf_find_talDB_anatindex(sbeamfile, 'Left Cerebellum');
        [~, coords_tal_4]=subf_find_talDB_anatindex(sbeamfile, 'Right Cerebellum');
        coords_tal = unique([coords_tal_1;coords_tal_2;coords_tal_3;coords_tal_4],'rows');
        maskname = 'AALmaskCerebrumCBellum';
    case 'COR'  
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Frontal Lobe');
        [~, coords_tal_2]=subf_find_talDB_anatindex(sbeamfile, 'Parietal Lobe');
        [~, coords_tal_3]=subf_find_talDB_anatindex(sbeamfile, 'Temporal Lobe');
        [~, coords_tal_4]=subf_find_talDB_anatindex(sbeamfile, 'Occipital Lobe');
        coords_tal = unique([coords_tal_1;coords_tal_2;coords_tal_3;coords_tal_4],'rows');
        maskname = 'AALmaskCorticalLobes';
    case 'BGS'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Caudate');
        [~, coords_tal_2]=subf_find_talDB_anatindex(sbeamfile, 'Putamen');
        [~, coords_tal_3]=subf_find_talDB_anatindex(sbeamfile, 'Medial Globus Pallidus');
        [~, coords_tal_4]=subf_find_talDB_anatindex(sbeamfile, 'Lateral Globus Pallidus');
        coords_tal = unique([coords_tal_1;coords_tal_2;coords_tal_3;coords_tal_4],'rows');
        maskname = 'AALmaskBasalGanglia';
    case 'ACC'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Anterior Cingulate');
        coords_tal = unique(coords_tal_1,'rows');
        maskname = 'AALmaskACC';
    case 'Lim'
        [~, coords_tal_1]=subf_find_talDB_anatindex(sbeamfile, 'Limbic Lobe');
        coords_tal = unique(coords_tal_1,'rows');
        maskname = 'AALmaskLimbicLabel';
    case 'MNI'
        vox=beam.voxels;
        %enter coordinate limits in x, y, z min/max format
        definput={'-75'; '75'; '-115'; '85'; '-45'; '90'};
        prompt={'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', 'Zmax'};
        dlgtitle = 'Enter numeric VOI boundaries'; dims = [1 35]; opts.Resize='on';
        mnibounds=inputdlg(prompt, dlgtitle, dims, definput, opts);  
        %find coords in beam.voxels that satisfy limits
        coords=str2double(mnibounds);
        indx=find(coords(1) < vox(:, 1) & vox(:, 1) < coords(2));
        indy=find(coords(3) < vox(:, 2) & vox(:, 2) < coords(4));
        indz=find(coords(5) < vox(:, 3) & vox(:, 3) < coords(6));
        indvox=intersect(indz, intersect(indx, indy)); %index for mask in each row of MNI       
        coords_meg = vox(indvox);
        maskname = 'MNImasked';
    otherwise
        disp('using default VOI')
        vox=beam.voxels;
        definput={'-75'; '75'; '-115'; '85'; '-45'; '90'};
        coords=str2double(definput);
        indx=find(coords(1) < vox(:, 1) & vox(:, 1) < coords(2));
        indy=find(coords(3) < vox(:, 2) & vox(:, 2) < coords(4));
        indz=find(coords(5) < vox(:, 3) & vox(:, 3) < coords(6));
        indvox=intersect(indz, intersect(indx, indy));       
        coords_meg = vox(indvox);
        maskname = 'MNIdefaultmask';
end
if strcmp(masktype, 'MNI') == 0
    global ndefaults
    ndefaults.mni2tal='brett';
    coords_mni = nut_tal2mni(coords_tal);
    coords_mri = nut_mni2mri(coords_mni);
    coords_meg = nut_mri2meg(coords_mri);
end
if exist('timepts') ~= 1
    timepts=beam.timepts;
    timewindow=beam.timewindow;
    srate=beam.srate;
    bands=beam.bands;
    s=beam.s;
    voxels=beam.voxels;
    voxelsize=beam.voxelsize;
    coreg=beam.coreg;
    sinfo=beam.sinfo; %doesn't exist in group avg data
    snpm=beam.snpm;
    MNIgoodvoxels=beam.MNIgoodvoxels;
end
%%%%%%% Round Tal-identified coords to match beam.voxels resolution %%%%%
%%%%%%% Then triangulate each Tal coord to closest beam.voxel coord %%%%%
for ii = 1:3
     voxelsrounded(:,ii) = voxelsize(ii)*round(coords_meg(:,ii)/voxelsize(ii)); % group mni/mri/meg same coordsys
end
uniquevoxels = unique(voxelsrounded,'rows');
[k, d] = dsearchn(uniquevoxels,beam.voxels);
% Value of k(i) is voxel index of beam.voxels that is 
% closest to uniquevoxels(i) with distance d(i) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(s,2)  
    for timewin = 1:size(s{i},2) %iterate time windows, earliest first
        for freqbin = 1:size(s{i},3) %iterate frequency bins, lowest first
            for ii = 1:size(k,1)
                %%%% the labeled voxel is already inflated from 2x2x2 to
                %%%% to 5x5x5 or 8x8x8, so don't allow more label expansion
                if d(ii) <= voxrez            %  lower = stricter inclusion
                    voxel_number(ii) = k(ii); % k(ii) is the index number of the VOI
                elseif d(ii) > voxrez   
                    voxel_number(ii) = 0;
                end
            end
            voxel_mask(1:size(k,1)) = 0;
            for ii = 1:size(k,1)
                if voxel_number(ii) ~= 0
                    voxel_mask(ii) = 1;
                end
            end
            voxel_mask = voxel_mask';
            masked_beam = voxel_mask.*s{i}(:,timewin,freqbin);
            masked_beam_snpm_p_uncorr_pos = voxel_mask.*snpm.p_uncorr_pos(:,timewin,freqbin);
            masked_beam_snpm_p_uncorr_neg = voxel_mask.*snpm.p_uncorr_neg(:,timewin,freqbin);
            index_zero_p_uncorr_pos = find(masked_beam_snpm_p_uncorr_pos == 0);
            index_zero_p_uncorr_neg = find(masked_beam_snpm_p_uncorr_neg == 0);
            masked_beam_snpm_p_uncorr_pos(index_zero_p_uncorr_pos) = 1;
            masked_beam_snpm_p_uncorr_neg(index_zero_p_uncorr_neg) = 1;
            snpm.p_corr_pos(index_zero_p_uncorr_pos,timewin,freqbin)=1; %
            snpm.p_corr_neg(index_zero_p_uncorr_neg,timewin,freqbin)=1; %
            %Create beam file with restricted voxels and uncorr pvalues
            voi_beam{i}(:,timewin,freqbin)=masked_beam;
            voi_mask{i}(:,timewin,freqbin)=voxel_mask;
            if isfield(snpm, 'T')            
                masked_beam_snpm_T = voxel_mask.*snpm.T(:,timewin,freqbin);
                snpm.T(:,timewin,freqbin)=masked_beam_snpm_T;
            elseif isfield(snpm, 'R')
                masked_beam_snpm_T = voxel_mask.*snpm.R(:,timewin,freqbin);
                snpm.R(:,timewin,freqbin)=masked_beam_snpm_T;
            else 
                error('Check that T or R stat values exist in snpm.')
            end
            snpm.p_uncorr_pos(:,timewin,freqbin)=masked_beam_snpm_p_uncorr_pos;
            voi_mask_snpm.p_uncorr_pos(:,timewin,freqbin)=voxel_mask;
            snpm.p_uncorr_neg(:,timewin,freqbin)=masked_beam_snpm_p_uncorr_neg;
            voi_mask_snpm.p_uncorr_neg(:,timewin,freqbin)=voxel_mask;
            voxel_mask=0; % reset mask for next
        end
    end
end
%Create s_beam with restricted voxels and .... define cluster correction
beam.s=voi_beam; %masked s structure with 0 values in out of bounds voxels
s=beam.s;
P=snpm.p_uncorr_pos;
N=snpm.p_uncorr_neg; % masked snpm.p_uncorr with p=1 in masked positions
R=snpm.T;
CP=nut_clusterstats_hard(R, P, beam, clusternum, hardthresh);
CN=nut_clusterstats_hard(R, N, beam, clusternum, hardthresh);
% save clustersize used and fill clustercorr menu item with new P values
snpm.clusterthres_pos=clusternum;
snpm.clusterthres_neg=clusternum;
snpm.p_cluster_corr_pos=CP;
snpm.p_cluster_corr_neg=CN;
clear beam CN CP N P R
if exist('mnibounds', 'var')
    snpm.mnimask=mnibounds;
end
beamname=[sbeamfile(1:length(sbeamfile)-4), '_cluster', num2str(clusternum), '_p', num2str(hardthresh), '_', maskname, '.mat'];
save(beamname, 'timepts','timewindow','bands','srate','voxelsize','voxels','s','coreg','MNIgoodvoxels','snpm','sinfo','voi_mask');

%%%%%%%%% Grabs the coords_tal associated with the label %%%%%%%%%%%%
function [index_anat_str,coords_tal] = subf_find_talDB_anatindex(sbeam,anat_str)
global rivets
for i = 1:5
    for j = 1:1105
        if strcmp(rivets.TalDB.labels{i,j},anat_str) == 1
            k(i,j)=1;
        else
            k(i,j)=0;
        end
    end
end
%Find index number for region match; can only match one string
for i = 1:5
index_anat_str{i}=find(k(i,:)); % Find which cell the matches are in (1 through 5)
end
j=0;
for i = 1:5
    if isempty(index_anat_str{i}) == 0 %Should only be one cell which is nonzero
        index_anat_str_nonempty = index_anat_str{i};
        index_anat_str_nonempty_rownumber = i;
        j=j+1; 
    end
end
texmess=sprintf('Found %i entry matching %s.', j, anat_str); % this should always be 1
disp(texmess);
index_anat_str = index_anat_str_nonempty;
lia=ismember(rivets.TalDB.data,index_anat_str);
voi_index = find(lia); %Find voxel locations (coords) for anat_str
voi_tal_coords = rivets.TalDB.coords(voi_index,:);
coords_tal = voi_tal_coords;
%%%%%%%%%%%%%%%%%
