function [meg,zs] = nut_check_channel(meg,z_thd,coordfile)

% NUT_CHANNEL_SELECT  detect "bad" channels: channels having either
%   strong deviation of their mean correlation with close channels or 
%   of their variance. 
%
%   Usage:
%   meg =nut_epoch_select(meg,[z_thd],[coordfile])
%
%   meg         meg field of nuts structure (nuts.meg)
%   z_thd       threshold to reject channels that have any properties above the 
%               given standard deviation (optional, default is 3)
%   coordfile   filename of electrode coordinate file. Can be spherical or
%               MNI template (optional).
%
% Output is the meg field with bad channels removed.
%
% Modified after code in EEGLAB's FASTER toolbox.


if nargin<2 || isempty(z_thd)
    z_thd=[3 3]; 
elseif isscalar(z_thd)
    z_thd=z_thd*ones(1,2);
end

if (size(meg.data,2) > length(meg.goodchannels))
    meg.data=meg.data(:,goodchannels,:);
end

[lentrial,numchan,numtrials] = size(meg.data);

list_properties=zeros(numchan,2);

% 1 Mean correlation between each channel and all other channels
if isfinite(z_thd(1))
%     if ~strcmpi(meg.reference,'AVG')
%         oldmeg=meg;
%         meg=nut_eegref(meg,'AVG');
%     end

    % Calculate correlations
    % compute euclidian distance of each channels with all others
    %[x y z] = textread('/Users/amoz/Documents/MATLAB/spm8/EEGtemplates/biosemi128.xyz','%f%f%f%*s','headerlines',1);
    if nargin>2
        elec=nut_read_eeg_coords(2,coordfile);
        c=[elec.x(meg.goodchannels) elec.y(meg.goodchannels) elec.z(meg.goodchannels)];
    elseif isfield(meg,'sensorCoord')
        c=meg.sensorCoord;
        if size(c,1)>size(meg.data,2)
            c=c(meg.goodchannels,:);
        end
    else
        elec=nut_read_eeg_coords(2);
        c=[elec.x(meg.goodchannels) elec.y(meg.goodchannels) elec.z(meg.goodchannels)];
    end
    C=zeros(numchan,numchan);
    for k=1:numchan
    %     euclidian distance
        C(k,:)=nut_rownorm(c-repmat(c(k,:),[numchan 1]));
    end

    mcorrs=zeros(numchan,numtrials);
    for epo=1:numtrials
        for k=1:numchan
        % find closest neighbours of electrode
            nbrs = (C(k,:)>0 & C(k,:)<30);
            % mean correlation with close electrodes
            mcorrs(k,epo) = nanmean(abs(corr(meg.data(:,k,epo),meg.data(:,nbrs,epo))));
        end
    end
    % mean correlation of all epochs
    mcorrs=nanmean(mcorrs,2);

    if ~strcmpi(meg.reference,'AVG')
        closetoref = (C(meg.goodchannels==meg.referenceidx,:)<30);
        % reference electrodes and neighbours to 1
        mcorrs(closetoref)=1;
    end

    list_properties(:,1) = -mcorrs;

%     if exist('oldmeg','var')
%         meg=oldmeg;
%         clear oldmeg
%     end
end
    
% 2 Variance of the channels
if isfinite(z_thd(2))
    vars = zeros(numtrials,numchan);
    for epo=1:numtrials
        v = var(meg.data(:,:,epo));
        v(~isfinite(v))=mean(v(isfinite(v)));
        vars(epo,:)=v;
    end
    vars=mean(vars,1);

    list_properties(:,2) = vars;
end

for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
	list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end

zs=list_properties-repmat(mean(list_properties,1),size(list_properties,1),1);
zs=zs./repmat(std(zs,[],1),size(list_properties,1),1);
zs(isnan(zs))=0;

all_l = zs > repmat(z_thd,[size(zs,1) 1]);
bad_channels_idx = any(all_l,2);

if any(bad_channels_idx)
    lab=meg.sensor_labels(meg.goodchannels);
    disp('The following channels are bad:')
    disp(lab(bad_channels_idx));
else
    disp('No bad channels detected.')
end

meg.data(:,bad_channels_idx,:)=[];
meg.goodchannels(bad_channels_idx)=[];
