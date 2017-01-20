function bad_channels_idx = nut_channel_select(data,z_thd)

% NUT_CHANNEL_SELECT  detect "bad" channels: flat channels or channels having either
%   strong deviation of their mean correlation with close channels or 
%   of their variance. !!! Specific for biosemi !!
% 
%
%   nut_epoch_select(data,z_thd)
%       data: nutmeg eeg/meg data matrix (nsamples,nchannels, ntrials)
%       z_thd: threshold to reject channels that have 
%       any properties above the given standard deviation
%       bad_channels_idx: logical vector indicating the rejected channels as true



measure = 1;
numchan = size(data,2);
numtrials = size(data,3);
lentrial = size(data,1);
% TEMPORAL PROPERTIES

% 1 Mean correlation between each channel and all other channels

eeg_chans = 1:numchan;




% Calculate correlations
% compute euclidian distance of each channels with all others
[x y z] = textread('/Users/amoz/Documents/MATLAB/spm8/EEGtemplates/biosemi128.xyz','%f%f%f%*s','headerlines',1);
c=[x y z];
for k=1:numchan
%     euclidian distance
    C(k,:)=nut_rownorm(c-repmat(c(k,:),[numchan 1]));
end

mcorrs=zeros(numchan,numtrials);
for epo=1:numtrials
%     ignore first channel since it is the reference
    for k=2:numchan
    % find closest neighbours of electrode
        nbrs = (C(k,:)>0 & C(k,:)<30);
        % mean correlation with close electrodes
        mcorrs(k,epo) = nanmean(abs(corr(data(:,k,epo),data(:,nbrs,epo))));
    end
end
% mean correlation of all epochs
mcorrs=nanmean(mcorrs,2);
closetoref = (C(1,:)<30);
% reference electrodes and neighbours to 1
mcorrs(closetoref)=1;

% calc_indices=setdiff(eeg_chans,ignore);
% corrs = abs(corrcoef(data(:,calc_indices)));
% mcorrs=zeros(size(eeg_chans));
% for u=1:length(calc_indices)
%     mcorrs(calc_indices(u))=mean(corrs(u,:));
% end
% mcorrs(ignore)=mean(mcorrs(calc_indices));

list_properties(:,measure) = -mcorrs

measure = measure + 1;



% if length(size(data)) > 2
%     %%%% concatenate all trials 
%     data = permute(data,[2,1,3]);
%     data = reshape(data, size(data,1),size(data,2)*size(data,3));
%     data=data';
% end
% 3 Variance of the channels
for epo=1:numtrials
    v = var(data(:,:,epo));
    v(~isfinite(v))=mean(v(isfinite(v)));
    vars(epo,:)=v;
end
vars=mean(vars,1);

list_properties(:,measure) = vars;

measure = measure + 1;

%%% 4 Hurst exponent
% for u=1:length(eeg_chans)
%     list_properties(u,measure) = -hurst_exponent(data(:,eeg_chans(u))');
% end

for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
	list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end





zs=list_properties-repmat(mean(list_properties,1),size(list_properties,1),1);
zs=zs./repmat(std(zs,[],1),size(list_properties,1),1);
zs(isnan(zs))=0;



all_l = zs > repmat(z_thd,size(list_properties));


bad_channels_idx = any(all_l,2);


