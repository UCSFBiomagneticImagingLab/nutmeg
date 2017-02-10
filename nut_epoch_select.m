function data = nut_epoch_select(data,N)

% NUT_EPOCH_SELECT  keep best trials based on: less deviation from channel
% mean, less global variance and less amplitude difference between max and
% min values
%
%   data = nut_epoch_select(data,N)
%       data: nutmeg eeg/meg data matrix (nsamples,nchannels,ntrials)
%       N: Number of desired data trials
%       


dims=size(data);
T = permute(data,[2 1 3]);
T = reshape(T,[dims(2) dims(1)*dims(3)]);
means = mean(T,2)'; %clear T

% 1 Epoch's mean deviation from channel means.
for u = 1:size(data,3)
	list_properties(u,1) = mean(abs(squeeze(mean(data(:,:,u),1)) - means));
end

% 2 Epoch variance
list_properties(:,2) = mean(squeeze(var(data,0,1)));


% 3 Max amplitude difference
%for t = 1:size(data,2)
%     for u = 1:size(data,3)
ampdiffs(:,:) = max(data(:,:,:),[],1) - min(data(:,:,:),[],1);
%     end
%end
list_properties(:,3) = mean(ampdiffs,1);

for v = 1:size(list_properties,2)
	list_properties(:,v) = list_properties(:,v) - median(list_properties(:,v));
end

%     rejection_options.measure=ones(1,size(list_properties,2));
%     rejection_options.z=std_thd*ones(1,size(list_properties,2));

% rejection_options.measure=logical(rejection_options.measure);
zs=list_properties-repmat(mean(list_properties,1),size(list_properties,1),1);
zs=zs./repmat(std(zs,[],1),size(list_properties,1),1);
zs(isnan(zs))=0;
[~,idx]=sort(mean(zs,2));
data = data(:,:,sort(idx(1:N)));


