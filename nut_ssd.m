function [outdata,W,A,K] = nut_ssd(data,lo,hi,noiseband,srate,numcomps,flag)
% NUT_SSD performs spatio-spectral decomposition
%
% [outdata,W,A,K] = nut_ssd(data,lo,hi,noiseband,srate,numcomponent,flag)
%
%   data            input data (time*channel*epochs)
%   lo, hi          cutoff frequencies in Hz for bandpass filter
%   noiseband       typically 3 Hz around the passband defined by lo/hi
%   srate           sampling rate of data in Hz
%   numcomponent    (optional) number of SSD components to use. By default,
%                   Haufe's heuristic is used to determine automatically.
%   flag            (optional) by default, low level factorization is given
%                   as outdata (i.e., the output data has the same dimension
%                   as input, but noise removed). If flag is set to
%                   'dimred', the SSD components are given instead.
%
% References:
% - Nikulin et al. NeuroImage 2011; 55: 1528
% - Haufe et al. NeuroImage 2014; 101: 583

dim = size(data);
if length(dim)<3, dim(3)=1; end

S = nut_filter2(data,'butter','bp',4,lo,hi,srate,1);
if (lo-noiseband>0)
    N = nut_filter2(data,'butter','bp',4,lo-noiseband,hi+noiseband,srate,1);
else
    N = nut_filter2(data,'butter','',4,0,hi+noiseband,srate,1);
end
% try
%     N = nut_filter2(N,'butter','notch',4,lo-1,hi+1,srate,1);
% catch   % out of memory
%     disp 'out of memory'
    N = N-S;
% end

C = nut_cov(S,0);
K = nut_cov(N,0); clear N

[W,D] = eig(C,K);                     % equation 12 of Nikulin et al.
if (nargin<6 || isempty(numcomps))
    E = diag(D);
    Q = prctile(E,[25 75]);
    thres = Q(2) + (Q(2) -  Q(1));    % equation 10 of Haufe et al.
    numcomps = sum(E>thres)
end
W = W(:,dim(2):-1:dim(2)-numcomps+1);
A = K * W;                            % equation 21 of Nikulin et al.

if nargin>6 && strcmpi(flag,'dimred')
    outdata = zeros(dim(1),numcomps,dim(3));
    for k=1:dim(3)
        outdata(:,:,k) = S(:,:,k) * W;    % equation 1 in Haufe et al.
    end
else
    outdata = zeros(dim);
    for k=1:dim(3)
        outdata(:,:,k) = (A * W.' * S(:,:,k).').';  % equation 2 in Haufe et al.
    end
end