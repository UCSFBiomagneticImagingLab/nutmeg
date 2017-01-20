function coreg = nut_correct_mripath(coreg,newbase)
% NUT_CORRECT_MRIPATH   tries to update MRI file location.
%
% coreg = nut_correct_mripath(coreg)
%
% 2015 Adrian G. Guggisberg

global ndefaults

if exist(coreg.mripath,'file'), return, end

fprintf('MRI path "%s" does not exist. Trying to find correct file location.\n',coreg.mripath)

if nargin<2
    if isfield(ndefaults,'mribase')
        if exist(ndefaults.mribase,'dir')~=7, error('The MRI base directory you specified does not exist.'), end
        newbase = ndefaults.mribase;    % MRI base set with nut_set_mribase
    else
        newbase = pwd;   % Trying to find it close to current directory.       
    end
end

% Check if MRI is SPM template. In this case we should find it easily.
[dum,fi,ext]=fileparts(coreg.mripath);
fi=[fi ext];
if strcmp(fi,'avg152T1.nii')
    whichdir=which(fi);
    if ~isempty(whichdir)
        coreg.mripath=whichdir;
        coreg.norm_mripath=whichdir;
        fprintf('Updated successfully. New path is %s.\n',coreg.mripath)
        return
    end
end

% Create cell array with all potential new directories
lim=find(newbase==filesep);
L=[[1;lim'+1] [lim'-1;length(newbase)]];
for k=1:size(L,1)
    D{k}=newbase(L(k,1):L(k,2));
end

[p,F,ext] = fileparts(coreg.mripath);
if ~isempty(p)
% Create cell array with directories of wrong path
    lim=find(p==filesep);
    L=[[1;lim'+1] [lim'-1;length(p)]];
    for k=1:size(L,1)
        E{k}=p(L(k,1):L(k,2));
    end

% Try to figure out new path
    comp = '';
    k=length(E);    
    while( exist(comp,'dir')~=7 && k>0 )
        C = [D E(k:end)];
        comp = fullfile(C{:});
        k = k-1;
    end
end
if isempty(p) || (k<1)
    C=D;
    comp=newbase;
end

% Find MRI in one of the directories of the new path
P = '';
j = length(C);
while( exist(P,'file')~=2 && j>0 )
    comp = fullfile(C{1:j});
    P = fullfile(comp,[F ext]);
    j=j-1;    
end
if j>0
    %P = strrep(P,[filesep filesep],filesep);
    coreg.mripath=P;
    if isfield(coreg,'norm_mripath')
        [lost,F,ext]=fileparts(coreg.norm_mripath);
        coreg.norm_mripath = fullfile(comp,[F ext]);
    end
    if isfield(coreg,'volfile')
        coreg=rmfield(coreg,'volfile');
    end
    fprintf('Updated successfully. New path is %s.\n',coreg.mripath)
else
    disp('Unable to update MRI path. You can set new path with NUT_SET_MRIBASE.')        
end

