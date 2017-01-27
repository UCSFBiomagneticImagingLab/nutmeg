function nut_opensession(savedfile,partial)
% NUT_OPENSESSION  opens NUTMEG session file.

global nuts

if nargin<2, partial=false; end

figh = findobj('tag','nutmegfig');
isgui = ~isempty(figh);

if ~exist('savedfile','var')
    [nutfilename, nutpathname]=uigetfile('*.mat','Open NUTMEG session...');
    if isequal(nutfilename,0)|isequal(nutpathname,0)
        return;
    end 
    savedfile=fullfile(nutpathname,nutfilename);
end

if partial
    warning('off','MATLAB:load:variableNotFound');
    nuts=load(savedfile,'coreg','voxels','voxelsize','Lp','voxor');
    warning('on','MATLAB:load:variableNotFound');
else
    nuts=load(savedfile);
end
if isfield(nuts,'nuts')
    nuts=nuts.nuts;
end

if isgui
    nuts.fig = figh; % update nutmeg fig handle
elseif isfield(nuts,'fig')
    nuts=rmfield(nuts,'fig'); 
end  

if(isfield(nuts,'megfile'))
    set(findobj('Tag','nut_megfile'),'String',nuts.meg.filename);
end

if( isfield(nuts,'preprocessing') && isfield(nuts.preprocessing,'beamformertype') && (nuts.preprocessing.beamformertype == 1) )  % legacy support... old way had numbers instead of strings
    nuts.preprocessing.beamformertype = 'Default Beamformer';
end

if (isfield(nuts,'coreg') && isfield(nuts.coreg,'mripath') && exist(nuts.coreg.mripath,'file')~=2 )
    nuts.coreg = nut_correct_mripath(nuts.coreg);
end

if isgui
    nut_refresh_image;  % load image into SPM
    nut_enabler;
    % nut_defaults;
end
