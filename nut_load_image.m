function nut_load_image(headname)
% NUT_LOAD_IMAGE

global st nuts
if exist('headname','var')
    nuts.coreg.mripath=headname;
else
	[mrifile, mripath]=uigetfile('*.img;*.nii','Please Select Image Volume...');
	if isequal(mrifile,0)|isequal(nuts.coreg.mripath,0)
        return;
	end
	nuts.coreg.mripath = fullfile(mripath,mrifile);
end
disp(['File: ' nuts.coreg.mripath]);
if length(nuts.coreg.mripath)>35
    imagestring=['...' nuts.coreg.mripath(end-30:end)];
    set(findobj('Tag','nut_ImageNameText'),'String',imagestring);
else
    set(findobj('Tag','nut_ImageNameText'),'String',nuts.coreg.mripath);
end

% fiducial and lsc pts should be cleared and tfm reset if loading new MRI
if(isfield(nuts,'meg') && isfield(nuts.meg,'lsc'))
    nuts.meg=rmfield(nuts.meg,'lsc');
end
if isfield(nuts.coreg,'fiducials_mri_mm')
    nuts.coreg=rmfield(nuts.coreg,'fiducials_mri_mm');
end
nuts.coreg.meg2mri_tfm=eye(4);

% loading image
spm_image('init',nuts.coreg.mripath);
nut_spmfig_setup;

if(strcmp(spm('ver'),'SPM2'))
    type = st.vols{1}.dim(4);
    if(st.vols{1}.private.hdr.dime.dim(1) < 1 | st.vols{1}.private.hdr.dime.dim(1) > 15)
        warning('looks like you''s guys got you''s endian all flipped around n'' whatnot');
    end;
else
    type = st.vols{1}.dt(1);
end


nuts.coreg.orientation = get(findobj('Tag','nut_orientation_menu'),'Value');
%set(findobj('Tag','nut_close_coreg'),'Enable','on')
if(~isempty(gcbo))  % only run nut_coreg_enabler if called from Coreg GUI
    nut_coreg_enabler
end
