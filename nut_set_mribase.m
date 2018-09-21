function nut_set_mribase(pat)
% NUT_SET_MRIBASE  Sets new MRI file paths.
% When changing to computers or modifying path names, NUTMEG can produce errors
% when sessions or s_beam files are opened, if it cannot find the MRIs
% anymore. You can resolve this with this function.
%
% nut_set_mribase yourprojectpath
% 
% yourprojectpath   must be a string indicating the complete path to the directory
%                   containing all MRI data. If you have MRIs in different
%                   subdirectories, indicate the parent directory which is
%                   common to all files. E.g.:
%                   'G:\Research\FantasticProject'
%                   The MRIs can then be in subdirectories, e.g.,
%                   'G:\Research\FantasticProject\Subject01\MRI\T1.nii'
%                   'G:\Research\FantasticProject\Subject02\MRI\T1.nii'
%                   If only the path to the parent directory, but not the
%                   name of the subject specific subdirectory changed,
%                   NUTMEG will figure out the new position automatically.

global ndefaults
if ~isstruct(ndefaults)
    nut_defaults;
end
ndefaults.mribase = pat;
