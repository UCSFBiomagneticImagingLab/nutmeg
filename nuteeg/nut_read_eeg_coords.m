function [elec,scalefactor] = nut_read_eeg_coords(coordtype,filename)

% Import sensor locations and headshape
if nargin<2
    if coordtype==2 && strncmp(spm('ver'),'SPM8',4), cwd=pwd; cd([fileparts(which('spm.m')) filesep 'EEGtemplates']); end
    [sensor_filename, sensor_path]=uigetfile( ...
               {'*.bvct;*.bvef;*.dat;*.elp;*.sfp;*.elc;*.xyz' ,'Supported Formats'; ...
                '*.elp'             ,'EMSE Locator'; ...
                '*.bvct;*.bvef'     ,'BrainVision'; ...
                '*.dat'             ,'NeuroScan 3dd ASCII Export'; ...
                '*.sfp'             ,'BESA'; ...
                '*.elc'             ,'ASA'; ...
                '*.xyz'             ,'Cartool electrode coordinates'}, ...
                'Load Sensor Coordinates...');
    if isequal(sensor_filename,0) || isequal(sensor_path,0)
        return
    end
    filename=fullfile(sensor_path,sensor_filename);
    if coordtype==2 && strncmp(spm('ver'),'SPM8',4), cd(cwd), clear cwd, end
end

[dum,dum,ext] = fileparts(filename);
switch lower(ext)
case '.dat'
    elec = elec_load_scan_3ddasc(filename);
    scalefactor = 10;
case '.elp'
    elec = elec_emse2matlab(filename);
    elec.hsp = [];
    scalefactor = 1000;
case {'.bvct' '.bvef'}
    elec = in_channel_brainvision(filename);
    elec.hsp = []; elec.ref=[];
    scalefactor = 1000;
case {'.sfp' '.elc'}
    if ~exist('ft_read_sens.m','file'), errordlg('You need FieldTrip in your Matlab path to import this electrode coordinate format.'), return, end
    elec = ft_read_sens(filename);
    if coordtype==1
        nasidx = [find(~cellfun(@isempty,strfind(lower(elec.label(1:3)),'nas'))) find(~cellfun(@isempty,strfind(lower(elec.label(1:3)),'nz')))];
        if isempty(nasidx), elec.nasion=[]; 
        elseif isfield(elec,'pnt'), elec.nasion = elec.pnt(nasidx,:); elec.pnt(nasidx,:)=[]; elec.label(nasidx)=[]; clear nasidx
        else elec.nasion = elec.chanpos(nasidx,:); elec.chanpos(nasidx,:)=[]; elec.label(nasidx)=[]; clear nasidx
        end
        lpaidx = find(~cellfun(@isempty,strfind(lower(elec.label(1:3)),'lpa')));
        if isempty(lpaidx), elec.lpa=[];
        elseif isfield(elec,'pnt'), elec.lpa = elec.pnt(lpaidx,:); elec.pnt(lpaidx,:)=[]; elec.label(lpaidx)=[]; clear lpaidx
        else elec.lpa = elec.chanpos(lpaidx,:); elec.chanpos(lpaidx,:)=[]; elec.label(lpaidx)=[]; clear lpaidx
        end
        rpaidx = find(~cellfun(@isempty,strfind(lower(elec.label(1:3)),'rpa')));
        if isempty(rpaidx), elec.rpa=[];
        elseif isfield(elec,'pnt'), elec.rpa = elec.pnt(rpaidx,:); elec.pnt(rpaidx,:)=[]; elec.label(rpaidx)=[]; clear rpaidx
        else elec.rpa = elec.chanpos(rpaidx,:); elec.chanpos(rpaidx,:)=[]; elec.label(rpaidx)=[]; clear rpaidx
        end
    else
        elec.nasion = []; elec.lpa=[]; elec.rpa=[];
    end
    if coordtype==2
        if isfield(elec,'pnt'), elec.pnt = nut_mni2mri(elec.pnt);
        else elec.chanpos = nut_mni2mri(elec.chanpos);
        end
    end
    if isfield(elec,'pnt'), elec.x = elec.pnt(:,1); elec.y = elec.pnt(:,2); elec.z = elec.pnt(:,3);
    else elec.x = elec.chanpos(:,1); elec.y = elec.chanpos(:,2); elec.z = elec.chanpos(:,3);
    end
    elec.hsp = []; elec.ref = [];
    if isfield(elec,'pnt'); elec = rmfield(elec,'pnt'); else elec = rmfield(elec,'chanpos'); end;
    scalefactor = 1; % Fieldtrip seems to convert to mm units if necessary.
case '.xyz'
    [elec.x elec.y elec.z elec.label] = textread(filename,'%f%f%f%s','headerlines',1);
    elec.hsp=[]; elec.ref=[]; elec.nasion=[]; elec.lpa=[]; elec.rpa=[];
    if coordtype==2
        xyz = nut_mni2mri([elec.x elec.y elec.z]);
        elec.x=xyz(:,1); elec.y=xyz(:,2); elec.z=xyz(:,3);
    end
    scalefactor = 1; 
end
