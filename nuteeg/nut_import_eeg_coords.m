function nut_import_eeg_coords(coordtype)
% Imports coordinates of EEG sensors and any headshape points. These
% coordinates are adjusted so that the origin is [0 0 0], with the nasion
% being [+X 0 0], lpa being [0 +Y 0], and rpa being [0 -Y 0], as is used in
% the MEG coordinate system. 
%
% Optional input:
%   coordtype       1=digitized electrode coordinates (e.g., with Polhemus),
%                     (default)
%                   2=template coordinates aligned to MNI template brain
%                   3=coordinates previously aligned to individual MRI 
%                     (e.g., in Cartool)
%                   
% The following fields are stored in the nuts structure:
%   nuts.meg.sensorCoord - the coordinates of the sensors in mm [nSensorsx3].
%                          Each sensor corresponds to a sensor label in
%                          nuts.meg.sensor_labels.
%   nuts.meg.nasion - coordinates of the nasion in mm, if available.
%   nuts.meg.lpa - coordinates of the lpa in mm, if available.
%   nuts.meg.rpa - coordinates of the rpa in mm, if available.
%
% @authors Daniel D.E. Wong and Adrian G. Guggisberg


global nuts

if nargin<1, coordtype=1; end
if nargin<2, shortcut=false; end

[elec,scalefactor]=nut_read_eeg_coords(coordtype);

if coordtype==1
    if isempty(elec.lpa)
        [idx,~] = listdlg('ListString',elec.label,'SelectionMode','Single','PromptString','Select LPA');
        elec.lpa = [elec.x(idx) elec.y(idx) elec.z(idx)];
    end
    
    if isempty(elec.rpa)
        [idx,~] = listdlg('ListString',elec.label,'SelectionMode','Single','PromptString','Select RPA');
        elec.rpa = [elec.x(idx) elec.y(idx) elec.z(idx)];
    end
    
    if isempty(elec.nasion)
        [idx,~] = listdlg('ListString',elec.label,'SelectionMode','Single','PromptString','Select Nasion');
        elec.nasion = [elec.x(idx) elec.y(idx) elec.z(idx)];
    end
            
    % Find fiducial plane normal
    pa_vec = elec.lpa-elec.rpa; pa_vec = pa_vec/norm(pa_vec);
    naspa_vec = elec.nasion - elec.rpa; naspa_vec = naspa_vec/norm(naspa_vec);
    fp_norm = cross(pa_vec,naspa_vec); fp_norm = fp_norm/norm(fp_norm);

    % Nasion vector (from nasion to origin)
    nasion_vec = cross(fp_norm,pa_vec);

    % Intersection between nasion and pa lines ([x y z] = t[a b c]+[x0 y0 z0])
    for i = 1:3; if nasion_vec(i); idx0 = i; break; end; end;   % index of non-zero nasion vector element
    for i = 1:3; if i ~= idx0 & pa_vec(i); idx1 = i; break; end; end;   % index of non-zero pa vector element
    t = (nasion_vec(idx1)*(elec.lpa(idx0)-elec.nasion(idx0))/nasion_vec(idx0) + elec.nasion(idx1) - elec.lpa(idx1)) / (pa_vec(idx1) - nasion_vec(idx1)*pa_vec(idx0)/nasion_vec(idx0));
    origin = t*pa_vec + elec.lpa;

    % Translate coordinates so that origin is [0 0 0]
    elec.nasion = elec.nasion - origin;
    elec.lpa = elec.lpa - origin;
    elec.rpa = elec.rpa - origin;
    elec.x = elec.x - origin(1);
    elec.y = elec.y - origin(2);
    elec.z = elec.z - origin(3);
    if ~isempty(elec.hsp); elec.hsp = elec.hsp - repmat(origin,size(elec.hsp,1),1); end;

    % Rotate coordinates so that nasion is [+X 0 0], lpa is [0 +Y 0], and rpa
    % is [0 -Y 0] (i.e. MEG coordinates)

    % Rotate nasion
    if elec.nasion(3)
        th = atan(elec.nasion(3)/elec.nasion(1)); % Y rotation (z->0)
        rot = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
        elec = transform_scan3dd(rot,elec);
    end
    if elec.nasion(2)
        th = atan(elec.nasion(2)/elec.nasion(1)); % Z rotation (y->0)
        rot = [cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
        elec = transform_scan3dd(rot,elec);
    end

    % Rotate PA
    if elec.lpa(3)
        th = atan(elec.lpa(3)/elec.lpa(2));   % X rotation (z->0)
        rot = [1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];
        elec = transform_scan3dd(rot,elec);
    end

    % Ensure orientations are correct
    if elec.nasion(1) < 0    % Ensure correct X orientation
        rot = [cos(pi) sin(pi) 0; -sin(pi) cos(pi) 0; 0 0 1];   % Z rotation 180 deg
        elec = transform_scan3dd(rot,elec);
    end
    if elec.lpa(2) < 0       % Ensure correct Y orientation
        rot = [1 0 0; 0 cos(pi) sin(pi); 0 -sin(pi) cos(pi)];   % X rotation 180 deg
        elec = transform_scan3dd(rot,elec);
    end
    
else
    
    [elec.x,elec.y,elec.z] = nut_mri2meg(elec.x,elec.y,elec.z);
    
end

% Store coordinates

% map sensorCoord entries to channels...
if isfield(nuts,'meg') && isfield(nuts.meg,'sensor_labels')
    indexmap = nan(1,length(nuts.meg.sensor_labels));
    for i=1:length(nuts.meg.sensor_labels)
        tmp = strmatch(upper(nuts.meg.sensor_labels{i}),upper(elec.label),'exact');
        if ~isempty(tmp)
            indexmap(i)=tmp;
        end
    end
    if all(isnan(indexmap))
        disp('Labels in data and electrode location files no not match. Assuming 1:1 correspondance.')
        indexmap = 1:length(nuts.meg.sensor_labels);
    elseif any(isnan(indexmap))
        warndlg(['No location found for the following electrode(s): ' sprintf('%s ',nuts.meg.sensor_labels{isnan(indexmap)})]);
    end
        
else
    warning('Cannot check for correct labels because no sensor labels found.')
    indexmap = 1:length(elec.label);
end
nuts.meg.sensorCoord = zeros(length(indexmap),3);
for i = find(isfinite(indexmap))
    nuts.meg.sensorCoord(i,:) = [elec.x(indexmap(i)), elec.y(indexmap(i)), elec.z(indexmap(i))]*scalefactor;
end

% If we're using VEO as a proxy for inion (since Scan3dd doesn't let me set the inion)
veoidx = strmatch('VEO',upper(elec.label));
if ~isempty(veoidx)
    answer = questdlg('Set inion to VEO?','Import EEG Coordinates','Yes','No','No');
    if strcmp(answer,'Yes')
        nuts.meg.inion = [elec.x(veoidx), elec.y(veoidx), elec.z(veoidx)]*scalefactor;
    end
end

if ~isempty(elec.nasion), nuts.meg.nasion = elec.nasion*scalefactor; end
if ~isempty(elec.lpa), nuts.meg.lpa = elec.lpa*scalefactor; end
if ~isempty(elec.rpa), nuts.meg.rpa = elec.rpa*scalefactor; end
if ~isempty(elec.hsp), nuts.meg.hsp = elec.hsp*scalefactor; end  % Store head shape points for multisphere fitting

nut_show_eegmralign;

%---------------------------------
% Performs a transformation operation on elec vectors.
% @param A the transformation matrix.
% @param elec the elec structure.
% @return the elec structure with transformed vectors.
%
function elec = transform_scan3dd(A,elec)
elec.nasion = (A*elec.nasion')';
elec.lpa = (A*elec.lpa')';
elec.rpa = (A*elec.rpa')';
xyz = (A*[elec.x elec.y elec.z]')';
elec.x = xyz(:,1); elec.y = xyz(:,2); elec.z = xyz(:,3);
if ~isempty(elec.hsp), elec.hsp = (A*elec.hsp')'; end
if ~isempty(elec.ref), elec.ref = (A*elec.ref')'; end