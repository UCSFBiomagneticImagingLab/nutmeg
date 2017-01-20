% warp_via_NFT() - Main warping function
%
% Usage:
%   >> [Ptm, ind, Cscalp_w,Cskull_w,CCSF_w,W,A,e,LMm2] =
%   warping_main_function(Cscalp, Escalp, Cskull, CCSF,LMm, Fm, Fd, pos);
%
% Inputs:
%   Cscalp, Escalp - scalp mesh
%   Cskull         - coordinates of skull mesh
%   CCSF           - coordinates of CSF mesh
%   LMm            - landmarks on the template mesh
%   Fm             - fiducials of the template mesh
%   Fd             - fiducials from the digitizer data
%   pos            - electrode locations
%
% Outputs:
%   Ptm      - electrode locations on the mesh
%   ind      - index of the electrodes on the mesh
%   Cscalp_w - warped scalp mesh coordinates
%   Cskull_w - warped skull mesh coordinates
%   CCSF_w   - warped CSF mesh coordinates
%   W, A, e  - warping transform parameters
%   LMm2     - warped landmarks 
%
%
% Author: Zeynep Akalin Acar, SCCN, 2008

% Copyright (C) 2007 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vol = warp_via_NFT(Fd, elo)

elo = [Fd;elo];
%[d, elo] = warping_distafterwarping([0 0 0 0 0 90], elo, elo); % arrange orientation ??? check!

load Warping_MNIdata4L.mat

a1 = max(elo) - min(elo);
a2 = max(Cscalp(:,2:4)) - min(Cscalp(:,2:4));
rat = mean(a2./a1);
% make the same scale with the mesh
if rat>500
    elo = elo * 1000; 
elseif rat>50
    elo = elo * 100; 
elseif rat>5
    elo = elo * 10;
end

[pos, Fd] = initial_registration(elo, Fd, Cscalp, Fm);

%[elox, dm] = warping_distmeshafterwarping([0 0 0 0 0 0], pos, Cscalp, Escalp);
%mdm = median(dm); sdm = std(dm);
index_kdm = 1:size(pos,1); %find((dm < 2*mdm)); % kdm gives the index of the electrodes close to the scalp
%index = [1:length(pos)]; rejected = setdiff(index,index_kdm)

figure;
eeglab_plotmesh(Escalp(:,2:4), Cscalp(:,2:4),[],1); hold on; axis on;
axis([-100 100 -200 100 -100 150])
view(165, 10)
plot3(pos(:,1), pos(:,2), pos(:,3), 'b.')
%plot3(pos(rejected,1),pos(rejected,2),pos(rejected,3),'ro')
% index_kdm is the index of the electrodes that are close to the scalp, to
% do the warping

% find warping parameters
[W,A,e,LMm2, Pt, back] = find_warping(Cscalp, Escalp, LMm, Fm, Fd, pos,index_kdm);

% warp the mesh
Cscalp_w = warped_mesh(Cscalp,A,W,LMm2);
Cskull_w = warped_mesh(Cskull,A,W,LMm2);
CCSF_w = warped_mesh(CCSF,A,W,LMm2);
Cbrain_w = warped_mesh(Cbrain,A,W,LMm2);

%%%%%%%
of = [pwd filesep];

% generate a file for StepSc2.txt for final improvement
f=fopen(sprintf('%sStepSc2.txt',of), 'w');
fprintf(f, 'correct 2\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 5\n');
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'improve 2 0.3 0.2\n');
fprintf(f, 'correct 5\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

conf = nft_get_config;
Mesh_WriteSMF(of, 'temp.smf', Cscalp_w, Escalp);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[Cscalp_w,Escalp] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

Mesh_WriteSMF(of, 'temp.smf', Cskull_w, Eskull);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[Cskull_w,Eskull] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

Mesh_WriteSMF(of, 'temp.smf', CCSF_w, ECSF);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[CCSF_w,ECSF] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

Mesh_WriteSMF(of, 'temp.smf', Cbrain_w, Ebrain);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[Cbrain_w,Ebrain] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

[so2, k1,k2] = mesh_check_intersection(Cbrain_w(:,2:4), CCSF_w, ECSF);
Cbrain_w(:,2:4) = so2;
[so2, k1,k2] = mesh_check_intersection(CCSF_w(:,2:4), Cskull_w, Eskull);
CCSF_w(:,2:4) = so2;
[so2, k1,k2] = mesh_check_intersection(Cskull_w(:,2:4), Cscalp_w, Escalp);
Cskull_w(:,2:4) = so2;


%%%%%%%
%Output

vol.bnd(1).vertices = Cscalp_w(:,2:4);
vol.bnd(1).faces = Escalp(:,2:4);
vol.bnd(2).vertices = Cskull_w(:,2:4);
vol.bnd(2).faces = Eskull(:,2:4);
vol.bnd(3).vertices = Cbrain_w(:,2:4);
vol.bnd(3).faces = Ebrain(:,2:4);

% show result
% figure;
% eeglab_plotmesh(Escalp(:,2:4), Cscalp_w(:,2:4),[],1); hold; axis on;
% plot3(Ptm(:,1), Ptm(:,2), Ptm(:,3), 'b.')
% axis([-100 100 -200 100 -100 150])
% view(165, 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mesh_WriteSMF(of, name, Coord, Elem);
nnp = size(Coord,1); 
nel = size(Elem,1);
fid = fopen([of name], 'w');
fprintf(fid,'v %f %f %f \r\n',Coord(:,2:4)');
fprintf(fid,'t %d %d %d \r\n',Elem(:,2:4)');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Coordw = warped_mesh(Coord,A,W,p);

r = Coord(:,2:4);
rw = warp_lm(r,A,W,p) + r;
Coordw = Coord;
Coordw(:,2:4) = rw;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rw] = warp_lm(r,A,W,p)
rw = r * A(1:3,1:3) + repmat(A(4,:), size(r,1), 1);
for i = 1 : size(p,1)
    U = sqrt(sum((r - repmat(p(i,:), size(r,1),1)).^2, 2));  
    rw = rw + U * W(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, A, e, LMm2, Ptbu, back] = find_warping(Coord, Elem, LMm, Fm, Fd, pos,index_kdm);
% Coord, Elem is the mesh that will be warped 
% LMm : Landmarks on the mesh (n2 x 3)
% Fm : fiducials on the mesh (n1 x 3)
% Fd : fiducials on the digitizer data (n1 x 3)
% pos : digitizer data (ne x 3)

ne = size(pos, 1); % number of electrodes
n1 = size(Fm, 1);  % number of fiducials
n2 = size(LMm, 1); % number of landmarks

options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun',1e-6);
Xo = [0 0 0 0 0 0];

% find translation and rotation for fiducials
X = fminsearch(@(X) funrstPP(X, Fd, Fm), Xo, options);
X = fminsearch(@(X) funrstPP(X, Fd, Fm), X, options);

% find translated digitizer fiducials
[d, Fdt] = warping_distafterwarping(X, Fd, Fm); 

% find translated and rotated digitizer locations
[d, Pt] = warping_distafterwarping(X, pos, ones(ne,3));

% Ptm are the rotated and translated digitizer locations 
% and moved to the closest point on the mesh 

Ptbu = Pt;
Pt = Pt(index_kdm,:);

[Ptm, dmi] = funrstp2(X, pos(index_kdm,:), Coord, Elem);
%[Ptm, dmi] = funrstp2(X, pos, Coord, Elem);

% find the index and distance between minimum distance Ptm and LMm  
for i = 1 : n2
    K = Ptm - ones(size(Ptm,1),1) * LMm(i,:);
    L = sum(K.*K,2);
    [k,l] = min(L);
    indm(i) = l;
    minm(i) = k;
end

% find the landmarks on the translated and rotated digitizer locations
LMd = Pt(indm,:); % on digitizer

LMm2 = Ptm(indm,:); % on model

% calculate the warping transformation
[W,A,e] = warp_transform(LMm2, LMd);

% warp_back parameters
[W2,A2,e2] = warp_transform(LMd, LMm2);
back.W = W2;
back.A = A2;
back.e = e2;
back.LMd = LMd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, A, e] = warp_transform(p, q)
% [W,A,e]=warp_transform(p,q)
% calculates nonlinear transformatin coefficents see Ermer's Thesis
% p... = Landmarks in system 1
% q... = landmarks in system 2
% e =warp energy
K = zeros(size(p,1),size(p,1));
for i = 1:size(K,1)
    for j = 1:size(K,2)
        K(i,j) = norm(p(i,:)-p(j,:));    
    end
end
P = [p ones(size(p,1),1)];
L = [K P;P' zeros(4,4)];
D = [q-p;zeros(4,3)];
H = L\D;
W = H(1:size(p,1),:);
A = H(size(p,1)+1:end,:);
e = sum(diag(W'*K*W));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = funrstPP(X, F, Fe);
% F is the point set of the digitizer
% d is the distance between translated and rotated F and Fe
% X is the vector of translation and rotation parameters

tx = X(1);
ty = X(2);
tz = X(3);
alpx = X(4)*pi/180;
alpy = X(5)*pi/180;
alpz = X(6)*pi/180;

x = F(:,1);
y = F(:,2);
z = F(:,3);

% rotation around x-axis
x1 = x;
y1 = y * cos(alpx) - z * sin(alpx);
z1 = y * sin(alpx) + z * cos(alpx);

% rotation around y-axis
x2 = z1 * sin(alpy) + x1 * cos(alpy);
y2 = y1;
z2 = z1 * cos(alpy) - x1 * sin(alpy);

% rotation around z-axis
x3 = x2 * cos(alpz) - y2 * sin(alpz);
y3 = x2 * sin(alpz) + y2 * cos(alpz);
z3 = z2;

% translation
x4 = x3 + tx;
y4 = y3 + ty;
z4 = z3 + tz;

N = size(Fe,1);

Ma = Fe - [x4 y4 z4];
d = sum(sqrt(sum(Ma.*Ma,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F2, dmi] = funrstp2(X, F, Coord, Elem);
% F : digitizer points
% F2 : are the rotated and translated version of F according to X
% they are moved to the closest points on the mesh.

tx = X(1);
ty = X(2);
tz = X(3);
alpx = X(4) * pi / 180;
alpy = X(5) * pi / 180;
alpz = X(6) * pi / 180;

x2 = F(:,1);
y2 = F(:,2);
z2 = F(:,3);

% rotation around x-axis
x3 = x2;
y3 = y2 * cos(alpx) - z2 * sin(alpx);
z3 = y2 * sin(alpx) + z2 * cos(alpx);

% rotation around y-axis
x4 = z3 * sin(alpy) + x3 * cos(alpy);
y4 = y3;
z4 = z3 * cos(alpy) - x3 * sin(alpy);

% rotation around z-axis
x5 = x4 * cos(alpz) - y4 * sin(alpz);
y5 = x4 * sin(alpz) + y4 * cos(alpz);
z5 = z4;

% translation
x1 = x5 + tx;
y1 = y5 + ty;
z1 = z5 + tz;

F2 = [x1 y1 z1];

hh = waitbar(0,'calculating the distance between digitizer and mesh');
for i = 1 : length(F2);
    waitbar(i/length(F2));
    [dm, Pm] = warping_distmeshpoint(F2(i,:), Coord, Elem);
    dmi(i) = dm;    Pmi(i,:) = Pm;
end; 
close(hh);
F2 = Pmi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos, Fd] = initial_registration(elo, F, Cscalp, Fm);

ne = size(elo,1);
[P1e, P2e] = find_new_points_for_reg(elo, F);
[P1m, P2m] = find_new_points_for_reg(Cscalp(:,2:4), Fm);

% find coarse scaling (x) using F2-F3
sx = (F(2,:)-F(3,:))/(Fm(2,:)-Fm(3,:)); sx2=abs(1-sx);
% find coarse scaling (y) using P1-F1
sy = (P1e-F(1,:))/(P1m-Fm(1,:)); sy2=abs(1-sy);
% find coarse scaling (z) using P1-P2
sz = (P1e-P2e)/(P1m-P2m); sz2=abs(1-sz);
% find the one closest to 1, smallest scaling factor
[k,l] = min([sx2 sy2 sz2]);
a = [sx sy sz]; min_sc=a(l);
% after scaling
elo2 = elo;
elo2(:,1) = elo2(:,1) / min_sc;
elo2(:,2) = elo2(:,2) / min_sc;
elo2(:,3) = elo2(:,3) / min_sc;
F2 = elo2(1:3,:);

[P1e2, P2e2] = find_new_points_for_reg(elo2, F2);

% find coarse translation using P1e, P1m
tr = P1e2 - P1m;
elo3 = elo2 - ones(length(elo2),1) * tr;
F3 = F2 - ones(3,1)*tr;

[P1e3, P2e3] = find_new_points_for_reg(elo3, F3);

% find the rotation using fiducials and P2
options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun',1e-6);
Xo = [0 0 0];
X = fminsearch(@(X) funrstPP_rot(X, [F3; P2e3], [Fm;P2m]), Xo, options);
X2=[0 0 0 X];
% find the rotated digitizer locations
[d, elo4] = warping_distafterwarping(X2, elo3, elo3);

[P1x, P2x] = find_new_points_for_reg(elo4, elo4(1:3,:));
%pos = elo4(4:ne,:); % don't take fiducials
pos = elo4; % save with fiducials
Fd = elo4(1:3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P1, P2] = find_new_points_for_reg(elo,F);
% elo is the digitizer location (ne x 3)
% F is the fiducials (3 x 3)
%     1st row nasion
%     2nd row LPA
%     3rd row RPA
% P1 is the mean for the ear fiducials
% P2 is the upper point of the line that is perpendicular to the F1-F2-F3 plane
%       that intersects the digitizer locations
 
ne = length(elo); % number of electrodes

F1 = F(1,:); % nasion
F2 = F(2,:); % LPA
F3 = F(3,:); % RPA


% P1 is the mean point of ear fiducials
P1 = (F2 + F3) / 2;

% plane equation for F1-F2-F3
AB = F2 - F1;
AC = F3 - F1;
n = cross(AB,AC); nz=n(3);
% find the line equation perperdicular to F1-F2-F3 plane r(t)
max_d = max(elo(:,3))/nz*2;
min_d = min(elo(:,3))/nz;
incr = (max(elo(:,3))-min(elo(:,3)))/nz/100;
t = min_d:incr:max_d;
r = ones(length(t),1)*P1 + t'*n; %in terms of t

% find the closest electrode point to r
for i=1:ne
    p1 = elo(i,:);
    M = r - ones(length(t),1)*p1;
    M = sqrt(sum(M.*M,2));
    [k,l] = min(M);
    dis(i,1) = k; % minimum distance
    dis(i,2) = l; % index of r
end

[k,l] = min(dis(:,1));

rm = dis(l,2); %the index of closest point on r
P2 = r(rm,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = funrstPP_rot(X, F, Fe);
% F is the point set of the digitizer
% d is the distance between translated and rotated F and Fe
% X is the vector of rotation parameters

alpx = X(1)*pi/180;
alpy = X(2)*pi/180;
alpz = X(3)*pi/180;

x = F(:,1);
y = F(:,2);
z = F(:,3);

% rotation around x-axis
x1 = x;
y1 = y * cos(alpx) - z * sin(alpx);
z1 = y * sin(alpx) + z * cos(alpx);

% rotation around y-axis
x2 = z1 * sin(alpy) + x1 * cos(alpy);
y2 = y1;
z2 = z1 * cos(alpy) - x1 * sin(alpy);

% rotation around z-axis
x3 = x2 * cos(alpz) - y2 * sin(alpz);
y3 = x2 * sin(alpz) + y2 * cos(alpz);
z3 = z2;

N = size(Fe,1);

Ma = Fe - [x3 y3 z3];
d = sum(sqrt(sum(Ma.*Ma,2)));

