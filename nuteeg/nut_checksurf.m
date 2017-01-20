function vol = nut_checksurf(vol)


for ii = 1:length(vol.bnd)-1
    % Despite what the instructions for surfboolean says, surfaces should
    % be ordered from inside-out!!
    [newnode, newelem] = surfboolean(vol.bnd(ii+1).vertices,vol.bnd(ii+1).faces,'decouple',vol.bnd(ii).vertices,vol.bnd(ii).faces);
    vol.bnd(ii+1).faces = newelem(newelem(:,4)==2,1:3) - size(vol.bnd(ii+1).vertices,1);
    vol.bnd(ii+1).vertices = newnode(newnode(:,4)==2,1:3);
end

for ii = 1:length(vol.bnd)
    [vol.bnd(ii).vertices, vol.bnd(ii).faces] = surfreorient(vol.bnd(ii).vertices, vol.bnd(ii).faces);
    vol.bnd(ii).faces = vol.bnd(ii).faces(:,[3 2 1]);
end



