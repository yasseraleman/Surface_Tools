function [Surfa] = Surf_Ext_Corr(Surf);

% Syntax :
% [Surfa] = Surf_Ext_Corr(Surf);
%
% Removes some topological errors from surface files
% 1 Refills deleted faces.
% 2 Remove repeated faces.
% 3 Remove repeated points.
%
% Input Parameters:
% Surf        : Matlab structure variable containing surface information.
%
%
% Output Parameters:
% Surfa        : Corrected Matlab structure variable containing surface.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf
% __________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% July 1st 2008
% Version $1.0

warning off
fclose('all');
% =====================Checking Input Parameters===========================%
if nargin ==0
    errordlg('Surface information is needed');
    return;
end

% =========================================================================%
%
%=========================Main program====================================%
fv.faces = Surf.SurfData.faces;
[fa,indtt] =sort(fv.faces');
[ford,indt2] =sortrows(fa');ford = double(ford);
temp = ford(1:end-1,:) - ford(2:end,:);temp = sum(abs(temp'));
fnd = find(temp == 0);
fv.faces(indt2(fnd),:)=[];
Surf.SurfData.faces = fv.faces;
Npoints = size(Surf.SurfData.vertices,1);
Nfaces = size(Surf.SurfData.faces,1);
[Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
Temp = sum(Tri);
Tri(:,Temp==0) = [];
Surf.Tri = Tri;
ind = find(Surf.Tri(:,2)<=2);
while ~isempty(ind)
    dfac = unique(nonzeros(Surf.Tri(ind,3:end)));
    Surf.SurfData.faces(dfac,:) = [];
    Surf.SurfData.vertices(ind,:) = [];
    if isfield(Surf.SurfData,'VertexNormals')
        Surf.SurfData.VertexNormals(ind,:) = [];
    end
    if isfield(Surf,'Is')
        Surf.Is(ind,:) = [];
    end
    if isfield(Surf.SurfData,'FaceVertexCData')
        Surf.SurfData.FaceVertexCData(ind,:) = [];
    end
    cont = 0;
    Mat = Surf.SurfData.faces;
    for i =1:size(ind,1)
        cont= cont+1;
        dvert = find(Surf.SurfData.faces >ind(i));
        Mat(dvert) = Surf.SurfData.faces(dvert) - cont;
    end
    Surf.SurfData.faces = Mat;
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri;
    ind = find(Surf.Tri(:,2)<=2);
end
try [Surf] = Corr_face(Surf); end;
Surfa = Surf;clear Surf;

%=========================End of Main Program==============================%
return;
%
%=========================Internal Functions==============================%
% Refilling empty faces.
function [Surf] = Corr_face(Surf);
for j = 1:size(Surf.SurfData.vertices,1)
    vert = Surf.SurfData.faces(nonzeros(Surf.Tri(j,3:end)),:);vert = vert(:);vert(vert == j) = [];
    c = accumarray(double(vert),ones(size(vert,1),1));
    ind = find(c == 1);
    if ~isempty(ind)&size(ind,1) ==2
        Surf.SurfData.faces(end+1,:) = [j ind'];
        [Tri] = Vert_Neib(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
        Temp = sum(Tri);
        Tri(:,Temp==0) = [];
        Surf.Tri = Tri;
    elseif ~isempty(ind)&size(ind,1) ==1
        fac = nonzeros(Surf.Tri(j,3:end));
        indtemp = find(Surf.SurfData.faces(fac,:) == ind);
        [rw,col] = ind2sub(size(Surf.SurfData.faces(fac,:)),indtemp);
        Surf.SurfData.faces(fac(rw),:) = [];
        [Tri] = Vert_Neib(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
        Temp = sum(Tri);
        Tri(:,Temp==0) = [];
        Surf.Tri = Tri;
    end
end
return;