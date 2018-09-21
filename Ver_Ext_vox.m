function [fv] = Ver_Ext_vox(dat,NPoints);
% (subfunction)
% Extract a mesh from a binary mask using matlab script isosurface.
%
% Author: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1 2006
% Version $1.0

%=========================Main program====================================%
% dim = V.dim(1:3)';
% siz = sqrt(sum(V.mat(1:3,1:3).^2));
if nargin == 1
    NPoints = 10000;
end
S = size(dat);
ind = find(dat);
[x,y,z] = ind2sub(S,ind);clear ind;
xc = nonzeros(single([min(x)-1:max(x)+1]'));clear x;
yc = nonzeros(single([min(y)-1:max(y)+1]'));clear y;
zc = nonzeros(single([min(z)-1:max(z)+1]'));clear z;
[meshx,meshy,meshz]=meshgrid(xc,yc,zc);
dat=permute(dat,[2 1 3]);
tic;fv=isosurface(meshx,meshy,meshz,dat(min(yc):max(yc),min(xc):max(xc),min(zc):max(zc)),0,'verbose');toc;
clear meshx meshy meshz;
[fa,indtt] =sort(fv.faces');
[ford,indt2] =sortrows(fa');ford = double(ford);
temp = ford(1:end-1,:) - ford(2:end,:);temp = sum(abs(temp'));
fnd = find(temp == 0);
fv.faces(indt2(fnd),:)=[];
fv.faces = uint32(fv.faces);
fv.vertices = double(fv.vertices);
if (NPoints<=length(fv.vertices))&(NPoints ~=0)
    disp('   ');
    disp('Reducing Surface...');
    factor = 1/(length(fv.vertices)/NPoints);
    tic;fv=reducepatch(fv,factor);toc;
elseif (NPoints>length(fv.vertices))&(NPoints ~=0)
    disp('The number of points is greater than the number of vertex in the original surface. No reduction');
end
[fa,indtt] =sort(fv.faces');
[ford,indt2] =sortrows(fa');ford = double(ford);
temp = ford(1:end-1,:) - ford(2:end,:);temp = sum(abs(temp'));
fnd = find(temp == 0);
fv.faces(indt2(fnd),:)=[];
fv.faces = uint32(fv.faces);
fv.vertices = double(fv.vertices)-1;
h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);
norma = sqrt(sum((normals').^2));
%fv.normals = normals./repmat(norma',[1 3]);
%========================End of main program==============================%
return;