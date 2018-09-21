function CM_Atrophy_Modeler
%%% Anotation average to labels

% mri_annotation2label --subject COMP_avesubj_thmax05 --hemi lh --labelbase $SUBJECTS_DIR/COMP_avesubj_thmax05/label/aparc-lh
% %%% Label average to native
% mri_label2label --srclabel /media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/label/aparc-lh-000.label --srcsubject COMP_avesubj_thmax05 --trglabel /media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_00002__MAD-20090606-T1-MAD-1/label/Testtttttt.label --trgsubject COMP_00002__MAD-20090606-T1-MAD-1 --regmethod surface --hemi lh
% fid = fopen('Testtttttt.label','r');
% lin = fgetl(fid);
% Npoints = str2num(fgetl(fid));
% Temp = reshape(textread('Testtttttt.label','%f',Npoints*5,'headerlines',2),5,Npoints)';
alpha = 0.33;
beta = 0.33;
scf = 1;
porc = 20;
direc = -1;
scalef = 2;
N = 10;
if isunix
    dira = '/media/COSAS/Test/freesurfer/';
    [aparc,ctab] = read_cfiles([dira 'ch2/label/lh.aparc.annot']);
    [thickness,ctab] = read_cfiles([ dira 'ch2/surf/lh.thickness']);
    [Sid, Sname, colx, coly, colz, offset] = textread([ dira 'ch2/label/aparc.annot.ctab'],'%f%s%u%u%u%u',36);
    [vertices, faces] = freesurfer_read_surf([dira 'ch2/surf/lh.white']);Surfw.SurfData.vertices = vertices;Surfw.SurfData.faces = faces;Surfw.Name = 'White';
    [vertices, faces] = freesurfer_read_surf([dira 'ch2/surf/lh.pial']);Surfp.SurfData.vertices = vertices;Surfp.SurfData.faces = faces;Surfp.Name = 'Pial';
else

    dira = 'E:\Test\freesurfer\'
    [aparc,ctab] = read_cfiles([dira 'ch2\label\lh.aparc.annot']);
    [thickness,ctab] = read_cfiles([ dira 'ch2\surf\lh.thickness']);
    [Sid, Sname, colx, coly, colz, offset] = textread([ dira 'ch2\label\aparc.annot.ctab'],'%f%s%u%u%u%u',36);

    [vertices, faces] = freesurfer_read_surf([dira 'ch2\surf\lh.white']);Surfw.SurfData.vertices = vertices;Surfw.SurfData.faces = faces;Surfw.Name = 'White';
    [vertices, faces] = freesurfer_read_surf([dira 'ch2\surf\lh.pial']);Surfp.SurfData.vertices = vertices;Surfp.SurfData.faces = faces;Surfp.Name = 'Pial';

end
labels = colx + coly.*2^8 + colz*2^16 + offset*2^24;%labels(5)=[];

%[OutFiles, SurfF] = Exp_Surf([dira 'COMP_00002__MAD-20090606-T1-MAD-1/surf/lh.white'], '0', '','', 'imp','n');Surfw = SurfF{1};
%[OutFiles, SurfF] = Exp_Surf([dira 'COMP_00002__MAD-20090606-T1-MAD-1/surf/lh.pial'], '0', '','', 'imp','n');Surfp = SurfF{1};
%ind = [1:20000];
ind = find(aparc==labels(2));
a = ismember(Surfw.SurfData.faces,ind);ind2 = find(sum(a')==0);
Surfwt =Surfw; Surfwt.SurfData.faces(ind2,:) = [];        Surfpt =Surfp; Surfpt.SurfData.faces(ind2,:) = [];
indf = unique(Surfwt.SurfData.faces(:));
%%%%%%%%%%%%%%%%%%%%%%%Corriegiendo Superficies para quedarme solo
%%%%%%%%%%%%%%%%%%%%%%%con la parte que me interesa
Surfwtemp = Surfwt;Surfptemp = Surfpt;
Surfwtemp.SurfData.vertices = zeros(size(indf,1),3);
%Surfwtemp.Tri = zeros(size(indf,1),size(Surfwt.Tri,2));
Surfptemp.SurfData.vertices = zeros(size(indf,1),3);
Surfwtemp.SurfData.faces = 0*Surfwt.SurfData.faces;
for i =1:size(indf,1)
    Surfwtemp.SurfData.vertices(i,:) = Surfwt.SurfData.vertices(indf(i),:);
    Surfptemp.SurfData.vertices(i,:) = Surfpt.SurfData.vertices(indf(i),:);
    %Surfwtemp.Tri(i,:) = Surfwt.Tri(indf(i),:);
    indn = find(Surfwt.SurfData.faces ==indf(i));Surfwtemp.SurfData.faces(indn) = i;
end
[Tri] = Vert_Neib(double(Surfwtemp.SurfData.faces),size(Surfwtemp.SurfData.vertices,1),size(Surfwtemp.SurfData.faces,1));
Temp = sum(Tri);
Tri(:,Temp==0) = [];
Surfptemp.Tri = Tri;Surfwtemp.Tri = Tri;
Nlevels = 6;
%Vertnlevels = Neib_level(Surfwtemp,Nlevels);
%Surfwtemp.Is(Vertnlevels(:,1),1) =Vertnlevels(:,2);

%Surfptemp.Tri = Surfwtemp.Tri;
Surfptemp.SurfData.faces = Surfwtemp.SurfData.faces; %  Surfwtemp.SurfData.faces(20:end,:) = []; Surfptemp.SurfData.faces(20:end,:) = [];
Surfa{1,1} = Surfptemp;Surfa{2,1} = Surfwtemp;

%%%%%%%%%%%%%%%%%%%Fast Marching



% vertex = Surfptemp.SurfData.vertices;
% faces = Surfptemp.SurfData.faces;
% 
% 
% options.nb_iter_max = Inf;
% 
% start_points = 23 ;
% [D,S,Q] = perform_fast_marching_mesh(vertex, faces, start_points, options);
% 
% npaths = 10; nverts = size(vertex,2);
% % select some points that are far enough from the starting point
% [tmp,I] = sort( D(:) ); I = I(end:-1:1); I = I(1:round(nverts*1));
% end_points = floor( rand(npaths,1)*(length(I)-1) )+1;
% end_points = I(end_points);
% % precompute some usefull information about the mesh
% options.v2v = compute_vertex_ring(faces);
% options.e2f = compute_edge_face_ring(faces);
% % extract the geodesics
% options.method = 'discrete';
% options.verb = 0;
% paths = compute_geodesic_mesh(D,vertex,faces, end_points, options);



%Ns = max(unique(Surfwtemp.Is)); 
%vf = sigmoid([0 Ns/2 Ns/4],double(Surfwtemp.Is));


%%
xw = Surfwtemp.SurfData.vertices(:,1); yw = Surfwtemp.SurfData.vertices(:,2); zw = Surfwtemp.SurfData.vertices(:,3);distw=sqrt(xw.^2+yw.^2+zw.^2);
xp = Surfptemp.SurfData.vertices(:,1); yp = Surfptemp.SurfData.vertices(:,2); zp = Surfptemp.SurfData.vertices(:,3);distp=sqrt(xp.^2+yp.^2+zp.^2);
norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);
u = (xp-xw)./norms; v = (yp-yw)./norms; w = (zp-zw)./norms;

%%%%%%%%%Direccion de la lesion


 vf =ones(size(norms));
 Surfo = Surfa;
 if direc == 1
     [Ca] = [xw yw zw]+[repmat(norms.*vf*porc/100,[ 1 3]).*[u v w]];distat=sqrt(Ca(:,1).^2+Ca(:,2).^2+Ca(:,3).^2);
     distatro = sqrt((Ca(:,1)-xw).^2+(Ca(:,2)-yw).^2+(Ca(:,3)-zw).^2);
 else direc == -1
     [Ca] = [xp yp zp]-[repmat(norms.*vf*porc/100,[ 1 3]).*[u v w]];distat=sqrt(Ca(:,1).^2+Ca(:,2).^2+Ca(:,3).^2);
     distatro = sqrt((xp-Ca(:,1)).^2+(yp-Ca(:,2)).^2+(zp-Ca(:,3)).^2);
 end
 Surft.SurfData.vertices = Ca;
 Surft.SurfData.faces = Surfptemp.SurfData.faces;
 Surfa{3,1} = Surft;
 N = 50; %dista = distatro/Nt;
% xw = Surfwtemp.SurfData.vertices(:,1); yw = Surfwtemp.SurfData.vertices(:,2); zw = Surfwtemp.SurfData.vertices(:,3);
% xp = Surfptemp.SurfData.vertices(:,1); yp = Surfptemp.SurfData.vertices(:,2); zp = Surfptemp.SurfData.vertices(:,3);
% norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);
% 
%%%%%%%%%%%%%%%%%%%%%%

% Cp = [xp yp zp]+[repmat(scalef*distatro,[ 1 3]).*[u v w]];
% Cw = [xw yw zw]-[repmat(scalef*distatro,[ 1 3]).*[u v w]];
% dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));
% steps = [0:1/N:1];dista = repmat(steps,[size(dista,1) 1]).*repmat(dista,[1 N+1]);
% CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(u,[1 size(dista,2)]).*dista ]';siz = size(CintX);%CintX = CintX(:);
% CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(v,[1 size(dista,2)]).*dista ]';%CintY = CintY(:);
% CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(w,[1 size(dista,2)]).*dista ]';%CintZ = CintZ(:);

%%%%%Quitar esta parte si no funciona
%------------------------------------
Cp = [xp yp zp];
Cw = [xw yw zw];
if direc == 1
    Cw = Ca;
elseif direc == -1
    Cp = Ca;
end

N=5;
%dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));
steps = [0:1/N:1];dista = repmat(steps,[size(distatro,1) 1]).*repmat(distatro,[1 N+1]);
CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(u,[1 size(dista,2)]).*dista ]';siz = size(CintX);%CintX = CintX(:);
CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(v,[1 size(dista,2)]).*dista ]';%CintY = CintY(:);
CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(w,[1 size(dista,2)]).*dista ]';%CintZ = CintZ(:);
%---------------------------------

%%%%%%%%% Leyendo el centro de coordenadas RAS
if isunix
    cras = textread([dira 'ch2/mri/transforms/talairach.lta'],'%s',5,'headerlines',20);
else
    cras = textread([dira 'ch2\mri\transforms\talairach.lta'],'%s',5,'headerlines',20);
end
cras = char(cras);cras = str2num(cras(3:end,:))';

V = spm_vol([dira 'ch2.nii']); It =spm_read_vols(V);       I =0*It;
V1 = spm_vol([dira 'p1ch2.nii']); It1 =spm_read_vols(V1); 
V2 = spm_vol([dira 'p2ch2.nii']); It2 =spm_read_vols(V2); 
V3 = spm_vol([dira 'p3ch2.nii']); It3 =spm_read_vols(V3); 
[pth, nm,ext] = fileparts(V.fname);
% 

voxs = (inv(V.mat)*[CintX(:)+cras(1) CintY(:)+cras(2) CintZ(:)+cras(3) ones(size(CintX(:),1),1)]')';voxX = reshape(voxs(:,1),[siz]);voxY = reshape(voxs(:,2),[siz]);voxZ = reshape(voxs(:,3),[siz]);


%% Probando cosas nuevas

% [X,Y,Z] = meshgrid([floor(min(voxX(:))):ceil(max(voxX(:)))],[floor(min(voxY(:))):ceil(max(voxY(:)))],[floor(min(voxZ(:))):ceil(max(voxZ(:)))]);X =X(:);Y =Y(:);Z =Z(:);
% indt = sub2ind(size(It),X,Y,Z);
% close all;
% figure;  colordef black;
% Is = zeros(size(voxY));Is1 = zeros(size(voxY));Is2 = zeros(size(voxY));Is3 = zeros(size(voxY));
%subplot(1,2,2);
%strsurf=patch(Surft.SurfData,'edgecolor','black','facecolor',[1 0 0],'tag', 'patch','facelighting','gouraud');camlight;axis equal;
% for j = 1:size(voxX,2)
%     a = floor([voxX(:,j) voxY(:,j) voxZ(:,j)]);b = ceil([voxX(:,j) voxY(:,j) voxZ(:,j)]);c = unique([a;b],'rows');indt = sub2ind(size(It),c(:,1),c(:,2),c(:,3));
%     indt = sub2ind(size(It),c(:,1),c(:,2),c(:,3));Vn = It(indt);Vn1 = It1(indt);Vn2 = It2(indt);Vn3 = It3(indt);
%     dMt = sqrt((X-voxs(j,1)).^2+(Y-voxs(j,2)).^2+(Z-voxs(j,3)).^2);
%     dMt = dist([[voxX(:,j) voxY(:,j) voxZ(:,j)]; c]'); 
%     Mi = dMt(1:size(voxX,1),size(voxX,1)+1:end);
%     Mi = logical(Mi)./Mi;
%     Mi = Mi./repmat(sum(Mi')',[1 size(c,1)]);
%     Is(:,j) = (Mi*Vn)';
%     Is1(:,j) = (Mi*Vn1)';
%     Is2(:,j) = (Mi*Vn2)';
%     Is3(:,j) = (Mi*Vn3)';
%     XAju =sqrt((CintX(:,j)-CintX(1,j)).^2+(CintY(:,j)-CintY(1,j)).^2+(CintZ(:,j)-CintZ(1,j)).^2)-scalef*distatro(j);
%     indcth = find((XAju>=0)&(XAju<=norms(j)));
%     if direc == 1
%         indatro = find(((XAju)>=0)&(XAju<=distatro(j)));
%     elseif direc == -1
%         indatro = find((XAju>=(norms(j)-distatro(j)))&(XAju<=norms(j)));
%     end
%     YAju = Is(:,j);
%     YAju1 = Is1(:,j);YAju2 = Is2(:,j);YAju3 = Is3(:,j);
%    % I(round(voxX(:,j)),round(voxY(:,j)),round(voxZ(:,j))) = 1;
%     %I(round(voxX(indcth,j)),round(voxY(indcth,j)),round(voxZ(indcth,j))) = 2;
%    % I(round(voxX(indatro,j)),round(voxY(indatro,j)),round(voxZ(indatro,j))) = 3;
%     subplot(2,2,1); h = gca;cla(h,'reset');
%     plot(XAju,YAju,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju(indcth),'or');
%     plot(XAju(indatro),YAju(indatro),'oy');
%     drawnow;
%     subplot(2,2,2);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju1,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju1(indcth),'or');
%     plot(XAju(indatro),YAju1(indatro),'oy');
%     drawnow;
%     subplot(2,2,3);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju2,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju2(indcth),'or');
%     plot(XAju(indatro),YAju2(indatro),'oy');
%     drawnow;
%     subplot(2,2,4);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju3,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju3(indcth),'or');
%     plot(XAju(indatro),YAju3(indatro),'oy');
%     drawnow;
%     %hold on;
%     %plot3(Surft.SurfData.vertices(1:j,1),Surft.SurfData.vertices(1:j,2),Surft.SurfData.vertices(1:j,3),'.g','Markersize',10);
%     
% end

%%%%%%%%%%%%%%


%% Probando cosas REnuevas

%[X,Y,Z] = meshgrid([floor(min(voxX(:))):ceil(max(voxX(:)))],[floor(min(voxY(:))):ceil(max(voxY(:)))],[floor(min(voxZ(:))):ceil(max(voxZ(:)))]);X =X(:);Y =Y(:);Z =Z(:);
a =unique(round([voxX(:) voxY(:) voxZ(:)]),'rows');
indt = sub2ind(size(It),a(:,1),a(:,2),a(:,3));
IsO = spm_sample_vol(It,voxX(:),voxY(:),voxZ(:),7);
%close all;
%figure;  colordef black;
%Is = zeros(size(voxY));Is1 = zeros(size(voxY));Is2 = zeros(size(voxY));Is3 = zeros(size(voxY));
%subplot(1,2,2);
%strsurf=patch(Surft.SurfData,'edgecolor','black','facecolor',[1 0 0],'tag', 'patch','facelighting','gouraud');camlight;axis equal;
Mat = zeros(size(indt,1),size(voxX(:),1));
th = 2; %Threshold para la vecindad
for j = 1:size(indt,1)
    %dMt = dist([[a(j) a(j) a(j)]; voxX(:) voxY(:) voxZ(:)]'); 
    dMt = sqrt((voxX(:)-a(j,1)).^2+(voxY(:)-a(j,2)).^2+(voxZ(:)-a(j,3)).^2);
    ind = find(dMt<=th);
    Mat(j,ind) = (dMt(ind))./sum(dMt(ind));
    
    % a1 = floor([voxX(:,j) voxY(:,j) voxZ(:,j)]);b1 = ceil([voxX(:,j) voxY(:,j) voxZ(:,j)]);c1 = unique([a1;b1],'rows');indt = sub2ind(size(It),c1(:,1),c1(:,2),c1(:,3));
    % indt = sub2ind(size(It),c1(:,1),c1(:,2),c1(:,3));Vn = It(indt);%Vn1 = It1(indt);Vn2 = It2(indt);Vn3 = It3(indt);
     
     %dMt1 = dist([[voxX(:,j) voxY(:,j) voxZ(:,j)]; c1]'); 
    % Mi = dMt1(1:size(voxX,1),size(voxX,1)+1:end);
    % Mi = logical(Mi)./Mi;
   %  Mi = Mi./repmat(sum(Mi')',[1 size(c1,1)]);
   %  Is(:,j) = (Mi*Vn)';

end



[X,Y,Z] = meshgrid([floor(min(voxX(:))):ceil(max(voxX(:)))],[floor(min(voxY(:))):ceil(max(voxY(:)))],[floor(min(voxZ(:))):ceil(max(voxZ(:)))]);X =X(:);Y =Y(:);Z =Z(:);
% indt = sub2ind(size(It),X,Y,Z);
% close all;
% figure;  colordef black;
Is = zeros(size(voxY));Is1 = zeros(size(voxY));Is2 = zeros(size(voxY));Is3 = zeros(size(voxY));
%subplot(1,2,2);
%strsurf=patch(Surft.SurfData,'edgecolor','black','facecolor',[1 0 0],'tag', 'patch','facelighting','gouraud');camlight;axis equal;
 for j = 1:size(voxX,1)
     a = floor([voxX(j) voxY(j) voxZ(j)]);b = ceil([voxX(:,j) voxY(:,j) voxZ(:,j)]);c = unique([a;b],'rows');indt = sub2ind(size(It),c(:,1),c(:,2),c(:,3));

     indt = sub2ind(size(It),c(:,1),c(:,2),c(:,3));Vn = It(indt);Vn1 = It1(indt);Vn2 = It2(indt);Vn3 = It3(indt);
     dMt = sqrt((X-voxX(:,j)).^2+(Y-voxY(:,j)).^2+(Z-voxZ(:,j)).^2);
     ind = find(dMt<=th);
%     dMt = dist([[voxX(:,j) voxY(:,j) voxZ(:,j)]; c]'); 
%     Mi = dMt(1:size(voxX,1),size(voxX,1)+1:end);
     Mi = logical(dMt(ind))./dMt(ind);
     Mi = Mi./repmat(sum(Mi')',[1 size(c,1)]);
%     Is(:,j) = (Mi*Vn)';
%     Is1(:,j) = (Mi*Vn1)';
%     Is2(:,j) = (Mi*Vn2)';
%     Is3(:,j) = (Mi*Vn3)';
%     XAju =sqrt((CintX(:,j)-CintX(1,j)).^2+(CintY(:,j)-CintY(1,j)).^2+(CintZ(:,j)-CintZ(1,j)).^2)-scalef*distatro(j);
%     indcth = find((XAju>=0)&(XAju<=norms(j)));
%     if direc == 1
%         indatro = find(((XAju)>=0)&(XAju<=distatro(j)));
%     elseif direc == -1
%         indatro = find((XAju>=(norms(j)-distatro(j)))&(XAju<=norms(j)));
%     end
%     YAju = Is(:,j);
%     YAju1 = Is1(:,j);YAju2 = Is2(:,j);YAju3 = Is3(:,j);
%    % I(round(voxX(:,j)),round(voxY(:,j)),round(voxZ(:,j))) = 1;
%     %I(round(voxX(indcth,j)),round(voxY(indcth,j)),round(voxZ(indcth,j))) = 2;
%    % I(round(voxX(indatro,j)),round(voxY(indatro,j)),round(voxZ(indatro,j))) = 3;
%     subplot(2,2,1); h = gca;cla(h,'reset');
%     plot(XAju,YAju,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju(indcth),'or');
%     plot(XAju(indatro),YAju(indatro),'oy');
%     drawnow;
%     subplot(2,2,2);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju1,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju1(indcth),'or');
%     plot(XAju(indatro),YAju1(indatro),'oy');
%     drawnow;
%     subplot(2,2,3);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju2,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju2(indcth),'or');
%     plot(XAju(indatro),YAju2(indatro),'oy');
%     drawnow;
%     subplot(2,2,4);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju3,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju3(indcth),'or');
%     plot(XAju(indatro),YAju3(indatro),'oy');
%     drawnow;
%     %hold on;
%     %plot3(Surft.SurfData.vertices(1:j,1),Surft.SurfData.vertices(1:j,2),Surft.SurfData.vertices(1:j,3),'.g','Markersize',10);
%     
 end








%It(indt) = Mat*
    
    
%     f = (Mat'*Mat)

%     Is1(:,j) = (Mi*Vn1)';
%     Is2(:,j) = (Mi*Vn2)';
%     Is3(:,j) = (Mi*Vn3)';
%     XAju =sqrt((CintX(:,j)-CintX(1,j)).^2+(CintY(:,j)-CintY(1,j)).^2+(CintZ(:,j)-CintZ(1,j)).^2)-scalef*distatro(j);
%     indcth = find((XAju>=0)&(XAju<=norms(j)));
%     if direc == 1
%         indatro = find(((XAju)>=0)&(XAju<=distatro(j)));
%     elseif direc == -1
%         indatro = find((XAju>=(norms(j)-distatro(j)))&(XAju<=norms(j)));
%     end
%     YAju = Is(:,j);
%     YAju1 = Is1(:,j);YAju2 = Is2(:,j);YAju3 = Is3(:,j);
%    % I(round(voxX(:,j)),round(voxY(:,j)),round(voxZ(:,j))) = 1;
%     %I(round(voxX(indcth,j)),round(voxY(indcth,j)),round(voxZ(indcth,j))) = 2;
%    % I(round(voxX(indatro,j)),round(voxY(indatro,j)),round(voxZ(indatro,j))) = 3;
%     subplot(2,2,1); h = gca;cla(h,'reset');
%     plot(XAju,YAju,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju(indcth),'or');
%     plot(XAju(indatro),YAju(indatro),'oy');
%     drawnow;
%     subplot(2,2,2);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju1,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju1(indcth),'or');
%     plot(XAju(indatro),YAju1(indatro),'oy');
%     drawnow;
%     subplot(2,2,3);
%     h = gca;cla(h,'reset');
%     plot(XAju,YAju2,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju2(indcth),'or');
%     plot(XAju(indatro),YAju2(indatro),'oy');
%     drawnow;
%     subplot(2,2,4);
% %     h = gca;cla(h,'reset');
%     plot(XAju,YAju3,'bo');%drawnow;
%     hold on;
%     plot(XAju(indcth),YAju3(indcth),'or');
%     plot(XAju(indatro),YAju3(indatro),'oy');
%     drawnow;
    %hold on;
    %plot3(Surft.SurfData.vertices(1:j,1),Surft.SurfData.vertices(1:j,2),Surft.SurfData.vertices(1:j,3),'.g','Markersize',10);
    
%end



%%


Is = reshape(spm_sample_vol(It,voxs(:,1),voxs(:,2),voxs(:,3),7),[siz]);
Isl = reshape(spm_sample_vol(It,voxs(:,1),voxs(:,2),voxs(:,3),7),[siz]);
%f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
close all;
figure; 
for j = 1:size(Is,2)
    j
    XAju =sqrt((CintX(:,j)-CintX(1,j)).^2+(CintY(:,j)-CintY(1,j)).^2+(CintZ(:,j)-CintZ(1,j)).^2)-scalef*distatro(j);
    indcth = find((XAju>=0)&(XAju<=norms(j)));
    if direc == 1
        indatro = find(((XAju)>=0)&(XAju<=distatro(j)));
    elseif direc == -1
        indatro = find((XAju>=(norms(j)-distatro(j)))&(XAju<=norms(j)));
    end
    YAju = Is(:,j);
    YAjul = Isl(:,j);
    %XAju = [1:length(YAju)]';
    %indmin = find(Is(:,j)==min(Is(:,j)));indmax = find(Is(:,j)==max(Is(:,j)));

    %p = nlinfit(XAju(indmin-1:indmax+1),YAju(indmin-1:indmax+1),f,[0 20 50 5]);
    %p = nlinfit(XAju,YAju,f,[0 20 50 5]);
    %plot(XAju(indmin-1:indmax+1),YAju(indmin-1:indmax+1),'bo');
    %line(XAju(indmin-1:indmax+1),f(p,XAju(indmin-1:indmax+1)),'color','r');
    h = gca;cla(h,'reset');
   % I(round(voxX(:,j)),round(voxY(:,j)),round(voxZ(:,j))) = 1;
    %I(round(voxX(indcth,j)),round(voxY(indcth,j)),round(voxZ(indcth,j))) = 2;
   % I(round(voxX(indatro,j)),round(voxY(indatro,j)),round(voxZ(indatro,j))) = 3;
    
    plot(XAju,YAju,'bo');%drawnow;
    hold on;
    plot(XAju(indcth),YAju(indcth),'or');
    plot(XAju(indatro),YAju(indatro),'oy');
    drawnow;
    
    %
    %line(XAju,f(p,XAju),'color','r');
end
V.fname = [pth filesep nm '_pitopitocolorito.nii']; 
    write_atlas_vol(V,I);

%%%%%%%%%%%Reduccion de corteza



%% 


% 
% %%%%%%%%%%%%%%%%%%%%%%
% xw = Surfwtemp.SurfData.vertices(:,1); yw = Surfwtemp.SurfData.vertices(:,2); zw = Surfwtemp.SurfData.vertices(:,3);
% xp = Surfptemp.SurfData.vertices(:,1); yp = Surfptemp.SurfData.vertices(:,2); zp = Surfptemp.SurfData.vertices(:,3);
% norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);
% u = (xp-xw)./norms; v = (yp-yw)./norms; w = (zp-zw)./norms;
% %%%%%%%%%%%Reduccion de corteza
% Surfo = Surfa;
% if direc == 1
%     Surfptemp.SurfData.vertices = Surfwtemp.SurfData.vertices+[repmat(norms.*v*porc/100,[ 1 3]).*[u v w]];Surfa{3,1} = Surfptemp;
% else direc == -1
%     Surfwtemp.SurfData.vertices = Surfptemp.SurfData.vertices-[repmat(norms.*v*porc/100,[ 1 3]).*[u v w]];Surfa{3,1} = Surfwtemp;
% end
% %%%%%%%%%%%%%%%%
% Plot_Surf(Surfa);
% 
% xw = Surfwtemp.SurfData.vertices(:,1); yw = Surfwtemp.SurfData.vertices(:,2); zw = Surfwtemp.SurfData.vertices(:,3);
% xp = Surfptemp.SurfData.vertices(:,1); yp = Surfptemp.SurfData.vertices(:,2); zp = Surfptemp.SurfData.vertices(:,3);
% norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);
% u = (xp-xw)./norms; v = (yp-yw)./norms; w = (zp-zw)./norms;
% 
% 
% steps = [0:1/N:1];steps([1 end]) = [];dista = repmat(steps,[size(norms,1) 1]).*repmat(norms,[1 N-1]);
% %     Nsurf = zeros(size(xw,1),length(steps)*3);
% %     Nsurf(:,1:3:end) = repmat(xw,[1 size(dista,2)])+repmat(u,[1 size(dista,2)]).*dista;
% %     Nsurf(:,2:3:end) = repmat(yw,[1 size(dista,2)])+repmat(v,[1 size(dista,2)]).*dista;
% %     Nsurf(:,3:3:end) = repmat(zw,[1 size(dista,2)])+repmat(w,[1 size(dista,2)]).*dista;
% Surft = Surfwtemp;
% for i = 1:N-1
%     Surft.SurfData.faces = [Surft.SurfData.faces; Surfwtemp.SurfData.faces+ size(Surft.SurfData.vertices,1)];
%     Surft.SurfData.vertices= [ Surft.SurfData.vertices;[xw+u.*dista(:,i) yw+v.*dista(:,i) zw+w.*dista(:,i)]];
% 
% end
% Surft.SurfData.faces = [Surft.SurfData.faces;;Surfptemp.SurfData.faces+size(Surft.SurfData.vertices,1)];
% Surft.SurfData.vertices= [ Surft.SurfData.vertices;Surfptemp.SurfData.vertices];
% Surfo{3,1} = Surft;Plot_Surf(Surfo);
% %     if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
% %     delx = diff([min(x(:)) max(x(:))])/n;
% %     dely = diff([min(y(:)) max(y(:))])/m;
% %     delz = diff([min(z(:)) max(y(:))])/max(m,n);
% %     del = sqrt(delx.^2 + dely.^2 + delz.^2);
% %     if del>0
% %         len = sqrt((u/del).^2 + (v/del).^2 + (w/del).^2);
% %         maxlen = max(len(:));
% %     else
% %         maxlen = 0;
% %     end
% %
% %     if maxlen>0
% %         scf = scf*0.9 / maxlen;
% %     else
% %         scf = scf*0.9;
% %     end
% %plot(norms-thickness(indf),'-b');
% if isunix
% cras = textread([dira 'ch2/mri/transforms/talairach.lta'],'%s',5,'headerlines',20);
% else
%     cras = textread([dira 'ch2\mri\transforms\talairach.lta'],'%s',5,'headerlines',20);
% end
% cras = char(cras);cras = str2num(cras(3:end,:))';
% 
% Surft.SurfData.vertices = Surft.SurfData.vertices + repmat(cras,[size(Surft.SurfData.vertices,1) 1]);
% % [Tri] = Vert_Neib(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
% % Temp = sum(Tri);
% % Tri(:,Temp==0) = [];
% % Surft.Tri = Tri;
% 
% V = spm_vol([dira 'ch2.nii']); It =spm_read_vols(V);I = 0*It;
% voxs = round(inv(V.mat)*[Surft.SurfData.vertices ones(size(Surft.SurfData.vertices,1),1)]')';
% ind = sub2ind(size(I),voxs(:,1), voxs(:,2), voxs(:,3));I(ind) = 1;I = imfill(I,'holes');
% [pth, nm,ext] = fileparts(V.fname);V.fname = [pth nm '_test.nii'];
% write_atlas_vol(V,I);
% 
% %[OutF, Surfa] = Exp_Surf(deblank(OutputFiles(i,:)), V(i).fname, '','frb', 'exp','y');
% freesurfer_write_surf([dira 'chd_lh.pial.frb'], Surft.SurfData.vertices, Surft.SurfData.faces);
% 
% 
% % col = abs([u v w]);
% % col = col./repmat(sqrt(sum((col').^2))',[1 3]);
% % u = u*scf; v = v*scf; w = w*scf;
% % beta = beta * sqrt(u.*u + v.*v + w.*w) ./ (sqrt(u.*u + v.*v) + eps);
% % hu = [xp-alpha*(u+beta.*(v+eps)) xp xp-alpha*(u-beta.*(v+eps))];
% % hv = [yp-alpha*(v-beta.*(u+eps)) yp yp-alpha*(v+beta.*(u+eps))];
% % hw = [zp-alpha*w zp zp-alpha*w];
% % for i = 1:size(u(:),1)
% %     line([xw(i) xp(i)], [yw(i) yp(i)],[zw(i) zp(i)],'Color',col(i,:),'LineWidth',1);
% %     line([hu(i,1) hu(i,2) hu(i,3)], [hv(i,1) hv(i,2) hv(i,3)],[hw(i,1) hw(i,2) hw(i,3)],'Color',col(i,:),'LineWidth',1);
% % end
% % inda
% % for i = 1:size(inda,1)
% %     line([xw(inda(i)) xp(inda(i))], [yw(inda(i)) yp(inda(i))],[zw(inda(i)) zp(inda(i))],'Color',col(inda(i),:),'LineWidth',1);
% %     line([hu(inda(i),1) hu(inda(i),2) hu(inda(i),3)], [hv(inda(i),1) hv(inda(i),2) hv(inda(i),3)],[hw(inda(i),1) hw(inda(i),2) hw(inda(i),3)],'Color',col(inda(i),:),'LineWidth',1);
% % end
% % return;


%%%%%%%%%%%%%
return

function v = sigmoid(params,range)

% Sigmoid creates a Sigmoid function using parameters in PARAMS and the 
% variable range.
% 
% V = SIGMOID(PARAMS,RANGE)
%
% PARAMS: a 3-vector, the entries of which are (in this order):
% amplitude value 
% 
% phase
% slope

amplitude = params(1);
Phase=params(2);
Slope=params(3);
%a=params(4);
v=1./(1+double(Phase)*exp(-double(Slope)*(range)));


function Vertnlevels = Neib_level(Stest,Nlevels);
if nargin ==1
    Nlevels =1;
end
col = [0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1); 
re = floor(Nlevels/Ncolor); col = repmat(col,[re+1 1]);
[Tri] = Vert_Neib(double(Stest.SurfData.faces),size(Stest.SurfData.vertices,1),size(Stest.SurfData.faces,1));
Temp = sum(Tri);
Tri(:,Temp==0) = [];
Stest.Tri =Tri;

[fa,indtt] =sort(Stest.SurfData.faces');
fa = fa';
%[ford,indt2] =sortrows(fa');ford = double(ford);
%temp = ford(1:end-1,:) - ford(2:end,:);temp = sum(abs(temp'));
%fnd = find(temp == 0);
%Stest.SurfData.faces(indt2(fnd),:)=[];
Stest.SurfData.faces = uint32(Stest.SurfData.faces);
bound =[0];
% Plot_Surf(Stest);
for i = 1:size(Stest.SurfData.vertices,1)
    if (sum(ismember(bound,i)) == 0)|i~=size(Stest.SurfData.vertices,1)
        Neib = Stest.Tri(i,3:end); Neib = unique(nonzeros(Neib(:)));A = unique(Stest.SurfData.faces(Neib,:));A(A==i) = [];
        for j = 1:size(A,1)
            if sum(ismember(bound,A(j))) == 0 
                X = sort([i A(j)]);
                a = [ismember(fa(:,1:2),X,'rows') ismember(fa(:,2:3),X,'rows') ismember(fa(:,[1 3]),X,'rows')];
                if sum(a(:)) ==1;
                    bound = [bound [i A(j)]];
                    
                else
                end
            end
        end
    end
end
bound(bound ==0) =[]; bound = unique(bound);bound = bound(:);
Vertnlevels = [bound  bound*0+1];
Vertnlevelst = Vertnlevels;

% hold on;plot3(Stest.SurfData.vertices(Vertnlevels(:,1),1),Stest.SurfData.vertices(Vertnlevels(:,1),2),Stest.SurfData.vertices(Vertnlevels(:,1),3),'.','Markersize',10,'Color',[0 1 0])
 
% if Nlevels > 1
    %for i = 2:Nlevels
    i =1;
    while sum(ismember([1:size(Stest.SurfData.vertices,1)]',Vertnlevels(:,1)))~=size(Stest.SurfData.vertices,1)
        i= i+1;
        boundt = 0;
        for j = 1:size(Vertnlevelst,1)
            Neib = Stest.Tri(Vertnlevelst(j),3:end); Neib = unique(nonzeros(Neib(:)));A = unique(Stest.SurfData.faces(Neib,:));A(A==Vertnlevelst(j)) = [];
            boundt =[boundt(:);A];
        end
        boundt(boundt ==0) =[]; boundt = unique(boundt);
        indd = ismember(boundt,Vertnlevels(:,1));boundt(indd) = [];
        
        %hold on;plot3(Stest.SurfData.vertices(boundt,1),Stest.SurfData.vertices(boundt,2),Stest.SurfData.vertices(boundt,3),'.','Markersize',10,'Color',col(i,:))
        
        Vertnlevelst = [boundt boundt*0+i]; 
        Vertnlevels = [Vertnlevels; Vertnlevelst];
    end
% end
% for i = 1:Nlevels
%     ind = find(Vertnlevels(:,2) == i);
%     hold on;plot3(Stest.SurfData.vertices(Vertnlevels(ind,1),1),Stest.SurfData.vertices(Vertnlevels(ind,1),2),Stest.SurfData.vertices(Vertnlevels(ind,1),3),'.','Markersize',10,'Color',col(i,:))
% end

return;



function freesurfer_fwrite3(fid, val)
% freesurfer_fwrite3 - FreeSurfer function to write a 3 byte integer to a file
%
% freesurfer_fwrite3(fid, val)
%
% see also freesurfer_read3, freesurfer_read_surf, freesurfer_write_surf
 if(nargin ~= 2)
  fprintf('USAGE: freesurfer_fwrite3(fid, val)\n');
  return;
end

%fwrite(fid, val, '3*uchar') ;
b1 = bitand(bitshift(val, -16), 255) ;
b2 = bitand(bitshift(val, -8), 255) ;
b3 = bitand(val, 255) ; 
fwrite(fid, b1, 'uchar') ;
fwrite(fid, b2, 'uchar') ;
fwrite(fid, b3, 'uchar') ;
 
return

function freesurfer_write_surf(fname, vert, face)

% freesurfer_write_surf - FreeSurfer I/O function to write a surface file
%
% freesurfer_write_surf(fname, vert, face)
%
% writes a surface triangulation into a binary file
% fname - name of file to write
% vert  - Nx3 matrix of vertex coordinates
% face  - Mx3 matrix of face triangulation indices
%
% The face matrix here must be matlab compatible
% (no zero indices).  It is converted to FreeSurfer
% indices that start at zero.
%
% See also freesurfer_read_surf, freesurfer_write_curv, freesurfer_write_wfile

if(nargin ~= 3)
    fprintf('USAGE: freesurfer_write_surf(fname, vert, face)\n');
    return;
end

if size(vert,2) ~= 3,
    error('vert must be Nx3 matrix');
end

if size(face,2) ~= 3,
    error('face must be Mx3 matrix');
end

fprintf('...subtracting 1 from face indices for FreeSurfer compatibility.\n');
face = face - 1;

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;

TRIANGLE_FILE_MAGIC_NUMBER = 16777214 ;
freesurfer_fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);

vnum = size(vert,1) ;  % number of vertices
fnum = size(face,1) ;  % number of faces

% Ouput a couple of text lines with creation date
str = sprintf('created from matalb on %s\n',datestr(now));
fwrite(fid, str,'char');
fwrite(fid, str,'char');

fwrite(fid, vnum,'int32');
fwrite(fid, fnum,'int32');

% reshape vert into column array and write
vert = reshape(vert',size(vert,1)*size(vert,2),1);
fwrite(fid, vert,'float32');

% reshape face into column array and write
face = reshape(face',size(face,1)*size(face,2),1);
fwrite(fid, face,'int32');

fclose(fid) ;

return


function [vertices, faces] = freesurfer_read_surf(fname)

% freesurfer_read_surf - FreeSurfer I/O function to read a surface file
%
% [vertices, faces] = freesurfer_read_surf(fname)
%
% Reads the vertex coordinates (mm) and face lists from a surface file.
%
% Surface files are stored as either triangulations or quadrangulations.
% That is, for a triangulation, each face is defined by 3 vertices.  For a
% quadrangulation, each face is defined by 4 vertices.  The rows of 'faces'
% contain indices into the rows of 'vertices', the latter holds the XYZ
% coordinates of each vertex.
%
% The freesurfer faces index the vertices in counter-clockwise order (when
% viewed from the outside of the surface).  This is consistent with a
% right-hand rule.  If we have vertices
%
% C           B
%
%
%       A
%
% Then we can calculate an edge vector from A to B (ie, AB = B - A) and
% another edge vector from A to C (ie, AC = C - A).  If you form a "gun"
% with your thumb and forefinger of the right hand, then align your thumb
% with the AB vector and your forefinger with the AC vector, your palm is
% facing out of the screen and extending your middle finger in the
% orthogonal direction to the plane of the screen will give the outward
% surface normal of the triangle ABC.  (If you lookup "triangle" on
% Wolfram's mathworld, you can see that AB is referred to as c and AC is
% referred to as b.)
%
% However, if this surface is read into matlab, it will give INWARD surface
% normals in the matlab patch command.  For some reason, matlab is not
% following the right hand rule.  To get OUTWARD normals with the matlab
% patch command, use faces(:,[1 3 2]) (see below).
%
% The vertex coordinates are in mm.  The FreeSurfer coordinate
% system for surfaces is quite simple, but relating to their MRI
% cor-??? files is too confusing to explain here; see the FreeSurfer
% homepage or google the documentation by Graham Wideman.  For the
% surfaces, at least, the origin is somewhere in the center of the
% head, and the vertex XYZ coordinates are oriented such that +X is
% right, +Y is anterior and +Z is superior (this is the
% FreeSurfer RAS coordinate system).
%
% Note that reading the faces of a quad file can take a long
% time due to their compact storage format.  In this case, the return of
% vertices can be faster if the face output variable is not specified; in
% this case, the faces are not read.
%
% Try this to visualize the surface:
% Hp = patch('vertices',vertices,'faces',faces(:,[1 3 2]),...
%       'facecolor',[.5 .5 .5],'edgecolor','none')
% camlight('headlight','infinite')
% vertnormals = get(Hp,'vertexnormals');
%
% See also freesurfer_write_surf, freesurfer_read_curv,
%          freesurfer_read_wfile
%

% $Revision: 1.1.2.1 $ $Date: 2007/12/09 22:50:25 $

% Copyright (C) 2000  Darren L. Weber
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  08/2000, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1.2.1 $ $Date: 2007/12/09 22:50:25 $';
fprintf('FREESURFER_READ_SURF [v %s]\n',ver(11:15));

if(nargin < 1)
    help freesurfer_read_surf;
    return;
end

%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER  =  16777214;
QUAD_FILE_MAGIC_NUMBER      =  16777215;


% open it as a big-endian file
fid = fopen(fname, 'rb', 'b');
if (fid < 0),
    str = sprintf('could not open surface file %s.', fname);
    error(str);
end

fprintf('...reading surface file: %s\n', fname);
tic;

magic = freesurfer_fread3(fid);

if (magic == QUAD_FILE_MAGIC_NUMBER),
    Nvertices = freesurfer_fread3(fid);
    Nfaces = freesurfer_fread3(fid);
    fprintf('...reading %d quad file vertices\n',Nvertices);
    vertices = fread(fid, Nvertices*3, 'int16') ./ 100 ;
    if (nargout > 1),
        fprintf('...reading %d quad file faces (please wait)\n',Nfaces);
        faces = zeros(Nfaces,4);
        for iface = 1:Nfaces,
            for n=1:4,
                faces(iface,n) = freesurfer_fread3(fid) ;
            end
            if(~rem(iface, 10000)), fprintf(' %7.0f',iface); end
            if(~rem(iface,100000)), fprintf('\n'); end
        end
    end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER),
    fprintf('...reading triangle file\n');
    tline = fgets(fid); % read creation date text line
    tline = fgets(fid); % read info text line

    Nvertices = fread(fid, 1, 'int32'); % number of vertices
    Nfaces = fread(fid, 1, 'int32'); % number of faces

    % vertices are read in column format and reshaped below
    vertices = fread(fid, Nvertices*3, 'float32');

    % faces are read in column format and reshaped
    faces = fread(fid, Nfaces*3, 'int32');
    faces = reshape(faces, 3, Nfaces)';
else
    str = sprintf('unknown magic number in surface file %s.', fname);
    error(str);
end

vertices = reshape(vertices, 3, Nvertices)';
fclose(fid);

fprintf('...adding 1 to face indices for matlab compatibility.\n');
faces = faces + 1;

t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return
















