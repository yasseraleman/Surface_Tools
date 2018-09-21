function  [Dist,Int,Icvar] = Surface_sphere_interpolation(Surf,cvar,Surft,flag,itype)
%Interpolacion Lineal sobre una esfera y determinacion de intersecciones
%entre una superficie y una linea

%Surf = read_meta_surf(deblank(Surfa(1,:)));
if nargin<4
    flag = 'indsurf';
else
    flag = 'sphere';
end
%% Reading Surfaces
wh = whos('Surf');
if(strcmp(wh.class,'cell'))
    Surftemp = Surf{1}; clear Surf;Surf = Surftemp;
elseif ischar(Surf);
    [OutFiles, SurfF] = Exp_Surf(Surf, '0', '','', 'imp','n');Surf = SurfF{1};
end

wh = whos('Surft');
if(strcmp(wh.class,'cell'))
    Surftemp = Surft{1}; clear Surft;Surft = Surftemp;
elseif ischar(Surft);
    [OutFiles, SurfF] = Exp_Surf(Surft, '0', '','', 'imp','n');Surft = SurfF{1};
end

%% Reading Characteristic
wh = whos('cvar');
if ischar(cvar);
    [cvar,ctab] = read_cfiles(cvar);
end

%% Line Part
%Surft = read_meta_surf(deblank(Surfa(2,:)));
Vert = Surf.SurfData.vertices;Faces = Surf.SurfData.faces;
n = cross(Vert(Faces(:,2),:)-Vert(Faces(:,1),:),Vert(Faces(:,3),:)-Vert(Faces(:,1),:),2); D = -1*dot(n,Vert(Faces(:,1),:),2);
%[Surft] = Reorg_Surf(Surft);
fv.vertices = Surft.SurfData.vertices;
fv.faces = Surft.SurfData.faces;
if ~isfield(Surft.SurfData,'VertexNormals')
h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
else
    normals = Surft.SurfData.VertexNormals;
end


 norma = sqrt(sum((normals').^2));
 normals = normals./repmat(norma',[1 3]);
% temp = sum(normals')';
% ind = find(isnan(temp) ==1);
% if ~isempty(ind)
%     indt = ind;
%     for i = 1:size(ind,1)
%         vert = Surf.SurfData.faces(nonzeros(Surf.Tri(ind(i),3:end)));ind2 = ismember(vert,indt);vert(ind2) = [];
%         tempN = sum([normals(vert,1) normals(vert,2) normals(vert,3)])/size(vert,1);
%         normals(ind(i),:) = tempN;
%         norma(ind(i)) = sqrt(sum((tempN').^2));
%         indt(indt == ind(i)) = [];
%     end
% end
% Surf.SurfData.VertexNormals = normals./repmat(norma',[1 3]);
% normals = Surf.SurfData.VertexNormals;




%P1 = Surft.SurfData.vertices+100*normals; 
P1 = Surft.SurfData.vertices*0; 
P2 = Surft.SurfData.vertices+100*normals;
%P2 = Surft.SurfData.vertices-100*normals;
Ns = size(Surft.SurfData.vertices,1);
 H = waitbar(0,['Computing Distance '],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');
  if ~isempty(strcmp(flag,'sphere'))
        [Th,Phi,R] = cart2sph(Surf.SurfData.vertices(:,1),Surf.SurfData.vertices(:,2),Surf.SurfData.vertices(:,3)); R = double(unique(int32(R)));
  else
      
  end
for i =1:Ns
    Num = dot(n,repmat(P1(i,:),[size(n,1) 1]),2)+D; Den = dot(n,repmat(P2(i,:)-P1(i,:),[size(n,1) 1]),2);
    %Num = single(n(:,1)*P1(i,1)'+n(:,2)*P1(i,2)'+n(:,3)*P1(i,3)'+repmat(D,[1 size(P1,1)])); Den = single(n(:,1)*(P2(i,1)-P1(i,1))'+n(:,2)*(P2(i,2)-P1(i,2))'+n(:,3)*(P2(i,3)-P1(i,3))'+eps);
    t = -1*Num./Den;clear Num Den;
    xint = single(repmat(P1(i,1)',[size(Faces,1) 1])+t.*(repmat((P2(i,1)-P1(i,1))',[size(Faces,1) 1])));
    yint = single(repmat(P1(i,2)',[size(Faces,1) 1])+t.*(repmat((P2(i,2)-P1(i,2))',[size(Faces,1) 1])));
    zint = single(repmat(P1(i,3)',[size(Faces,1) 1])+t.*(repmat((P2(i,3)-P1(i,3))',[size(Faces,1) 1])));clear t;
    
    
     waitbar(i/Ns,H,['Computing Distance: Point ' num2str(i) ' of '  num2str(Ns)]);
    PpP1 =  Vert(Faces(:,1),:)-[xint yint zint];
    PpP2 =  Vert(Faces(:,2),:)-[xint yint zint];
    PpP3 =  Vert(Faces(:,3),:)-[xint yint zint];
    angP2p3 = (acos(dot(PpP2,PpP3,2)./(normm(PpP2).*normm(PpP3))))*180/pi;
    angP1p3 = (acos(dot(PpP1,PpP3,2)./(normm(PpP1).*normm(PpP3))))*180/pi;
    angP1p2 = (acos(dot(PpP1,PpP2,2)./(normm(PpP1).*normm(PpP2))))*180/pi;
    ind = find(round(angP2p3+angP1p3+angP1p2) == 360);
    [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord);
    Dmat = dist([Surft.SurfData.vertices(i,:);inte]');Dmat = Dmat(1,:);Dmat(Dmat==0) = [];
    %angs = acos(dot(Vertexint,tempinte,2)./(normm(Vertexint).*normm(tempinte)))
    if isempty(Dmat)
        Dmat = 0;
    else
        Dist(i) = min(Dmat);Int{i} = inte(Dmat==min(Dmat),:);Faceint =ind(Dmat==min(Dmat));
        Facet = unique(Faces(Faceint(:),:));
        Faceint =ind(Dmat==min(Dmat));Vertexint = Vert(Facet(:),:);
        
    end
    if ~isempty(strcmp(flag,'sphere'))
        Int{1} = Int{1} + sign(Int{1}).*abs((R-normm(repmat(Int{1},[3 1]))).*(dot([R 0 0;0 R 0; 0 0 R],repmat(Int{1},[3 1]),2))./(normm([R 0 0;0 R 0; 0 0 R]).*normm(repmat(Int{1},[3 1]))))';
    else
    end
    tempinte =repmat(Int{i},[size(Vertexint,1) 1]);
    Geodist = R*acos(dot(Vertexint,tempinte,2)./(normm(Vertexint).*normm(tempinte)));
    [Swei] = Geodesic_Smoothing(Geodist,'lineal');
    if strcmp(itype,'lineal')
        Icvar(i)=dot(cvar(Facet(:)),Swei)./normm(Swei);
    elseif strcmp(itype,'nearest')
        ind = find(Swei == min(Swei));
        Icvar(i)=cvar(Facet(ind(1)));
    end
end
 close(H)
return;

function norma = normm(M)
norma = sqrt(sum((M').^2))';
return