function Plot_Surf_Norm(SFile,Ptrl,sa);
%
% Syntax :
% Plot_Surf_Norm(SFile,Ptrl,sa);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa       : Surfaces files.
%               : An atlas surface file is considered as a single surface.
%   Ptrl        : Point index specifying which surface point will be plotted.
%   sa          : Boolean variable(sa = y, plot surfaces in the same window,
%                 sa = n, plot surfaces in diferent windows)
%
% Output Parameters:
%
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf quiver3
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 30th 2008
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin==0
    [SFile,sts] = spm_select([1],'any','Selecting Surface File','',cd);
end
if nargin < 2
    Ptrl = 'all';
end
if nargin < 3
    sa ='y';
end
%=========================================================================%

%=========================Main Program====================================%
alpha = 0.33;
beta = 0.33;
scf = 1;
wh = whos('SFile');
if (strcmp(wh.class,'struct'))
    Surfa = SFile;
    Ns = size(Surfa,1);
elseif(strcmp(wh.class,'cell'))
    Ns = size(SFile,1);
    Surfa = SFile;
elseif ischar(SFile(1,:));
    [OutFiles, Surf] = Exp_Surf(SFile, '0', '','', 'imp','n');
    Surfa = Surf{1};
    Ns = size(Surfa,1);
end
wh = whos('Surfa');
for j = 1:Ns
    if (strcmp(wh.class,'cell'))
        Surf = SFile{j};
    elseif (strcmp(wh.class,'struct'))
        Surf = Surfa(j);
    end
    indf = Sel_Structs(Ptrl,[1:size(Surf.SurfData.vertices,1)]');indf = indf(:);
    if max(indf)>size(Surf.SurfData.vertices,1)
        errordlg('Specified index are greater than the number of points in the surface');
        return;
    end
    fv.vertices = Surf.SurfData.vertices;
    fv.faces = Surf.SurfData.faces;
    if ~isfield(Surf.SurfData,'VertexNormals')
        h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
    else
        normals = Surf.SurfData.VertexNormals;
    end
    norma = sqrt(sum((normals').^2));
    normals = normals./repmat(norma',[1 3]);
    temp = sum(normals')';
    ind = find(isnan(temp) ==1);
    if ~isempty(ind)
        indt = ind;
        for i = 1:size(ind,1)
            vert = Surf.SurfData.faces(nonzeros(Surf.Tri(ind(i),3:end)));ind2 = ismember(vert,indt);vert(ind2) = [];
            tempN = sum([normals(vert,1) normals(vert,2) normals(vert,3)])/size(vert,1);
            normals(ind(i),:) = tempN;
            norma(ind(i)) = sqrt(sum((tempN').^2));
            indt(indt == ind(i)) = [];
        end
    end
    Surf.SurfData.VertexNormals = normals./repmat(norma',[1 3]);
    fv.vertices = Surf.SurfData.vertices;
    indfaces = ismember(Surf.SurfData.faces,indf);
    fv.faces = Surf.SurfData.faces((sum(indfaces')>0)',:);

    if strcmp(sa,'y')&(j==1)
        nm = Surf(1).Name;
        colordef black;h = figure('numbertitle','off','name',['IBASPM Surface Ploting...:  Plot ' nm]);hold on;  grid on; box on;
    elseif strcmp(sa,'n')
        if size(Surf,2)~=1
            nm = Surf(1).Name;
        else
            nm = Surf.Name;
        end
        colordef black;h = figure('numbertitle','off','name',['IBASPM Surface Ploting...:  Plot ' nm]);hold on;  grid on; box on;
    end
    strsurf=patch(fv,'facecolor',[1 0 0],'edgecolor','black','tag','model0','facelighting','gouraud');
    x = Surf.SurfData.vertices(indf,1); y = Surf.SurfData.vertices(indf,2); z = Surf.SurfData.vertices(indf,3);
    u = Surf.SurfData.VertexNormals(indf,1); v = Surf.SurfData.VertexNormals(indf,2); w = Surf.SurfData.VertexNormals(indf,3);
    if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
    delx = diff([min(x(:)) max(x(:))])/n;
    dely = diff([min(y(:)) max(y(:))])/m;
    delz = diff([min(z(:)) max(y(:))])/max(m,n);
    del = sqrt(delx.^2 + dely.^2 + delz.^2);
    if del>0
        len = sqrt((u/del).^2 + (v/del).^2 + (w/del).^2);
        maxlen = max(len(:));
    else
        maxlen = 0;
    end

    if maxlen>0
        scf = scf*0.9 / maxlen;
    else
        scf = scf*0.9;
    end
    col = abs([u v w]);
    col = col./repmat(sqrt(sum((col').^2))',[1 3]);
    u = u*scf; v = v*scf; w = w*scf;
    beta = beta * sqrt(u.*u + v.*v + w.*w) ./ (sqrt(u.*u + v.*v) + eps);
    hu = [x+u-alpha*(u+beta.*(v+eps)) x+u x+u-alpha*(u-beta.*(v+eps))];
    hv = [y+v-alpha*(v-beta.*(u+eps)) y+v y+v-alpha*(v+beta.*(u+eps))];
    hw = [z+w-alpha*w z+w z+w-alpha*w];
    for i = 1:size(u(:),1)
        line([x(i) x(i)+u(i)], [y(i) y(i)+v(i)],[z(i) z(i)+w(i)],'Color',col(i,:),'LineWidth',1);
        line([hu(i,1) hu(i,2) hu(i,3)], [hv(i,1) hv(i,2) hv(i,3)],[hw(i,1) hw(i,2) hw(i,3)],'Color',col(i,:),'LineWidth',1);
    end
    if ~strcmp(sa,'y')|(j==1)
        axis image;
        view(3);
        camlight;
    end
end
if strcmp(sa,'y')
    axis image;
    view(3);
end
%========================End of main program==============================%
return;