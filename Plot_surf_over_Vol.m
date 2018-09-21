function Plot_surf_over_Vol(Image,Surfa,CFiles,slices,ptsurf,cl,sw)
%
% Syntax :
% show_surf_over_image(Image,Surfa,slices,ptsurf)
%
% This function plots the surfaces contained in Surfa over the volume Image 
% for selected slices. 
%
% Input Parameters:
%   Image       : Image Filename.
%   Surfa       : Surfaces files.
%   CFiles      : Grosor or atlas text files.
%   slices      : Selected Slices to visualize surface intersection. Slices
%                 must be separated by , or ;(ie. 0 0 90,90 0 0;100 50 0,...)
%   ptsurf      : Variable to set surface view. If 'line', only the intersection 
%                 between the selected slices and surfaces will be plotted. 
%                 If 'patch', a 3D surface patch will be
%                 plotted over the selected slices.
%                 as a line. 
%    cl         : Colormap used to see the results.
%    sw         : Boolean variable(sa = 1, plot surfaces in the same window, 
%                 sa = 0, plot surfaces in diferent windows)
%
% Output Parameters:
%  
%
% Related references:
% 
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf
% Exp_Surf Plot_Surf
%__________________________________________________
% Authors: Pedro Valdes Hernandez and Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% April 30th 2007
% Version $1.0

%=====================Checking Input Parameters===========================%
if ~exist('cl','var')|(isempty(cl))
    cl = 'jet';
end
if ~exist('sw','var')|(isempty(sw))
    sw = 'y';
end
if ~exist('ptsurf','var')|(isempty(ptsurf))
    ptsurf = 'line';
end
if ~exist('Image','var')|(isempty(Image))
    [Image,sts] = spm_select([1],'image','Selecting Image File','',cd);
end
if ~exist('Surfa','var')|(isempty(Surfa))
    [Surfa,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
end
if ~exist('slices','var')|(isempty(slices))
    slices = input('Please enter the slices (ie, [ 89 12 65; 90 0 90; ...]):  ','s');
    ind=strfind(slices,',');
    if ~isempty(ind)
        slices(ind)=';';
    end
end
if ~exist('CFiles','var')|(isempty(CFiles))
    CFiles ='';
end
if exist('CFiles','var')&(~isempty(CFiles))
    Ns = size(Surfa,1);
    Nc = size(CFiles,1);
    if Ns~=Nc
        errordlg('Please Select as much text files as surfaces');
        return;
    end
end

%=========================================================================%
%=========================Main program====================================%
global defaults;
spm_defaults
warning off;
Osys=computer;
ind=strfind(slices,',');
if ~isempty(ind)
    slices(ind)=';';
end
slices = eval(['[' slices ']']);
ns = size(slices,1);
nr = ceil(sqrt(ns));
nc = nr;
V = spm_vol(Image);
I = spm_read_vols(V);
if isfield(defaults,'analyze')
    if (defaults.analyze.flip ==1)
        I = flipdim(I,1);
    end
end
I = permute(I,[2 1 3]); %I = zeros(size(I));
v = 90*eye(3);
%hw = waitbar(0,'printing overlaid surface...');
[pth,nm,ext] = fileparts(deblank(Image));
if strcmp(sw,'y')
    colordef black;h = figure('numbertitle','off','name',['IBASPM Plotting Surfaces over ' nm]);
end
for i = 1:ns
    cont=0;
    color = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    if strcmp(sw,'y')
        subplot(nr,nc,i);
    elseif strcmp(sw,'n')
        colordef black;h = figure('numbertitle','off','name',['IBASPM Plotting Surfaces over ' nm ', Slice ' num2str(slices(i,:))]);
    end
    hold on
    set(h,'Color',[0 0 0]);
    slicei = slices(i,:);
    fig = slice(I,slicei(1),slicei(2),slicei(3));
    slicei = logical(slicei);
    fign = fig(~slicei(1:3));
    figy = fig(slicei(1:3));
    if ~isempty(fign)&size(fign,1)==2
        set(fign(1),'Visible','off');
        set(fign(2),'Visible','off');
    elseif ~isempty(fign)&size(fign,1)==1
        set(fign(1),'Visible','off');
    end
    set(figy,'LineStyle','none');
    colormap gray; axis equal; hold on; 
    if sum(slicei,2)==1
        view(v(slicei,[1 3])); 
    else
        view(3)
    end
    h1 =gca;
    set(gca,'XTicklabel',num2str((str2num(get(gca,'XTicklabel')))-abs(V.mat(1,4))));set(h1,'XColor',[1 1 1]);xlabel(h1,' X axis(mm)');
    set(gca,'YTicklabel',num2str((str2num(get(gca,'YTicklabel')))-abs(V.mat(2,4))));set(h1,'YColor',[1 1 1]);ylabel(h1,' Y axis(mm)');
    set(gca,'ZTicklabel',num2str((str2num(get(gca,'ZTicklabel')))-abs(V.mat(3,4))));set(h1,'ZColor',[1 1 1]);zlabel(h1,' Z axis(mm)');
    set(gca,'Color',[0 0 0]);
    title(['Slice ' num2str(slices(i,:)) ' vox, ' num2str(slices(i,:) - [abs(V.mat(1,4)) abs(V.mat(2,4)) abs(V.mat(3,4))]) ' mm']);
    plane = double(slicei);
    wh = whos('Surfa');
    for s = 1:size(Surfa,1)
        if ischar(Surfa(s,:))
            [OutFiles, Surf] = Exp_Surf(deblank(Surfa(s,:)), Image, '','', 'imp','n');
            Surf = Surfa{1};
        elseif (strcmp(wh.class,'cell'))
            Surf = Surfa{s};
        elseif (strcmp(wh.class,'struct'))
            Surf = Surfa(s);
        end
        for j = 1:size(Surf,2)
            if ~isempty(CFiles)
                Npoints = size(Surf(j).SurfData.vertices,1);
                [pth,nm,ext] = fileparts(deblank(CFiles(s,:)));
                CFile = [pth filesep nm ext(1:4)];
                Surf(j).Is = single((textread(CFile,'%f',Npoints)));
            end
            cont=cont+1;
            if cont>size(color,1);color=repmat(color,[2 1]);end
            Surf(j).SurfData.vertices = [Surf(j).SurfData.vertices ones(size(Surf(j).SurfData.vertices,1),1)]*inv(V.mat)';
            Surf(j).SurfData.vertices(:,4) = [];
            if strcmp(ptsurf,'patch')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isfield(Surf(j),'Is')&~isfield(Surf(j).SurfData,'FaceVertexCData')
                    strsurf=patch(Surf(j).SurfData,'facecolor',color(cont,:),'edgecolor','none','tag', ...
                        'model0','facelighting','gouraud');
                elseif isfield(Surf(j),'Is')&~isfield(Surf(j).SurfData,'FaceVertexCData')
                    %% Plot using vertices
                    [Colors] = Surf_Color(Surf(j),cl);
                    Surf(j).SurfData.FaceVertexCData = Colors;
                    Surf(j).SurfData.FaceColor = 'interp';
                    strsurf=patch(Surf(j).SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
                    Surf(j).SurfData = rmfield(Surf(j).SurfData,'FaceVertexCData');
                elseif ~isfield(Surf(j),'Is')&(isfield(Surf(j).SurfData,'FaceVertexCData'))
                    Surf(j).SurfData.FaceColor = 'interp';
                    strsurf=patch(Surf(j).SurfData,'edgecolor','none','tag', 'patch','facelighting','gouraud');
                elseif isfield(Surf(j),'Is')&isfield(Surf(j).SurfData,'FaceVertexCData')
                    Surf(j).SurfData.FaceColor = 'interp';
                    strsurf=patch(Surf(j).SurfData,'edgecolor','none','tag', 'patch','facelighting','gouraud');
                end
            elseif strcmp(ptsurf,'line')
                ind = find(plane == 1);
                for k = 1:size(ind,2)
                    planet = zeros(1,4);planet(ind(k))=1;
                    planet(4) = -slices(i,ind(k));
                    [Xline,Yline,Zline] = cut_surf(planet,Surf(j).SurfData);
                    line(Xline,Yline,Zline,'LineWidth',1.5,'Color',color(cont,:));
                end
            end
            names{cont}=Surf(j).Name;
        end
        
    end
    %        waitbar(i/ns,hw)
    if strcmp(ptsurf,'patch')
        axis image;
        camlight;
    elseif strcmp(ptsurf,'line')
        axis image;
    end
end
%close(hw);
warning off;
%========================End of main program==============================%
return

%=========================================================================%
%=======================Internal functions================================%

function [Xline,Yline,Zline] = cut_surf(plane,fv)

% This function cuts a surface (defined by vertices and faces as represented 
% by MATLAB patches) by a plane, leading to a curve in 3D space. The
% resulting curve is represented by a set of contigous lines in the space
% 
% Syntax:
% [Xline,Yline,Zline] = cut_surf(plane,SurfData)
% 
% INPUTS:
% plane : A 4-length vector with the parameters of the plane. If plane = [A
% B C D] then every 3D point P = (x,y,z) belonging to the plane satisfies
% plane*[P; 1]' = A*x + B*y + C*z + D = 0
% SurfData : surface structure as represented in MATLAB by patches:
%            SurfData.vertices
%            SurfData.faces
%
% OUTPUTS:
% Xline,Yline,Zline: Matrices with the line coordinates.
% The entire curve can be plotted by simply typing: 
% line(Xline,Yline,Zline,'Properties1',Value1,...);
%
% Pedro Antonio Vald�s Hern�ndez
% 
% October 29 2008

warning off %#ok
Xline = []; Yline = []; Zline = [];
oo = ones(size(fv.vertices,1),1);
maxdist = sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
vertx = find(abs(dot([fv.vertices oo],plane(oo,:),2)/norm(plane(1:3)))<maxdist);
indf = ismember(fv.faces,vertx);
[rindf,cindf] = find(indf); %#ok
rindf = unique(rindf);
Nf = length(rindf);
% h = waitbar(0,'cutting surface...');
for i = 1:Nf
   verts = fv.vertices(fv.faces(rindf(i),:),:);
   verts(:,4) = 1;
   diffv(1,:) = diff(verts([1 2],:));
   diffv(2,:) = diff(verts([2 3],:));
   diffv(3,:) = diff(verts([3 1],:));
   alpha = -verts*plane'./(diffv*plane');
   % NaN   : contains
   % < 0   : not contains down
   % -Inf  : parallel down
   % > 1   : not contains up
   % +Inf  : parallel up
   ind = find((alpha<1 & alpha >=0) | (alpha<=1 & alpha >0))  ;
   if ~isempty(ind) && length(ind)==2
       points = verts(ind,1:3) + alpha(ind,[1 1 1]).*diffv(ind,1:3);
       Xline = [Xline points(:,1)]; %#ok
       Yline = [Yline points(:,2)]; %#ok
       Zline = [Zline points(:,3)]; %#ok
   end
%    waitbar(i/Nf,h);
end
% close(h);
warning on %#ok
%=========================================================================%
