function Plot_Surf(Surfa,tr,sa);
%
% Syntax :
% Plot_Surf(Surfa,tr,sa);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa       : Surfaces files.
%               : An atlas surface file is considered as a single surface.
%   tr          : Transparency vector(values have to be between 0 and 1).
%   sa          : Boolean variable(sa = 1, plot surfaces in the same window,
%                 sa = 0, plot surfaces in diferent windows)
%
% Output Parameters:
%
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin==0
    [Surfa,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
    tr = input('Please select the transparence factor:   ');
    if (max(tr)>1)|(max(tr)<0)
        errordlg('Please select correct transparency values');
        return;
    end
    sa = input('Plot surfaces in the same window ? (y/n):   ','s');
    if (lower(sa)~='y')&(lower(sa)~='n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
end
if nargin < 2
    tr =1;
end
if nargin < 3
    sa ='y';
end
%=========================================================================%

%=========================Main Program====================================%
switch lower(sa)
    case 'y'
        sa = 1;
    case 'n'
        sa = 0;
end
[tr] = sort(tr);
wh = whos('Surfa');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = size(Surfa,1);
elseif ischar(Surfa(1,:));
    Ns = size(Surfa,1);
end
tr((size(tr,2)+1):Ns) = max(tr);
tr = tr(1:Ns);
[tr,Surfa] = Reord_Tr(Surfa,tr,Ns,wh);

%col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
cont = 0;
Osys = computer;
for j = 1:Ns
    Surf = Surfa{j,1};
    re = floor(length(Surf)/Ncolor); col = repmat(col,[re+1 1]);
    if ~isfield(Surf,'SurfData');
        errordlg('This file does not contain surface information');
        continue;
    end
    if (sa ==1)&(j==1)
        if ~isfield(Surf(1),'Name')
            Surf(1).Name = 'Unknown';
        end
        nm = Surf(1).Name;
        colordef black;h = figure('numbertitle','off','Color', 'black','name',['IBASPM Surface Ploting...:  Plot ' nm]);hold on;  grid on; box on;
    elseif (sa ==0)
        if ~isfield(Surf(1),'Name')
            Surf(1).Name = 'Unknown';
        end
        nm = Surf(1).Name;
        colordef black;h = figure('numbertitle','off','Color','black','name',['IBASPM Surface Ploting...:  Plot ' nm]);hold on;  grid on; box on;
    end
    if size(Surf,2)~=1
        for i = 1:size(Surf,2);
            if isfield(Surf(i).SurfData,'VertexNormals');
                Surf(i).SurfData = rmfield(Surf(i).SurfData,'VertexNormals');
            end
            strsurf=patch(Surf(i).SurfData,'facecolor',col(i,:),'edgecolor','none','tag', ...
                'model0','facelighting','gouraud');
            set(strsurf,'FaceAlpha',tr(j));
        end
    else
        if isfield(Surf.SurfData,'VertexNormals');
            Surf.SurfData = rmfield(Surf.SurfData,'VertexNormals');
        end
        if ~isfield(Surf,'Is')&~isfield(Surf.SurfData,'FaceVertexCData')
            strsurf=patch(Surf.SurfData,'facecolor',col(j,:),'edgecolor','none','tag', ...
                'model0','facelighting','gouraud');
        elseif isfield(Surf,'Is')&~isfield(Surf.SurfData,'FaceVertexCData')
            %% Plot using vertices
            [Colors] = Surf_Color(Surf,'jet');
            Surf.SurfData.FaceVertexCData = Colors;
            Surf.SurfData.FaceColor = 'interp';
            strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
            if sum(Surf.Is-floor(Surf.Is)) ~=0
                h = colorbar;
                mat = max(Surf(1).Is);
                set(h,'YTickLabel',num2str(str2num(get(h,'YTickLabel'))*mat));
            end
            Surf.SurfData = rmfield(Surf.SurfData,'FaceVertexCData');
        elseif ~isfield(Surf,'Is')&(isfield(Surf.SurfData,'FaceVertexCData'))
            Surf.SurfData.FaceColor = 'interp';
            strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','gouraud');
        end
    end
    if isfield(Surf(1),'Is')&isfield(Surf(1).SurfData,'FaceVertexCData')
        cont = cont + 1;
        Surf.SurfData.FaceColor = 'interp';
        strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','gouraud');
        if sum(Surf(1).Is-floor(Surf(1).Is)) ~=0
            h = colorbar;
            mat = max(Surf(1).Is);
            set(h,'YTickLabel',num2str(str2num(get(h,'YTickLabel'))*mat));
        end
    end
    if (sa == 0)|(j==1)
        axis image;
        view(3);
        if ~strcmp(Osys,'PCWIN')|(cont == 0)
            camlight;
        end
    end
    set(strsurf,'FaceAlpha',tr(j));
end
if (sa == 1)
    axis image;
    view(3);
end
%========================End of main program==============================%
return

%=======================Internal functions================================%

function [trn,Surfa] = Reord_Tr(SurfF,tr,Ns,wh);

if strcmp(wh.class,'struct')
    Surfa{1,1} = SurfF;
elseif strcmp(wh.class,'cell')
    Surfa = SurfF;
end

for j = 1:Ns
    if ischar(SurfF(1,:));
        [ pth nm ext] = fileparts(SurfF(j,:));
        if strcmp(ext(1:4),'.mat');
            Surf = load('-mat',[pth filesep nm ext(1:4)]);
            Surfa{j,1} = Surf.Surf;
        else
            [OutFiles, Surf] = Exp_Surf(SurfF(j,:), '0', '','', 'imp','n');
            Surfa{j,1} = Surf{1};
        end
    end
    Surf = Surfa{j,1};
    for i = 1:size(Surf,2);
nam(j) =abs(min(Surf(i).SurfData.vertices(:,1)));
    end
end
[b,t] = sort(nam);
trn(t(end:-1:1)) = tr;

return
%=========================================================================%
%=========================Internal Functions===============================





