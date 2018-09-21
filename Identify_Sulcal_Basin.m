function varargout = Identify_Sulcal_Basin(varargin);
%
% Syntax :
% hlist = Identify_Sulcal_Basin(Surf);
%
% This function plots each surface patch for identifying sulcal basin. The
% patch label must be specified in the .Is field.
%
% Input Parameters:
%   Surf       : Surfaces files.
%
% Output Parameters:
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

%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surf = varargin{1};

    Surf = Surface_Checking(Surf);
    if isfield(Surf.SurfData,'FaceVertexCData')
        Surf.SurfData = rmfield(Surf.SurfData,'FaceVertexCData');
    end
    % Parameters
    pauseTime = 1;
    
    % Figure Properties
    param.figcolor = [0 0 0]; % Color
    param.figurevisible = 'on'; % Visible
    param.newfigure = 1    ;      % Boolean variable to create new figure
    
    set(0,'units','centimeters');
    cm_screen = get(0,'screensize');
        
    figPosition = [1 1 cm_screen(3)-3 cm_screen(4)-3*cm_screen(3)/cm_screen(4)]; % Position
    figUnits = 'centimeters'; % Units
end
if nargin==2
    pauseTime = varargin{2};
end
if nargin >2 
    error('To many outputs');
    return;
end
%% ================== End of Checking input parameters ================================= %

%% ================================== Main Program ===================================== %

structs = unique(nonzeros(Surf.Is));
Nsulc = length(structs);


% New Figure
FigID = figure('numbertitle','off','name','New Figure','Color',[0 0 0], 'units',figUnits,'Position',figPosition,'InvertHardcopy','off');

for i = 1:Nsulc % For each sulcal spacelabSulbasin
    ind = find(Surf.Is == structs(i)); % Points belonging to a specified region
    Surftemp = Surf;
    Surftemp.Is = zeros(size(Surftemp.SurfData.vertices,1),1);
    Surftemp.Is(ind) = 1;
    
    axIds1 = subplot(1,2,1);
    h = title(['Sulcal Basin Label: ' num2str(structs(i))]);set(h,'Color',[[1 1 1]-param.figcolor]);
    set(axIds1,'Color',param.figcolor);
    set(axIds1,'XColor',[[1 1 1]-param.figcolor]*.7);
    set(axIds1,'YColor',[[1 1 1]-param.figcolor]*.7);
    set(axIds1,'ZColor',[[1 1 1]-param.figcolor]*.7);
    [hSurf] = Plot_Surf(Surftemp,'FigID',FigID);grid on;
    view([270 0 ]);
    
    axIds2 = subplot(1,2,2);
    h = title(['Sulcal Basin Label: ' num2str(structs(i))]);set(h,'Color',[[1 1 1]-param.figcolor]);
    set(axIds2,'Color',param.figcolor);
    set(axIds2,'XColor',[[1 1 1]-param.figcolor]*.7);
    set(axIds2,'YColor',[[1 1 1]-param.figcolor]*.7);
    set(axIds2,'ZColor',[[1 1 1]-param.figcolor]*.7);
    [hSurf] = Plot_Surf(Surftemp,'FigID',FigID);grid on;
    view([90 0 ]);
    pause(pauseTime);
    cla(axIds1,'reset');
    cla(axIds2,'reset');
end
%% ================================== Main Program ===================================== %
% Outputs


%% =========================== End of Main Program ===================================== %
return;