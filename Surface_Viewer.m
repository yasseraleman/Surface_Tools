function varargout = test_Surface_Viewer(Surfa,varargin);
%
% Syntax :
% hlist = Plot_Surf(Surfa,transpVal,colMap);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa       : Surfaces files.
%               : An atlas surface file is considered as a single surface.
%   transpVal          : Transparency vector(values have to be between 0 and 1).
%   colMap          : Colormap used to see the results.
%
% Output Parameters:
%   hlist       : Handles list.
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
global hlink;
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    
    % Parameters
    transpVal = 1; % Opacity
    colMap = 'jet'; % ColorMap
    newFig = 'y'; %
    
    % Figure Properties
    param.figcolor = [0 0 0]; % Color
    param.figurevisible = 'on'; % Visible
    param.newfigure = 1    ;      % Boolean variable to create new figure
    
    set(0,'units','centimeters');
    cm_screen = get(0,'screensize');
    
    figPosition = [1 1 cm_screen(3)-3 cm_screen(4)-3*cm_screen(3)/cm_screen(4)]; % Position
    figUnits = 'centimeters'; % Units
end
%% =================== End of checking input parameters ================= %


% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'transpVal' % Group Separation
                    transpVal=varargin{2};
                case 'colMap' % Jitter compression percent
                    colMap=varargin{2};
                case 'figcolor'
                    param.figcolor=varargin{2};
                case 'figurevisible'
                    param.figurevisible=varargin{2};
                case 'figUnits'
                    figUnits=varargin{2};
                case 'figPosition'
                    figPosition=varargin{2};
                case 'figUnits'
                    figUnits=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% =================== End of checking input parameters ================= %

%% ========================= Main Program =============================== %

% New Figure
FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
menuId = uimenu(FigID,'Label','Surface Editing'); % Creating Menu  for figure.


Nsurf = length(Surfa);
% if Nsurf>6
%     Nsurf = 6;
% end

nc = ceil(sqrt(Nsurf));
nr = ceil(Nsurf/nc);
clmaps = colormaps_colors;
for i = 1:Nsurf
    temp = Surfa{i,1};
    clmap = deblank(clmaps(i,:));
    axIds(i) = subplot(nr,nc,i);
    set(axIds(i),'Color',param.figcolor);
    set(axIds(i),'XColor',[[1 1 1]-param.figcolor]*.7);
    set(axIds(i),'YColor',[[1 1 1]-param.figcolor]*.7);
    set(axIds(i),'ZColor',[[1 1 1]-param.figcolor]*.7);
    grid on;
    if ischar(temp)
        if exist(deblank(temp(1,:)),'file')
            for j = 1:size(temp,1);
                tempSurf = deblank(temp(j,:));
                Surf(j) = Surface_Checking(tempSurf);
            end
            [hSurf, jSurf, surfRealMap] = test_Plot_Surf(temp,'FigID',FigID,'transpVal',1,'colMap',clmap);grid on;
        else
            error('Wrong surface filename');
            return
        end
    else
        [hSurf, jSurf, surfRealMap] = test_Plot_Surf(temp,'FigID',FigID,'transpVal',1,'colMap',clmap);grid on;
    end
    
    hSurfs{i} = hSurf; % All handles
    jSurfs{i} = jSurf; % All Surfaces
    surfRealMaps{i} = surfRealMap; % Index vector of all real values for each surface map
    
    % Detecting Colormap
    bb = findobj(hSurf,'Type','patch');
    cont = 0;
    for j = 1:length(bb)
        if isempty(bb(j).FaceVertexCData)
            cont = cont + 1;
        end
    end
    if cont == length(bb)
        clmap = 'Custom color...';
    end
    %     axIds(i) = gca;
    axPropert.PropFig = '';
    axPropert.ColorbarFig = '';
    axPropert.Colormap = clmap;
    axPropert.colorbaron = 1;
    axPropert.gridon = 1;
    axPropert.unlinkon = 0;
    axPropert.colnormalizeon = 0;
    axPropert.Limits = [axIds(i).XLim(1) axIds(i).XLim(2)...
        axIds(i).YLim(1) axIds(i).YLim(2)...
        axIds(i).ZLim(1) axIds(i).ZLim(2)];
    
    %% Title Properties
    axIds(i).Title.String = ['Axis ' num2str(i)];
    axIds(i).Title.Color = [1 1 1] - param.figcolor;
    axPropert.Title.titleon = 0;
    axPropert.Title.String = ['Introduce New Title for Axis ' num2str(i)] ;
    axPropert.Title.Enable = 'off';
    % Font Properties Title
    axPropert.Title.FontName = axIds(i).Title.FontName;
    axPropert.Title.FontSize = axIds(i).Title.FontSize;
    axPropert.Title.FontWeight = axIds(i).Title.FontWeight;
    axPropert.Title.FontColor = axIds(i).Title.Color;
    axPropert.Title.FontAngle = axIds(i).Title.FontWeight;
    
    
    %% XLabel Properties
    axPropert.XLabel.titleon = 0;
    axPropert.XLabel.String = 'Introduce New X Axis Label';
    axPropert.XLabel.Enable = 'off';
    % Font Properties XLabel
    axPropert.XLabel.FontName = axIds(i).XLabel.FontName;
    axPropert.XLabel.FontSize = axIds(i).XLabel.FontSize;
    axPropert.XLabel.FontWeight = axIds(i).XLabel.FontWeight;
    axPropert.XLabel.FontColor = axIds(i).XLabel.Color;
    axPropert.XLabel.FontAngle = axIds(i).XLabel.FontWeight;
    axPropert.XAxis.gridon = axIds(i).XGrid;
    axPropert.XAxis.Color = axIds(i).XColor;
    axPropert.XAxis.numXTicks = length(axIds(i).XTick);
    
    
    %% YLabel Properties
    axPropert.YLabel.titleon = 0;
    axPropert.YLabel.String = 'Introduce New Y Axis Label';
    axPropert.YLabel.Enable = 'off';
    % Font Properties YLabel
    axPropert.YLabel.FontName = axIds(i).YLabel.FontName;
    axPropert.YLabel.FontSize = axIds(i).YLabel.FontSize;
    axPropert.YLabel.FontWeight = axIds(i).YLabel.FontWeight;
    axPropert.YLabel.FontColor = axIds(i).YLabel.Color;
    axPropert.YLabel.FontAngle = axIds(i).YLabel.FontWeight;
    axPropert.YAxis.gridon = axIds(i).YGrid;
    axPropert.YAxis.Color = axIds(i).YColor;
    axPropert.YAxis.numXTicks = length(axIds(i).YTick);
    
    
    %% ZLabel Properties
    axPropert.ZLabel.titleon = 0;
    axPropert.ZLabel.String = 'Introduce New Z Axis Label';
    axPropert.ZLabel.Enable = 'off';
    % Font Properties ZLabel
    axPropert.ZLabel.FontName = axIds(i).ZLabel.FontName;
    axPropert.ZLabel.FontSize = axIds(i).ZLabel.FontSize;
    axPropert.ZLabel.FontWeight = axIds(i).ZLabel.FontWeight;
    axPropert.ZLabel.FontColor = axIds(i).ZLabel.Color;
    axPropert.ZLabel.FontAngle = axIds(i).ZLabel.FontWeight;
    axPropert.ZAxis.gridon = axIds(i).ZGrid;
    axPropert.ZAxis.Color = axIds(i).ZColor;
    axPropert.ZAxis.numXTicks = length(axIds(i).ZTick);
    
    axIds(i).UserData = axPropert;
    
    cbar = findobj(hSurfs{i},'Tag','Colorbar');
    if ~isempty(cbar)
        
        colorbarProper.Title.titleon = 1;
        colorbarProper.Title.String = 'Colorbar Title';
        tempHand = title(cbar, colorbarProper.Title.String);
        colorbarProper.Title.Handle = tempHand;
        colorbarProper.Title.Handle.Color = [1 1 1] - param.figcolor;
        colorbarProper.Title.Enable = 'on';
        colorbarProper.Ticks.number = 11;
        colorbarProper.FontName = cbar.FontName;
        colorbarProper.FontWeight = cbar.FontWeight;
        colorbarProper.FontSize = cbar.FontSize;
        colorbarProper.FontAngle = cbar.FontAngle;
        colorbarProper.FontColor = cbar.Color;
        cbar.UserData = colorbarProper;
    end
    
    
    %submenuId(i) = uimenu(menuId,'Label',['Axis: ' num2str(i)],'Callback',{@change_axis_properties,axIds(i),hSurfs(i)}); % Creating submenus
    % submenuId(i) = uimenu(menuId,'Label',['Axis: ' num2str(i)],'Callback',{@change_axis_properties,axIds(i),hSurfs(i),jSurfs(i),surfRealMaps(i)}); % Creating submenus
    
    
    
    % %     contextmenuId = uicontextmenu(FigID); % Creating Context Menu
    % %     axIds(i).UIContextMenu = contextmenuId;
    % %     submenuIdc(i) = uimenu(contextmenuId,'Label',['Axis: ' num2str(i)],'Callback',{@change_axis_properties,axIds(i),hSurf}); % Creating submenus
    
end
if Nsurf > 1
    for i = 1:Nsurf
        submenuId(i) = uimenu(menuId,'Label',['Axis: ' num2str(i)],'Callback',{@change_axes_properties,axIds,hSurfs,jSurfs,surfRealMaps,i}); % Creating submenus for editing all axes at the same time
    end
    axNumb = 1:Nsurf;
    submenuId(i+1) = uimenu(menuId,'Label','All Axes','Callback',{@change_axes_properties,axIds,hSurfs,jSurfs,surfRealMaps,axNumb}); % Creating submenus for editing all axes at the same time
    
else
    axNumb = 1;
    submenuId(1) = uimenu(menuId,'Label',['All Axes'],'Callback',{@change_axes_properties,axIds,hSurfs,jSurfs,surfRealMaps, axNumb}); % Creating submenus for editing all axes at the same time
end
hlink = linkprop(axIds(:)',{'CameraPosition','CameraUpVector'});
%% ======================== End of Main Program ========================= %
% Outputs
varargout{1} = hlink;
varargout{2} = axIds;
return;


function FigpropID = change_axes_properties(h, callbackdata, axIds, hSurfs, jSurfs, surfRealMaps, axNumb);
% close all;
axPropert = axIds(axNumb(1)).UserData; % Restoring axis properties. It allows to keep configuration in case the the user closes the properties window

figID = axIds(axNumb(1)).Parent;

figUnits = 'centimeters';
figColor = [0.3 0.3 0.3];
set(0,'units',figUnits);
cm_screen = get(0,'screensize');

% % Creating Figure Dimensions
scaleFactorWidth = 3; % Figure Width scale factor
scaleFactorHeight = 2; % Figure Height scale factor

figWidth  = (cm_screen(3)-cm_screen(1))/scaleFactorWidth; % Figure width
figHeight = (cm_screen(4)-3*cm_screen(3)/cm_screen(4))/scaleFactorHeight; % Figure height
figCornerX = cm_screen(3)/2-(figWidth/2); % Figure Corner (X coordinate)
figCornerY = cm_screen(4)/2-(figHeight/2); % Figure Corner (Y coordinate)

figPosition = [figCornerX figCornerY figWidth figHeight]; % Figure Position
if ~isempty(axPropert.PropFig)&&ishandle(axPropert.PropFig) % Detecting and closing old figures related with this axis
    close(axPropert.PropFig);
end

% Creating figure
hand.FigpropID = figure('numbertitle',    'off',...
    'name',           'Modifying all axes at the same time ...',...
    'Color',          figColor,...
    'units',          figUnits,...
    'Position',       figPosition,...
    'InvertHardcopy', 'off',...
    'MenuBar',        'None',...
    'Resize',         'off');

% Saving Figure Handle in the axPropert variable
axPropert.PropFig =  hand.FigpropID;
axIds(axNumb(1)).UserData = axPropert;

%% %%%%%%%%%%%%%%%%% Defining Left Panel Objects %%%%%%%%%%%%%%%%%%%%%%%% %

% Creating Left Panel Properties
hand.panelID(1) = uipanel(hand.FigpropID,...
    'Position',            [.01 .16 .48 .74],...
    'Units',               'normalized',...
    'Title',               'Axis properties',...
    'BackgroundColor',     figColor,...
    'HighlightColor',      figColor/4,...
    'ForegroundColor',     [1 1 1]);
axNames = '';
for itemp = 1:length(axIds)
    axNames = strvcat(axNames, ['Modifying Axis: ' num2str(itemp)]);
end
if length(axIds) >1
    axNames = strvcat(axNames, ['Modifying all axes']);
end

hand.popupID = uicontrol(hand.FigpropID,...
    'Style',               'popup',...
    'Tag',                 'axispopup',...
    'String',              cellstr(axNames),...
    'Value',               axNumb(1),...
    'Units',               'normalized',...
    'Position',            [.02 .87 .96 .11],...
    'Callback',            {@getaxis,axIds, hSurfs, jSurfs, surfRealMaps});

%% ====================== Defining Colormap ============================= %

% Colormaps Text
hand.textID(1) = uicontrol(hand.panelID(1),...
    'Style',               'text',...
    'String',              'Select Colormap',...
    'Units',               'normalized',...
    'Position',            [.025 .92 .44 .045],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0]);

% Creating colormaps list
defColormaps = strvcat('Custom color...' ,colormaps_colors);
colormaps = cellstr(defColormaps);

% Detecting posible colormaps
indposclmap = find(ismember(defColormaps(:,1:length(axPropert.Colormap)),axPropert.Colormap,'rows'));

% Colormaps Popup
% % hand.popupID = uicontrol(hand.panelID(1),...
% %               'Style',               'popup',...
% %               'String',              colormaps,...
% %               'Value',               indposclmap,...
% %               'Units',               'normalized',...
% %               'Position',            [.025 .83 .95 .07],...
% %               'Callback',            {@setmap,axId,hSurfs});

hand.popupID = uicontrol(hand.panelID(1),...
    'Style',               'popup',...
    'Tag',                 'colormappopup',...
    'String',              colormaps,...
    'Value',               indposclmap,...
    'Units',               'normalized',...
    'Position',            [.025 .83 .95 .07],...
    'Callback',            {@setmap,axIds, hSurfs, jSurfs, surfRealMaps});



set(hand.popupID,'FontName','Courier','BackgroundColor',[1 1 1]);
set(hand.popupID,'Units','characters');
menuPos = get(hand.popupID,'Position');
set(hand.popupID,'Units','pixels');

allLength = cellfun(@numel,colormaps);
maxLength = max(allLength);
cmapHTML = [];
Npopcolor = 8;
for i = 1:numel(colormaps)
    arrow = [repmat('-',1,maxLength-allLength(i)+1) '>'];
    if i== 1
        
        %cmapFun = str2func(['@(x) ' lower(colormaps{i}) '(x)']);
        cData = repmat([1 1 1],[Npopcolor 1]);
    else
        cData = abs(colormaps_colors(deblank(lower(colormaps{i})),Npopcolor));
        %         cData = cmapFun(16);
    end
    curHTML = ['<HTML>' colormaps{i} '<FONT color="#FFFFFF">' arrow '<>'];
    for nc = 1:Npopcolor
        HEX = rgb2hex([cData(nc,1) cData(nc,2) cData(nc,3)]);
        curHTML = [curHTML '<FONT BGCOLOR="' HEX '" COLOR="' HEX '">__'];
    end
    
    
    curHTML = curHTML(1:end-2);
    curHTML = [curHTML '</FONT></HTML>'];
    cmapHTML = [cmapHTML; {curHTML}];
end
set(hand.popupID,'String',cmapHTML);

%% ==================== End of Defining Colormap ======================== %

%% ======================== Defining Views  ============================= %

% Camera View Panel Angle Text
hand.panelID(2) = uipanel(hand.panelID(1),...
    'Position',            [.025 .5 .95 .31],...
    'Title',               'Camera Point of View',...
    'BackgroundColor',     figColor,...
    'HighlightColor',      figColor/4,...
    'ForegroundColor',     [1 1 1]);

% Camera View text
hand.textID(2) = uicontrol(hand.panelID(1),...
    'Style',               'text',...
    'String',              'Camera Angles ',...
    'Units',               'normalized',...
    'Position',            [0.05    0.63    0.5    0.05],...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0],...
    'HorizontalAlignment', 'left');

% Azimuthal Angle Text
hand.textID(3) = uicontrol(hand.panelID(1),...
    'Style',               'text',...
    'String',              'Azimuthal',...
    'Units',               'normalized',...
    'Position',            [0.435    0.69    0.25    0.05],...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [1 1 1],...
    'HorizontalAlignment', 'left');

% Elevation Angle Text
hand.textID(4) = uicontrol(hand.panelID(1),...
    'Style',               'text',...
    'String',              'Elevation',...
    'Units',               'normalized',...
    'Position',            [0.71    0.69    0.23    0.05],...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [1 1 1],...
    'HorizontalAlignment', 'left');

% Obtain current axis view
[azAngle,elevAngle] = view(axIds(axNumb(1)));

% Azimuthal Angle Edit
hand.editID(1) = uicontrol(hand.panelID(1),...
    'Tag',                 'azimuthaledit',...
    'Style',               'edit',...
    'String',              num2str(azAngle),...
    'Units',               'normalized',...
    'Position',            [0.455    0.63    0.20    .055],...
    'Callback',            {@setazimuthangle,axIds});

% Elevation Angle Edit
hand.editID(2) = uicontrol(hand.panelID(1),...
    'Tag',                 'elevationedit',...
    'Style',               'edit',...
    'String',              num2str(elevAngle),...
    'Units',               'normalized',...
    'Position',            [0.72    0.63    0.20   .055],...
    'Callback',            {@setelevationangle,axIds});

% Reset Axis View Button
hand.buttonID(1) = uicontrol(hand.panelID(1),...
    'Style',               'pushbutton',...
    'Value',               0,...
    'String',              'Reset Camera View',...
    'Units',               'normalized',...
    'Position',            [0.06 0.53 0.87 0.07],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0],...
    'Callback',            {@resetcameraview, axIds});
% %
%% ==================== End of Defining Views =========================== %
% %
% %
%% =================== Colorbar properties checkboxes =================== %
% Colorbar Panel
hand.panelID(3) = uipanel(hand.panelID(1),...
    'Position',        [.025 .20 .95 .27],...
    'Title',           'Colorbar Properties',...
    'BackgroundColor', figColor,...
    'HighlightColor',  figColor/4,...
    'ForegroundColor', [1 1 1]);

% Display Colorbar Checkbox
hand.checkID(1) = uicontrol(hand.panelID(1),...
    'Tag',                 'showcolorbarcheck',...
    'Style',               'checkbox',...
    'Value',               axPropert.colorbaron,...
    'String',              'Display colorbar',...
    'Units',               'normalized',...
    'Position',            [.07 .33 .45 .07],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0],...
    'Callback',            {@switchcolobar, axIds, hSurfs});

% Edit ColorBar Properties
hand.buttonID(2) = uicontrol(hand.panelID(1),...
    'Style',               'pushbutton',...
    'Value',               0,...
    'String',              'Change Colorbar Properties',...
    'Units',               'normalized',...
    'Position',            [0.06 0.23 0.87 0.07],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0],...
    'Callback',            {@changecolorbar, axIds, hSurfs});



%% ===================== End of Colorbar properties ===================== %

%% ======================== Some checkboxes ============================= %
% hand.panelID(4) = uipanel(hand.panelID(1),'Position',[.025 .02 .95 .15],'Title','Some Checkboxes','BackgroundColor',figColor,'HighlightColor',figColor/4, 'ForegroundColor',[1 1 1]);
% Display Axis Grid Checkbox
hand.checkID(2) = uicontrol(hand.panelID(1),...
    'Tag',                 'showaxisgridcheck',...
    'Style',               'checkbox',...
    'Value',               axPropert.gridon,...
    'String',              'Display axis grid',...
    'Units',               'normalized',...
    'Position',            [.04 .03 .46 .07],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0],...
    'Callback',            {@switchaxisgrid, axIds});

% Unlink Axis Checkbox
hand.checkID(3) = uicontrol(hand.panelID(1),...
    'Tag',                 'unlinkpropercheck',...
    'Style',               'checkbox',...
    'Value',               axPropert.unlinkon,...
    'String',              'Unlink axis',...
    'Units',               'normalized',...
    'Position',            [.04 .1 .44 .07],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [0 1 1],...
    'Callback',            {@unlinkaxis, axIds});


if length(axNumb) > 1
    axPropert.colnormalizevisible = 'on';
else
    axPropert.colnormalizevisible = 'off';
end
axIds(axNumb(1)).UserData = axPropert;
hand.checkID(3) = uicontrol(hand.panelID(1),...
    'Tag',                 'colnormalizecheck',...
    'Style',               'checkbox',...
    'Visible',             axPropert.colnormalizevisible,...
    'Value',               axPropert.colnormalizeon,...
    'String',              'Normalize colormaps',...
    'Units',               'normalized',...
    'Position',            [.49 .1 .44 .07],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [0 1 1],...
    'Callback',            {@normalizecolormaps, axIds, hSurfs, jSurfs, surfRealMaps});
%% ==================== End of Some checkboxes ============================= %

%% ========================== Export Figure ============================= %
% Export button panel
hand.panelID(4) = uipanel(hand.FigpropID, ...
    'Position',        [.01 .02 .48 .13],...
    'Units',           'normalized',...
    'BackgroundColor', figColor,...
    'ForegroundColor', [1 0 0]);

% Export Button
hand.buttonID(3) = uicontrol(hand.panelID(4), ...
    'Style',           'pushbutton',...
    'Value',           0,...
    'String',          'Export Figure',...
    'FontName',        'Helvetica',...
    'FontSize',        10,...
    'FontUnits',       'points',...
    'FontWeight',      'bold',...
    'FontAngle',       'normal',...
    'Units',           'normalized',...
    'Position',        [.001 .001 .998 .998],...
    'BackgroundColor', [0.6902    0.1647    0.1647],...
    'ForegroundColor', [1 1 1],...
    'Callback',        {@exportfigure, figID});

%% ==================== End of Export Figure ============================ %


%% %%%%%%%%%%%%%%%%% End of Left Panel Objects %%%%%%%%%%%%%%%%%%%%%%%%%% %

%% %%%%%%%%%%%%%%%%% Defining Right Panel Objects %%%%%%%%%%%%%%%%%%%%%%% %

% Creating Uipanel Properties
hand.panelID(5) = uipanel(hand.FigpropID,...
    'Position',        [.51 .01 .48 .89],...
    'Units',           'normalized',...
    'Title',           'Axis Appearance',...
    'BackgroundColor', figColor,...
    'HighlightColor',  figColor/4,...
    'ForegroundColor', [1 1 1]);

% Creating Tab Group
s = warning('off', 'MATLAB:uitabgroup:OldVersion');
hand.hTabGroup(1) = uitabgroup('Parent',    hand.panelID(5),...
    'Position',  [.0225 .0225 .96 .955],...
    'Units', 'normalized');

%% ==================== Creating Tab 1 (Axis Tab) ======================= %
hTabs(1) = uitab('Parent',            hand.hTabGroup(1),...
    'Title',             'Axis   ',...
    'BackgroundColor',   figColor,...
    'ForegroundColor',   figColor);

% Axis Title
hand.checkID(4) = uicontrol('Parent',              hTabs(1),...
    'Tag',                 'showtitlecheck',...
    'Style',               'checkbox',...
    'Value',               axPropert.Title.titleon,...
    'String',              'Display Axis Title',...
    'Units',               'normalized',...
    'Position',            [.045 0.92 0.7 0.045],...
    'HorizontalAlignment', 'left',...
    'BackgroundColor',     figColor,...
    'ForegroundColor',     [.81 .81 0],...
    'Callback',            {@setnewaxistitle,axIds});

hand.editID(3) = uicontrol('Parent',               hTabs(1),...
    'Tag',                  'showtitleedit',...
    'Style',                'edit',...
    'Enable',               axPropert.Title.Enable,...
    'String',               axPropert.Title.String,...
    'Units',                'normalized',...
    'Position',             [.045 .84 .91 .06],...
    'Callback',             {@introduceaxistitle,axIds});

% Axis Limits
hand.panelID(6) = uipanel(hTabs(1),...
    'Position',              [.02   0.300    0.9600    0.5000],...
    'Title',                 'Axis Limits',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

hand.textID(5) = uicontrol(hTabs(1),...
    'Style',                 'text',...
    'String',                'Min',...
    'Units',                 'normalized',...
    'Position',              [.52 .71 .18 .055],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [1 1 1]);

hand.textID(6) = uicontrol(hTabs(1),...
    'Style',                 'text',...
    'String',                'Max',...
    'Units',                 'normalized',...
    'Position',              [.76 .71 .18 .055],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [1 1 1]);

minXlim = axIds(axNumb(1)).XLim(1);
maxXlim = axIds(axNumb(1)).XLim(2);
hand.textID(7) = uicontrol(hTabs(1),...
    'Style',                 'text',...
    'String',                'X Axis Limits',...
    'Units',                 'normalized',...
    'Position',              [0.025    0.645    0.5    0.0450],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);

hand.editID(4) = uicontrol(hTabs(1),...
    'Tag',                   'xMinEdit',...
    'Style',                 'edit',...
    'String',                num2str(minXlim),...
    'Units',                 'normalized',...
    'Position',              [.52 .64 .18 .055],...
    'Callback',              {@setaxminlim,axIds});

hand.editID(5) = uicontrol(hTabs(1),...
    'Tag',                   'xMaxEdit',...
    'Style',                 'edit',...
    'String',                num2str(maxXlim),...
    'Units',                 'normalized',...
    'Position',              [.76 .64 .18 .055],...
    'Callback',              {@setaxmaxlim,axIds});


minYlim = axIds(axNumb(1)).YLim(1);
maxYlim = axIds(axNumb(1)).YLim(2);
hand.textID(8) = uicontrol(hTabs(1),...
    'Style',                 'text',...
    'String',                'Y Axis Limits',...
    'Units',                 'normalized',...
    'Position',              [0.025    0.545    0.5    0.0450],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);

hand.editID(6) = uicontrol(hTabs(1),...
    'Tag',                   'yMinEdit',...
    'Style',                 'edit',...
    'String',                num2str(minYlim),...
    'Units',                 'normalized',...
    'Position',              [.52 .54 .18 .055],...
    'Callback',              {@setayminlim,axIds});

hand.editID(7) = uicontrol(hTabs(1),...
    'Tag',                   'yMaxEdit',...
    'Style',                 'edit',...
    'String',                num2str(maxYlim),...
    'Units',                 'normalized',...
    'Position',              [.76 .54 .18 .055],...
    'Callback',              {@setaymaxlim,axIds});


minZlim = axIds(axNumb(1)).ZLim(1);
maxZlim = axIds(axNumb(1)).ZLim(2);

hand.textID(9) = uicontrol(hTabs(1),...
    'Style',                 'text',...
    'String',                'Z Axis Limits',...
    'Units',                 'normalized',...
    'Position',              [0.025    0.445    0.5    0.0450],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);

hand.editID(10) = uicontrol(hTabs(1),...
    'Tag',                   'zMinEdit',...
    'Style',                 'edit',...
    'String',                num2str(minZlim),...
    'Units',                 'normalized',...
    'Position',              [.52 .44 .18 .055],...
    'Callback',              {@setazminlim,axIds});

hand.editID(11) = uicontrol(hTabs(1),...
    'Tag',                   'zMaxEdit',...
    'Style',                 'edit',...
    'String',                num2str(maxZlim),...
    'Units',                 'normalized',...
    'Position',              [.76 .44 .18 .055],...
    'Callback',              {@setazmaxlim,axIds});

hand.buttonID(4) = uicontrol(hTabs(1),...
    'Style',                 'pushbutton',...
    'Value',                 0,...
    'String',                'Reset Axis Limits',...
    'Units',                 'normalized',...
    'Position',              [0.06 0.33 0.87 0.085],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@resetaxislimits, axIds});

% Font Properties
hand.panelID(7) = uipanel(hTabs(1),...
    'Position',              [.02   0.02    0.9600    0.26],...
    'Title',                 'Title Font Properties',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

hand.buttonID(5) = uicontrol(hTabs(1),...
    'Style',                 'pushbutton',...
    'Tag',                   'titlefonttype',...
    'Value',                 0,...
    'String',                [axPropert.Title.FontName ', '  axPropert.Title.FontWeight ', '  axPropert.Title.FontAngle ' (' num2str(axPropert.Title.FontSize) ' points)'],...
    'Units',                 'normalized',...
    'Position',              [0.042 0.14 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@changefonttitleprop, axIds});

hand.buttonID(6) = uicontrol(hTabs(1),...
    'Style',                 'pushbutton',...
    'Tag',                   'titlefontcolorbutton',...
    'Value',                 0,...
    'String',                'Change Font Color...',...
    'Units',                 'normalized',...
    'Position',              [0.042 0.04 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       axPropert.Title.FontColor,...
    'ForegroundColor',       [1 1 1] - axPropert.Title.FontColor,...
    'Callback',              {@changefonttitlecolor, axIds});


%% ================== End of Creating Tab 1 (Axis Tab) ================== %


%% ==================== Creating Tab 2 (X Axis Tab) ===================== %
hTabs(2) = uitab('Parent',            hand.hTabGroup(1),...
    'Title',             'X Axis   ',...
    'BackgroundColor',   figColor,...
    'ForegroundColor',   figColor);

hand.checkID(5) = uicontrol('Parent',                hTabs(2),...
    'Tag',                   'showxlabelcheck',...
    'Style',                 'checkbox',...
    'Value',                 axPropert.XLabel.titleon,...
    'String',                'Display X Axis Label',...
    'Units',                 'normalized',...
    'Position',              [.045 0.92    0.7    0.045],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@setnewaxisxlabel,axIds});


hand.editID(10) = uicontrol('Parent',                hTabs(2),...
    'Tag',                   'showxlabeledit',...
    'Style',                 'edit',...
    'Enable',                axPropert.XLabel.Enable,...
    'String',                axPropert.XLabel.String,...
    'Units',                 'normalized',...
    'Position',              [.045 .84 .91 .06],...
    'Callback',              {@introduceaxisxlabel,axIds});

% X Axis appearance
hand.panelID(8) = uipanel(hTabs(2),...
    'Position',              [.02   0.31    0.9600    0.5],...
    'Title',                 'X Axis appearance',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

switch axPropert.XAxis.gridon
    case 'on'
        tempValue = 1;
    case 'off'
        tempValue = 0;
end
hand.checkID(5) = uicontrol('Parent',                hTabs(2),...
    'Tag',                   'showxaxisgrid',...
    'Style',                 'checkbox',...
    'Value',                 tempValue,...
    'String',                'Display X Axis Grid',...
    'Units',                 'normalized',...
    'Position',              [.07 .68 .55 .07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@showxaxisgrid,axIds});

hand.textID(9) = uicontrol(hTabs(2),...
    'Style',                 'text',...
    'String',                'Number of Ticks in X',...
    'Units',                 'normalized',...
    'Position',              [.07    0.58   .55    0.0450],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);

hand.editID(10) = uicontrol(hTabs(2),...
    'Tag',                   'xTickNumb',...
    'Style',                 'edit',...
    'String',                num2str(axPropert.XAxis.numXTicks),...
    'Units',                 'normalized',...
    'Position',              [0.67 .58 .2 .055],...
    'Callback',              {@setxticksnumber,axIds});

hand.buttonID(8) = uicontrol(hTabs(2),...
    'Style',                 'pushbutton',...
    'Tag',                   'xaxiscolorbutton',...
    'Value',                 0,...
    'String',                'Change X Axis Color...',...
    'Units',                 'normalized',...
    'Position',              [0.042 0.04 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       axPropert.XLabel.FontColor,...
    'ForegroundColor',       [1 1 1] - axPropert.XLabel.FontColor,...
    'Callback',              {@changefontxlabelcolor, axIds});

% Font Properties
hand.panelID(8) = uipanel(hTabs(2),...
    'Position',              [.02   0.02    0.9600    0.26],...
    'Title',                 'XLabel Font Properties',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

hand.buttonID(7) = uicontrol(hTabs(2),...
    'Style',                 'pushbutton',...
    'Tag',                   'xlabelfonttype',...
    'Value',                 0,...
    'String',                [axPropert.XLabel.FontName ', '  axPropert.XLabel.FontWeight ', '  axPropert.XLabel.FontAngle ' (' num2str(axPropert.XLabel.FontSize) ' points)'],...
    'Units',                 'normalized',...
    'Position',              [0.042 0.14 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@changefontxlabelprop, axIds});

hand.buttonID(8) = uicontrol(hTabs(2),...
    'Style',                 'pushbutton',...
    'Tag',                   'xlabelfontcolorbutton',...
    'Value',                 0,...
    'String',                'Change Font Color...',...
    'Units',                 'normalized',...
    'Position',              [0.042 0.04 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       [1 1 1] - axPropert.XLabel.FontColor,...
    'ForegroundColor',       axPropert.XLabel.FontColor,...
    'Callback',              {@changefontxlabelcolor, axIds});

%% ================ End of Creating Tab 2 (X Axis Tab) ================== %


%% ==================== Creating Tab 3 (Y Axis Tab) ===================== %
hTabs(3) = uitab('Parent',            hand.hTabGroup(1),...
    'Title',             'Y Axis   ',...
    'BackgroundColor',   figColor,...
    'ForegroundColor',   figColor);

hand.checkID(6) = uicontrol('Parent',                hTabs(3),...
    'Tag',                   'showylabelcheck',...
    'Style',                 'checkbox',...
    'Value',                 axPropert.YLabel.titleon,...
    'String',                'Display Y Axis Label',...
    'Units',                 'normalized',...
    'Position',              [.045 0.92    0.7    0.045],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@setnewaxisylabel,axIds});


hand.editID(11) = uicontrol('Parent',                hTabs(3),...
    'Tag',                   'showylabeledit',...
    'Style',                 'edit',...
    'Enable',                axPropert.YLabel.Enable,...
    'String',                axPropert.YLabel.String,...
    'Units',                 'normalized',...
    'Position',              [.045 .84 .91 .06],...
    'Callback',              {@introduceaxisylabel,axIds});

% Y Axis appearance
hand.panelID(8) = uipanel(hTabs(3),...
    'Position',              [.02   0.31    0.9600    0.5],...
    'Title',                 'Y Axis appearance',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

switch axPropert.YAxis.gridon
    case 'on'
        tempValue = 1;
    case 'off'
        tempValue = 0;
end
hand.checkID(5) = uicontrol('Parent',                hTabs(3),...
    'Tag',                   'showyaxisgrid',...
    'Style',                 'checkbox',...
    'Value',                 tempValue,...
    'String',                'Display Y Axis Grid',...
    'Units',                 'normalized',...
    'Position',              [.07 .68 .55 .07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@showyaxisgrid,axIds});



hand.textID(9) = uicontrol(hTabs(3),...
    'Style',                 'text',...
    'String',                'Number of Ticks in Y',...
    'Units',                 'normalized',...
    'Position',              [.07    0.58   .55    0.0450],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);

hand.editID(10) = uicontrol(hTabs(3),...
    'Tag',                   'yTickNumb',...
    'Style',                 'edit',...
    'String',                num2str(axPropert.YAxis.numXTicks),...
    'Units',                 'normalized',...
    'Position',              [0.67 .58 .2 .055],...
    'Callback',              {@setyticksnumber,axIds});

% Font Properties
hand.panelID(9) = uipanel(hTabs(3),...
    'Position',              [.02   0.02    0.9600    0.26],...
    'Title',                 'YLabel Font Properties',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

hand.buttonID(9) = uicontrol(hTabs(3),...
    'Style',                 'pushbutton',...
    'Tag',                   'ylabelfonttype',...
    'Value',                 0,...
    'String',                [axPropert.YLabel.FontName ', '  axPropert.YLabel.FontWeight ', '  axPropert.YLabel.FontAngle ' (' num2str(axPropert.YLabel.FontSize) ' points)'],...
    'Units',                 'normalized',...
    'Position',              [0.042 0.14 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@changefontylabelprop, axIds});

hand.buttonID(10) = uicontrol(hTabs(3),...
    'Style',                 'pushbutton',...
    'Tag',                   'ylabelfontcolorbutton',...
    'Value',                 0,...
    'String',                'Change Font Color...',...
    'Units',                 'normalized',...
    'Position',              [0.042 0.04 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       axPropert.YLabel.FontColor,...
    'ForegroundColor',       [1 1 1] - axPropert.YLabel.FontColor,...
    'Callback',              {@changefontylabelcolor, axIds});
%% ================== End of Creating Tab Y (Axis Tab) ================== %


%% ==================== Creating Tab 4 (Z Axis Tab) ===================== %
hTabs(4) = uitab('Parent',            hand.hTabGroup(1),...
    'Title',             'Z Axis   ',...
    'BackgroundColor',   figColor,...
    'ForegroundColor',   figColor);

hand.checkID(7) = uicontrol('Parent',                hTabs(4),...
    'Tag',                   'showzlabelcheck',...
    'Style',                 'checkbox',...
    'Value',                 axPropert.ZLabel.titleon,...
    'String',                'Display Z Axis Label',...
    'Units',                 'normalized',...
    'Position',              [.045 0.92    0.7    0.045],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@setnewaxiszlabel,axIds});


hand.editID(12) = uicontrol('Parent',                hTabs(4),...
    'Tag',                   'showzlabeledit',...
    'Style',                 'edit',...
    'Enable',                axPropert.ZLabel.Enable,...
    'String',                axPropert.ZLabel.String,...
    'Units',                 'normalized',...
    'Position',              [.045 .84 .91 .06],...
    'Callback',              {@introduceaxiszlabel,axIds});

% Y Axis appearance
hand.panelID(8) = uipanel(hTabs(4),...
    'Position',              [.02   0.31    0.9600    0.5],...
    'Title',                 'Z Axis appearance',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

switch axPropert.ZAxis.gridon
    case 'on'
        tempValue = 1;
    case 'off'
        tempValue = 0;
end
hand.checkID(5) = uicontrol('Parent',                hTabs(4),...
    'Tag',                   'showzaxisgrid',...
    'Style',                 'checkbox',...
    'Value',                 tempValue,...
    'String',                'Display Z Axis Grid',...
    'Units',                 'normalized',...
    'Position',              [.07 .68 .55 .07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@showzaxisgrid,axIds});

hand.textID(9) = uicontrol(hTabs(4),...
    'Style',                 'text',...
    'String',                'Number of Ticks in Z',...
    'Units',                 'normalized',...
    'Position',              [.07    0.58   .55    0.0450],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);

hand.editID(10) = uicontrol(hTabs(4),...
    'Tag',                   'zTickNumb',...
    'Style',                 'edit',...
    'String',                num2str(axPropert.ZAxis.numXTicks),...
    'Units',                 'normalized',...
    'Position',              [0.67 .58 .2 .055],...
    'Callback',              {@setzticksnumber,axIds});




% Font Properties
hand.panelID(10) = uipanel(hTabs(4),...
    'Position',              [.02   0.02    0.9600    0.26],...
    'Title',                 'ZLabel Font Properties',...
    'BackgroundColor',       figColor,...
    'HighlightColor',        figColor/4,...
    'ForegroundColor',       [1 1 1]);

hand.buttonID(11) = uicontrol(hTabs(4),...
    'Style',                 'pushbutton',...
    'Tag',                   'zlabelfonttype',...
    'Value',                 0,...
    'String',                [axPropert.ZLabel.FontName ', '  axPropert.ZLabel.FontWeight ', '  axPropert.ZLabel.FontAngle ' (' num2str(axPropert.ZLabel.FontSize) ' points)'],...
    'Units',                 'normalized',...
    'Position',              [0.042 0.14 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@changefontzlabelprop, axIds});

hand.buttonID(12) = uicontrol(hTabs(4),...
    'Style',                 'pushbutton',...
    'Tag',                   'zlabelfontcolorbutton',...
    'Value',                 0,...
    'String',                'Change Font Color...',...
    'Units',                 'normalized',...
    'Position',              [0.042 0.04 0.9 0.07],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       axPropert.ZLabel.FontColor,...
    'ForegroundColor',       [1 1 1] - axPropert.ZLabel.FontColor,...
    'Callback',              {@changefontzlabelcolor, axIds});
%% ================== End of Creating Tab Z (Axis Tab) ================== %



guidata(hand.FigpropID,hand);
movegui('center')
return;


%% =================== Setting Individual Axis limits =================== %
% 1. ------------------- Minimum Limits
% Minimum X limit
function setaxminlim(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    valtext = source.String;
    val = str2num(valtext);
    if ~isempty(val)
        if val < axIds(j).XLim(2)
            axIds(j).XLim(1) = val;
        end
    else
        warning(['Wrong minimun limit in the X axis( Axis handle: ' num2str(axIds(j)) '). The text could not be converted to number. Default limits are keeped']);
    end
end
return

% Minimum Y limit
function setayminlim(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    valtext = source.String;
    val = str2num(valtext);
    if ~isempty(val)
        if val < axIds(j).YLim(2)
            axIds(j).YLim(1) = val;
        end
    else
        warning(['Wrong minimun limit in the Y axis( Axis handle: ' num2str(axIds(j)) '). The text could not be converted to number. Default limits are keeped']);
    end
end
return

% Minimum Z limit
function setazminlim(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    valtext = source.String;
    val = str2num(valtext);
    if ~isempty(val)
        if val < axIds(j).ZLim(2)
            axIds(j).ZLim(1) = val;
        end
    else
        warning(['Wrong minimun limit in the Z axis( Axis handle: ' num2str(axIds(j)) '). The text could not be converted to number. Default limits are keeped']);
    end
end
return

% 2. ------------------- Maximum Limits
% Maximum Y limit
function setaxmaxlim(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
hSurfs = hSurfs(Naxis);

for j = 1:length(axIds)
    valtext = source.String;
    val = str2num(valtext);
    if ~isempty(val)
        if val > axIds(j).XLim(1)
            axIds(j).XLim(2) = val;
        end
    else
        warning(['Wrong maximun limit in the X axis( Axis handle: ' num2str(axIds(j)) '). The text could not be converted to number. Default limits are keeped']);
    end
end
return

% Maximum Y limit
function setaymaxlim(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    valtext = source.String;
    val = str2num(valtext);
    if ~isempty(val)
        if val > axIds(j).YLim(1)
            axIds(j).YLim(2) = val;
        end
    else
        warning(['Wrong maximun limit in the Y axis( Axis handle: ' num2str(axIds(j)) '). The text could not be converted to number. Default limits are keeped']);
    end
end
return

% Maximum Z limit
function setazmaxlim(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    valtext = source.String;
    val = str2num(valtext);
    if ~isempty(val)
        if val > axIds(j).ZLim(1)
            axIds(j).ZLim(2) = val;
        end
    else
        warning(['Wrong maximun limit in the Z axis( Axis handle: ' num2str(axIds(j)) '). The text could not be converted to number. Default limits are keeped']);
    end
end
return

%% =============== End of setting Individual Axis limits ================ %
% % function setmap(source,callbackdata, axIds, hSurfs);
function setmap(source,callbackdata, axIds, hSurfs, jSurfs, surfRealMaps);
val = source.Value;

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
hSurfs = hSurfs(Naxis);
jSurfs = jSurfs(Naxis);
surfRealMaps = surfRealMaps(Naxis);



if val > 1
    for j = 1:length(axIds)
        hSurf = hSurfs{j};
        jsurf = jSurfs{j};
        surfRealMap = surfRealMaps{j};
        
        maps = cellstr(strvcat('Custom color...' ,colormaps_colors));
        newmap = maps{val};
        clmap = newmap;
        
        
        if surfRealMap
            % For R2014a and earlier:
            % val = get(source,'Value');
            % maps = get(source,'String');
            bb = findobj(axIds(j),'Type','patch');
            Surf.Is = jsurf.Is(surfRealMap);
            [Colors] = Surf_Color(Surf,newmap);
            bb.FaceVertexCData = Colors;
            bb.FaceColor = 'interp';
            
            %             for i = 1:length(bb)
            %                 if ~isempty(bb(i).FaceVertexCData)
            %                     colors = bb(i).FaceVertexCData*255;
            %                     Surf.Is = colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16;
            %                     [Colors] = Surf_Color(Surf,newmap);
            %
            %                 end
            %             end
            
            cent_col=(0-min(Surf.Is))/(max(Surf.Is)-min(Surf.Is));
            colObject = findobj(hSurf,'Tag','Colorbar');
            if ~isempty(colObject)
                set(colObject,'Visible','on');
            end
            colormap(axIds(j),colormaps_colors(newmap,128,cent_col));
        else
            bb = findobj(axIds(j),'Type','patch');
            bb.FaceColor = 'interp';
        end
        axPropert = axIds(j).UserData;
        axPropert.Colormap = clmap;
        axIds(j).UserData = axPropert;
    end
else
    color = uisetcolor;
    clmap = 'Custom color...';
    for j = 1:length(axIds)
        hSurf = hSurfs{j};
        bb = findobj(axIds(j),'Type','patch');
        bb.FaceColor = color;
        colObject = findobj(hSurf,'Tag','Colorbar');
        if ~isempty(colObject)
            set(colObject,'Visible','off');
        end
        
        axPropert = axIds(j).UserData;
        axPropert.Colormap = clmap;
        axIds(j).UserData = axPropert;
    end
end
return;

function setazimuthangle(source,callbackdata,axIds);
global hlink hand;


tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

% Detecting Elevation Angle
editobj = findobj(source.Parent,'Tag','elevationedit');
elevvaltext = editobj.String;
elev = str2num(elevvaltext);
if isempty(elev)
    error('Wrong Elevation Angle value');
    return;
end

% Selecting AzimuthalAngle
valtext = source.String;
azith = str2num(valtext);
if isempty(azith)
    error('Wrong Azimuthal Angle value');
    return;
end
for j = 1:length(axIds)
    removetarget(hlink,axIds(j));
    view(axIds(j),[azith elev]);
end
checkobj = findobj(source.Parent,'Tag','unlinkpropercheck'); % Unlink axis
set(checkobj,'Value',1);
return;


function setelevationangle(source,callbackdata,axIds);
global hlink;

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

% Detecting Elevation Angle
editobj = findobj(source.Parent,'Tag','azimuthaledit');
azithvaltext = editobj.String;
azith = str2num(azithvaltext);
if isempty(azith)
    error('Wrong Azimuthal Angle value');
    return;
end

valtext = source.String;
elev = str2num(valtext);
if isempty(elev)
    error(['Wrong Elevation Angle value']);
    return;
end
for j = 1:length(axIds)
    removetarget(hlink,axIds(j));
    view(axIds(j),[azith elev]);
end
checkobj = findobj(source.Parent,'Tag','unlinkpropercheck'); % Unlink axis
set(checkobj,'Value',1);
return;

function resetcameraview(source,callbackdata, axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

azithobj = findobj(source.Parent,'Tag','azimuthaledit');
elevobj = findobj(source.Parent,'Tag','elevationedit');
for j = 1:length(axIds)
    [azith, elev] = view(axIds(j),3);
    azithobj.String = num2str(azith);
    elevobj.String = num2str(elev);
end
return;

function switchcolobar(source,callbackdata,axIds, hSurfs);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
hSurfs = hSurfs(Naxis);

value = source.Value;
for j = 1:length(axIds)
    hSurf = hSurfs{j};
    axPropert = axIds(j).UserData;
    colObject = findobj(hSurf,'Tag','Colorbar');
    if ~isempty(colObject)
        if value == 1
            set(colObject,'Visible','on');
        else
            set(colObject,'Visible','off');
        end
    end
    axPropert.colorbaron = value;
    axIds(j).UserData = axPropert;
end
return

function switchaxisgrid(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    if value == 1
        axIds(j).XColor = [1 1 1] - axIds(j).Color;
        axIds(j).XGrid = 'on';
        axIds(j).YColor = [1 1 1] - axIds(j).Color;
        axIds(j).YGrid = 'on';
        axIds(j).ZColor = [1 1 1] - axIds(j).Color;
        axIds(j).ZGrid = 'on';
        tempobj = findobj('Tag','showxaxisgrid');
        tempobj.Value = 1;
        tempobj = findobj('Tag','showyaxisgrid');
        tempobj.Value = 1;
        tempobj = findobj('Tag','showzaxisgrid');
        tempobj.Value = 1;
        
    else
        axIds(j).XColor = axIds(j).Color;
        axIds(j).XGrid = 'off';
        axIds(j).YColor = axIds(j).Color;
        axIds(j).YGrid = 'off';
        axIds(j).ZColor = axIds(j).Color;
        axIds(j).ZGrid = 'off';
        tempobj = findobj('Tag','showxaxisgrid');
        tempobj.Value = 0;
        tempobj = findobj('Tag','showxlabelcheck');
        tempobj.Value = 0;
        tempobj = findobj('Tag','showxlabeledit');
        tempobj.Enable = 'off'
        tempobj = findobj('Tag','showyaxisgrid');
        tempobj.Value = 0;
        tempobj = findobj('Tag','showylabelcheck');
        tempobj.Value = 0;
        tempobj = findobj('Tag','showylabeledit');
        tempobj.Enable = 'off'
        tempobj = findobj('Tag','showzaxisgrid');
        tempobj.Value = 0;
        tempobj = findobj('Tag','showzlabelcheck');
        tempobj.Value = 0;
        tempobj = findobj('Tag','showzlabeledit');
        tempobj.Enable = 'off'
    end
    axPropert.XAxis.gridon = axIds(j).XGrid;
    axPropert.XAxis.Color = axIds(j).XColor;
    axPropert.YAxis.gridon = axIds(j).YGrid;
    axPropert.YAxis.Color = axIds(j).YColor;
    axPropert.ZAxis.gridon = axIds(j).ZGrid;
    axPropert.ZAxis.Color = axIds(j).ZColor;
    axPropert.gridon = value;
    axIds(j).UserData = axPropert;
end

function unlinkaxis(source,callbackdata,axIds);
global hlink

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    if value == 1
        removetarget(hlink,axIds(j));
    else
        addtarget(hlink,axIds(j));
    end
    axPropert.unlinkon = value;
    axIds(j).UserData = axPropert;
end
return;

% Change Axis Title
function setnewaxistitle(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;
titleobj = findobj(source.Parent,'Tag','showtitleedit'); % Unlink axis
butcolor = findobj(source.Parent,'Tag','titlefontcolorbutton');

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axColor = axIds(j).Title.Color;
    if value == 1
        set(titleobj,'Enable','on');
        axIds(j).Title.String = titleobj.String;
        axIds(j).Title.Visible = 'on';
        if axIds(j).Title.Color == axIds(j).Color;
            axIds(j).Title.Color = [1 1 1] - axColor;
            butcolor.BackgroundColor =  [1 1 1] - axColor;
            butcolor.ForegroundColor = axColor;
        else
            butcolor.BackgroundColor =  axColor;
            butcolor.ForegroundColor = [1 1 1] - axColor;
        end
    else
        set(titleobj,'Enable','off');
        axIds(j).Title.Visible = 'off';
    end
    axPropert.Title.titleon = value;
    axPropert.Title.String = titleobj.String;
    axPropert.Title.Enable = axIds(j).Title.Visible;
    axIds(j).UserData = axPropert;
end


return

function introduceaxistitle(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).Title.String = source.String;
    axIds(j).Title.Visible = 'on';
    
    axPropert.Title.titleon = 1;
    axPropert.Title.String = source.String;
    axIds(j).UserData = axPropert;
end
return

% Change Axis XLABEL
function setnewaxisxlabel(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;
xlabelobj = findobj(source.Parent,'Tag','showxlabeledit'); % Unlink axis
butcolor = findobj(source.Parent,'Tag','xlabelfontcolorbutton');

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axColor = axIds(j).XLabel.Color;
    if value == 1
        set(xlabelobj,'Enable','on');
        axIds(j).XLabel.String = xlabelobj.String;
        axIds(j).XLabel.Visible = 'on';
        if axIds(j).XLabel.Color == axIds(j).Color;
            axIds(j).XLabel.Color = [1 1 1] - axColor;
            butcolor.BackgroundColor =  [1 1 1] - axColor;
            butcolor.ForegroundColor = axColor;
        else
            butcolor.BackgroundColor =  axColor;
            butcolor.ForegroundColor = [1 1 1] - axColor;
        end
        
    else
        set(xlabelobj,'Enable','off');
        axIds(j).XLabel.Visible = 'off';
    end
    axPropert.XLabel.titleon = value;
    axPropert.XLabel.String = xlabelobj.String;
    axPropert.XLabel.Enable = axIds(j).XLabel.Visible;
    axIds(j).UserData = axPropert;
end

return

function introduceaxisxlabel(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
% hSurfs = hSurfs(Naxis);

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).XLabel.String = source.String;
    axIds(j).XLabel.Visible = 'on';
    
    axPropert.XLabel.titleon = 1;
    axPropert.XLabel.String = source.String;
    axIds(j).UserData = axPropert;
end
return

% Change Axis YLABEL
function setnewaxisylabel(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;
ylabelobj = findobj(source.Parent,'Tag','showylabeledit'); % Unlink axis
butcolor = findobj(source.Parent,'Tag','ylabelfontcolorbutton');

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axColor = axIds(j).YLabel.Color;
    if value == 1
        set(ylabelobj,'Enable','on');
        axIds(j).YLabel.String = ylabelobj.String;
        axIds(j).YLabel.Visible = 'on';
        if axIds(j).YLabel.Color == axIds(j).Color;
            axIds(j).YLabel.Color = [1 1 1] - axColor;
            butcolor.BackgroundColor =  [1 1 1] - axColor;
            butcolor.ForegroundColor = axColor;
        else
            butcolor.BackgroundColor =  axColor;
            butcolor.ForegroundColor = [1 1 1] - axColor;
        end
    else
        set(ylabelobj,'Enable','off');
        axIds(j).YLabel.Visible = 'off';
    end
    axPropert.YLabel.titleon = value;
    axPropert.YLabel.String = ylabelobj.String;
    axPropert.YLabel.Enable = axIds(j).YLabel.Visible;
    axIds(j).UserData = axPropert;
end

return

function introduceaxisylabel(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).YLabel.String = source.String;
    axIds(j).YLabel.Visible = 'on';
    
    axPropert.YLabel.titleon = 1;
    axPropert.YLabel.String = source.String;
    axIds(j).UserData = axPropert;
end
return

% Change Axis ZLABEL
function setnewaxiszlabel(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;
zlabelobj = findobj(source.Parent,'Tag','showzlabeledit'); % Unlink axis
butcolor = findobj(source.Parent,'Tag','zlabelfontcolorbutton');

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axColor = axIds(j).ZLabel.Color;
    if value == 1
        set(zlabelobj,'Enable','on');
        axIds(j).ZLabel.String = zlabelobj.String;
        axIds(j).ZLabel.Visible = 'on';
        if axIds(j).ZLabel.Color == axIds(j).Color;
            axIds(j).ZLabel.Color = [1 1 1] - axColor;
            butcolor.BackgroundColor =  [1 1 1] - axColor;
            butcolor.ForegroundColor = axColor;
        else
            butcolor.BackgroundColor =  axColor;
            butcolor.ForegroundColor = [1 1 1] - axColor;
        end
    else
        set(zlabelobj,'Enable','off');
        axIds(j).ZLabel.Visible = 'off';
    end
    axPropert.ZLabel.titleon = value;
    axPropert.ZLabel.String = zlabelobj.String;
    axPropert.ZLabel.Enable = axIds(j).ZLabel.Visible;
    axIds(j).UserData = axPropert;
end

return

function introduceaxiszlabel(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).ZLabel.String = source.String;
    axIds(j).ZLabel.Visible = 'on';
    
    axPropert.ZLabel.titleon = 1;
    axPropert.ZLabel.String = source.String;
    axIds(j).UserData = axPropert;
end
return


function resetaxislimits(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

xMinobj = findobj(source.Parent,'Tag','xMinEdit');
xMaxobj = findobj(source.Parent,'Tag','xMaxEdit');

yMinobj = findobj(source.Parent,'Tag','yMinEdit');
yMaxobj = findobj(source.Parent,'Tag','yMaxEdit');

zMinobj = findobj(source.Parent,'Tag','zMinEdit');
zMaxobj = findobj(source.Parent,'Tag','zMaxEdit');

for j = 1:length(axIds)
    axis(axIds(j),'image');
    xMinobj.String = num2str(axIds(j).XLim(1));
    xMaxobj.String = num2str(axIds(j).XLim(2));
    yMinobj.String = num2str(axIds(j).YLim(1));
    yMaxobj.String = num2str(axIds(j).YLim(2));
    zMinobj.String = num2str(axIds(j).ZLim(1));
    zMaxobj.String = num2str(axIds(j).ZLim(2));
end
return;


% Change Title Font properties
function changefonttitleprop(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontprop = uisetfont('Selecting New Font for Axis Title');
source.String = [fontprop.FontName ', '  fontprop.FontWeight ', '  fontprop.FontAngle  '(' num2str(fontprop.FontSize) ' points)'];
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).Title.FontName = fontprop.FontName;
    axIds(j).Title.FontWeight = fontprop.FontWeight;
    axIds(j).Title.FontSize = fontprop.FontSize;
    axIds(j).Title.FontAngle = fontprop.FontAngle;
    
    axPropert.Title.FontName = fontprop.FontName;
    axPropert.Title.FontSize = fontprop.FontSize;
    axPropert.Title.FontWeight = fontprop.FontWeight;
    axPropert.Title.FontAngle = fontprop.FontAngle;
    axIds(j).UserData = axPropert;
end
return;
% Change Title Font Color
function changefonttitlecolor(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontcolor = uisetcolor('Selecting New Font Color');
source.BackgroundColor = fontcolor;
source.ForegroundColor = [1 1 1] - fontcolor;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).Title.Color = fontcolor;
    axPropert.Title.FontColor = fontcolor;
    axIds(j).UserData = axPropert;
end
return;

% Change XLabel Font properties
function changefontxlabelprop(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontprop = uisetfont('Selecting New Font for Axis Title');
source.String = [fontprop.FontName ', '  fontprop.FontWeight ', '  fontprop.FontAngle  '(' num2str(fontprop.FontSize) ' points)'];
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).XLabel.FontName = fontprop.FontName;
    axIds(j).XLabel.FontWeight = fontprop.FontWeight;
    axIds(j).XLabel.FontSize = fontprop.FontSize;
    axIds(j).XLabel.FontAngle = fontprop.FontAngle;
    
    axPropert.XLabel.FontName = fontprop.FontName;
    axPropert.XLabel.FontSize = fontprop.FontSize;
    axPropert.XLabel.FontWeight = fontprop.FontWeight;
    axPropert.XLabel.FontAngle = fontprop.FontAngle;
    axIds(j).UserData = axPropert;
end
return;
% Change XLabel Font Color
function changefontxlabelcolor(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontcolor = uisetcolor('Selecting New Font Color');
source.BackgroundColor = fontcolor;
source.ForegroundColor = [1 1 1] - fontcolor;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).XLabel.Color = fontcolor;
    axPropert.XLabel.FontColor = fontcolor;
    axIds(j).UserData = axPropert;
end
return;

% Change YLabel Font properties
function changefontylabelprop(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontprop = uisetfont('Selecting New Font for Axis Title');
source.String = [fontprop.FontName ', '  fontprop.FontWeight ', '  fontprop.FontAngle  '(' num2str(fontprop.FontSize) ' points)'];
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).YLabel.FontName = fontprop.FontName;
    axIds(j).YLabel.FontWeight = fontprop.FontWeight;
    axIds(j).YLabel.FontSize = fontprop.FontSize;
    axIds(j).YLabel.FontAngle = fontprop.FontAngle;
    
    axPropert.YLabel.FontName = fontprop.FontName;
    axPropert.YLabel.FontSize = fontprop.FontSize;
    axPropert.YLabel.FontWeight = fontprop.FontWeight;
    axPropert.YLabel.FontAngle = fontprop.FontAngle;
    axIds(j).UserData = axPropert;
end
return;

% Change YLabel Font Color
function changefontylabelcolor(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontcolor = uisetcolor('Selecting New Font Color');
source.BackgroundColor = fontcolor;
source.ForegroundColor = [1 1 1] - fontcolor;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).YLabel.Color = fontcolor;
    axPropert.YLabel.FontColor = fontcolor;
    axIds(j).UserData = axPropert;
end
return;

% Change ZLabel Font Properties
function changefontzlabelprop(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontprop = uisetfont('Selecting New Font for Axis Title');
source.String = [fontprop.FontName ', '  fontprop.FontWeight ', '  fontprop.FontAngle  '(' num2str(fontprop.FontSize) ' points)'];
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).ZLabel.FontName = fontprop.FontName;
    axIds(j).ZLabel.FontWeight = fontprop.FontWeight;
    axIds(j).ZLabel.FontSize = fontprop.FontSize;
    axIds(j).ZLabel.FontAngle = fontprop.FontAngle;
    
    axPropert.ZLabel.FontName = fontprop.FontName;
    axPropert.ZLabel.FontSize = fontprop.FontSize;
    axPropert.ZLabel.FontWeight = fontprop.FontWeight;
    axPropert.ZLabel.FontAngle = fontprop.FontAngle;
    axIds(j).UserData = axPropert;
end
return;

% Change ZLabel Font Color
function changefontzlabelcolor(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

fontcolor = uisetcolor('Selecting New Font Color');
source.BackgroundColor = fontcolor;
source.ForegroundColor = [1 1 1] - fontcolor;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    axIds(j).ZLabel.Color = fontcolor;
    axPropert.ZLabel.FontColor = fontcolor;
    axIds(j).UserData = axPropert;
end
return;


function showxaxisgrid(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

value = source.Value;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    if value == 1
        axIds(j).XColor = [1 1 1] - axIds(j).Color;
        axIds(j).XGrid = 'on';
        axPropert.XAxis.Color = axIds(j).XColor;
        axPropert.XAxis.gridon = 'on';
    else
        axIds(j).XColor = axIds(j).Color;
        axPropert.XAxis.gridon = 'off';
        axPropert.XAxis.Color = axIds(j).XColor;
        axPropert.XAxis.gridon = 'off';
    end
    axIds(j).UserData = axPropert;
end
return;

function showyaxisgrid(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
value = source.Value;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    if value == 1
        axIds(j).XColor = [1 1 1] - axIds(j).Color;
        axIds(j).YGrid = 'on';
        axPropert.YAxis.Color = axIds(j).XColor;
        axPropert.YAxis.gridon = 'on';
    else
        axIds(j).YColor = axIds(j).Color;
        axPropert.YAxis.gridon = 'off';
        axPropert.YAxis.Color = axIds(j).XColor;
        axPropert.YAxis.gridon = 'off';
    end
    axIds(j).UserData = axPropert;
end
return;

function showzaxisgrid(source,callbackdata,axIds);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
value = source.Value;

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    if value == 1
        axIds(j).ZColor = [1 1 1] - axIds(j).Color;
        axIds(j).ZGrid = 'on';
        axPropert.ZAxis.Color = axIds(j).XColor;
        axPropert.ZAxis.gridon = 'on';
    else
        axIds(j).ZColor = axIds(j).Color;
        axPropert.ZAxis.gridon = 'off';
        axPropert.ZAxis.Color = axIds(j).XColor;
        axPropert.ZAxis.gridon = 'off';
    end
    axIds(j).UserData = axPropert;
end
return;

function  setxticksnumber(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);

numTicks = round(str2num(source.String));
for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    tickValues = axIds(j).XTick;
    minVal  = min(tickValues);
    maxVal  = max(tickValues);
    range = maxVal-minVal;
    if numTicks<=1
        values = [minVal maxVal];
    else
        values =  minVal:range/(numTicks-1):maxVal;
    end
    cbar.Limits = [min(values) max(values)];
    axIds(j).XTick = values(:);
    axIds(j).XTickLabel = cellstr(num2str(values(:)));
    axPropert.XAxis.numXTicks = numTicks;
    axIds(j).UserData = axPropert;
end
return

function  setyticksnumber(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
numTicks = round(str2num(source.String));

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    tickValues = axIds(j).YTick;
    minVal  = min(tickValues);
    maxVal  = max(tickValues);
    range = maxVal-minVal;
    if numTicks<=1
        values = [minVal maxVal];
    else
        values =  minVal:range/(numTicks-1):maxVal;
    end
    cbar.Limits = [min(values) max(values)];
    axIds(j).YTick = values(:);
    axIds(j).YTickLabel = cellstr(num2str(values(:)));
    axPropert.YAxis.numXTicks = numTicks;
    axIds(j).UserData = axPropert;
end
return

function  setzticksnumber(source,callbackdata,axIds);
tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
numTicks = round(str2num(source.String));

for j = 1:length(axIds)
    axPropert = axIds(j).UserData;
    tickValues = axIds(j).ZTick;
    minVal  = min(tickValues);
    maxVal  = max(tickValues);
    range = maxVal-minVal;
    if numTicks<=1
        values = [minVal maxVal];
    else
        values =  minVal:range/(numTicks-1):maxVal;
    end
    cbar.Limits = [min(values) max(values)];
    axIds(j).ZTick = values(:);
    axIds(j).ZTickLabel = cellstr(num2str(values(:)));
    axPropert.YAxis.numXTicks = numTicks;
    axIds(j).UserData = axPropert;
end
return


function changecolorbar(source,callbackdata, axIds, hSurfs);

tempobj = findobj('Tag','axispopup');
if length(tempobj.String) > 1
    if tempobj.Value == length(tempobj.String)
        Naxis = [1:length(tempobj.String)-1];
    else
        Naxis = tempobj.Value;
    end
else
    Naxis = 1;
end
axIds = axIds(Naxis);
hSurfs = hSurfs(Naxis);

figUnits = 'centimeters';
figColor = [0.3 0.3 0.3];
set(0,'units',figUnits);
cm_screen = get(0,'screensize');

scaleFactorWidth = 7;                                                     % Figure Width scale factor
scaleFactorHeight = 5;                                                    % Figure Height scale factor

figWidth  = (cm_screen(3)-cm_screen(1))/scaleFactorWidth;                 % Figure width
figHeight = (cm_screen(4)-3*cm_screen(3)/cm_screen(4))/scaleFactorHeight; % Figure height
figCornerX = cm_screen(3)/2-(figWidth/2);                                 % Figure Corner (X coordinate)
figCornerY = cm_screen(4)/2-(figHeight/2);                                % Figure Corner (Y coordinate)
figPosition = [figCornerX figCornerY figWidth figHeight];                 % Figure Position
% if ~isempty(axPropert.ColorbarFig)&&ishandle(axPropert.ColorbarFig) % Detecting and closing old figures related with this axis
%     close(axPropert.ColorbarFig);
% end
% Creating figure


cont = 0;
colorBars{1} = '';
for j = 1:length(axIds)
    hSurf = hSurfs{j};
    colObject = findobj(hSurf,'Tag','Colorbar');
    if ~isempty(colObject)
        cont = cont + 1;
        colorBars{cont} = colObject;
    end
end

if ~isempty(colorBars)
    tete = colorBars{1};
else
end


axPropert.ColorbarFig = figure('numbertitle',    'off',...
    'name',           'Change colorbar properties ...',...
    'Color',          figColor,...
    'units',          figUnits,...
    'Position',       figPosition,...
    'InvertHardcopy', 'off',...
    'MenuBar',        'None',...
    'Resize',         'off');

% Creating Left Panel Properties
hand.panelID(1) = uipanel(axPropert.ColorbarFig,...
    'Position',            [.02 scaleFactorWidth/scaleFactorHeight*.02 .96 .96],...
    'Units',               'normalized',...
    'Title',               'Colorbar properties',...
    'BackgroundColor',     figColor,...
    'HighlightColor',      figColor/4,...
    'ForegroundColor',     [1 1 1]);

cbar = colorBars{1};
colorbarProper = cbar.UserData;
hand.checkID(7) = uicontrol('Parent',                hand.panelID(1),...
    'Tag',                   'showcolorbartitle',...
    'Style',                 'checkbox',...
    'Value',                 colorbarProper.Title.titleon,...
    'String',                colorbarProper.Title.String,...
    'Units',                 'normalized',...
    'Position',              [.045 0.87    0.7    0.09],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@showcolorbartitle,colorBars});


hand.editID(12) = uicontrol('Parent',                hand.panelID(1),...
    'Tag',                   'introduceacolorbartitle',...
    'Style',                 'edit',...
    'Enable',                colorbarProper.Title.Enable,...
    'String',                colorbarProper.Title.String,...
    'Units',                 'normalized',...
    'Position',              [.045 .72 .91 .11],...
    'Callback',              {@introduceacolorbartitle,colorBars});



hand.textID(9) = uicontrol(hand.panelID(1),...
    'Style',                 'text',...
    'String',                'Number of Colorbar Ticks',...
    'Units',                 'normalized',...
    'Position',              [.02    0.537   .7    0.09],...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);


hand.editID(10) = uicontrol(hand.panelID(1),...
    'Tag',                   'ColorbTickNumb',...
    'Style',                 'edit',...
    'String',                num2str(colorbarProper.Ticks.number),...
    'Units',                 'normalized',...
    'Position',              [0.7 .53 .2 .11],...
    'Callback',              {@setcolorbticknumber,colorBars});

% Font Properties

hand.buttonID(11) = uicontrol(hand.panelID(1),...
    'Style',                 'pushbutton',...
    'Value',                 0,...
    'String',               [colorbarProper.FontName ', '  colorbarProper.FontWeight ', '  colorbarProper.FontAngle ' (' num2str(colorbarProper.FontSize) ' points)'],...
    'Units',                 'normalized',...
    'Position',              [0.042 0.3 0.91 0.15],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0],...
    'Callback',              {@changecolorbarfont, colorBars});
% 'String',                ,...


hand.buttonID(12) = uicontrol(hand.panelID(1),...
    'Style',                 'pushbutton',...
    'Tag',                   'zlabelfontcolorbutton',...
    'Value',                 0,...
    'String',                'Change Font Color...',...
    'Units',                 'normalized',...
    'Position',              [0.042 0.07 0.91 0.15],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       colorbarProper.FontColor,...
    'ForegroundColor',       [1 1 1] - colorbarProper.FontColor,...
    'Callback',              {@changecolorbarfontcolor, colorBars});
%                             'BackgroundColor',       axPropert.ZLabel.FontColor,...
%                             'ForegroundColor',       [1 1 1] - axPropert.ZLabel.FontColor,...


%% %%%%%%%%%%%%%%%%% Defining Left Panel Objects %%%%%%%%%%%%%%%%%%%%%%%% %


return;


function showcolorbartitle(source,callbackdata, colorBars)
cbarFIG = source.Parent;
editobj = findobj(source.Parent,'Tag','introduceacolorbartitle');
if ~isempty(colorBars)
    origFIG = get(colorBars{1},'Parent');
    for j = 1:length(colorBars)
        cbar = colorBars{j};
        colorbarProper = cbar.UserData;
        tempHand = title(cbar, source.String);
        if  source.Value
            colorbarProper.Title.titleon = 1;
            if tempHand.Color == origFIG.Color;
                tempHand.Color = [1 1 1] - origFIG.Color;
            end
            set(tempHand,'Visible','on');
            set(editobj,'Enable','on');
        else
            colorbarProper.Title.titleon = 0;
            set(editobj,'Enable','off');
            set(tempHand,'Visible','off');
        end
        colorbarProper.Title.String = source.String;
        colorbarProper.Title.Color = tempHand.Color;
        colorbarProper.Title.Enable = editobj.Enable;
        cbar.UserData = colorbarProper;
    end
end
return;

function introduceacolorbartitle(source,callbackdata, colorBars);
if ~isempty(colorBars)
    origFIG = get(colorBars{1},'Parent');
    for j = 1:length(colorBars)
        cbar = colorBars{j};
        colorbarProper = cbar.UserData;
        tempHand = title(cbar, source.String,'Tag','colorbartitle');
        if tempHand.Color == origFIG.Color;
            tempHand.Color = [1 1 1] - origFIG.Color;
        end
        colorbarProper.Title.Handle = tempHand;
        colorbarProper.Title.String = source.String;
        colorbarProper.Title.Color = tempHand.Color;
        cbar.UserData = colorbarProper;
    end
end
return

function setcolorbticknumber(source,callbackdata, colorBars);

for j = 1:length(colorBars)
    cbar = colorBars{j};
    colorbarProper = cbar.UserData;
    tickValues = str2num(char(cbar.TickLabels));
    
    minVal  = min(tickValues);
    maxVal  = max(tickValues);
    range = maxVal-minVal;
    numTicks = round(str2num(source.String));
    if numTicks<=1
        values = [minVal maxVal];
    else
        values =  minVal:range/(numTicks-1):maxVal;
    end
    cbar.Limits = [min(values) max(values)];
    cbar.Ticks = linspace(cbar.Limits(1),cbar.Limits(2),length(values));
    cbar.TickLabels = sprintf('%.2f\n',values(:));
    colorbarProper.Ticks.number = numTicks;
    cbar.UserData = colorbarProper;
end
return

function changecolorbarfont(source,callbackdata, colorBars);
fontprop = uisetfont('Selecting New Font for Colorbar ');
source.String = [fontprop.FontName ', '  fontprop.FontWeight ', '  fontprop.FontAngle  '(' num2str(fontprop.FontSize) ' points)'];
for j = 1:length(colorBars)
    cbar = colorBars{j};
    colorbarProper = cbar.UserData;
    cbar.FontName = fontprop.FontName;
    cbar.FontWeight = fontprop.FontWeight;
    cbar.FontSize = fontprop.FontSize;
    cbar.FontAngle = fontprop.FontAngle;
    
    
    colorbarProper.FontName = fontprop.FontName;
    colorbarProper.FontSize = fontprop.FontSize;
    colorbarProper.FontWeight = fontprop.FontWeight;
    colorbarProper.FontAngle = fontprop.FontAngle;
    cbar.UserData = colorbarProper;
end
return;

function changecolorbarfontcolor(source,callbackdata, colorBars);
fontcolor = uisetcolor('Selecting New Font Color');
source.BackgroundColor = fontcolor;
source.ForegroundColor = [1 1 1] - fontcolor;

for j = 1:length(colorBars)
    cbar = colorBars{j};
    colorbarProper = cbar.UserData;
    cbar.Color = fontcolor;
    colorbarProper.Title.Handle.Color = fontcolor;
    colorbarProper.FontColor = fontcolor;
    colorbarProper.Title.Color = fontcolor;
    cbar.UserData = colorbarProper;
end
return;


function normalizecolormaps(source,callbackdata, axIds, hSurfs, jSurfs, surfRealMaps);

tempobj = findobj('Tag','axispopup');

if tempobj.Value == length(tempobj.String)& tempobj.Value ~= 1
    val = source.Value; % Checkbox Value
    popobj = findobj('Tag','colormappopup');
    popval = popobj.Value;
    source.Enable = 'on';
    maps = cellstr(strvcat('Custom color...' ,colormaps_colors));
    newmap = maps{popval};
    clmap = newmap;
    if val
        allIsReal = '';
        if popval > 1
            cont = 0;
            for j = 1:length(axIds)
                jsurf = jSurfs{j};
                surfRealMap = surfRealMaps{j};
                
                if surfRealMap
                    Surf.Is = jsurf.Is(surfRealMap);
                    cont = cont + 1;
                    if cont == 1
                        allIsReal = Surf.Is;
                        tempVar(cont,:)= [j 1 length(Surf.Is)];
                    else
                        tempVar(cont,:)= [j length(allIsReal)+1 length(allIsReal)+length(Surf.Is)];
                        allIsReal = [allIsReal;Surf.Is];
                    end
                end
            end
            if ~isempty(allIsReal)
                Surf.Is = allIsReal;
                
                [Colors] = Surf_Color(Surf,newmap);
                minVal = min(allIsReal);
                maxVal = max(allIsReal);
                cent_col=(0-minVal)/(maxVal-minVal);
                
                for j = 1:cont
                    hSurf = hSurfs{tempVar(j,1)};
                    jsurf = jSurfs{tempVar(j,1)};
                    surfRealMap = surfRealMaps{tempVar(j,1)};
                    
                    
                    if surfRealMap
                        bb = findobj(axIds(tempVar(j,1)),'Type','patch');
                        bb.FaceVertexCData(surfRealMap,:) = Colors(tempVar(j,2):tempVar(j,3),:);
                        bb.FaceColor = 'interp';
                        
                        
                        colObject = findobj(hSurf,'Tag','Colorbar');
                        if ~isempty(colObject)
                            set(colObject,'Visible','on','Limits',[ minVal maxVal]);
                        end
                        set(axIds(tempVar(j,1)),'CLim',[minVal,maxVal]);
                        colormap(axIds(tempVar(j,1)),colormaps_colors(newmap,128,cent_col));
                    end
                    
                    axPropert = axIds(tempVar(j,1)).UserData;
                    axPropert.colnormalizeon = 1;
                    axPropert.Colormap = clmap;
                    axIds(tempVar(j,1)).UserData = axPropert;
                    
                end
            end
        else
            color = uisetcolor;
            for j = 1:length(axIds)
                bb = findobj(axIds(j),'Type','patch');
                bb.FaceColor = color;
                colObject = findobj(hSurf,'Tag','Colorbar');
                if ~isempty(colObject)
                    set(colObject,'Visible','off');
                end
            end
            clmap = 'Custom color...';
            axPropert = axIds(j).UserData;
            axPropert.colnormalizeon = 1;
            axPropert.Colormap = clmap;
            axIds(j).UserData = axPropert;
        end
    else
        for j = 1:length(axIds)
            hSurf = hSurfs{j};
            jsurf = jSurfs{j};
            surfRealMap = surfRealMaps{j};
            maps = cellstr(strvcat('Custom color...' ,colormaps_colors));
            newmap = maps{popval};
            clmap = newmap;
            if popval > 1
                
                if surfRealMap
                    % For R2014a and earlier:
                    % val = get(source,'Value');
                    % maps = get(source,'String');
                    bb = findobj(axIds(j),'Type','patch');
                    Surf.Is = jsurf.Is(surfRealMap);
                    [Colors] = Surf_Color(Surf,newmap);
                    bb.FaceVertexCData = Colors;
                    bb.FaceColor = 'interp';
                    
                    minVal = min(Surf.Is);
                    maxVal = max(Surf.Is);
                    cent_col=(0-minVal)/(maxVal-minVal);
                    colObject = findobj(hSurf,'Tag','Colorbar');
                    if ~isempty(colObject)
                        set(colObject,'Visible','on');
                        set(colObject,'Visible','on','Limits',[ minVal maxVal]);
                    end
                    set(axIds(j),'CLim',[minVal,maxVal]);
                    colormap(axIds(j),colormaps_colors(newmap,128,cent_col));
                end
            else
                color = uisetcolor;
                bb = findobj(axIds(j),'Type','patch');
                bb.FaceColor = color;
                colObject = findobj(hSurf,'Tag','Colorbar');
                if ~isempty(colObject)
                    set(colObject,'Visible','off');
                end
                clmap = 'Custom color...';
            end
            axPropert = axIds(j).UserData;
            axPropert.colnormalizeon = 0;
            axPropert.Colormap = clmap;
            axIds(j).UserData = axPropert;
        end
    end
else
    source.Enable = 'off';
end
return;

function outVal = getaxis(source,callbackdata,axIds,hSurfs,jSurfs,surfRealMaps);
outVal = source.Value;

if length(source.String) > 1
    if outVal == length(source.String)
        Naxis = [1:length(source.String)-1];
    else
        Naxis = outVal;
    end
else
    Naxis = 1;
end


for j = 1:length(Naxis);
    axes(axIds(Naxis(j)));
    
    axPropert = axIds(Naxis(j)).UserData;
    
    
    % Creating colormaps list
    defColormaps = strvcat('Custom color...' ,colormaps_colors);
    
    % Detecting posible colormaps
    indposclmap = find(ismember(defColormaps(:,1:length(axPropert.Colormap)),axPropert.Colormap,'rows'));
    
    tempobj = findobj('Tag','colormappopup');
    tempobj.Value = indposclmap;
    objParent = tempobj.Parent;
    figure(objParent.Parent);
    
    %% ==================== End of Defining Colormap ======================== %
    
    %% ======================== Defining Views  ============================= %
    
    % Obtain current axis view
    [azAngle,elevAngle] = view(axIds(Naxis(j)));
    
    
    tempobj = findobj('Tag','azimuthaledit');
    tempobj.String = num2str(azAngle);
    
    tempobj = findobj('Tag','elevationedit');
    tempobj.String = num2str(elevAngle);
    
    % %
    %% ==================== End of Defining Views =========================== %
    % %
    % %
    %% =================== Colorbar properties checkboxes =================== %
    
    
    
    tempobj = findobj('Tag','showcolorbarcheck');
    tempobj.Value = axPropert.colorbaron;
    
    %% ===================== End of Colorbar properties ===================== %
    
    %% ======================== Some checkboxes ============================= %
    
    tempobj = findobj('Tag','showaxisgridcheck');
    tempobj.Value = axPropert.gridon;
    
    tempobj = findobj('Tag','unlinkpropercheck');
    tempobj.Value = axPropert.unlinkon;
    
    
    tempobj = findobj('Tag','colnormalizecheck');
    tempobj.Value = axPropert.colnormalizeon;
    
    if length(Naxis) > 1
        axPropert.colnormalizevisible = 'on';
    else
        axPropert.colnormalizevisible = 'off';
    end
     tempobj.Visible = axPropert.colnormalizevisible;
    
    %% ==================== End of Some checkboxes ============================= %
    
    
    %% %%%%%%%%%%%%%%%%% End of Left Panel Objects %%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    %% %%%%%%%%%%%%%%%%% Defining Right Panel Objects %%%%%%%%%%%%%%%%%%%%%%% %
    
    %% ==================== Creating Tab 1 (Title Tab) ====================== %
    tempobj = findobj('Tag','showtitleedit');
    tempobj.Enable = axPropert.Title.Enable;
    tempobj.String = axPropert.Title.String;
    tempobj = findobj('Tag','showtitlecheck');
    tempobj.Value = axPropert.Title.titleon;
    
    minXlim = axIds(Naxis(j)).XLim(1);
    maxXlim = axIds(Naxis(j)).XLim(2);
    tempobj = findobj('Tag','xMinEdit');
    tempobj.String = num2str(minXlim);
    tempobj = findobj('Tag','xMaxEdit');
    tempobj.String = num2str(maxXlim);
    
    minYlim = axIds(Naxis(j)).YLim(1);
    maxYlim = axIds(Naxis(j)).YLim(2);
    tempobj = findobj('Tag','yMinEdit');
    tempobj.String = num2str(minYlim);
    tempobj = findobj('Tag','yMaxEdit');
    tempobj.String = num2str(maxYlim);
    
    minZlim = axIds(Naxis(j)).ZLim(1);
    maxZlim = axIds(Naxis(j)).ZLim(2);
    tempobj = findobj('Tag','zMinEdit');
    tempobj.String = num2str(minZlim);
    tempobj = findobj('Tag','zMaxEdit');
    tempobj.String = num2str(maxZlim);
    
    
    
    tempobj = findobj('Tag','titlefonttype');
    tempobj.String = [axPropert.Title.FontName ', '  axPropert.Title.FontWeight ', '  axPropert.Title.FontAngle ' (' num2str(axPropert.Title.FontSize) ' points)'];
    
    tempobj = findobj('Tag','titlefontcolorbutton');
    tempobj.BackgroundColor = axPropert.Title.FontColor;
    tempobj.ForegroundColor = [1 1 1] - axPropert.Title.FontColor;
    %% ================== End of Creating Tab 1 (Title Tab) ================= %
    
    
    %% ==================== Creating Tab 2 (X Axis Tab) ===================== %
    tempobj = findobj('Tag','showxlabelcheck');
    tempobj.Value = axPropert.XLabel.titleon;
    
    
    tempobj = findobj('Tag','showxlabeledit');
    tempobj.Enable = axPropert.XLabel.Enable;
    tempobj.String = axPropert.XLabel.String;
    
    switch axPropert.XAxis.gridon
        case 'on'
            tempValue = 1;
        case 'off'
            tempValue = 0;
    end
    tempobj = findobj('Tag','showxaxisgrid');
    tempobj.Value = tempValue;
    
    tempobj = findobj('Tag','xTickNumb');
    tempobj.String =  num2str(axPropert.YAxis.numXTicks);
    
    tempobj = findobj('Tag','xlabelfonttype');
    tempobj.String = [axPropert.XLabel.FontName ', '  axPropert.XLabel.FontWeight ', '  axPropert.XLabel.FontAngle ' (' num2str(axPropert.XLabel.FontSize) ' points)'];
    
    tempobj = findobj('Tag','xlabelfontcolorbutton');
    tempobj.BackgroundColor = axPropert.XLabel.FontColor;
    tempobj.ForegroundColor = [1 1 1] - axPropert.XLabel.FontColor;
    %% ================== End of Creating Tab 2 (X Axis Tab) ================ %
    
    
    
    %% ==================== Creating Tab 3 (Y Axis Tab) ===================== %
    tempobj = findobj('Tag','showylabelcheck');
    tempobj.Value = axPropert.YLabel.titleon;
    
    
    tempobj = findobj('Tag','showylabeledit');
    tempobj.Enable = axPropert.YLabel.Enable;
    tempobj.String = axPropert.YLabel.String;
    
    switch axPropert.YAxis.gridon
        case 'on'
            tempValue = 1;
        case 'off'
            tempValue = 0;
    end
    tempobj = findobj('Tag','showyaxisgrid');
    tempobj.Value = tempValue;
    
    tempobj = findobj('Tag','yTickNumb');
    tempobj.String =  num2str(axPropert.YAxis.numXTicks);
    
    tempobj = findobj('Tag','ylabelfonttype');
    tempobj.String = [axPropert.YLabel.FontName ', '  axPropert.YLabel.FontWeight ', '  axPropert.YLabel.FontAngle ' (' num2str(axPropert.YLabel.FontSize) ' points)'];
    
    tempobj = findobj('Tag','ylabelfontcolorbutton');
    tempobj.BackgroundColor = axPropert.YLabel.FontColor;
    tempobj.ForegroundColor = [1 1 1] - axPropert.YLabel.FontColor;
    
    %% ================== End of Creating Tab 3 (Axis Tab) ================== %
    
    %% ==================== Creating Tab 4 (Z Axis Tab) ===================== %
    tempobj = findobj('Tag','showzlabelcheck');
    tempobj.Value = axPropert.ZLabel.titleon;
    
    
    tempobj = findobj('Tag','showzlabeledit');
    tempobj.Enable = axPropert.ZLabel.Enable;
    tempobj.String = axPropert.ZLabel.String;
    
    switch axPropert.ZAxis.gridon
        case 'on'
            tempValue = 1;
        case 'off'
            tempValue = 0;
    end
    tempobj = findobj('Tag','showzaxisgrid');
    tempobj.Value = tempValue;
    
    tempobj = findobj('Tag','zTickNumb');
    tempobj.String =  num2str(axPropert.ZAxis.numXTicks);
    
    tempobj = findobj('Tag','zlabelfonttype');
    tempobj.String = [axPropert.ZLabel.FontName ', '  axPropert.ZLabel.FontWeight ', '  axPropert.ZLabel.FontAngle ' (' num2str(axPropert.ZLabel.FontSize) ' points)'];
    
    tempobj = findobj('Tag','zlabelfontcolorbutton');
    tempobj.BackgroundColor = axPropert.ZLabel.FontColor;
    tempobj.ForegroundColor = [1 1 1] - axPropert.ZLabel.FontColor;
    
    %% ================ End of Creating Tab 4 (Z Axis Tab) ================== %
    
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %
    % % % % % %         end
    % % % % % %         source.UserData = outVal;
    % % % % % %
    % % % % % %
    % % % % % %         axIds(i) = gca;
    % % % % % %         axPropert.PropFig = '';
    % % % % % %         axPropert.ColorbarFig = '';
    % % % % % %         axPropert.Colormap = clmap;
    % % % % % %         axPropert.colorbaron = 1;
    % % % % % %         axPropert.gridon = 1;
    % % % % % %         axPropert.unlinkon = 0;
    % % % % % %         axPropert.colnormalizeon = 0;
    % % % % % %         axPropert.Limits = [axIds(i).XLim(1) axIds(i).XLim(2)...
    % % % % % %             axIds(i).YLim(1) axIds(i).YLim(2)...
    % % % % % %             axIds(i).ZLim(1) axIds(i).ZLim(2)];
    % % % % % %
    % % % % % %         %% Title Properties
    % % % % % %         axIds(i).Title.String = ['Axis ' num2str(i)];
    % % % % % %     axIds(i).Title.Color = [1 1 1] - param.figcolor;
    % % % % % %     axPropert.Title.titleon = 0;
    % % % % % %     axPropert.Title.String = ['Introduce New Title for Axis ' num2str(i)] ;
    % % % % % %     axPropert.Title.Enable = 'off';
    % % % % % %     % Font Properties Title
    % % % % % %     axPropert.Title.FontName = axIds(i).Title.FontName;
    % % % % % %     axPropert.Title.FontSize = axIds(i).Title.FontSize;
    % % % % % %     axPropert.Title.FontWeight = axIds(i).Title.FontWeight;
    % % % % % %     axPropert.Title.FontColor = axIds(i).Title.Color;
    % % % % % %     axPropert.Title.FontAngle = axIds(i).Title.FontWeight;
    % % % % % %
    % % % % % %
    % % % % % %     %% XLabel Properties
    % % % % % %     axPropert.XLabel.titleon = 0;
    % % % % % %     axPropert.XLabel.String = 'Introduce New X Axis Label';
    % % % % % %     axPropert.XLabel.Enable = 'off';
    % % % % % %     % Font Properties XLabel
    % % % % % %     axPropert.XLabel.FontName = axIds(i).XLabel.FontName;
    % % % % % %     axPropert.XLabel.FontSize = axIds(i).XLabel.FontSize;
    % % % % % %     axPropert.XLabel.FontWeight = axIds(i).XLabel.FontWeight;
    % % % % % %     axPropert.XLabel.FontColor = axIds(i).XLabel.Color;
    % % % % % %     axPropert.XLabel.FontAngle = axIds(i).XLabel.FontWeight;
    % % % % % %     axPropert.XAxis.gridon = axIds(i).XGrid;
    % % % % % %     axPropert.XAxis.Color = axIds(i).XColor;
    % % % % % %     axPropert.XAxis.numXTicks = length(axIds(i).XTick);
    % % % % % %
    % % % % % %
    % % % % % %     %% YLabel Properties
    % % % % % %     axPropert.YLabel.titleon = 0;
    % % % % % %     axPropert.YLabel.String = 'Introduce New Y Axis Label';
    % % % % % %     axPropert.YLabel.Enable = 'off';
    % % % % % %     % Font Properties YLabel
    % % % % % %     axPropert.YLabel.FontName = axIds(i).YLabel.FontName;
    % % % % % %     axPropert.YLabel.FontSize = axIds(i).YLabel.FontSize;
    % % % % % %     axPropert.YLabel.FontWeight = axIds(i).YLabel.FontWeight;
    % % % % % %     axPropert.YLabel.FontColor = axIds(i).YLabel.Color;
    % % % % % %     axPropert.YLabel.FontAngle = axIds(i).YLabel.FontWeight;
    % % % % % %     axPropert.YAxis.gridon = axIds(i).YGrid;
    % % % % % %     axPropert.YAxis.Color = axIds(i).YColor;
    % % % % % %     axPropert.YAxis.numXTicks = length(axIds(i).YTick);
    % % % % % %
    % % % % % %
    % % % % % %     %% ZLabel Properties
    % % % % % %     axPropert.ZLabel.titleon = 0;
    % % % % % %     axPropert.ZLabel.String = 'Introduce New Z Axis Label';
    % % % % % %     axPropert.ZLabel.Enable = 'off';
    % % % % % %     % Font Properties ZLabel
    % % % % % %     axPropert.ZLabel.FontName = axIds(i).ZLabel.FontName;
    % % % % % %     axPropert.ZLabel.FontSize = axIds(i).ZLabel.FontSize;
    % % % % % %     axPropert.ZLabel.FontWeight = axIds(i).ZLabel.FontWeight;
    % % % % % %     axPropert.ZLabel.FontColor = axIds(i).ZLabel.Color;
    % % % % % %     axPropert.ZLabel.FontAngle = axIds(i).ZLabel.FontWeight;
    % % % % % %     axPropert.ZAxis.gridon = axIds(i).ZGrid;
    % % % % % %     axPropert.ZAxis.Color = axIds(i).ZColor;
    % % % % % %     axPropert.ZAxis.numXTicks = length(axIds(i).ZTick);
    % % % % % %
    % % % % % %     axIds(i).UserData = axPropert;
    % % % % % %
    % % % % % %     cbar = findobj(hSurfs{i},'Tag','Colorbar');
    % % % % % %     if ~isempty(cbar)
    % % % % % %
    % % % % % %         colorbarProper.Title.titleon = 1;
    % % % % % %         colorbarProper.Title.String = 'Colorbar Title';
    % % % % % %         tempHand = title(cbar, colorbarProper.Title.String);
    % % % % % %         colorbarProper.Title.Handle = tempHand;
    % % % % % %         colorbarProper.Title.Handle.Color = [1 1 1] - param.figcolor;
    % % % % % %         colorbarProper.Title.Enable = 'on';
    % % % % % %         colorbarProper.Ticks.number = 11;
    % % % % % %         colorbarProper.FontName = cbar.FontName;
    % % % % % %         colorbarProper.FontWeight = cbar.FontWeight;
    % % % % % %         colorbarProper.FontSize = cbar.FontSize;
    % % % % % %         colorbarProper.FontAngle = cbar.FontAngle;
    % % % % % %         colorbarProper.FontColor = cbar.Color;
    % % % % % %         cbar.UserData = colorbarProper;
end
return;

function exportfigure(source,callbackdata, figID);

figUnits = 'centimeters';
figColor = [0.3 0.3 0.3];
set(0,'units',figUnits);
cm_screen = get(0,'screensize');

scaleFactorWidth = 7;                                                     % Figure Width scale factor
scaleFactorHeight = 5;                                                    % Figure Height scale factor

figWidth  = (cm_screen(3)-cm_screen(1))/scaleFactorWidth;                 % Figure width
figHeight = (cm_screen(4)-3*cm_screen(3)/cm_screen(4))/scaleFactorHeight; % Figure height
figCornerX = cm_screen(3)/2-(figWidth/2);                                 % Figure Corner (X coordinate)
figCornerY = cm_screen(4)/2-(figHeight/2);                                % Figure Corner (Y coordinate)
figPosition = [figCornerX figCornerY figWidth figHeight];                 % Figure Position
% if ~isempty(axPropert.ColorbarFig)&&ishandle(axPropert.ColorbarFig) % Detecting and closing old figures related with this axis
%     close(axPropert.ColorbarFig);
% end
% Creating figure



axPropert.ExportFig = figure('numbertitle',    'off',...
    'name',           'Export Figure ...',...
    'Color',          figColor,...
    'units',          figUnits,...
    'Position',       figPosition,...
    'InvertHardcopy', 'off',...
    'MenuBar',        'None',...
    'Resize',         'off');

% Creating Left Panel Properties
hand.panelID(1) = uipanel(axPropert.ExportFig,...
    'Position',            [.02 scaleFactorWidth/scaleFactorHeight*.02 .96 .96],...
    'Units',               'normalized',...
    'Title',               'Export properties',...
    'BackgroundColor',     figColor,...
    'HighlightColor',      figColor/4,...
    'ForegroundColor',     [1 1 1]);

hand.checkID(7) = uicontrol('Parent',                hand.panelID(1),...
    'Style',                 'text',...
    'String',                'Figure Filename',...
    'Units',                 'normalized',...
    'Position',              [.045 0.87    0.7    0.09],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);


hand.saveeditID = uicontrol('Parent',                hand.panelID(1),...
    'Tag',                   'intfigurefilename',...
    'Style',                 'edit',...
    'Enable',                'on',...
    'String',                'Introduce Figure Filename',...
    'Units',                 'normalized',...
    'Position',              [.045 .72 .71 .11]);


hand.buttonID(11) = uicontrol(hand.panelID(1),...
    'Style',                 'pushbutton',...
    'Value',                 0,...
    'String',                '...',...
    'Units',                 'normalized',...
    'Position',              [0.78 .715 0.2 0.13],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [1 1 1]-figColor,...
    'FontSize',              12,... 
    'FontWeight',            'bold',...
    'Callback',              {@changefigfilename});



hand.checkID(7) = uicontrol('Parent',                hand.panelID(1),...
    'Style',                 'text',...
    'String',                'Resolution (dpi)',...
    'Units',                 'normalized',...
    'Position',              [.045 0.535    0.4    0.09],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);


hand.resoeditID = uicontrol('Parent',                hand.panelID(1),...
    'Tag',                   'spatresoedit',...
    'Style',                 'edit',...
    'Enable',                'on',...
    'String',                '300',...
    'Units',                 'normalized',...
    'Position',              [.45 .532 .3 .11]);


hand.checkID(7) = uicontrol('Parent',                hand.panelID(1),...
    'Style',                 'text',...
    'String',                'Image Format',...
    'Units',                 'normalized',...
    'Position',              [.045 0.39    0.4    0.09],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       figColor,...
    'ForegroundColor',       [.81 .81 0]);



% Image Formats
imFormatsString = cellstr(strvcat(...
    '-djpeg100  (JPEG image, quality level of 100)',...
    '-dtiff     (TIFF with packbits compression)',...
    '-dps       (PostScript for black and white printers)',...
    '-dpsc      (PostScript for color printers)',...
    '-dps2      (Level 2 PostScript for black and white printers)',...
    '-dpsc2     (Level 2 PostScript for color printers)',...
    '-deps      (Encapsulated PostScript)',...
    '-depsc     (Encapsulated Color PostScript)',...
    '-deps2     (Encapsulated Level 2 PostScript)',...
    '-depsc2    (Encapsulated Level 2 Color PostScript)',...
    '-dpdf      (Color PDF file format)',...
    '-dsvg      (Scalable Vector Graphics)',...
    '-dpng      (Portable Network Graphic 24-bit truecolor image)',...
    '-dbmpmono  (Monochrome .BMP file format)',...
    '-dbmp256   (8-bit (256-color) .BMP file format)',...
    '-dbmp16m   (24-bit .BMP file format)',...
    '-dpcxmono  (Monochrome PCX file format)',...
    '-dpcx16    (Older color PCX file format (EGA/VGA, 16-color))',...
    '-dpcx256   (Newer color PCX file format (256-color))',...
    '-dpcx24b   (24-bit color PCX file format, 3 8-bit planes)'));

hand.popupID = uicontrol(hand.panelID(1),...
    'Style',               'popup',...
    'Tag',                 'imageformatpopup',...
    'String',              imFormatsString,...
    'Value',               1,...
    'Units',               'normalized',...
    'Position',            [.045 .29 .91 .07]);

hand.buttonID(12) = uicontrol(hand.panelID(1),...
    'Style',                 'pushbutton',...
    'Tag',                   'figsavebutton',...
    'Value',                 0,...
    'String',                'Save Figure',...
    'Units',                 'normalized',...
    'Position',              [0.045 0.01 0.91 0.18],...
    'HorizontalAlignment',   'left',...
    'BackgroundColor',       [0.6902    0.1647    0.1647],...
    'ForegroundColor',       [1 1 1],...
    'FontSize',              10,...
    'FontWeight',            'bold',...
    'Callback',              {@savefigure, figID});
return;

function changefigfilename(source,callback);

tempobj = findobj('Tag','intfigurefilename');
[filename, pathname, filterindex] = uiputfile( ...
    {'*.jpg',   'JPEG image (djpeg100)';...
    '*.tiff',  'TIFF image (dtiff)'; ...
    '*.fig',   'Matlab figure'; ...
    '*.ps',    'PostScript image (dps)';...
    '*.ps',    'PostScript image (dpsc)';...
    '*.ps',    'PostScript image (dps2)';...
    '*.ps',    'PostScript image (dpsc2)';...
    '*.ps',    'PostScript image (deps)';...
    '*.ps',    'PostScript image (depsc)';...
    '*.ps',    'PostScript image (deps2)';...
    '*.ps',    'PostScript image (depsc2)';...
    '*.pdf',   'PDF image (dpdf)'; ...
    '*.svg',   'SVG image (dsvg)'; ...
    '*.png',   'PNG image (dpng)'; ...
    '*.bmp',   'Bitmap image (dbmpmono)';...
    '*.bmp',   'Bitmap image (dbmp256)';...
    '*.bmp',   'Bitmap image (dbmp16m)';...
    '*.pcx',   'PCX image (dpcxmono)'; ...
    '*.pcx',   'PCX image (dpcx16)';...
    '*.pcx',   'PCX image (dpcx256)';...
    '*.pcx',   'PCX image (dpcx24b)';...
    '*.*',  'All Files (*.*)'}, ...
    'Save as');
switch filterindex
    case 1
        ext = '.jpg';
        imformat = '-djpeg100';
    case 2
        ext = '.tiff';
        imformat = '-dtiff';
    case 3
        ext = '.fig';
        imformat = '-dfig';
    case 4
        ext = '.ps';
        imformat = '-dps';
    case 5
        ext = '.ps';
        imformat = '-dpsc';
    case 6
        ext = '.ps';
        imformat = '-dps2';
    case 7
        ext = '.ps';
        imformat = '-dpsc2';
    case 8
        ext = '.ps';
        imformat = '-deps';
    case 9
        ext = '.ps';
        imformat = '-depsc';
    case 10
        ext = '.ps';
        imformat = '-deps2';
    case 11
        ext = '.ps';
        imformat = '-depsc2';
    case 12
        ext = '.pdf';
        imformat = '-dpdf';
    case 13
        ext = '.svg';
        imformat = '-dsvg';
    case 14
        ext = '.png';
        imformat = '-dpng';
    case 15
        ext = '.bmp';
        imformat = '-dbmpmono';
    case 16
        ext = '.bmp';
        imformat = '-dbmp256';
    case 17
        ext = '.bmp';
        imformat = '-dbmp16m';
    case 18
        ext = '.pcx';
        imformat = '-dpcxmono';
    case 19
        ext = '.pcx';
        imformat = '-dpcx16';
    case 20
        ext = '.pcx';
        imformat = '-dpcx256';
    case 21
        ext = '.pcx';
        imformat = '-dpcx24b';
    otherwise
        [pth,nm, ext] = fileparts(filename);
        switch deblank(ext)
            case '.jpg'
                imformat = '-djpeg100';
            case '.tiff'
                imformat = '-dtiff';
            case '.fig'
                imformat = '-dfig';
            case '.bmp'
                imformat = '-dpcx24b';
            case '.png'
                imformat = '-dpng';
            case '.pcx'
                imformat = '-dbmp16m';
            case '.pdf'
                imformat = '-dpdf';
            case '.ps'
                imformat = '-dpsc2';
            case '.svg'
                imformat = '-dsvg';
            otherwise
                error('Unrecognized image format');
                return
        end
end
if isempty(ext)
    ext = 'jpg';
    imformat = '-djpeg100';
end
if strfind(filename,ext)
    tempobj.String = [pathname filename];
else
    tempobj.String = [pathname filename ext];
end
tempobj.UserData = imformat;
return

function savefigure(source, callback, figId);
tempobj = findobj('Tag','intfigurefilename');
figFilename = tempobj.String;
imformat = tempobj.UserData

tempobj = findobj('Tag','spatresoedit');
spatres = tempobj.String;
set(0,'DefaultFigureColor','remove');

% if exist(which('export_fig'));
%     export_fig(figFilename,imformat, figId, '-nocrop','-transparent','-opengl',['-r' num2str(spatres)]);
%     
% else
    print( figId,imformat,['-r' num2str(spatres)],figFilename);
% end
return;

% % warning off;mkdir(DirPics);
% % filename = 'Left-CurvatureMap_SulcalLines-lateral.tiff';
% % Figurename = [DirPics filesep Id '-' filename];
% % export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r500' );

%% ====================== End of Internal functions ===================== %