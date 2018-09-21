function varargout = test_Plot_Surf(Surfa,varargin)
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
%   hlist        : Handles list.
%   Surf         : Joined Surface
%   rangeColReal : Surface Indexes with real values in the surface map. It
%                  can be used to unify colormaps.
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
    Surfa = Surface_Checking(Surfa);
    % Parameters
    transpVal = 1; % Opacity
    colMap = 'jet'; % ColorMap
    newFig = 'y'; %
    lightIntensity = .25; % Light Intensity
    
    % Figure Properties
    param.figcolor = [0 0 0]; % Color
    param.figurevisible = 'on'; % Visible
    param.newfigure = 1    ;      % Boolean variable to create new figure
    
    set(0,'units','centimeters');
    cm_screen = get(0,'screensize');
    
    
    figPosition = [1 1 cm_screen(3)-3 cm_screen(4)-3*cm_screen(3)/cm_screen(4)]; % Position
    figUnits = 'centimeters'; % Units
end


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
                case 'lightIntensity'
                    lightIntensity = varargin{2};
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
% Detect figure and axis
if ~exist('FigID','var')
    FigID = figure('numbertitle','off','name','New Figure','Color', param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
    axesID = axes('Color',param.figcolor); grid on;
else
    if ~isempty(FigID)
        figure(FigID);
    else
        FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
        axesID = axes('Color',param.figcolor); grid on;
    end
end
%% =================== End of checking input parameters ================= %

%% ======================== Main Program ================================ %
[transpVal] = sort(transpVal);
wh = whos('Surfa');
if (strcmp(wh.class,'struct'))
    Ns = 1;
elseif (strcmp(wh.class,'cell'))
    Surfa = Surfa(:);
    Ns = size(Surfa,1);
elseif ischar(Surfa(1,:));
    Ns = size(Surfa,1);
end
transpVal((size(transpVal,2)+1):Ns) = max(transpVal);
transpVal = transpVal(1:Ns);
[transpVal,Surfa] = Reord_Tr(Surfa,transpVal,Ns,wh);
% [Surfa] = Reord_Tr(Surfa,Ns,wh);
transpVal = repmat(transpVal,[Ns 1]);
%col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
contVals = 0;
Osys = computer;
conth = 0; % handles counter
contnone = 0; % Counter for surfaces that do not contain Is maps.
contreal = 0; % Counter for surfaces that contains real values  in the Is maps.
continteger = 0; % Counter for surfaces that contains integer values in the Is maps.
context = 0;  % Counter along all surfaces 
for j = 1:Ns
    Surf = Surfa{j,1};
    re = floor(length(Surf)/Ncolor); col = repmat(col,[re+1 1]);
    if ~isfield(Surf(1),'SurfData');
        errordlg('This file does not contain surface information');
        continue;
    end
    %% ====================== Joining Surfaces ========================== %

        for i = 1:length(Surf);
            context = context + 1; % Counter along all surfaces
            
            % Trying to stablish a unique range for all surfaces in the
            % same axis
            if context == 1        % Initialization 
                Surftemp.SurfData.vertices = Surf(i).SurfData.vertices; % New Surface Vertices
                Surftemp.SurfData.faces = Surf(i).SurfData.faces; % New Surface Faces
                Colors = repmat(col(context,:),[size(Surf(i).SurfData.vertices,1) 1]);
                tempFaceVertexCData= Colors; % Final color data
                if isfield(Surf(i), 'Is')&~isfield(Surf(i).SurfData, 'FaceVertexCData'); % if the Is map exist then detect if it is integer(ie. atlas) or real(ie. thickness)
                    Surftemp.Is = Surf(i).Is;
                    if sum(Surf(i).Is - floor(Surf(i).Is))== 0 % Atlas condition
                        continteger = continteger + 1;
                        allIsInteger = Surf(i).Is;
                        surfLimitsInteger(continteger,:) = [1 length(Surf(i).Is)];
                    else % Real map condition
                        contreal = contreal + 1;
                        allIsReal = Surf(i).Is;
                        surfLimitsReal(contreal,:) = [1 length(Surf(i).Is)];
                    end
                    
                elseif ~isfield(Surf(i), 'Is')&~isfield(Surf(i).SurfData, 'FaceVertexCData') % if the Is map do no exist then apply a unique color
                    contnone = contnone + 1;
                    Surftemp.Is = ones(size(Surf(i).SurfData.vertices,1),1)*contnone;
                    allIsNone = ones(size(Surf(i).SurfData.vertices,1),1)*contnone;
                    surfLimitsNone(contnone,:) = [1 size(Surf(i).SurfData.vertices,1)];
                    
                elseif ~isfield(Surf(i), 'Is')&isfield(Surf(i).SurfData, 'FaceVertexCData') % if previous colordata exist and the Is map do no exist then use the old colordata
                    Surftemp.Is = ones(size(Surf(i).SurfData.vertices,1),1)*i;
                    tempFaceVertexCData = Surf(i).SurfData.FaceVertexCData;
                    
                elseif isfield(Surf(i), 'Is')&isfield(Surf(i).SurfData, 'FaceVertexCData') % if previous colordata a nd the Is map exist then use the old colordata and add the Is map to create a unified colordata
                    Surftemp.Is = Surf(i).Is;
                    tempFaceVertexCData = Surf(i).SurfData.FaceVertexCData;
                    if sum(Surf(i).Is - floor(Surf(i).Is))~= 0  % Real Data condition
                        contreal = contreal + 1;
                        allIsReal = Surf(i).Is;
                        surfLimitsReal(contreal,:) = [1 length(Surf(i).Is)];
                    end
                end
                
            else % Adding Values 
                if isfield(Surf(i), 'Is')&~isfield(Surf(i).SurfData, 'FaceVertexCData');
                    Surftemp.Is = [Surftemp.Is;Surf(i).Is];
                    if sum(Surf(i).Is - floor(Surf(i).Is))== 0
                        continteger = continteger + 1;
                        if continteger == 1
                            allIsInteger = Surf(i).Is;
                            surfLimitsInteger(continteger,:) = [1 length(Surf(i).Is)];
                        else
                            allIsInteger = [allIsInteger;Surf(i).Is];
                            surfLimitsInteger(continteger,:) = [size(Surftemp.SurfData.vertices,1)+1 size(Surftemp.SurfData.vertices,1)+length(Surf(i).Is)];
                        end
                    else
                        contreal = contreal + 1;
                        if contreal == 1
                            allIsReal = Surf(i).Is;
                            surfLimitsReal(contreal,:) = [1 length(Surf(i).Is)];
                        else
                            allIsReal = [allIsReal;Surf(i).Is];
                            surfLimitsReal(contreal,:) = [size(Surftemp.SurfData.vertices,1)+1 size(Surftemp.SurfData.vertices,1)+length(Surf(i).Is)];
                        end
                    end
                    Colors = Surf_Color(Surf(i),colMap);
                    tempFaceVertexCData = [tempFaceVertexCData;Colors];
                    
                   
                elseif ~isfield(Surf(i), 'Is')&~isfield(Surf(i).SurfData, 'FaceVertexCData')
                    contnone = contnone + 1;
                    Surftemp.Is = [Surftemp.Is;ones(size(Surf(i).SurfData.vertices,1),1)*contnone];
                    if contnone == 1
                        allIsNone = [ones(size(Surf(i).SurfData.vertices,1),1)*contnone];
                        surfLimitsNone(contnone,:) = [1 size(Surf(i).SurfData.vertices,1)];                   
                    else
                        allIsNone = [allIsNone;ones(size(Surf(i).SurfData.vertices,1),1)*contnone];
                        surfLimitsNone(contnone,:) = [size(Surftemp.SurfData.vertices,1)+1 size(Surftemp.SurfData.vertices,1)+size(Surf(i).SurfData.vertices,1)];
                    end
                    Colors = repmat(col(i,:),[size(Surf(i).SurfData.vertices,1) 1]);
                    tempFaceVertexCData = [tempFaceVertexCData;Colors];
                    
                elseif ~isfield(Surf(i), 'Is')&isfield(Surf(i).SurfData, 'FaceVertexCData')
                    tempFaceVertexCData = [tempFaceVertexCData;Surf(i).SurfData.FaceVertexCData];
                    
                elseif isfield(Surf(i), 'Is')&isfield(Surf(i).SurfData, 'FaceVertexCData')
                    Surftemp.Is = [Surftemp.Is;Surf(i).Is];
                    tempFaceVertexCData = [tempFaceVertexCData;Surf(i).SurfData.FaceVertexCData];
                    if sum(Surf(i).Is - floor(Surf(i).Is))~= 0
                        contreal = contreal + 1;
                        
                        if contreal == 1
                            allIsReal = Surf(i).Is;
                            surfLimitsReal(contreal,:) = [1 length(Surf(i).Is)];
                        else
                            allIsReal = [allIsReal;Surf(i).Is];
                            surfLimitsReal(contreal,:) = [size(Surftemp.SurfData.vertices,1)+1 size(Surftemp.SurfData.vertices,1)+length(Surf(i).Is)];
                        end
                    end
                end
                Surftemp.SurfData.vertices = [Surftemp.SurfData.vertices; Surf(i).SurfData.vertices]; % New Surface Vertices
                Surftemp.SurfData.faces = [Surftemp.SurfData.faces; Surf(i).SurfData.faces+max(Surftemp.SurfData.faces(:))]; % New Surface Faces
            end
            
            
            %% ================ End of Joining Surfaces ========================= %
        end
end
rangeColReal = 0;
if exist('surfLimitsReal','var')
    tempVar = [1;surfLimitsReal(:,2) - surfLimitsReal(:,1)];
    
    for k= 1:size(surfLimitsReal,1)
        rangeColReal = [rangeColReal eval([num2str(surfLimitsReal(k,1)) ':' num2str(surfLimitsReal(k,2))])];
    end
    Surft.Is = allIsReal;
    Colors = Surf_Color(Surft,colMap);
    rangeColReal(1) = [];
    tempFaceVertexCData(rangeColReal,:) = Colors;
end
rangeCol = 0;
if exist('surfLimitsInteger','var')
    for k= 1:size(surfLimitsInteger,1)
        rangeCol = [rangeCol eval([num2str(surfLimitsInteger(k,1)) ':' num2str(surfLimitsInteger(k,2))])];
    end
    Surft.Is = allIsInteger;
    Colors = Surf_Color(Surft,colMap);
    rangeCol(1) = [];
    tempFaceVertexCData(rangeCol,:) = Colors;
end
Surftemp.SurfData.FaceVertexCData = tempFaceVertexCData;

Surf = Surftemp; clear Surftemp;
%% ==================== Plotting Surfaces =========================== %

if size(Surf.SurfData.faces,2) == 3
    Surf.SurfData.FaceColor = 'interp';
    strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
    conth = conth + 1;
    hlist(conth) = strsurf;
    set(strsurf,'SpecularExponent',60);
elseif size(Surf.SurfData.faces,2) == 2
    poss = Surf.SurfData.faces;
    strlines = line([Surf.SurfData.vertices(poss(:,1),1) Surf.SurfData.vertices(poss(:,2),1)]',[Surf.SurfData.vertices(poss(:,1),2) Surf.SurfData.vertices(poss(:,2),2)]',[Surf.SurfData.vertices(poss(:,1),3) Surf.SurfData.vertices(poss(:,2),3)]','Color',col(j,:));
    conth = conth + 1;
    hlist(conth) = strlines(1);
else
    strlines = line(Surf.SurfData.vertices(:,1),Surf.SurfData.vertices(:,2),Surf.SurfData.vertices(:,3),'Color',col(j,:));
    conth = conth + 1;
    hlist(conth) = strlines(1);
end
% Transparent surfaces
if exist('strsurf','var')
    set(strsurf,'FaceAlpha',transpVal(j));
    set(strsurf,'SpecularExponent',60);
end
%% ==================== End of Plotting Surfaces ==================== %
%% ====================== Creating colorbar ============================= %
if contreal
    cent_col=(0-min(allIsReal))/(max(allIsReal)-min(allIsReal));
    
range = max(allIsReal)-min(allIsReal);
    values =  min(allIsReal):range/4:max(allIsReal);
    set(gca,'CLim',[min(values),max(values)]);
    colormap(gca,colormaps_colors(colMap,128,cent_col));
    
    hc = colorbar('peer',gca,'Limits',[min(values) max(values)]);
   % set(hc,'Ticks',values,'TickLabels',cellstr(num2str(values(:))));
   
    hc.Ticks = linspace(hc.Limits(1),hc.Limits(2),length(values));
    set(hc,'YTickLabel',sprintf('%.2f\n',values(:)));
    
    figColor = get(FigID,'Color');
    set(hc,'Color',([1 1 1] - figColor)*.8)
end
%% ====================== End of Creating colorbar ====================== %

%% ====================== Camera Angle and Lights ======================= %
axis image;
view(3);
% hc = camlight;
% set(hc, 'Color',[1 1 1]*lightIntensity);

[x,y, z] = sph2cart(35,0,1);
hl1 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);

[x,y, z] = sph2cart(125,0,1);
hl2 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);

[x,y, z] = sph2cart(215,0,1);
hl3 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);

[x,y, z] = sph2cart(305,0,1);
hl4 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);

[x,y, z] = sph2cart(0,90,1);
hl5 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);

[x,y, z] = sph2cart(0,270,1);
hl6 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
if contreal
    hlist = [hlist(:);[hl1 hl2 hl3 hl4 hl5 hl6 hc]'];
else
    hlist = [hlist(:);[hl1 hl2 hl3 hl4 hl5 hl6]'];
end
%% ===================== End of Camera Angle and Lights ================= %
%========================End of main program==============================%
% Outputs
varargout{1} = hlist;
varargout{2} = Surf;
varargout{3} = rangeColReal;
return

%=======================Internal functions================================%

function [trn,Surfa] = Reord_Tr(SurfF,transpVal,Ns,wh);
%function [Surfa] = Reord_Tr(SurfF,Ns,wh);
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
trn(t(end:-1:1)) = transpVal;

return
%=========================================================================%
%=========================Internal Functions===============================