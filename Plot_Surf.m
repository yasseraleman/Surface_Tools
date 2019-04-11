function varargout = Plot_Surf(Surfa,varargin)
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
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surfa = Surface_Checking(Surfa);
    % Parameters
    transpVal = 1; % Opacity
    colMap = 'jet'; % ColorMap
    newFig = 'y'; %
    lightIntensity = .2; % Light Intensity
    noLight = 0; % Boolean variable to activate Lights
    boolFlag = 0; % Boolean variable to detect if the field FaceVertexCData exist
    
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
                    noLight = 1;
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
    FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
    axesID = axes('Color',param.figcolor); grid on;
else
    if ishandle(FigID)
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
col = [1 1 1; 1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
contVals = 0;
Osys = computer;
conth = 0; % handles counter
re = floor(length(Surfa)/Ncolor); col = repmat(col,[re+1 1]);
for j = 1:Ns
    Surf = Surfa{j,1};
    re = floor(length(Surf)/Ncolor); col = repmat(col,[re+1 1]);
    if ~isfield(Surf,'SurfData');
        errordlg('This file does not contain surface information');
        continue;
    end
    if length(Surf)~=1
        AllIs = 0;
        for i = 1:length(Surf);
            if isfield(Surf(i).SurfData,'VertexNormals');
                Surf(i).SurfData = rmfield(Surf(i).SurfData,'VertexNormals');
            end
            if size(Surf(i).SurfData.faces,2) == 3
                if isfield(Surf(i).SurfData,'FaceVertexCData')
                    Surf(i).SurfData.FaceColor = 'interp';
                    strsurf=patch(Surf(i).SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
                    conth = conth + 1;
                    hlist(conth) = strsurf;
                    set(strsurf,'SpecularExponent',60);
                elseif ~isfield(Surf(i).SurfData,'FaceVertexCData')&isfield(Surf(i),'Is')
                    if length(unique(Surf(i).Is))>1
                        Surf(i).SurfData.FaceColor = 'interp';
                        [Colors] = Surf_Color(Surf(i),colMap);
                        Surf(i).SurfData.FaceVertexCData = Colors;
                        strsurf=patch(Surf(i).SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
                        conth = conth + 1;
                        hlist(conth) = strsurf;
                        set(strsurf,'SpecularExponent',90);
                        if sum(Surf(i).Is-floor(Surf(i).Is)) ~=0
                            contVals = contVals  + 1;
                            maxVals(contVals) = max(Surf(i).Is);
                            minVals(contVals) = min(Surf(i).Is);
                        end
                        AllIs = [AllIs;Surf(i).Is];
                    else
                        if isfield(Surf(i),'Color')
                            col(i,:) = Surf(i).Color;
                        end
                        Surf(i).SurfData.FaceVertexCData = repmat(col(i,:),[length(Surf(i).SurfData.vertices) 1]);
                        strsurf=patch(Surf(i).SurfData,'facecolor',col(i,:),'edgecolor','none','tag', ...
                            'model0','facelighting','phong');
                        conth = conth + 1;
                        hlist(conth) = strsurf;
                        set(strsurf,'SpecularExponent',90);
                    end
                else
                    if isfield(Surf(i),'Color')
                        col(i,:) = Surf(i).Color;
                    end
                    Surf(i).SurfData.FaceVertexCData = repmat(col(i,:),[length(Surf(i).SurfData.vertices) 1]);
                    strsurf=patch(Surf(i).SurfData,'facecolor',col(i,:),'edgecolor','none','tag', ...
                        'model0','facelighting','phong');
                    conth = conth + 1;
                    hlist(conth) = strsurf;
                    set(strsurf,'SpecularExponent',90);
                end
                set(strsurf,'FaceAlpha',transpVal(j));
            elseif size(Surf(i).SurfData.faces,2) == 2
                poss = Surf(i).SurfData.faces;
                strlines = line([Surf(i).SurfData.vertices(poss(:,1),1) Surf(i).SurfData.vertices(poss(:,2),1)]',[Surf(i).SurfData.vertices(poss(:,1),2) Surf(i).SurfData.vertices(poss(:,2),2)]',[Surf(i).SurfData.vertices(poss(:,1),3) Surf(i).SurfData.vertices(poss(:,2),3)]','Color',col(i,:));
                conth = conth + 1;
                hlist(conth) = strlines(1);
            else
                strlines = line(Surf(i).SurfData.vertices(:,1),Surf(i).SurfData.vertices(:,2),Surf(i).SurfData.vertices(:,3),'Color',col(i,:));
                conth = conth + 1;
                hlist(conth) = strlines(1);
            end
            
            if exist('strsurf','var')
                if isfield(Surf(i),'Is')
                    dcmH = customDataCursor(strsurf, [Surf(i).SurfData.FaceVertexCData*255 Surf(i).Is]);
                else
                    dcmH = customDataCursor(strsurf, [Surf(i).SurfData.FaceVertexCData*255 zeros(length(Surf(i).SurfData.vertices),1)]);
                end
            end
            
        end
        AllIs(1,:) = [];
        [uord,ord] = unique(AllIs);
        
        colMap = Val2colors(AllIs(ord),colMap);
    else
        if isfield(Surf,'Color')
            col(j,:) = Surf.Color;
            Surf.SurfData.FaceVertexCData = repmat(Surf.Color,[size(Surf.SurfData.vertices,1) 1]);
        end
        if isfield(Surf.SurfData,'VertexNormals');
            Surf.SurfData = rmfield(Surf.SurfData,'VertexNormals');
        end
        if ~isfield(Surf,'Is')&~isfield(Surf.SurfData,'FaceVertexCData')
            if size(Surf.SurfData.faces,2) == 3
                Surf.SurfData.FaceVertexCData = repmat(col(j,:),[length(Surf.SurfData.vertices) 1]);
                strsurf=patch(Surf.SurfData,'facecolor',col(j,:),'edgecolor','none','tag', ...
                    'model0','facelighting','phong');
                conth = conth + 1;
                hlist(conth) = strsurf;
                set(strsurf,'SpecularExponent',90);
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
        elseif isfield(Surf,'Is')&~isfield(Surf.SurfData,'FaceVertexCData')
            %% Plot using vertices
            [Colors] = Surf_Color(Surf,colMap);
            
            Surf.SurfData.FaceVertexCData = Colors;
            if sum(Surf.Is-floor(Surf.Is)) ==0
                ind = find(Surf.Is == 0);
                Surf.SurfData.FaceVertexCData(ind,:) = repmat([1 1 1],[length(ind) 1]);
            end
            Surf.SurfData.FaceColor = 'interp';
            if size(Surf.SurfData.faces,2) == 3
                strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','phong');
                conth = conth + 1;
                hlist(conth) = strsurf;
                set(strsurf,'SpecularExponent',90);
            elseif size(Surf.SurfData.faces,2) == 2
                poss = Surf.SurfData.faces;
                ColorsLine = Surf_Color(Surf,colMap);
                ColorsLine = [ColorsLine;ColorsLine(1,:)];
                ColorsLine = (ColorsLine(1:end-1,:) + ColorsLine(2:end,:))/2;
                hold on;
                for k = 1:size(ColorsLine,1);
                    strlines(k) = plot3(Surf.SurfData.vertices(poss(k,:),1),Surf.SurfData.vertices(poss(k,:),2),Surf.SurfData.vertices(poss(k,:),3),'Color',ColorsLine(k,:),'Linewidth',7);
                end
                conth = conth + 1;
                hlist(conth) = strlines(1);
            else
                [ColorsLine] = Surf_Color(Surf,colMap);
                ColorsLine = Val2colors(Surf.Is,colMap);
                ColorsLine = [ColorsLine;ColorsLine(1,:)];
                ColorsLine = (ColorsLine(1:end-1,:) + ColorsLine(2:end,:))/2;
                poss = [[1:length(Surf.Is)-1]' [2:length(Surf.Is)]';[length(Surf.Is) 1]];
                hold on;
                for k = 1:size(ColorsLine,1);
                    strlines(k) = plot3(Surf.SurfData.vertices(poss(k,:),1),Surf.SurfData.vertices(poss(k,:),2),Surf.SurfData.vertices(poss(k,:),3),'Color',ColorsLine(k,:),'Linewidth',7);
                end
                conth = conth + 1;
                hlist(conth) = strlines(1);
            end
            if sum(Surf.Is-floor(Surf.Is)) ~=0
                contVals = contVals  + 1;
                maxVals(contVals) = max(Surf.Is);
                minVals(contVals) = min(Surf.Is);
                [uord,ord] = unique(Surf.Is);
                colMap = Surf.SurfData.FaceVertexCData(ord,:);
                
                %                 h = colorbar;
                %                 colormap(colormaps_colors(colMap,64))
                %                 range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10: max(Surf(1).Is);
                %                 set(h,'YTickLabel',num2str(values'));
            end
%             Surf.SurfData = rmfield(Surf.SurfData,'FaceVertexCData');
            boolFlag = 1;
        elseif isfield(Surf,'Is')&isfield(Surf.SurfData,'FaceVertexCData')&boolFlag
            if isfield(Surf,'Is')
                if sum(Surf.Is-floor(Surf.Is)) ~=0
                    contVals = contVals  + 1;
                    maxVals(contVals) = max(Surf.Is);
                    minVals(contVals) = min(Surf.Is);
                    [uord,ord] = unique(Surf.Is);
                    colMap = Surf.SurfData.FaceVertexCData(ord,:);
                end
            end
            
            Surf.SurfData.FaceColor = 'interp';
            if size(Surf.SurfData.faces,2) == 3
                strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
                conth = conth + 1;
                hlist(conth) = strsurf;
                set(strsurf,'SpecularExponent',90);
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
        elseif ~isfield(Surf,'Is')&isfield(Surf.SurfData,'FaceVertexCData')
                        
            Surf.SurfData.FaceColor = 'interp';
            if size(Surf.SurfData.faces,2) == 3
                strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
                conth = conth + 1;
                hlist(conth) = strsurf;
                set(strsurf,'SpecularExponent',90);
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
            
        end
         if exist('strsurf','var')
                if isfield(Surf,'Is')
                    dcmH = customDataCursor(strsurf, [Surf.SurfData.FaceVertexCData*255 Surf.Is]);
                else
                    dcmH = customDataCursor(strsurf, [Surf.SurfData.FaceVertexCData*255 zeros(length(Surf.SurfData.vertices),1)]);
                end
            end
        
    end
    if exist('strsurf','var')
        set(strsurf,'FaceAlpha',transpVal(j));
        set(strsurf,'SpecularExponent',90);
    end
end
if contVals
    hc = colorbar('peer',gca);
    if isstr(colMap)
        colormap(gca,colormaps_colors(colMap,64));
    else
        colormap(gca,colMap);
    end
    
    if exist('ord','var')
        if length(ord) < 11
            nTicks = length(ord);
        else
            nTicks = 11;
        end
        interValues = floor(linspace(1,length(ord),nTicks));
        
    else
        nTicks = 11;
        interValues = floor(linspace(1,11,nTicks));
    end
    
    values = uord(interValues);
    tempc = colMap(interValues,:);
    defvalues = linspace(values(1),values(end),nTicks); % Definitive Tick Values
    
    if length(defvalues) >1
        col = interp1(values,tempc,defvalues); % Definitive ColorValues
        range = max(defvalues)-min(defvalues);
        valuesInterp =  min(minVals):range/255:max(maxVals);
        colMap = interp1(defvalues,col,valuesInterp);
    else
        colMap = tempc;
    end
    colormap(gca,colMap);
    hc.Ticks = linspace(hc.Limits(1),hc.Limits(2),nTicks);
    set(hc,'YTickLabel',num2str(defvalues(:)));
    %
    %     set(hc,'YTickLabel',num2str(values'));
    
    figColor = get(FigID,'Color');
    set(hc,'Color',([1 1 1] - figColor)*.8)
end
axis image;
axis off;

% hc = camlight;
% set(hc, 'Color',[1 1 1]*lightIntensity);
AxID = gca;
AxID.Clipping = 'off';
if isempty(findall(AxID.Children,'Type','light'))
    
    view(3);
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
    % if contVals
    %     hlist = [hlist(:);[hl1 hl2 hl3 hl4 hl5 hl6 hc]'];
    % else
    hlist = [hlist(:);[hl1 hl2 hl3 hl4 hl5 hl6]'];
    % end
end
%========================End of main program==============================%
% Outputs;
varargout{1} = hlist;

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

function dcmH = customDataCursor(h, datalabels)
%  Modified from customDataCursor.m created by Dave Van Tol.
%
% customDataCursor allows a user to label each point in a data series with
% custom labels. When the data cursor mode is enabled (DATACURSORMODE ON),
% clicking on a data point will display the custom label, then the x and y
% values of the data point. If more than one data point of the data series
% is in the same location, the data cursor will display all of the labels
% first, then the x and y locations.
%
% Usage: customDataCursor(H,DATALABELS) applies the custom labels found in
% cell array DATALABELS to the line with handle H. DATALABELS must be a
% cell array with the same length as 'XDATA' and 'YDATA' in the line; in
% other words, there must be a label for each data point in the line.
%
% dcmH = customDataCursor(H,DATALABELS) returns the handle to the data
% cursor in dcmH.
%
% Note: CUSTOMDATACURSOR uses the 'UserData' property of a line to store
% and retrieve the labels. If that property is used by another portion of a
% user's program, this function may be broken, or this function may break
% the user's program.
% 
% Example:
%   % generate some data and chart it
%   N = 20;
%   x = rand(N,1);
%   y = rand(N,1);
%   h = plot(x,y,'ko');
%   % make up some labels and apply them to the chart
%   labels = cell(N,1);
%   for ii = 1:N
%       labels{ii} = ii;
%   end
%   customDataCursor(h,labels)
%
% Copyright (c) 2008, Dave Van Tol
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% Check input arguments
if ~exist('h','var')||~exist('datalabels','var')
    error('Improper inputs to function customDataCursor')
elseif ~ishandle(h)
    error('Nonexistent handle passed to function customDataCursor')
elseif length(datalabels) ~= length(get(h,'Vertices'))
    error(['Error in input to function customDataCursor: '...
        'number of labels is different than the number of data points in the line'])
end

% Put the labels in the 'userdata' property of the line
set(h,'userdata',datalabels)
% find the handle for the data cursor; set the update function 
dcmH = datacursormode(gcf); 
set(dcmH,'UpdateFcn',@cursorUpdateFcn)

function output_txt = cursorUpdateFcn(obj,cursor_info)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

% read in the labels from 'UserData' of the line
labels = get(get(cursor_info,'Target'),'UserData');
% read the x and y data
vertCoords = get(get(cursor_info,'Target'),'Vertices');
datapoint= find(ismembertol(vertCoords,cursor_info.Position,0.01,'ByRows',true));

% now figure out which data point was selected
% datapoint = find( (xvals==pos(1))&(yvals==pos(2)) );

% create the text to be displayed
output_txt = { ['Label: ' num2str(labels(datapoint(1),4))];...
    ['Color: ' num2str(round(labels(datapoint(1),1))) ' ' num2str(round(labels(datapoint(1),2))) ' ' num2str(round(labels(datapoint(1),3)))];...
    ['X: ',num2str(cursor_info.Position(1),4)];...
    ['Y: ',num2str(cursor_info.Position(2),4)];... 
    ['Z: ',num2str(cursor_info.Position(3),4)]};
% If there is a Z-coordinate in the position, display it as well