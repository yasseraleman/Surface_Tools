function  varargout = Create_Surface_Maps_Isolines(varargin);
%
% Syntax :
%     isoLines = Create_Surface_Maps_Isolines(Surf,distanceMap,numbIsolines);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surfa                : Surface file.
%   distanceMap          : Distance map to create the Isolines.
%   numbIsolines         : Number of Isolines (default: numbIsolines = 20).
%
% Output Parameters:
%   isoLines    : Isolines (Nx5 variable). First 3 columns are X, Y and Z
%                 coordinates. Column 4 is the number of clusters in the
%                 same iso line and the fifth column is the isoline number.
%   refIso      : Reference Isoline Number (1 for positive or negative distance maps and 
%                 different from one if the distance map contains positive and negative values).       
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
if nargin < 1 % the indispensable input arguments are not provided
    error('One input is mandatory');
    return;
elseif nargin < 2
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
    if isfield(Surf,'Is')
            distanceMap = Surf.Is;
    else
        error('There is not distance map');
        return;
    end
    numbIsolines =  20;
end
if nargin == 2
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
    distanceMap = varargin{2};
    numbIsolines =  20;
end
if nargin == 3
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
    distanceMap  = varargin{2};
    numbIsolines = varargin{3};
end
if nargin > 3
    error('To many inputs');
    return;
end

if nargout > 2
    error('To many outputs');
    return;
end

if length(distanceMap)~=size(Surf.SurfData.vertices,1)
    error('Different sizes between distance map and surfac vertices ');
    return;
end
%% ================== End of Checking input parameters ===================== %

%% ================== Main Program ===================== %
min_distance = min(distanceMap);
max_distance = max(distanceMap);
isolines_step = (max_distance-min_distance)/numbIsolines;
isolines_distances = min_distance: isolines_step : max_distance;
refIso = 1;

%  Changing Reference Isoline if distance map contains not only positive
%  positive values but also negative values.
if ~isempty(find(isolines_distances<0))&~isempty(find(isolines_distances>=0))
    posIsos = ceil(numbIsolines*max_distance/(max_distance-min_distance));
    posRange = linspace(0,max_distance,posIsos);
    negRange = linspace(min_distance,0,numbIsolines-posIsos);
    isolines_distances = [negRange(1:end-1) posRange];
    refIso = find(isolines_distances == 0);
    isolines_distances(refIso) = 0.05;
end


isolines_distances(1) = isolines_distances(1) + (isolines_distances(2)-isolines_distances(1))/3;
isolines_distances(end) = isolines_distances(end) + (isolines_distances(end)-isolines_distances(end-1))/3;

[LS,LD,I] = isolines(Surf.SurfData.vertices,Surf.SurfData.faces,distanceMap,isolines_distances);

% Isolines Labels
isoLines = [0 0 0 0 0];
sts = nonzeros(unique(I));
Niso = length(sts);
for i = 1:Niso
    
    % Selecting all lines belonging to the same isoline
    ind = find(I == sts(i)); 
    X = [LS(ind,1) LD(ind,1)]';
    Y = [LS(ind,2) LD(ind,2)]';
    Z = [LS(ind,3) LD(ind,3)]';
    
    % Minimum lines segments length
    distEdge = min(sqrt((X(:,2)-X(:,1)).^2 + (Y(:,2)-Y(:,1)).^2 +(Z(:,2)-Z(:,1)).^2)); 
    
    % Reorganicing lines segments
    [intCurve] = Reorganize_Interception_Line([X(:) Y(:) Z(:)]);
    
    % Distance between continuos points
    disTemp = sqrt(sum((intCurve(2:end,1:3)- intCurve(1:end-1,1:3)).^2,2)); 
    ind2del = find(disTemp < distEdge/10);
    intCurve(ind2del,:) = [];
    
 
    % [X Y Z ClustLabels IsoLabel]
    isoLines = [isoLines;intCurve sts(i)*ones(size(intCurve,1),1)];
end
isoLines(1,:) = [];
%% ================== End of Main Program ===================== %
% Outputs
varargout{1} = isoLines;
varargout{2} = refIso;
return;
