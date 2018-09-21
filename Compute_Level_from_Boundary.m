function varargout = Compute_Level_from_Boundary(varargin);
%
% Syntax :
% bLevel = Compute_Level_from_Boundary(Surf);
%
% This function computes the neighbor levels for a surface mesh. The
% starting level (level = 1) are boundary points
%
% Input Parameters:
%   Surf       : Surface.
%
% Output Parameters:
%   bLevel      : Output level labels. Points with Level = 1 are boundary points.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% February 17th 2007
% Version $1.0

%% ============================= Checking Inputs ======================= %%
if nargin < 1
    error('One Input is mandatory');
    return
end
Surf = varargin{1};

% Surface Checking
Surf = Surface_Checking(Surf); 
if ~isfield(Surf,'Is')
    Surf.Is = ones(size(Surf.SurfData.vertices,1),1);
end

if nargout > 1
    errordlg('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%

%% ============================ Main Program =========================== %%

% ------------------ Computing Boundary Points -------------------------- %
% % % Remove Edges belonging to the gyral crown
% % allEdges = [Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,2); Surf.SurfData.faces(:,2) Surf.SurfData.faces(:,3);Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,3)] ; % Edges from the intersection faces
% % allEdges = sort(allEdges')'; % Sorting edges according to its point indexes
% % [reordEdges,b,order] = unique(allEdges,'rows');  % Unifying edges to remove repeated edges
% % tempVar = accumarray(order,ones(length(order),1)); % Detecting how many faces containg each edge. Boundary edges belong to a single face!
% % indboundEdges = find(tempVar == 1); % Boundary edges
% % gCrownEdges = reordEdges(indboundEdges,:); % Boundary edges (Level = 1)

gCrownEdges = Compute_Boundary_Edges(Surf);

% ------------------ End of Computing Boundary Points ------------------- %

% ------------------------- Computing Levels ---------------------------- %
% Neighborbood
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);

% Detecting other levels
bLevel = zeros(size(Surf.SurfData.vertices,1),1);
bLevel([gCrownEdges(:,1) gCrownEdges(:,2)]) = 1;

indn = find(bLevel == 0); % Unlabeled vertices
cont = 1;
while ~isempty(indn)
    cont = cont + 1;
    indlab = find(bLevel); % Labeled vertices
    neigh = temp(indlab,:);
    neigh = unique(nonzeros(neigh(:))); % Neighbors
    neigh(ismember(neigh,indlab)) = []; % Removing Labeled vertices
    bLevel(neigh) = cont;
    indn = find(bLevel == 0); % Unlabeled vertices
end

%% ========================= End of Main Program ======================= %%

% Outputs
varargout{1} = bLevel;
return;