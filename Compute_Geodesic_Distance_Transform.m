function varargout = Compute_Geodesic_Distance_Transform(varargin)
%
% Syntax :
% distTransf = Compute_Geodesic_Distance_Transform(Surf, startpoints);
%
% This function computes the geodesic distance transform for a surface.
%
% Input Parameters:
%   Surf            : Surfaces file in Matlab Format.
%                     - Mandatory fields
%                       Surf.SurfData.vertices
%                       Surf.SurfData.faces
%   startpoints     : Selected source points.
%           
% Output Parameters:
%   distTransf       : Distance Transform vector.
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
    return;
else
    Surf =varargin{1};
    Surf = Surface_Checking(Surf);
    method = 'exact';
end
if nargin == 2
     bound_indexes = varargin{2};
     method = 'exact';
end
if nargin>2 
    method = varargin{3};
    switch method
        case 'exact'
            
        case 'manhattan'
            
        otherwise
            error('Wrong method');
            return;
    
    end
end

if nargout > 2
    error('To many outputs');
    return;
end

if nargout > 3
    error('Only three outputs are allowed');
    return;
end

%% =================== End of Checking input parameters ================= %

%% ======================== Main Program =============================== %%
if nargin <2
    bLevel = Compute_Level_from_Boundary(Surf); % Detecting vertex level (1: boundary vertices).
    bound_indexes = find(bLevel == 1);
else
    bound_indexes = varargin{2};
end
graphSurf = surface2graph(Surf);            % Creating Graph from surface triangles

switch method
    case 'exact'
        global geodesic_library;
        geodesic_library = 'geodesic_release';      %"release" is faster and "debug" does additional checks
        rand('state', 0);                         %comment this statement if you want to produce random mesh every time
        
        
        vertices = Surf.SurfData.vertices;
        faces = Surf.SurfData.faces;
        N = size(vertices,1);
        
        mesh = geodesic_new_mesh(vertices,faces);         %initilize new mesh
        algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm
        source_points = geodesic_create_surface_point('vertex',bound_indexes,vertices(bound_indexes,:));
        geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)
        
        
        clear vertex_id
        points_indexes = [1:N];
        destination = geodesic_create_surface_point('vertex',points_indexes,vertices(points_indexes,:));
        path = geodesic_trace_back(algorithm, destination);     %find a shortest path from source to destination
        distances = zeros(N,1);              %find distances to all vertices of the mesh (actual pathes are not computed)
        [source_id, distTransf] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1

        
    case 'manhattan'
        geoDist = graphallshortestpaths(graphSurf); % Geodesic Distance Matrix
        geoDist(logical(eye(size(geoDist)))) = Inf; %
        distTransf = min(geoDist(:,bound_indexes),[],2);      % Distance Transform
end

%% ==================== End of Main Program ============================ %%

% Outputs
varargout{1} = distTransf;

return;