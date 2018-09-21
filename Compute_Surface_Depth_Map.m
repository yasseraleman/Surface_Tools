function varargout = Compute_Surface_Depth_Map(varargin);
%
% Syntax :
% depthMap = Compute_Surface_Depth_Map(Surf);
%
% This function computes the surface depth map.
%
% Input Parameters:
%   Surf            : Surfaces file in Matlab Format.
%                     - Mandatory fields
%                       Surf.SurfData.vertices
%                       Surf.SurfData.faces
%   method          : Selected method ('euclidean' or geodesic).
%           
% Output Parameters:
%   depthMap        : Depth Map.
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
    method = 'euclidean';
end
if nargin == 2
     method = varargin{2};
end

if nargin > 2
    error('To many inputs');
    return;
end

if nargout > 1
    error('To many outputs');
    return;
end


%% =================== End of Checking input parameters ================= %

%% ======================== Main Program =============================== %%

switch method
    case 'euclidean'
        % Computing Hull Surface
        HullSurfMat = Compute_Hull_from_Surface(Surf,.4);
        opts.verb = 0;opts.nsub = 2;
        [vertex,faces] = perform_mesh_subdivision(HullSurfMat.SurfData.vertices',HullSurfMat.SurfData.faces',opts.nsub,opts);
        HullSurfMat.SurfData.faces = faces';
        HullSurfMat.SurfData.vertices = vertex';
        
        % Computing Temporal Depth Map
        indeNear = dsearchn(HullSurfMat.SurfData.vertices,Surf.SurfData.vertices);
        depthMap = sqrt(sum((HullSurfMat.SurfData.vertices(indeNear,:)- Surf.SurfData.vertices).^2,2));
    case 'geodesic'

end

%% ==================== End of Main Program ============================ %%

% Outputs
varargout{1} = depthMap;