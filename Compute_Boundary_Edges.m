function varargout = Compute_Boundary_Edges(varargin);
%
% Syntax :
% boundEdges = Compute_Boundary_Edges(Surf);
%
% This function computes the bound edges of an open surface.
%
% Input Parameters:
%   Surf        : Surfaces files.
%       
% Output Parameters:
%   boundEdges   : Boudary edges list.
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
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
end
if nargout > 1
   error('To many outputs');
end
%% ====================== End of Checking input parameters ===================== %

%% =============================== Main Program ======================================== %


% Find all edges in mesh, note internal edges are repeated
E = sort([Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,2); Surf.SurfData.faces(:,2) Surf.SurfData.faces(:,3); Surf.SurfData.faces(:,3) Surf.SurfData.faces(:,1)]')';
% determine uniqueness of edges
[E,m,n] = unique(E,'rows');
% determine counts for each unique edge
counts = accumarray(n(:), 1);
% extract edges that only occurred once
boundEdges = E(counts==1,:);
%% ========================== End of Main Program ====================================== %
% Outputs
varargout{1} = boundEdges;

return