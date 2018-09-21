function varargout = Compute_Geodesic_Center_of_Mass(varargin)
%
% Syntax :
% [centIndex,centCoord] = Compute_Geodesic_Center_of_Mass(Surf);
%
% This function computes the geodesic center of mass for a surface.
%
% Input Parameters:
%   Surf            : Surfaces file in Matlab Format.
%                     - Mandatory fields
%                       Surf.SurfData.vertices
%                       Surf.SurfData.faces
%           
% Output Parameters:
%   centIndex       : Index of the Center of Mass.
%   centCoord       : Center of Mass coordinates.
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
end
if nargout >2
    error('Only two outputs are allowed');
    return;
end
%% =================== End of Checking input parameters ================= %

%% ======================== Main Program =============================== %%
graphSurf = surface2graph(Surf);            % Creating Graph from surface triangles
geoDist = graphallshortestpaths(graphSurf); % Geodesic Distance Matrix
[~,centIndex] = min(sum(geoDist,2));        % Index of the Center of Mass   
centCoord = Surf.SurfData.vertices(centIndex,:);    % Center of Mass coordinates 
%% ==================== End of Main Program ============================ %%

% Outputs
varargout{1} = centIndex;
varargout{2} = centCoord;

return;