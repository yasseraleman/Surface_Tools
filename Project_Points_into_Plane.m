function varargout = Project_Points_into_Plane(varargin);
%
% Syntax :
%  projCoords = Project_Points_into_Plane(Coords, planEq, pointCoord);
%
% This script projects a poit cloud onto a plane. The plane equation must be
% specified as a 1x4 vector ([A B C D]);
%
% Input Parameters:
%       Coords                  : Points coordinates (Nx3 array).    
%       planEq                  : Plane equation must be specified as a 
%                                 1x4 vector ([A B C D]);
%       pointCoord              : Coordinates to center of the plane ([Xc
%                                 Yc Zc]). If the coordinate center is not 
%                                 specified the the origin will be set as 
%                                 coordinate center.
%
% Output Parameters:
%       projCoords              : Projected points.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%

if nargin < 3
    error('The plane equation and points coordinates are mandatory');
    return
elseif nargin == 3
    Coords = varargin{1};
    planEq = varargin{2};
    pointCoord = varargin{3};
elseif nargin > 3
    error('Please enter only two inputs as maximum');
    return
end

%% ======================= End of Input parameters  ======================%

%% ============================== Main Program  ==========================%

 N = planEq(1:3)/norm(planEq(1:3)); % <-- do this if N is not normalized
 N2 = N.'*N;
 projCoords = Coords*(eye(3)-N2)+repmat(pointCoord*N2,size(Coords,1),1);

%% =========================== End of Main Program   =====================%
% Output
varargout{1} = projCoords;
return;