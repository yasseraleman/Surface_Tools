function [varargout] = Separate_Independent_Surface(varargin);
%
% Syntax :
%   [Surfout] = Separate_Independent_Surface(Surf);
%
% This script divides a Surface matlab variable into individual surfaces
% according to its connectivity. Each cluster will be saved as a separated
% surface.
%
% Input Parameters:
%        Surf             :  Matlab Surface Variable
%
% Output Parameters:
%       Surfout           :  Separated Surfaces
%
% Related references:
%
%
% See also: 
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Cuban Neuroscience Center
% August 14th 2007
% Version $1.0

%% ====================== Checking input parameters ===================== %

if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
end
%% ========================= Checking input parameters ================================= %

%% =========================== Main Program ============================================ %

graphSu = surface2graph(Surf);
LabNet = Label_Graph_Components(graphSu);
txt = max(LabNet)';

sts = unique(txt);
Ns = length(sts);

for i = 1:Ns
    ind = find(txt==sts(i));

    % Surface Vertices
    Surft = Surf;

    % Selecting Faces
    Faces = Surf.SurfData.faces;
    a = ismember(Surf.SurfData.faces,ind);
    ind2 = find(sum(a')==0);
    Surft.SurfData.faces(ind2,:) = [];
    finalInd = Surft.SurfData.faces(:);

    % Reorganicing Surface
    [Surft] = Reorg_Surf(Surft);
    Surfout(i) = Surft;
end
%% ====================== End of Main Program ========================================== %

% Outputs
varargout{1} = Surfout;

return