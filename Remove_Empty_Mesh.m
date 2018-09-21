function varargout = Remove_Empty_Mesh(varargin);
%
% Syntax :
% outMesh = Remove_Empty_Mesh(inMesh, outMesh);
% 
% This function removes empty nodes from a specified mesh file. The
% corrected surface filename can be fixed in a second argument. Otherwise
% the corrected surface will be saved in the folder of the input mesh.
%
% Input Parameters:
%   inMesh       : Surface to correct.
%   outMesh      : Corrected surface.
%
% Output Parameters:
%   outMesh      : Output level labels. Points with Level = 1 are  starting points.
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

%% ============================= Checking Inputs ====================================== %%
if nargin < 1
    error('One input is mandatory');
    return;
elseif nargin == 1
    inMesh = varargin{1};
    [pth,nm,ext] = fileparts(inMesh);
    outMesh = [pth filesep nm '_mod.mesh'];
elseif nargin == 2
     inMesh = varargin{1};
     outMesh = varargin{2};
end
Surf = Surface_Checking(inMesh);
%% ========================== End of Checking Inputs ================================== %%


%% ================================ Main Program ====================================== %%

% Detecting empty nodes
contdelsurf = 0;
scont = 0;
indel = 0;
for j = 1:length(Surf)
    if isempty(Surf(j).Tri)
        contdelsurf = contdelsurf + 1;
        indel(contdelsurf) = j;
    else
        scont = scont + 1;
        sizes(scont) = size(Surf(j).Tri,1);
    end
end
if sum(indel) > 0
    % Removing empty Surfaces
    Surf(indel) = [];
end

% Saving the new mesh surface
save_mesh(Surf, outMesh);
%% ================================ Main Program ====================================== %%
% Outputs
varargout{1} = outMesh;
return;