function [varargout] = Struct_Area_Comp(varargin);
%
% Syntax :
% [Surf] = Struct_Area_Comp(Surf);
%
% This function computes the parcelated structures area for a given surface.
%
% Input Parameters:
%   Surf       : Surface.
%
% Output Parameters:
%   Surf      : Output Surfaces. The area values for each structure is stored 
%  in  StructS field inside Surf struct.
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

uni = sort(unique(Surf.Is));
[At, Af] = Area_Comp(Surf);
Nstruct = length(uni);
Surf.StructS = zeros(Nstruct,3);

for k =1:Nstruct
    indvert = find(Surf.Is == uni(k));
    Surf.StructS(k,1) = uni(k);
    
    indexes = ismember(Surf.SurfData.faces, indvert);
    ind3 = find(sum(indexes,2) == 3);
    ind2 = find(sum(indexes,2) == 2);
    ind1 = find(sum(indexes,2) == 1);
    Surf.StructS(k,2) = sum(Af(ind3)) + 2/3*sum(Af(ind2)) + 1/3*sum(Af(ind1));
    Surf.StructS(k,3) = length(indvert) ;
end
%========================End of main program==============================%

% Outputs
varargout{1} = Surf;
return
