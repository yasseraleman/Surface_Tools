function varargout = Surface_Vox2mm(varargin)
%
% Syntax :
% Surf = Surface_Vox2mm(Surf,Vr,convcad);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surf        : Surface file (Matlab Format).
%   Vr          : Volume or image filename.
%   convcad     : Conversion string ('vox2mm'or 'mm2vox').
%
% Output Parameters:
%   hlist       : Handles list.
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
Surf = varargin{1};
Surf = Surface_Checking(Surf);

Vr =  varargin{2}; %
if ischar(Vr)
    if ~exist(Vr,'file')
        error('Please enter a correct image volume');
        return;
    end
    Vr = spm_vol(Vr);
elseif isstruct(Vr)
    Vr = Vr;
else
    error('Please enter a correct image volume');
    return;
end
convcad = varargin{3};
%% ================== End of Checking input parameters ================== %

%% ================================ Main Program ======================== %

switch convcad
    case 'mm2vox'
        tempVar = [inv(Vr.mat)*[Surf.SurfData.vertices ones(size(Surf.SurfData.vertices,1),1)]']';
    case 'vox2mm'
        tempVar = [Vr.mat*[Surf.SurfData.vertices ones(size(Surf.SurfData.vertices,1),1)]']';
end
Surf.SurfData.vertices = tempVar(:,1:3);

%% ======================== End of Main Program ========================= %
% Outputs
varargout{1} = Surf;

return;