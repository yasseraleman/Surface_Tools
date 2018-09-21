function varargout = Create_Region_ConvexHull_From_Surface(varargin);
%
% Syntax :
%   HullSurfMat = Create_Region_ConvexHull_From_Surface(Surf);
%
% This script creates hemispheric convex hull from surface using the
% boundary matlab function.
%
% Input Parameters:
%      Surf           :  Surfaces file.
%
% Output Parameters:
%     HullSurfMat     :  Hemispheric Convex Hull Surfaces
%
% Related references:
%
%
% See also: 
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% May 12th 2012
% Version $1.0

if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
end


%% ==================== ENTRY VERIFICATION ===============================%
% Surf = Read_Surface('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/surf/lh.pial','aparc.annot');


HullSurfMat = Compute_Hull_from_Surface(Surf,.4);

HullSurf.SurfData.vertices = HullSurfMat.SurfData.vertices;
HullSurf.SurfData.faces    = HullSurfMat.SurfData.faces;
clear HullSurfMat;
HullSurfMat = HullSurf;

if isfield(Surf,'Is')
    
    opts.verb = 0;opts.nsub = 2;
    [vertex,faces] = perform_mesh_subdivision(HullSurfMat.SurfData.vertices',HullSurfMat.SurfData.faces',opts.nsub,opts);
    HullSurfMat.SurfData.faces = faces';
    HullSurfMat.SurfData.vertices = vertex';

    indclosest = dsearchn(Surf.SurfData.vertices,HullSurfMat.SurfData.vertices);
    HullSurfMat.Is = Surf.Is(indclosest);
end
varargout{1} = HullSurfMat;
