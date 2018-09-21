function varargout = Geodesic_Map_Smoothing(varargin);
%
% Syntax :
%     smoothMap = Geodesic_Map_Smoothing(Surf, cortMap);
%
% This function smoothes cortical maps acording to its 
%
% Input Parameters:
%       Surf                    : Surface variable (file, struct or cellarray).
%       cortMap                 : Cortical map (file or Nx1 vector, where N is the 
%                                 number of points in its corresponding surface).
%       opts (optional input)   : Options:
%                                 opts.nSmooth - Smooth Iterations
%                                 opts.sMethod - Method (1: 'linear' 
%                                                        2: 'weighted'(1/r^2)
%                                                        3: 'gaussian'
%                                                        4: 'diffusion').
%
% Output Parameters:
%       smoothMap               : Smoothed cortical map (N x 1)
%
%
% See also: Surface_Checking Gyral_Crowns_Extraction Extracting_Sulcal_Line 
% Recursively_Remove_Branches Reordering_Gyral_Crowns Merge_Watershed_Regions 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0


%% ============================= Checking Inputs ======================= %%
if nargin < 2
    error('Two Inputs are mandatory');
    return
end
Surf = varargin{1};
cortMap = varargin{2};

% Surface Checking
Surf = Surface_Checking(Surf); 

if ischar(cortMap) % Verifying if the cortical map is a file
    if exist(cortMap,'file') % Verifying if the cortMap exists
        try
            [cortMap] = read_cfiles(cortMap); % Reading surface
        catch
            error('Unrecognized format');
            return;
        end
    end
end
try
    if (size(cortMap,1) ~= size(Surf.SurfData.vertices,1))
        error('Different sizes between cortical map and surface');
        return;
    end
catch
    error('Different sizes between cortical map and surface');
    return;
end

if nargin < 3
    opts.nSmooth = 1;        % Number of smoothing iterations
    opts.sMethod = 1 ;        % Smoothing method
elseif nargin == 3
    opts = varargin{3};
    if isstruct(opts)
        if ~isfield(opts,'nSmooth')
            opts.nSmooth = 1;        % Number of smoothing iterations
        end
        if ~isfield(opts,'sMethod')
            opts.sMethod = 1 ;        % Smoothing method
        end
    else
        error('Unrecognized options');
        return;
    end
end

if nargin > 3
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%
disp(' ')
disp('Smoothing Cortical Map...');
tic;
%% ========================= Main Program ============================== %%
% ---------------- Smoothing Curvature Map
% parameters for the operator
laplacian_type = 'distance';
options.symmetrize = 1;
options.normalize = 1; % it must be normalized for filtering
options.verb = 0;
%W = compute_mesh_weight(Surf.SurfData.vertices,Surf.SurfData.faces,laplacian_type,options);
A = triangulation2adjacency(Surf.SurfData.faces);
[i,j,s] = find(sparse(A));

switch opts.sMethod
    case 1
        d = sqrt(sum( (Surf.SurfData.vertices(i,:) - Surf.SurfData.vertices(j,:)).^2, 2));
    case 2
        d = sum( (Surf.SurfData.vertices(i,:) - Surf.SurfData.vertices(j,:)).^2, 2);
    case 3
        
    case 4        
        d = abs(cortMap(i) - cortMap(j));
end

W = sparse(i,j,d);
W(W>0) = 1./W(W>0);
W = (W+W')/2;
W = diag(sum(W,2).^(-1)) * W;

smoothMap = cortMap;
options.face_vertex_color = [];
 for it=1:opts.nSmooth
  %  curvmapd = (W*(W*curvmapd));
    smoothMap = (W*(W*smoothMap));
end
%% ======================= End of Main Program ========================= %%
% Outputs
varargout{1} = smoothMap;
return;