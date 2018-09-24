function varargout = Erode_Surface_Label(varargin);
%
% Syntax :
%     [dilLabels] = Erode_Surface_Label(Surf, surfLabels, opts);
%
% This function erodes labels throughout a specified surface.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%        surfLabels                     : Surface labels (N x 1 postive and integer vector ).
%        opts (optional input)          : Options:
%                                         opts.iterations - Number of
%                                         iterations for the erotion procedure.
%                                         opts.fill - Boolean variable to
%                                         fill region wholes. (0: no
%                                         refill, 1: refill).
%                                         opts.reminter - Boolean variable
%                                         to remove labels from boundary
%                                         points.
% Output Parameters:
%         dilLabels                     : Dilated labels
%
%
% See also:
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
surfLabels = varargin{2};

% Surface Checking
Surf = Surface_Checking(Surf);

% Checking Labels Map
if ischar(surfLabels) % Verifying if the labels variable is a file
    if exist(surfLabels,'file') % Verifying if the labels map exists
        try
            [surfLabels] = read_cfiles(surfLabels); % Reading label
        catch
            error('Unrecognized curvature format');
            return;
        end
    end
end
try
    if (size(surfLabels,1) ~= size(Surf.SurfData.vertices,1))
        error('Different sizes between labels map and surface');
        return;
    end
catch
    error('Different sizes between labels map and surface');
    return;
end

if ~sum(surfLabels-floor(surfLabels)) ==0
    error('Labels Map must be a postive and integer Npoints x 1 vector');
    return;
end
if sum(surfLabels)==0
    error('There are not initial regions in the specified label map');
    return;
end

% Checking opts
if nargin < 3
    opts.iterations  = 1; % Number of iterations
    opts.fill = 1;       % Refilling
    opts.reminter =1;    % Remove boundary labels
    opts.verbose = 0;    % Verbose boolean variable
elseif nargin == 3
    opts = varargin{3};
    if isstruct(opts)
        if ~isfield(opts,'iterations')
            opts.iterations = 1;        % Number of iterations
        end
        if ~isfield(opts,'fill')
            opts.fill = 1;        % Refilling
        end
        if ~isfield(opts,'reminter')
            opts.reminter = 1;        % Mantaining connectivity
        end
        if ~isfield(opts,'verbose')
            opts.verbose = 0;        % Verbose boolean variable
        end
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

%% ============================= Main Program ========================== %%
if opts.verbose
    disp(' ')
    disp('Eroding Surface Labels...');
    tic;
end
%% ========================= Removing Medial Wall points =============== %%
indmwall = find(surfLabels == 4000);
surfLabels(indmwall) = 0;

%% =================== Computing Neighbor points ======================= %%
%[Surft, LabBranches] = Branch_Labelling(Surf, Slines);
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;


% Eroding clusters
for it = 1:opts.iterations
    if opts.verbose
        disp(['Iteration: ' num2str(it)]);
    end
    
    
    tempVar = repmat(surfLabels,[1 size(temp,2)]) - surfLabels(temp);
    tempVar(indz) = 0;
    vert_to_erode = find(sum(tempVar,2)~=0);
    
    surfLabels(vert_to_erode) = 0;
     
    
    
    
%     vert_to_grow = find(sum(Labels_Mat,2) ~= 0); % Finding Labeled Neighbors
%     vert_to_grow(ismember(vert_to_grow,find(dilLabels~=0))) = [];% Vertex to grow
%     
%     
%     % Labeled with the label of the majority of its labeled neighbors
%     dilLabels(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
%     
%     dilLabels(vert_to_grow) = arrayfun( @(vert_to_grow)mode(nonzeros(Labels_Mat(vert_to_grow,:))), vert_to_grow);
%     
%     % Updating the labels matrix
%     Labels_Mat = dilLabels(temp);
%     Labels_Mat(indz) = 0;
%     
%     % Detecting if someone is missing for labelling
%     indn = find(dilLabels == 0);
%     if isempty(indn)&(it ~= opts.iterations)
%         break;
%     end
end

% if opts.fill
%     % Labels Refilling
%     Surf.Is = dilLabels;
%     indzeros = find(Surf.Is == 0);
%     Surf.Is(indzeros) = max(Surf.Is) + 1;
%     [Surf] = Surf_Corr(Surf);
%     indzeros = find(Surf.Is == max(Surf.Is));
%     Surf.Is(indzeros) = 0;
%     dilLabels = Surf.Is;
% end


% % ------------------ Removing Interceptions ----------------------------- %
% if opts.reminter
%     indcrown = find(dilLabels);
%     indcrownfaces = find(sum(ismember(Surf.SurfData.faces,indcrown),2) ~= 0);
%     crownfaces = Surf.SurfData.faces(indcrownfaces,:);
%     crossEdges = [crownfaces(:,1) crownfaces(:,2); crownfaces(:,2) crownfaces(:,3);crownfaces(:,1) crownfaces(:,3)] ; % Edges from the intersection faces
%     labCrossEdges = dilLabels(crossEdges);
%     [X, Y] = find(labCrossEdges == 0);
%     crossEdges(X,:)  = [];
%     ind2rem = find(dilLabels(crossEdges(:,1)) - dilLabels(crossEdges(:,2)) ~=0);
%     vert2rem = crossEdges(ind2rem,:);
%     dilLabels(vert2rem(:)) = 0;
% end
% % ------------------- End of  Removing Interceptions -------------------- %
% 
% 
% 
% % ---------------- Parcellating Medial Wall
% dilLabels(indmwall) = 4000;
% if isfield(opts,'mwallind')
%     dilLabels(indmwall) = 4000;
% end
% if opts.verbose
%     toc;
% end
%% ====================== End of Main Program ========================== %%

% Outputs
varargout{1} = surfLabels;
return;