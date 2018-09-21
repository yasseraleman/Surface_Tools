function varargout = Surface_Region_Growing(varargin);
%
% Syntax :
%     [dilLabels] = Surface_Region_Growing(Surf, surfLabels);
%
% This function dilates labels throughout a specified surface.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%        surfLabels                     : Surface labels (N x 1 postive and integer vector ).
%
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

% if ~sum(surfLabels-floor(surfLabels)) ==0
%     error('Labels Map must be a postive and integer Npoints x 1 vector');
%     return;
% end
if sum(surfLabels)==0
    error('There are not initial regions in the specified label map');
    return;
end

if nargin > 2
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%
%load('/media/COSAS/scripts/Gyral_Crowns_and_Sulcal_Lines_Extraction/matlab.mat');
Surf = varargin{1};
surfLabels = varargin{2};

%% ============================= Main Program ========================== %%
disp(' ')
disp('Growing Surface Labels...');
tic;

%% =================== Computing Neighbor points ======================= %%
%[Surft, LabBranches] = Branch_Labelling(Surf, Slines);
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

% ================== Creating the Labels Matrix ===================== %
% This matrix contains information about the
% labels inside each point neighborhood

Labels_Mat = surfLabels(temp);
Labels_Mat(indz) = 0;
% =============== End of Creating the Labels Matrix ================= %

% Dilating clusters
dilLabels = surfLabels;
indn = find(dilLabels == 0);

dilLabelsant = dilLabels*0;
it = 0;
while ~isempty(indn)&(sum(dilLabelsant-dilLabels)~=0)
    dilLabelsant = dilLabels;
    it = it + 1;
    disp(['Iteration: ' num2str(it)]);
    vert_to_grow = find(sum(Labels_Mat,2) ~= 0); % Finding Labeled Neighbors
    vert_to_grow(ismember(vert_to_grow,find(dilLabels~=0))) = [];% Vertex to grow
    % vert_to_grow(find(sum(Labels_Mat(vert_to_grow,:),2) == Inf)) = []; % Restricting growing to noninfinite values
    
    
    % Labeled with the label of the majority of its labeled neighbors
    dilLabels(vert_to_grow) = arrayfun( @(vert_to_grow)mode(nonzeros(Labels_Mat(vert_to_grow,:))), vert_to_grow);
    
    
    % Updating the labels matrix
    Labels_Mat = dilLabels(temp);
    Labels_Mat(indz) = 0;
    
    % Detecting if someone is missing for labelling
    indn = find(dilLabels == 0);
end

% Labels Refilling
Surf.Is = dilLabels;
toc;
%% ====================== End of Main Program ========================== %%

% Outputs
varargout{1} = dilLabels;
return;