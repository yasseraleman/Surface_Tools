function varargout = Compound_Surf_Objects(Surfa,varargin)
%
% Syntax :
% Surfout = Compound_Surf_Objects(Surfa);
%
% This function employs mesh_boolean (GPToolbox) to create compound objects by peforming boolean operations.
%
% Input Parameters:
%   Surfa       : Surfaces files.
%               : An atlas surface file is considered as a single surface.
%
% Output Parameters:
%   Surfj       : Output Surface.
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

%% ============================ Checking input parameters ============================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surfa = Surface_Checking(Surfa);
%     Parameters
    smoothIter = 0; % Smooth Iterations
    boolOperation = 'union'    ;      % Boolean variable to create new figure

end
% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'boolOperation' % Boolean Operation
                    boolOperation=varargin{2};
                    switch lower(boolOperation)
                        case 'union'
                            boolOperation = 'union';
                        case 'intersect'
                            boolOperation = 'intersect';
                        case 'minus'
                            boolOperation = 'minus';
                        case 'xor'
                            boolOperation = 'xor';
                        case 'resolve'
                            boolOperation = 'resolve';
                        otherwise
                            error('Wrong Boolean operation');
                    end
                    
                case 'smoothIter' % Smooth Iterations for boundary regions
                    smoothIter=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

% % % % % % % %% =================== End of checking input parameters ================= %


% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/testboolean.mat');
% Surfa = SulcSurfMat;

%% ======================== Main Program ================================ %
wh = whos('Surfa');

if (strcmp(wh.class,'struct'))
    Ns = 1;
elseif (strcmp(wh.class,'cell'))
    Surfa = Surfa(:);
    Ns = size(Surfa,1);
elseif ischar(Surfa(1,:));
    Ns = size(Surfa,1);
end
cont = 0;
allIds = 0;
allVertex = [0 0 0];
for j = 1:Ns
    
    if (strcmp(wh.class,'struct'))
        Surf = Surfa;
    elseif (strcmp(wh.class,'cell'))
        Surf = Surfa{j,1};
    elseif ischar(Surfa(1,:))
        Surf = Read_Surface(deblank(Surfa(j,:)));
    end
    for i = 1:length(Surf);
        Surft = Reorg_Surf(Surf(i));
        cont = cont + 1;
        if ~isfield(Surft, 'Is');
            Surft.Is = cont*ones(size(Surft.SurfData.vertices,1),1);
        end
        allIds = [allIds; Surft.Is ];
        allVertex = [allVertex;Surft.SurfData.vertices];
        if cont == 1
            Surfj.SurfData.vertices = Surft.SurfData.vertices;
            Surfj.SurfData.faces = Surft.SurfData.faces;
            % ================== Adding Scalar Values ======================= %
            if isfield(Surft, 'Is');
                Surfj.Is = Surft.Is;
            else
                Surfj.Is = cont*zeros(size(Surft.SurfData.vertices,1),1);
            end
            % ================== Adding Colors Field ======================== %
%             if isfield(Surf(i).SurfData, 'FaceVertexCData');
%                 Surfj.SurfData.FaceVertexCData = Surft.SurfData.FaceVertexCData;
%             else
%                 Surfj.SurfData.FaceVertexCData = ones(size(Surft.SurfData.vertices,1),3);
%             end
        else
            switch boolOperation
                case 'union'
                    [joinVertex,joinFaces, faceCorresp] = mesh_boolean(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces),double(Surft.SurfData.vertices),double(Surft.SurfData.faces),'union');
                    if isempty(joinFaces)
                        [joinVertex,joinFaces, faceCorresp] = mesh_boolean(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces),double(Surft.SurfData.vertices),double(Surft.SurfData.faces),'resolve');
                        [joinVertex,joinFaces,flip] = outer_hull(joinVertex,joinFaces);
                    end
                case 'intersect'
                    [joinVertex,joinFaces, faceCorresp] = mesh_boolean(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces),double(Surft.SurfData.vertices),double(Surft.SurfData.faces),'intersect');
                    
                case 'minus'
                    [joinVertex,joinFaces, faceCorresp] = mesh_boolean(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces),double(Surft.SurfData.vertices),double(Surft.SurfData.faces),'minus');
                case 'xor'
                    [joinVertex,joinFaces, faceCorresp] = mesh_boolean(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces),double(Surft.SurfData.vertices),double(Surft.SurfData.faces),'xor');
                case 'resolve'
                    [joinVertex,joinFaces, faceCorresp] = mesh_boolean(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces),double(Surft.SurfData.vertices),double(Surft.SurfData.faces),'resolve');
            end
            Surfj.SurfData.vertices = joinVertex;
            Surfj.SurfData.faces    = joinFaces;
        end
    end
end
% Surfj = Compound_Surf(Surfa);
% [joinVertex,joinFaces,flip] = outer_hull(double(Surfj.SurfData.vertices),double(Surfj.SurfData.faces));
%             Surfj.SurfData.vertices = joinVertex;
%             Surfj.SurfData.faces    = joinFaces;

allIds(1,:) = [];
allVertex(1,:) = [];

Surfj.Is = zeros(size(Surfj.SurfData.vertices,1),1);
[ind,ord] = ismember(Surfj.SurfData.vertices,allVertex,'rows');
Surfj.Is(find(ind)) =  allIds(ord(find(ind)));

if smoothIter
    FV=smoothpatch(Surfj.SurfData,1,smoothIter);
    ind = find(Surfj.Is ==0 );
    Surfj.SurfData.vertices(ind,:) = FV.vertices(ind,:);
end
Surfj = Surf_Corr(Surfj);

if sum(Surfj.Is-floor(Surfj.Is)) == 0
    %- Defining Neighborhoods %
    [Trip] = Vert_Neibp(double(Surfj.SurfData.faces),size(Surfj.SurfData.vertices,1),size(Surfj.SurfData.faces,1));
    Temp = sum(Trip);
    Trip(:,Temp==0) = [];
    temp = Trip(:,3:end);
    indz = find(temp == 0);
    temp(indz) = 1;
    
    
    % Defining Boundary points
    temp1 = Surfj.Is(temp);
    temp1(indz) =  max(temp1(:))+1;
    NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
    NewMat(indz) = 0;
    
    % Setting Boundary points = 0
    Surfj.Is(find(logical(sum(logical(NewMat)')'))) = 0;
    Surfj = Surf_Corr(Surfj);
end

%========================End of main program==============================%
% Outputs;
varargout{1} = Surfj;

return;