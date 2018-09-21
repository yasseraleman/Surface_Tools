function varargout = surface2graph(varargin);
%
% Syntax :
% [Graph] = surface2graph(Surf, bool);
%
% This scripts creates a sparse Graph from a surface matlab struct
%
% Input Parameters:
%       Surf                    : Matlab surface variable. 
%       bool                    : Boolean variable to create a logical
%                                 graph. All edges weights will be 1.
%
% Output Parameters:
%      Graph                    : Output Graph.
%
% See also:  
%______________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =================== Checking Input Parameters ======================= %%
if nargin == 0
    errormsg('Please enter a correct Surface struct');
    return;
end
Surf = varargin{1};
%% =================== End of Checking Input Parameters ================ %%
%% ========================== Main Program ============================= %%

Nv = length(Surf.SurfData.vertices); % Graph dimension
Graph = sparse(Nv,Nv); % Creating an empty sparse matrix dimension

Edges = [Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,2); Surf.SurfData.faces(:,2) Surf.SurfData.faces(:,3);Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,3)] ; % Edges
ind = sub2ind(size(Graph),[Edges(:,1);Edges(:,2)],[Edges(:,2);Edges(:,1)]); % Edges indexes

% Computing Neighborhood
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

% Creating Edges weights
if nargin == 2
    mapVector = varargin{2};
    % Logical graph
    if  length(mapVector)== 1
        Graph(ind) = 1;
    elseif length(mapVector)== size(Surf.SurfData.vertices,1)
        
        % Mean value between each point and its neighbors
        tempMat = ((mapVector(temp) +  [repmat(mapVector,[1 size(temp,2)])])/2).^-1 + eps;
        
        tempMat(indz)=0;
        tempMat = tempMat(:);
        
        % Creating weighted Graph
        temp(indz) = 0;
        Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)];
        ind = find(sum(logical(Coord),2) == 2);
        Coord = Coord(ind,:);
        tempMat = tempMat(ind);
        ind = sub2ind(size(Graph),Coord(:,1),Coord(:,2));
        Graph(ind) = tempMat; % Grafo de Distancia
    end
    
else
    % Edges weights are the euclidean distances between edges points
    
   
     % Computing Distance
    X = Surf.SurfData.vertices(temp,1);
    X = (reshape(X,size(temp)) - repmat(Surf.SurfData.vertices(Trip(:,1),1),[1 size(temp,2)])).^2;
    
    Y = Surf.SurfData.vertices(temp,2);
    Y = (reshape(Y,size(temp)) - repmat(Surf.SurfData.vertices(Trip(:,1),2),[1 size(temp,2)])).^2;
    
    Z = Surf.SurfData.vertices(temp,3);
    Z = (reshape(Z,size(temp)) - repmat(Surf.SurfData.vertices(Trip(:,1),3),[1 size(temp,2)])).^2;
    
    Distance=sqrt(X+Y+Z);
    Distance(indz)=0;
    Distance = Distance(:);
    
    % Creating weighted Graph
    temp(indz) = 0;
    Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)];
    ind = find(sum(logical(Coord),2) == 2);
    Coord = Coord(ind,:);
    Distance = Distance(ind);
    ind = sub2ind(size(Graph),Coord(:,1),Coord(:,2));
    Graph(ind) = Distance; % Grafo de Distancia
end
varargout{1} = Graph;
%% ==========================End of Main Program ============================= %%
return