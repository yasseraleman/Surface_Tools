function [varargout] = Cluster_Separation(varargin);
%
% Syntax :
%     [Surf] = Cluster_Separation(Surf);
%
% This function detects separated clusters inside the same structure and
% assigns a new label to each of them.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%                                         This structure variable most
%                                         contains a .Is field. This fields
%                                         storages the surface parcellation
%
% Output Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%                                         This structure variable contains a new .Is field. 
%
%
% See also: Surf_Sep
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 13th 2016
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
    error('The surface variable does not contain the labels field (.Is)');
    return;
end
% Is checking
try
    if (size(Surf.Is,1) ~= size(Surf.SurfData.vertices,1))
        error('Different sizes between Labels map and surface');
        return;
    end
catch
    error('Different sizes between Labels map and surface');
    return;
end

if nargin > 1
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%


% load('/media/COSAS/scripts/matlab.mat');
% Surf.Is = newparc;

%% ======================= Main Program ================================ %%

strl = unique(Surf.Is);strl(strl==0) = [];

Nc = length(strl);
newparc = zeros(size(Surf.SurfData.vertices,1),1);
for i = 1:Nc
    ind = find(Surf.Is == strl(i));
    relatFaces = Surf.SurfData.faces(find(sum(ismember(Surf.SurfData.faces,ind),2)~=0),:);
    
    relatEdges = [relatFaces(:,1) relatFaces(:,2);relatFaces(:,2) relatFaces(:,3);relatFaces(:,1) relatFaces(:,3)];
       
    [ordVertices,~,auxOrder] = unique(relatEdges(:));   % Locating the order in the original edge matrix
    Nedge = size(relatEdges,1);               % Number o edges
    reorgrelatEdges = reshape(auxOrder,[Nedge 2]);  % Reorganicing edges to avoid memory problems in the graph creatin step
    
    orderBranch=accumarray(reorgrelatEdges(:),reorgrelatEdges(:)*0+1);
    endpoints = find(orderBranch == 1);
    
    graphMat = sparse(max(reorgrelatEdges(:)),max(reorgrelatEdges(:))); % Creating empty Graph Matrix
    ind = sub2ind(size(graphMat),[reorgrelatEdges(:,1) reorgrelatEdges(:,2)],[reorgrelatEdges(:,2) reorgrelatEdges(:,1)]); % Establishing connections
    graphMat(ind) = 1;
    
    % Detecting Clusters
    graphMat = Label_Graph_Components(graphMat);
    newparc(ordVertices) = max(newparc) + max(graphMat);
    
    
end
%% ====================== End of Main Program ========================== %%
varargout{1} = Surf;
return