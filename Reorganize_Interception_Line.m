function varargout = Reorganize_Interception_Line(varargin);
%
% Syntax :
%     [reordCoords, reordintEdges] = Reorganize_Interception_Line(intCoords);
%
% This script reorganizes the interception points between surfaces.
%
% Input Parameters:
%       intCoords              : Interception points coordinates.
%
% Output Parameters:
%       reordCoords             : Reordered Coordinates. The third
%                                row shows the cluster Id for each
%                                point.
%       reordintEdges           : Reordered Edges. The third
%                                row shows the cluster Id for each
%                                point.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
% % load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test2.mat');

if nargin < 1
    error('There are missing input variables');
    return;
end

%% =========================== Main Program ============================= %
intCoords  = varargin{1}; % Interception LIne

eucDist = dist(intCoords');
eucDist(logical(eye(size(eucDist)))) = Inf;

Npoints = length(eucDist);
tempVar = [1:2:Npoints-1]';
intEdges = [tempVar tempVar+1];

[~,loc] = min(eucDist,[],2);
tempVar = [1:length(eucDist)]';

intEdges = [intEdges;tempVar loc];
intEdges =  unique(sort([intEdges]')','rows');
Graph = zeros(Npoints,Npoints);
ind2graph = unique(sub2ind(size(eucDist),[intEdges(:,1) ;intEdges(:,2)],[intEdges(:,2) ;intEdges(:,1)]));
Graph(ind2graph)= eucDist(ind2graph)+eps;



% Detecting if the sulcus is cutted in two pieces by the hull surface
LabNet = Label_Graph_Components(Graph);

idClusters = nonzeros(unique(LabNet(:))); % Number of possible 
Nc = length(idClusters); % Number of clusters
edgeId = max(full(LabNet));% Labelling each edge according to its cluster number
edgeId = edgeId(:);

%% ================= Reordering the interception edges ================== %
neworder = 0;
reordEdgeId = 0;
reordintEdges = [0 0 0];
cont = 0;
for i  = 1:Nc
    ind = find(edgeId == idClusters(i)); % Selecting the cluster
    subGraph = Graph(ind,ind); % Creating a subgraph for each cluster
    indextrems = find(sum(subGraph,2)==1);
    if isempty(indextrems)
        [Xrem,Yrem] = find(subGraph);
        start_point = Xrem(1); % Selecting the start and the endpoint
        end_point = Yrem(1);
        subGraph(start_point,end_point) = 0; % Disconecting the cluster
        subGraph(end_point,start_point) = 0;
        [DIST, PATHS]=graphshortestpath(sparse(subGraph),start_point,end_point); % Detecting the order between points
        PATHS = [PATHS(end) PATHS(1:end-1)];
        PATHS = PATHS(:);
        neworder = [neworder;ind(PATHS(:))]; % New order
        tempVar = ind(PATHS(:));
        reordintEdges  = [reordintEdges;[tempVar(1:end-1) tempVar(2:end);tempVar(end) tempVar(1)] i*ones(length(tempVar),1)];
    else
        cont = cont + 1;
%         distances=graphallshortestpaths(sparse(subGraph));
%         distances(find(distances==Inf))=0;
%         [start_point,end_point] = find(distances == max(distances(:)));
        [DIST, PATHS]=graphshortestpath(sparse(subGraph),indextrems(1),indextrems(2)); % Detecting the order between points
        
        if cont == 1
            extPairs(cont,:) = [1 length(PATHS) idClusters(i)];
        else
            extPairs(cont,:) = [length(neworder) length(neworder)+length(PATHS)-1 idClusters(i)];
        end
        neworder = [neworder;ind(PATHS(:))]; % New order
        tempVar = ind(PATHS(:));
        reordintEdges  = [reordintEdges;tempVar(1:end-1) tempVar(2:end) i*ones(length(tempVar)-1,1)];
    end
    reordEdgeId = [reordEdgeId;edgeId(ind(PATHS(:)))];
end
neworder(1) = [];
reordEdgeId(1) = [];
reordintEdges(1,:) = [];
reordCoords = [intCoords(neworder,:) reordEdgeId];


if exist('extPairs','var')
    if size(extPairs,1) >1
        
%         Graph = surface2graph(intSurf);
        distances = dist(reordCoords');
        distances=graphallshortestpaths(Graph);
        distances(find(distances==Inf))=0;
        [extremes,fillZero] = unique(extPairs(:,1:2));
        fillZero = reshape(fillZero,[2 size(extPairs,1)])';
        Nlen = length(extremes);
        indfillZero = sub2ind([Nlen Nlen],[fillZero(:,1);fillZero(:,2)],[fillZero(:,2);fillZero(:,1)]);
        temDist = distances(reordintEdges(extremes),reordintEdges(extremes));
        [Xextr,Yextr] = find(temDist == max(temDist(:)));
        
        
        [X,Y] = meshgrid(extremes,extremes);X = X(:);Y = Y(:);
        tempOrd = unique(sort([X Y]')','rows');
        ord2rem = find(ismember(tempOrd,[[extremes extremes];extPairs(:,1:2)],'rows'));
        tempOrd(ord2rem,:) =[];
        [a,b] = ismember(tempOrd,extremes);
        ind = sub2ind(size(temDist),b(:,1),b(:,2));
        
        distVect = [extremes(b) temDist(ind)];
        
        tempexT = extPairs(:,1:2);
        tempVal = ismember(tempexT,extremes(Xextr(1)));
        neighPos = find(tempVal - repmat(sum(tempVal,2),[1 2]) == -1);
        stNeigh = tempexT(neighPos);
        reordExtremes = [extremes(Xextr(1)) stNeigh 1];
        tempexT = [tempexT;flipdim(tempexT,2)];
        for i = 1:Nc-1
            oldstNeigh = stNeigh;
            [X,Y] = find(distVect(:,1:2) == oldstNeigh);
            [~,pos] = min(distVect(X,3));
            stNeight = distVect(X(pos),1:2);
            stNeight(stNeight == oldstNeigh) = [];
            
            [X,Y] = find(tempexT(:,1) == stNeight);
            stNeigh = tempexT(X,2);
            reordExtremes = [reordExtremes;stNeight stNeigh i+1];
            %         extPairs == Xextr
        end
        %
        flipBool = (reordExtremes(:,1) - reordExtremes(:,2));
        definitCoords(1,:) = [0 0 0];
        definitEdges(1,:) = [0 0 0];
        for i = 1:Nc
            partExtremes = reordExtremes(i,:);
            indrow = find(sum(ismember(extPairs(:,1:2),partExtremes(:)),2)== 2);
            clustId = extPairs(indrow,3);
            indTemp = find(reordEdgeId(:,3) == clustId);
            IntercCoordsTemp = reordCoords(indTemp,:);
            IntercEdgeTemp   = reordintEdges(indTemp,:);

            if flipBool(i)>0
                IntercCoordsTemp = flipdim(IntercCoordsTemp,1);
                IntercEdgeTemp   = flipdim(IntercEdgeTemp,1);

            end
            definitCoords = [definitCoords;IntercCoordsTemp(:,1:2) ones(size(IntercCoordsTemp,1),1)*i];
            definitEdges = [definitEdges;IntercEdgeTemp(:,1:2) ones(size(IntercEdgeTemp,1),1)*i];

        end
        definitCoords(1,:) = [];
        definitEdges(1,:) = [];
        reordCoords = definitCoords;
        reordintEdges  = definitEdges;

    end
end
%% =========================== End of Main Program ====================== %

varargout{1} = reordCoords;
varargout{2} = reordintEdges;

