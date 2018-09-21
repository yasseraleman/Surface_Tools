function varargout = Reorganize_Edges(varargin);




SulcSurfMat = varargin{1};
Edges = varargin{2};

allEdges = [SulcSurfMat.SurfData.faces(:,1) SulcSurfMat.SurfData.faces(:,2); SulcSurfMat.SurfData.faces(:,2) SulcSurfMat.SurfData.faces(:,3);SulcSurfMat.SurfData.faces(:,1) SulcSurfMat.SurfData.faces(:,3)] ; % Edges from the intersection faces

Nfaces = size(SulcSurfMat.SurfData.faces,1);
allEdges = sort(allEdges')';  % Unifying edges to remove repeated edges
allEdges = [allEdges repmat([1:Nfaces]',[3 1])];

[indkeep,b] = ismember(Edges,allEdges(:,1:2),'rows');
faces = allEdges(b,3);
Nedg = size(Edges,1);
faces = zeros(Nedg,2);
for i  = 1:Nedg
    ind2 = find(ismember(allEdges(:,1:2),Edges(i,:),'rows'));
    faces(i,:) = allEdges(ind2,3)';
end

sts = unique(faces);
Nf = length(sts);
Graph = sparse(Nedg,Nedg);
for i = 1:Nf
    [X,Y] = find(faces == sts(i));
    if length(X)==2
        Graph(X(1),X(2)) = 1;
        Graph(X(2),X(1)) = 1;
    end
end

LabNet = Label_Graph_Components(Graph);
edgeIds = max(LabNet,[],2);
clustIds = nonzeros(unique(LabNet));
Nc = length(clustIds);
reordEdges = [0 0 0];
for i = 1:Nc
    ind = find(edgeIds == clustIds(i));
    subGraph = Graph(ind,ind);
    [Xel,Yel] = find(subGraph);
    subGraph(Xel(1),Yel(1)) = 0;
    subGraph(Yel(1),Xel(1)) = 0;
    start_ed = Xel(1);
    end_ed = Yel(1);
    [DIST, PATHS]=graphshortestpath(subGraph,start_ed,end_ed);
    reordEdges = [reordEdges;Edges(ind(PATHS),:) ones(length(PATHS),1)*i];
end
reordEdges(1,:) = [];
varargout{1} = reordEdges;
return;




