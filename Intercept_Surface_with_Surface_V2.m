function varargout = Intercept_Surface_with_Surface_V2(varargin);
%
% Syntax :
%       [intCurve, intNormals, edgeLabels] = Intercept_Surface_with_Surface(estSurface,  intSurface);
%
% This function computes the interception line between two surfaces
%
%
% Input Parameters:
%       estSurface              : Static Surface
%       intSurface              : Intercepting Surface.
%
% Output Parameters:
%      intCurve                 : Interception Curve.
%      intNormals               : Normals at interception curve.
%      edgeLabels               : Interception edges labeled according to
%                                 the its endpoints labels the interception
%                                 surface.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
% % if nargin < 2
% %     error('Three Inputs are needed');
% %     return
% % end
estSurface = varargin{1};
intSurface = varargin{2};
epsilon = 10^-5;

number_of_isolines = 10;
% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_single.mat');
%% =========================== Main Program ============================= %
% Surface planes
if ~isfield(estSurface.SurfData,'VertexNormals')
    estSurface = Compute_Surface_Normals(estSurface);
end
if ~isfield(intSurface.SurfData,'VertexNormals')
    intSurface = Compute_Surface_Normals(intSurface);
end

intFaces = intSurface.SurfData.faces;
allEdges = [intFaces(:,1) intFaces(:,2); intFaces(:,2) intFaces(:,3);intFaces(:,1) intFaces(:,3)] ; % Edges from the intersection faces

% Mean edge length
meanEdgeLength = mean(sqrt(sum((intSurface.SurfData.vertices(allEdges(:,2),:) - intSurface.SurfData.vertices(allEdges(:,1),:)).^2,2))) ;


[Surfout] = Separate_Independent_Surface(intSurface);
Nn = length(Surfout);
for k = 1:Nn
    parcIds = Surfout(k).Is;
    [fF] = intersect_other(estSurface.SurfData.vertices,estSurface.SurfData.faces, Surfout(k).SurfData.vertices, double(Surfout(k).SurfData.faces));
    if ~isempty(fF)
        Surft1.SurfData.vertices = Surfout(k).SurfData.vertices;
        Surft1.Is = Surfout(k).Is;
        Surft1.SurfData.faces = Surfout(k).SurfData.faces(fF(:,2),:);
        Surft1 = Reorg_Surf(Surft1);
%         Surfa{1,1} = Surft1;
%         Surft1 = Surf_Corr(Surft1);
%         Surfa{2,1} = Surft1;
        
        
        Surfh.SurfData.vertices = estSurface.SurfData.vertices;
        Surfh.SurfData.faces = estSurface.SurfData.faces(fF(:,1),:);
        Surfh.SurfData.VertexNormals = estSurface.SurfData.VertexNormals;
        Surfh = Reorg_Surf(Surfh);
        
        
        % Bottom surface of the sulcus (Part below the hull surface)
        [W,H, J] = mesh_boolean(estSurface.SurfData.vertices,estSurface.SurfData.faces,Surfout(k).SurfData.vertices,Surfout(k).SurfData.faces, 'intersect');
        
        % Top surface of the sulcus (Part above the hull surface)
        [W1,H1, J] = mesh_boolean(Surfout(k).SurfData.vertices,Surfout(k).SurfData.faces,estSurface.SurfData.vertices,estSurface.SurfData.faces, 'minus');
        
        % Interception between both parts of the sulci
        intercVert = ismember(W,W1,'rows');
        
        % Detecting common faces from both surfaces
        ind = find(intercVert);
        infacesIntercep = find(sum(ismember(H,ind),2) == 3);
        
        % Detecting faces from below part of the surface
        tempVar = ismember(H,ind);
        infacesBelow = find(sum(tempVar,2) ~= 3);
        
        clear Surftest;
        
        %% Common Surface
        Surftest.SurfData.vertices = W;
        Surftest.SurfData.faces = H(infacesIntercep,:);
        Surftest = Reorg_Surf(Surftest);
        
        %% Below Surface
        Surfbelow.SurfData.vertices = W;
        Surfbelow.SurfData.faces = H(infacesBelow,:);
        Surfbelow = Reorg_Surf(Surfbelow);
        distancesBelow = Compute_Geodesic_Distance_Transform(Surfbelow);
        
        
        %% Above Surface
        intercVert = ismember(W1,W,'rows');
        
        % Detecting common faces from both surfaces
        ind = find(intercVert);
        infacesIntercep = find(sum(ismember(H1,ind),2) == 3);
        
        % Detecting faces from below part of the surface
        tempVar = ismember(H1,ind);
        infacesAbove = find(sum(tempVar,2) ~= 3 & sum(tempVar,2) == 0);
        
        % Below Surface
        Surfabove.SurfData.vertices = W1;
        Surfabove.SurfData.faces = H1(infacesAbove,:);
        Surfabove = Reorg_Surf(Surfabove);
        distancesAbove = Compute_Geodesic_Distance_Transform(Surfabove);
        
        
        %% Refilling Labelling
        allDist = [-1*distancesAbove;distancesBelow];
        allVerts = [Surfabove.SurfData.vertices;Surfbelow.SurfData.vertices];
        
        locbottom = dsearchn(allVerts, Surfout(k).SurfData.vertices);
        Surfout(k).Is = allDist(locbottom);
        
        %% Creating Surface Isolines
        Niso = 2*ceil((max(Surfout(k).Is)-min(Surfout(k).Is))/meanEdgeLength);
        [isoLines, zeroIso] = Create_Surface_Maps_Isolines(Surfout(k), Surfout(k).Is, Niso);
        
% %         %% Moving reference Surface
% %         stsIsos = unique(isoLines(:,5));
% %         
% %         [col] = colorGradient([255 255 255],[255 255 0],length(stsIsos));
% %         col(zeroIso,:) = [0 1 1];
% %         Surfout(k).SurfData.FaceVertexCData = Surf_Color(Surfout(k));
% %         Plot_Surf(Surfout(k))
% %         for isoIndex = 1:length(stsIsos)
% %             indIso = find(isoLines(:,5) == stsIsos(isoIndex));
% %             isoCoords = isoLines(indIso,:);
% %             for k = 1:max(isoCoords(:,4))
% %                 indClust = find(isoCoords(:,4) == k);
% %                 hold on;
% %                 plot3(isoCoords(indClust,1), isoCoords(indClust,2),isoCoords(indClust,3),'.','Color',col(isoIndex,:), 'LineWidth', 3);
% %             end
% %         end
        
        % Labelling Surface
        locbottom = dsearchn(Surfout(1).SurfData.vertices,Surftest.SurfData.vertices);
        Surftest.Is = parcIds(locbottom);
        
        [Trip] = Vert_Neibp(double(Surftest.SurfData.faces),size(Surftest.SurfData.vertices,1),size(Surftest.SurfData.faces,1));
        Temp = sum(Trip);
        Trip(:,Temp==0) = [];
        temp = Trip(:,3:end);
        indz = find(temp == 0);
        temp(indz) = 1;
        
        temp1 = Surftest.Is(temp);
        temp1(indz) =  max(temp1(:))+1;
        NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
        NewMat(indz) = 0;
        a = logical(sum(logical(NewMat)')');
        indc = find(a);
        Surftest.Is(indc) = 0;
        Surftest = Surf_Corr(Surftest);
        
        
        % Subdividing Interception mesh
        [vertex,faces] = perform_mesh_subdivision(Surftest.SurfData.vertices',Surftest.SurfData.faces',1);
        Surftsub.SurfData.faces = faces';
        Surftsub.SurfData.vertices = vertex';
        
        % Refilling Labelling
        locbottom = dsearchn(Surfout(k).SurfData.vertices, Surftsub.SurfData.vertices);
        Surftsub.Is = Surfout(k).Is(locbottom);
        
        % Distance Transform
        distances = Compute_Geodesic_Distance_Transform(Surftsub);
        Surftsub.Is = distances + eps;
        
        % Sulcal Topology
        [skelEdges, Graph] = Extracting_Sulcal_Topology_V2(Surftsub);
        Surfa{k} = Surftsub;
        smeanLines{k} = skelEdges;
        
    end
end
varargout{1} = Surfa;
varargout{2} = smeanLines;