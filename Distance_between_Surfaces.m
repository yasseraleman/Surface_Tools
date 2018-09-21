function varargout = Distance_between_Surfaces(varargin);
%
% Syntax :
%     [distMap] = Distance_between_Surfaces(Surf, Surfref, indexes, opts);
%
% This function computes the distance between two surfaces. The second
% input surface will be taken as the reference surface.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%        Surfref                        : Reference Surface variable (file, struct or cellarray).
%        indexes                        : The distance between surfaces
%                                         will be computed only for the
%                                         selected indexes.
%        opts (optional input)          : Options:
%                                         opts.metric - Distance metric (closest or normal)
%                                               closest - locate the
%                                               closest point in the
%                                               reference surface measuring
%                                               the Euclidean distance
%                                               between points.
%                                               normal  - intercept point
%                                               normal with the reference
%                                               surface.
%
% Output Parameters:
%         distMap                       : Distance Map.
%
%
% See also: Surface_Checking 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0


%% ============================= Checking Inputs ======================= %%
if nargin < 2
    error('The program needs two surface to compute the distance between them');
    return
end

Surf = varargin{1};
% Surface Checking
Surf = Surface_Checking(Surf);

Surfref = varargin{2};
% Surface Checking
Surfref = Surface_Checking(Surfref);

if nargin >= 3
    indexes = varargin{3};
end

if nargin < 4
    opts.metric  = 'closest';       % Distance metric
elseif nargin == 4
    opts = varargin{4};
    if isstruct(opts)
        if ~isfield(opts,'metric')
            opts.metric = 'closest';        % Distance metric
        end
    else
        error('Unrecognized options');
        return;
    end
end
% if nargin > 4
%     error('To Many Input Parameters');
%     return;
% end
% if nargout > 1
%     error('To Many Output Parameters');
%     return;
% end
%% ========================= End of Checking Inputs ==================== %%

%% ============== Computing Distance between surfaces  ================== %
disp(' ');
disp('Computing Distance between surface and its outher surface ... ');
tic;


VertP = Surfref.SurfData.vertices;FacesP = Surfref.SurfData.faces;NormalsP = Surfref.SurfData.VertexNormals;
Verts = Surf.SurfData.vertices;Normals = Surf.SurfData.VertexNormals;
if nargin >= 3
    Verts = Verts(indexes,:);
    Normals = Normals(indexes,:);
end

Npoints = size(Verts,1);
distMap = zeros(size(Verts,1),1);

% Computing Distances
switch opts.metric
    case 'normal' % intercept point normal with the reference surface.
%
        np = cross(VertP(FacesP(:,1),:)-VertP(FacesP(:,2),:),VertP(FacesP(:,3),:)-VertP(FacesP(:,2),:),2); Dp = -1*dot(np,VertP(FacesP(:,1),:),2);
        P1 = Verts+100*Normals; P2 = Verts-100*Normals; % Creating perpendicular lines
        
        
        for po =1:Npoints
            Num = dot(np,repmat(P1(po,:),[size(np,1) 1]),2)+Dp; Den = dot(np,repmat(P2(po,:)-P1(po,:),[size(np,1) 1]),2); % Creating Lines
            t = -1*Num./Den;clear Num Den; % Lines parameters
            xint = single(repmat(P1(po,1)',[size(FacesP,1) 1])+t.*(repmat((P2(po,1)-P1(po,1))',[size(FacesP,1) 1]))); % Line parametrization
            yint = single(repmat(P1(po,2)',[size(FacesP,1) 1])+t.*(repmat((P2(po,2)-P1(po,2))',[size(FacesP,1) 1])));
            zint = single(repmat(P1(po,3)',[size(FacesP,1) 1])+t.*(repmat((P2(po,3)-P1(po,3))',[size(FacesP,1) 1])));clear t;
            
            PpP1 =  VertP(FacesP(:,1),:)-[xint yint zint];
            PpP2 =  VertP(FacesP(:,2),:)-[xint yint zint];
            PpP3 =  VertP(FacesP(:,3),:)-[xint yint zint];
            angP2p3 = (acos(dot(PpP2,PpP3,2)./(normm(PpP2).*normm(PpP3))))*180/pi; % Angles between each face point and the interception point inpial surface
            angP1p3 = (acos(dot(PpP1,PpP3,2)./(normm(PpP1).*normm(PpP3))))*180/pi;
            angP1p2 = (acos(dot(PpP1,PpP2,2)./(normm(PpP1).*normm(PpP2))))*180/pi;
            
            ind = find(round(angP2p3+angP1p3+angP1p2) == 360);
            [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with the reference surface
            if ~isempty(inte)
                %     orientsign = sign((dot(inte -repmat(Surf2Process.SurfData.vertices(po,:),[size(inte,1) 1]),repmat(Surf2Process.SurfData.VertexNormals(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
                distTemp = sqrt(sum((inte - repmat(Verts(po,:),[size(inte,1) 1])).^2,2)); % Computing distance to the reference surface
                [valMin,posit] = min(distTemp);
                distMap(po) = valMin;
            end
        end
    case 'closest' % locate the closest point in the reference surface measuring the Euclidean distance between points.
        Nref = size(VertP,1);
        for po =1:Npoints
            distTemp = sqrt(sum((VertP - repmat(Verts(po,:),[Nref 1])).^2,2));
            [valMin,posit] = min(distTemp);
            distMap(po) = valMin;
        end
end
%
if nargin >= 3
    distMap = distMap(indexes);
end
%% ========== End of Computing Distance between surfaces  =============== %

% Outputs
varargout{1} = distMap;
return;