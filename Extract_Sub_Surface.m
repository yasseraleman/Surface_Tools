function varargout = Extract_Sub_Surface(varargin);
%
% Syntax :
%    [subSurf, subsurfLabels, subsurfIndexes] = Extract_Sub_Surface(Surf, stindex, opts);
%
% This scripts extracts sulcal skeleton from the curvature map computed for
% a surface.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%        stindex                        : Structures labels
%        opts                           : Options:
%                                         - opts.keep (Number of surface clusters (0: Keep all surface clusters))
%
% Output Parameters:
%        subSurf                        : Matlab surface variable.
%        subsurfLabels                  : Surface labels (vector). It contains the
%                                       same labels as stindex. If any
%                                       point in the surface contains a
%                                       specified label this value is
%                                       removed from the output labels 
%                                       vector.
%        subsurfIndexes                 : Sub surface indexes. Cell array
%                                       containing the indexes in the 
%                                       original surface for each extracted
%                                       subsurface.
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez and Javier Santoja
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 3
    opts.keep  = 0;            % Number of surface clusters (0: Keep all surface clusters)
else
    opts = varargin{3};
    if ~isfield(opts,'keep')
        opts.keep = 0;      % Number of surface clusters (0: Keep all surface clusters)
    end
end

Surf = varargin{1};

% Surface Checking
Surf = Surface_Checking(Surf);

if nargin >= 2
    stindex = varargin{2};
else
    stindex = unique(Surf.Is);
end

if nargin < 1
    error('One Inputs is mandatory');
    return
end

% Checking Options
if nargin > 3
    error('To Many Input Parameters');
    return;
end
if nargout > 3
    error('To Many Output Parameters');
    return;
end

%% ====================== End of Input parameters  =======================%

%% ============================= Main Program ========================== %%
cont = 0;
Surft = '';
subsurfLabels = '';
subsurfIndexes = '';
for i = 1:length(stindex)
    ind = find(Surf.Is == stindex(i)); % Points belonging to a specified region
    if ~isempty(ind)
        cont = cont +1;
        Surft = Surf;
        
        % Selecting only the faces in the sulci
        %indneg = find(subSurf.Is <0);
        a = ismember(Surft.SurfData.faces,ind);
        faces2rem = find(sum(a,2) ~= 3);
        Surft.SurfData.faces(faces2rem,:) = [];
        a = ismember(ind,Surft.SurfData.faces(:));
        ind(find(a == 0)) = [];
        
        % Detecting posible surface clusters.
        allEdges = [Surft.SurfData.faces(:,1) Surft.SurfData.faces(:,2); Surft.SurfData.faces(:,2) Surft.SurfData.faces(:,3);Surft.SurfData.faces(:,1) Surft.SurfData.faces(:,3)] ; % Edges from the intersection faces
        %       allEdges = allEdges(find(sum(ismember(allEdges,ind),2)==2),:);
        [~,b] = ismember(allEdges,ind);
        
        temp =  [ [b(:,1) ;b(:,2)] [b(:,2);b(:,1)] ];
        [C,iak,~] = unique(temp,'rows','last');
        Mat = sparse(C(:,1),C(:,2),1);
        tempMat = Label_Graph_Components(Mat);
        
        if opts.keep
            acc = accumarray(nonzeros(tempMat(:)),nonzeros(tempMat(:))*0+1);
            [reordacc, ord] = sort(acc,'descend');
            if length(ord) < opts.keep
                opts.keep = length(ord);
            end
            loc = ord(1:opts.keep);

            tempacc = [1:length(acc)]';
            tempacc(loc) = [];
            if ~isempty(tempacc)
                [X2rem,Y2rem] = find(ismember(tempMat,tempacc));
                ind2rem = unique([X2rem;Y2rem]);
                ind2rem = find(ismember(b,ind2rem));
                outPoints = allEdges(ind2rem);
                outPoints = unique(outPoints);
                a = ismember(ind,outPoints);
                ind(find(a == 1)) = [];
                
                a = ismember(Surft.SurfData.faces,ind);
                faces2rem = find(sum(a,2) ~= 3);
                Surft.SurfData.faces(faces2rem,:) = [];
            end
        end
        Surft = Reorg_Surf(Surft);
        subSurf(cont) = Surft;
        subsurfLabels(cont,1) = stindex(i);
        subsurfIndexes{cont} = ind;
    end
    
end
%% ====================== End of Main Program ========================== %%
% Outputs
if ~exist('subSurf','var')
    subSurf = Surf;
    subsurfIndexes{1} = find(Surf.SurfData.faces(:));
end

varargout{1} = subSurf;
varargout{2} = subsurfLabels;
varargout{3} = subsurfIndexes;

return;