function varargout = Compound_Surf(Surfa,varargin)
%
% Syntax :
% Surfout = Compound_Surf(Surfa);
%
% This function joints surfaces from different SurfFiles.
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

%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surfa = Surface_Checking(Surfa);
end

%% =================== End of checking input parameters ================= %

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
        if cont == 1
            Surfj.SurfData.vertices = Surft.SurfData.vertices;
            Surfj.SurfData.faces = Surft.SurfData.faces;
            % ================== Adding Scalar Values ======================= %
            if isfield(Surft, 'Is');
                Surfj.Is = Surft.Is;
            else
                Surfj.Is = zeros(size(Surft.SurfData.vertices,1),1);
            end
            % ================== Adding Colors Field ======================== %
            if isfield(Surf(i).SurfData, 'FaceVertexCData');
                Surfj.SurfData.FaceVertexCData = Surft.SurfData.FaceVertexCData;
            else
                Surfj.SurfData.FaceVertexCData = ones(size(Surft.SurfData.vertices,1),3);
            end
        else
            Surfj.SurfData.vertices = [Surfj.SurfData.vertices;Surft.SurfData.vertices];
            Surfj.SurfData.faces =    [Surfj.SurfData.faces;Surft.SurfData.faces+max(Surfj.SurfData.faces(:))];
            % ================== Adding Scalar Values ======================= %
            if isfield(Surft, 'Is');
                Surfj.Is = [Surfj.Is; Surft.Is];
            else
                Surfj.Is = [Surfj.Is; zeros(size(Surft.SurfData.vertices,1),1)];
            end
            % ================== Adding Colors Field ======================== %
            
            if isfield(Surf(i).SurfData, 'FaceVertexCData');
                Surfj.SurfData.FaceVertexCData = [Surfj.SurfData.FaceVertexCData; Surft.SurfData.FaceVertexCData];
            else
                Surfj.SurfData.FaceVertexCData = [Surfj.SurfData.FaceVertexCData; ones(size(Surft.SurfData.vertices,1),3)];
            end
        end
        
        
    end
    
end

%========================End of main program==============================%
% Outputs;
varargout{1} = Surfj;

return;