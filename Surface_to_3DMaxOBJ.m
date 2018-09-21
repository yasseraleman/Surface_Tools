function varargout = Surface_to_3DMaxOBJ(varargin);
%
% Syntax :
% OBJFile = Surface_to_3DMaxOBJ(Surf, OBJFile, colors);
%
% This function exports Matlab Surfaces into WaveFront OBJ format.
% Materials file will be also created.
%
% Input Parameters:
%   Surf         : Surface Variable (Matlab Format)
%    OBJFile     : Output Object Filename
%   colortable   : Struct variable similar to freesurfer's colortable. It
%                  includes 2 fields:
%                  1. struct_names: Cellarray containing 
%                  structures names.
%                  2. table. Nstructures X 6 Matrix containing some
%                  material parameters (Colors, Transparency ( Tr), Ids)
%                  Table order (R G B 0 Id Tr)

%
% Output Parameters:
%  OBJFile      : Output Object Filename
%
% Related references:
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% September 12th 2012
% Version $1.0

%% ============================= Checking Inputs ======================= %%
if nargin < 2
    error('Two Inputs are mandatory');
    return
end
Surf   = varargin{1};
OBJFile   = varargin{2};
[pth,nm,ext] = fileparts(OBJFile);



% Checking Gyral Crown matrix
if ~exist(pth,'dir')
    error(['Unknown Directory: ' pth]);
    return;
end

if nargin  == 3
    colors = varargin{3}
end

if nargin > 3
    errordlg('To Many Input Parameters');
    return;
end
if nargout > 1
    errordlg('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%


%% ========================= Main Program ============================== %%
fid = fopen(OBJFile,'wt');
fprintf(fid,'%s\n','# Max2Obj Version 4.0 Mar 10th, 2001');
fprintf(fid,'%s\n','#');
ind = strfind(nm,'.');nm(ind) = [];
filename = [deblank(pth) filesep deblank(nm) '.mtl' ];
fidmlt = fopen(filename,'wt');
fprintf(fidmlt,'%s\n','# Max2Mtl Version 4.0 Mar 10th, 2001');
fprintf(fidmlt,'%s\n','#');
fprintf(fid, '%s\n',['mtllib ./' deblank(nm) '.mtl']);
if iscell(Surf); % Verifying if the surface is a matlab cell
    cont = 0;
    for i = 1:length(Surf)
        Surfttemp = Surf{i};
        % Surface Checking
        Surfttemp = Surface_Checking(Surfttemp);
        Ntemp = length(Surfttemp);
        for j = 1:Ntemp
            cont = cont + 1;
            Surfttempt.SurfData.vertices = Surfttemp(j).SurfData.vertices;
            Surfttempt.SurfData.faces = Surfttemp(j).SurfData.faces;
            Surfout(cont) = Surfttempt;
        end
    end
elseif isstruct(Surf); % Verifying if the surface is a matlab structure
    Surfout = Surface_Checking(Surf);
elseif ischar(Surf)
    cont = 0;
    for i = 1:size(Surf,1)
        % Surface Checking
        Surfttemp = Surface_Checking(deblank(Surf(i,:)));
        cont = cont + 1;
        Surfttempt.SurfData.vertices = Surfttemp.SurfData.vertices;
        Surfttempt.SurfData.faces = Surfttemp.SurfData.faces;
        Surfout(cont) = Surfttempt;
    end
end

if exist('Surfout','var')
    if isfield(Surfout,'SurfData')% Verifying important fields
        if ~isfield(Surfout(1).SurfData,'faces')|~isfield(Surfout(1).SurfData,'vertices'); % Verifying important fields
            error('Important fields from the Surface matlab are Missing');
            return
        else
            
            Ns = length(Surfout);
            if nargin < 3
                col = [240,163,255;0,117,220;153,63,0;76,0,92;0,92,49;43,206,72;...
                    255,204,153;148,255,181;143,124,0;157,204,0;194,0,136;...
                    0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
                    224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5]/255;
                
                Ncolor = size(col,1);
                re = floor(Ns/Ncolor); colors = repmat(col,[re+1 1]);
            end
            Nfaces = 0;
            for i = 1:Ns
                Surft = Surfout(i);
                Matname = sprintf('%.4d',i);
                
                Color = colors(i,1:3);
                disp(['Adding ' Matname ' Surface']);
               
                %% Saving Vertices
                fprintf(fid, '\n');
                fprintf(fid, '%s\n','g');
                fprintf(fid, '%s\n',['# object ' Matname ' to come ...']);
                fprintf(fid, '%s\n','#');
                Mat = Surft.SurfData.vertices';
                fprintf(fid,'v %.6f %.6f %.6f\n', Mat(:));
                fprintf(fid, '%s\n',['# ' num2str(size(Surft.SurfData.vertices,1)) ' vertices']);
                fprintf(fid, '\n');
                fprintf(fid, '%s\n',['g ' Matname]);
                %% Saving Structure Faces
                
                Mat = Surft.SurfData.faces'+Nfaces;
                Newfaces = size(Surft.SurfData.vertices,1);
                Nfaces = Nfaces+max(Surft.SurfData.faces(:));
                fprintf(fid, '%s\n',['usemtl ' Matname]);
                fprintf(fid,'%s\n', ['s 2']);
                fprintf(fid,'f %u %u %u\n', Mat(:));
                fprintf(fid, '%s\n',['# ' num2str(size(Surft.SurfData.faces,1)) ' faces']);
                % -----------------------------------------
                %% Creating Material Structure
                fprintf(fidmlt,'%s\n',['newmtl ' Matname]);
                fprintf(fidmlt,'Ka  %.1f %.1f %.1f\n',Color(1),Color(2),Color(3));
                fprintf(fidmlt,'Kd  %.1f %.1f %.1f\n',Color(1),Color(2),Color(3));
                fprintf(fidmlt,'Ks  %.1f %.1f %.1f\n',1,1,1);
                fprintf(fidmlt,'%s\n',['d ' sprintf('%.1f',1)]);
                fprintf(fidmlt,'%s\n',['Ns 8.0']);
                fprintf(fidmlt,'%s\n',['illum 2']);
                fprintf(fidmlt,'%s\n',['#']);
                % ----------------------------------------
            end
            fprintf(fid, '\n');
            fprintf(fid, '%s','g');
            fprintf(fidmlt,'%s\n','# EOF');
            fclose(fid);
            fclose(fidmlt);
        end
    else
        error('Unrecognized Surface format');
        return;
    end
end
%% ===================== End of Main Program =========================== %%
%  Outputs 
varargout{1} = OBJFile;

return

function [Surft] = Reorg_Surf(Surf);
%This scripts reorganize Surfaces in case of deleted faces. 

indf = unique(Surf.SurfData.faces(:));
%%%%%%%%%%%%%%%%%%%%%%%Corriegiendo Superficies para quedarme solo
%%%%%%%%%%%%%%%%%%%%%%%con la parte que me interesa
Surft = Surf;
Surft.SurfData.vertices = zeros(size(indf,1),3);
Surft.Is = zeros(size(indf,1),1);
Surft.SurfData.VertexNormals = zeros(size(indf,1),3);
Surft.SurfData.FaceVertexCData = zeros(size(indf,1),3);
Surft.SurfData.faces = 0*Surft.SurfData.faces;
for i =1:size(indf,1)
    Surft.SurfData.vertices(i,:) = Surf.SurfData.vertices(indf(i),:);
    if isfield(Surf.SurfData,'VertexNormals')
        Surft.SurfData.VertexNormals(i,:) = Surf.SurfData.VertexNormals(indf(i),:);
    end
    if isfield(Surf.SurfData,'FaceVertexCData')
        Surft.SurfData.FaceVertexCData(i,:) = Surf.SurfData.FaceVertexCData(indf(i),:);
    end
    if isfield(Surf.SurfData,'Is')
        Surft.Is(i,:) = Surf.Is(indf(i),:);
    end
    indn = find(Surf.SurfData.faces ==indf(i));Surft.SurfData.faces(indn) = i;
end