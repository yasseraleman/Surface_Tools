function Plot_oversurf(SurfFile,CFiles,cl,tr,sp,sw);
% Syntax :
% Plot_oversurf(SurfFile,CFiles,cl,tr,sp,sw);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   SurfFile    : Surface files.
%   CFiles      : Grosor or atlas text files
%   cl          : Colormap used to see the results.
%   tr          : Transparency vector(values have to be between 0 and 1).
%   sp          : Boolean variable(sp = 1, plot results over the sphere,
%                 sa = 0, do not plot results over the sphere).
%   sw          : Boolean variable(sw = 1, plot surfaces in the same window,
%                 sw = 0, plot surfaces in diferent windows)
%
% Output Parameters:
%
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_Surf Atlas_Surf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0

%=====================Checking Input Parameters===========================%
if ~exist('cl','var')|(isempty(cl))
    cl = 'jet';
end
if ~exist('tr','var')|(isempty(tr))
    tr = 1;
end
if ~exist('sp','var')|(isempty(sp))
    sp = 0;
end
if ~exist('sw','var')|(isempty(sw))
    sw = 'y';
end
if sp == 0
    if ~exist('SurfFile','var')|(isempty(SurfFile))
        [SurfFile,sts] = spm_select([1 2],'any','Selecting Surface Files','',cd);
    end
    wh = whos('SurfFile');
    if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'));
        Surfa = SurfFile;
    elseif ischar(SurfFile);
        for i = 1:size(SurfFile,1)
            [pth nm ext] = fileparts(SurfFile(i,:));
            if ext(1:4) == '.mat';
                Surf = load('-mat',[pth filesep nm ext(1:4)]);
                if isfield(Surf,'Surf')
                    Surfa{i,1} = Surf.Surf;
                else
                    Surfa{i,1} = Surf;
                end
            else
                [OutFiles, Surf] = Exp_Surf(SurfFile(i,:), '0', '','', 'imp','n');
                Surfa{i,1} = Surf{1};
            end
        end
    end
elseif sp == 1
    temp = which('Sphere.mat');
    Surf = load('-mat',temp);
    if isfield(Surf,'Surf')
        Surfa{1,1} = Surf.Surf;
    else
        Surfa{1,1} = Surf;
    end
end
if ~exist('CFiles','var')|(isempty(CFiles))
    [CFiles,sts] = spm_select([1 2],'any','Selecting Thickness/Curvature Files','',cd);
end
mat = 0;
if sp == 1
    Surfa = repmat(Surfa,[size(CFiles,1) 1]);
end
Ns = size(Surfa,1);
Nc = size(CFiles,1);
if Ns~=Nc
    errordlg('Please Select as much text files as surfaces');
    return;
end
%=========================================================================%
%=========================Main Program====================================%
for i = 1:Ns
    Surf = Surfa{i,1};
    Npoints = size(Surf.SurfData.vertices,1);
    [pth,nm,ext] = fileparts(deblank(CFiles(i,:)));
    CFile = [pth filesep nm deblank(ext)];
    if isfield(Surf.SurfData,'FaceVertexCData')
        rmfield(Surf.SurfData,'FaceVertexCData');
    end
    [txt,ctab] = read_cfiles(CFile);
    %% ===================== Painting Surface ============================== %%
    if ctab.table == 0
        [Colors] = Surf_Color(Surf,cl);
        Surf.SurfData.FaceVertexCData = Colors;
    else
        sts = unique(txt);
        Nst = length(sts);
        Surf.SurfData.FaceVertexCData = zeros(size(Surf.SurfData.vertices,1),3);
        for j = 1:Nst
            ind = find(txt==sts(j));
            indc = find(ctab.table(:,5)==sts(j));
            if isempty(indc)
                Matname = 'unknown_structure';
                Color =   [1 1 1];      % Color
            else
                Matname = char(ctab.struct_names{indc});
                Color =  ctab.table(indc,1:3)/255;       % Color
            end
            Surf.SurfData.FaceVertexCData(ind,:) = repmat(Color,[length(ind) 1 ]);
        end
    end
    
    %% ===================== End of Painting Surface ======================= %%
    Surf.Is = txt;
%     [Colors] = Surf_Color(Surf,cl);
%     Surf.SurfData.FaceVertexCData = Colors;
    Surfa{i,1} = Surf;
end
Plot_Surf(Surfa,tr,sw,cl);
%========================End of main program==============================%
return;


function [curv, fnum] = read_char(fname);

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
    str = sprintf('could not open file %s.', fname) ;
    error(str) ;
end
% vnum = fread3(fid) ;
b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vals_per_vertex = fread(fid, 1, 'int32') ;
    curv = fread(fid, vnum, 'float') ;

    fclose(fid) ;
else
    b1 = fread(fid, 1, 'uchar') ;
    b2 = fread(fid, 1, 'uchar') ;
    b3 = fread(fid, 1, 'uchar') ;
    vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
    curv = fread(fid, vnum, 'int16') ./ 100 ;
    fclose(fid) ;
end

function [retval] = fread3(fid)

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;




