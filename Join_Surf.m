function [OutputFiles, Surfj] = Join_Surf(SurfFiles,Out_File,formats,sa,pl,verb)
%
% Syntax :
% [OutputFiles Surfj] = Join_Surf(SurfFiles,Out_File,sa);
%
% This function joints surfaces from different SurfFiles.
%
% Input Parameters:
%   SurfFiles   : Individual Surfaces.
%   Out_File    : Output File.
%   formats     : Diferent Surface Formats to save the resulting surface
%   sa          : Boolean variable to save or not the resulting surface (y/n).
%   pl          : Boolean variable to plot or not the resulting surface (y/n).
%   surface.
%   verb         : verb = 'verbose' Show status progress, '' do not show.
% Output Parameters:
%   OutputFiles  : Atlased Surfaces Files.
%   Surfj            : Matlab variable containing the resulting surface.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin==0
    [SurfFiles,sts] = spm_select([1 Inf],'any','Selecting Surfaces Files','',cd);
    contt = 0;ind = '';
    warning on;
    for i = 1:size(SurfFiles,1)
        [pth nm ext] = fileparts(SurfFiles(i,:));ext = lower(ext);
        if ~strcmp(ext(1:4),'.mat')&~strcmp(ext(1:4),'.obj')&~strcmp(ext(1:4),'.dfs')&~strcmp(ext(1:4),'.mes')&~strcmp(ext(1:4),'.off')&~strcmp(ext(1:4),'.srx')&~strcmp(ext(1:4),'.txt')~strcmp(ext(1:4),'.asc')
            warning(['The ' ext(1:4) ' file format is not suported. Please import the surface ' [nm ext(1:4)] ' file to our supported formats']);
            contt = contt + 1;
            ind(contt) = i;
        end
    end
    warning off;
    if ~isempty(ind)
        SurfF(ind,:) = [];
    end
    sa = input('Do you want to save the surfaces(y/n):   ','s');
    if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
    if strcmp(lower(sa),'y')
        [FileName, Out_dir] = uiputfile('Saving Resulting Surface...');
        Out_File = [Out_dir 'Joint_' FileName];
    else
        Out_File = '';
    end
    pl = 'y';
    verb='verbose';
else
    if ~exist('SurfFiles','var')|isempty(SurfFiles)
        [SurfFiles,sts] = spm_select([1 Inf],'mat','Selecting Surfaces Files','',cd);
        contt = 0;ind = '';
        warning on;
        for i = 1:size(SurfFiles,1)
            [pth nm ext] = fileparts(SurfFiles(i,:));ext = lower(ext);
            if ~strcmp(ext(1:4),'.mat')&~strcmp(ext(1:4),'.obj')&~strcmp(ext(1:4),'.dfs')&~strcmp(ext(1:4),'.mes')&~strcmp(ext(1:4),'.off')&~strcmp(ext(1:4),'.srx')&~strcmp(ext(1:4),'.txt')
                warning(['The ' ext(1:4) ' file format is not suported. Please import the surface ' [nm ext(1:4)] ' file to our supported formats']);
                contt = contt + 1;
                ind(contt) = i;
            end
        end
        warning off;
        if ~isempty(ind)
            SurfF(ind,:) = [];
        end
    end;
    if ~exist('verb','var')
        verb = '';
    end
    if ~exist('sa','var')|isempty(sa)
        sa = input('Do you want to save the surfaces(y/n):   ','s');
        if (lower(sa)~='y')&(lower(sa)~='n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
        if strcmp(lower(sa),'y')
            [FileName, Out_dir] = uiputfile('Saving Resulting Surface...');
            Out_File = [Out_dir 'Joint_' FileName];
        else
            Out_File = '';
        end
    end
end
warning off;
%=========================================================================%
%
%=========================Main program=====================================
wh = whos('SurfFiles');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = size(SurfFiles,1);
    SurfF = SurfFiles;
elseif ischar(SurfFiles(1,:));
    Ns = size(SurfFiles,1);
end
wh = whos('SurfF');
OutputFiles = '';
Surfj.Area = 0;
for i =1:Ns
    [OutputFiles, SurfF] = Exp_Surf(SurfFiles(i,:), '0', '','', 'imp','n','joining');
    Surf = SurfF{1};
    N = size(Surf,1);
    for j = 1:N
        if strcmp(verb,'verbose')
            disp(['Case ===> ' num2str(i) ' Adding ' Surf(j).Name]);
        end;
        %
        if j*i>=2
            Surfj.SurfData.vertices = [Surfj.SurfData.vertices; Surf(j).SurfData.vertices];
            Surfj.SurfData.faces = [Surfj.SurfData.faces; Surf(j).SurfData.faces+max(Surfj.SurfData.faces(:))];
            Surfj.SurfData.VertexNormals = [Surfj.SurfData.VertexNormals; Surf(j).SurfData.VertexNormals];
            if isfield(Surf(j), 'Is');
                Surfj.Is = [Surfj.Is; Surf(j).Is];
            else
                Surfj.Is = [Surfj.Is; zeros(size(Surf(j).SurfData.vertices,1),1)];
            end
        else
            Surfj.Name = ['Joint_' Surf(1).Name];
            Surfj.SurfData.vertices =Surf(j).SurfData.vertices;
            Surfj.SurfData.faces = Surf(j).SurfData.faces;
            Surfj.SurfData.VertexNormals = Surf(j).SurfData.VertexNormals;
            if isfield(Surf(j), 'Is')
                Surfj.Is =Surf(j).Is;
            else
                Surfj.Is = zeros(size(Surf(j).SurfData.vertices,1),1);
            end
        end
        Surfj.Area = Surfj.Area+ Surf(j).Area;
    end
end
[Tri] = Vert_Neib(double(Surfj.SurfData.faces),size(Surfj.SurfData.vertices,1),size(Surfj.SurfData.faces,1));
Temp = sum(Tri);
Tri(:,Temp==0) = [];
Surfj.Tri = Tri;
Surfj.Orig = Surf.Orig;
Surfj.VoxSize = Surf.VoxSize;
Surfj.Dim = Surf.Dim;
Surfj.Code = '';
Surfj.Type = 'Joint';
Surf = Surfj;
Surf.Imp='Joint';
if strcmp(pl,'y') 
    Plot_Surf(Surfj,1,'y');
end
%=====Saving Resulting Surface ============
Osys=computer;
if strcmp(sa,'y')
    [Out_dir,nm,ext] = fileparts(Out_File(1,:));
    if ~isempty(Out_dir)
        ind = strfind(Out_dir, filesep);
        if ind(end)~=length((Out_dir))
            fold = lower(Out_dir(ind(end)+1:end));
            if strcmp(fold,'joint_surf')
                Output_dir = deblank(Out_dir(1:ind(end)-1));
            else
                Output_dir = deblank(Out_dir);
            end
        elseif (size(ind)==1)&(ind(end)~=length((Out_dir)))
            fold = lower(Out_dir(ind(end-1)+1:ind(end)-1));
            if strcmp(fold,'joint_surf')
                Output_dir = deblank(Out_dir(1:ind(end-1)-1));
            else
                Output_dir = deblank(Out_dir);
            end
        else
            Output_dir = deblank(Out_dir);
        end
        ind = strfind(Output_dir, filesep);
        if ind(end)==size(Output_dir,2)
            Output_dir = Output_dir(1,1:ind(end)-1);
        end
    end
    mkdir(Output_dir,'Joint_Surf'); npth = [Output_dir filesep 'Joint_Surf'];
    temp = strfind(nm,'Joint_');
    if isempty(temp)
        nm = ['Joint_' nm];
    end
    for k = 1:size(formats,1)
        form = deblank(formats(k,:));
        disp(['==>  Exporting Surface ' deblank(nm) ' to .' form ' format <==']);
        [OutFiles, Surfa] = Exp_Surf(Surf, '0', npth, form, 'exp', 'y','joining');
        OutputFiles = strvcat(OutputFiles,OutFiles);
    end
    Surfj = '';
elseif strcmp(sa,'n')
    OutputFiles = '';
end
%========================End of main program==============================%
return;
