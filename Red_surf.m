function [OutFiles,rSurf] = Red_Surf(SurfF,Out_dir,Npoints,sa,verb);
%
% Syntax :
% [OutFiles,rSurf] = Red_Surf(SurfF,Out_dir,Npoints,sa,verb);
%
% This function computes a reduction process to the surface contained in the Surf variable.
% The resulting surface will have the number of points specified by the
% NPoints variable.
%
% Input Parameters:
%   SurfF        : A matlab structure, cellarray(1xNumber of Subjects) or filename
%                 containing the information about the surface that will be reduced.
%   Out_dir      : Output directory for reduced surfaces.
%   NPoints     : Number of points in the resulting surface, if the number of
%                 points given by NPoints is smaller than the number of
%                 points.
%   verb         : verb = 'verbose' Show reduction progress, 0 do not show.
%    sa          : Variable to save or not the matlab surfaces.
%
% Output Parameters:
%   rSurf          : A matlab structure containing the reduced surface.
%  OutFiles        : Reduced surfaces filenames.
%
% Related references:
%
%
% See also: Red_Surf Surf_Comp Plot_Surf Plot_oversurf Exp_Surf
% Atlas_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 3st 2006
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin==0
    [SurfF,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
    Npoints = input('Please enter the number of points:   ');
    if Npoints<0
        errordlg('Please select a non negative numer of points');
        return
    end
    sa = input('Do you want to save the surfaces in a new folder(y/n):   ','s');
    if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
    if strcmp(sa,'y')
        [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
    elseif strcmp(sa,'n')
        Out_dir = '';
    end
    verb = 'verbose';
end
if ~exist('SurfF','var')|(isempty(SurfF))
    [SurfF,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
end
if ~exist('sa','var')|isempty(sa)
        sa = input('Do you want to save the surfaces(y/n):   ','s');
        if (lower(sa)~='y')&(lower(sa)~='n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
        if strcmp(lower(sa),'y')
            [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
        else
            Out_dir = '';
        end
    end
if ~exist('verb','var')
   verb = '';
end
if ~exist('Npoints','var')|(isempty(Npoints))
    Npoints = input('Please enter the number of points:   ');
    if Npoints<0
        errordlg('Please select a non negative numer of points');
        return
    end
end
OutFiles = '';
%=========================================================================%
%=======================Main Program======================================%
wh = whos('SurfF');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = size(SurfF,1);
elseif ischar(SurfF(1,:));
    Ns = size(SurfF,1);
    SurfFiles = SurfF;
    [OutFiles, SurfF] = Exp_Surf(SurfF, '0', '','', 'imp','n');
end
if (size(Out_dir,1)~=Ns)&(size(Out_dir,1)==1)
        Out_dir = repmat(Out_dir,[Ns 1]);
end
wh = whos('SurfF');
for t = 1:Ns
    if strcmp(wh.class,'struct')
        Surf = SurfF;
    elseif strcmp(wh.class,'cell')
        Surf = SurfF{t,1};
    end
    if ~isempty(Out_dir)
        ind = strfind(Out_dir(t,:), filesep);
        if ind(end)~=length((Out_dir(t,:)))
            fold = lower(Out_dir(t,ind(end)+1:end));
            if strcmp(fold,'reduced_surf')
                Output_dir = Out_dir(t,1:ind(end)-1);
            else
                Output_dir = Out_dir(t,:);
            end
        elseif (size(ind)==1)&(ind(end)~=length((Out_dir(t,:))))
            fold = lower(Out_dir(t,ind(end-1)+1:ind(end)-1));
            if strcmp(fold,'reduced_surf')
                Output_dir = Out_dir(t,1:ind(end-1)-1);
            else
                Output_dir = Out_dir(t,:);
            end
        else
            Output_dir = Out_dir(t,:);
        end
        ind = strfind(Output_dir, filesep);
        if ind(end)==size(Output_dir,2)
            Output_dir = Output_dir(1,1:ind(end)-1);
        end
    end
    if Npoints == 0
        rSurf = Surf;
        return
    end
    for k = 1:length(Surf)
        if strcmp(verb,'verbose')
            disp(['Case ===> ' num2str(t) ' Reducing ' Surf(k).Name]);
        end;
        if (Npoints<=size(Surf(k).SurfData.vertices,1))&(Npoints ~=0)
            fv.vertices = double(Surf(k).SurfData.vertices);
            fv.faces = Surf(k).SurfData.faces;
            factor = 1/(size(fv.vertices,1)/Npoints);
            fv=reducepatch(fv,factor);
            Nv = size(fv.vertices,1);
            Nf = size(fv.faces,1);
            [Tri] = Vert_Neib(double(fv.faces),Nv,Nf);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf(k).SurfData.vertices = fv.vertices;
            Surf(k).SurfData.faces = fv.faces;
            h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
            norma = sqrt(sum((normals').^2));
            Surf(k).SurfData.VertexNormals = normals./repmat(norma',[1 3]);
            Surf(k).Tri = Tri; clear fv;
        elseif (Npoints>size(Surf(k).SurfData.vertices,1))&(Npoints ~=0)
            disp('The number of points is greater than the number of vertex in the surface. No reduction');
        end
        try [Surfa] = Surf_Ext_Corr(Surf(k)); Surf(k) = Surfa; end
    end
    if strcmp(sa,'y')
        if exist('Out_dir','var')
            pth = Output_dir;
        end
        mkdir(pth,'Reduced_Surf'); npth = [pth filesep 'Reduced_Surf'];
        [ptht,nm,ext] = fileparts(SurfFiles(t,:));
        if strcmp(ext(1:4),'.mat')
            temp = strfind(nm,'Red_');
            if isempty(temp)
                nm = ['Red_' nm];
            end
            [Outfile] = savematfile([npth filesep nm '.mat'], Surf);
            OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
            rSurf = '';
        elseif strcmp(ext(1:4),'.mes')|strcmp(ext(1:4),'.frb')|strcmp(ext(1:4),'.obj')|strcmp(ext(1:4),'.dfs')|strcmp(ext(1:4),'.asc')...
                |strcmp(ext(1:4),'.srx')|strcmp(ext(1:4),'.txt')|strcmp(ext(1:4),'.off')|strcmp(ext(1:4),'.vtk')
            temp = strfind(Surf.Name,'Red_');
            if isempty(temp)
                Surf.Name = ['Red_' Surf.Name];
            end
            [OutputFiles, SurfF] = Exp_Surf(Surf, '0', npth, ext(2:4), 'exp','y');
            OutFiles = strvcat(OutFiles,OutputFiles);
            rSurf = '';
        else
            temp = strfind(Surf.Name,'Red_');
            if isempty(temp)
                Surf.Name = ['Red_' Surf.Name];
            end
            [OutputFiles, SurfF] = Exp_Surf(Surf, '0', npth, 'frb', 'exp','y');
            OutFiles = strvcat(OutFiles,OutputFiles);
            rSurf = '';
        end
    elseif strcmp(sa,'n')
        OutFiles = '';
        rSurf{t,1} = Surf;
    end
end
%========================End of main program==============================%
return;








