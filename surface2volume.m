function [OutFiles,sSurf] = Smooth_surf(SurfF,Out_dir,Iter,sa,verb);
%
% Syntax :
% [OutFiles,sSurf] = Smooth_surf(SurfF,Out_dir,Iter,sa,verb);
%
% This function computes a smoothing process to the surface contained in the Surf variable.
% The smoothing process will be applied the number of times given by the input variable Iter.
%
% Input Parameters:
%   SurfF        : A matlab structure, cellarray(1xNumber of Subjects) or filename
%                 containing the information about the surface that will be smoothed.
%   Out_dir      : Output directory for smoothed surfaces.
%   Iter         : Smooth Iterations. If Iter = 0, no smoothing.
%   verb         : verb = 'verbose' Show reduction progress, 0 do not show.
%    sa          : Variable to save or not the matlab surfaces.
%
% Output Parameters:
%   sSurf          : A matlab structure containing the smoothed surface.
%  OutFiles        : Smoothed surfaces filenames.
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
    Iter = input('Smooth Iterations:   ');
    if Iter<0
        errordlg('Please select a non negative numer of iterations')
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
if ~exist('Out_dir','var')
    [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
end
if ~exist('sa','var')|(isempty(sa))
   sa = 'y';
end
if ~exist('verb','var')
   verb = '';
end
if ~exist('Iter','var')|(isempty(Iter))
    Iter = input('Smooth Iterations:   ');
    if Iter<0
        errordlg('Please select a non negative numer of iterations')
        return
    end
end
if nargin==2
    verb = 0;
end
SurfFiles = '';
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


% parameters for the operator
laplacian_type = 'distance';
options.symmetrize = 1;
options.normalize = 1; % it must be normalized for filtering
options.verb = 0;


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
            if strcmp(fold,'smoothed_surf')
                Output_dir = Out_dir(t,1:ind(end)-1);
            else
                Output_dir = Out_dir(t,:);
            end
        elseif (size(ind)==1)&(ind(end)~=length((Out_dir(t,:))))
            fold = lower(Out_dir(t,ind(end-1)+1:ind(end)-1));
            if strcmp(fold,'smoothed_surf')
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
    if Iter == 0
        sSurf = Surf;
        return
    end
    for k = 1:length(Surf)
        if strcmp(verb,'verbose')
            disp(['Case ===> ' num2str(t) ' Smoothing ' Surf(k).Name]);
        end;
        W = compute_mesh_weight(Surf(k).SurfData.vertices,Surf(k).SurfData.faces,laplacian_type,options);
        % This is the corresponding laplacian
        %L = compute_mesh_laplacian(vertex,faces,laplacian_type,options);
        vertex2 = Surf(k).SurfData.vertices;
        % clf;
        options.face_vertex_color = [];
        for j=1:Iter
            vertex2 = (W*(W*vertex2));
        end
        Surf(k).SurfData.vertices = vertex2;
    end
    
    if strcmp(sa,'y')
        if exist('Out_dir','var')
            pth = Output_dir;
        end
        mkdir(pth,'Smoothed_Surf'); npth = [pth filesep 'Smoothed_Surf'];
        [ptht,nm,ext] = fileparts(SurfFiles(t,:));
        if strcmp(ext(1:4),'.mat')
            temp = strfind(nm,'Smooth_');
            if isempty(temp)
                nm = ['Smooth_' nm];
            end
            [Out] = savematfile([npth filesep nm '.mat'], Surf);
            OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
            sSurf = '';
        elseif strcmp(ext(1:4),'.mes')
            temp = strfind(Surf.Name,'Smooth_');
            if isempty(temp)
                Surf.Name = ['Smooth_' Surf.Name];
            end
            [OutFiles, SurfF] = Exp_Surf(Surf, '0', npth, 'mesh', 'exp','y');
            sSurf = '';
        elseif strcmp(deblank(ext),'.inflated')|strcmp(deblank(ext),'.pial')|strcmp(deblank(ext),'.sphere')|strcmp(deblank(ext),'.white')
            temp = strfind(Surf.Name,'Smooth_');
            if isempty(temp)
                Surf.Name = ['Smooth_' Surf.Name];
            end
            [OutFiles, SurfF] = Exp_Surf(Surf, '0', npth, deblank(ext(2:end)), 'exp','y');
            sSurf = '';
        else
            
            temp = strfind(Surf.Name,'Smooth_');
            if isempty(temp)
                Surf.Name = ['Smooth_' Surf.Name];
            end
            [OutFiles, SurfF] = Exp_Surf(Surf, '0', npth, ext(2:4), 'exp','y');
            sSurf = '';
        end
    elseif strcmp(sa,'n')
        OutFiles = '';
        sSurf{t,1} = Surf;
    end
end
%========================End of main program==============================%
return;

%========================Internal Functions===============================%
% Correcting NaN normal components
function [Surf] = Norn_Corr(Surf);
normals = Surf.SurfData.VertexNormals;
norma = sqrt(sum((normals').^2));
normals = normals./repmat(norma',[1 3]);
temp = sum(normals')';
ind = find(isnan(temp) ==1);
if ~isempty(ind)
    indt = ind;
    for i = 1:size(ind,1)
        vert = Surf.SurfData.faces(nonzeros(Surf.Tri(ind(i),3:end)));ind2 = ismember(vert,indt);vert(ind2) = [];
        tempN = sum([normals(vert,1) normals(vert,2) normals(vert,3)])/size(vert,1);
        normals(ind(i),:) = tempN;
        norma(ind(i)) = sqrt(sum((tempN').^2));
        indt(indt == ind(i)) = [];
    end
end
Surf.SurfData.VertexNormals = normals./repmat(norma',[1 3]);
return;