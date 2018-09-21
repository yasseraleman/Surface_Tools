function [OutputFiles,Surfa,At] = Surf_Comp(AImages,CodeFile,Out_dir,StList,NPoints,pl,tr,Iter,comp_a,form, sa);
%
% Syntax :
% [At] =  Surf_Comp(MImages,AImages,CodeFile,Out_file,Out_dir,StList,NPoints,pl,ch,tr,Iter,comp_a);
%
% Computes the surface extraction for an atlas file or a mask file. If the
% imput images are atlas files, it extracts the surfaces for the structures
% specified in the structure list(StList).If the imput images are binary or
% mask images it will perform the surface extraction for all this images.
%
% Input Parameters:
%   AImages     : Individual Atlas files.
%   CodeFile    : Text File with the relationship between structures codes
%                 and structures names
%   Out_dir     : Output directory for surface files.
%                 If the user doesn't change the output directory, the
%                 resulting files are saved in the same address than the
%                individual Mask/Atlas files.
%   StList      : List of structures codes.
%   NPoints     : Number of vertex in the surface.
%   pl          : Boolean variable for plotting surfaces(0 do not plot surfaces, 1 plot surfaces).
%   tr          : Set the surface transparence.
%   Iter        : Number of time
%   comp_a      : Boolean variable(0,do not compute the surfaces area; 1,compute the surfaces areas).
%   form        : Format to export the computed surface(srx & txt: Imagic format,
%   obj: CLASP surface format, off: FSL Betall format, mesh: Brainvisa surface format,
%   dfs: BrainSuite surface format).
%   sa          : Variable to save or not the matlab surfaces.
%
% Output Parameters:
%  OutputFiles    : List with the computed surfaces.
%  Surfa        : Cell Array with surfaces in matlab variables format.
%   At          : N(subjets)xM(structures/mask) surfaces areas matrix.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2006
% Version $1.0
warning off
fclose('all');
%=====================Checking Input Parameters===========================%
if nargin ==0
    sa ='y';
    tr=1;pl =1;
    [AImages,sts] = spm_select([1 Inf],'image','Selecting Atlased Images','',cd);
    [CodeFile,sts] = spm_select([1],'cod','Reading Code File ...','',cd);
    StList = input('Please enter the structures codes:  ','s');
    [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
    NPoints = input('Please enter the number of points in the surface (ie: NPoints = 0, No Reduction ):  ');
    Iter = input('Please enter the smoothing iterations (ie: Iter = 0, No smoothing ):  ');
    ar = input('Compute Surface Area (ie: y/n ):  ','s');
    if strcmp(ar,'y')&strcmp(ar,'n')
        errordlg('Please answer y or n ');
    end
    if strcmp(ar,'y')
        comp_a = 1
    elseif strcmp(ar,'n')
        comp_a =0;
    end
    form = input('Please enter the format( ie:  off,obj,dfs,mesh,srx,non):   ','s');
    if ~strcmp(form,'off')&~strcmp(form,'obj')&~strcmp(form,'dfs')&~strcmp(form,'mesh')&~strcmp(form,'srx')&~strcmp(form,'txt')&~strcmp(form,'none')
        errordlg('Please enter a correct format( ie: off,obj,dfs,mesh,srx,txt,none) ');
        return;
    end;
end

%=========================================================================%
%=========================Main program====================================%
OutputFiles = '';
fid = fopen(CodeFile);
if fid >0
    cont = 0;
    StructNames = '';
    while 1
        cont = cont + 1;
        line = fgetl(fid);
        if ~ischar(line),   break,   end
        ind = strfind(line,'=');
        StructNames = strvcat(StructNames,line(ind+1:end));
        StructCodest(1,cont) = str2num(deblank(line(1:ind-1)));
    end
    fclose(fid);
    clear tline cont;

    if strcmp(lower(StList),'all')
        StrList = StructCodest;
    else
        StrList = Sel_Structs(StList,StructCodest);
    end
end
V = spm_vol(AImages);
Ns  = length(V);
if (size(Out_dir,1)~=Ns)&(size(Out_dir,1)==1)
    Out_dir = repmat(Out_dir,[Ns 1]);
end
%H = waitbar(0,['Creating Surface ' num2str(Ns)],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');
for i = 1:Ns
    [pth nm ext] = fileparts(V(i).fname);
    I = uint16(spm_read_vols(V(i))); I = imfill(I,'holes');
    if (fid<0)|((fid>0)&(sum(ismember(StrList,StructCodest)) ~= length(StrList)))
        StructCodes = unique(I(:)');StructCodes(StructCodes == 0) = [];
        StructNames = num2str(StructCodes');
        StrList = Sel_Structs(lower(StList),StructCodes);
    end
    if ((fid>0)&(sum(ismember(StrList,StructCodest)) == length(StrList))&(~isempty(StrList)))
        StructCodes = unique(I(:)');StructCodes(StructCodes == 0) = [];
        StrList = Sel_Structs(lower(StList),StructCodes);
        ind = ismember(StrList,StructCodes);StrList(ind == 0) = [];
        ind = ismember(StructCodest,StrList);
        StructCodes = StructCodest(ind);
        StructNamest = deblank(StructNames(ind,:));
    end
    if (length(unique(I(:)))-1 ==1)
        NStruct = 1;
        I = logical(I);
        Itype = 'Mask';
        StructNames = nm;
        StructCodes =1;
        StrList = 1;
    else
        NStruct = length(StrList);
        Itype = 'Atlas';
    end
    
    ind = strfind(Out_dir(i,:), filesep);
    if ind(end)~=length((Out_dir(i,:)))
        fold = lower(Out_dir(i,ind(end)+1:end));
        if strcmp(fold,'atlased')|strcmp(fold,'struct_surf')
            Output_dir = Out_dir(i,1:ind(end)-1);
        else
            Output_dir = Out_dir(i,:);
        end
    elseif (size(ind)==1)&(ind(end)~=length((Out_dir(i,:))))
        fold = lower(Out_dir(i,ind(end-1)+1:ind(end)-1));
        if strcmp(fold,'atlased')|strcmp(fold,'struct_surf')
            Output_dir = Out_dir(i,1:ind(end-1)-1);
        else
            Output_dir = Out_dir(i,:);
        end
    else
        Output_dir = Out_dir(i,:);
    end
    ind = strfind(Output_dir, filesep);
    if ind(end)==size(Output_dir,2)
        Output_dir = Output_dir(1,1:ind(end)-1);
    end
    voxsize = sqrt(sum(V(i).mat(1:3,1:3).^2));
    Surf = struct('Name','','Area','','SurfData','','Tri','','Orig','','Dim','','VoxSize','','Code','');
    for j = 1:NStruct
        disp('  ');
        inds = find(StructCodes==StrList(j));
        if length(StructCodes)~=1
            disp(['Case --> ' num2str(i) ' ====>> Creating Surface for ' num2str(StructCodes(1,inds)) '_' StructNamest(inds,:) ' <<====']);
        else
            disp(['Case --> ' num2str(i) ' ====>> Creating Surface for ' nm ' <<====']);
        end
        if NPoints ==0
            %                 if j==1
            %                     %waitbar(i/Ns,H,['Creating Surface ' num2str(i) ' of ' num2str(Ns) ': No Reduction']);
            %                 end
            disp(['No Reduction...']);
        else
            %                 if j==1
            %                     %waitbar(i/Ns,H,['Creating Surface ' num2str(i) ' of ' num2str(Ns) ': Number of Points == ' num2str(NPoints)]);
            %                 end
            disp(['The surface will be reduced to ' num2str(NPoints) ' points']);
        end
        if Iter ==0
            disp(['No Smoothing...']);
        else
            disp(['Smooth Iterations == ' num2str(Iter)]);
        end
        if comp_a ==0
            disp(['No Area Computation...']);
        end
        dat = logical(zeros(size(I)));
        ind  = find(I ==StrList(j));
        dat(ind) = 1;
        bw = bwlabeln(double(dat));
        ind1 = find(dat);
        c = accumarray(bw(ind1),ones(length(bw(ind)),1));
        indpos = find(c == max(c));
        dat(bw~= indpos) = 0;clear bw;
        It = logical(zeros(size(I,1)+4,size(I,2)+4,size(I,3)+4));It(3:size(I,1)+2,3:size(I,2)+2,3:size(I,3)+2) = dat;
        [fv] = Ver_Ext(V(i),NPoints,It);
        clear dat;
        disp('  ');
        disp('Computing Neighbor Points .....');
        Surf(j).SurfData.faces = fv.faces;
        Surf(j).SurfData.vertices = fv.vertices;
        Surf(j).SurfData.VertexNormals = fv.normals;clear fv;
        Npoints = size(Surf(j).SurfData.vertices,1);
        Nfaces = size(Surf(j).SurfData.faces,1);
        tic;[Tri] = Vert_Neib(double(Surf(j).SurfData.faces),Npoints,Nfaces);toc;
        Temp = sum(Tri);
        Tri(:,Temp==0) = [];
        Surf(j).Tri = Tri;
        if strcmp(Itype,'Atlas')
            Surf(j).Name = deblank([nm '_' num2str(StructCodes(inds)) '_' StructNamest(inds,:)]);
            Surf(j).Code = StrList(j);
        elseif strcmp(Itype,'Mask')
            Surf(j).Name = nm;
            Surf(j).Code = 1;
        end
        try [Surfa] = Surf_Ext_Corr(Surf(j)); Surf(j) = Surfa; end
        if comp_a ==1
            disp('   ');
            disp('Computing Surface Area...');
            tic;[At] = Area_Comp(Surf(j).SurfData);
            disp(['Surf_Area  ==  ' num2str(At) ' cm^2']);toc;
        else
            At = 0;
        end
        Surf(j).Area = At;
        Surf(j).Orig = abs(V(i).mat(1:3,4)');
        Surf(j).Dim = abs(V(i).dim(1:3));
        Surf(j).VoxSize = voxsize;
        Surf(j).Imp='Ext';
        Surf(1).Type = Itype;
        if Iter~=0
            disp('   ');
            disp('Smoothing... ');
            try tic;[OutFiles,sSurf] = Smooth_surf(Surf(j),'',Iter,'n','');toc; end
            Surf(j) = sSurf{1};
        end
    end
    clear dat Tri It;
    if strcmp(sa,'y')
        mkdir(Output_dir,'Struct_Surf'); npth = [Output_dir filesep 'Struct_Surf'];
        if length(Surf)==1
            OutputFiles = strvcat(OutputFiles,[npth filesep Surf.Name '_Struct_Surf.mat']);
        else
            OutputFiles = strvcat(OutputFiles,[npth filesep nm '_Struct_Surf.mat']);
        end
        [Outfile] = savematfile(deblank(OutputFiles(i,:)), Surf);
        Surfa = '';
    elseif strcmp(sa,'n')
        OutFiles = '';
        Surfa{i,1} = Surf;
    end
    if ~strcmp(form,'none')
        disp('   ');
        disp(['Exporting Surface to the specified formats... ']);
        tic;[OutF, Surfa] = Exp_Surf(deblank(OutputFiles(i,:)), V(i).fname, '',form, 'exp','y');toc;
    end
    clear Surf;

end
%close(H);
if pl==1
    if strcmp(sa,'y')
        Plot_Surf(OutputFiles(1,:),tr,sa);
    elseif strcmp(sa,'n')
        Plot_Surf(Surfa{1,1},tr,sa);
    end
end
%========================End of main program==============================%
return;

%========================Internal Functions===============================%

function [V, P] = crop_st(P,VA, borders);
S = size(P);
borders = 2*borders;
ind = find(P);
[x,y,z] = ind2sub(S,ind);
P = P(min(x):max(x),min(y):max(y),min(z):max(z));
S = size(P);
P(end + borders(1),end + borders(2),end + borders(3)) = 0;
d = round((size(P)/2) - (S/2));
P = translateImageN0(P,d(1),d(2),d(3));
[pathstr,name,ext,vers] = fileparts(VA.fname);
V = VA;
V.dim(1:3) = size(P);
% vec = [min(x)-1 min(y)-1 min(z)-1] - borders/2; vec = V.mat(1:3,1:3)*vec';
% V.mat(1:3,4) = V.mat(1:3,4) + vec;
return;

%==========================================================================