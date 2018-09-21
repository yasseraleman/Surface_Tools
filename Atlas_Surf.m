function [OutFiles,Surfa,models] = Atlas_Surf(SurfF,SVol,SMask,StrL,Out_dir,sa,ua,ma);
%
% Syntax :
% [OutFiles,Surfa] = Atlas_Surf(SurfF,SVol,SMask,StrL,Out_dir,sa,ua);
%
% This function says which points of the surface belongs to each anatomical
% structures taking into an account a given reference atlas.  The inter hemispheres
% points are labeled 300.
%
% Input Parameters:
%   SurfF     : Individual Surfaces.
%   SVol      : Reference Atlas
%  SMask    : White Matter Mask.
%  StrL        :  Cotains a list with the structures codes that will not
%                   be used in the labelling process.
%    - 'None'        : Use all the structures in the labelling process.
% i.e - '1,78-92,6' : The following structures will not be used. Structure
%                 with code equal to 1, structures with codes from 78 to 92 and
%                  structure with code equal to 6.
%                  It offers the chance to choose 2 or more structures with
%                  codes that aren't one after the other, just separating the
%                  structures codes with a coma(,).
%                  If you want to select a range of codes, you have to put the first
%                  structure code, then '-' and next the last structure code,
%                  i.e, 78-92. This expresion removes the structures with codes
%                  from 78 to 92 from the labelling process.
% Out_dir   : Output Directory.
%   sa         : Boolean variable to save or not the atlased surfaces (y/n).
%   ua         : Using a reference atlas(y/n)
% Output Parameters:
%   OutFiles  : Atlased Surfaces Files.
%   Surfa     : Matlab variable containing surface information.
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
% December 1st 2006
% Version $1.0

%%
%=====================Checking Input Parameters===========================%
if nargin==0
    [SurfF,sts] = spm_select([1 Inf],'any','Selecting Surfaces Files','',cd);
    contt = 0;ind = '';
    warning on;
    for i = 1:size(SurfF,1)
        [pth nm ext] = fileparts(SurfF(i,:));ext = lower(ext);
        if ~strcmp(ext(1:4),'.mat')&~strcmp(ext(1:4),'.obj')&~strcmp(ext(1:4),'.dfs')&~strcmp(ext(1:4),'.mes')&~strcmp(ext(1:4),'.off')&~strcmp(ext(1:4),'.srx')&~strcmp(ext(1:4),'.txt')&~strcmp(ext(1:4),'.asc')&~strcmp(ext(1:4),'.frb')
            warning(['The ' ext(1:4) ' file format is not suported. Please import the surface ' [nm ext(1:4)] ' file to our supported formats']);
            contt = contt + 1;
            ind(contt) = i;
        end
    end
    warning off;
    if ~isempty(ind)
        SurfF(ind,:) = [];
    end
    ua = input('Do you want to use a reference atlas?(y/n):   ','s');
    if strcmp(lower(ua),'y')&strcmp(lower(ua),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    elseif strcmp(lower(ua),'y')
        [SVol,sts] = spm_select([1],'image','Selecting Reference Atlas File','',cd);
        cdir = which('ibaspm'); [pth,nm,ext] = fileparts(cdir);
        SMask =[cdir filesep 'Hemi_Mask.img'] ;
        SVol = repmat(SVol, [size(SurfF,1) 1]);
        SMask = repmat(SMask, [size(SurfF,1) 1]);
    elseif strcmp(lower(ua),'n')
        [SVol,sts] = spm_select([1 Inf],'image','Selecting Individual Atlases','',cd);
        [SMask,sts] = spm_select([1 Inf],'image','Selecting White Matter Masks','',cd);
        if ~isempty(ind)
            SMask = SMask(~ind,:);
            SVol = SVol(~ind,:);
        end
    end
    ma = input('Do you want to create models matrix for BMA analysis?(y/n):   ','s');
    if strcmp(lower(ma),'y')&strcmp(lower(ma),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
    Na = size(SurfF,1);
    Ns = size(SVol);
    if Na ~= Ns
        errordlg('Please select the same number of files for surfaces and individual atlas files');
        return
    end
    sa = input('Do you want to save the surfaces(y/n):   ','s');
    if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
    StrL = input('Please enter the structures codes:  ','s');
    if strcmp(lower(sa),'y')
        [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
    else
        Out_dir = '';
    end

else
    if ~exist('SurfF','var')|isempty(SurfF)
        [SurfF,sts] = spm_select([1 Inf],'mat','Selecting Surfaces Files','',cd);
        contt = 0;ind = '';
        warning on;
        for i = 1:size(SurfF,1)
            [pth nm ext] = fileparts(SurfF(i,:));ext = lower(ext);
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
    if ~exist('ua','var')|isempty(ua)
        ua = input('Do you want to use a reference atlas?(y/n):   ','s');
        if (lower(ua)~='y')&(lower(ua)~='n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
    end
    if ~exist('ma','var')|isempty(ma)
        ma = input('Do you want to create models matrix for BMA analysis?(y/n):   ','s');
        if strcmp(lower(ma),'y')&strcmp(lower(ma),'n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
    end
    if ~exist('SVol','var')|isempty(SVol)
        if strcmp(lower(ua),'y')
            [SVol,sts] = spm_select([1],'image','Selecting Reference Atlas File','',cd);
            SMask =[cdir filesep 'Hemi_Mask.img'] ;
            SVol = repmat(SVol, [size(SurfF,1) 1]);
            SMask = repmat(SMask, [size(SurfF,1) 1]);
        elseif strcmp(lower(ua),'n')
            [SVol,sts] = spm_select([1 inf],'image','Selecting Individual Atlases','',cd);
        end
        Na = size(SurfF,1);
        Ns = size(SVol,1);
        if Na ~= Ns
            errordlg('Please select the same number of files for surfaces and individual atlas files or use a reference atlas');
            return
        end
    end
    if ~exist('SMask','var')|isempty(SMask)
        if strcmp(lower(ua),'n')
            [SMask,sts] = spm_select([1],'image','Selecting White Matter Masks','',cd);
        end
        Na = size(SurfF,1);
        Ns = size(SMask,1);
        if Na ~= Ns
            errordlg('Please select the same number of files for surfaces and individual white matter masks or use a reference atlas');
            return
        end
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
    if ~exist('StrL','var')|isempty(StrL)
        StrL = input('Please enter the structures codes:  ','s');
    end
end
if ~strcmp(lower(StrL),'none')
    StrList = Sel_Structs(StrL);
end
Na = size(SurfF,1);
Ns = size(SVol);
if Na ~= Ns
    errordlg('Please select the same number of files for surfaces and individual atlas files or use a reference atlas');
    return
end
warning off;
%=========================================================================%
%
%=========================Main program====================================%
Vt = spm_vol(SVol);
Vm = spm_vol(SMask);
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
Ns = size(SurfF,1);
wh = whos('SurfF');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = size(SurfF,1);
elseif ischar(SurfF(1,:));
    Ns = size(SurfF,1);
end
OutFiles = '';
Osys=computer;
for i = 1:Ns
    if ischar(SurfF(1,:));
        [OutputFiles, Surfa] = Exp_Surf(SurfF(i,:), SVol, '','', 'imp','n');
    end
    [pths,nms,exts] = fileparts(SurfF(i,:));
    [ptha,nma,exta] = fileparts(SVol(i,:));
    disp([' =======> Atlasing surface ' num2str(i) ' __ ' deblank(nms) deblank(exts) ' using '  deblank(nma) deblank(exta)  ' Atlas File <============']);
    It = spm_read_vols(Vt(i));
    satlas = length(unique(It(:)))-1;
    if ~strcmp(lower(StrL),'none')
        rinds = ismember(It(:),StrList');
        It(rinds) = 0;
    end
    Im = spm_read_vols(Vm(i));
    if strcmp(lower(ua),'y')
        ind = find((Im~=0)&(It == 0));
        Temp = zeros(size(Im));
        Temp(ind) = Im(ind);
        Temp((Temp ==16)|(Temp ==26)|(Temp ==36)|(Temp ==46)|(Temp ==56)) = 6;
        Temp((Temp ==13)|(Temp ==23)|(Temp ==33)|(Temp ==43)|(Temp ==53)) = 3;
        Temp((Temp ==7)|(Temp ==10)|(Temp ==1)|(Temp ==2)) = 0;
        ind = find(Temp~=0);
        [x,y,z] = ind2sub(size(Temp),ind);
        xt = x-1; indt = sub2ind(size(Temp),xt,y,z);
        it = find((Temp(indt)~=0)&(Temp(indt)~=Temp(ind)));
        [x,y,z] = ind2sub(size(Temp),ind(it));
        xmi = min(x);xmax = max(x);
        XTemp = zeros(size(Temp));
        XTemp(ind(it)) = 1;
        Nhood = ones(3,3,3);Nhood(1,:,:) =zeros(3,3);Nhood(3,:,:) =zeros(3,3);
        XTemp = imerode(XTemp,strel(Nhood));
        bw = bwlabeln(XTemp);
        ind = find(XTemp);
        c = accumarray(bw(ind),ones(length(bw(ind)),1));
        indpos = find(c == max(c));
        XTemp(bw~= indpos) = 0;
        Xtemp = imdilate(XTemp,strel(Nhood));
        for j = -5:5
            indt = sub2ind(size(XTemp),x-j,y,z);
            XTemp(indt) =1;
        end
        for j = 1:size(XTemp,1)
            if sum(sum(squeeze(XTemp(j,:,:)))) ~=0
                XTemp(j,:,:) = imfill(squeeze(XTemp(j,:,:)),'holes');
            end
        end
        Im = XTemp; clear XTemp Temp;
        bw = bwlabeln(Im);
        ind = find(Im);
        c = accumarray(bw(ind),ones(length(bw(ind)),1));
        indpos = find(c == max(c));
        Im(bw~= indpos) = 0;
        It(logical(Im)) = 300;
    elseif strcmp(lower(ua),'n')
        ind = find((Im~=0)&(It == 0));
        Temp = zeros(size(Im));
        Temp(ind) = Im(ind);
        Temp((Temp ==7)|(Temp ==10)|(Temp ==2)|(Temp ==1)) = 0;
        ind = find(Temp~=0);
        [x,y,z] = ind2sub(size(Temp),ind);
        xt = x-1; indt = sub2ind(size(Temp),xt,y,z);
        it = find((Temp(indt)~=0)&(Temp(indt)~=Temp(ind)));
        [x,y,z] = ind2sub(size(Temp),ind(it));
        xmi = min(x);xmax = max(x);
        XTemp = zeros(size(Temp));
        XTemp(ind(it)) = 1;
        for j = -5:5
            indt = sub2ind(size(XTemp),x+j,y,z);
            XTemp(indt) =1;
        end
        for j = 1:size(XTemp,1)
            if sum(sum(squeeze(XTemp(j,:,:)))) ~=0
                XTemp(j,:,:) = imfill(squeeze(XTemp(j,:,:)),'holes');
            end
        end
        Im = XTemp; clear XTemp Temp;
        bw = bwlabeln(Im);
        ind = find(Im);
        c = accumarray(bw(ind),ones(length(bw(ind)),1));
        indpos = find(c == max(c));
        Im(bw~= indpos) = 0;
        It(logical(Im)) = 300;
    end
    Surf = Surfa{1,1};
    if strcmp(wh.class,'struct')
        if strcmp(Surf(1).Type,'Atlas')
            errordlg('This variable contains an atlas structures surfaces');
            break ;
        end
    end
    if ~isempty(Out_dir)
        ind = strfind(Out_dir(i,:), filesep);
        if ind(end)~=length((Out_dir(i,:)))
            fold = lower(Out_dir(i,ind(end)+1:end));
            if strcmp(fold,'atlased_surf')
                Output_dir = Out_dir(i,1:ind(end)-1);
            else
                Output_dir = Out_dir(i,:);
            end
        elseif (size(ind)==1)&(ind(end)~=length((Out_dir(i,:))))
            fold = lower(Out_dir(i,ind(end-1)+1:ind(end)-1));
            if strcmp(fold,'atlased_surf')
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
    end
    if Vt(i).mat(1,1) <0
        Vt(i).mat(1,1) = -1*Vt(i).mat(1,1);
        Vt(i).mat(1,4) = -1*Vt(i).mat(1,4);
    end
    Surftemp=Surf.SurfData.vertices;
    Surf.SurfData.vertices = [Surf.SurfData.vertices ones(size(Surf.SurfData.vertices,1),1)]*inv(Vt(i).mat)';
    Surf.SurfData.vertices(:,4) = [];
    %Ism = spm_sample_vol(Im,double(Xi),double(Yi),double(Zi),0);%clear Im;
    It = imdilate(It,strel([1 1 1; 1 1 1; 1 1 1]));
    Is = spm_sample_vol(It,double(Surf.SurfData.vertices(:,1)),double(Surf.SurfData.vertices(:,2)),double(Surf.SurfData.vertices(:,3)),0);
    Surf.SurfData.vertices = Surftemp; clear Surftemp;
    Is = Recur(Is, Surf,3);
    ind = find(Is == 0);
    if ~isempty(ind)
        Is = Recur(Is, Surf,2);
    end
    ind = find(Is == 0);
    if ~isempty(ind)
        Is = Recur(Is, Surf,1);
    end
    uni = sort(unique(Is));
    re = floor(length(uni)/Ncolor); col = repmat(col,[re+1 1]);
    Surf.Is = int16(Is);
    [Surf] = Surf_Corr(Surf);
    [Surf] = Struct_Area_Comp(Surf);
    if strcmp(ma,'y')
       [models] = Create_BMAmodels(Surf,satlas);
    else
    end
    %%%%%%%%%%%%%%%%%%%
    if strcmp(sa,'y')
        if exist('Out_dir','var')
            pth = Output_dir;
        end
        mkdir(deblank(pth),'Atlased_Surf'); npth = [deblank(pth) filesep 'Atlased_Surf'];
        [ptht,nm,ext] = fileparts(SurfF(i,:));
        if strcmp(ext(1:4),'.mat')
            temp = strfind(nm,'Atlas_');
            if isempty(temp)
                nm = ['Atlas_' nm];
            end
            [Temp] = Save_Surf(Surf, SurfF(i,:));
            OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
            [out] = savematfile([npth filesep nm '.mat'], Surf);
            Surfa = '';
        else
            Tempext = '.txt';
            temp = strfind(nm,'Atlas_');
            if isempty(temp)
                nm = ['Atlas_' nm];
            end
            [Temp] = Save_Surf(Surf, SurfF(i,:));
            mkdir(deblank(pth),'Atlased_Surf'); npth = [deblank(pth) filesep 'Atlased_Surf'];
            OutFiles = strvcat(OutFiles,[npth filesep nm Tempext]);
            fid = fopen(deblank(OutFiles(i,:)),'w');
            fprintf(fid,'%4u\n',Surf.Is);
            fclose(fid);
            Surfa = '';
        end
        if strcmp(lower(ma),'y')
            [pthm,nmm,extm] = fileparts(deblank(OutFiles(i,:)));
            ModelFile = [pthm filesep nmm '_models.mat'];
            [ModelFile] = savematfile(ModelFile, models);
        else
            models = '';
        end
    elseif strcmp(sa,'n')
        OutFiles = '';
        Surfa{i,1} = Surf;
        if strcmp(lower(ma),'n')
            models = '';
        end
    end
end
%=========================End of main program=============================%
return



%==================Internal Functions======================================
function [models] = Create_BMAmodels(Surf,satlas)

[cd ,nm ]=fileparts(which('IBASPM'));
if satlas == 116
    AtlasFileCod = [cd filesep 'atlas116.cod'];
elseif abs(71-satlas) <=3 
    AtlasFileCod = [cd filesep 'atlas71.cod'];
elseif satlas == 16
    AtlasFileCod = [cd filesep 'Hemi_Mask.cod'];
elseif satlas == 19
    AtlasFileCod = [cd filesep 'Atlas116_Lobes.cod'];
else
    AtlasFileCod = '';
end
if ~isempty(AtlasFileCod)
    fid = fopen(AtlasFileCod);
    cont = 0;
    while 1
        cont = cont + 1;
        line = fgetl(fid);
        if ~ischar(line),   break,   end
        ind = strfind(line,'=');
        StructNames{1, cont} = line(ind+1:end);
        StructCodes(1,cont) = str2num(deblank(line(1:ind-1)));
    end
    fclose(fid);
    clear tline cont;
else
    uni = unique(Surf.Is);
    uni(uni == 0) = [];
    StructNames = num2cell(uni)';
    StructCodes = uni;
end
uni = unique(Surf.Is);
uni(uni==0)=[];
Nstruct = size(uni,1);
for j =1:Nstruct
    ind = find(Surf.Is == uni(j));
    indt = find(StructCodes ==uni(j));
    if uni(j)~=300
        if ~isempty(AtlasFileCod)
            models(j).name = StructNames{1,indt};
        else
        models(j).name = num2str(StructNames{1,indt});
        end
        models(j).indices = ind;
    else
        models(j).name = 'Corpus Callosum';
        models(j).indices = ind;
    end
end
models(1).Npoints = size(Surf.SurfData.vertices,1);
models(1).Normals = Surf.SurfData.VertexNormals;
return

%==========================================================================