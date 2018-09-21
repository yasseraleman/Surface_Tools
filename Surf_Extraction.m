function [varargout] = Surf_Extraction(varargin);
%
% Syntax :
% [Surfa] = Surf_Comp(AImages,StructsId,Space);
%
% Computes the surface extraction for an atlas file or a mask file. If the
% imput images are atlas files, it extracts the surfaces for the structures
% specified in the structure list(StList).If the imput images are binary or
% mask images it will perform the surface extraction for all this images.
%
% Input Parameters:
%   AImages     : Individual Atlas files.
%   StructsId   : Structures Ids; 
%    Space      : Surface Space ('vox' or 'mm')
% Output Parameters:
%  Surfa        : Cell Array with surfaces in matlab variables format.
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


%% ===================== Checking Input Parameter======================= %%
if nargin == 0
   error('One Input is mandatory');
   return; 
end
if nargin > 3
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
AImages = varargin{1};
if ischar(AImages) % Verifying if the image exists
    if exist(AImages,'file')
        try
            V = spm_vol(AImages);
            I = uint16(spm_read_vols(V));
            voxsize = sqrt(sum(V.mat(1:3,1:3).^2));
        catch
            error('Unrecognized Image format or compressed image.');
            return;
        end
    end
elseif length(size(AImages)) == 3
    I = AImages;
    voxsize = [1 1 1];
else
    error('Unrecognized Image format or compressed image.');
    return;
end
clear AImages;
if nargin == 1
    Space = 'vox';
    StructsId = 1;
end
if nargin == 2
    StructsId = varargin{2};
    Space = 'vox';
end
if nargin == 3
    StructsId = varargin{2};
    Space = varargin{3};
end

%% =================== End of Checking Input Parameters ================ %%



%% ========================Main program================================= %%
Surf = '';
if nargin == 1
    I = logical(I);
end
if nargin >= 2
    indss = ismember(I(:),StructsId);
    indt = find(indss == 0);
    I(indt) = 0;
end
sts = nonzeros(unique(I(:)));
I = imfill(I,'holes');
Nst = length(sts);
Surf = struct('Name','','Area','','SurfData','','Tri','','Orig','','Dim','','VoxSize','','Code','');
for j = 1:Nst
    dat = int16(I*0);
    ind = find(I == sts(j));
    dat(ind) = 1;
    dat = logical(dat);
    
    
    bw = bwlabeln(double(dat));
    ind1 = find(dat);
    c = accumarray(bw(ind1),ones(length(bw(ind1)),1));
    indpos = find(c == max(c));
    dat(bw~= indpos) = 0;clear bw;
    It = logical(zeros(size(I,1)+4,size(I,2)+4,size(I,3)+4));It(3:size(I,1)+2,3:size(I,2)+2,3:size(I,3)+2) = dat;
    
    for i = 1:4
        It = imdilate(It,strel(ones(3,3,3)));
    end
    for i = 1:4
        It = imerode(It,strel(ones(3,3,3)));
    end
    
    switch Space
        case 'mm'
            [fv] = Ver_Ext(V(1),2000,It);
            fv=smoothpatch(fv,0,2);
        case 'vox'
            [fv] = Ver_Ext_vox(It,600);
            fv=smoothpatch(fv,0,4);
    end
    
    clear dat;
    disp('  ');
    disp('Computing Neighbor Points .....');
    Surf(j).SurfData.faces = fv.faces;
    Surf(j).SurfData.vertices = fv.vertices;
%    Surf(j).SurfData.VertexNormals = fv.normals;clear fv;
    Npoints = size(Surf(j).SurfData.vertices,1);
    Nfaces = size(Surf(j).SurfData.faces,1);
    tic;[Tri] = Vert_Neib(double(Surf(j).SurfData.faces),Npoints,Nfaces);toc;
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf(j).Tri = Tri;
    try [Surfa] = Surf_Ext_Corr(Surf(j)); Surf(j) = Surfa; end
    Surf(j).Is =  ones(size(Surf(j).SurfData.vertices,1),1)*double(sts(j));
    disp('   ');
    disp('Smoothing... ');
    %
%      try tic;[OutFiles,sSurf] = Smooth_surf(Surf(j),'',2,'n','');toc; end
%      Surf(j)= sSurf{1};
end
%% =======================End of main program=========================== %%

% Outputs
varargout{1} = Surf;
return;