function varargout = Slice2Surf(bgImage,slices,varargin);
%
% Syntax :
% Surfout = Slice2Surf(bgImage,slices,colMap);
%
% This function creates a surface from a specified image slices.
%
% Input Parameters:
%   bgImage     : Image filename.
%   slices      : Slices (ie [0 0 128], [128 128 0])
%   colMap      : Colormap used to see the results.
%
% Output Parameters:
%   Surfout     : Surface of the slices.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% March 30th 2017
% Version $1.0

%% ====================== Checking input parameters ===================== %
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    % Reading ODF File or Variable
    temp = whos('bgImage');
    switch temp.class
        case 'char'
            if ~exist(bgImage,'file')
                error('The Image file does not exist');
                return
            else
                [pthbg,nmbg,extbg] = fileparts(bgImage);
                switch deblank(extbg)
                    case '.gz'
                        boolzipBG = 1;
                        tempBG = unzip_nifti(bgImage);
                    case '.nii'
                        boolzipBG = 0;
                        tempBG = bgImage;
                    otherwise
                        error('Unrecognized image format.');
                        return;
                end
                Vbg = spm_vol(tempBG);
                
                param.Xlim = [1 Vbg(1).dim(1)];
                param.Ylim = [1 Vbg(1).dim(2)];
                param.Zlim = [1 Vbg(1).dim(3)];
                voxsize = sqrt(sum(Vbg(1).mat(1:3,1:3).^2));
            end
        otherwise
            boolzipBG = 0;
    end
    
    % Parameters
    colMap = 'jet'; % ColorMap
        
end


% deal with the input arguments
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'colMap' % Colormap
                    colMap=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% =================== End of checking input parameters ================= %


%% ============================ Main Program ============================ %
cont = 0;
for t = 1:size(slices,1)
    slice = slices(t,:);
    ind = find(slice);
    for i = 1:length(ind)
        cont = cont + 1;
        
        switch temp.class
            case {'double','single'}
                
                if ind(i) == 3
                    Ibg = squeeze(bgImage(:,:,slice(3)));
                    sliceVal =  slice(3)-1;
                elseif ind(i) == 2
                    Ibg = squeeze(bgImage(:,slice(2),:));
                    sliceVal =  slice(2)-1;
                elseif ind(i) == 1
                    Ibg = squeeze(bgImage(slice(1),:,:));
                    sliceVal =  slice(1)-1;
                end
            otherwise
                %% ============ Detecting Slice ========================== %%
                if ind(i) == 3
                    DIM = Vbg(1).dim([1 2]);
                    C=[1 0 0 0;0 1 0 0;0 0 1 slice(3);0 0 0 1]; % Axial Transformation Matrix
                    Xlim = [param.Xlim(1):param.Xlim(2)]; % Axial Limits
                    Ylim = [param.Ylim(1):param.Ylim(2)];
                    sliceVal =  slice(3)-1;
                elseif ind(i) == 2
                    DIM = Vbg(1).dim([1 3]);
                    C=[1 0 0 0;0 0 1 slice(2);0 1 0 0;0 0 0 1]; % Coronal Transformation Matrix
                    Xlim = [param.Xlim(1):param.Xlim(2)]; % Coronal Limits
                    Ylim = [param.Zlim(1):param.Zlim(2)];
                    sliceVal =  slice(2)-1;
                elseif ind(i) == 1
                    DIM = Vbg(1).dim([2 3]);
                    C=[0 0 1 slice(1);1 0 0 0;0 1 0 0;0 0 0 1]; % Sagital Transformation Matrix
                    Xlim = [param.Ylim(1):param.Ylim(2)]; % Sagital Limits
                    Ylim = [param.Zlim(1):param.Zlim(2)];
                    sliceVal =  slice(1)-1;
                end
                
                %% ================== Creating Slice Surface ================== %%
                
                Ibg = spm_slice_vol(Vbg,C,DIM,0);
                Ibg = Ibg(Xlim,Ylim);
        end

        Ibg=rot90(Ibg);
        Ibg = flip(Ibg,1);
        dims = size(Ibg);
        [X,Y] = meshgrid([1:dims(1)],[1:dims(2)]);
        
        if ind(i) == 3
            Vertices = [X(:) Y(:) X(:)*0+sliceVal];
            Ibg = Ibg';

            Colors = Val2colors(Ibg(:),colMap);
        elseif ind(i) == 2
            Vertices = [X(:)*0+sliceVal X(:) Y(:)];
            Ibg = Ibg';
            Colors = Val2colors(Ibg(:),colMap);
        elseif ind(i) == 1
            Vertices = [X(:)  X(:)*0+sliceVal Y(:)];
            Ibg = Ibg';
            Colors = Val2colors(Ibg(:),colMap);
        end
        tri = delaunay(X,Y);
        Surfout(cont).SurfData.vertices = Vertices;
        Surfout(cont).SurfData.faces = tri;
        Surfout(cont).SurfData.FaceVertexCData = Colors;

    end
    %% ================ End of Creating Slice Surface ================== %%
    
end

cont = 0;
for i = 1:length(Surfout);
    cont = cont + 1;
    if cont == 1
        Surfj.SurfData.vertices = Surfout(i).SurfData.vertices;
        Surfj.SurfData.faces = Surfout(i).SurfData.faces;
        Surfj.SurfData.FaceVertexCData = Surfout(i).SurfData.FaceVertexCData;

    else
        Surfj.SurfData.vertices = [Surfj.SurfData.vertices;Surfout(i).SurfData.vertices];
        Surfj.SurfData.faces =    [Surfj.SurfData.faces;Surfout(i).SurfData.faces+max(Surfj.SurfData.faces(:))];
        Surfj.SurfData.FaceVertexCData = [Surfj.SurfData.FaceVertexCData; Surfout(i).SurfData.FaceVertexCData];
    end
end
clear Surfout;
Surfout = Surfj;
clear Surfj;

%% ====================== Compressing Images =========================== %%
if  (boolzipBG == 1)&(exist('Vbg','var'))
    zip_nifti(Vbg(1).fname);
    delete(Vbg(1).fname);
    remove_niimat(Vbg(1).fname);
end
%% ====================== End of Compressing Images ==================== %%

%% ======================== End of Main Program ============================ %

% Outputs
varargout{1} = Surfout;
return;
