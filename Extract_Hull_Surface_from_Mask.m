function varargout = Extract_Hull_Surface_from_Mask(varargin);
%
% Syntax :
% maskImage = Extract_Hull_Surface_from_Mask(maskImage);
%
% This script a hull surface from a mask image.
%
% Input Parameters:
%       maskImage               : Mask image
%
% Output Parameters:
%       Surf                    : Hull Surface in matlab format.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

warning off;
if nargin == 0
    error('A mask image must be specified');
    return;
else
    maskImage = varargin{1};
    if ischar(maskImage);
        maskImage = spm_read_vols(spm_vol(maskImage));
    end
end
%% Mask Refilling
for z = 1:size(maskImage,3)
    Temp = squeeze(maskImage(:,:,z));
    %     Temp = imfill(Temp,'holes');
    Temp = bwmorph(Temp,'fill',Inf);
    maskImage(:,:,z) = Temp;
end;
for z = 1:size(maskImage,1)
    Temp = squeeze(maskImage(z,:,:));
    %     Temp = imfill(Temp,'holes');
    Temp = bwmorph(Temp,'fill',Inf);
    
    maskImage(z,:,:) = Temp;
end;
for z = 1:size(maskImage,2)
    Temp = squeeze(maskImage(:,z,:));
    %     Temp = imfill(Temp,'holes');
    Temp = bwmorph(Temp,'fill',Inf);
    
    maskImage(:,z,:) = Temp;
end;
bw = bwlabeln(maskImage);
ind = find(maskImage);
c = accumarray(bw(ind),ones(length(bw(ind)),1));
indpos = find(c == max(c));
maskImage(bw~= indpos) = 0;clear bw;

Iold = maskImage;
maskImage = imclose(Iold,strel(ones(3,3,3)));
while sum(Iold(:)-maskImage(:))~=0
    Iold = maskImage;
    for z = 1:size(Iold,3)
        Temp = squeeze(Iold(:,:,z));
        %         Temp = imfill(Temp,'holes');
        Temp = bwmorph(Temp,'fill',Inf);
        
        Iold(:,:,z) = Temp;
    end;
    for z = 1:size(Iold,1)
        Temp = squeeze(Iold(z,:,:));
        %         Temp = imfill(Temp,'holes');
        Temp = bwmorph(Temp,'fill',Inf);
        
        Iold(z,:,:) = Temp;
    end;
    for z = 1:size(Iold,2)
        Temp = squeeze(Iold(:,z,:));
        %         Temp = imfill(Temp,'holes');
        Temp = bwmorph(Temp,'fill',Inf);
        Iold(:,z,:) = Temp;
    end;
    bw = bwlabeln(Iold);
    ind = find(Iold);
    c = accumarray(bw(ind),ones(length(bw(ind)),1));
    indpos = find(c == max(c));
    Iold(bw~= indpos) = 0;clear bw;
    maskImage = imclose(Iold,strel(ones(3,3,3)));
end
maskImage = imfill(maskImage,'holes');

Surf = Surf_Extraction(maskImage);
Surf.Is = zeros(size(Surf.SurfData.vertices,1),1);

[labid] = Recur_Corr(Surf,0,zeros(size(Surf.Is)),1);
Surf.Is = labid;

index = accumarray(labid,labid*0+1);
[~,indmax]  = max(index);
ind2del = find(labid~=indmax);
face2del = find(sum(ismember(Surf.SurfData.faces,ind2del),2) ~= 0);
Surf.SurfData.faces(face2del,:) = [];
Surf = Reorg_Surf(Surf);

faces =   boundary(Surf.SurfData.vertices(:,1),Surf.SurfData.vertices(:,2),Surf.SurfData.vertices(:,3),.98);
Surf.SurfData.faces = faces;
Surf = Reorg_Surf(Surf);





% % % % % % % % disp(['Left Hemisphere: Refilling Sulcal Space ... Z']);
% % % % % % % % tic;
% % % % % % % % [X,Y] = meshgrid(1:size(maskImage,1),1:size(maskImage,2));X =X(:);Y =Y(:);
% % % % % % % % for i = 1:size(maskImage,3)
% % % % % % % %     A = maskImage(:,:,i);
% % % % % % % %     T = A*0;
% % % % % % % %     if sum(A(:))
% % % % % % % %         T = Hull_from_Mask_Slice(A,Y,X);
% % % % % % % %         % % %         [B,L] = bwboundaries(logical(A),'noholes');
% % % % % % % %         % % %         for j = 1:length(B)
% % % % % % % %         % % %             boundary = B{j};
% % % % % % % %         % % %             k=LineCurvature2D(boundary);
% % % % % % % %         % % %             ind = find(k >0);
% % % % % % % %         % % %             boundary(ind,:) = [];
% % % % % % % %         % % %             indold = 0;
% % % % % % % %         % % %             while (isequal(ind,indold) ==0)|~isempty(ind)
% % % % % % % %         % % %                 indold = ind;
% % % % % % % %         % % %                 k=LineCurvature2D(boundary);
% % % % % % % %         % % %                 indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
% % % % % % % %         % % %                 ind = find((k >0));
% % % % % % % %         % % %                 boundary(ind,:) = [];
% % % % % % % %         % % %             end
% % % % % % % %         % % %             boundary = [boundary;boundary(1,:)];
% % % % % % % %         % % %             indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
% % % % % % % %         % % %             indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
% % % % % % % %         % % %             T(indin) = j;
% % % % % % % %         % % %         end
% % % % % % % %     end
% % % % % % % %     Temp(:,:,i) = T;
% % % % % % % % end
% % % % % % % % 
% % % % % % % % disp(['Left Hemisphere: Refilling Sulcal Space ... Y']);
% % % % % % % % tic;
% % % % % % % % [X,Y] = meshgrid(1:size(Temp,1),1:size(Temp,3));X =X(:);Y =Y(:);
% % % % % % % % for i = 1:size(Temp,2)
% % % % % % % %     A = squeeze(Temp(:,i,:));
% % % % % % % %     T = A*0;
% % % % % % % %     if sum(A(:))
% % % % % % % %         T = Hull_from_Mask_Slice(A,Y,X);
% % % % % % % %         % % %         [B,L] = bwboundaries(logical(A),'noholes');
% % % % % % % %         % % %         for j = 1:length(B)
% % % % % % % %         % % %             boundary = B{j};
% % % % % % % %         % % %             k=LineCurvature2D(boundary);
% % % % % % % %         % % %             ind = find(k >0);
% % % % % % % %         % % %             boundary(ind,:) = [];
% % % % % % % %         % % %             indold = 0;
% % % % % % % %         % % %             while (isequal(ind,indold) ==0)|~isempty(ind)
% % % % % % % %         % % %                 indold = ind;
% % % % % % % %         % % %                 k=LineCurvature2D(boundary);
% % % % % % % %         % % %                 indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
% % % % % % % %         % % %                 ind = find(k >0);
% % % % % % % %         % % %                 boundary(ind,:) = [];
% % % % % % % %         % % %             end
% % % % % % % %         % % %             boundary = [boundary;boundary(1,:)];
% % % % % % % %         % % %             indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
% % % % % % % %         % % %             indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
% % % % % % % %         % % %             T(indin) = j;
% % % % % % % %         % % %         end
% % % % % % % %     end
% % % % % % % %     Temp(:,i,:) = T;
% % % % % % % % end
% % % % % % % % 
% % % % % % % % % % disp(['Left Hemisphere: Refilling Sulcal Space ... X']);
% % % % % % % % % % tic;
% % % % % % % % % % [X,Y] = meshgrid(1:size(Temp,2),1:size(Temp,3));X =X(:);Y =Y(:);
% % % % % % % % % % for i = 1:size(Temp,2)
% % % % % % % % % %     A = squeeze(Temp(i,:,:));
% % % % % % % % % %     T = A*0;
% % % % % % % % % %     if sum(A(:))
% % % % % % % % % %         [B,L] = bwboundaries(logical(A),'noholes');
% % % % % % % % % %         for j = 1:length(B)
% % % % % % % % % %             boundary = B{j};
% % % % % % % % % %             k=LineCurvature2D(boundary);
% % % % % % % % % %             ind = find(k >0);
% % % % % % % % % %             boundary(ind,:) = [];
% % % % % % % % % %             indold = 0;
% % % % % % % % % %             while (isequal(ind,indold) ==0)|~isempty(ind)
% % % % % % % % % %                 indold = ind;
% % % % % % % % % %                 k=LineCurvature2D(boundary);
% % % % % % % % % %                 indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
% % % % % % % % % %                 %ind = find(k >0);
% % % % % % % % % %                 indr = ismember(A(indid),[3000 3001 3014 3000 3006 3016 3021 7000 3032 3033 3034 3035]);
% % % % % % % % % %                 ind = find((k >0)&indr==0);
% % % % % % % % % %                 boundary(ind,:) = [];
% % % % % % % % % %             end
% % % % % % % % % %             boundary = [boundary;boundary(1,:)];
% % % % % % % % % %             indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
% % % % % % % % % %             indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
% % % % % % % % % %             T(indin) = j;
% % % % % % % % % %         end
% % % % % % % % % %     end
% % % % % % % % % %     Temp(i,:,:) = T;
% % % % % % % % % % end
% % % % % % % % % % Temp = imdilate(Temp,strel(ones(3,3,3)));
% % % % % % % % [fv] = Ver_Ext_vox('',0,logical(Temp));
% % % % % % % % Surf.SurfData.vertices = fv.vertices;
% % % % % % % % Surf.SurfData.faces = fv.faces;
varargout{1} = Surf;
return

function T = Hull_from_Mask_Slice(A,Y,X);
T = A*0;
Told = A;
while sum(Told(:) - T(:))~=0
    Told = T;
    [B,L] = bwboundaries(logical(A),'noholes');
    for j = 1:length(B)
        boundary = B{j};
        k=LineCurvature2D(boundary);
        ind = find(k >0);
        boundary(ind,:) = [];
        indold = 0;
        while (isequal(ind,indold) ==0)|~isempty(ind)
            indold = ind;
            k=LineCurvature2D(boundary);
            indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
            ind = find((k >0));
            boundary(ind,:) = [];
        end
        boundary = [boundary;boundary(1,:)];
        indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
        indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
        T(indin) = j;
    end
    A = T;
end


return;