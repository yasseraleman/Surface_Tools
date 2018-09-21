function Surf_Plot_CaliStatResults(name);

InputDir = '/home/yaleman/Stat_Results';
Ydir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Qdec/Ys'
srootdir = '/home/yaleman/Stat_Results/Pics';
mkdir('/home/yaleman/Stat_Results','Pics')



% nameall = strvcat('*.thickness.fwhm0.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm0.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% 
%  nameall = strvcat(nameall,'*.thickness.fwhm5.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm5.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% 
%  nameall = strvcat(nameall,'*.thickness.fwhm10.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm10.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% 
%  nameall = strvcat(nameall,'*.thickness.fwhm15.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% %
%  nameall = strvcat(nameall,'*.thickness.fwhm20.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% 
%  nameall = strvcat(nameall,'*.thickness.fwhm25.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');

  nameall = strvcat('*.thickness.fwhm0.COMP_avesubj_thmax05_BCN-vs-ZGZ_RM.mgh','*.thickness.fwhm0.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-ZGZ_RM.mgh');
 % nameall = strvcat(nameall,'*.thickness.fwhm5.COMP_avesubj_thmax05_BCN-vs-ZGZ_RM.mgh','*.thickness.fwhm5.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-ZGZ_RM.mgh');
  nameall = strvcat(nameall,'*.thickness.fwhm10.COMP_avesubj_thmax05_BCN-vs-ZGZ_RM.mgh','*.thickness.fwhm10.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-ZGZ_RM.mgh');
 % nameall = strvcat(nameall,'*.thickness.fwhm15.COMP_avesubj_thmax05_BCN-vs-ZGZ_RM.mgh','*.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-ZGZ_RM.mgh');
   nameall = strvcat('*.thickness.fwhm20.COMP_avesubj_thmax05_BCN-vs-ZGZ_RM.mgh','*.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-ZGZ_RM.mgh');
%  % nameall = strvcat(nameall,'*.thickness.fwhm25.COMP_avesubj_thmax05_BCN-vs-ZGZ_RM.mgh','*.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-ZGZ_RM.mgh');
%   
% % 
   nameall = strvcat(nameall,'*.thickness.fwhm0.COMP_avesubj_thmax05_BCN-vs-VAL_RM.mgh','*.thickness.fwhm0.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-VAL_RM.mgh');
%   nameall = strvcat(nameall,'*.thickness.fwhm5.COMP_avesubj_thmax05_BCN-vs-VAL_RM.mgh','*.thickness.fwhm5.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-VAL_RM.mgh');
   nameall = strvcat(nameall,'*.thickness.fwhm10.COMP_avesubj_thmax05_BCN-vs-VAL_RM.mgh','*.thickness.fwhm10.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-VAL_RM.mgh');
%   nameall = strvcat(nameall,'*.thickness.fwhm15.COMP_avesubj_thmax05_BCN-vs-VAL_RM.mgh','*.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-VAL_RM.mgh');
   nameall = strvcat(nameall,'*.thickness.fwhm20.COMP_avesubj_thmax05_BCN-vs-VAL_RM.mgh','*.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-VAL_RM.mgh');
%   nameall = strvcat(nameall,'*.thickness.fwhm25.COMP_avesubj_thmax05_BCN-vs-VAL_RM.mgh','*.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_BCN-vs-VAL_RM.mgh');
  
%   nameall = strvcat(nameall,'*.thickness.fwhm15.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% % %
%   nameall = strvcat(nameall,'*.thickness.fwhm20.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');
% % 
%   nameall =  strvcat(nameall,'*.thickness.fwhm25.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh','*.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh');



% name = '*.thickness.fwhm5.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh'
% name = '*.thickness.fwhm5.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh'
%
% name = '*.thickness.fwhm10.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh'
% name = '*.thickness.fwhm10.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh'
%
% name = '*.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh'
% name = '*.thickness.fwhm15.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh'
%
% name = '*.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4_Contrast_FTest_RM.mgh'
% name = '*.thickness.fwhm20.COMP_avesubj_thmax05_Contrast_FTest_RM.mgh'







%if strfind(name,'lh')
[labelsl,ctab] = read_cfiles('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/label/lh.aparc.annot');
[OutFiles, SurfF] = Exp_Surf('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/surf/lh.inflated', '0', '','', 'imp','n');Surfl = SurfF{1};
%elseif strfind(name,'rh')
[labelsr,ctab] = read_cfiles('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/label/rh.aparc.annot');
[OutFiles, SurfF] = Exp_Surf('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/surf/rh.inflated', '0', '','', 'imp','n');Surfr = SurfF{1};
for i =1:size(nameall,1)
    name = deblank(nameall(i,:));
    CFilest = sel_files(InputDir, name);


    [cvarl,var] = read_cfiles(deblank(CFilest(1,:)));
    p = power(repmat(10,[size(cvarl,1) 1]),-1*cvarl);p(p>1)=1;[pID,pN] = FDR(p,0.05);if ~isempty(pID);ind = find(p>pID);ind = find(p>pID);else pID = 0.05;end
    %end
    p(ind) = 0;
    Surfl.Is = -1*p;
    [Colors] = Surf_Color(Surfl,'hot'); Colors(ind,:) = repmat([0.37 0.37 0.37],[length(ind) 1]);Surfl.SurfData.FaceVertexCData = Colors;

%     if strfind(deblank(CFilest(1,:)),'-vs-')
%         [cols] = detectcols(deblank(CFilest(1,:)));
%         ind = strfind(name,'thmax05');nname = name(3:ind+6);
%         yname = [Ydir filesep 'Y_lh.' nname '.mgh'];
% 
%         [txt, M, mr_parms, volsz] = load_mgh(yname);
%         new = reshape(squeeze(txt),[163842 6 6]);
%         meanss = squeeze(mean(new,3));
%         intermeans = meanss(:,cols);meanmaq = intermeans(:,1)-intermeans(:,2); Surfl.Is = abs(meanmaq);
%         stdsss = [std(squeeze(new(:,:,1))')' std(squeeze(new(:,:,2))')' std(squeeze(new(:,:,3))')' std(squeeze(new(:,:,4))')' std(squeeze(new(:,:,5))')' std(squeeze(new(:,:,6))')'];
%         interstdss = stdsss(:,cols);
%         %factors = abs(2*(intermeans(:,1)-intermeans(:,2))./(interstdss(:,1)+interstdss(:,2)+eps));factors(factors>10) = 0;indnan = isnan(factors);factors(indnan)=0;
%         [Colors] = Surf_Color(Surfl,'jet');% Colors(ind,:) = repmat([0.37 0.37 0.37],[length(ind) 1]);
%         Surfl.SurfData.FaceVertexCData = Colors;
%     end

    colordef black;hf = figure('numbertitle','off','Color','black','name',['IBASPM Surface Ploting...:  Contrast  T' name(4:end-5)]);
    subplot(2,2,1)
    %strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','gouraud');
    Plot_Surf(Surfl,'hot');
    axis off;axis tight;axis equal;
    h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
    % if strfind(name,'lh')
    view([270 0]);
    % elseif strfind(name,'rh')
    %     view([90 0]);
    %end
    camlight;


    subplot(2,2,2)
    Plot_Surf(Surfl,'hot');
    axis off;axis tight;axis equal;
    h=title(['Medial View']);set(h,'FontSize',15,'FontName','Arial');
    % if strfind(name,'lh')
    view([90 0]);
    % elseif strfind(name,'rh')
    %     view([270 0]);
    % end
    camlight;




    %%%%%%%%%%%%%%%%%%%%Right Hemisphere
    [cvarr,var] = read_cfiles(deblank(CFilest(2,:)));
    p = power(repmat(10,[size(cvarr,1) 1]),-1*cvarr);p(p>1)=1;[pID,pN] = FDR(p,0.05);if ~isempty(pID);ind = find(p>pID);else pID = 0.05;end
    p(ind) = 0;
    Surfr.Is = -1*p;
    [Colors] = Surf_Color(Surfr,'hot'); Colors(ind,:) = repmat([0.37 0.37 0.37],[length(ind) 1]);Surfr.SurfData.FaceVertexCData = Colors;
%     if strfind(deblank(CFilest(2,:)),'-vs-')
%         [cols] = detectcols(deblank(CFilest(1,:)));
%         ind = strfind(name,'thmax05');nname = name(3:ind+6);
%         yname = [Ydir filesep 'Y_rh.' nname '.mgh'];
% 
%         [txt, M, mr_parms, volsz] = load_mgh(yname);
%         new = reshape(squeeze(txt),[163842 6 6]);
%         meanss = squeeze(mean(new,3));
%         intermeans = meanss(:,cols);meanmaq = intermeans(:,1)-intermeans(:,2);Surfr.Is = abs(meanmaq);
%         stdsss = [std(squeeze(new(:,:,1))')' std(squeeze(new(:,:,2))')' std(squeeze(new(:,:,3))')' std(squeeze(new(:,:,4))')' std(squeeze(new(:,:,5))')' std(squeeze(new(:,:,6))')'];
%         interstdss = stdsss(:,cols);
%         %factors = abs(2*(intermeans(:,1)-intermeans(:,2))./(interstdss(:,1)+interstdss(:,2)+eps));factors(factors>10) = 0;indnan = isnan(factors);factors(indnan)=0;
%         [Colors] = Surf_Color(Surfr,'jet');% Colors(ind,:) = repmat([0.37 0.37 0.37],[length(ind) 1]);
%         Surfr.SurfData.FaceVertexCData = Colors;
%     end
    
    
    
    subplot(2,2,3)
    %strsurf=patch(Surf.SurfData,'edgecolor','none','tag', 'patch','facelighting','gouraud');
    Plot_Surf(Surfr,'hot');
    axis off;axis tight;axis equal;
    %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
    % if strfind(name,'lh')
    %view([270 0]);
    % elseif strfind(name,'rh')
    view([90 0]);
    %end
    camlight;


    subplot(2,2,4)
    Plot_Surf(Surfr,'hot');
    axis off;axis tight;axis equal;
    %h=title(['Medial View']);set(h,'FontSize',15,'FontName','Arial');
    % if strfind(name,'lh')
    %    view([90 0]);
    % elseif strfind(name,'rh')
    view([270 0]);
    % end
    camlight;
    %if sum(Surf(1).Is-floor(Surfr(1).Is)) ~=0
%         h = colorbar;
         mat = max([Surfr(1).Is;Surfl.Is]);%mat=mat-mat/2;
%         set(h,'YTickLabel',[0 mat/2 mat]);
    %end
    disp([name  '               '    num2str(mat)])    
    
    
    
    
    imname = [srootdir filesep name(3:end-4)  '.jpg'];
    print( hf, '-djpeg','-r200', imname);
    close(hf);
    
    
    
end









function [cols] = detectcols(contrast);

if strfind(contrast,'MAD-vs-OVI')
    cols= [1 2];
elseif strfind(contrast,'MAD-vs-BCN')
    cols= [1 3];
elseif strfind(contrast,'MAD-vs-ZGZ')
    cols= [1 4];
elseif strfind(contrast,'MAD-vs-VAL')
    cols= [1 5];
elseif strfind(contrast,'MAD-vs-VIT')
    cols= [1 6];
elseif strfind(contrast,'OVI-vs-BCN')
    cols= [2 3];
elseif strfind(contrast,'OVI-vs-ZGZ')
    cols= [2 4];
elseif strfind(contrast,'OVI-vs-VAL')
    cols= [2 5];
elseif strfind(contrast,'OVI-vs-VIT')
    cols= [2 6];
elseif strfind(contrast,'BCN-vs-ZGZ')
    cols= [3 4];
elseif strfind(contrast,'BCN-vs-VAL')
    cols= [3 5];
elseif strfind(contrast,'BCN-vs-VIT')
    cols= [3 6];
elseif strfind(contrast,'ZGZ-vs-VAL')
    cols= [4 5];
elseif strfind(contrast,'ZGZ-vs-VIT')
    cols= [4 6];
elseif strfind(contrast,'VAL-vs-VIT')
    cols= [5 6];
end
return;



















function custom_plotsurf(Surf,cl);
if isfield(Surf.SurfData,'VertexNormals');
    Surf.SurfData = rmfield(Surf.SurfData,'VertexNormals');
end
if nargin ==1
    cl = 'jet';
end
if ~isfield(Surf.SurfData,'FaceVertexCData')
    [Colors] = Surf_Color(Surf,cl);
    Surf.SurfData.FaceVertexCData = Colors;
    Surf.SurfData.FaceColor = 'interp';
end
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
axis off;axis tight;axis equal;
return