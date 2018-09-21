function AnoMat = Generate_Deviation;
clc
hemi = 'lh';
InputDir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Qdec';
CFilest = sel_files(InputDir,['Y_' hemi '*']);
if hemi == 'lh'
    [cvar,ctab] = read_cfiles('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/label/lh.aparc.annot');
    [OutFiles, SurfF] = Exp_Surf('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/surf/lh.inflated', '0', '','', 'imp','n');Surf = SurfF{1};
elseif hemi == 'rh'
    [cvar,ctab] = read_cfiles('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/label/rh.aparc.annot');
    [OutFiles, SurfF] = Exp_Surf('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_avesubj_thmax05/surf/rh.inflated', '0', '','', 'imp','n');Surf = SurfF{1};
end

% Pruebas = strvcat('lh.thickness.COMP_avesubj_thmax05.mgh','lh.thickness.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'lh.thickness.fwhm0.COMP_avesubj_thmax05.mgh','lh.thickness.fwhm0.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'lh.thickness.fwhm5.COMP_avesubj_thmax05.mgh','lh.thickness.fwhm5.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'lh.thickness.fwhm10.COMP_avesubj_thmax05.mgh','lh.thickness.fwhm10.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'lh.thickness.fwhm15.COMP_avesubj_thmax05.mgh','lh.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'lh.thickness.fwhm20.COMP_avesubj_thmax05.mgh','lh.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'lh.thickness.fwhm25.COMP_avesubj_thmax05.mgh','lh.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh');
% Pruebas = strvcat(Pruebas,'rh.thickness.COMP_avesubj_thmax05.mgh','rh.thickness.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'rh.thickness.fwhm0.COMP_avesubj_thmax05.mgh','rh.thickness.fwhm0.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'rh.thickness.fwhm5.COMP_avesubj_thmax05.mgh','rh.thickness.fwhm5.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'rh.thickness.fwhm10.COMP_avesubj_thmax05.mgh','rh.thickness.fwhm10.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'rh.thickness.fwhm15.COMP_avesubj_thmax05.mgh','rh.thickness.fwhm15.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'rh.thickness.fwhm20.COMP_avesubj_thmax05.mgh','rh.thickness.fwhm20.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh',...
%                   'rh.thickness.fwhm25.COMP_avesubj_thmax05.mgh','rh.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh');  
%               Pruebas = 'rh.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh';
        Contrasts = strvcat('Contrast_FTest','MAD-vs-OVI','OVI-vs-BCN','BCN-vs-ZGZ','ZGZ-vs-VAL','VAL-vs-VIT');
          

for i  =1:size(CFilest,1)
        [txt, M, mr_parms, volsz] = load_mgh(deblank(CFilest(i,:)));
        new = reshape(squeeze(txt),[163842 6 6]);
        meanss = mean(new,2);
        stdsss = [std(squeeze(new(:,:,1))')' std(squeeze(new(:,:,2))')' std(squeeze(new(:,:,3))')' std(squeeze(new(:,:,4))')' std(squeeze(new(:,:,5))')' std(squeeze(new(:,:,6))')'];
        
end
    newname = [InputDir filesep 'Qdec' filesep 'Y_' deblank(Pruebas(j,:))];
    save_mgh(txt, newname, M, mr_parms);
%     tic;
%     for k=1:size(Contrasts,1)
%         cad=['mri_glmfit --glmdir Results --y ' newname ' --fsgd qdec.fsgd --C ' deblank(Contrasts(k,:)) '.mtx'];
%         %system(cad);
%         ind =strfind(Pruebas(j,:),'.');names=deblank(Pruebas(j,1:ind(end)-1));
%         %movefile([InputDir filesep 'Results' filesep deblank(Contrasts(k,:)) filesep 'sig.mgh'],[InputDir filesep 'Results' filesep 'Sig_' names '_' deblank(Contrasts(k,:)) '.mgh']);
%         %disp([InputDir filesep 'Results' filesep deblank(Contrasts(k,:)) filesep 'sig.mgh' '  =================>>> ' InputDir filesep 'Results' filesep 'Sig_' names '_' deblank(Contrasts(k,:)) '.mgh']);
%     end
%     toc;


%AnoMat = reshape(AnoMat,[Nsubjects Nsites]);
return