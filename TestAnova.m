function AnoMat = GenerateYs
clc


 InputDir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Qdec';
 Ydir = [InputDir filesep 'Ys'];
 Cdir = [InputDir filesep 'Contrasts'];
 Home = '/home/yaleman/Stat_Results';
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
              %Pruebas = 'rh.thickness.fwhm25.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_6_M4.mgh';
              
               Contrasts= strvcat('Contrast_FTest_RM','MAD-vs-OVI_RM','MAD-vs-BCN_RM','MAD-vs-ZGZ_RM','MAD-vs-VAL_RM','MAD-vs-VIT_RM','OVI-vs-BCN_RM','OVI-vs-ZGZ_RM','OVI-vs-VAL_RM','OVI-vs-VIT_RM',...
               'BCN-vs-ZGZ_RM','BCN-vs-VAL_RM','BCN-vs-VIT_RM','ZGZ-vs-VAL_RM','ZGZ-vs-VIT_RM','VAL-vs-VIT_RM');
                %Contrasts ='VAL-vs-VIT_RM';
              
         
%[txt, M, mr_parms, volsz] = load_mgh('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Qdec/y.mgh');
Filest = sel_files(Ydir,'Y_rh*');
cont = 0;Ns= size(Filest,1)*size(Contrasts,1);
H = waitbar(0,['Number of points in the surface  ' num2str(Ns)],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');

for j = 1:size(Filest,1)
%     CFilest = sel_files(InputDir,deblank(Pruebas(j,:)));
%     CFilest = Sort_Data(CFilest,'-1');
%     for i  =1:size(CFilest,1)
%         [cvar,ctab] = read_cfiles(deblank(CFilest(i,:)));
%         txt(:,1,1,i) = cvar;
%     end
%     newname = [InputDir filesep 'Qdec' filesep 'Y_' deblank(Pruebas(j,:))];
%     save_mgh(txt, newname, M, mr_parms);
    tic;
    for k=1:size(Contrasts,1)
        waitbar(i/Ns,H,['Correcting Point  ' num2str(cont) ' of ' num2str(Ns)]);
        cont = cont+1
        cad=['mri_glmfit --glmdir /media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Qdec/Results --y ' deblank(Filest(j,:)) ' --fsgd ' InputDir filesep 'qdec2.fsgd doss --C ' Cdir filesep deblank(Contrasts(k,:)) '.mtx'];
        system(cad);
         ind =strfind(Filest(j,:),filesep);names=deblank(Filest(j,ind(end)+1:end)); ind = strfind(names,'.');names = names(1:ind(end)-1);
        copyfile([InputDir filesep 'Results' filesep deblank(Contrasts(k,:)) filesep 'sig.mgh'],[Home filesep 'Sig_' names '_' deblank(Contrasts(k,:)) '.mgh']);
        disp([InputDir filesep 'Results' filesep deblank(Contrasts(k,:)) filesep 'sig.mgh' '  =================>>> ' InputDir filesep 'Results' filesep 'Sig_' names '_' deblank(Contrasts(k,:)) '.mgh']);
        %p = power(repmat(10,[size(cvar1,1) 1]),-1*cvar1);[pID,pN] = FDR(p,0.05);ind = find(p>pID);pt = p;pt(ind) = 0;Surf.Is = pt;Plot_Surf(Surf)
    end
    toc;
    
end
close(H)
%AnoMat = reshape(AnoMat,[Nsubjects Nsites]);
return