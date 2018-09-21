Outdir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Calibration_Results/Calibration_Results_6sujetos';
hemis = strvcat('lh','rh');
filters = strvcat('COMP_avesubj_thmax05','fwhm0.COMP_avesubj_thmax05','fwhm5.COMP_avesubj_thmax05','fwhm10.COMP_avesubj_thmax05','fwhm15.COMP_avesubj_thmax05',...
    'fwhm20.COMP_avesubj_thmax05','fwhm25.COMP_avesubj_thmax05','COMP_avesubj_thmax10',...
'fwhm0.COMP_avesubj_thmax10','fwhm5.COMP_avesubj_thmax10','fwhm10.COMP_avesubj_thmax10','fwhm15.COMP_avesubj_thmax10','fwhm20.COMP_avesubj_thmax10',...
'fwhm25.COMP_avesubj_thmax10');
for i = 1:size(hemis,1)
    for j = 1:size(filters,1)
        disp(['Processing  Hemi: ' upper(hemis(i,:)) '   . Filter:  ' filters(j,:) ])
        Mean_cort_results([deblank(filters(j,:)) '.'],hemis(i,:));
    end
end
