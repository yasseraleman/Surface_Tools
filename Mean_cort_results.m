function Mean_cort_results(tempfile,hemi);
% tempfile = 'fwhm0.COMP_avesubj_thmax05.';
% hemi = 'lh';
%ext = ['.' tempfile];
method =4;
InputDir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer';
ind = strfind(tempfile,'.'); 
if length(ind)>1
    ntempfile = tempfile(ind(1)+1:ind(end)-1);
else
    ntempfile = tempfile(1:ind-1);
end
[cvar,ctab] = read_cfiles(['/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/' ntempfile '/label/lh.aparc.annot']);
[Sid, Sname, colx, coly, colz, offset] = textread('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/COMP_00002__MAD-20090606-T1-MAD-1/label/aparc.annot.ctab','%f%s%u%u%u%u',36);
    labels = colx + coly.*2^8 + colz*2^16 + offset*2^24;labels(5)=[];
CFilest = sel_files(InputDir,[hemi '.thickness.' tempfile 'mgh']);
if (method ==1)|(method ==2)
    load([InputDir filesep 'Calibration_Results' filesep hemi '.thickness.' tempfile 'mat']);
    Num = size(ResultsTV.NewNcth,2);
    
elseif (method ==3)|(method ==4)
    CFilest = Sort_Data(CFilest,'-1');
    load([InputDir filesep 'Calibration_Results' filesep 'Calibration_Results_6sujetos' filesep hemi '.thickness.' tempfile 'mat']);
    Num = size(ResultsTV.NewNcth,2);
end

Ns = size(labels,1);
for i = 1:size(CFilest,1)
    [txt, M, mr_parms, volsz] = load_mgh(deblank(CFilest(i,:)));
    ind = strfind(deblank(CFilest(i,:)),filesep);
    Statfile =[deblank(CFilest(i,1:ind(end-1))) 'stats' filesep hemi '.aparc.stats'];ind = strfind(deblank(CFilest(i,:)),'.');
    if strfind(deblank(CFilest(i,:)),'COMP_00002__MAD');
        col =1;
        if ~isempty(strfind(deblank(CFilest(i,:)),'-MAD-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-OVI-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-BCN-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-ZGZ-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VAL-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VIT-2'));
            col =2*col;
        end
    elseif strfind(deblank(CFilest(i,:)),'COMP_00006__OVI');
        col = 2;
        if ~isempty(strfind(deblank(CFilest(i,:)),'-MAD-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-OVI-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-BCN-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-ZGZ-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VAL-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VIT-2'));
            col =2*col;
        end
    elseif strfind(deblank(CFilest(i,:)),'COMP_00005__BCN');
        col = 3;
        if ~isempty(strfind(deblank(CFilest(i,:)),'-MAD-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-OVI-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-BCN-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-ZGZ-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VAL-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VIT-2'));
            col =2*col;
        end
    elseif strfind(deblank(CFilest(i,:)),'COMP_00003__ZGZ');
        col = 4;
        if ~isempty(strfind(deblank(CFilest(i,:)),'-MAD-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-OVI-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-BCN-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-ZGZ-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VAL-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VIT-2'));
            col =2*col;
        end
    elseif strfind(deblank(CFilest(i,:)),'COMP_00007__VAL');
        col = 5;
        if ~isempty(strfind(deblank(CFilest(i,:)),'-MAD-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-OVI-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-BCN-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-ZGZ-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VAL-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VIT-2'));
            col =2*col;
        end
    elseif strfind(deblank(CFilest(i,:)),'COMP_00004__VIT');
        col = 6;
        if  ~isempty(strfind(deblank(CFilest(i,:)),'-MAD-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-OVI-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-BCN-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-ZGZ-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VAL-2'))|~isempty(strfind(deblank(CFilest(i,:)),'-VIT-2'));
            col =2*col;
        end
    end
    if method ==1
         newname = [deblank(CFilest(i,1:ind(end)-1)) '_Correct_' ntempfile '_' num2str(Num) '_M1.mgh'];
         save_mgh(ResultsTV.NewNcth(:,col), newname, M, mr_parms);
         
         for i = 1:Ns
             ind = find(cvar ==labels(i));
             meancth(i,1)= mean(ResultsTV.NewNcth(ind,col));
             meancth(i,2)= std(ResultsTV.NewNcth(ind,col));
         end
    read_write_statfile(Statfile,Num,tempfile,meancth,method);
    elseif method==2
        newname = [deblank(CFilest(i,1:ind(end)-1)) '_Correct_' ntempfile '_' num2str(Num) '_M2.mgh'];
        if strfind(deblank(CFilest(i,:)),'-MAD-')|strfind(deblank(CFilest(i,:)),'-MAD-2');
            colt =1;
        elseif strfind(deblank(CFilest(i,:)),'-OVI-')|strfind(deblank(CFilest(i,:)),'-OVI-2');
            colt = 2;
        elseif strfind(deblank(CFilest(i,:)),'-BCN-')|strfind(deblank(CFilest(i,:)),'-BCN-2');
            colt = 3;
        elseif strfind(deblank(CFilest(i,:)),'-ZGZ-')|strfind(deblank(CFilest(i,:)),'-ZGZ-2');
            colt = 4;
        elseif strfind(deblank(CFilest(i,:)),'-VAL-')|strfind(deblank(CFilest(i,:)),'-VAL-2');
            colt = 5;
        elseif strfind(deblank(CFilest(i,:)),'-VIT-')|strfind(deblank(CFilest(i,:)),'-VIT-2');
            colt = 6;
        end
        cthcorr = (txt-ResultsTV.Ct(:,colt))./(ResultsTV.Bt(:,colt)+eps);
        save_mgh(cthcorr, newname, M, mr_parms);
        meancth = [0 0];
        for i = 2:Ns
             ind = find(cvar ==labels(i));
              ind2 = find(abs(cthcorr(ind))>10);ind(ind2)=[];
             meancth(i,1)= mean(cthcorr(ind));
             ind2 = find(abs(cthcorr(ind))>3*meancth(i,1));ind(ind2)=[];
             meancth(i,2)= std(cthcorr(ind));
         end
        read_write_statfile(Statfile,Num,tempfile,abs(meancth),method);
    elseif method==3
        newname = [deblank(CFilest(i,1:ind(end)-1)) '_Correct_' ntempfile '_' num2str(Num) '_M3.mgh'];
        save_mgh(ResultsTV.NewNcth(:,col), newname, M, mr_parms);
        for i = 1:Ns
             ind = find(cvar ==labels(i));
             meancth(i,1)= mean(ResultsTV.NewNcth(ind,col));
             meancth(i,2)= std(ResultsTV.NewNcth(ind,col));
         end
    read_write_statfile(Statfile,Num,tempfile,meancth,method);
    elseif method==4
        newname = [deblank(CFilest(i,1:ind(end)-1)) '_Correct_' ntempfile '_' num2str(Num) '_M4.mgh'];
        cthcorr = (txt-ResultsTV.Ct(:,col))./(ResultsTV.Bt(:,col)+eps);
        save_mgh(cthcorr, newname, M, mr_parms);
        meancth = [0 0];
        for i = 2:Ns
             ind = find(cvar ==labels(i));
              ind2 = find(abs(cthcorr(ind))>10);ind(ind2)=[];
             meancth(i,1)= mean(cthcorr(ind));
             ind2 = find(abs(cthcorr(ind))>3*meancth(i,1));ind(ind2)=[];
             meancth(i,2)= std(cthcorr(ind));
         end
        read_write_statfile(Statfile,Num,tempfile,abs(meancth),method);
    end
end
return;