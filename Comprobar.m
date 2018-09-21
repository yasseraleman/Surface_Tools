function Comprobar(hemi);
%% Comprobar si los datos estan corregidos correctamente.
%%
hemi = 'lh';
% %filter = '.pial';
%afilter = '.aparc.annot';
estnumber = ['-1'];
InputDir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer';
mkdir('/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Calibration_Results/Calibration_Results_6sujetos','Pics')
srootdir = '/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Calibration_Results/Calibration_Results_6sujetos/Pics'
%SurfFiles = sel_files(InputDir,[hemi filter]);

%%  Sujeto 1
cfilter = '*thickness.COMP_avesubj_thmax05.mgh';
Califile = ['/media/Data/ibas/Compatibilidad_Proyecto_PEPs/Procesamiento_FreeSurfer/Calibration_Results/Calibration_Results_6sujetos' filesep hemi '.' cfilter(2:end-4) '.mat'];
load(Califile);
CFilest = sel_files(InputDir,[hemi cfilter]);
CFiles = Sort_Data(CFilest,estnumber);%CFiles2 = Sort_Data(CFilest,['-2' filesep]);
Ns = size(CFiles,1);
for i = 1:Ns
    [cvar,ctab] = read_cfiles(deblank(CFiles(i,:)));
    mcth(i,1:length(cvar)) = cvar;
end
%%
%%  Sujeto 2
cfilter1 = '*.thickness.COMP_avesubj_thmax05_Correct_COMP_avesubj_thmax05_12_M2.mgh';
CFilest = sel_files(InputDir,[hemi cfilter1]);
CFiles = Sort_Data(CFilest,estnumber);%CFiles2 = Sort_Data(CFilest,['-2' filesep]);
Ns = size(CFiles,1);
for i = 1:Ns
    [cvar,ctab] = read_cfiles(deblank(CFiles(i,:)));
    mcth1(i,1:length(cvar)) = cvar;
end
%%

ind = strfind(cfilter,'.');cads = cfilter(ind(1)+1:ind(2)-1);
tempfile = sel_files([InputDir filesep cads],[hemi '.aparc.annot']);
[vertices, cvar1, colortable] = read_annotation(deblank(tempfile));
strs = unique(cvar1);strs(strs==0) =[];
vals = zeros(36,size(strs,1)); names = '';
for i = 1:size(strs,1)
    ind = find(cvar1 == strs(i));
    indn= find(colortable.table(:,5) == strs(i));
    names =strvcat(names,colortable.struct_names{indn});
    Temp = mcth(:,ind); valso(:,i) = mean(Temp,2);
    Temp = mcth1(:,ind); valsco(:,i) = mean(Temp,2);
    Temp = ResultsTV.Bt(ind,:); B(i,:) = mean(Temp,1);
    Temp = ResultsTV.Ct(ind,:); C(i,:) = mean(Temp,1);
end

%% Removing unknown and Corpus callosum
ind = ismember(names,strvcat('corpuscallosum','unknown'),'rows');
names(ind,:) = [];
valso(:,ind) = [];
valsco(:,ind) = [];
B(ind,:) = [];
C(ind,:) = [];


% Lobulo Frontal
FL = strvcat('superiorfrontal','rostralmiddlefrontal','caudalmiddlefrontal','parsopercularis','parstriangularis','parsorbitalis','lateralorbitofrontal','medialorbitofrontal','frontalpole','precentral','paracentral');
ind = ismember(names,FL,'rows');
names = strvcat(names,'frontal-lobe');
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

% Lobulo Parietal
PL = strvcat('postcentral','supramarginal','superiorparietal','inferiorparietal','precuneus');
ind = ismember(names,PL,'rows');
names = strvcat(names,'parietal-lobe');
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

% Lobulo Temporal Medial
TML = strvcat('entorhinal','parahippocampal','temporalpole','fusiform');
ind = ismember(names,TML,'rows');
names = strvcat(names,'medial-temporal-lobe');
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

% Lobulo Temporal Lateral
TLL = strvcat('superiortemporal','bankssts','middletemporal','inferiortemporal','transversetemporal');
ind = ismember(names,TLL,'rows');
names = strvcat(names,'lateral-temporal-lobe');
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

% Lobulo Occipital
OL = strvcat('lingual','pericalcarine','cuneus','lateraloccipital');
ind = ismember(names,OL,'rows');
names = strvcat(names,'occipital-lobe');
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

% Corteza Cingulada
CC = strvcat('rostralanteriorcingulate','caudalanteriorcingulate','posteriorcingulate','isthmuscingulate');
ind = ismember(names,CC,'rows');
names = strvcat(names,'cingulate-cortex');
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

% Hemisphere
HeM = strvcat(FL,PL,TML,TLL,OL,CC);
ind = ismember(names,HeM,'rows');
names = strvcat(names,[hemi(1) '-hemisphere']);
valso(:,end+1) = mean(valso(:,ind),2);
valsco(:,end+1) = mean(valsco(:,ind),2);
B(end+1,:) = mean(B(ind,:));
C(end+1,:) = mean(C(ind,:));

%%
m = strvcat('x','*','+','o','s','d');
for i = 1:size(valsco,2)
    hf = figure('Color',[0 0 0]);
    colordef black;hold on
    for j = 0:5
        plot(1,valso(j*6+1,i),'Color',[1 0 0],'Marker',m(j+1),'Markersize',10);
        plot(1.5,valso(j*6+2,i),'Color',[0 1 0],'Marker',m(j+1),'Markersize',10);
        plot(2,valso(j*6+3,i),'Color',[0 0 1],'Marker',m(j+1),'Markersize',10);
        plot(2.5,valso(j*6+4,i),'Color',[1 1 0],'Marker',m(j+1),'Markersize',10);
        plot(3,valso(j*6+5,i),'Color',[0 1 1],'Marker',m(j+1),'Markersize',10);
        plot(3.5,valso(j*6+6,i),'Color',[1 1 1],'Marker',m(j+1),'Markersize',10);
        
        plot(7,valsco(j*6+1,i),'Color',[1 0 0],'Marker',m(j+1),'Markersize',10);
        plot(7.5,valsco(j*6+2,i),'Color',[0 1 0],'Marker',m(j+1),'Markersize',10);
        plot(8,valsco(j*6+3,i),'Color',[0 0 1],'Marker',m(j+1),'Markersize',10);
        plot(8.5,valsco(j*6+4,i),'Color',[1 1 0],'Marker',m(j+1),'Markersize',10);
        plot(9,valsco(j*6+5,i),'Color',[0 1 1],'Marker',m(j+1),'Markersize',10);
        plot(9.5,valsco(j*6+6,i),'Color',[1 1 1],'Marker',m(j+1),'Markersize',10);
    end
    
    set(gca,'XLim',[0.5 10]);
    set(gca,'YLim',[0 5]);
    set(gca,'XTick',[0:0.5:10]);
    set(gca,'XTickLabel',['   ';'   ';'MAD';'OVI';'BCN';'ZGZ';'VAL';'VIT';'   ';'   ';'   ';'   ';'   ';'   ';'MAD';'OVI';'BCN';'ZGZ';'VAL';'VIT';'   ']);
    
    
    h=Title(['Uncorrected  vs Corrected. Structure name:  ' upper(names(i,:))]);
    YLabel('Cortical Thickness (mm)');
    XLabel('Scanner');
    %%============== Correction Factors ======================================%
    %% Madrid
    text(7.5,4.5,'Corrected Data FWHM = 0');
    text(7,0.9,'Corr-Fac');
    text(7,0.8,['B = ' num2str(round(B(i,1)*10)/10)]);
    text(7,0.7,['C = ' num2str(round(C(i,1)*10)/10)]);
    %% Oviedo
    text(7.5,0.9,'Corr-Fac');
    text(7.5,0.8,['B = ' num2str(round(B(i,2)*10)/10)]);
    text(7.5,0.7,['C = ' num2str(round(C(i,2)*10)/10)]);
    %% Barcelona
    text(8,0.9,'Corr-Fac');
    text(8,0.8,['B = ' num2str(round(B(i,3)*10)/10)]);
    text(8,0.7,['C = ' num2str(round(C(i,3)*10)/10)]);
    %% Zaragoza
    text(8.5,0.9,'Corr-Fac');
    text(8.5,0.8,['B = ' num2str(round(B(i,4)*10)/10)]);
    text(8.5,0.7,['C = ' num2str(round(C(i,4)*10)/10)]);
    %% Valencia
    text(9,0.9,'Corr-Fac');
    text(9,0.8,['B = ' num2str(round(B(i,5)*10)/10)]);
    text(9,0.7,['C = ' num2str(round(C(i,5)*10)/10)]);
    %% Vitoria
    text(9.5,0.9,'Corr-Fac');
    text(9.5,0.8,['B = ' num2str(round(B(i,6)*10)/10)]);
    text(9.5,0.7,['C = ' num2str(round(C(i,6)*10)/10)]);
    %%============== Correction Factors ======================================%
    
    %%============== Factores Descriptivos======================================%
    %% Madrid
    text(6.8,4.3,'Desc-Var');
    text(6.8,4.2,['Mean = ' num2str(round(mean(valsco(1:6:end,i))*10)/10)]);
    text(6.8,4.1,['Std = ' num2str(round(std(valsco(1:6:end,i))*10)/10)]);
    %% Oviedo
    text(7.3,4.3,'Desc-Var');
    text(7.3,4.2,['Mean = ' num2str(round(mean(valsco(2:6:end,i))*10)/10)]);
    text(7.3,4.1,['Std = ' num2str(round(std(valsco(2:6:end,i))*10)/10)]);
    %% Barcelona
    text(7.8,4.3,'Desc-Var');
    text(7.8,4.2,['Mean = ' num2str(round(mean(valsco(3:6:end,i))*10)/10)]);
    text(7.8,4.1,['Std = ' num2str(round(std(valsco(3:6:end,i))*10)/10)]);
    %% Zaragoza
    text(8.3,4.3,'Desc-Var');
    text(8.3,4.2,['Mean = ' num2str(round(mean(valsco(4:6:end,i))*10)/10)]);
    text(8.3,4.1,['Std = ' num2str(round(mean(valsco(4:6:end,i))*10)/10)]);
    %% Valencia
    text(9,4.3,'Desc-Var');
    text(9,4.2,['Mean = ' num2str(round(mean(valsco(5:6:end,i))*10)/10)]);
    text(9,4.1,['Std = ' num2str(round(std(valsco(5:6:end,i))*10)/10)]);
    %% Vitoria
    text(9.5,4.3,'Desc-Var');
    text(9.5,4.2,['Mean = ' num2str(round(mean(valsco(6:6:end,i))*10)/10)]);
    text(9.5,4.1,['Std = ' num2str(round(std(valsco(6:6:end,i))*10)/10)]);
    %%============== Factores Descriptivos ===============================%
    
     %%============== Factores Descriptivos======================================%
     text(1.5,4.5,'Uncorrected Data FWHM = 0');
    %% Madrid
    text(0.8,4.3,'Desc-Var');
    text(0.8,4.2,['Mean = ' num2str(round(mean(valso(1:6:end,i))*10)/10)]);
    text(0.8,4.1,['Std = ' num2str(round(std(valso(1:6:end,i))*10)/10)]);
    %% Oviedo
    text(1.3,4.3,'Desc-Var');
    text(1.3,4.2,['Mean = ' num2str(round(mean(valso(2:6:end,i))*10)/10)]);
    text(1.3,4.1,['Std = ' num2str(round(std(valso(2:6:end,i))*10)/10)]);
    %% Barcelona
    text(1.8,4.3,'Desc-Var');
    text(1.8,4.2,['Mean = ' num2str(round(mean(valso(3:6:end,i))*10)/10)]);
    text(1.8,4.1,['Std = ' num2str(round(std(valso(3:6:end,i))*10)/10)]);
    %% Zaragoza
    text(2.3,4.3,'Desc-Var');
    text(2.3,4.2,['Mean = ' num2str(round(mean(valso(4:6:end,i))*10)/10)]);
    text(2.3,4.1,['Std = ' num2str(round(mean(valso(4:6:end,i))*10)/10)]);
    %% Valencia
    text(2.8,4.3,'Desc-Var');
    text(2.8,4.2,['Mean = ' num2str(round(mean(valso(5:6:end,i))*10)/10)]);
    text(2.8,4.1,['Std = ' num2str(round(std(valso(5:6:end,i))*10)/10)]);
    %% Vitoria
    text(3.3,4.3,'Desc-Var');
    text(3.3,4.2,['Mean = ' num2str(round(mean(valso(6:6:end,i))*10)/10)]);
    text(3.3,4.1,['Std = ' num2str(round(std(valso(6:6:end,i))*10)/10)]);
    %%============== Factores Descriptivos ===============================%
    
    
    imname = [srootdir filesep upper(hemi) '-' deblank(names(i,:))  '.jpg'];
    print( hf, '-djpeg','-r200', imname);
    close(hf);
    % Subjects = legend(' Subject 1',' Subject 2',' Subject 3',' Subject 4',' Subject 5',' Subject 6',' Subject 7',' Subject 8',' Subject 9',' Subject 10',' Subject 11',' Subject 12',' Subject 13');
    
    
end
return;










a = 1;

function SurfFiles = Sort_Data(SurfFiles,estnumber,both);
contmad = 0;contvit = 0;contbcn = 0;contzgz = 0;contovi = 0;contval = 0;
contmad2 = 0;contvit2 = 0;contbcn2 = 0;contzgz2 = 0;contovi2 = 0;contval2 = 0;
for i =1:size(SurfFiles,1)
    indmad = strfind(deblank(SurfFiles(i,:)),['-MAD' estnumber]);
    if ~isempty(indmad)
        contmad = contmad+1;
        indxmad(contmad) =i;
    end
    
    indmad2 = strfind(deblank(SurfFiles(i,:)),['-MAD-2']);
    if ~isempty(indmad2)
        contmad2 = contmad2+1;
        indxmad2(contmad2) =i;
    end
    
    indvit = strfind(deblank(SurfFiles(i,:)),['-VIT' estnumber]);
    if ~isempty(indvit)
        contvit = contvit+1;
        indxvit(contvit) =i;
    end
    
    indvit2 = strfind(deblank(SurfFiles(i,:)),['-VIT-2']);
    if ~isempty(indvit2)
        contvit2 = contvit2+1;
        indxvit2(contvit2) =i;
    end
    
    indbcn = strfind(deblank(SurfFiles(i,:)),['-BCN' estnumber]);
    if ~isempty(indbcn)
        contbcn = contbcn+1;
        indxbcn(contbcn) =i;
    end
    
    indbcn2 = strfind(deblank(SurfFiles(i,:)),['-BCN-2']);
    if ~isempty(indbcn2)
        contbcn2 = contbcn2+1;
        indxbcn2(contbcn2) =i;
    end
    
    
    indzgz = strfind(deblank(SurfFiles(i,:)),['-ZGZ' estnumber]);
    if ~isempty(indzgz)
        contzgz = contzgz+1;
        indxzgz(contzgz) =i;
    end
    
    indzgz2 = strfind(deblank(SurfFiles(i,:)),['-ZGZ-2']);
    if ~isempty(indzgz2)
        contzgz2 = contzgz2+1;
        indxzgz2(contzgz2) =i;
    end
    
    indval = strfind(deblank(SurfFiles(i,:)),['-VAL' estnumber]);
    if ~isempty(indval)
        contval = contval+1;
        indxval(contval) =i;
    end
    
    indval2 = strfind(deblank(SurfFiles(i,:)),['-VAL-2']);
    if ~isempty(indval2)
        contval2 = contval2+1;
        indxval2(contval2) =i;
    end
    
    indovi = strfind(deblank(SurfFiles(i,:)),['-OVI' estnumber]);
    if ~isempty(indovi)
        contovi = contovi+1;
        indxovi(contovi) =i;
    end
    
    indovi2 = strfind(deblank(SurfFiles(i,:)),['-OVI-2']);
    if ~isempty(indovi2)
        contovi2 = contovi2+1;
        indxovi2(contovi2) =i;
    end
    
end

if exist('both')
    if strcmpi(both,'b')
        SurfFilest = SurfFiles;
        SurfFilest(1:12:(length(indxmad))*12,:) = SurfFiles(indxmad,:);
        SurfFilest(2:12:(length(indxmad2))*12,:) = SurfFiles(indxmad2,:);
        SurfFilest(3:12:(length(indxovi))*12,:) = SurfFiles(indxovi,:);
        SurfFilest(4:12:(length(indxovi2))*12,:) = SurfFiles(indxovi2,:);
        SurfFilest(5:12:(length(indxbcn))*12,:) = SurfFiles(indxbcn,:);
        SurfFilest(6:12:(length(indxbcn2))*12,:) = SurfFiles(indxbcn2,:);
        SurfFilest(7:12:(length(indxzgz))*12,:) = SurfFiles(indxzgz,:);
        SurfFilest(8:12:(length(indxzgz2))*12,:) = SurfFiles(indxzgz2,:);
        SurfFilest(9:12:(length(indxval))*12,:) = SurfFiles(indxval,:);
        SurfFilest(10:12:(length(indxval2))*12,:) = SurfFiles(indxval2,:);
        SurfFilest(11:12:(length(indxvit))*12,:) = SurfFiles(indxvit,:);
        SurfFilest(12:12:(length(indxvit2))*12,:) = SurfFiles(indxvit2,:);
    else
        %SurfFilest = SurfFiles(1:36,:);
        SurfFilest(1:6:(length(indxmad))*6,:) = SurfFiles(indxmad,:);
        SurfFilest(2:6:(length(indxovi))*6,:) = SurfFiles(indxovi,:);
        SurfFilest(3:6:(length(indxbcn))*6,:) = SurfFiles(indxbcn,:);
        SurfFilest(4:6:(length(indxzgz))*6,:) = SurfFiles(indxzgz,:);
        SurfFilest(5:6:(length(indxval))*6,:) = SurfFiles(indxval,:);
        %if strcmp(estnumber,['-2' filesep])
        SurfFilest(6:6:(length(indxvit))*6,:) = SurfFiles(indxvit,:);
        % else
        % SurfFilest(6:6:end) = SurfFiles(indxvit,:);
        % end
    end
else
    %SurfFilest = SurfFiles(1:36,:);
    SurfFilest(1:6:(length(indxmad))*6,:) = SurfFiles(indxmad,:);
    SurfFilest(2:6:(length(indxovi))*6,:) = SurfFiles(indxovi,:);
    SurfFilest(3:6:(length(indxbcn))*6,:) = SurfFiles(indxbcn,:);
    SurfFilest(4:6:(length(indxzgz))*6,:) = SurfFiles(indxzgz,:);
    SurfFilest(5:6:(length(indxval))*6,:) = SurfFiles(indxval,:);
    %if strcmp(estnumber,['-2' filesep])
    SurfFilest(6:6:(length(indxvit))*6,:) = SurfFiles(indxvit,:);
end
SurfFiles = SurfFilest;
% SurfFilest(indx,:) = [];
% Nsites = 6;
% % Stemp = SurfFiles;
% %  for i =1:Nsites
%
%  end
%     ind = strfind(deblank(SurfFiles(i,:)),estnumber);
%     if isempty(ind)
%         cont = cont+1;
%        indx(cont) =cont;
%     end
% end
return;
