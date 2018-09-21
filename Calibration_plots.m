function Calibration_plots(Id, hemi);
Id = ['COMP_00002__MAD-20090606-T1-MAD-1'];%'COMP_00004__VIT-20091017-T1-MAD-2'];
hemi = 'lh';
atlas ='COMP_avesubj_thmax10';
ext = ['.' atlas];
InputDir = 'X:\NEURO\ibas\Compatibilidad_Proyecto_PEPs\Procesamiento_FreeSurfer';
Ns = size(Id,1);
if Ns == 1
    %% Sin corregir
    surfaces = strvcat([InputDir filesep Id filesep 'surf' filesep hemi '.sphere'],[InputDir filesep Id filesep 'surf' filesep hemi '.inflated'],[InputDir filesep Id filesep 'surf' filesep hemi '.pial']);
    [OutFiles, SurfF] = Exp_Surf(surfaces, '0', '','', 'imp','n');
    thick_nat_var = read_cfiles([InputDir filesep Id filesep 'surf' filesep hemi '.thickness']);
    annot_nat_var = read_cfiles([InputDir filesep Id filesep  'label' filesep hemi '.aparc.annot']);
    curv_nat_var = read_cfiles([InputDir filesep Id filesep 'surf' filesep hemi '.curv']);

    
    colordef black;hf = figure('numbertitle','off','Color','black','name','IBASPM Surface Ploting...:  One Subject Stuffs ');
    %%%%%%%%%%%%%%%%%%%%Pial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(3,3,1);Surf = SurfF{3};Surf.Is = curv_nat_var;
        custom_plotsurf(Surf,'hot');h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        axis tight;axis equal; %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(3,3,2);Surf = SurfF{3};Surf.Is = annot_nat_var;
        custom_plotsurf(Surf);h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([90 0]);axis off;
        axis tight;axis equal; %h=title(['Medial View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(3,3,3);Surf = SurfF{3};Surf.Is = thick_nat_var;
        custom_plotsurf(Surf);h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([0 90]);axis off;
        axis tight;axis equal; %h=title(['Top View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        %%%%%%%%%%%%%%%%%%%%%%%%%%Pial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%Inflated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         subplot(3,3,4);Surf = SurfF{2};Surf.Is = curv_nat_var;
        custom_plotsurf(Surf,'hot');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        axis tight;axis equal; %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(3,3,5);Surf = SurfF{2};Surf.Is = annot_nat_var;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([90 0]);axis off;
        axis tight;axis equal; %h=title(['Medial View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(3,3,6);Surf = SurfF{2};Surf.Is = thick_nat_var;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([0 90]);axis off;
        axis tight;axis equal; %h=title(['Top View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        
        %%%%%%%%%%%%%%%%%%%%Inflated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%Sphere%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(3,3,7);Surf = SurfF{1};Surf.Is = curv_nat_var;
        custom_plotsurf(Surf,'hot');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([180 270]);axis off;
        axis tight;axis equal; %h=title(['Bottom View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(3,3,8);Surf = SurfF{1};Surf.Is = annot_nat_var;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([0 90]);axis off;
        axis tight;axis equal; %h=title(['Top View']);camlight;set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(3,3,9);Surf = SurfF{1};Surf.Is = thick_nat_var;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([180 270]);axis off;
        axis tight;axis equal; %h=title(['Bottom View']);camlight;set(h,'FontSize',15,'FontName','Arial');
        camlight;
        
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            mat = max(Surf(1).Is);
            set(h,'YTickLabel',num2str(str2num(get(h,'YTickLabel'))*mat));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%Sphere%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    surfaces = strvcat([InputDir filesep atlas filesep 'surf' filesep hemi '.sphere.reg.avg'],[InputDir filesep atlas filesep 'surf' filesep hemi '.inflated_avg'],[InputDir filesep atlas filesep 'surf' filesep hemi '.pial_avg']);
    [OutFiles, SurfF] = Exp_Surf(surfaces, '0', '','', 'imp','n');
    
    colordef black;hf = figure('numbertitle','off','Color','black','name','IBASPM Surface Ploting...:  One Subject Changes ');
    %%  Sin corregir
        
    thick_reg_var0 =  read_cfiles([InputDir filesep Id filesep 'surf' filesep 'lh.thickness.fwhm0' ext '.mgh']);
    thick_reg_var10 =  read_cfiles([InputDir filesep Id filesep 'surf' filesep 'lh.thickness.fwhm10' ext '.mgh']);
    thick_reg_var20 =  read_cfiles([InputDir filesep Id filesep 'surf' filesep 'lh.thickness.fwhm20' ext '.mgh']);
    %%
    %%Corregido
    load([InputDir filesep 'Calibration_Results' filesep 'lh.thickness.fwhm0' ext '.mat']);
    num = Id2num(Id,size(ResultsTV.NewNcth));thick_reg_var0_corr = ResultsTV.NewNcth(:,num);
    
    load([InputDir filesep 'Calibration_Results' filesep 'lh.thickness.fwhm10' ext '.mat']);
    num = Id2num(Id,size(ResultsTV.NewNcth));thick_reg_var10_corr = ResultsTV.NewNcth(:,num);
    
    load([InputDir filesep 'Calibration_Results' filesep 'lh.thickness.fwhm20' ext '.mat']);
    num = Id2num(Id,size(ResultsTV.NewNcth));thick_reg_var20_corr = ResultsTV.NewNcth(:,num);
    
    %%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%Inflated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,3,1);Surf = SurfF{2};Surf.Is = thick_reg_var0;
        custom_plotsurf(Surf); h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        axis tight;axis equal; %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(2,3,2);Surf = SurfF{2};Surf.Is = thick_reg_var10;
        custom_plotsurf(Surf);h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([90 0]);axis off;
        axis tight;axis equal; %h=title(['Medial View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(2,3,3);Surf = SurfF{2};Surf.Is = thick_reg_var20;
        custom_plotsurf(Surf);h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([0 90]);axis off;
        axis tight;axis equal; %h=title(['Top View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        %%%%%%%%%%%%%%%%%%%%Inflated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%Inflated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         subplot(2,3,4);Surf = SurfF{2};Surf.Is = thick_reg_var0_corr;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        axis tight;axis equal; %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(2,3,5);Surf = SurfF{2};Surf.Is = thick_reg_var10_corr;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([90 0]);axis off;
        axis tight;axis equal; %h=title(['Medial View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        subplot(2,3,6);Surf = SurfF{2};Surf.Is = thick_reg_var20_corr;
        custom_plotsurf(Surf);
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([0 90]);axis off;
        axis tight;axis equal; %h=title(['Top View']);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        
        %%%%%%%%%%%%%%%%%%%%Inflated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%

elseif Ns > 1
    colordef black;hf = figure('numbertitle','off','Color','black','name','IBASPM Surface Ploting...:  Two Subject Changes ');cont =0;
    for i =1:Ns
        surfaces = strvcat([InputDir filesep  deblank(Id(i,:)) filesep 'surf' filesep hemi '.sphere.reg']);%,[InputDir filesep  deblank(Id(i,:)) filesep 'surf' filesep hemi '.inflated'],[InputDir filesep  deblank(Id(i,:)) filesep 'surf' filesep hemi '.pial']);
        thick_nat_var = read_cfiles([InputDir filesep deblank(Id(i,:)) filesep 'surf' filesep hemi '.thickness']);
        annot_nat_var = read_cfiles([InputDir filesep deblank(Id(i,:)) filesep 'label' filesep hemi '.aparc.annot']);
        curv_nat_var = read_cfiles([InputDir filesep deblank(Id(i,:)) filesep 'surf' filesep hemi '.curv']);
        [OutFiles, SurfF] = Exp_Surf(surfaces, '0', '','', 'imp','n');
        cont = cont+1;
        subplot(Ns,3,cont);Surf = SurfF{1};Surf.Is = curv_nat_var;
        custom_plotsurf(Surf,'hot');if i ==1;h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');end
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        axis tight;axis equal; 
        camlight;cont = cont+1;
        subplot(Ns,3,cont);Surf = SurfF{1};Surf.Is = annot_nat_var;
        custom_plotsurf(Surf);if i ==1;h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');end
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        axis tight;axis equal; 
        camlight;cont = cont+1;
        subplot(Ns,3,cont);Surf = SurfF{1};Surf.Is = thick_nat_var;
        custom_plotsurf(Surf); if i ==1;h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');end
        if hemi == 'lh'
            view([270 0]);
        elseif hemi == 'rh'
            view([90 0]);
        end
        axis off;
        %view([0 90]);axis off;
        axis tight;axis equal;
        camlight;

    end
    
end
% if sum(Surf.Is-floor(Surf.Is)) ~=0
%     h = colorbar;
%     mat = max(Surf(1).Is);
%     set(h,'YTickLabel',num2str(str2num(get(h,'YTickLabel'))*mat));
% end
function num = Id2num(Id,N)
if strfind(Id,'COMP_00002__MAD')
    if N == 12
        num = [1 2];
    else
        num = 1;
    end
elseif strfind(Id,'COMP_00003__ZGZ')
     if N == 12
        num = [3 4];
    else
        num = 3;
    end
elseif strfind(Id,'COMP_00004__VIT')
    if N == 12
        num = [5 6];
    else
        num = 5;
    end
elseif strfind(Id,'COMP_00005__BCN')
     if N == 12
        num = [7 8];
    else
        num = 7;
    end
elseif strfind(Id,'COMP_00006__OVI')
    if N == 12
        num = [9 10];
    else
        num = 9;
    end
elseif strfind(Id,'COMP_00007__VAL')
    if N == 12
        num = [11 12];
    else
        num = 11;
    end
end
return

function custom_plotsurf(Surf,cl);
if nargin ==1
    cl = 'jet';
end
[Colors] = Surf_Color(Surf,cl);
Surf.SurfData.FaceVertexCData = Colors;
Surf.SurfData.FaceColor = 'interp';
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
return