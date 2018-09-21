function Ord10_Compress_MultisCale_Connectivity_Matrix;


connectDir = '/media/Data/PROCESSING_RESULTS/HCP/7-connectome';
IdsFile = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/HCPrest_Ids.txt';

%% ================== Subjects Processing ============================== %%
Ids = char(textread(IdsFile,'%s'));
%% ==== End of Creating Multiscale Atlas for template subject ====== %%
Nsubj = size(Ids,1);
errormessa = '';
for j= 1:Nsubj
    subjId = deblank(Ids(j,:));
    disp(['Processing =======>  Subject ID: ' Id ' . ---  ' num2str(j) ' of ' num2str(Nsubj)]);
    
    opts.pipe.subjId = subjId;
    opts.anat.atlas.type = 'kmeans 832';
    opts.anat.atlas.filename = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/label/rh.aparc_scale05_Nzones0410.annot';
    opts.diff.tracking.method = 'graphmodel';
    opts.diff.tracking.ttype = 'odf';
    
    
    [aparcl,aparcctabl] = read_cfiles('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/label/lh.aparc.annot');
    [aparcr,aparcctabr] = read_cfiles('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/label/rh.aparc.annot');
    
    aparcctabl.struct_names([1 5]) = [];
    aparcctabr.struct_names([1 5]) = [];
    Nl = length(aparcctabl.struct_names);
    Nr = length(aparcctabr.struct_names);
    
    % Connect = Read_Connectivity_Matrix([connectDir filesep 'tractres' filesep 'rumba' filesep subjId '-Connectivity_Matrix-aparc+aseg_RUMBA.mat']);
    load([connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Connectivity_Matrix-aparc+aseg_RUMBA.mat']);
    Connect.StructNames = Connect.Names;
    Connect.StructCodes = Connect.Codes;
    
    load([connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Distance_Matrix-aparc+aseg_RUMBA.mat'])
    Distance.StructNames = Distance.Names;
    Distance.StructCodes = Distance.Codes;
    
    NamesCort = strvcat([repmat('ctx-lh-',[Nl 1]) char(aparcctabl.struct_names)],[repmat('ctx-rh-',[Nr 1]) char(aparcctabr.struct_names)]);
    NamesO= strvcat(Connect.StructNames(1:14,:),NamesCort);
    Names= strvcat(NamesO(1:14,:),NamesCort);
    
    
    %% ====================== Freesurfer parcellation ========================%
    % Subcortical Structures
    SubCL = [10:13 17:18 26];
    SubCR = [ 49:54 58];
    
    % Selecting Frontal regions
    FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032 1026 1002]);
    FroIdsR = FroIds + 1000;
    
    % Selecting Temporal regions
    TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
    TempIdsR = TempIds + 1000;
    
    % Selecting Parietal regions
    ParIds = sort([1029 1008 1031 1022 1025 1010 1023]);
    ParIdsR = ParIds + 1000;
    
    % Selecting Occipital regions
    OccIds = sort([1011 1013 1005 1021]);
    OccIdsR = OccIds + 1000;
    
    % Selecting Insula regions
    InsIds = [1035];
    InsIdsR = [2035];
    OldOrg = [SubCL SubCR sort([FroIds FroIdsR ParIds ParIdsR TempIds TempIdsR OccIds OccIdsR InsIds InsIdsR])];
    Norg = [SubCL FroIds ParIds TempIds OccIds InsIds SubCR FroIdsR ParIdsR TempIdsR OccIdsR InsIdsR];
    
    [~,it] = ismember(Norg,OldOrg);
    Names = Names(it,:);
    NamesO = NamesO(it,:);
    %% ====================== End of Freesurfer parcellation =================%
    
    
    % Reorganicing Connectivity Matrix
    Ns = size(Names,1);
    AllCodes = 0;
    AllNames ='';
    AllInds = 0;
    
    for i = 1:Ns
        nm2find = deblank(Names(i,:));
        ind = find(ismember(Connect.StructNames(:,1:length(nm2find)),nm2find,'rows'));
        if length(ind) >1
            [reordNames,order] = sortrows(Connect.StructNames(ind,:));
            Codes= Connect.StructCodes(ind,:);
            reordCodes = Codes(order);
            ind2reorg = ind(order);
        else
            reordNames = Connect.StructNames(ind,:);
            reordCodes = Connect.StructCodes(ind,:);
            ind2reorg = ind;
        end
        AllCodes = [AllCodes;reordCodes];
        AllNames = [AllNames;reordNames];
        AllInds = [AllInds;ind2reorg];
        
    end
    AllCodes(1) = [];
    AllInds(1) = [];
    
    t = Connect.Matrix(:,AllInds);
    t = t(AllInds,:);
    reordConnect.Matrix = t./max(t(:));
    reordConnect.StructNames = AllNames;
    reordConnect.StructCodes = AllCodes;
    
    t = Distance.Matrix(:,AllInds);
    t = t(AllInds,:);
    reordDistance.Matrix = t;
    reordDistance.StructNames = AllNames;
    reordDistance.StructCodes = AllCodes;
    
    
    opts.diff.connect.matrix = [connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Connectivity_Matrix-aparc+aseg_RUMBA_reord_832Struct.txt'];
    Save_Connectivity_Matrix(opts.diff.connect.matrix,reordConnect,opts);
    
    opts.diff.distance.matrix = [connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Distance_Matrix-aparc+aseg_RUMBA_reord_832Struct.txt'];
    Save_Connectivity_Matrix(opts.diff.distance.matrix,reordDistance,opts);
    
    %  Connect =  reordConnect;
    %  Distance = reordDistance;
    
    
    Nsubcort = 14;
    
    Nscales = size(lhregions,1);
    
    refConnectMatrix = [connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Connectivity_Matrix-aparc+aseg_RUMBA_reord_832Struct.txt'];
    refDistanceMatrix = [connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Distance_Matrix-aparc+aseg_RUMBA_reord_832Struct.txt'];
    Connect = Read_Connectivity_Matrix(refConnectMatrix);
    Distance = Read_Connectivity_Matrix(refDistanceMatrix);
    
    
    
    
    for i = 1:size(AllNames,1)
        
        temp = strread(AllNames(i,:),'%s','delimiter','_');
        scalBelong(i,1:2) = [i length(temp)];
        if i == 1
            Nscales = length(temp);
        else
            Nscales = max([Nscales length(temp)]);
        end
    end
    tempScaleNames = AllNames;
    varTempNames{Nscales} = AllNames;
    varTempCodes{Nscales} = AllCodes;
    for i = Nscales-1:-1:1
        tempVar = tempScaleNames;
        ind = find(scalBelong(:,2) == i+1);
        %     AllNames(ind,:)
        for j = 1:length(ind)
            tempName = tempVar(ind(j),:);
            ind2del = strfind(tempName,'_');
            tempName(ind2del(i):end) = repmat(' ',[1 length(tempName)-ind2del(i)+1]);
            tempVar(ind(j),:) = tempName;
        end
        scalBelong(ind,2) = scalBelong(ind,2) -1;
        [a,b]=unique(tempVar,'rows');
        varTempNames{i} = tempVar(sort(b),:);
        varTempCodes{i} = AllCodes(sort(b),:);
        
        tempScaleNames = tempVar;
    end
    
    for k = 1:Nscales-1
        
        
        origCodes = varTempCodes{k};
        Names = varTempNames{k};
        
        
        ConnectR.StructNames = Names;
        ConnectR.StructCodes = origCodes;
        
        DistanceR.StructNames = Names;
        DistanceR.StructCodes = origCodes;
        
        
        
        % Reorganicing Connectivity Matrix
        Ns = size(ConnectR.StructNames,1);
        AllNames = Connect.StructNames;
        for i = 1:Ns
            nm2find = deblank(Names(i,:));
            ind = find(ismember(AllNames(:,1:length(nm2find)),nm2find,'rows'));
            grouping{i} = ind;
        end
        Ns = length(grouping);
        Mat = zeros(Ns,Ns);
        Matd = zeros(Ns,Ns);
        for i = 1:Ns-1
            xindex = grouping{i};
            for j = i+1:Ns
                jindex = grouping{j};
                [X,Y] = meshgrid(xindex,jindex); X = X(:);Y = Y(:);
                indmat = sub2ind(size(Connect.Matrix),X,Y);
                Mat(i,j) = max(Connect.Matrix(indmat));
                Mat(j,i) = max(Connect.Matrix(indmat));
                
                Matd(i,j) = mean(Distance.Matrix(indmat));
                Matd(j,i) = mean(Distance.Matrix(indmat));
            end
        end
        ConnectR.Matrix = Mat./max(Mat(:));
        DistanceR.Matrix = Matd;
        
        
        %     opts.pipe.subjId = 'HCP_100307-20150807-T1rest';
        opts.anat.atlas.type = ['kmeans ' num2str(Ns)];
        opts.anat.atlas.filename = deblank(lhregions(k,:));
        %     opts.diff.tracking.method = 'probtrack';
        %     opts.diff.tracking.ttype = 'odf';
        
        
        refConnectMatrix = [connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Connectivity_Matrix-aparc+aseg_RUMBA_reord_' num2str(Ns) 'Struct.txt'];
        opts.diff.connect.matrix = refConnectMatrix;
        Save_Connectivity_Matrix(opts.diff.connect.matrix,ConnectR,opts);
        
        refDistanceMatrix = [connectDir filesep subjId filesep 'tractres' filesep 'rumba' filesep subjId '-Distance_Matrix-aparc+aseg_RUMBA_reord_' num2str(Ns) 'Struct.txt'];
        opts.diff.distance.matrix = refDistanceMatrix;
        Save_Connectivity_Matrix(opts.diff.distance.matrix,DistanceR,opts);
        
        clear grouping;
        
    end
end
return;




