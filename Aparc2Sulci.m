function varargout = Aparc2Sulci(varargin);
%
% Syntax :
%     OutAnnotFile = Aparc2Sulci(Surf, AnnotFile, OutAnnotFile);
%
% This function creates sulci labeles according using the aparc boundaries.
% See
%
% Input Parameters:
%       Surf                   : Surface variable (file, struct or cellarray).
%       AnnotFile              : Annotation File
%       OutAnnotFile           : Saving the sulci parcellation in an
%                                annotation file
%
% Output Parameters:
%        OutAnnotFile           : Saving the sulci parcellation in an
%                                 annotation file. If nargin <3
%                                 OutAnnotFile is a vector file containing
%                                 sulci parcellation.
%
% See also: save_annotfiles
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0


nargin = 2;
varargin{1} = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/fsaverage/surf/lh.inflated';
varargin{2} = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/fsaverage/label/lh.aparc.annot';

%% ============================= Checking Inputs ======================= %%
if nargin < 2
    error('Two Inputs are mandatory');
    return
end
Surf = varargin{1};
AnnotFile = varargin{2};

% Surface Checking
Surf = Surface_Checking(Surf); 

if ~exist(AnnotFile,'file')
    errordlg('The Annotation File Does not exist');
    return;
else
    try
        % Reading Parcellation
        [txt,colortable] = read_cfiles(AnnotFile);
        if ~isfield(colortable,'table')
            errordlg('The File is not an annotation file');
            return;
        end
    catch
        errordlg('Error reading Annotation File');
        return;
    end
end
if nargin == 3
    OutAnnotFile = varargin{3};
end
if nargin > 3
    errordlg('To Many Input Parameters');
    return;
end
if nargout > 1
    errordlg('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%

tabfile = which('aparc.annot.klein2012.sulci.ctab');
if ~exist(tabfile,'file')
    error('Wrong CTab file');
    return;
end

%% ========================= Main Program =============================== %

% --------------------------- Creating Colortable ----------------------- %
% Definning Sulci Names
% 1stsegmentposttempsulcus = ”1st segment of post. sup. temporal sulcus / primary intermediate sulcus”
% 1sttransversetemporalsulcus = “1st transverse temporal sulcus and Heschl’s sulcus_ftts”
% sulcinames = {'frontomarginalsulcus_fms';'superiorfrontalsulcus_sfrs';'inferiorfrontalsulcus_ifrs';'precentralsulcus_prcs';'centralsulcus_cs';'postcentralsulcus_pocs';...
%     'intraparietalsulcus_itps';'1stsegmentposttempsulcus_csts1';'sylvianfissure_ls';'lateraloccipitalsulcus_locs';'anterioroccipitalsulcus_aocs';'superiortemporalsulcus_sts';...
%     'inferiortemporalsulcus_its';'circularsulcus_crs';'cingulatesulcus_cgs'; 'paracentralsulcus_pcs';'parietooccipitalfissure_pos';...
%     'calcarine fissure_ccs';'superiorrostralsulcus_sros';'callosalsulcus_cas';'lateralHshapedorbitalsulcus_lhos';'olfactorysulcus_olfs';'occipitotemporalsulcus_ots';'collateralsulcus_cos'};
% 
% % Definning Sulci Colors
% colors = [240,163,255;...
%     0,117,220;...
%     153,63,0;...
%     76,0,92;...
%     0,92,49;...
%     43,206,72;...
%     255,204,153;...
%     148,255,181;...
%     143,124,0;...
%     157,204,0;...
%     194,0,136;...
%     210 113 26;...
%     0,51,128;...
%     255,164,5;...
%     255,168,187;...
%     66,102,0;...
%     255,0,16;...
%     94,241,242;...
%     0,153,143;...
%     224,255,102;...
%     116,10,255;...
%     153,0,0;...
%     255,255,128;...
%     255,80,5];

[sulcinames, r,g,b,tr] = textread(tabfile,'%s%u%u%u%u','delimiter',' ');
colors = [r g b];

%%
% 
% 
colortab = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];

% Creating Colortable
colortable.numEntries = size(colortab,1);
colortable.orig_tab = 'Sulci Colortable';
colortable.struct_names = sulcinames;
colortable.table = colortab;         
% ---------------------- End of Creating Colortable --------------------- %

% --------------------- Reading Input Annot File ------------------------ %
[txt,ctab] = read_cfiles(AnnotFile);
% ctab.table = [ctab.table [0:size(ctab.table,1)-1]' ];
% txt_relabel = txt*0;
% for i = 1:length(ctab.table(:,5))
%     ind = find(txt == ctab.table(i,5));
%     txt_relabel(ind) = i-1;
% end
[a,txt_relabel] = ismember(txt,ctab.table(:,5));
txt_relabel = txt_relabel-1;

% ------------------ End of Reading Input Annot File -------------------- %

% ------------------ Computing Regional Boundaries ---------------------- %
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;
temp1 = txt_relabel(temp);
temp1(indz) =  max(temp1(:))+1;
NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')');
indc = find(a);
gCrowns = zeros(length(txt_relabel),1);
gCrowns(indc) = 1;
gCrowns = txt_relabel.*gCrowns;
% ----------------------- End of Extracting Crown points ---------------- %


% ---------------- Computing Crossing Edges ----------------------------- %
indcrown = find(gCrowns);
indcrownfaces = find(sum(ismember(Surf.SurfData.faces,indcrown),2) ~= 0);
crownfaces = Surf.SurfData.faces(indcrownfaces,:);
crossEdges = [crownfaces(:,1) crownfaces(:,2); crownfaces(:,2) crownfaces(:,3);crownfaces(:,1) crownfaces(:,3)] ; % Edges from the intersection faces

% ------------------- End of Computing Crossing Edges ------------------- %

% ----------------------- Labelling Crossing Edges ---------------------- %
labCrossEdges = sort(txt_relabel(crossEdges)')';
% -------------------- End of Labelling Crossing Edges ------------------ %

% --------------------------- Sulci definition  ------------------------- %            
sulcpairmat{1} = [28 32];
sulcpairmat{2} = [3 28;27 28];
sulcpairmat{3} = [3 18 ;3 19; 3 20; 18 27;19 27;20 27];
sulcpairmat{4} = [24 28;3 24;24 27;18 24;19 24;20 24];
sulcpairmat{5} = [22 24];
sulcpairmat{6} = [22 29 ;22 31];
sulcpairmat{7} = [29 31;8 29];
sulcpairmat{8} = [8 31];
sulcpairmat{9} = [30 31];
sulcpairmat{10} = [8 11;11 29];
sulcpairmat{11} = [11 15; 9 11];
sulcpairmat{12} = [15 30;1 30];
sulcpairmat{13} = [9 15];
sulcpairmat{14} = [[12,35]; [30,35]; [34,35]; [2,35];[10,35];[23,35];[26,35];[22,35]; [24,35]; [31,35]; [18,35]; [19,35]; [20,35]];
% sulcpairmat{15} = [30 34];
sulcpairmat{15} = [[2,14]; [10,14]; [14,23];  [2,28]; [10,28]; [23,28];  [2,17]; [10,17]; [17,23]; [17,26]; [17,25]];%[26,28];[14,26]; 
sulcpairmat{16} = [17 28];
sulcpairmat{17} = [5 25;[5 13];[21 25]];
sulcpairmat{18} = [[13,25]; [2,13]; [10,13]; [13,23]; [13,26];[13 21]];
sulcpairmat{19} = [14 28;[14,26];28 32];% Remover [14,26];
sulcpairmat{20} = [[2,4]; [4,10]; [4,23]; [4,26]];
sulcpairmat{21} = [[3,12]; [12,27]; [12,18]; [12,19]; [12,20]];
sulcpairmat{22} = [12 14];
sulcpairmat{23} = [[7,9]; [7,11]];
sulcpairmat{24} = [[6,7]; [7,16]; [7,13]];
% ---------------------- End of Sulci definition  ----------------------- % 

% --------------------------- Sulci Estimation  ------------------------- %
sulcAnnot = zeros(length(txt),1);
ind2rem = '';
cont = 0;
for  i = 1:length(sulcpairmat)
    %     indsulci = find(sum(ismember(labCrossEdges,sulcpairmat{i}),2) == 2);
    indsulci = find(ismember(labCrossEdges,sulcpairmat{i},'rows'));
    if isempty(indsulci)
        cont = cont + 1;
        ind2rem(cont) = i;
    end
    temp = crossEdges(indsulci);
    sulcAnnot(temp(:)) = colortable.table(i,5);
end
colortable.struct_names(ind2rem) = [];
colortable.table(ind2rem,:) = [];

% Adding 0
colortable.table = [0 0 0 0 0;colortable.table];
colortable.struct_names = [{'unknown'};colortable.struct_names];
colortable.numEntries = size(colortable.table,1);

[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;
% ================== Creating the Labels Matrix ===================== %
% This matrix contains information about the
% labels inside each point neighborhood

Labels_Mat = sulcAnnot(temp);
Labels_Mat(indz) = 0;
% =============== End of Creating the Labels Matrix ================= %

% Dilating clusters 
dilLabels = sulcAnnot;
for it = 1:2
    vert_to_grow = find(sum(Labels_Mat,2) ~= 0); % Finding Labeled Neighbors
    vert_to_grow(ismember(vert_to_grow,find(dilLabels~=0))) = [];% Vertex to grow
    
   
    % Labeled with the label of the majority of its labeled neighbors
    dilLabels(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
    
    dilLabels(vert_to_grow) = arrayfun( @(vert_to_grow)mode(nonzeros(Labels_Mat(vert_to_grow,:))), vert_to_grow);
    
    % Updating the labels matrix
    Labels_Mat = sulcAnnot(temp);
    Labels_Mat(indz) = 0;
    
    % Detecting if someone is missing for labelling
    indn = find(dilLabels == 0);
    if isempty(indn)&(it ~= opts.iterations)
        break;
    end
end
sulcAnnot = dilLabels;
% ------------------------- End of Sulci Estimation  -------------------- % 

% ---------------- Computing Crossing Edges ----------------------------- %
indcrown = find(sulcAnnot);
indcrownfaces = find(sum(ismember(Surf.SurfData.faces,indcrown),2) ~= 0);
crownfaces = Surf.SurfData.faces(indcrownfaces,:);
crossEdges = [crownfaces(:,1) crownfaces(:,2); crownfaces(:,2) crownfaces(:,3);crownfaces(:,1) crownfaces(:,3)] ; % Edges from the intersection faces
labCrossEdges = sulcAnnot(crossEdges);
[X, Y] = find(labCrossEdges == 0);
crossEdges(X,:)  = [];
ind2rem = find(sulcAnnot(crossEdges(:,1)) - sulcAnnot(crossEdges(:,2)) ~=0);
vert2rem = crossEdges(ind2rem,:);
sulcAnnot(vert2rem(:)) = 0; 

% ------------------- End of Computing Crossing Edges ------------------- %

% Medial Wall Indexes
% indmwall = find(txt == 0);
% sulcAnnot(indmwall) = colortable.table(1,5);

%% ========================= End of Main Program ======================== %
if nargin <= 2;
    [a, sulcAnnot] = ismember(sulcAnnot,colortable.table(:,5));
    varargout{1} = sulcAnnot;
else
    OutAnnotFile = save_annotfiles(sulcAnnot,OutAnnotFile,colortable);
    varargout{1} = OutAnnotFile;
end
return
    

 
