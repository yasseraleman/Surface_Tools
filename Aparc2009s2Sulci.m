function varargout = Aparc2009s2Sulci(varargin);
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


% % nargin = 2;
% % varargin{1} = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/fsaverage/surf/lh.inflated';
% % varargin{2} = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/fsaverage/label/lh.aparc.a2009s.annot';

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
tabfile = which('aparc.annot.a2009s.sulci.ctab');
if ~exist(tabfile,'file')
    error('Wrong aparc2009s CTab file');
    return;
end


%% ========================= Main Program =============================== %

% --------------------------- Creating Colortable ----------------------- %
% Definning Sulci Names
% 1stsegmentposttempsulcus = ”1st segment of post. sup. temporal sulcus / primary intermediate sulcus”
% 1sttransversetemporalsulcus = “1st transverse temporal sulcus and Heschl’s sulcus_ftts”

[sulcinames, r,g,b,tr] = textread(tabfile,'%s%u%u%u%u','delimiter',' ');


% Definning Sulci Colors
colors = [r g b];
colortab = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];

% Creating Colortable
colortable.numEntries = size(colortab,1);
colortable.orig_tab = 'Sulci aparc2009s Colortable';
colortable.struct_names = sulcinames;
colortable.table = colortab;         
% ---------------------- End of Creating Colortable --------------------- %

% --------------------- Reading Input Annot File ------------------------ %
[sulcAnnot,ctab] = read_cfiles(AnnotFile);

[a,txt_relabel] = ismember(sulcAnnot,colortab(:,5));
ind2rem = find(a == 0);
sulcAnnot(ind2rem) = 0;

% Joinning Cingulate
ind = find(ismember(sulcAnnot,colortable.table(1:4,5)));
sulcAnnot(ind) = colortable.table(1,5);
colortable.struct_names{1} = 'S_cingulate';

% Joinning Intraparietal
ind = find(ismember(sulcAnnot,colortable.table(5:6,5)));
sulcAnnot(ind) = colortable.table(5,5);
colortable.struct_names{5} = 'S_intraparietal';

% Joinning Precentral
ind = find(ismember(sulcAnnot,colortable.table(7:8,5)));
sulcAnnot(ind) = colortable.table(7,5);
colortable.struct_names{7} = 'S_precentral';

% Removing labels and names
colortable.table([2 3 4 6 8],:) = [];
colortable.struct_names([2 3 4 6 8]) = [];

% Adding 0
colortable.table = [0 0 0 0 0;colortable.table];
colortable.struct_names = [{'unknown'};colortable.struct_names];

ind2rem = find(ismember(colortable.table(:,5),sulcAnnot) == 0);
colortable.table(ind2rem,:) = [];
colortable.struct_names(ind2rem) = [];

colortable.numEntries = size(colortable.table,1);

% Labels Refilling
Surf.Is = sulcAnnot;
indzeros = find(Surf.Is == 0);
Surf.Is(indzeros) = max(Surf.Is) + 1;
[Surf] = Surf_Corr(Surf);
indzeros = find(Surf.Is == max(Surf.Is));
Surf.Is(indzeros) = 0;
sulcAnnot = Surf.Is;
%% ========================= End of Main Program ======================== %
if nargin <= 2;
    [a, sulcAnnot] = ismember(sulcAnnot,colortable.table(:,5));
    varargout{1} = sulcAnnot;
else
    OutAnnotFile = save_annotfiles(sulcAnnot,OutAnnotFile,colortable);
    varargout{1} = OutAnnotFile;
end
return
    

 
