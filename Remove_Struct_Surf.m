function [OutputFiles,remSurf] = Remove_Struct_Surf(SurfFiles, Out_dir, StList, Rem, Ref, sa);
%
% Syntax:
%  [OutputFiles] = Remove_Struct_Surf(SurfFiles, Out_dir, StList, Rem, Ref);
%
% This function removes the specified structures from an atlased image.
% The strutures codes must be specified at StList variable.
%
% Input Parameters:
%   SurfFiles    : Individual Atlased Surfaces.
%   Out_dir      : Output Directory to save the resulting surface.
%   StList         : Structures code list. It cotains a list with the structures codes
% i.e - '1,78-92,6' : Removes/lefts structures with the specified codes. Structure
%                  with codes iqual to 1, structures with codes from 78 to 92 and
%                  structure with code iqual to 6.
%                  It offers the chance to choose 2 or more structures with
%                  codes that aren't one after the other, just separating the
%                  structures codes with a coma(,).
%                  If you want to select a range of codes, you have to put the first
%                  structure code, then '-' and next the last structure code,
%                  i.e, 78-92. This expresion removes/lefts for all the
%                  structures with codes from 78 to 92.
%   Rem       : Boolean variable to remove or left the specified
%                  structures. 1(remove) o 0(left);
%  Ref         : Boolean variable to refill the resulting surface. 1(refill) o 0(no refill);
%    sa          : Variable to save or not the matlab surfaces.
% Note: A correct structure selection is very IMPORTANT!!! to select the specified structures.
% 1.- If you want to remove/left individual structures, you must
% separate the structures codes with a coma (ie. 34,67,89,90).
% 2.- If you want to remove/left structures for a range of structures, you must
%  separate the first and the last code with a '-'(i.e, 67-90).
%
% Output Parameters:
%    OutputFiles   : Surfaces Filenames without/with the structures.
%
%
% See also: Atlasing Auto_Labelling
%__________________________________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2006
% Version $1.0

warning off
fclose('all');
%=====================Checking Input Parameters===========================%
if nargin==0
    [SurfFiles,sts] = spm_select([1 Inf],'any','Selecting Atlased Surfaces','',cd);
    sa = input('Do you want to save the surfaces in a new folder(y/n):   ','s');
    if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
    if strcmp(sa,'y')
        [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
    elseif strcmp(sa,'n')
        Out_dir = '';
    end
    StList = input('Please enter the structures codes:  ','s');
    Rem = input(' Select 1 if you want to remove or 0 if you want to left the structures:  ','f');
    Ref = input(' Select 1 if you want to refill the resulting surface or 0 if you do not:  ','f');
else
    if isempty(SurfFiles)
        [SurfFiles,sts] = spm_select([1 Inf],'any','Selecting Atlased Surfaces','',cd);
    end;
    if isempty(StList)
        StList = input('Please enter the structures codes:  ','s');
    end;
    if isempty(Rem)
        Rem = input(' Select 1 if you want to remove or 0 if you want to left the structures:  ','f');
    end;
    if isempty(Ref)
        Ref = input(' Select 1 if you want to refill the resulting surface or 0 if you do not:  ','f');
    end;
    if ~exist('sa','var')|isempty(sa)
        sa = input('Do you want to save the surfaces(y/n):   ','s');
        if (lower(sa)~='y')&(lower(sa)~='n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
        if strcmp(lower(sa),'y')
            [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
        else
            Out_dir = '';
        end
    end
end
%=========================================================================%

%=========================Main program====================================%
Ns = size(SurfFiles,1);
if (size(Out_dir,1)~=Ns)&(size(Out_dir,1)==1)
    Out_dir = repmat(Out_dir,[Ns 1]);
end
StrList = Sel_Struct(StList);
OutputFiles = '';
for i = 1:Ns
    disp(['Case --> ' num2str(i) ' ====>> Removing structures from  surface  '  num2str(i) ' <<====']);
    [pth,nm,ext] = fileparts(SurfFiles(i,:));
    if strcmp(ext(1:4),'.mat')
        Surf = load('-mat',[pth filesep nm ext(1:4)]);
        try Surf = Surf.Surf;
        catch
            warndlg([pthm filesep nmm extm(1:4) '  is does not contain surface information '])
        end
    else
        try
            [OutFiles, SurfF] = Exp_Surf(SurfFiles(i,:), '0', '','', 'imp','n');
            Surf = SurfF{1};
            if exist([pth filesep 'Atlas_' nm '.txt'])
                txt = load([pth filesep 'Atlas_' nm '.txt']);txt = txt(:);
            elseif exist([pth filesep 'Atlased_Surf' filesep 'Atlas_' nm '.txt'])
                txt = load([pth filesep 'Atlased_Surf' filesep 'Atlas_' nm '.txt']);txt = txt(:);
            end
            if size(Surf.SurfData.vertices,1) == size(txt,1)
                Surf.Is =txt;
            else
                errordlg(['The file Atlas_' nm '.txt ' 'is not the atlas file for ' nm ' surface' ]);
                return;
            end
        catch
            warndlg([pthm filesep nmm extm(1:4) '  is does not contain surface information ']);
        end

    end
    I = Surf(i).Is;
    Temp = ismember(I(:)',StrList);
    if Rem ==1
        I(Temp') = 0;
    elseif Rem ==0
        I(~Temp') = 0;
    else
        errordlg('Please select a correct boolean value for removing/lefting structures');
    end
    if Ref ==1
        I = Recur(I,Surf,2);
    elseif (Ref ~=0)&(Ref ~=1)
        errordlg('Please select a correct boolean value for removing/lefting structures');
    end
    if strcmp(sa,'y')
        ind = strfind(Out_dir(i,:), filesep);
        if ind(end)~=length((Out_dir(i,:)))
            fold = lower(Out_dir(i,ind(end)+1:end));
            if strcmp(fold,'joint_surf')|strcmp(fold,'rem_surf')||strcmp(fold,'atlas_surf')
                Output_dir = Out_dir(i,1:ind(end)-1);
            else
                Output_dir = Out_dir(i,:);
            end
         elseif (size(ind)==1)&(ind(end)~=length((Out_dir(i,:))))
            fold = lower(Out_dir(i,ind(end-1)+1:ind(end)-1));
            if strcmp(fold,'joint_surf')|strcmp(fold,'rem_surf')||strcmp(fold,'atlas_surf')
                Output_dir = Out_dir(i,1:ind(end-1)-1);
            else
                Output_dir = Out_dir(i,:);
            end
        else
            Output_dir = Out_dir(i,:);
        end
        mkdir(Output_dir,'Rem_Surf')
        temp = strfind(nm,'Rem_');
        if isempty(temp)
            nm = ['Rem_' nm];
        end
        if strcmp(ext(1:4),'.mat')
            Surf.Is = I;
            [Outfile] = savematfile([Output_dir filesep 'Rem_Surf' filesep nm ext(1:4)], Surf);
            OutputFiles = strvcat(OutputFiles,[Output_dir filesep 'Rem_Surf' filesep nm ext(1:4)]);
            remSurf = '';
        else
            OutputFiles = strvcat(OutputFiles,[Output_dir filesep 'Rem_Surf' filesep nm '.txt']);
            fid = fopen(deblank(OutputFiles(i,:)),'w');
            fprintf(fid,'%4u\n',I);
            fclose(fid);
            remSurf = '';
        end
    elseif strcmp(sa,'n')
        remSurf{i,1} = Surf;
        OutputFiles = '';
    else
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
end
%========================End of main program==============================%
return;