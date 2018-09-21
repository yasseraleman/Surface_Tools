function Multi_Run_Euclidean_Distance_betweenPial_and_Outer(Idfile);

% Idfile = strvcat('HCP_201111-20140807-T1wMPR1');
Idfile = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/HCPrest_Ids.txt';
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
Ids = 'HCP_857263-20150807-T1rest';


Ns = size(Ids,1);
failed = '';
FreeSDir = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing';
cad = '';
Ido = '';
for  i = 1:Ns
    
    Id = deblank(Ids(i,:));
    disp(['Processing Subject ' Id '. Subject ' num2str(i) ' of ' num2str(Ns)]);
    try
        delete([ FreeSDir filesep Id filesep  'surf' filesep Id 'SULCLINES_FAILED.txt']);
    end
    %     dirl = dir([FreeSDir filesep Id filesep 'surf' filesep 'lh.lines.crowns.depth.mat']);
    %     dirr = dir([FreeSDir filesep Id filesep 'surf' filesep 'rh.lines.crowns.depth.mat']);
    %     try
    %         if (dirl.datenum < 736174)|(dirr.datenum < 736174);
    try
        %% ======================== Left Hemisphere ============================= %
        pialfile = [FreeSDir filesep Id filesep 'surf' filesep 'lh.pial'];
        pialfileouter = [FreeSDir filesep Id filesep 'surf' filesep 'lh.pial-outer-smoothed'];
        Surf = Read_Surface(pialfile);
        [distMap] = Distance_between_Surfaces(pialfile, pialfileouter);
        fname = [FreeSDir filesep Id filesep 'surf' filesep 'lh.pial_eucdepth' ];
        save_char(fname, distMap, size(Surf.SurfData.faces,1));
        
        %% ======================== Right Hemisphere ============================ %
        
        pialfile = [FreeSDir filesep Id filesep 'surf' filesep 'rh.pial'];
        pialfileouter = [FreeSDir filesep Id filesep 'surf' filesep 'rh.pial-outer-smoothed'];
        Surf = Read_Surface(pialfile);
        [distMap] = Distance_between_Surfaces(pialfile, pialfileouter);
        fname = [FreeSDir filesep Id filesep 'surf' filesep 'rh.pial_eucdepth' ];
        save_char(fname, distMap, size(Surf.SurfData.faces,1));
    catch
        fid = fopen([ FreeSDir filesep Id filesep  'surf' filesep Id 'SULCLINES_FAILED.txt'],'wt');
        fclose(fid);
    end
    %         end
    %     catch
    %         outfiles = FreeSurfer_Gyral_Sulcal_Depth_Lines_Extraction(FreeSDir, Id);
    %     end
    disp(' ');
    
end
return;


