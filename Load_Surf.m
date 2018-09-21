function [Surfa] = Load_Surf(SurfFiles);
%
% Syntax :
% [Surfa] = Load_Surf(SurfFiles);
%
% Loads surfaces files into a Matlab variable SurfF
%
% Input Parameters:
%   SurfFiles   : Surfaces filenames.
%
% Output Parameters:
%  Surfa        : Output Matlab variable

% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0
warning off;
%=====================Checking Input Parameters===========================%
if nargin==0
    [SurfFiles,sts] = spm_select([1],'any','Selecting Surface Files','',cd);
end
if ~exist('SurfFiles','var')|isempty(SurfFiles)
    [SurfFiles,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
end
%=========================================================================%
%=========================Main program====================================%
Ns = size(SurfFiles,1);
for i = 1:Ns
    [pth,nm,ext] = fileparts(deblank(SurfFiles(i,:)));
    if ~isempty(pth)
        Sfile = [pth filesep nm ext(1:4)];
    else
        Sfile = [nm ext(1:4)];
    end
    form = lower(ext(2:4));
    switch form
        case 'srx'
            Surf.Imp = 'srx';
            Surf.Name = nm;
            Surf.Area = 0;
            fid = fopen(Sfile, 'r');
            Npoints = fread(fid, 1, 'int32');
            Nfaces = fread(fid, 1, 'int32');

            % Reading Vertices
            xyzn = fread(fid, 3*Npoints, 'float32');

            % Reading Faces
            face = fread(fid, 3*Nfaces, 'int32');
            fclose(fid);
            vert = reshape(xyzn, [3 Npoints])';Surf.SurfData.vertices=vert;
            face = reshape(face, [3 Nfaces])';face = face+1;
            Surf.SurfData.faces = face;clear face;
            dx = ceil(max(vert(:,1))-min(vert(:,1)))+2;
            dy = ceil(max(vert(:,2))-min(vert(:,2)))+2;
            dz = ceil(max(vert(:,3))-min(vert(:,3)))+2;
            Surf.Orig = [abs(min(vert(:,1)))+ abs(dx/2)+1 abs(min(vert(:,2)))+ abs(dy/2)+1 abs(min(vert(:,3)))+ abs(dz/2)+1];
            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
            AT = [1 0 0 1; 0 0 -1 dy; 0 1 0 1; 0 0 0 1];
            vert = [vert ones(Npoints,1)]*AT'; vert(:,4) = [];
            Surf.VoxSize = [1 1 1];
            Surf.SurfData.vertices(:,1) = [];
            [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf.Tri = Tri; clear Tri;
            h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
            norma = sqrt(sum((Normals').^2));
            Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
            Surf.Type = 'Mask';
            Surfa{i,1} = Surf;
        case 'txt'
            Surf.Imp = 'txt';
            Surf.Name = nm;
            Surf.Area = 0;
            fid = fopen(Sfile, 'r');

            % Reading Vertices
            Npoints = textread(Sfile,'%f',1,'headerlines',0);
            vert = reshape(textread(Sfile,'%f',Npoints*4,'headerlines',1),4,Npoints)';
            Surf.SurfData.vertices=vert;
            % Reading Faces
            Nfaces = textread(Sfile,'%f',1,'headerlines',Npoints+1);
            face = reshape(textread(Sfile,'%f',Nfaces*4,'headerlines',Npoints+2),4,Nfaces)'+1;
            face(:,1) = [];vert(:,1) = [];
            Surf.SurfData.faces = face;clear face;
            dx = ceil(max(vert(:,1))-min(vert(:,1)))+2;
            dy = ceil(max(vert(:,2))-min(vert(:,2)))+2;
            dz = ceil(max(vert(:,3))-min(vert(:,3)))+2;
            Surf.Orig = [abs(min(vert(:,1)))+ abs(dx/2)+1 abs(min(vert(:,2)))+ abs(dy/2)+1 abs(min(vert(:,3)))+ abs(dz/2)+1];
            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
            Surf.VoxSize = [1 1 1];
            [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
            Surf.SurfData.vertices(:,1) = [];
            % Computing Normals
            h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
            norma = sqrt(sum((Normals').^2));
            Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf.Tri = Tri; clear Tri;
            Surf.Type = 'Mask';
            Surfa{i,1} = Surf;
        case 'off'
            Surf.Imp = 'off';
            fid = fopen(deblank(Sfile),'r');
            lin = fgetl(fid);
            lin = fgetl(fid);
            [Npoints,Nfaces,Car] = strread(lin,'%f%f%f','delimiter',' ');

            % Reading Vertices
            Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',2),3,Npoints)';
            xyzn=Surf.SurfData.vertices;
            % Reading faces
            Surf.SurfData.faces = uint32(reshape(textread(Sfile,'%f',Nfaces*4,'headerlines',Npoints+2),4,Nfaces)')+1;
            Surf.SurfData.faces(:,1) = [];
            dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
            dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
            dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
            Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
            Surf.VoxSize = [1 1 1];
            [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);

            % Computing Normals
            h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
            norma = sqrt(sum((Normals').^2));
            Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf.Tri = Tri;
            clear Tri;
            Surf.Name = nm;
            Surf.Area = 0;
            Surf.VoxSize = [1 1 1];
            Surf.Type = 'Mask';
            Surfa{i,1} = Surf;
        case 'obj'
            Surf.Imp = 'obj';
            fid = fopen(deblank(SurfF(i,:)),'r');
            lin = fgetl(fid);
            [typ,ac,dif,spr,spc,tr,Npoints] = strread(lin,'%s%n%n%n%u%n%u','delimiter',' ');
            if lower(char(typ)) == 'p';

                % Reading Vertices
                Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',1),3,Npoints)';

                % Reading Normals
                Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',Npoints+2),3,Npoints)';
                Surf.SurfData.FaceColor = 'interp';

                % Reading Faces
                Nfaces = textread(Sfile,'%f',1,'headerlines',Npoints*2+3);
                cl_flag=int32(textread(Sfile,'%f',1,'headerlines',Npoints*2+4));
                if cl_flag == 0
                    cl=int32(textread(Sfile,'%f',5,'headerlines',Npoints*2+4));
                    end_ind=textread(Sfile,'%f',Nfaces,'headerlines',Npoints*2+6);
                    temps = repmat(Nfaces,[1 8]);modp = mod(temps,[2:9]);
                    ind = find(modp==0);
                    if ~isempty(ind)
                        Nl = max(ind)+1;
                    else
                        pr = primes(Nfaces);
                        temps = repmat(Nfaces,[1 size(pr,1)]);
                        modp = mod(temps,pr);
                        ind = find(modp==0);Nl = pr(max(ind));
                    end
                    leng = max(end_ind);
                    ind=textread(Sfile,'%f',double(leng),'headerlines',Npoints*2+6+Nfaces/Nl);
                    ind = ind+1;
                    Surf.SurfData.FaceColor = [1 1 1];
                elseif cl_flag == 1
                elseif cl_flag == 2
                    temp = textread(Sfile,'%f',Npoints*4+1,'headerlines',Npoints*2+4); temp(1) = [];temp(4:4:end) = [];
                    Surf.SurfData.FaceVertexCData = reshape(temp,3,Npoints)';clear temp;
                    end_ind=textread(Sfile,'%f',Nfaces,'headerlines',Npoints*3+5);
                    temps = repmat(Nfaces,[1 8]);modp = mod(temps,[2:9]);
                    ind = find(modp==0);
                    if ~isempty(ind)
                        Nl = max(ind)+1;
                    else
                        pr = primes(Nfaces);
                        temps = repmat(Nfaces,[1 size(pr,1)]);
                        modp = mod(temps,pr);
                        ind = find(modp==0);Nl = pr(min(ind));
                    end
                    leng = max(end_ind);
                    ind=textread(Sfile,'%f',double(leng),'headerlines',Npoints*3+5+Nfaces/Nl);
                    ind = ind+1;
                end
                faces =  [ind(end_ind-2) ind(end_ind-1) ind(end_ind)];
                dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
                dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
                dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
                Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
                Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                Surf.VoxSize = [1 1 1];
                Surf.SurfData.faces = faces; clear faces;
                clear ind end_ind;
                [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                Temp = sum(Tri);
                Tri(:,Temp==0) = [];
                Surf.Tri = Tri;
                clear Tri;
                Surf.Name = nm;
                Surf.Area = 0;
                Surf.Type = 'Mask';
                Surfa{i,1} = Surf;
            end
        case 'dfs'
            Surf.Imp = 'dfs';
            Sfile = [pth filesep nm ext(1:4)];
            fid=fopen(deblank(Sfile),'rb','ieee-le');
            magic=char(fread(fid,8,'char'));
            version=fread(fid,4,'char');
            hdrsize=fread(fid,1,'int32');
            mdoffset=fread(fid,1,'int32');
            pdoffset=fread(fid,1,'int32');
            Nfaces=fread(fid,1,'int32');
            Npoints=fread(fid,1,'int32');
            nStrips=fread(fid,1,'int32');
            stripSize=fread(fid,1,'int32');
            Normals=fread(fid,1,'int32');
            uvStart=fread(fid,1,'int32');
            vcoffset=fread(fid,1,'int32');
            precision=fread(fid,1,'int32');
            orientation=fread(fid,[4 4],'float64');
            fseek(fid,hdrsize,-1);

            % Reading Faces
            faces = fread(fid,[3 Nfaces],'int32')+1;

            % Reading Vertices
            vert=fread(fid,[3 Npoints],'float32');
            Surf.SurfData.faces = faces' ;clear faces;
            Surf.SurfData.vertices = vert' ;
            dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
            dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
            dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
            Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
            Surf.VoxSize = [1 1 1];
            [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf.Tri = Tri;
            Surf.Name = nm;
            Surf.Area = 0;
            Surf.VoxSize = [1 1 1];
            Surf.Type = 'Mask';

            % Reading Normals
            if (Normals>0)
                fseek(fid,Normals,-1);
                Surf.SurfData.VertexNormals = fread(fid,[3 Npoints],'float32')';
            end;
            Surfa{i,1} = Surf;
        case 'mes'
            Surf.Imp = 'mesh';
            Sfile = [pth filesep nm ext(1:5)];
            fid = fopen(deblank(Sfile),'rb');
            Type = char(fread(fid, 5, 'uchar'));  %- 'ascii' or 'binar'
            if strcmp(Type','binar');
                [byte_swapping, COUNT]     = fread(fid, 1, 'uint32'); %- 'ABCD' or 'DCBA'
                ff = strcmp(dec2hex(byte_swapping),'41424344');
                if ~ff
                    [fn, pm, mf] = fopen(1); %- machine format
                    fclose(fid);
                    if strmatch(mf,'ieee-le');
                        fid = fopen(deblank(Sfile),'r','ieee-be');
                    else
                        fid = fopen(deblank(Sfile),'r','ieee-le');
                    end
                    [file_format, COUNT]   = fread(fid, 5, 'uchar');
                    [byte_swapping, COUNT] = fread(fid, 1, 'uint32');
                end
                Val = fread(fid,1,'uint32');
                Text = char(fread(fid,4,'char'));
                Pol = fread(fid,1,'uint32');
                Tsteps = fread(fid,1,'uint32');
                for t=1:Tsteps
                    Tinst = fread(fid,1,'uint32');

                    % Reading vertices
                    Npoints = fread(fid,1,'uint32');
                    vert = fread(fid,Npoints*3,'float32')';
                    vert = reshape(vert,[3,Npoints])';
                    Surf(t).SurfData.vertices =vert;
                    % Reading normals
                    T = fread(fid,1,'uint32');
                    normals = fread(fid,Npoints*3,'float32');
                    normals = reshape(normals,[3,Npoints])';
                    normals = Mat*[normals ones(Npoints,1)]';
                    normals = normals';normals(:,4) = [];

                    % Reading Faces
                    T = fread(fid,1,'uint32');
                    Nfaces = fread(fid,1,'uint32');
                    faces = fread(fid,Nfaces*Pol,'uint32');
                    faces = reshape(faces,[Pol,Nfaces])'+1;

                    Surf(t).SurfData.VertexNormals = normals;clear normals
                    Surf(t).SurfData.faces = faces;clear faces vert;
                    dx = ceil(max(Surf(1).SurfData.vertices(:,1))-min(Surf(1).SurfData.vertices(:,1)))+2;
                    dy = ceil(max(Surf(1).SurfData.vertices(:,2))-min(Surf(1).SurfData.vertices(:,2)))+2;
                    dz = ceil(max(Surf(1).SurfData.vertices(:,3))-min(Surf(1).SurfData.vertices(:,3)))+2;
                    Surf(t).Orig = [abs(min(Surf(1).SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf(1).SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf(1).SurfData.vertices(:,3)))+ abs(dz/2)+1];
                    Surf(t).Dim = [abs(dx) abs(dy) abs(dz)];
                    Surf(t).VoxSize = [1 1 1];

                    [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf(t).Tri = Tri;
                    if Tsteps~=1
                        Surf(t).Name = [nm '_' sprintf('%.3d',t)];
                    else
                        Surf(t).Name = nm;
                    end
                    Surf(t).Area = 0;
                    Surf(t).Type = 'Mask';
                end
            elseif strcmp(Type(1:5,1)','ascii');
                Surf(t).Name = nm;
                Text = fscanf(fid,'%s',1);
                Pol = fscanf(fid,'%d',1);
                Tsteps = fscanf(fid,'%d',1);
                for t=1:Tsteps
                    ms = fscanf(fid,'\n%d',1);

                    % Reading vertices
                    Npoints = fscanf(fid,'\n%d\n',1);
                    vert = fscanf(fid,'(%f ,%f ,%f) ',3*vertex_number);
                    vert = reshape(vert,[3,Npoints])';

                    % Reading Normals
                    T = fscanf(fid,'\n%d\n',1);
                    normals = fscanf(fid,'(%f ,%f ,%f) ',3*Npoints);
                    normals = reshape(normals,[3,Npoints])';
                    normals = Mat*[normals ones(Npoints,1)]';
                    normals = normals';normals(:,4) = [];
                    T = fscanf(fid,'\n%d\n',1);

                    % Reading Faces
                    Nfaces = fscanf(fid,'\n%d\n',1);
                    faces = fscanf(fid,'(%d ,%d ,%d) ',Pol*Nfaces);
                    faces = reshape(faces,[Pol,Nfaces])'+1;
                    Surf(t).SurfData.vertices = vert;
                    Surf(t).SurfData.VertexNormals = normals;clear normals
                    Surf(t).SurfData.faces = faces;clear faces vert;
                    dx = ceil(max(Surf(1).SurfData.vertices(:,1))-min(Surf(1).SurfData.vertices(:,1)))+2;
                    dy = ceil(max(Surf(1).SurfData.vertices(:,2))-min(Surf(1).SurfData.vertices(:,2)))+2;
                    dz = ceil(max(Surf(1).SurfData.vertices(:,3))-min(Surf(1).SurfData.vertices(:,3)))+2;
                    Surf(t).Orig = [abs(min(Surf(1).SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf(1).SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf(1).SurfData.vertices(:,3)))+ abs(dz/2)+1];
                    Surf(t).Dim = [abs(dx) abs(dy) abs(dz)];
                    Surf(t).VoxSize = [1 1 1];
                    [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf(t).Tri = Tri;
                    Surf(t).Area = 0;
                    Surf(t).Type = 'Mask';
                end
            end
            Surfa{i,1} = Surf;
            fclose(fid);
        case 'asc'
            Surf.Name = nm;
            fid = fopen(Sfile, 'rt');
            line = fgetl(fid);
            line = fgetl(fid);
            [Npoints,Nfaces] = strread(line,'%u%u','delimiter',' ');

            % Reading Vertices
            Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*4,'headerlines',2),4,Npoints)';
            Surf.SurfData.vertices(:,4) = [];

            % Reading Faces
            Surf.SurfData.faces = reshape(textread(Sfile,'%u',Nfaces*4,'headerlines',2+Npoints),4,Nfaces)'+1;
            Surf.SurfData.faces(:,4) = [];
            fclose(fid);
            dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
            dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
            dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
            Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
            Surf.VoxSize = [1 1 1];
            [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf.Tri = Tri; clear Tri;
            h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
            norma = sqrt(sum((Normals').^2));
            Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
            Surf.Area = 0;
            Surf.Imp = 'asc';
            Surf.Type = 'Mask';
            Surfa{i,1} = Surf;
        case 'vtk'
            Surf.Name = nm;
            fid = fopen(Sfile, 'rt');
            line = fgetl(fid);line = fgetl(fid);line = fgetl(fid);
            line = fgetl(fid);
            [Info,typ] = strread(line,'%s%s','delimiter',' ');
            if ~strcmp(lower(typ),'polydata')
                errordlg('Please select a correct surface format');
                return;
            end
            line = fgetl(fid);
            % Reading Vertices
            [txt,Npoints,typ] = strread(line,'%s%n%s','delimiter',' ');
            Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',5),3,Npoints)';
            dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
            dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
            dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
            Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
            Surf.VoxSize = [1 1 1];
            % Reading Faces
            [line,Nfaces] = textread(Sfile,'%s %u',1,'headerlines',5+Npoints);
            pol = textread(Sfile,'%u',1,'headerlines',6+Npoints);pol = pol+1;
            Surf.SurfData.faces = reshape(textread(Sfile,'%f',Nfaces*pol,'headerlines',6+Npoints),pol,Nfaces)';
            Surf.SurfData.faces(:,1)=[]; Surf.SurfData.faces= Surf.SurfData.faces+1;
            Surf.Is = textread(Sfile,'%f',Npoints,'headerlines',9+Npoints+Nfaces);
            if ~isempty(Surf.Is)
                Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',10+2*Npoints+Nfaces),3,Npoints)';
            else
                Surf =rmfield(Surf,'Is');
            end
            [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf.Tri = Tri; clear Tri;
            Surf.Area = 0;
            Surf.Type = 'Mask';
            Surf.Imp = 'vtk';
            Surfa{i,1} = Surf;
        otherwise
            Tmn =  16777214 ;
            Qmn =  16777215 ;
            [pth,nm,ext] = fileparts(deblank(SurfF(i,:)));
            fid = fopen(deblank(SurfF(i,:)),'rb','b');
            if fid >0
                MNumber = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
                if(MNumber == Qmn)
                    Surf.Imp = 'frb';
                    Npoints = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
                    Nfaces = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
                    Surf.SurfData.vertices = fread(fid, Npoints*3, 'int16') ./ 100 ;
                    for t=1:Npoints
                        for k=1:4
                            Surf.SurfData.faces(t,k) =  bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar'); ;
                        end
                    end
                    Surf.SurfData.faces = Surf.SurfData.faces+1;
                elseif (MNumber == Tmn)
                    Surf.Imp = 'frb';
                    tline = fgets(fid);
                    tline = fgets(fid);
                    Np = fread(fid,1,'int32');
                    Nf = fread(fid,1,'int32');
                    vertices = fread(fid,Np*3,'float32');
                    Surf.Name = [nm ext];
                    Surf.SurfData.vertices = reshape(vertices,3,Np)';
                    faces = fread(fid,Nf*3,'int32');
                    Surf.SurfData.faces = reshape(faces,3,Nf)'+1;
                    fclose(fid);
                else
                    errordlg('This file is not a known Surface File. Please try again');
                    return
                end
                dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
                dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
                dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
                Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
                Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                Surf.VoxSize = [1 1 1];
                [Tri] = Vert_Neib(double(Surf.SurfData.faces),Np,Nf);
                Temp = sum(Tri);
                Tri(:,Temp==0) = [];
                Surf.Tri = Tri;
                Surf.Name = nm;
                Surf.Area = 0;
                Surf.Type = 'Mask';
                Surfa{i,1} = Surf;
            else
                errordlg('This file is not a known Surface File. Please try again');
                return
            end
    end
end
%=========================End of main program=============================%
return;

