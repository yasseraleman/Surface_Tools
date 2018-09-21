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