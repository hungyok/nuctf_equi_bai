%clear;
function [NDR_X]=ndr_cut(i)
    load /nuctf_equi_bai/ndr_call/yy1_lee.mat;
    %A=1; A1=2; B=3; C=4; D=5; E=6; F=7; G=8;
    NDR_X=cell(16,1); 
<<<<<<< HEAD
    SD=0.0678; NDRcut=0.8-(i-1)*SD;   
=======
    sd=0.0678; NDRcut=0.8-(i-1)*sd;   
>>>>>>> ccacfc09b980e21049c52b9ef6d21469c3985309
    for chr=1:16
        mposx=zeros(1,3); 
        x1=x1_lee{chr}; yy1=y1_lee{chr}; 
        flgin=1; cnt=0;
        for j=1:length(x1)
            if yy1(j)<NDRcut
               if flgin>0
                  jin=j; flgin=0;
               end
               cnt=cnt+1;
            else
               jout=j; 
               if cnt>0
                  dx=x1(jout,1)-1-(x1(jin,1)-1); 
                  if dx>4
                     mpos(1,1)=x1(jin,1)-1; mpos(1,2)=x1(jout,1)-1; mpos(1,3)=NDRcut;
                     mposx=cat(1,mposx,mpos);          
                  end
               end
               cnt=0; flgin=1;
            end        
        end 
        NDR_X{chr,1}=mposx(2:end,:);
        fprintf('i...%d...chr...%d \n',i,chr); 
    end 
    %fname=sprintf('NDR_%d.mat',i);
    %fpath=strcat(pathx,fname);
    %save(fpath,'NDR_X','-v7.3');
end


