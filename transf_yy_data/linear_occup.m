%clear;
function [x1_lee,y1_lee]=linear_occup(S2,c,g,d_lee,d_indx)
 j=1; 
 for i=1:length(S2)  % Lee
     x1(j,1)=S2(i,1); 
     y1(j,1)=S2(i,2);
     if y1(j,1)>1.26 
        j=j-1;
     end
     j=j+1;
 end
 y1s=smooth(y1,5);
 y1=y1s;
 y2=c*exp(g*y1); %y2mean=mean(y2); 
 j00=1; x1_lee=cell(16,1); y1_lee=cell(16,1);
  for chr=1:16
      J=d_indx(chr,1); x1=d_lee(j00:J,1); yy1=y2(j00:J); 
      x1_lee{chr,1}=x1;
      y1_lee{chr,1}=yy1;
      j00=J+1;
  end
end
 %save yourpath\transf_yy_data\para_cg.txt cg -ascii;
    
