%clear;
function [c,g]=transf_yy_data(S2)
j=1; cg=zeros(1,2); para=zeros(0,2);
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
y2=zeros(length(y1),1);

 c=0.82;
 while c<0.852
       g=0.11;
       while g<0.18
             cnt=0; i=1;
             while i<=length(y1)
                   y2(i,1)=c*exp(g*y1(i,1));
                   if y2(i,1)<1
                      cnt=cnt+1;
                   else
                      i=length(y1);
                   end
                   i=i+1; 
             end
             if cnt==length(y1)
                if abs(mean(y2)-0.8)<=0.0001 
                   g1=g; c1=c; cg(1)=c; cg(2)=g;
                   para=cat(1,para,cg);
                   g=2; fprintf('g1 %f c1 %f cnt %d mean(y2) %f \n',g1,c1,cnt,mean(y2));
                end
             end
             g=g+0.0001; 
       end 
       c=c+0.0001; fprintf('c %f \n',c);
 end
 [x,I]=max(para(:,2));
 c=para(I,1); g=para(I,2);
end
 %save /nuctf_equi_bai/transf_yy_data/para_cg.txt cg -ascii;
    
