function [xptfa]=tf_cluster(tfn,pos_octf)
TFlist=importdata('D:\Cell_protocol\tf_energy_all\listbai_all.txt'); % list of TF name and size 
TF_list1=TFlist.data; clear TFlist;
load D:\Cell_protocol\tf_energy_all\tfindx.txt; % ranked listbai_all.txt indices
TF_list=TF_list1(tfindx(1:tfn),1);
CL1=cell(1,1); flgnn=0; rng('shuffle');
for i=1:length(pos_octf{1,1})
    cls=i; j=i+1; flg=1; 
    x1=pos_octf{1,3}(i,2:end); x=TF_list(x1(x1>0)); xa=pos_octf{1,1}(i)+floor(max(x)/2);
    while flg>0 && j<(length(pos_octf{1,1})+1)
          if i<length(pos_octf{1,1})
             x1=pos_octf{1,3}(j,2:end); x=TF_list(x1(x1>0)); xb=pos_octf{1,1}(j)-floor(max(x)/2);
          else
             xb=pos_octf{1,1}(i)+500;
          end
          dx=xb-xa;
          if dx<1
             cls=cat(1,cls,j);
          else
             flg=0;
          end
          j=j+1;
    end
    if length(cls)>1
       CL1=cat(1,CL1,cls); flgnn=1; %cluster
    else
       if flgnn==0
          cls=i; CL1=cat(1,CL1,cls); %isolation
       else
          flgnn=0;
       end
    end
end
CL=CL1(2:end,1); xptf1=zeros(1,2); xp=zeros(1,2); flchge=0; 
xptfa1=zeros(1,2); % remodeling position
for i=1:length(CL)
    n1=length(CL{i}); sump=0; CLp1=zeros(1,2); oc1=zeros(1,1); flgnn1=1;
    for j=1:n1
        x1=pos_octf{1,2}(CL{i}(j),2:end); x=x1(x1>0); oc1=cat(2,oc1,x); 
        sump=sump+pos_octf{1,2}(CL{i}(j),1);
        for j1=1:pos_octf{1,2}(CL{i}(j),1) 
            cx=[CL{i}(j),j1];
            CLp1=cat(1,CLp1,cx);
        end
    end
    CLp=CLp1(2:end,:); pmx=max(oc1(2:end)); nn1=zeros(sump,2);
    for j=1:sump
        nn1(j,1)=0; nn1(j,2)=j;
    end
    if sump>1
       nn2=nchoosek(1:sump,2);
       nn=cat(1,nn1,nn2); 
    else
       nn=nn1; 
    end
    if rand<=pmx
       if flchge>0 % compare & filter
          nn1=nn; 
          %fprintf('1..chr...%d...i...%d...n1...%d...sump...%d...length nn...%d \n',chr,i,n1,sump,length(nn(:,1)));
          for j=1:2
              if xp(1,j)>0
                 k=1;
                 while k<(length(nn(:,1))+1)
                       if nn(k,1)>0
                          if xp(1,j)==CLp(nn(k,1),1) || xp(1,j)==CLp(nn(k,2),1)
                             nn1(k,1)=0; nn1(k,2)=0;
                          end
                       else
                          if xp(1,j)==CLp(nn(k,2),1)
                             nn1(k,1)=0; nn1(k,2)=0;
                          end
                       end
                       k=k+1;
                 end
              end
          end
          if sum(nn1(:,2))>0
             nn=nn1(nn1(:,2)>0,:); flgnn1=1; 
             %fprintf('2..chr...%d...i...%d...n1...%d...sump...%d...length nn...%d \n',chr,i,n1,sump,length(nn(:,1))); 
          else
             flgnn1=0;
          end
       end
       %fprintf('3..chr...%d...i...%d...n1...%d...sump...%d...length nn...%d \n',chr,i,n1,sump,length(nn(:,1)));
       if flgnn1>0
          j=randi(length(nn(:,1))); flchge=1; 
          if nn(j,1)==0
             xp(1,1)=0; xp(1,2)=CLp(j,1); %single 
             xptf1=cat(1,xptf1,[xp(1,2),CLp(j,2)]);
             x1=pos_octf{1,2}(xp(1,2),2:end); x=x1(x1>0); pt=x(CLp(j,2));
             if pmx*rand<=pt
                xptfa1=cat(1,xptfa1,[xp(1,2),CLp(j,2)]);
             end
          else
             xp(1,1)=CLp(nn(j,1),1); xp(1,2)=CLp(nn(j,2),1); %double 
             xptf1=cat(1,xptf1,[xp(1,1),CLp(nn(j,1),2)]); 
             x1=pos_octf{1,2}(xp(1,1),2:end); x=x1(x1>0); pt=x(CLp(nn(j,1),2));
             if pmx*rand<=pt
                xptfa1=cat(1,xptfa1,[xp(1,1),CLp(nn(j,1),2)]);
             end
             xptf1=cat(1,xptf1,[xp(1,2),CLp(nn(j,2),2)]);
             x1=pos_octf{1,2}(xp(1,2),2:end); x=x1(x1>0); pt=x(CLp(nn(j,2),2));
             if pmx*rand<=pt
                xptfa1=cat(1,xptfa1,[xp(1,2),CLp(nn(j,2),2)]);
             end
         end
       end
    else
       flchge=0;
    end
    %fprintf('CL i...%d\n',i);
end
xptfa=xptfa1(2:end,:); 
end
