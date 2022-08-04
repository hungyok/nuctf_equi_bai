% occupancy, remodeling, avg_occ, sigma, Pndr, Auc
function [Y]=occupR(remod,xfit)
load /nuctf_equi_bai/nuc_energy/E_Em.mat;
tfnx=floor(length(xfit)/2)-1;
xfit(1)=remod(1); xfit(2)=remod(2);

TFlist=importdata('/nuctf_equi_bai/tf_energy_all/listbai_all.txt'); % list of TF name and size 
TF_list1=TFlist.data; clear TFlist;
load /nuctf_equi_bai/tf_energy_all/tfindx.txt; % ranked listbai_all.txt indices
TF_list=TF_list1(tfindx(1:tfnx),1);
load /nuctf_equi_bai/NucRemod/simplexM_remod/input/pos_octf/pos_octf.mat;
pos_octf1=pos_octf;
pos_octf=cell(1,3);
nx=100; ht=remod(3)/0.212; sig=remod(4); A=ht; B=1/(2*sig^2); CHR=16; occup_nuc=cell(CHR,1);
for ix=1:nx % # realization
    tic;
    for chr=1:CHR % indie_test 2
        E1=E{chr};
        pos_octf{1,1}=pos_octf1{chr,1};
        pos_octf{1,2}=pos_octf1{chr,2};
        pos_octf{1,3}=pos_octf1{chr,3};
        [xptfa]=tf_cluster(tfnx,pos_octf); xptf=xptfa;
        lx=length(E{chr}); xx=1:lx; x=xx';
        Esum=zeros(lx,1);
        for jx=1:length(xptf(:,1))
            x1=pos_octf{1,3}(xptf(jx,1),2:end); x2=TF_list(x1(x1>0));
            xo=pos_octf{1,1}(xptf(jx,1),1)+floor(x2(xptf(jx,2))/2)+1;
            x1=x-xo; xo1=pos_octf{1,1}(xptf(jx,1),1)-(147+floor(x2(xptf(jx,2))/2)); x2=x-xo1;
            if xo1>0
               yr=A*exp(-(x1.^2)*B); 
               yr(1:(xo-1),1)=zeros((xo-1),1);
               yl=A*exp(-(x2.^2)*B); 
               yl(xo1+1:lx,1)=zeros((lx-xo1),1);
               ym=zeros(lx,1); ym((xo1+1):(xo-1))=A;
               y=yl+yr+ym;
            else
               yr=A*exp(-(x1.^2)*B);
               yr(1:(xo-1),1)=zeros((xo-1),1);
               yl=zeros(lx,1);
               ym=zeros(lx,1); ym(1:(xo-1))=A;
               y=yl+yr+ym;
            end
            Esum=Esum+y;
        end
        E1=E1-Esum;
        occup_nuc{chr,1}=occup_nucs(xfit,E1,Em,chr);
    end
    if ix<2
       occup_nuc0=occup_nuc;  
    else
       for chr=1:CHR
           occup_nuc0{chr,1}=occup_nuc0{chr,1}+occup_nuc{chr,1};
       end
    end
    fprintf('ix...%d \n',ix); toc;
end
load /nuctf_equi_bai/ndr_call/yy3A_lee.mat;
ONLx=cell(CHR,1);
for chr=1:CHR
    occup_nuc{chr,1}=occup_nuc0{chr,1}*(1/nx); ONLx{chr,1}=occup_nuc{chr,1}(x1_lee{chr,1}(100:end-100));
    if chr==1
       Om=occup_nuc{1,1}; ONL=ONLx{chr,1}; yL=y1_lee{chr,1}(100:end-100);
    else
       Om=cat(1,Om,occup_nuc{chr,1}); ONL=cat(1,ONL,ONLx{chr,1}); yL=cat(1,yL,y1_lee{chr,1}(100:end-100));
    end
end
avg_occ=mean(Om);
Oer=ONL-yL;
sigma=sqrt((Oer'*Oer)/length(Oer));
load /nuctf_equi_bai/ndr_call/ndrpos_chrA.mat;
ndrcnt=0; NDRcut=0.6643; sumNDRx=0;
for chr=1:CHR
    for j=1:length(ndr_chr{chr,1})
        ndrp=floor((ndr_chr{chr,1}(j,1)+ndr_chr{chr,1}(j,2))/2);
        ndr=ndrp+39; i=ndrp-39; cnt=0;
        while i<=ndr
              if occup_nuc{chr,1}(i,1)<=NDRcut
                 cnt=cnt+1;
              end
              i=i+1;
        end
        if cnt>0
           ndrcnt=ndrcnt+1;
        end
    end
    fprintf('NDR chr %d \n',chr); sumNDRx=sumNDRx+length(ndr_chr{chr,1});
end
Pndr=ndrcnt/sumNDRx;
%Roc_auc calculation
cut_off=0.04:0.01:0.96;
TPR=zeros(length(cut_off),1); FPR=zeros(length(cut_off),1);
for cut=1:length(cut_off)
    NDRcut=cut_off(cut); TP=0; FN=0; FP=0; TN=0;
    for chr=1:CHR
        mposx=zeros(1,3);
        x1=x1_lee{chr,1}; yy1=y1_lee{chr,1};
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
                     mpos(1,1)=x1(jin,1)-1; mpos(1,2)=x1(jout,1)-1; mpos(1,3)=0;
                     mposx=cat(1,mposx,mpos);
                  end
               end
               cnt=0; flgin=1;
            end
        end
        if length(mposx(:,1))>2
           lx=length(mposx(2:end,1));
           bipos1a=zeros(lx+lx-1,3);
           bipos1a(1,:)=mposx(2,:); j=2;
           for i=1:(lx-1)
               mpos(1,1)=mposx(i+1,2); mpos(1,2)=mposx(i+2,1); mpos(1,3)=1;
               bipos1a(j,:)=mpos;
               bipos1a(j+1,:)=mposx(i+2,:); j=j+2;
           end
           for i=1:length(bipos1a)
               xx=floor((bipos1a(i,2)+bipos1a(i,1))/2);
               biposb=occup_nuc{chr,1}(xx,1);
               if biposb<NDRcut
                  nuc=0;
               else
                  nuc=1;
               end
               if nuc==0 && bipos1a(i,3)==0
                  TP=TP+1;
               elseif nuc==1 && bipos1a(i,3)==0
                  FN=FN+1;
               elseif nuc==0 && bipos1a(i,3)==1
                  FP=FP+1;
               else
                  TN=TN+1;
               end
           end
        end
    end
    fprintf('cut....%d \n',cut);
    if (FP+TN)==0
       FPR(cut,1)=0;
    else
       FPR(cut,1)=FP/(FP+TN);
    end
    if (TP+FN)==0
       TPR(cut,1)=0;
    else
       TPR(cut,1)=TP/(TP+FN);
    end
end
Auc=0;
for i=2:length(cut_off)
    Auc=((TPR(i,1)+TPR(i-1,1))/2)*abs(FPR(i,1)-FPR(i-1,1))+Auc;
end
fprintf('ndrcnt...%d...sumNDRx...%d...avg_occ...%.6f...sigma...%.6f...Pndr...%.6f...Auc...%.6f \n',ndrcnt,sumNDRx,avg_occ,sigma,Pndr,Auc); 
% important quality parameters: mean occupancy(avg_occ), sigma(RMSD), Pndr(NDR count), Auc(ROC area-under-curve).
Y=occup_nuc;
%save /nuctf_equi_bai/NucRemod/occup_profile/Y.mat Y -v7.3;
end
