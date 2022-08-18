%clear;  
function [Y]=occupx(TF,Eseq,xfit)
path1='/nuctf_equi_bai/ndr_call/';
path2='/nuctf_equi_bai/tf_energy_all/Etf_allmat_chr/';
fpath= strcat(path2,'Emtfall.mat');
load(fpath); % mean TF energy
load /nuctf_equi_bai/nuc_energy/E_Em.mat; % Nuc energy
tfn=30;
TFlist=importdata('/nuctf_equi_bai/tf_energy_all/listbai_all.txt'); % list of 104 TF names and sizes
TF_list1=TFlist.data; clear TFlist;
load /nuctf_equi_bai/tf_energy_all/tfindx.txt; % tfindx is a set (size 104) of ranked (or rearranged) TF indices of listbai_all.txt. You can change to any set of TFs
TF_list=TF_list1(tfindx(1:tfn),1); % tfn=5, means considering only top 5 TFs.
l=147; dLb=2000; L=5*dLb;    
Ei=-7; CHR=16; Etfmul=cell(1,CHR); % CHR no. of chromosomes
for chr=1:CHR
    fname = sprintf('Etf_chr%d.mat',chr); % TF energy
    fpath= strcat(path2,fname);
    load(fpath);
    Etfmul{1,chr}=Etf(:,tfindx(1:tfn));
end                      
   a=Emtf(tfindx(1:tfn),1)';
   E2=a;
   for i=1:(dLb-1)
       E2=cat(1,E2,a);
   end
   E1=zeros(dLb,1)+Em; clear a;

   CN=xfit(1,1); gamaN=xfit(2,1);   
   C=xfit(3:(3+tfn-1),1); gama=xfit((3+tfn):((tfn+1)*2),1);   
   % Eseq=1 if sequence effect present, otherwise Eseq=0; TF=1 if TF present, otherwise TF=0;  
   if TF<1
      C=zeros(length(xfit(3:(3+tfn-1),1)),1);
   end
   if Eseq<1
      for chr=1:CHR
          E{chr}=-Ei*ones(length(E{chr}),1); CN=0.001; gamaN=1; 
      end      
      Em=-Ei; E1=zeros(dLb,1)+Em; 
   end
   nor=1e140; Y=cell(CHR,1); Pngtf=cell(tfn,1); pLmcnt=0; tic; 
for chr=1:CHR
    Lgnome=length(E{1,chr});
    m=mod(Lgnome-3*dLb,dLb);
    if m==0
       pLm=1+(Lgnome-3*dLb)/dLb; EN=E{1,chr}; ETF=Etfmul{1,chr}; 
    else
       pLm=2+(Lgnome-3*dLb-m)/dLb; 
       EN=cat(1,E{1,chr},(zeros(dLb-m,1)+Em));       
       ETF=cat(1,Etfmul{1,chr},E2(1:dLb-m,1:tfn));     
    end   
    i1=0; pLmcnt=pLmcnt+pLm;
    while i1<=(pLm-1)
        Zf=zeros(L,1); Zb=zeros(L,1); ks=min(TF_list);
        Es=cat(1,E1,EN((1+dLb*i1):(dLb*i1+3*dLb),1),E1);  
        Zf(1:ks-1,1)=ones(ks-1,1)/nor;
        Zb(L-ks+2:L,1)=ones(ks-1,1)/nor;
        Estf=cat(1,E2,ETF((1+dLb*i1):(dLb*i1+3*dLb),1:tfn),E2);
        for i=1:(L-ks+1) % DNA length
            Zbsum=0; Zfsum=0;
            jf=i+(ks-1)-TF_list(1:tfn)+1; % position of TF, nuc
            jbx=L-ks+2-i+TF_list(1:tfn);       
            for t=1:tfn %T #transcription factor                           
                if jf(t)>=1                      % forward TF  
                   if jf(t)==1
                     Zfor=C(t,1)*exp(gama(t,1)*Estf(jf(t),t))/nor; 
                   else
                     Zfor=Zf(jf(t)-1)*C(t,1)*exp(gama(t,1)*Estf(jf(t),t));  
                   end               
                else
                   Zfor=0;
                end
                jb=L-ks+2-i;                    % backward TF
                if jbx(t)<=(L+1)
                   if jbx(t)==(L+1)
                      Zbor=C(t,1)*exp(gama(t,1)*Estf(jb,t))/nor;
                   else
                      Zbor=Zb(jbx(t))*C(t,1)*exp(gama(t,1)*Estf(jb,t));
                   end
                else
                   Zbor=0;
                end
                Zfsum=Zfsum+Zfor; 
                Zbsum=Zbsum+Zbor;
            end
            if (i+(ks-1))>=l                  % forwards nucleosome
               if (i+(ks-1))==l
                  Zfor=CN*exp(gamaN*Es((i+(ks-1)-l+1),1))/nor;
               else
                 Zfor=Zf(i+(ks-1)-l)*CN*exp(gamaN*Es((i+(ks-1)-l+1),1));
               end
            else
               Zfor=0;
            end
            if (ks-1+i)>=l                    % backwards nucleosome
               if (L-ks+2-i+l)>L
                  Zbor=CN*exp(gamaN*Es((L-ks+2-i),1))/nor;
               else
                  Zbor=Zb(L-ks+2-i+l)*CN*exp(gamaN*Es((L-ks+2-i),1));
               end
            else
               Zbor=0;
            end
            Zfsum=Zfsum+Zfor;         
            Zbsum=Zbsum+Zbor;  
            Zf(i+(ks-1))=Zf(i+(ks-1)-1) + Zfsum; % forward sums
            Zb(L-ks+2-i)=Zb(L-ks+2-i+1) + Zbsum; % backward sums
        end
       clear jf jbx; 
       Pn=zeros(L,1);
       Pn(1,1)=(1/nor)*CN*exp(gamaN*Es(1,1))*Zb(1+l)/Zf(L);
       Pn(L-l+1,1)=Zf(L-l)*CN*exp(gamaN*Es(L-l+1,1))*(1/nor)/Zf(L);
       Pn(2:L-l,1)=(((Zf((2:L-l)-1)*CN).*exp(gamaN*Es(2:L-l,1))).*Zb((2:L-l)+l))/Zf(L);
       if i1==0
          Png=Pn(dLb+1:3*dLb,1);
       elseif i1==(pLm-1)
          Pni=Pn(2*dLb+1:4*dLb,1);
          Png=cat(1,Png,Pni);
       else
          Pni=Pn(2*dLb+1:3*dLb,1);
          Png=cat(1,Png,Pni);
       end
       Pntf=cell(tfn,1); Pnitf=cell(tfn,1);
       for t=1:tfn
           Pntf{t,1}=zeros(L,1); ltf=TF_list(t,1);
           Pntf{t,1}(1,1)=(1/nor)*C(t,1)*exp(gama(t,1)*Estf(1,t))*Zb(1+ltf)/Zf(L);
           Pntf{t,1}(L-ltf+1,1)=Zf(L-ltf)*C(t,1)*exp(gama(t,1)*Estf(L-ltf+1,t))*(1/nor)/Zf(L);
           Pntf{t,1}(2:L-ltf,1)=(((Zf((2:L-ltf)-1)*C(t,1)).*exp(gama(t,1)*Estf(2:L-ltf,t))).*Zb((2:L-ltf)+ltf))/Zf(L);
           if i1==0
              Pngtf{t,1}=Pntf{t,1}(dLb+1:3*dLb,1);
           elseif i1==(pLm-1)
              Pnitf{t,1}=Pntf{t,1}(2*dLb+1:4*dLb,1);
              Pngtf{t,1}=cat(1,Pngtf{t,1},Pnitf{t,1});
           else
              Pnitf{t,1}=Pntf{t,1}(2*dLb+1:3*dLb,1);
              Pngtf{t,1}=cat(1,Pngtf{t,1},Pnitf{t,1});
           end
       end
       if Zf(L)==Inf
          nor=nor*1e10; i1=-1; jc=1; fprintf('inf...\n');
       elseif (Zf(L)/nor)<1e-300
          nor=nor/1e10; i1=-1; jc=1; fprintf('zero...\n');
       end       
       i1=i1+1; clear Pn Pni Es Estf; 
    end
    clear EN ETF; O=zeros(Lgnome,1); ks=l; %ks=l nucleosome
    for i=1:Lgnome 
        if i<ks
           O(i)=sum(Png(1:i));
        elseif i>(Lgnome-ks+1)   
           O(i)=sum(Png(i-ks+1:Lgnome-ks+1));
        else
           O(i)=sum(Png((i-ks+1):i));
        end
    end
    ON=nor*O; clear Png; 
    Y{chr,1}=ON; Ytfx=zeros(Lgnome,0);
    for t=1:tfn
        O=zeros(Lgnome,1); ks=TF_list(t,1);
        for i=1:Lgnome 
            if i<ks
               O(i)=sum(Pngtf{t,1}(1:i));
            elseif i>(Lgnome-ks+1)   
               O(i)=sum(Pngtf{t,1}(i-ks+1:Lgnome-ks+1));
            else
               O(i)=sum(Pngtf{t,1}((i-ks+1):i));
            end
        end
        ON=nor*O; clear O; 
        Ytfx=cat(2,Ytfx,ON);
    end
    Y{chr,1}=cat(2,Y{chr,1},Ytfx);
    fprintf('chr %d \n',chr); 
end 
toc;
fpath= strcat(path1,'yy3A_lee.mat');
load(fpath);
ONLx=cell(CHR,1); occup_nuc=cell(CHR,1);
for chr=1:CHR
    occup_nuc{chr,1}=Y{chr,1}(:,1); ONLx{chr,1}=occup_nuc{chr,1}(x1_lee{chr,1}(100:end-100));
    if chr==1
       Om=occup_nuc{1,1}; ONL=ONLx{chr,1}; yL=y1_lee{chr,1}(100:end-100);
    else
       Om=cat(1,Om,occup_nuc{chr,1}); ONL=cat(1,ONL,ONLx{chr,1}); yL=cat(1,yL,y1_lee{chr,1}(100:end-100));
    end
end
avg_occ=mean(Om);
Oer=ONL-yL;
sigma=sqrt((Oer'*Oer)/length(Oer));
fpath= strcat(path1,'ndrpos_chrA.mat');
load(fpath);
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
end
