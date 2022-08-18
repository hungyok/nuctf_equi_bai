% occupancy
function [occupx]=occup_nucs(xfit,E1,Em,chr1) 
path1 = '/nuctf_equi_bai/tf_energy_all/Etf_allmat_chr/';
fpath= strcat(path1,'Emtfall.mat');
load(fpath); % mean TF energy
TFlist=importdata('/nuctf_equi_bai/tf_energy_all/listbai_all.txt'); % list of TF name and size 
TF_list1=TFlist.data; clear TFlist;
load /nuctf_equi_bai/tf_energy_all/tfindx.txt; % ranked listbai_all.txt indices: 1st or top (high NDR predictor) to 104th or bottom (low NDR predictor) TF
l=147; dLb=2000; L=5*dLb;    
tfn=floor(length(xfit)/2)-1; E=E1;
fname = sprintf('Etf_chr%d.mat',chr1);
fpath= strcat(path1,fname);
load(fpath);
Etfmul=Etf(:,tfindx(1:tfn));                   
TF_list=TF_list1(tfindx(1:tfn),1);
a=Emtf(tfindx(1:tfn),1)';
E2=a;
for i=1:(dLb-1)
    E2=cat(1,E2,a);
end
E1a=zeros(dLb,1)+Em; clear a;
CN=xfit(1,1); gamaN=xfit(2,1);   
C=xfit(3:(3+tfn-1),1); gama=xfit((3+tfn):((tfn+1)*2),1);
nor=1e150; 
Lgnome=length(E);
m=mod(Lgnome-3*dLb,dLb);
if m==0
   pLm=1+(Lgnome-3*dLb)/dLb; EN=E; ETF=Etfmul; 
else
   pLm=2+(Lgnome-3*dLb-m)/dLb; 
   EN=cat(1,E,(zeros(dLb-m,1)+Em));       
   ETF=cat(1,Etfmul,E2(1:dLb-m,1:tfn));     
end   
i1=0;
while i1<=(pLm-1)
      Zf=zeros(L,1); Zb=zeros(L,1); ks=min(TF_list);
      Es=cat(1,E1a,EN((1+dLb*i1):(dLb*i1+3*dLb),1),E1a);  
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
      if Zf(L)==Inf
         nor=nor*1e10; i1=-1; fprintf('inf...\n');
      elseif (Zf(L)/nor)<1e-300
         nor=nor/1e10; i1=-1; fprintf('zero...\n');
      end       
      i1=i1+1; %clear Pn Pni Es Estf; 
end
O=zeros(Lgnome,1); ks=l; % clear EN ETF;
for i=1:Lgnome %ks=l nucleosome
    if i<ks
       O(i)=sum(Png(1:i));
    elseif i>(Lgnome-ks+1)   
       O(i)=sum(Png(i-ks+1:Lgnome-ks+1));
    else
       O(i)=sum(Png((i-ks+1):i));
    end
end
ON=nor*O;  
occupx=ON; %fprintf('chr1...%d\n',chr1);
end
