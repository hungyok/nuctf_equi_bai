 %clear;  
 function [sigma]=optimbfunc_ver3a(xi,E,Em,Etfmul,Emtf)
 load D:\Cell_protocol\NucTF\simplexM_tf30\input\yy3A_lee.mat;
 load D:\Cell_protocol\NucTF\simplexM_tf30\input\rand_genomeB.mat;
 tfn=30; 
 TFlist=importdata('D:\Cell_protocol\tf_energy_all\listbai_all.txt');
 TF_list1=TFlist.data; clear TFlist;
 load D:\Cell_protocol\tf_energy_all\tfindx.txt;
 TF_list=TF_list1(tfindx(1:tfn),1); 
 l=147; dLb=3500; dLx=750; L=3*dLb+2*dLx; 
 invld=0; 
 if xi(1,1)>10 || xi(1,1)<0
    invld=1;
 end
 if xi(2,1)>5 || xi(2,1)<0
    invld=1;
 end
 for i7=3:(3+tfn-1)
     if xi(i7,1)>0.5 || xi(i7,1)<0
        invld=1;
     end
 end
 for i7=(3+tfn):((tfn+1)*2)
     if xi(i7,1)>5 || xi(i7,1)<0
        invld=1;
     end
 end
 if invld==1
    sigma=5;
 else
    a=Emtf(tfindx(1:tfn),1)';
    E2=a;
    for i=1:(dLb-1)
        E2=cat(1,E2,a);
    end
    E1=zeros(dLx,1)+Em;
    %C=0.0*ones(tfn,1); gama=ones(tfn,1);
    C=xi(3:(3+tfn-1),1); gama=xi((3+tfn):((tfn+1)*2),1); CN=xi(1,1); gamaN=xi(2,1);
    ONLx=cell(16,1); pLmsum=0; nor=1e120; tic;
    for chr=1:chrindx(1,1)
        Lgnome=length(E{1,chrlist{1,1}(chr)});
        m=mod(Lgnome-3*dLb,dLb);
        if m==0
           pLm=1+(Lgnome-3*dLb)/dLb; EN=E{1,chrlist{1,1}(chr)}; ETF=Etfmul{1,chrlist{1,1}(chr)};
        else
           pLm=2+(Lgnome-3*dLb-m)/dLb;
           EN=cat(1,E{1,chrlist{1,1}(chr)},(zeros(dLb-m,1)+Em));
           ETF=cat(1,Etfmul{1,chrlist{1,1}(chr)},E2(1:dLb-m,1:tfn));
        end
        i1=0; pLmsum=pLmsum+pLm;
        while i1<=(pLm-1)
              Zf=zeros(L,1); Zb=zeros(L,1); ks=min(TF_list);
              Zf(1:ks-1,1)=ones(ks-1,1)/nor;
              Zb(L-ks+2:L,1)=ones(ks-1,1)/nor;
              Es=cat(1,E1,EN((1+dLb*i1):(dLb*i1+3*dLb),1),E1);
              Estf=cat(1,E2(1:dLx,1:tfn),ETF((1+dLb*i1):(dLb*i1+3*dLb),1:tfn),E2(1:dLx,1:tfn));
              for i=1:(L-ks+1) % DNA length
                  Zbsum=0; Zfsum=0;
                  jf=i+(ks-1)-TF_list(1:tfn)+1; % position of TF, nuc
                  jbx=L-ks+2-i+TF_list(1:tfn); jb=L-ks+2-i;      
                  for t=1:tfn %T #transcription factor                           
                      if jf(t)>=1                      % forward TF  
                         if jf(t)==1
                            Zfor=C(t,1)*exp(gama(t,1)*Estf(jf(t),t))/nor; 
                         else
                            Zfor=Zf(jf(t)-1)*C(t,1)*exp(gama(t,1)*Estf(jf(t),t));  
                         end   %fprintf('bbbb C(t,1) %f exp(gama(t,1)*Estf(jf,1)) %f \n',C(t,1),exp(gama(t,1)*Estf(jf,1)));             
                      else
                         Zfor=0;
                      end                                  
                      if jbx(t)<=(L+1) % backward TF
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
                  %fprintf('i1 %d Zfor %f Zf %f \n',i1,Zfor,Zf(i+(ks-1)));
              end
              Pn=zeros(L,1);
              Pn(1,1)=(1/nor)*CN*exp(gamaN*Es(1,1))*Zb(1+l)/Zf(L);
              Pn(L-l+1,1)=Zf(L-l)*CN*exp(gamaN*Es(L-l+1,1))*(1/nor)/Zf(L);
              Pn(2:L-l,1)=(((Zf((2:L-l)-1)*CN).*exp(gamaN*Es(2:L-l,1))).*Zb((2:L-l)+l))/Zf(L);
              if i1==0
                 Png=Pn(dLx+1:(2*dLb+dLx),1);
              elseif i1==(pLm-1)
                 Pni=Pn(dLb+dLx+1:(3*dLb+dLx),1);
                 Png=cat(1,Png,Pni);
              else
                 Pni=Pn(dLb+dLx+1:(2*dLb+dLx),1);
                 Png=cat(1,Png,Pni);
              end
              if Zf(L)==Inf
                 nor=nor*1e10; i1=-1; fprintf('inf...\n');
              elseif (Zf(L)/nor)<1e-300  
                 nor=nor/1e10; i1=-1; fprintf('zero...\n');
              end       
              i1=i1+1; %fprintf('pLm %d i1 %d Zf(L) %f \n',pLm,i1,Zf(L));
        end
        s=size(Png(1:Lgnome,1)); Lg=s(1,1); O=zeros(Lg,1); ks=l;
        for i=1:Lg %ks=l nucleosome
            if i<ks
               O(i)=sum(Png(1:i));
            elseif i>(Lg-ks+1)
               O(i)=sum(Png(i-ks+1:Lg-ks+1));
            else
               O(i)=sum(Png((i-ks+1):i));
            end
        end
        ON=nor*O;
        ONLx{chrlist{1,1}(chr),1}=ON(x1_lee{chrlist{1,1}(chr),1}(100:end-100));
        %fprintf('chr %d pLm %d \n',chr,pLm); 
    end
    fprintf('function optimbfunc_ver3a...tfn...%d...pLmsum=%d\n',tfn,pLmsum);
    toc;
    for chr=1:chrindx(1,1)
        if chr==1
           ONL=ONLx{chrlist{1,1}(chr),1}; yL=y1_lee{chrlist{1,1}(chr),1}(100:end-100);
        else
           ONL=cat(1,ONL,ONLx{chrlist{1,1}(chr),1}); yL=cat(1,yL,y1_lee{chrlist{1,1}(chr),1}(100:end-100)); 
        end
    end
    Oer=ONL-yL;
    sigma=sqrt((Oer'*Oer)/length(Oer)); 
 end
end
