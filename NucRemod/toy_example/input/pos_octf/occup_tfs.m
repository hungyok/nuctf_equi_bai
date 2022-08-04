 % To find TF occupancy and position at occupancy tfcut=0.0022 
 load /nuctf_equi_bai/nuc_energy/E_Em.mat;
 load /nuctf_equi_bai/tf_energy_all/Etf_allmat_chr/Emtfall.mat;
 load /nuctf_equi_bai/tf_energy_all/tfindx.txt;
 TFlist=importdata('/nuctf_equi_bai/tf_energy_all/listbai_all.txt');
 TF_list1=TFlist.data; clear TFlist;
 l=147; dLb=2000; L=5*dLb; 
 pth = '/nuctf_equi_bai/tf_energy_all/Etf_allmat_chr/';
 load /nuctf_equi_bai/NucTF/toy_example/output/simplex_xval.txt;
 xfit=simplex_xval(:,1); tfn=floor(length(xfit)/2)-1; tfcut=zeros(tfn,1); tfcut(:,1)=0.0022;
 Etfmul=cell(1,16); Pngtf=cell(tfn,1); pos_octf=cell(16,3); 
for k =1:16
    filename = sprintf('Etf_chr%d.mat',k);
    fpath= strcat(pth,filename);
    load(fpath);
    Etfmul{1,k}=Etf(:,tfindx(1:tfn));
end                   
   TF_list=TF_list1(tfindx(1:tfn),1);
   a=Emtf(tfindx(1:tfn),1)';
   E2=a;
   for i=1:(dLb-1)
       E2=cat(1,E2,a);
   end
   E1=zeros(dLb,1)+Em; clear a;
   CN=xfit(1,1); gamaN=xfit(2,1);   %Ea=0.9392; 
   C=xfit(3:(3+tfn-1),1); gama=xfit((3+tfn):((tfn+1)*2),1);
   nor=1e140; pLmcnt=0; tic;
for chr=1:16
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
           % fprintf('i1 %d Zfor %f Zf %f \n',i1,Zfor,Zf(i+(ks-1)));
        end
        clear jf jbx; Pntf=cell(tfn,1); Pnitf=cell(tfn,1);
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
        i1=i1+1; clear Es Estf; %fprintf('chr...%d...pLm...%d...i1...%d \n',chr,pLm,i1);
    end
    clear EN ETF; %ON1=zeros(Lgnome,1);  
    Of=zeros(Lgnome,1); pos_tf1=zeros(Lgnome,1); pos_oc1=zeros(Lgnome,1);
    for t=1:tfn
        Of1=zeros(Lgnome,1);
        s=size(Pngtf{t,1}(1:Lgnome,1)); Lg=s(1,1); O=zeros(Lg,1); ks=TF_list(t,1);
        for i=1:Lg %ks=l nucleosome
            if i<ks
               O(i)=sum(Pngtf{t,1}(1:i));
            elseif i>(Lg-ks+1)   
               O(i)=sum(Pngtf{t,1}(i-ks+1:Lg-ks+1));
            else
               O(i)=sum(Pngtf{t,1}((i-ks+1):i));
            end
        end
        ON=nor*O; clear O; %ONL=cat(1,ONL,ON); 
        flgin=1; cnt=0; mposx=zeros(1,2);
        for i=1:length(ON)
            if ON(i)>tfcut(t)
               if flgin>0
                  jin=i; flgin=0;
               end
               cnt=cnt+1;
            else
               jout=i;
               if cnt>1
                  mpos(1,1)=floor((jin+jout)/2); mpos(1,2)=max(ON(jin:jout));
                  mposx=cat(1,mposx,mpos);
               end
               cnt=0; flgin=1;
            end 
        end
        Ot=mposx(2:end,:);
        if length(Ot)>0
           Of1(Ot(:,1))=Ot(:,2);
           pos_oc1=cat(2,pos_oc1,Of1);
           Of1(Ot(:,1))=1;
           pos_tf1=cat(2,pos_tf1,Of1.*(t*ones(Lgnome,1)));
           Of=Of+Of1;
        end
    end
    pos_oc=cat(2,Of,pos_oc1(:,2:end)); pos_tf=cat(2,Of,pos_tf1(:,2:end));
    pos_oc2=zeros(1,length(pos_oc(1,:))); pos_tf2=zeros(1,length(pos_tf(1,:))); posx=zeros(1,1);
    for i=1:length(Of)
        if Of(i)>0
           pos_oc2=cat(1,pos_oc2,pos_oc(i,:));
           pos_tf2=cat(1,pos_tf2,pos_tf(i,:));
           posx=cat(1,posx,i);
        end
    end
    fprintf('chr...%d \n',chr);
    pos_octf{chr,1}=posx(2:end,1);
    pos_octf{chr,2}=pos_oc2(2:end,:);
    pos_octf{chr,3}=pos_tf2(2:end,:);
end 
toc;
save /nuctf_equi_bai/NucRemod/toy_example/input/pos_octf/pos_octf.mat pos_octf -v7.3;

