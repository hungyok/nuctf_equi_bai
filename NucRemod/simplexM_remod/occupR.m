% Remodeling and rmsd
function [sigma]=occupR(xi)
load /nuctf_equi_bai/nuc_energy/E_Em.mat;
load /nuctf_equi_bai/NucTF/simplexM_tf30/input/rand_genomeB.mat;
load /nuctf_equi_bai/NucTF/simplexM_tf30/input/yy3A_lee.mat;
load /nuctf_equi_bai/NucTF/simplexM_tf30/output2/simplex_xval.txt;
xfit=simplex_xval(:,1); tfnx=floor(length(xfit)/2)-1;
xfit(1,1)=xi(1,1); xfit(2,1)=xi(2,1);
invld=0;
if xi(1,1)>100 || xi(1,1)<0
   invld=1;
end
if xi(2,1)>5 || xi(2,1)<0
   invld=1;
end
if xi(3,1)>25 || xi(3,1)<0
   invld=1;
end
if xi(4,1)>200 || xi(4,1)<1
   invld=1;
end
if invld==1
   sigma=5;
else  
   TFlist=importdata('/nuctf_equi_bai/tf_energy_all/listbai_all.txt'); % list of TF name and size 
   TF_list1=TFlist.data; clear TFlist;
   load /nuctf_equi_bai/tf_energy_all/tfindx.txt; % ranked listbai_all.txt indices: 1st or top (high NDR predictor) to 104th or bottom (low NDR predictor) TF.
   TF_list=TF_list1(tfindx(1:tfnx),1);
   load /nuctf_equi_bai/NucRemod/simplexM_remod/input/pos_octf/pos_octf.mat; % TF occupancy & position at occupancy tfcut=0.0022.
   pos_octf1=pos_octf;
   pos_octf=cell(1,3); idx=1;
   nx=10; ht=xi(3,1)/0.212; xd=xi(4,1)-floor(xi(4,1)); 
   if xd>0.5
      xds=floor(xi(4,1))+1;
   else
      xds=floor(xi(4,1));
   end
   sig=xds; A=ht; B=1/(2*sig^2); occup_nuc=cell(chrindx(idx,1),1); tic;
   for ix=1:nx % # realization
       for chr=1:chrindx(idx,1) % indie_test 2
           chr1=chrlist{idx,1}(chr);
           E1=E{chr1};
           pos_octf{1,1}=pos_octf1{chr1,1};
           pos_octf{1,2}=pos_octf1{chr1,2};
           pos_octf{1,3}=pos_octf1{chr1,3};
           [xptfa]=tf_cluster(tfnx,pos_octf); xptf=xptfa;
           lx=length(E{chr1}); xx=1:lx; x=xx';
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
           occup_nuc{chr,1}=occup_nucs(xfit,E1,Em,chr1);
       end
       if ix<2
          occup_nuc0=occup_nuc;  
       else
          for chr=1:chrindx(idx,1)
              occup_nuc0{chr,1}=occup_nuc0{chr,1}+occup_nuc{chr,1};
          end
       end
   end
   fprintf('function occupR toptf...%d...nx...%d\n',tfnx,nx); toc;
   for chr=1:chrindx(idx,1)
       chr1=chrlist{idx,1}(chr);
       occup_nuc{chr,1}=occup_nuc0{chr,1}*(1/nx);
       if chr==1
          ONL=occup_nuc{chr,1}(x1_lee{chr1,1}(100:end-100)); yL=y1_lee{chr1,1}(100:end-100);
       else
          ONL=cat(1,ONL,occup_nuc{chr,1}(x1_lee{chr1,1}(100:end-100))); yL=cat(1,yL,y1_lee{chr1,1}(100:end-100));
       end
   end
   Oer=ONL-yL;
   sigma=sqrt((Oer'*Oer)/length(Oer));
end
end
