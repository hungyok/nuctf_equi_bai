%clear;  
function [occup,ndr_chr]=ndr_pos_cal2(fname)
path1='nuctf_equi_bai/ndr_call/';
folder=fname;
fnx = fullfile(strcat(path1,folder),'yy1_lee.mat'); load(fnx);
%A=1; A1=2; B=3; C=4; D=5; E=6; F=7; G=8;
fnx = fullfile(strcat(path1,folder),'NDR_1.mat'); load(fnx); NDR_A=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_2.mat'); load(fnx); NDR_A1=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_3.mat'); load(fnx); NDR_B=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_4.mat'); load(fnx); NDR_C=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_5.mat'); load(fnx); NDR_D=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_6.mat'); load(fnx); NDR_E=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_7.mat'); load(fnx); NDR_F=NDR_X;
fnx = fullfile(strcat(path1,folder),'NDR_8.mat'); load(fnx); NDR_G=NDR_X;

fa=1; yytr=zeros(1,1); yyblu=zeros(1,1); SD=0.0678; ndrcnt=0; %fa=0.55;
AA=zeros(9,1); dxA_cf=100; dx_cf=50; dxkiA1_cf=125; dxkiB_cf=3*dxkiA1_cf; dxna_cf=110; ndra=0; ndrb=1025;
% fa=1 means occupancy at NDR is seto zero.
% dxA_cf is the slope cutoff which is the distance between cut-points at NDR_80 and at the lowest occupancy on the same side;
% dx_cf is the smallest NDR size cutoff below which is not an NDR; 
% dxkiA1_cf is the kink size (extreme ends of multiple kinks) at NDR_A1; dxkiB_cf is the kink size (extreme ends of multiple kinks) at NDR_B;
% dxna_cf is the distance between cut-points at NDR_A1; ndrb is the max NDR size;
% all these values are in bp.
ndr_chr=cell(16,1); dxkiA11=zeros(1,1); dxkiB1=zeros(1,1); dxna1=zeros(1,1); ndrw1=zeros(1,1); occup=cell(16,1);
for i=1:9
    AA(i,1)=0.8-(i-1)*SD;
end
for chr=1:16
    x1=x1_lee{chr}; yy1=y1_lee{chr}; yy1a=yy1; ndrj=0;
    for i=1:length(NDR_A{chr,1}(:,1))
        ndrcell=cell(7,1); cnt=0; cntA=0; jA=zeros(1,1); flgI=1;
        for j=1:7
            ndrcell{j,1}=zeros(1,1);
        end
        for j=1:length(NDR_A1{chr,1}(:,1))
            x=floor((NDR_A1{chr,1}(j,1)+NDR_A1{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0
           ndrcell{1,1}=jA(2:end,1); cnt=cnt+1; %....A1
        end
        cntA=0; jA=zeros(1,1); 
        for j=1:length(NDR_B{chr,1}(:,1))
            x=floor((NDR_B{chr,1}(j,1)+NDR_B{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0  
           ndrcell{2,1}=jA(2:end,1); cnt=cnt+1;    %....B
        end
        cntA=0; jA=zeros(1,1);
        for j=1:length(NDR_C{chr,1}(:,1))
            x=floor((NDR_C{chr,1}(j,1)+NDR_C{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0
           ndrcell{3,1}=jA(2:end,1); cnt=cnt+1;  %....C
        end
        cntA=0; jA=zeros(1,1);
        for j=1:length(NDR_D{chr,1}(:,1))
            x=floor((NDR_D{chr,1}(j,1)+NDR_D{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0
           ndrcell{4,1}=jA(2:end,1); cnt=cnt+1; %....D
        end
        cntA=0; jA=zeros(1,1);
        for j=1:length(NDR_E{chr,1}(:,1))
            x=floor((NDR_E{chr,1}(j,1)+NDR_E{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0
           ndrcell{5,1}=jA(2:end,1); cnt=cnt+1; %....E
        end
        cntA=0; jA=zeros(1,1);
        for j=1:length(NDR_F{chr,1}(:,1))
            x=floor((NDR_F{chr,1}(j,1)+NDR_F{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0
           ndrcell{6,1}=jA(2:end,1); cnt=cnt+1; %....F
        end
        cntA=0; jA=zeros(1,1);
        for j=1:length(NDR_G{chr,1}(:,1))
            x=floor((NDR_G{chr,1}(j,1)+NDR_G{chr,1}(j,2))/2);
            if x>=NDR_A{chr,1}(i,1) && x<=NDR_A{chr,1}(i,2)
               jA=cat(1,jA,j); cntA=cntA+1;
            end
        end
        if cntA>0
           ndrcell{7,1}=jA(2:end,1); cnt=cnt+1; %....G
        end 
        dx=NDR_A{chr,1}(i,2)-NDR_A{chr,1}(i,1); 
        if dx>dx_cf && cnt>0
           flgndr=1; j=ndrcell{1,1}(1,1);
           dxxAL=NDR_A1{chr,1}(j,1)-NDR_A{chr,1}(i,1); 
           j=ndrcell{1,1}(end,1);
           dxxAR=NDR_A{chr,1}(i,2)-NDR_A1{chr,1}(j,2); 
           if dxxAL<dxA_cf || dxxAR<dxA_cf
              if cnt==1
                 part=1;
                 j=ndrcell{1,1}(1,1);
                 xa=NDR_A1{chr,1}(j,1); 
                 j0=1;
                 while xa>x1(j0,1)
                       j0=j0+1;
                 end
                 ja=j0;
                 j=ndrcell{1,1}(end,1);
                 xb=NDR_A1{chr,1}(j,2);
                 j0=1;
                 while xb>x1(j0,1)
                       j0=j0+1;
                 end
                 jb=j0; Jxa=ja; Jxb=jb;
                 if length(ndrcell{1,1})>=2
                    jx1=ndrcell{1,1}(1,1); jx2=ndrcell{1,1}(end,1);
                    dx=NDR_A1{chr,1}(jx2,1)-NDR_A1{chr,1}(jx1,2); dxkiA11=cat(1,dxkiA11,dx);
                    if dx>dxkiA1_cf  % kink A ********
                       flgndr=0;              
                    end
                 end
                 if length(ndrcell{1,1})==1 
                    jx1=ndrcell{1,1}(1,1); 
                    dx=NDR_A1{chr,1}(jx1,2)-NDR_A1{chr,1}(jx1,1); dxna1=cat(1,dxna1,dx);
                    if dx<dxna_cf % narrow ********
                       flgndr=0;
                    end
                 end         
              else               
                 %left
                 part=2;
                 jcnt=2;
                 while jcnt<=cnt
                       j=ndrcell{jcnt,1}(1,1);
                       j0=ndrcell{jcnt-1,1}(1,1);
                       if jcnt>2
                          jn=ndrcell{jcnt-2,1}(1,1);
                       end
                       if jcnt==2
                          xa1=NDR_A{chr,1}(i,1); xa2=NDR_A1{chr,1}(j0,1); xa3=NDR_B{chr,1}(j,1); 
                       elseif jcnt==3
                          xa1=NDR_A1{chr,1}(jn,1); xa2=NDR_B{chr,1}(j0,1); xa3=NDR_C{chr,1}(j,1); 
                       elseif jcnt==4
                          xa1=NDR_B{chr,1}(jn,1); xa2=NDR_C{chr,1}(j0,1); xa3=NDR_D{chr,1}(j,1); 
                       elseif jcnt==5
                          xa1=NDR_C{chr,1}(jn,1); xa2=NDR_D{chr,1}(j0,1); xa3=NDR_E{chr,1}(j,1); 
                       elseif jcnt==6
                          xa1=NDR_D{chr,1}(jn,1); xa2=NDR_E{chr,1}(j0,1); xa3=NDR_F{chr,1}(j,1); 
                       else
                          xa1=NDR_E{chr,1}(jn,1); xa2=NDR_F{chr,1}(j0,1); xa3=NDR_G{chr,1}(j,1); 
                       end
                       dxL=xa2-xa1; 
                       mslope=SD/dxL;
                       x=(SD/mslope)+xa2;
                       if xa3<=(x+0)
                          dxL=xa3-xa2; 
                          mslope=SD/dxL;
                          x=(SD/mslope)+xa3;
                          if jcnt==cnt
                             xa=floor((xa3+x)/2); cntL=AA((jcnt+1)+1,1); jxa=xa3;
                          end
                       else
                          xa=floor((x+xa2)/2); cntL=AA(jcnt+1,1); jcnt=cnt; jxa=xa2;
                       end
                       jcnt=jcnt+1;
                 end
                 j0=1;
                 while xa>x1(j0,1)
                       j0=j0+1;
                 end
                 Ja=j0;
                 j0=1;
                 while jxa>x1(j0,1)
                       j0=j0+1;
                 end
                 Jxa=j0;
                 %right
                 jcnt=2;
                 while jcnt<=cnt
                       j=ndrcell{jcnt,1}(end,1);
                       j0=ndrcell{jcnt-1,1}(end,1);
                       if jcnt>2
                          jn=ndrcell{jcnt-2,1}(end,1);
                       end
                       if jcnt==2
                          xb1=NDR_A{chr,1}(i,2); xb2=NDR_A1{chr,1}(j0,2); xb3=NDR_B{chr,1}(j,2); 
                       elseif jcnt==3
                          xb1=NDR_A1{chr,1}(jn,2); xb2=NDR_B{chr,1}(j0,2); xb3=NDR_C{chr,1}(j,2); 
                       elseif jcnt==4
                          xb1=NDR_B{chr,1}(jn,2); xb2=NDR_C{chr,1}(j0,2); xb3=NDR_D{chr,1}(j,2); 
                       elseif jcnt==5
                          xb1=NDR_C{chr,1}(jn,2); xb2=NDR_D{chr,1}(j0,2); xb3=NDR_E{chr,1}(j,2); 
                       elseif jcnt==6
                          xb1=NDR_D{chr,1}(jn,2); xb2=NDR_E{chr,1}(j0,2); xb3=NDR_F{chr,1}(j,2); 
                       else
                          xb1=NDR_E{chr,1}(jn,2); xb2=NDR_F{chr,1}(j0,2); xb3=NDR_G{chr,1}(j,2); 
                       end
                       dxR=xb1-xb2; 
                       mslope=SD/dxR;
                       x=(-SD/mslope)+xb2;
                       if xb3>=(x-0)
                          dxR=xb2-xb3; 
                          mslope=SD/dxR;
                          x=(-SD/mslope)+xb3;
                          if jcnt==cnt
                             xb=floor((xb3+x)/2); cntR=AA((jcnt+1)+1,1); jxb=xb3;
                          end
                       else
                          xb=floor((x+xb2)/2); cntR=AA(jcnt+1,1); jcnt=cnt; jxb=xb2;
                       end
                       jcnt=jcnt+1;
                 end
                 j0=1;
                 while xb>x1(j0,1)
                       j0=j0+1;
                 end
                 Jb=j0; j0=1;
                 while jxb>x1(j0,1)
                       j0=j0+1;
                 end
                 Jxb=j0; 
                 if length(ndrcell{2,1})>=2
                    kx1=ndrcell{2,1}(1,1); kx2=ndrcell{2,1}(end,1);
                    dx=NDR_B{chr,1}(kx2,1)-NDR_B{chr,1}(kx1,2); dxkiB1=cat(1,dxkiB1,dx);
                    if dx>dxkiB_cf  % kink B ********
                       flgndr=0;   % no transformation
                    end
                 end
                 if length(ndrcell{1,1})==1 
                    jx1=ndrcell{1,1}(1,1); 
                    dx=NDR_A1{chr,1}(jx1,2)-NDR_A1{chr,1}(jx1,1); dxna1=cat(1,dxna1,dx);
                    if dx<dxna_cf % narrow ********
                       flgndr=0;
                    end
                 end
                 if length(ndrcell{1,1})>=2
                    part=11;
                    jx1=ndrcell{1,1}(1,1); jx2=ndrcell{1,1}(end,1);
                    dx=NDR_A1{chr,1}(jx2,1)-NDR_A1{chr,1}(jx1,2); %dxki1=cat(1,dxki1,dx);
                    if dx>dxkiB_cf  % kink B ********
                       flgndr=0;              
                    end
                    j0=1;
                    while NDR_A1{chr,1}(jx1,2)>x1(j0,1)
                          j0=j0+1;
                    end
                    ja=j0; j0=1;
                    while NDR_A1{chr,1}(jx2,1)>x1(j0,1)
                          j0=j0+1;
                    end
                    jb=j0;
                    if (0.8-max(yy1(ja:jb,1)))<SD/3 %..............below 0.8.........................
                       if dx<dxkiA1_cf  % kink A ********
                          flgndr=1; part=3;         
                       else 
                           part=4;
                          [m,I]=max(yy1(ja:jb,1)); Imx=ja+I-1; flgI=0;
                          for i1=1:2
                              if i1==1
                                 NDR_Aa=NDR_A{chr,1}(i,1); NDR_Ab=x1(Imx,1);
                              else
                                 NDR_Aa=x1(Imx,1); NDR_Ab=NDR_A{chr,1}(i,2);
                              end
                              cnt=0; cntA=0; jA=zeros(1,1); 
                              for j=1:7
                                  ndrcell{j,1}=zeros(1,1);
                              end
                              for j=1:length(NDR_A1{chr,1}(:,1))
                                  x=floor((NDR_A1{chr,1}(j,1)+NDR_A1{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0
                                 ndrcell{1,1}=jA(2:end,1); cnt=cnt+1; %....A1
                              end
                              cntA=0; jA=zeros(1,1); 
                              for j=1:length(NDR_B{chr,1}(:,1))
                                  x=floor((NDR_B{chr,1}(j,1)+NDR_B{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0  
                                 ndrcell{2,1}=jA(2:end,1); cnt=cnt+1;    %....B
                              end
                              cntA=0; jA=zeros(1,1);
                              for j=1:length(NDR_C{chr,1}(:,1))
                                  x=floor((NDR_C{chr,1}(j,1)+NDR_C{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0
                                 ndrcell{3,1}=jA(2:end,1); cnt=cnt+1;  %....C
                              end
                              cntA=0; jA=zeros(1,1);
                              for j=1:length(NDR_D{chr,1}(:,1))
                                  x=floor((NDR_D{chr,1}(j,1)+NDR_D{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0
                                 ndrcell{4,1}=jA(2:end,1); cnt=cnt+1; %....D
                              end
                              cntA=0; jA=zeros(1,1);
                              for j=1:length(NDR_E{chr,1}(:,1))
                                  x=floor((NDR_E{chr,1}(j,1)+NDR_E{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0
                                 ndrcell{5,1}=jA(2:end,1); cnt=cnt+1; %....E
                              end
                              cntA=0; jA=zeros(1,1);
                              for j=1:length(NDR_F{chr,1}(:,1))
                                  x=floor((NDR_F{chr,1}(j,1)+NDR_F{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0
                                 ndrcell{6,1}=jA(2:end,1); cnt=cnt+1; %....F
                              end
                              cntA=0; jA=zeros(1,1);
                              for j=1:length(NDR_G{chr,1}(:,1))
                                  x=floor((NDR_G{chr,1}(j,1)+NDR_G{chr,1}(j,2))/2);
                                  if x>=NDR_Aa && x<=NDR_Ab
                                     jA=cat(1,jA,j); cntA=cntA+1;
                                  end
                              end
                              if cntA>0
                                 ndrcell{7,1}=jA(2:end,1); cnt=cnt+1; %....G
                              end 
                              dx=NDR_Ab-NDR_Aa; 
                              if dx>dx_cf && cnt>0
                                 flgndr=1;
                                 j=ndrcell{1,1}(1,1);
                                 dxxAL=NDR_A1{chr,1}(j,1)-NDR_Aa; 
                                 j=ndrcell{1,1}(end,1);
                                 dxxAR=NDR_Ab-NDR_A1{chr,1}(j,2); 
                                 if dxxAL<dxA_cf || dxxAR<dxA_cf
                                    if cnt==1
                                       part=5;
                                       j=ndrcell{1,1}(1,1);
                                       xa=NDR_A1{chr,1}(j,1); 
                                       j0=1;
                                       while xa>x1(j0,1)
                                             j0=j0+1;
                                       end
                                       ja=j0;
                                       j=ndrcell{1,1}(end,1);
                                       xb=NDR_A1{chr,1}(j,2);
                                       j0=1;
                                       while xb>x1(j0,1)
                                             j0=j0+1;
                                       end
                                       jb=j0; Jxa=ja; Jxb=jb;
                                       if length(ndrcell{1,1})>=2
                                          jx1=ndrcell{1,1}(1,1); jx2=ndrcell{1,1}(end,1);
                                          dx=NDR_A1{chr,1}(jx2,1)-NDR_A1{chr,1}(jx1,2); dxkiA11=cat(1,dxkiA11,dx);
                                          if dx>dxkiA1_cf  % kink A********
                                             flgndr=0;              
                                          end
                                       end
                                       if length(ndrcell{1,1})==1 
                                          jx1=ndrcell{1,1}(1,1); 
                                          dx=NDR_A1{chr,1}(jx1,2)-NDR_A1{chr,1}(jx1,1); dxna1=cat(1,dxna1,dx);
                                          if dx<dxna_cf % narrow ********
                                             flgndr=0;
                                          end
                                       end
                                    else
                                       %left
                                       part=6;
                                       jcnt=2;
                                       while jcnt<=cnt
                                             j=ndrcell{jcnt,1}(1,1);
                                             j0=ndrcell{jcnt-1,1}(1,1);
                                             if jcnt>2
                                                jn=ndrcell{jcnt-2,1}(1,1);
                                             end
                                             if jcnt==2
                                                xa1=NDR_Aa; xa2=NDR_A1{chr,1}(j0,1); xa3=NDR_B{chr,1}(j,1); 
                                             elseif jcnt==3
                                                xa1=NDR_A1{chr,1}(jn,1); xa2=NDR_B{chr,1}(j0,1); xa3=NDR_C{chr,1}(j,1); 
                                             elseif jcnt==4
                                                xa1=NDR_B{chr,1}(jn,1); xa2=NDR_C{chr,1}(j0,1); xa3=NDR_D{chr,1}(j,1); 
                                             elseif jcnt==5
                                                xa1=NDR_C{chr,1}(jn,1); xa2=NDR_D{chr,1}(j0,1); xa3=NDR_E{chr,1}(j,1); 
                                             elseif jcnt==6
                                                xa1=NDR_D{chr,1}(jn,1); xa2=NDR_E{chr,1}(j0,1); xa3=NDR_F{chr,1}(j,1); 
                                             else
                                                xa1=NDR_E{chr,1}(jn,1); xa2=NDR_F{chr,1}(j0,1); xa3=NDR_G{chr,1}(j,1); 
                                             end
                                             dxL=xa2-xa1; 
                                             mslope=SD/dxL;
                                             x=(SD/mslope)+xa2;
                                             if xa3<=(x+0)
                                                dxL=xa3-xa2; 
                                                mslope=SD/dxL;
                                                x=(SD/mslope)+xa3;
                                                if jcnt==cnt
                                                   xa=floor((xa3+x)/2); cntL=AA((jcnt+1)+1,1); jxa=xa3;
                                                end
                                             else
                                                xa=floor((x+xa2)/2); cntL=AA(jcnt+1,1); jcnt=cnt; jxa=xa2;
                                             end
                                             jcnt=jcnt+1;
                                       end
                                       j0=1;
                                       while xa>x1(j0,1)
                                             j0=j0+1;
                                       end
                                       Ja=j0; j0=1;
                                       while jxa>x1(j0,1)
                                             j0=j0+1;
                                       end
                                       Jxa=j0;
                                       %right
                                       jcnt=2;
                                       while jcnt<=cnt
                                             j=ndrcell{jcnt,1}(end,1);
                                             j0=ndrcell{jcnt-1,1}(end,1);
                                             if jcnt>2
                                                jn=ndrcell{jcnt-2,1}(end,1);
                                             end
                                             if jcnt==2
                                                xb1=NDR_Ab; xb2=NDR_A1{chr,1}(j0,2); xb3=NDR_B{chr,1}(j,2); 
                                             elseif jcnt==3
                                                xb1=NDR_A1{chr,1}(jn,2); xb2=NDR_B{chr,1}(j0,2); xb3=NDR_C{chr,1}(j,2); 
                                             elseif jcnt==4
                                                xb1=NDR_B{chr,1}(jn,2); xb2=NDR_C{chr,1}(j0,2); xb3=NDR_D{chr,1}(j,2); 
                                             elseif jcnt==5
                                                xb1=NDR_C{chr,1}(jn,2); xb2=NDR_D{chr,1}(j0,2); xb3=NDR_E{chr,1}(j,2); 
                                             elseif jcnt==6
                                                xb1=NDR_D{chr,1}(jn,2); xb2=NDR_E{chr,1}(j0,2); xb3=NDR_F{chr,1}(j,2); 
                                             else
                                                xb1=NDR_E{chr,1}(jn,2); xb2=NDR_F{chr,1}(j0,2); xb3=NDR_G{chr,1}(j,2); 
                                             end
                                             dxR=xb1-xb2; 
                                             mslope=SD/dxR;
                                             x=(-SD/mslope)+xb2;
                                             if xb3>=(x-0)
                                                dxR=xb2-xb3; 
                                                mslope=SD/dxR;
                                                x=(-SD/mslope)+xb3;
                                                if jcnt==cnt
                                                   xb=floor((xb3+x)/2); cntR=AA((jcnt+1)+1,1); jxb=xb3;
                                                end
                                             else
                                                xb=floor((x+xb2)/2); cntR=AA(jcnt+1,1); jcnt=cnt; jxb=xb2;
                                             end
                                             jcnt=jcnt+1;
                                       end
                                       j0=1;
                                       while xb>x1(j0,1)
                                             j0=j0+1;
                                       end
                                       Jb=j0; j0=1;
                                       while jxb>x1(j0,1)
                                             j0=j0+1;
                                       end
                                       Jxb=j0; 
                                       if length(ndrcell{1,1})>=2
                                          jx1=ndrcell{1,1}(1,1); jx2=ndrcell{1,1}(end,1);
                                          dx=NDR_A1{chr,1}(jx2,1)-NDR_A1{chr,1}(jx1,2); %dxki1=cat(1,dxki1,dx);
                                          if dx>dxkiB_cf % kink B ********
                                             flgndr=0;              
                                          end
                                       end
                                       if length(ndrcell{2,1})>=2
                                          kx1=ndrcell{2,1}(1,1); kx2=ndrcell{2,1}(end,1);
                                          dx=NDR_B{chr,1}(kx2,1)-NDR_B{chr,1}(kx1,2); dxkiB1=cat(1,dxkiB1,dx);
                                          if dx>dxkiB_cf % kink B********
                                             flgndr=0;
                                          end
                                       end
                                       if length(ndrcell{1,1})==1 
                                          jx1=ndrcell{1,1}(1,1); 
                                          dx=NDR_A1{chr,1}(jx1,2)-NDR_A1{chr,1}(jx1,1); dxna1=cat(1,dxna1,dx);
                                          if dx<dxna_cf % narrow ********
                                             flgndr=0;
                                          end
                                       end
                                       ja=Ja; jb=Jb;
                                    end
                                    if flgndr>0
                                       if (x1(jb,1)-x1(ja,1))<ndrb && (x1(jb,1)-x1(ja,1))>ndra
                                          djxa=ja-Jxa; djxb=Jxb-jb;
                                          if djxa>1
                                             djxa=x1(ja,1)-x1(Jxa,1);
                                             mslope=-SD/djxa;
                                             yy1a(Jxa+1:ja-1,1)=mslope*(x1(Jxa+1:ja-1,1)-x1(Jxa,1))+cntL+SD;
                                          end
                                          if djxb>1
                                             djxb=x1(Jxb,1)-x1(jb,1);
                                             mslope=SD/djxb;
                                             yy1a(jb+1:Jxb-1,1)=mslope*(x1(jb+1:Jxb-1,1)-x1(Jxb,1))+cntR+SD;
                                          end
                                          yy1a(ja:jb,1)=(1-fa)*(1-yy1(ja:jb,1)); %apply nfr;
                                          ndrcnt=ndrcnt+1; ndrj=ndrj+1;
                                          ndr_chr{chr,1}(ndrj,1)=x1(ja,1);
                                          ndr_chr{chr,1}(ndrj,2)=x1(jb,1);
                                          ndr_chr{chr,1}(ndrj,3)=fa;
                                          if (ndrj-1)>0
                                             dxz=ndr_chr{chr,1}(ndrj,1)-ndr_chr{chr,1}(ndrj-1,2);
                                             if dxz<dxkiA1_cf
                                                yy1a(jb0:ja,1)=(1-fa)*(1-yy1(jb0:ja,1)); 
                                                ndrcnt=ndrcnt-1; ndrj=ndrj-1;
                                                ndr_chr{chr,1}(ndrj,1)=x1(ja0,1);
                                                ndr_chr{chr,1}(ndrj,2)=x1(jb,1);
                                                ndr_chr{chr,1}(ndrj,3)=fa;
                                                ja=ja0; 
                                                ndrw1=cat(1,ndrw1,x1(jb,1)-x1(ja,1)); 
                                             else
                                                ndrw1=cat(1,ndrw1,x1(jb,1)-x1(ja,1)); 
                                             end
                                          else
                                             ndrw1=cat(1,ndrw1,x1(jb,1)-x1(ja,1));  
                                          end
                                          ja0=ja; jb0=jb;
                                          %fprintf('below 0.8...what the hell...i...%d...i1...%d\n',i,i1);
                                       end   
                                    end
                                 end
                              end
                          end
                       end
                    end
                 end
                 ja=Ja; jb=Jb;
              end
              if flgndr>0 && flgI>0
                 if (x1(jb,1)-x1(ja,1))<ndrb && (x1(jb,1)-x1(ja,1))>ndra
                    djxa=ja-Jxa; djxb=Jxb-jb;
                    if djxa>1
                       djxa=x1(ja,1)-x1(Jxa,1);
                       mslope=-SD/djxa;
                       yy1a(Jxa+1:ja-1,1)=mslope*(x1(Jxa+1:ja-1,1)-x1(Jxa,1))+cntL+SD;
                    end
                    if djxb>1
                       djxb=x1(Jxb,1)-x1(jb,1);
                       mslope=SD/djxb;
                       yy1a(jb+1:Jxb-1,1)=mslope*(x1(jb+1:Jxb-1,1)-x1(Jxb,1))+cntR+SD;
                    end
                    yy1a(ja:jb,1)=(1-fa)*(1-yy1(ja:jb,1)); %apply nfr;
                    ndrcnt=ndrcnt+1; ndrj=ndrj+1;
                    ndr_chr{chr,1}(ndrj,1)=x1(ja,1);
                    ndr_chr{chr,1}(ndrj,2)=x1(jb,1);
                    ndr_chr{chr,1}(ndrj,3)=fa;
                    if (ndrj-1)>0
                       dxz=ndr_chr{chr,1}(ndrj,1)-ndr_chr{chr,1}(ndrj-1,2);
                       if dxz<dxkiA1_cf
                          yy1a(jb0:ja,1)=(1-fa)*(1-yy1(jb0:ja,1)); 
                          ndrcnt=ndrcnt-1; ndrj=ndrj-1;
                          ndr_chr{chr,1}(ndrj,1)=x1(ja0,1);
                          ndr_chr{chr,1}(ndrj,2)=x1(jb,1);
                          ndr_chr{chr,1}(ndrj,3)=fa;
                          ja=ja0; 
                          ndrw1=cat(1,ndrw1,x1(jb,1)-x1(ja,1)); 
                       else
                          ndrw1=cat(1,ndrw1,x1(jb,1)-x1(ja,1)); 
                       end
                    else
                       ndrw1=cat(1,ndrw1,x1(jb,1)-x1(ja,1));  
                    end
                    ja0=ja; jb0=jb;
                    %fprintf('what the hell...i....%d\n',i);
                 end
              end
           end 
        end
    end
    fprintf('chr....%d\n',chr);
    %yytr{chr,1}=yy1a;
    occup{chr,1}=yy1a;
    %yyblu=cat(1,yyblu,yy3);
end
%save /nuctf_equi_bai/ndr_call/yy3A.mat occup -v7.3;
%save /nuctf_equi_bai/ndr_call/ndrpos_chrA.mat ndr_chr -v7.3;
end
