%clear;
%function [Emtf,k]=tf_binding_pot(path1,path2)
%scan motifs and calculate binding energy
    TFdata=importdata('D:\Cell_protocol\tf_energy_all\wmsbai_data_all.txt');
    TFlist=importdata('D:\Cell_protocol\tf_energy_all\listbai_all.txt');
    path1='D:\Cell_protocol\tf_energy_all\sgd_genome';
    path2='D:\Cell_protocol\tf_energy_all\Etf_allmat_chr\';
    Etf1=cell(1,16);
for k=1:1
    fname = sprintf('sgd_chr%d.fa',k);
    Seqdata=importdata(fullfile(path1,fname)); chr=Seqdata{2,1};
    for i=3:length(Seqdata)
        chr=strcat(chr,Seqdata{i,1}); 
    end
    chrcp=seqcomplement(chr);
    L=length(chr);
    sumt=0; Etf=zeros(L,0);
    for t=1:length(TFlist.data) % T #transcription factor
        ks=TFlist.data(t); Etf0=zeros(L,3); 
        for i=1:(L-ks+1) % DNA length 
            w1=zeros(TFlist.data(t),1); w2=zeros(TFlist.data(t),1); w3=zeros(TFlist.data(t),1); w4=zeros(TFlist.data(t),1); 
            i1r=zeros(2,1);
            for i1=0:(TFlist.data(t)-1)
                %1strand
                if chr(i+i1)=='A'
                   w1(i1+1,1)=TFdata(i1+1+sumt,2);
                elseif chr(i+i1)=='C'
                   w1(i1+1,1)=TFdata(i1+1+sumt,3);
                elseif chr(i+i1)=='G'
                   w1(i1+1,1)=TFdata(i1+1+sumt,4);
                else   
                   w1(i1+1,1)=TFdata(i1+1+sumt,5);
                end
                %2strand
                if chrcp(i+i1)=='A'
                   w4(i1+1,1)=TFdata(1+(TFlist.data(t)-1)-i1+sumt,2);
                elseif chrcp(i+i1)=='C'
                   w4(i1+1,1)=TFdata(1+(TFlist.data(t)-1)-i1+sumt,3);
                elseif chrcp(i+i1)=='G'
                   w4(i1+1,1)=TFdata(1+(TFlist.data(t)-1)-i1+sumt,4);
                else   
                   w4(i1+1,1)=TFdata(1+(TFlist.data(t)-1)-i1+sumt,5);
                end
            end
            i1r(1,1)=sum(w1); i1r(2,1)=sum(w4); 
            i1=max(i1r); Etf0(i,1)=i1r(1,1); Etf0(i,2)=i1r(2,1); Etf0(i,3)=i1;
        end
        if t==2
           Etf2=Etf0; 
        end
        %Etf=cat(2,Etf,Etf0);
        sumt=sumt+TFlist.data(t); 
    end
    %Etf1{1,k}=Etf(:,2:105);
    fprintf('chr...%d \n',k);
end

%fname = sprintf('Emtfall.mat');
%fnx=strcat(path2,fname);
%save(fnx,'Emtf','-v7.3');
%end
