clear;
path1 ='yourpath\nuc_energy\'; 
path2 = 'yourpath\tf_energy_all\';    
path3 ='yourpath\ndr_call\';
load /NucTF/simplexM_top30/output2/simplex_xval.txt;
xfit = simplex_xval (:,1); TF = 1; Eseq = 1;
[O] = occupx(TF,Eseq,xfit,tfx,path1,path2,path3);
