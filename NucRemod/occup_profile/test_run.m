clear;
path1 ='yourpath\nuc_energy\'; 
path2 = 'yourpath\tf_energy_all\';    
path3 ='yourpath\NucRemod\simplexM_remod\';
path4 ='yourpath\ndr_call\';
load /NucRemod/simplexM_remod/output2/simplex_xval.txt;
remod = simplex_xval(:,1); % load the best-fit remodeling parameters
load /NucTF/simplexM_top30/output2/simplex_xval.txt;
xfit = simplex_xval(:,1); % load the best-fit TF parameters
[O] = occupR(remod,xfit,tfx,path1,path2,path3,path4); 
% save O.mat O;
