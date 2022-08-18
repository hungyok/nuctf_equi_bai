% a test run on how to load and run the code occupR.m
clear;
load /nuctf_equi_bai/NucRemod/simplexM_remod/output2/simplex_xval.txt;
remod=simplex_xval(:,1); 
load /nuctf_equi_bai/NucTF/simplexM_tf30/output2/simplex_xval.txt;
xfit=simplex_xval(:,1); 
Y=occupR(remod,xfit);
