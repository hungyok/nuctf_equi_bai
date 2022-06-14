clear;
load D:\Cell_protocol\NucRemod\simplexM_remod\output2\ismplx1\simplex_xval.txt;
remod=simplex_xval(:,1); 
load D:\Cell_protocol\NucTF\simplexM_tf30\output2\simplex_xval.txt;
xfit=simplex_xval(:,1); 
Y=occupR(remod,xfit);