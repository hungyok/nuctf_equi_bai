clear;
load D:\Cell_protocol\NucTF\simplexM_tf30\output2\simplex_xval.txt;
xfit=simplex_xval(:,1); TF=1; Eseq=1;
[Y]=occupx(TF,Eseq,xfit);