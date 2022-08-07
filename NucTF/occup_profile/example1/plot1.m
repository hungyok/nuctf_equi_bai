clear;
load /nuctf_equi_bai/NucTF/occup_profile/example1/O_noseq.mat;
y1=Y{1};
load /nuctf_equi_bai/NucTF/occup_profile/example1/O_seq.mat;
y2=Y{1};

clf;
x=1:length(y1); 
figure('Renderer', 'painters','units','centimeters', 'Position', [0, 0, 25, 8]);
plot (x,y1,'color','k','DisplayName','bkg','LineWidth',1);
hold on;
plot (x,y2,'color','r','DisplayName','bkg','LineWidth',1);
hold on
axis([100000 105000,0.6 1])
