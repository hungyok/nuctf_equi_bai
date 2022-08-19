clear;
load /nuctf_equi_bai/ndr_call/yy3A_lee.mat;
xlee=x1_lee{1};
ylee=y1_lee{1};
load /nuctf_equi_bai/NucTF/occup_profile/example2/O_good.mat;
y1=Y{1}(:,1);
y2=Y{1}(:,2);
y3=Y{1}(:,3);
y4=Y{1}(:,4);
y5=Y{1}(:,5);

clf;
x=1:length(y1); 
figure('Renderer', 'painters','units','centimeters', 'Position', [0, 0, 25, 8]);
plot (xlee,ylee,'color','k','DisplayName','bkg','LineWidth',1.5);
hold on;
plot (x,y1,'color','r','DisplayName','bkg','LineWidth',1.5);
hold on
plot (x,y2,'color','b','DisplayName','bkg','LineWidth',1.5);
hold on
plot (x,y3,'color','c','DisplayName','bkg','LineWidth',1.5);
hold on
plot (x,y4,'color','m','DisplayName','bkg','LineWidth',1.5);
hold on
plot (x,y5,'color','g','DisplayName','bkg','LineWidth',1.5);
hold on
axis([100000 105000,0 1])
