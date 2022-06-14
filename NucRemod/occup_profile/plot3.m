clear;
load D:\Cell_protocol\ndr_call\yy3A_lee.mat;
xlee=x1_lee{1};
ylee=y1_lee{1};
load D:\Cell_protocol\NucTF\occup_profile\example2\O_good.mat;
y1=Y{1}(:,1);
load D:\Cell_protocol\NucRemod\occup_profile\O.mat;
y2=occup_nuc{1}(:,1);
clf;
x=1:length(y1); 
figure('Renderer', 'painters','units','centimeters', 'Position', [0, 0, 25, 8]);
plot (xlee,ylee,'color','k','DisplayName','bkg','LineWidth',1.5);
hold on;
plot (x,y1,'color','r','DisplayName','bkg','LineWidth',1.5);
hold on;
plot (x,y2,'color','b','DisplayName','bkg','LineWidth',1.5);
hold on;

axis([100000 105000,0 1]);
