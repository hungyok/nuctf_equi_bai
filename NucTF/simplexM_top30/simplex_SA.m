% Optimization
function simplex_SA(path1,path2,path3)
foldername1='simplexM_top30'; foldername2='Etf_allmat_chr'; foldername3='input'; foldername4='output'; foldername5='initialp';
fnx = fullfile(strcat(path3,foldername1),foldername4,'simplex_xval.txt');
fnxbak = fullfile(strcat(path3,foldername1),foldername4,'simplex_xvalbak.txt');
fnf = fullfile(strcat(path3,foldername1),foldername4,'simplex_fval.txt');
fnt = fullfile(strcat(path3,foldername1),foldername4,'simplex_tval.txt');
fnT = fullfile(strcat(path3,foldername1),foldername4,'simplex_Temp.txt');
fileID = fopen(fnT,'w'); 
fpath= strcat(path1,'E_Em.mat');
load(fpath);
fpath= strcat(path2,'tfindx.txt');
load(fpath);
fpath=fullfile(strcat(path2,foldername2),'Emtfall.mat');
load(fpath);
Etfmul=cell(1,16); tfn=30; contd=0; % if contd=1 rerun with last output.
for k =1:16
    filename = sprintf('Etf_chr%d.mat',k);
    fpath=fullfile(strcat(path2,foldername2),filename);
    load(fpath);
    Etfmul{1,k}=Etf(:,tfindx(1:tfn));
end
n=(tfn+1)*2; Tmin=0.00001; alfa=0.9; err=0.000001; clear tfindx Etf;
x=zeros(n,n+1); xr=zeros(n,n+1); xnew=zeros(n,n+1); Fsig=zeros(n+1,1); Fsigxr=zeros(n+1,1); rng('shuffle'); 
h=zeros(1,n); e=eye(n);
if contd>0
   load(fnx);
   xh=simplex_xval(:,1);
   x(1:n,1)=xh; 
   for i=1:n
       if xh(i)<0.00001
          h(i)=0.00025;
       else
          h(i)=0.05*xh(i);
       end
   end
   for i=2:(n+1)
       x(:,i)=x(:,1)+h(i-1)*e(:,i-1);    
   end  
else
   fnx1 = fullfile(strcat(path3,foldername1),foldername3,foldername5,'top30_xhIIa.mat');
   load(fnx1);
   x(1:n,1)=xh(1,1:n); h(1:n)=xh(1,(n+1):2*n);
   for i=2:(n+1)
       x(:,i)=x(:,1)+h(i-1)*e(:,i-1);    
   end 
end
xi=x;
for i=1:(n+1) 
    Fsig(i,1)=optimbfunc_ver3a(xi(:,i),E,Em,Etfmul,Emtf,path2,path3); fprintf('print i....%d \n',i); 
end
[Fi,Fsig]=sortf(Fsig,n);
for i=1:n+1
    xnew(:,i)=x(:,Fi(i,1));
end
x=xnew;
ter=Fsig(n+1,1)-Fsig(1,1); Tmax=-ter/log(0.9); Temp=Tmax; fprintf('ter...%f Temp...%f Fsig1...%f \n',ter,Temp,Fsig(1,1));
if Temp>0.01
   M=100;
elseif Temp<=0.01 && Temp>0.001
   M=83;
elseif Temp<=0.001 && Temp>0.0001
   M=67;
elseif Temp<=0.0001 && Temp>0.00002
   M=50;
else
   M=33;
end 
fprintf(fileID,'ter...%f...Temp...%f...Fsig(1,1)...%f \n',ter,Temp,Fsig(1,1)); Tempt=[ter Temp Fsig(1,1)];
dlmwrite(fnx,x,'precision','%.6f'); dlmwrite(fnf,Fsig','precision','%.6f'); dlmwrite(fnt,Tempt,'precision','%.6f');
if ter<=err || Temp<Tmin
   flgx=0;
   fprintf('whileloop ter...%f Temp...%f Fsig1...%f \n',ter,Temp,Fsig(1,1));
else
   flgx=1;
   while flgx>0
       for mx=1:M
           k=1; flgx1=1; count=0;
           while flgx1>0
                 if k<=6 % k<=n for low n
                    sumx0=zeros(n,1);
                    for ki=1:(n-k+1)
                        sumx0=sumx0+x(:,ki);
                    end
                    x0=sumx0/(n-k+1); % x0=(x(:,1)+x(:,2))/2; %centriod
                    for ki=(n-k+2):(n+1)
                        ro=0.9+rand*(1.1-0.9); 
                        xr(:,ki)=x0+ro*(x0-x(:,ki)); %disp(xi);          
                        xi=xr(:,ki);
                        Fsigxr(ki,1)=optimbfunc_ver3a(xi,E,Em,Etfmul,Emtf,path2,path3);
                    end
                    fm=min(Fsigxr((n-k+2):(n+1),1));
                    if fm<Fsig(1,1)
                       x(:,(n-k+2):n+1)=xr(:,(n-k+2):n+1); Fsig((n-k+2):n+1,1)=Fsigxr((n-k+2):n+1,1); [Fi,Fsig]=sortf(Fsig,n);
                       for i=1:n+1
                           xnew(:,i)=x(:,Fi(i,1));
                       end
                       x=xnew;
                       flgx1=0;  fprintf('reflect accept1.... \n'); %accept;
                    else
                       if fm>1
                          Px=0; fprintf('reflect wow! high objfunction!.... \n');
                       else
                          Px=exp(-(fm-Fsig(1,1))/Temp);
                       end
                       if rand<Px
                          x(:,(n-k+2):n+1)=xr(:,(n-k+2):n+1); Fsig((n-k+2):n+1,1)=Fsigxr((n-k+2):n+1,1); [Fi,Fsig]=sortf(Fsig,n);
                          for i=1:n+1
                              xnew(:,i)=x(:,Fi(i,1));
                          end
                          x=xnew;
                          flgx1=0; fprintf('reflect accept2.... \n');   %accept;
                       else
                          k=k+1; 
                       end
                    end
                 else
                    for i=2:n+1 % shrink
                        x(:,i)=x(:,1)+0.5*(x(:,i)-x(:,1));
                    end
                    xi=x(:,2:n+1);
                    for i=2:n+1  
                        Fsig(i,1)=optimbfunc_ver3a(xi(:,i-1),E,Em,Etfmul,Emtf,path2,path3);
                    end
                    [Fi,Fsig]=sortf(Fsig,n); flgx1=0; fprintf('shrink.... \n');
                    for i=1:n+1
                        xnew(:,i)=x(:,Fi(i,1));
                    end
                    x=xnew; dlmwrite(fnx,x,'precision','%.6f'); 
                 end
                 count=count+1; fprintf('mx %d count %d \n',mx,count); 
           end
       end
       ter=Fsig(n+1,1)-Fsig(1,1); fprintf('whileloop ter...%f Temp...%f Fsig1...%f \n',ter,Temp,Fsig(1,1)); 
       fprintf(fileID,'ter...%f...Temp...%f...Fsig(1,1)...%f \n',ter,Temp,Fsig(1,1)); Tempt=[ter Temp Fsig(1,1)];
       dlmwrite(fnx,x,'precision','%.6f'); dlmwrite(fnxbak,x,'precision','%.6f'); dlmwrite(fnf,Fsig','precision','%.6f'); dlmwrite(fnt,Tempt,'precision','%.6f');
       if ter<=err || Temp<Tmin
          flgx=0;
       end
       Temp=alfa*Temp;
       if Temp>0.01
          M=100;
       elseif Temp<=0.01 && Temp>0.001
          M=83;
       elseif Temp<=0.001 && Temp>0.0001
          M=67;
       elseif Temp<=0.0001 && Temp>0.00002
          M=50;
       else
          M=33;
       end
   end
end
fclose(fileID);
end        
