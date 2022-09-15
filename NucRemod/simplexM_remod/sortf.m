% sorting
function [Fi,Fsig]=sortf(Fsig,n)
Fi=1:n+1; Fi=Fi'; 
i=1; 
while i<(n+1)
      for j=(i+1):(n+1)
          if Fsig(i,1)>Fsig(j,1)
             Fsm=Fsig(j,1); Fim=Fi(j,1);
             Fsig(j,1)=Fsig(i,1); Fi(j,1)=Fi(i,1);
             Fsig(i,1)=Fsm; Fi(i,1)=Fim;
           end
      end
      i=i+1;
end
