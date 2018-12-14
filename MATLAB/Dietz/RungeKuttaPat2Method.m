function [TOut,YOut] = RungeKuttaPat2Method(tBegin,tEnd,y,dt,RK)
global Reak Species
t=tBegin;
YOut=zeros(size(y,1),(tEnd-tBegin)/dt);
TOut=zeros(1,(tEnd-tBegin)/dt);
TOut(1)=t;
YOut(:,1)=y;

Y=zeros(size(y,1),RK.stage);
MY(RK.stage).M=sparse(size(y,1),size(y,1));

for i=1:(tEnd-tBegin)/dt
  Y(:,1)=y;
  for j=2:RK.stage
    MY(j-1).M=MHatMatrix(Reak,Species,t,Y(:,j-1));
    M=sparse(size(y,1),size(y,1));
    for k=1:j-1
      M=M+RK.a(j,k)*MY(k).M;
    end
    for k=1:size(y,1)
      M(:,k)=M(:,k)/max(Y(k,1),eps);
    end
    Y(:,j)=(eye(size(y,1))-dt*M)\y;
  end
  MY(RK.stage).M=MHatMatrix(Reak,Species,t,Y(:,RK.stage));
  M=sparse(size(y,1),size(y,1));
  for j=1:RK.stage
    M=M+RK.b(j)*MY(j).M;
  end
  for j=1:size(y,1)
      sigma=y(j)*(Y(j,2)/y(j))^(1/RK.a(2,1));
      M(:,j)=M(:,j)/max(sigma,eps);
   end
  y=(eye(size(y,1))-dt*M)\y;
  t=t+dt;
  TOut(i+1)=t;
  YOut(:,i+1)=y;
end
TOut=TOut';
YOut=YOut';
end

