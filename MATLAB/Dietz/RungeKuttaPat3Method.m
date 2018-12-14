function [TOut,YOut] = RungeKuttaPat3Method(tBegin,tEnd,y,dt,RK)
global Reak Species
t=tBegin;
YOut=zeros(size(y,1),(tEnd-tBegin)/dt);
TOut=zeros(1,(tEnd-tBegin)/dt);
TOut(1)=t;
YOut(:,1)=y;

Y=zeros(size(y,1),RK.stage);
MY(RK.stage).M=sparse(size(y,1),size(y,1));
% Sigma=zeros(size(y,1),1);
p=3*RK.a(2,1)*(RK.a(3,1)+RK.a(3,2))*RK.b(3);

for i=1:(tEnd-tBegin)/dt
    Y(:,1)=y;%Y1
    MY(1).M=MHatMatrix(Reak,Species,t,Y(:,1));
    M=sparse(size(y,1),size(y,1));
    for j=2:RK.stage-1%Y2
        for k=1:j-1
            M=M+RK.a(j,k)*MY(k).M;
        end
        for k=1:size(y,1)
            M(:,k)=M(:,k)/max(Y(k,1),eps);
        end
        Y(:,j)=(eye(size(y,1))-dt*M)\y;
    end
    MY(2).M=MHatMatrix(Reak,Species,t,Y(:,2));
    M=sparse(size(y,1),size(y,1));
    for j=3:RK.stage%Y3
        for k=1:j-1
            M=M+RK.a(j,k)*MY(k).M;
        end
        for k=1:size(y,1)
            rho=y(k)*(Y(k,2)/y(k))^(1/p);%zu Y(3)
            M(:,k)=M(:,k)/max(rho,eps);
        end
        Y(:,j)=(eye(size(y,1))-dt*M)\y;
    end
    MY(3).M=MHatMatrix(Reak,Species,t,Y(:,3));
    M=sparse(size(y,1),size(y,1));
    for j=1:RK.stage-1%SIGMA
        M=M+RK.beta(j)*MY(j).M;
    end
    for j=1:size(y,1)
        my=y(j)*(Y(j,2)/y(j))^(1/RK.a(2,1));
        M(:,j)=M(:,j)/max(my,eps);
    end
    Sigma=(eye(size(y,1))-dt*M)\y;
    
    M=sparse(size(y,1),size(y,1));
    for j=1:RK.stage%Y4
        M=M+RK.b(j)*MY(j).M;
    end
    for j=1:size(y,1)
        M(:,j)=M(:,j)/max(Sigma(j,1),eps);
    end
    y=(eye(size(y,1))-dt*M)\y;
    t=t+dt;
    TOut(i+1)=t;
    YOut(:,i+1)=y;
end
TOut=TOut';
YOut=YOut';
end

