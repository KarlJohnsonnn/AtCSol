function [RK] = RungeKutta(name)

switch name
    case 'RK1'
        RK.stage=1;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        RK.b(1)=1;
    case 'RK2a'
        RK.stage=2;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        RK.a(2,1)=0.5;
        RK.b(2)=1;
    case 'RK2b'
        RK.stage=2;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        RK.a(2,1)=1;
        RK.b(1)=0.5;
        RK.b(2)=0.5;
    case 'SSP1'
        RK.stage=3;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        gamma=.7;
        RK.a(2,1)=1;
        RK.a(3,1)=gamma;
        RK.a(3,2)=gamma;
        RK.b(1)=(5*gamma-1)/(6*gamma);
        RK.b(2)=1./6.;
        RK.b(3)=1./(6*gamma);
        RK.beta(2)=1/(2*RK.a(2,1));
        RK.beta(1)=1-RK.beta(2);
    case 'meister1'
        RK.stage=3;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        RK.a(2,1)=1;%a21
        RK.a(3,1)=1./4.;%a31
        RK.a(3,2)=1./4.;%a32
        RK.b(1)=1./6.;%b1
        RK.b(2)=1./6.;%b2
        RK.b(3)=2./3.;%b3
        RK.beta(2)=1/(2*RK.a(2,1));
        RK.beta(1)=1-RK.beta(2);
    case 'meister2'
        RK.stage=3;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        RK.a(2,1)=1./2.;%a21
        RK.a(3,1)=0;%a31
        RK.a(3,2)=3./4.;%a32
        RK.b(1)=2./9.;%b1
        RK.b(2)=1./3.;%b2
        RK.b(3)=4./9.;%b3
        RK.beta(2)=1/(2*RK.a(2,1));
        RK.beta(1)=1-RK.beta(2);
    case 'meister3'
        RK.stage=3;
        RK.a=zeros(RK.stage,RK.stage);
        RK.b=zeros(1,RK.stage);
        gamma=1./2.;
        RK.a(2,1)=2./3.;
        RK.a(3,1)=1./6.;
        RK.a(3,2)=1./2.;
        RK.b(1)=1./4.;
        RK.b(2)=1./4.;
        RK.b(3)=1./2.;
        RK.beta(2)=1/(2*RK.a(2,1));
        RK.beta(1)=1-RK.beta(2);
end
end

