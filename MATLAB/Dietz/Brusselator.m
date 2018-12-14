function [Reak,Species,tBegin,tEnd,dt]=Brusselator() 
% Brusselator aus Unconditionally positive and conservative third order
% modified Patankar-Runge-Kutta discretizations of
% production-destruction systems
% Stefan Kopecza and Andreas Meistera
% Page 24
tBegin      = 0.0d0 ;         % Startzeit Integration   (in [sec])
tEnd        = 6.0d0;     % Endzeit   Integration   (in [sec])
dt=0.05;
%1.Reaktion
% A1 = A5
Reak(2,1).con(1)=99999;
Reak(1).Left=1;
Reak(1).NameL(1).Name='A1';
Reak(1).KoeffL(1)=1;
Reak(1).Right=1;
Reak(1).NameR(1).Name='A5';
Reak(1).KoeffR(1)=1;
Reak(1).Type='CONST';
Reak(1).con(1)=1;
%2.Reaktion
% A2 + A5 = A3 + A6
Reak(2).Left=2;
Reak(2).NameL(1).Name='A2';
Reak(2).NameL(2).Name='A5';
Reak(2).KoeffL(1)=1;
Reak(2).KoeffL(2)=1;
Reak(2).Right=2;
Reak(2).NameR(1).Name='A3';
Reak(2).KoeffR(1)=1;
Reak(2).NameR(2).Name='A6';
Reak(2).KoeffR(2)=1;
Reak(2).Type='CONST';
Reak(2).con(1)=1;
%3.Reaktion
% A6 + A5 + A5 = 3A5
Reak(3).Left=3;
Reak(3).NameL(1).Name='A6';
Reak(3).NameL(2).Name='A5';
Reak(3).NameL(3).Name='A5';
Reak(3).KoeffL(1)=1;
Reak(3).KoeffL(2)=1;
Reak(3).KoeffL(3)=1;
Reak(3).Right=1;
Reak(3).NameR(1).Name='A5';
Reak(3).KoeffR(1)=3;
Reak(3).Type='CONST';
Reak(3).con(1)=1;
%4.Reaktion
% A5 = A4
Reak(4).Left=1;
Reak(4).NameL(1).Name='A5';
Reak(4).KoeffL(1)=1;
Reak(4).Right=1;
Reak(4).NameR(1).Name='A4';
Reak(4).KoeffR(1)=1;
Reak(4).Type='CONST';
Reak(4).con(1)=1;

Species(3,1).c=9999;

Species(1).Name='A1';
Species(1).Pos=1;
Species(1).MolMass=1;
Species(1).c= 10;
Species(2).Name='A2';
Species(2).Pos=2;
Species(2).MolMass=1;
Species(2).c= 10;
Species(3).Name='A3';
Species(3).Pos=3;
Species(3).MolMass=1;
Species(3).c=eps; %+epsilon?
Species(4).Name='A4';
Species(4).Pos=4;
Species(4).MolMass=1;
Species(4).c=eps; %+epsilon?
Species(5).Name='A5';
Species(5).Pos=5;
Species(5).MolMass=1;
Species(5).c=0.1;
Species(6).Name='A6';
Species(6).Pos=6;
Species(6).MolMass=1;
Species(6).c=0.1;
end


