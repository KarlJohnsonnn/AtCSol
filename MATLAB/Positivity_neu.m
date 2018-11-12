function Positivity_neu()
%02.11.18 Dietz
%kurzes Programm bezüglich positivity,  um die Ausgabe der Matrizen von AtcSol zu vergleichen
close all;
clear all; %#ok<*CLALL>
 

%load reactions,species,start and end time
  [Reak,Species,tBegin,tEnd]=SmallStratoKPP();  %SmallStratoKPP

%    [Reak,Species,tBegin,tEnd]=Robertson()


for i=1:size(Reak,1)
    Reak(i).M=0;
    for j=1:Reak(i).Left
        Pos=GetPosList(Reak(i).NameL(j).Name,Species);
        Reak(i).M=Reak(i).M+Species(Pos).MolMass*Reak(i).KoeffL(j)
    end
end

y=zeros(size(Species));
for i=1:size(Species,1)
    y(Species(i).Pos)=Species(i).c;
end

M=zeros(size(Species,1),size(Species,1));
for i=1:size(Reak,1)
    M=InsertReak(y,Reak(i),M,Species,tBegin) %#ok<*NOPRT>
end


global Mfinal;
% Mfinal=M*inv(diag(y))
 Mfinal=M/diag(y)

% Mfinal=M


 tEnd=tBegin+5000;
%   tEnd=tBegin+300;

%ODE Solver
tspan=[tBegin tEnd]
options = odeset('RelTol',1e-8,'AbsTol',1.e-8);
% options = odeset('RelTol',1e-4,'AbsTol',[1e-6 1e-10 1e-6]);
%[tt,yy] = ode45(@odefun,tspan,y,options);%geht nicht
%[tt,yy] = ode15s(@odefun,tspan,y,options);

[tt,yy] = ode23s(@odefun,tspan,y,options);



% yy(:,2) = 1e9*yy(:,2);% robertson


fun_plot(tt,yy,Species,y)

end

%%

function dydt = odefun(~,y)
% dydt = zeros(size(y));
global Mfinal;

dydt = Mfinal*y;
end

function fun_plot(tt,yy,Species,y)
%Plot function
for i=1:size(y)
    figure
      plot(tt,yy(:,i),'-o')
    legend({Species(i).Name},'FontSize',22);  
    title( 'Concentrations' ,'FontSize',28);
end
 end

% 
% function fun_plot(tt,yy,Species,y)
% %Plot function
% for i=1:size(y)
%      hold on
%      plot(tt,yy(:,i),'-o')
%     tmplegend(i)={Species(i).Name}; %#ok<AGROW>
% end
% legend(tmplegend,'FontSize',22);
% title( 'Concentrations' ,'FontSize',28);
% end
 %%


function M=InsertReak(y,Reak,MM,Species,tBegin)
M=MM;
% Computation of the ReactionRate
r=Rate(Reak,tBegin);
for i=1:Reak.Left
    Pos=GetPos(Reak.NameL(i).Name,Species);
    r=r*y(Pos);
end

%Insert Diagonal
for i=1:Reak.Left
    PosL=GetPos(Reak.NameL(i).Name,Species);
    M(PosL,PosL)=M(PosL,PosL)-r*Reak.KoeffL(i);
    for j=1:Reak.Right
        PosR=GetPos(Reak.NameR(j).Name,Species);
        M(PosR,PosL)=M(PosR,PosL)+Reak.KoeffR(j)*r*Species(PosL).MolMass/Reak.M;
    end
end
end

function Pos=GetPos(name,Species)
Pos=0;
for i=1:size(Species,1)
    if strcmp(name,Species(i).Name)
        Pos=Species(i).Pos;
    end
end
end

function Pos=GetPosList(name,Species)
Pos=0;
for i=1:size(Species,1)
    if strcmp(name,Species(i).Name)
        Pos=i;
    end
end
end


function f=Rate(Reak,tBegin)
switch Reak.Type
    
    case 'Photo'
        suntime=mod(tBegin/3600,24);% REAL(dp), PARAMETER :: SunRise=4.50_dp, SunSet=19.50_dp
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'Photo2'
        suntime=mod(tBegin/3600,24);
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'Photo3'
        suntime=mod(tBegin/3600,24);
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun*sun*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'Photo1'
        disp('type unknown')
    case 'Const'
        f=Reak.con(1);
end
end

function sun=updatesun(Tlocal)
%calculation sun for SmallStratoKPP
SunRise=4.5;
SunSet=19.5;
Ttmp = (2*Tlocal-SunRise-SunSet) / (SunSet-SunRise);
if Ttmp>0
    Ttmp=Ttmp*Ttmp;
else
    Ttmp=-Ttmp*Ttmp;
end
sun = (1+cos(pi*Ttmp)) * 0.5;
end

function [Reak,Species,tBegin,tEnd]=SmallStratoKPP()
%Aus SmallStratoKPP.run
% tBegin      = 0.0d0          % Startzeit Integration   (in [sec])
  tBegin      = 43200.0d0;          % Startzeit Integration   (in [sec])
%tBegin      = 46800.0d0          % Startzeit Integration   (in [sec])

tEnd        = 172800.0d0;     % Endzeit   Integration   (in [sec])
idate       = 010621;       % ---------''---------    (Datum: YYMMDD)
rlat        = 45.0d0;       % ---------''---------    (Breitengrad)
rlon        = 0.0d0;        % ---------''---------    (Längengrad)
Temperature0= 280.0d0;
% mAir=  2.2063583412639662E+019 % mAir = (N2+O2) [molec/cm3] * 298.15 [K] / Temperature [1/K] * 850 [hPa] / 1013.25 [hPa]
% O2 = 2.0 O   
% O + O2 = O3    
% O3 = O + O2           
% O + O3 = 2.0 O2    
% O3 = O1D + O2            
% O1D + M  = O + M  
% O1D + O3 = 2.0 O2      
% NO + O3 = NO2 + O2   
% NO2 + O = NO + O2    
% NO2 = NO + O               


% CLASS: GAS
% O2 = 2.0 O
% PHOTO3: A: 2.643e-10
% 
% CLASS: GAS
% O + O2 = O3
% CONST: A: 8.018e-17
% 
% CLASS: GAS
% O3 = O + O2
% PHOTO: A: 6.120E-04
% 
% #CLASS: GAS
% #O + O3 = 2.0 O2
% #CONST: A: 1.576E-15
% 
% #CLASS: GAS
% #O3 = O1D + O2
% #PHOTO2: A: 1.070E-03
% 
% #CLASS: GAS
% #O1D = O
% #CONST: A: 8.1200e+16
% #O1D + M  = O + M
% #CONST: A: 7.110E-11
% 
% #CLASS: GAS
% #O1D + O3 = 2.0 O2
% #CONST: A: 1.200E-10
% 
% #CLASS: GAS
% #NO + O3 = NO2 + O2
% #CONST: A: 6.062E-15
% 
% #CLASS: GAS
% #NO2 + O = NO + O2
% #CONST: A: 1.069E-11
% 
% #CLASS: GAS
% #NO2 = NO + O
% #PHOTO: A: 1.289E-02

%1.Reaktion
% O2 = 2.0 O                 
Reak(2,1).con(1)=99999;
Reak(1).Left=1;
Reak(1).NameL(1).Name='O2';
Reak(1).KoeffL(1)=1;
Reak(1).Right=1;
Reak(1).NameR(1).Name='O';
Reak(1).KoeffR(1)=2;
Reak(1).Type='Photo3'; 
Reak(1).con(1)=2.643e-10;
%2.Reaktion
% O + O2 = O3       
Reak(2).Left=2;
Reak(2).NameL(1).Name='O';
Reak(2).KoeffL(1)=1;
Reak(2).NameL(2).Name='O2';
Reak(2).KoeffL(2)=1;
Reak(2).Right=1;
Reak(2).NameR(1).Name='O3';
Reak(2).KoeffR(1)=1;
Reak(2).Type='Const';
Reak(2).con(1)=8.018e-17;
%3.Reaktion
% O3 = O + O2                
Reak(3).Left=1;
Reak(3).NameL(1).Name='O3';
Reak(3).KoeffL(1)=1;
Reak(3).Right=2;
Reak(3).NameR(1).Name='O';
Reak(3).KoeffR(1)=1;
Reak(3).NameR(2).Name='O2';
Reak(3).KoeffR(2)=1;
Reak(3).Type='Photo';
Reak(3).con(1)=6.120E-04;
%4.Reaktion
% O + O3 = 2.0 O2       
Reak(4).Left=2;
Reak(4).NameL(1).Name='O';
Reak(4).KoeffL(1)=1;
Reak(4).NameL(2).Name='O3';
Reak(4).KoeffL(2)=1;
Reak(4).Right=1;
Reak(4).NameR(1).Name='O2';
Reak(4).KoeffR(1)=2;
Reak(4).Type='Const';
Reak(4).con(1)=1.576E-15;
%5.Reaktion
% O3 = O1D + O2             
Reak(5).Left=1;
Reak(5).NameL(1).Name='O3';
Reak(5).KoeffL(1)=1;
Reak(5).Right=2;
Reak(5).NameR(1).Name='O1D';
Reak(5).KoeffR(1)=1;
Reak(5).NameR(2).Name='O2';
Reak(5).KoeffR(2)=1;
Reak(5).Type='Photo2';
Reak(5).con(1)=1.070E-03;
%6.Reaktion
% #O1D = O
Reak(6).Left=1;
Reak(6).NameL(1).Name='O1D';
Reak(6).KoeffL(1)=1;
Reak(6).Right=1;
Reak(6).NameR(1).Name='O';
Reak(6).KoeffR(1)=1;
Reak(6).Type='Const';
% Reak(6).con(1)=8.1200e+16;
Reak(6).con(1)=7.110E-11;
%7.Reaktion
% O1D + O3 = 2.0 O2      
Reak(7).Left=2;
Reak(7).NameL(1).Name='O1D';
Reak(7).KoeffL(1)=1;
Reak(7).NameL(2).Name='O3';
Reak(7).KoeffL(2)=1;
Reak(7).Right=1;
Reak(7).NameR(1).Name='O2';
Reak(7).KoeffR(1)=2;
Reak(7).Type='Const';
Reak(7).con(1)=1.200E-10;
%8.Reaktion
% NO + O3 = NO2 + O2   
Reak(8).Left=2;
Reak(8).NameL(1).Name='NO';
Reak(8).KoeffL(1)=1;
Reak(8).NameL(2).Name='O3';
Reak(8).KoeffL(2)=1;
Reak(8).Right=2;
Reak(8).NameR(1).Name='NO2';
Reak(8).KoeffR(1)=1;
Reak(8).NameR(2).Name='O2';
Reak(8).KoeffR(2)=1;
Reak(8).Type='Const';
Reak(8).con(1)=6.062E-15;
%9.Reaktion
% NO2 + O = NO + O2    
Reak(9).Left=2;
Reak(9).NameL(1).Name='NO2';
Reak(9).KoeffL(1)=1;
Reak(9).NameL(2).Name='O';
Reak(9).KoeffL(2)=1;
Reak(9).Right=2;
Reak(9).NameR(1).Name='NO';
Reak(9).KoeffR(1)=1;
Reak(9).NameR(2).Name='O2';
Reak(9).KoeffR(2)=1;
Reak(9).Type='Const';
Reak(9).con(1)=1.069E-11;
%10.Reaktion
% NO2 = NO + O              
Reak(10).Left=1;
Reak(10).NameL(1).Name='NO2';
Reak(10).KoeffL(1)=1;
Reak(10).Right=2;
Reak(10).NameR(1).Name='NO';
Reak(10).KoeffR(1)=1;
Reak(10).NameR(2).Name='O';
Reak(10).KoeffR(2)=1;
Reak(10).Type='Photo';
Reak(10).con(1)=1.289E-02;


Species(3,1).c=9999;

% Species(1).Name='O2';
% Species(1).Pos=1;
% Species(1).MolMass=32;
% Species(1).c=1.697E+16;
% Species(2).Name='O';
% Species(2).Pos=2;
% Species(2).MolMass=16;
% Species(2).c=6.624E+08  ;
% Species(3).Name='O3';
% Species(3).Pos=3;
% Species(3).MolMass=48;
% Species(3).c=5.326E+11;
% Species(4).Name='O1D';
% Species(4).Pos=4;
% Species(4).MolMass=16;
% Species(4).c=9.906E+01 ;
% Species(5).Name='NO';
% Species(5).Pos=5;
% Species(5).MolMass=30;
% Species(5).c=8.725E+08;
% Species(6).Name='NO2';
% Species(6).Pos=6;
% Species(6).MolMass=46;
% Species(6).c=2.240E+08;

%test
Species(1).Name='O2';
Species(1).Pos=5;
Species(1).MolMass=32;
Species(1).c=1.697E+16;
Species(2).Name='O';
Species(2).Pos=3;
Species(2).MolMass=16;
Species(2).c=6.624E+08  ;
Species(3).Name='O3';
Species(3).Pos=6;
Species(3).MolMass=48;
Species(3).c=5.326E+11;
Species(4).Name='O1D';
Species(4).Pos=4;
Species(4).MolMass=16;
Species(4).c=9.906E+01 ;
Species(5).Name='NO';
Species(5).Pos=1;
Species(5).MolMass=30;
Species(5).c=8.725E+08;
Species(6).Name='NO2';
Species(6).Pos=2;
Species(6).MolMass=46;
Species(6).c=2.240E+08;

end


function [Reak,Species,tBegin,tEnd]=Robertson() %#ok<*DEFNU>
%Aus Robertson.run
 tBegin      = 0.0d0 ;         % Startzeit Integration   (in [sec])
%  tBegin      = 43200.0d0;          % Startzeit Integration   (in [sec])
%tBegin      = 46800.0d0          % Startzeit Integration   (in [sec])

tEnd        = 172800.0d0;     % Endzeit   Integration   (in [sec])
idate       = 010621;       %#ok<*NASGU> % ---------''---------    (Datum: YYMMDD)
rlat        = 45.0d0;       % ---------''---------    (Breitengrad)
rlon        = 0.0d0;        % ---------''---------    (Längengrad)
Temperature0= 280.0d0;
% mAir=  2.2063583412639662E+019 % mAir = (N2+O2) [molec/cm3] * 298.15 [K] / Temperature [1/K] * 850 [hPa] / 1013.25 [hPa]
 
%    BEGIN_INITIAL                #  Initial Concentrations [molec/cm^3]
% A      1.0e0
%    END_INITIAL

% CLASS: GAS  
% A = B
% CONST:  A: 4e-2
% 
% CLASS: GAS
% 2.0 B = B + C 
% CONST:  A: 3e+7
% 
% CLASS: GAS
% B + C = A + C
% CONST:  A: 1e+4

%1.Reaktion
% A = B                
Reak(2,1).con(1)=99999;
Reak(1).Left=1;
Reak(1).NameL(1).Name='A';
Reak(1).KoeffL(1)=1;
Reak(1).Right=1;
Reak(1).NameR(1).Name='B';
Reak(1).KoeffR(1)=1;
Reak(1).Type='Const'; 
Reak(1).con(1)=4e-2;
%2.Reaktion
% 2.0 B = B + C       
Reak(2).Left=1;
Reak(2).NameL(1).Name='B';
Reak(2).KoeffL(1)=2;  
Reak(2).Right=2;
Reak(2).NameR(1).Name='B';
Reak(2).KoeffR(1)=1;
Reak(2).NameR(2).Name='C';
Reak(2).KoeffR(2)=1;
Reak(2).Type='Const';
Reak(2).con(1)=3e+7;
%3.Reaktion
% B + C = A + C               
Reak(3).Left=2;
Reak(3).NameL(1).Name='B';
Reak(3).NameL(2).Name='C';
Reak(3).KoeffL(1)=1;
Reak(3).KoeffL(2)=1;
Reak(3).Right=2;
Reak(3).NameR(1).Name='A';
Reak(3).KoeffR(1)=1;
Reak(3).NameR(2).Name='C';
Reak(3).KoeffR(2)=1;
Reak(3).Type='Const';
Reak(3).con(1)=1e+4;
 


Species(3,1).c=9999;

Species(1).Name='A';
Species(1).Pos=1;
Species(1).MolMass=32;
Species(1).c= 1.0e0;
Species(2).Name='B';
Species(2).Pos=2;
Species(2).MolMass=32;
Species(2).c=0+eps  ;%+epsilon?
Species(3).Name='C';
Species(3).Pos=3;
Species(3).MolMass=32;
Species(3).c=0+eps ;%+epsilon?
 
end
