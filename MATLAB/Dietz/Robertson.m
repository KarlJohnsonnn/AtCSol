function [Reak,Species,tBegin,tEnd,dt]=Robertson()
%Aus Robertson.run
tBegin      = 1e-6 ;         % Startzeit Integration   (in [sec])
tBegin      = 0;         % Startzeit Integration   (in [sec])

%  tBegin      = 43200.0d0;          % Startzeit Integration   (in [sec])
%tBegin      = 46800.0d0          % Startzeit Integration   (in [sec])

tEnd        = 1e+10;     % Endzeit   Integration   (in [sec])
idate       = 010621;       %#ok<*NASGU> % ---------''---------    (Datum: YYMMDD)
rlat        = 45.0d0;       % ---------''---------    (Breitengrad)
rlon        = 0.0d0;        % ---------''---------    (Längengrad)
Temperature0= 280.0d0;
dt=100;

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
Reak(1).Type='CONST';
Reak(1).con(1)=4e-2;
%2.Reaktion
% B + B = B + C
Reak(2).Left=2;
Reak(2).NameL(1).Name='B';
Reak(2).NameL(2).Name='B';
Reak(2).KoeffL(1)=1;
Reak(2).KoeffL(2)=1;
Reak(2).Right=2;
Reak(2).NameR(1).Name='B';
Reak(2).KoeffR(1)=1;
Reak(2).NameR(2).Name='C';
Reak(2).KoeffR(2)=1;
Reak(2).Type='CONST';
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
Reak(3).Type='CONST';
Reak(3).con(1)=1e+4;

Species(3,1).c=9999;

Species(1).Name='A';
Species(1).Pos=1;
Species(1).MolMass=32;
Species(1).c= 1.0e0;
Species(2).Name='B';
Species(2).Pos=2;
Species(2).MolMass=32;
Species(2).c=0+1.e-8;  %+epsilon?
Species(3).Name='C';
Species(3).Pos=3;
Species(3).MolMass=32;
Species(3).c=0+1.e-8; %+epsilon?
end