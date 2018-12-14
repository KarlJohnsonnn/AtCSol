function [Reak,Species,tBegin,tEnd,dt]=SmallStratoKPP()
%Aus SmallStratoKPP.run
%Für Vergleich mit Kpp plots:
%tBegin=12*3600;%43200
% tEnd=84*3600;%302400;
dt=1;

tBegin      = 43200.0d0;          % Startzeit Integration   (in [sec])
tEnd        = 46800.0d0;     % Endzeit   Integration   (in [sec])

% tEnd        = 302400.0d0;     % Endzeit   Integration   (in [sec])
idate       = 010621;       % ---------''---------    (Datum: YYMMDD)
rlat        = 45.0d0;       % ---------''---------    (Breitengrad)
rlon        = 0.0d0;        % ---------''---------    (Längengrad)
Temperature0= 280.0d0; %#ok<*NASGU>
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
Reak(10,1).con(1)=99999;
Reak(1).Left=1;
Reak(1).NameL(1).Name='O2';
Reak(1).KoeffL(1)=1;
Reak(1).Right=1;
Reak(1).NameR(1).Name='O';
Reak(1).KoeffR(1)=2;
Reak(1).Type='PHOTO3';
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
Reak(2).Type='CONST';
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
Reak(3).Type='PHOTO';
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
Reak(4).Type='CONST';
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
Reak(5).Type='PHOTO2';
Reak(5).con(1)=1.070E-03;
%6.Reaktion
% #O1D + M = O + M
Reak(6).Left=1;
Reak(6).NameL(1).Name='O1D';
Reak(6).KoeffL(1)=1;
Reak(6).Right=1;
Reak(6).NameR(1).Name='O';
Reak(6).KoeffR(1)=1;
Reak(6).Type='CONST';
% Reak(6).con(1)=7.110E-11;%alt
Reak(6).con(1)=7.110E-11*8.120E+16;
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
Reak(7).Type='CONST';
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
Reak(8).Type='CONST';
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
Reak(9).Type='CONST';
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
Reak(10).Type='PHOTO';
Reak(10).con(1)=1.289E-02;

Species(3,1).c=9999;

% %alte Reihenfolge
Species(1).Name='O2';
Species(1).Pos=1;
Species(1).MolMass=32;
Species(1).c=1.697E+16;
Species(2).Name='O';
Species(2).Pos=2;
Species(2).MolMass=16;
Species(2).c=6.624E+08  ;
Species(3).Name='O3';
Species(3).Pos=3;
Species(3).MolMass=48;
Species(3).c=5.326E+11;
Species(4).Name='O1D';
Species(4).Pos=4;
Species(4).MolMass=16;
Species(4).c=9.906E+01 ;
Species(5).Name='NO';
Species(5).Pos=5;
Species(5).MolMass=30;
Species(5).c=8.725E+08;
Species(6).Name='NO2';
Species(6).Pos=6;
Species(6).MolMass=46;
Species(6).c=2.240E+08;

%AtCSol Reihenfolge:
% Species(1).Name='O2';
% Species(1).Pos=5;
% Species(1).MolMass=32;
% Species(1).c=1.697E+16;
% Species(2).Name='O';
% Species(2).Pos=3;
% Species(2).MolMass=16;
% Species(2).c=6.624E+08  ;
% Species(3).Name='O3';
% Species(3).Pos=6;
% Species(3).MolMass=48;
% Species(3).c=5.326E+11;
% Species(4).Name='O1D';
% Species(4).Pos=4;
% Species(4).MolMass=16;
% Species(4).c=9.906E+01 ;
% Species(5).Name='NO';
% Species(5).Pos=1;
% Species(5).MolMass=30;
% Species(5).c=8.725E+08;
% Species(6).Name='NO2';
% Species(6).Pos=2;
% Species(6).MolMass=46;
% Species(6).c=2.240E+08;

end