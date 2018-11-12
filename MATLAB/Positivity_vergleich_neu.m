function Positivity_vergleich_neu()
%29.10.18 Dietz
%kurzes Programm bezüglich positivity,  um die Ausgabe der Matrizen von AtcSol zu vergleichen
close all;
 clear all;
%Liste wird von oben abgearbeitet. Eswird stufenweise eine weitere Reaktion %aktiviert
% O2 = 2.0 O                0.0000000000000000
% O + O2 = O3       901297607.04000008
% O3 = O + O2               0.0000000000000000
% O + O3 = 2.0 O2      556003.72224000003
% O3 = O1D + O2             0.0000000000000000
% O1D + M  = O + M  571905079.20000005
% O1D + O3 = 2.0 O2      6331.1227200000003
% NO + O3 = NO2 + O2  2816971.9970000004
% NO2 + O = NO + O2   1586156.5440000000
% NO2 = NO + O              0.0000000000000000

Reak(2,1).con(1)=99999;

Reak(1).Left=1;
Reak(1).NameL(1).Name='O2';
Reak(1).Right=1;
Reak(1).NameR(1).Name='O';
Reak(1).KoeffR(1)=2;
Reak(1).Type='Photo';


Reak(2).Left=2;
Reak(2).NameL(1).Name='O';
Reak(2).NameL(2).Name='O2';
Reak(2).Right=1;
Reak(2).NameR(1).Name='O3';
Reak(2).KoeffR(1)=1;
Reak(2).Type='Const';
Reak(2).con(1)=901297607.04000008;

Species(3,1).c=9999;

Species(1).Name='O2';
Species(1).Pos=1;
Species(1).MolMass=32;
Species(1).c=5;
Species(2).Name='O';
Species(2).Pos=2;
Species(2).MolMass=16;
Species(2).c=30;
Species(3).Name='O3';
Species(3).Pos=3;
Species(3).MolMass=48;
Species(3).c=70;

for i=1:size(Reak,1)
    Reak(i).M=0;
    for j=1:Reak(i).Left
        Pos=GetPosList(Reak(i).NameL(j).Name,Species);
        Reak(i).M=Reak(i).M+Species(Pos).MolMass;
    end
end

y=zeros(size(Species));
for i=1:size(Species,1)
    y(Species(i).Pos)=Species(i).c;
end

M=zeros(3,3);
M=InsertReak(y,Reak(2),M,Species);


% %%
% %example: 10 Strato reactions
% n=6;
% m=10;
% r=[1,1,1,1,1,1,1,1,1,1;];
% Matrix=pos_berechnung(n,m,r)


%%
%example: 1 Strato reactions
n=2;
m=1;
r=[1,0,0,0,0,0,0,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 2 Strato reactions
n=3;
m=2;
r=[1,1,0,0,0,0,0,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 3 Strato reactions
n=3;
m=3;
r=[1,1,1,0,0,0,0,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 4 Strato reactions
n=3;
m=4;
r=[1,1,1,1,0,0,0,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 5 Strato reactions
n=4;
m=5;
r=[1,1,1,1,1,0,0,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 6 Strato reactions
n=4;
m=6;
r=[1,1,1,1,1,1,0,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 7 Strato reactions
n=4;
m=7;
r=[1,1,1,1,1,1,1,0,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 8 Strato reactions
n=6;
m=8;
r=[1,1,1,1,1,1,1,1,0,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 9 Strato reactions
n=6;
m=9;
r=[1,1,1,1,1,1,1,1,1,0;];
Matrix=pos_berechnung(n,m,r)
%%
%example: 10 Strato reactions
n=6;
m=10;
r=[1,1,1,1,1,1,1,1,1,1;];
Matrix=pos_berechnung(n,m,r)
%%
end
function M=InsertReak(y,Reak,MM,Species)
M=MM;
% Computation of the ReactionRate
r=Rate(Reak);
for i=1:Reak.Left
    Pos=GetPos(Reak.NameL(i).Name,Species);
    r=r*y(Pos);
end
    
%Insert Diagonal
for i=1:Reak.Left
    PosL=GetPos(Reak.NameL(i).Name,Species);
    M(PosL,PosL)=M(PosL,PosL)-r;
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

function f=Rate(Reak)
switch Reak.Type
    
    case 'Photo'
        f=0;
    case 'Const'
        f=Reak.con(1);
end
end

function Matrix=pos_berechnung(n,m,r)
% m...reactions
% n...Species 
% r...Vektor der aktiven reactions
  
%Reaktionsgeschwindigkeiten aus ATcSol
rVec2 =[  0.0000000000000000;        901297607.04000008;  0.0000000000000000 ;       556003.72224000003 ;       0.0000000000000000 ;        8.0436720000000000E+018 ;       6331.1227200000003 ;       2816971.9970000004 ;       1586156.5440000000 ;       0.0000000000000000;];
%Startkonzentrationen aus ATcSol
yVec2=[8.725E+08,2.240E+08,6.624E+08,9.906E+01 ,1.697E+16 ,5.326E+11;];
%Vektor der Molekulargewichte aus ATcSol
mi2=[30,46,16,16,32,48;];
        Vergleich=0;


switch m
    
    case 1
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])
        
        alpha= [0,1;]
        
        beta=  [2,0;]
        
        yVec(1,:)=yVec2(3);
        yVec(2,:)=yVec2(5);
        mi(1)=mi2(3);
        mi(2)=mi2(5);
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        
    case 2
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])
        
        alpha= [0,1,0;
            1,1,0;]
        
        beta=  [2,0,0;
            0,0,1;]
         
           yVec(1,:)=yVec2(3);
         mi(1)=mi2(3);
    for i=2:1:n
            yVec(i,:)=yVec2(i+3);
            mi(i)=mi2(i+3);
    end
        
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        
    case 3
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,1,0;
    1,1,0;
    0,0,1;]

beta=  [2,0,0;
    0,0,1;
    1,1,0;             ]
         
         yVec(1,:)=yVec2(3);
         mi(1)=mi2(3);
    for i=2:1:n
            yVec(i,:)=yVec2(i+3);
            mi(i)=mi2(i+3);
    end
        
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        
    case 4
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,1,0;
    1,1,0;
    0,0,1;
    1,0,1;                     ]

beta=  [2,0,0;
    0,0,1;
    1,1,0;
    0,2,0;             ]
         
         yVec(1,:)=yVec2(3);
         mi(1)=mi2(3);
    for i=2:1:n
            yVec(i,:)=yVec2(i+3);
            mi(i)=mi2(i+3);
    end
        
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
         %Diag M aus AtcSol:
        Test=diag([            -1.3614939776000001
            -5.3111232000000006E-008
            -1.0439424000000001E-006;]);
                Vergleich=1;
           
     case 5
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,0,1,0;
    1,0,1,0;
    0,0,0,1;
    1,0,0,1;
    0,0,0,1;                 ]

beta=  [2,0,0,0;
    0,0,0,1;
    1,0,1,0;
    0,0,2,0;
    0,1,1,0;             ]
    for i=1:1:n
            yVec(i,:)=yVec2(i+2);
            mi(i)=mi2(i+2);
    end
        
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        %Diag M aus AtcSol:
        Test=diag([            -1.3614939776000001
            0.0000000000000000
            -5.3111232000000006E-008
            -1.0439424000000001E-006;]);
                Vergleich=1;
            
    case 6
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,0,1,0;
    1,0,1,0;
    0,0,0,1;
    1,0,0,1;
    0,0,0,1;
    0,1,0,0;             ]

beta=  [2,0,0,0;
    0,0,0,1;
    1,0,1,0;
    0,0,2,0;
    0,1,1,0;
    1,0,0,0;             ]
    for i=1:1:n
            yVec(i,:)=yVec2(i+2);
            mi(i)=mi2(i+2);
    end
        
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
         %Diag M aus AtcSol:
        Test=diag([            -1.3614939776000001
            -81200000000000000.
            -5.3111232000000006E-008
            -1.0439424000000001E-006;]);
                Vergleich=1;
                 
        case 7
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,0,1,0;
    1,0,1,0;
    0,0,0,1;
    1,0,0,1;
    0,0,0,1;
    0,1,0,0;
    0,1,0,1;         ]

beta=  [2,0,0,0;
    0,0,0,1;
    1,0,1,0;
    0,0,2,0;
    0,1,1,0;
    1,0,0,0;
    0,0,2,0;         ]
    for i=1:1:n
            yVec(i,:)=yVec2(i+2);
            mi(i)=mi2(i+2);
    end
        
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        %Diag M aus AtcSol:
        Test=diag([            -1.3638885376000001
            -81200000000000064.
            -5.3111232000000006E-008
            -1.0558296000000001E-006;]);
                Vergleich=1;
                  
    case 8
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,0,0,0,1,0;
    0,0,1,0,1,0;
    0,0,0,0,0,1;
    0,0,1,0,0,1;
    0,0,0,0,0,1;
    0,0,0,1,0,0;
    0,0,0,1,0,1;
    1,0,0,0,0,1;        ]

beta=  [0,0,2,0,0,0;
    0,0,0,0,0,1;
    0,0,1,0,1,0;
    0,0,0,0,2,0;
    0,0,0,1,1,0;
    0,0,1,0,0,0;
    0,0,0,0,2,0;
    0,1,0,0,1,0;        ]
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        yVec=yVec2;
         mi=mi2;
          %Diag M aus AtcSol:
        Test=diag([-3.2286212000000006E-003
            0.0000000000000000
            -1.3638885376000001
            -81200000000000064.
            -5.3111232000000006E-008
            -6.3449246000000012E-006;]);
                Vergleich=1;
              
    case 9
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

        
alpha= [0,0,0,0,1,0;
    0,0,1,0,1,0;
    0,0,0,0,0,1;
    0,0,1,0,0,1;
    0,0,0,0,0,1;
    0,0,0,1,0,0;
    0,0,0,1,0,1;
    1,0,0,0,0,1;
    0,1,1,0,0,0;    ]

beta=  [0,0,2,0,0,0;
    0,0,0,0,0,1;
    0,0,1,0,1,0;
    0,0,0,0,2,0;
    0,0,0,1,1,0;
    0,0,1,0,0,0;
    0,0,0,0,2,0;
    0,1,0,0,1,0;
    1,0,0,0,1,0;    ]
        for i=1:1:m
            rVec(i,:)=rVec2(i);
        end
        yVec=yVec2;
         mi=mi2;
           %Diag M aus AtcSol:
        Test=diag([-3.2286212000000006E-003
            -7.0810559999999996E-003
            -1.3638885376000001
            -81200000000000064.
            -5.3111232000000006E-008
            -6.3449246000000012E-006;]);
                Vergleich=1;

    case 10
        disp(['example: ', num2str(m),' reactions and ', num2str(n),' Species'])

alpha= [0,0,0,0,1,0;
    0,0,1,0,1,0;
    0,0,0,0,0,1;
    0,0,1,0,0,1;
    0,0,0,0,0,1;
    0,0,0,1,0,0;
    0,0,0,1,0,1;
    1,0,0,0,0,1;
    0,1,1,0,0,0;
    0,1,0,0,0,0;]

beta=  [0,0,2,0,0,0;
    0,0,0,0,0,1;
    0,0,1,0,1,0;
    0,0,0,0,2,0;
    0,0,0,1,1,0;
    0,0,1,0,0,0;
    0,0,0,0,2,0;
    0,1,0,0,1,0;
    1,0,0,0,1,0;
    1,0,1,0,0,0;]

        yVec=yVec2;
        rVec=rVec2; 
        mi=mi2;
        %Diag M aus AtcSol: 
        Test=diag([-3.2286212000000006E-003
            -7.0810559999999996E-003
            -1.3638885376000001
            -81200000000000064.
            -5.3111232000000006E-008
            -6.3449246000000012E-006;]);
                       Vergleich=1;
     %Diag M aus AtcSol: startzeit 432000                     
    Test=diag([                      -3.2286212000000006E-003
      -1.9971056000000001E-002
     -1.3638885376000001
     -81200000000000064.
    -5.3375532000000004E-008
      -1.6883449245999999E-003])
   
 
    otherwise
        warning('Unexpected reaction type. No calculation possible.')
end



%Vektor der Summe der molekülmasse pro Reaktion
MWVect=alpha*mi';
%Diagonalform von mi:
Dmi=diag(mi);

Dalpha=zeros(n,n);
for i=1:1:m
    Dalpha=Dalpha+diag(alpha(i,:))*rVec(i);
end
%Diagonalform von 1/y:
Dyy=diag(yVec);
Dy=inv(Dyy);
%D(r/MWVec):
rVec=diag(rVec);
MWVect=diag(MWVect);

%Berechnung der M-Matrix:
Matrix1= (beta'*rVec/MWVect*alpha*Dmi)*Dy;
Matrix2= Dalpha*Dy;
Matrix= Matrix1-Matrix2;
if Vergleich == 1 
%optionaler Test
Vergleich_der_Hauptdiagonale_mit_AtcSol=diag(Matrix)-diag(Test)
end
end