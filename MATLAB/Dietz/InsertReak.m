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
