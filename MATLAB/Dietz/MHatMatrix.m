function M=MHatMatrix(Reak,Species,tBegin,y)
for i=1:size(Reak,1)
    Reak(i).M=0;
    for j=1:Reak(i).Left
        Pos=GetPosList(Reak(i).NameL(j).Name,Species);
        Reak(i).M=Reak(i).M+Species(Pos).MolMass*Reak(i).KoeffL(j);
    end
end
M=zeros(size(Species,1),size(Species,1));
for i=1:size(Reak,1)
  M=InsertReak(y,Reak(i),M,Species,tBegin);
end
end

