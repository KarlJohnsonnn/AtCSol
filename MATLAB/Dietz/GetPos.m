function Pos=GetPos(name,Species)
Pos=0;
for i=1:size(Species,1)
    if strcmp(name,Species(i).Name)
        Pos=Species(i).Pos;
    end
end
end
