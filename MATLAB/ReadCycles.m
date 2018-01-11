%% reading list of cycles
function cyclic_set=ReadCycles(path)

fileID = fopen(path,'r');
iC = 1;
tline = fgetl(fileID);
while ischar(tline)
    cyclic_set(iC).len = 0;
    cyclic_set(iC).Species = [];
    cyclic_set(iC).Reactions = [];
    tline = fgetl(fileID);
    iC = iC + 1;
end
frewind(fileID);

iC = 1;
tline = fgetl(fileID);
while ischar(tline)
    cyclic_set(iC).Species = str2num(tline);
    cyclic_set(iC).len  = length(cyclic_set(iC).Species);
    tline = fgetl(fileID);
    iC = iC + 1;
end
fclose(fileID);
end 