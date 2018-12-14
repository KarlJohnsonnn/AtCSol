

function ChemOutPut = ReadChemFile2(Path)
ChemOutPut.Mechanism   = ReadMechanismName(Path);
ChemOutPut.SpcNumbers  = GetSpeciesNumbers(Path);
ChemOutPut.SpcNames    = ReadSpeciesNames(Path,ChemOutPut.SpcNumbers.nspc);
ChemOutPut.ReacNumbers = GetReactionNumbers(Path);

ChemOutPut.ReacCoefficients= GetReacCoefficients(Path);%new
end


%% subroutines

function [num] = GetReacCoefficients(Path)

fileID = fopen(Path,'r');
num = cell(10,1);%hier reaktionsanzahl
while ~feof(fileID)
    
    tline = strtrim(fgetl(fileID));
    found = contains(tline,'======================  Reactions');
   
    if found
        strtrim(fgetl(fileID)); % skip blank line
       
        
        for i=1:10
            
            strtrim(fgetl(fileID)); % skip blank line 
       strtrim(fgetl(fileID)); % skip blank line
       strtrim(fgetl(fileID)); % skip blank line
       strtrim(fgetl(fileID)); % skip blank line
        strtrim(fgetl(fileID)); % skip blank line
       strtrim(fgetl(fileID)); % skip blank line
       
%         num.coeff = sscanf( tline(27:end) , '%s ');
        num.coeff2(i)  = sscanf(strtrim(fgetl(fileID)), '%d');
                num.coeff(i)  = sscanf(strtrim(fgetl(fileID)), '%*d %f');
        end
        break;
        
        
%           for i = 1:ns
%             tline = strtrim(fgetl(fileID));
%             name{i} = sscanf( tline(2:end-1) , '%s ');
%           end
%         
    end
end
end




% gather the species names
function name = ReadMechanismName(Path)

fileID = fopen(Path,'r');

while ~feof(fileID)
    
    tline = strtrim(fgetl(fileID));
    found = contains(tline,'Chemical Mechanism:');
    
    if found
        name = sscanf( tline(27:end) , '%s ');
        break;
    end
end
end

% gather the numbers of species
function [num] = GetSpeciesNumbers(Path)

fileID = fopen(Path,'r');

while ~feof(fileID)
    
    found = contains(strtrim(fgetl(fileID)),'Numbers');
    
    if found
        strtrim(fgetl(fileID)); % skip blank line
        num.nspc  = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nsgas  = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nsaqua = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nspart = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nssoli = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nspass = sscanf(strtrim(fgetl(fileID)), '%d');
        break;
    end
end
end

% gather the species names
function [name] = ReadSpeciesNames(Path,ns)

fileID = fopen(Path,'r');
name = cell(ns,1);

while ~feof(fileID)
    
    found = contains(strtrim(fgetl(fileID)),'Species Names');
    
    if found
        strtrim(fgetl(fileID)); % skip blank line% read matrix dimension
        for i = 1:ns
            tline = strtrim(fgetl(fileID));
            name{i} = sscanf( tline(2:end-1) , '%s ');
        end
        break;
    end
end
end

% gather the numbers of species
function [num] = GetReactionNumbers(Path)

fileID = fopen(Path,'r');

while ~feof(fileID)
    
    found = contains(strtrim(fgetl(fileID)),'Description of Reactions');
    
    if found
        a=strtrim(fgetl(fileID)); % skip blank line
        num.nreac = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % gaseous reaction numbers
        num.nrgas = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_photo   = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_const   = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_temp    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_simp    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_lind    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_troe    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_spec    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrgas_special = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % henry reaction numbers
        num.nrhenry = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % diss reaction numbers
        num.nrdiss = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrdiss_dconst  = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrdiss_dtemp   = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrdiss_special = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % aqueous reaction numbers
        num.nraqua = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nraqua_photo   = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nraqua_const   = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nraqua_temp    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nraqua_spec    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nraqua_special = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % parti reaction numbers
        num.nrparti = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrparti_special = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % solid reaction numbers
        num.nrsolid = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrsolid_dtemp3  = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrsolid_equi    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrsolid_spec    = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrsolid_special = sscanf(strtrim(fgetl(fileID)), '%d');
        
        % mirco reaction numbers
        num.nrmicro = sscanf(strtrim(fgetl(fileID)), '%d');
        num.nrmicro_special = sscanf(strtrim(fgetl(fileID)), '%d');
        break;
    end
end
end
