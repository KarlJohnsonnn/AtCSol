%% Funktion zum Auslesen der in Fortran90 produzierten Matrizen
%
%

% function [SpeciesNames] = ReadSpeciesNames(Path)
function [SpeciesNames] = ReadSpeciesNames(Path)

fileID = fopen(Path,'r');

% skip headlines
for i=1:24
    tline = fgetl(fileID);
    if (i == 16)
        ns = sscanf( tline , '%d ' );   
    end
    if (i == 21 )
        nkat = sscanf( tline , '%d ' );
    end
end

SpeciesNames = cell(ns,1);

% read matrix dimension
for i = 1:ns
    tline = fgetl(fileID);
    SpeciesNames( i , :) = {sscanf( tline(3:end-1) , '%s ' )};
end
end

