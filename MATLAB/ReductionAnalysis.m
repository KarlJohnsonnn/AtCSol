


%% reduction analysis
close all
clear all

global b1 b2 len_f len_r mech

%% --------------> DEFINE PATHS <--------------

% declare path to AtCSol folder
AtCSol_path = '/Users/schimmel/Code/AtCSol/';

% declare mechanism (.chem, BSP, ini)
% chem_file_f = 'MCM32_CAPRAM40_full.chem';
% chem_file_r = 'MCM32+CAPRAM40_red.chem';
% mech     = 'MCM32+C40';             % name BSP in .run file
% ini_file = 'Urban2.ini';            % name of initial data set in .run file


% RACM + CAPRAM2.4 
chem_file_f = 'RACM+C24.chem';
chem_file_r = 'RACM+C24_red.chem';
mech     = 'RACM+C24';
ini_file = 'Urban.ini';





%% BEGIN MAIN

chem_path = [AtCSol_path,'CHEM/'];
ncdf_path = [AtCSol_path,'NetCDF/'];
fig_path  = [AtCSol_path,'OUTPUT/'];
ini_path  = [AtCSol_path,'INI/'];


% NetCDF files
m_full = [ncdf_path, mech, '_full.nc'];
m_red  = [ncdf_path, mech, '_red.nc'];

% Output file
outID = fopen([AtCSol_path,'OUTPUT/Statistics_',mech,'.txt'],'w');

% plot intervall
b1 = 10.0;   % 12 noon
b2 = 37.0;   % 12 noon next day


% read chem-file (numbers and species names)
Full_Mech    = ReadChemFile([chem_path,chem_file_f]);
Reduced_Mech = ReadChemFile([chem_path,chem_file_r]);

% print head
fprintf(outID,' \n');fprintf(outID,' \n');fprintf(outID,' \n');
fprintf(outID,'   MECHANISM ::    FULL ::  %-40s\n',Full_Mech.Mechanism);
fprintf(outID,'             :: REDUCED ::  %-40s\n',Reduced_Mech.Mechanism);

% analyse Numbers
Analyse_SpeciesNumers(outID,Full_Mech,Reduced_Mech);
Analyse_ReactionNumers(outID,Full_Mech,Reduced_Mech);

% gather diagnose species
Diag_Spc = GetDiagnoseSpecies([ini_path,ini_file],Reduced_Mech.SpcNames);

% gather time arrays
[t_f, t_r] = GetNetCDF_Data(m_full,m_red,'time');

% cut of some of the values
idx_f = cutOff(t_f,b1,b2);
idx_r = cutOff(t_r,b1,b2);
t_f = t_f(idx_f(1):idx_f(2));
t_r = t_r(idx_f(1):idx_f(2));
len_f = length(t_f);
len_r = length(t_r);

% gater species concentrations
n_diag = length(Diag_Spc);
c_f = cell(n_diag,1);
c_r = cell(n_diag,1);
dev_max  = cell(n_diag,1);
dev_mean = cell(n_diag,1);

% Plot
for iSpc = 1:n_diag
    [c_f{iSpc}, c_r{iSpc}] = GetNetCDF_Data(m_full,m_red,Diag_Spc{iSpc}.ncdf);
    
    c_f{iSpc} = c_f{iSpc}(idx_f(1):idx_f(2));
    c_r{iSpc} = c_r{iSpc}(idx_f(1):idx_f(2));
    
    % calculate deviations
    dev_max{iSpc} = max_deviation(c_f{iSpc},c_r{iSpc});
    dev_mean{iSpc} = mean_deviation(c_f{iSpc},c_r{iSpc});
    
%    if ( strcmp(Diag_Spc{iSpc}.name,'aHCL') || strcmp(Diag_Spc{iSpc}.name,'HCL') )
%     Plot_Concentrations(t_f,t_r,c_f{iSpc},c_r{iSpc},dev{iSpc},Diag_Spc{iSpc})
%    end
end

Analyse_SpeciesDeviation(outID,c_f,c_r,dev_max,Diag_Spc)

fclose(outID);

%% END MAIN




%% BEGIN SUBROUTINES

% Plot the statistic for species numbers
function n = Analyse_SpeciesNumers(id,full,red)

n_f = [ full.SpcNumbers.nspc,   full.SpcNumbers.nsgas,  full.SpcNumbers.nsaqua,...
    full.SpcNumbers.nspart, full.SpcNumbers.nssoli, full.SpcNumbers.nspass];
n_r = [ red.SpcNumbers.nspc,    red.SpcNumbers.nsgas,   red.SpcNumbers.nsaqua, ...
    red.SpcNumbers.nspart,  red.SpcNumbers.nssoli,  red.SpcNumbers.nspass];

n_perc = zeros(1,6);
for i=1:6
    if ( n_f(i) > 0 )
        n_perc(i) = 100.0 - n_r(i)/n_f(i)*100;
    end
end
n = [ n_f ; n_r ; n_perc]';


fprintf(id,' \n');fprintf(id,' \n');fprintf(id,' \n');
fprintf(id,'           +-----------------------------------------------------------------+ \n');
fprintf(id,'           |---------------->  Analysis of species numbers  <----------------| \n');
fprintf(id,'           |--+-----------------------------------------------------------+--| \n');
fprintf(id,'           +--+                                                           +--+ \n');
fprintf(id,' \n');
fprintf(id,'                                |      full     |    reduced   |  reduction-rate    \n');
fprintf(id,'         +----------------------+---------------+--------------+------------------+ \n');
fprintf(id,'         |   number of species  |    %6d     |    %6d    |    %8.4f [%%]  |\n',n(1,:));
fprintf(id,'         +----------------------+---------------+--------------+------------------+ \n');
fprintf(id,'                 - gaseous      |    %6d     |    %6d    |    %8.4f [%%]\n',n(2,:));
fprintf(id,'                 - aqueous      |    %6d     |    %6d    |    %8.4f [%%]\n',n(3,:));
fprintf(id,'                 - paricular    |    %6d     |    %6d    |    %8.4f [%%]\n',n(4,:));
fprintf(id,'                 - solid        |    %6d     |    %6d    |    %8.4f [%%]\n',n(5,:));
fprintf(id,'                 - non-reactive |    %6d     |    %6d    |    %8.4f [%%]\n',n(6,:));

fprintf(id,' \n');
end

% Plot the statistic for species numbers
function Analyse_SpeciesDeviation(id,c_f,c_r,dev,SpcList)

nSpc = length(SpcList);

exp_thresh = -18;

fprintf(id,' \n');fprintf(id,' \n');fprintf(id,' \n');
fprintf(id,'           +----------------------------------------------------------------------------------+ \n');
fprintf(id,'           |-------------------------> Analysis of species deviation <------------------------| \n');
fprintf(id,'           |-+------------------------------------------------------------------------------+-|\n');
fprintf(id,'           +-+             printing species with concentration values > 1.0e%3d             +-+ \n',exp_thresh);
fprintf(id,'    |-------------------+---------------------------------------+---------------------------------------------+ \n');
fprintf(id,'    |                   |              daily maxima             |                   |                         | \n');
fprintf(id,'    |       phase       |       full        |      reduced      |   max. deviation  |       species names     | \n');
fprintf(id,'    |-------------------+-------------------+-------------------+-------------------|-------------------------+ \n');


for iSpc = 1:nSpc
    % add legend with maxdev
    [pks_f,locs_f] = findpeaks(c_f{iSpc});
    [pks_r,locs_f] = findpeaks(c_r{iSpc});
    
    if isempty(pks_f)
        [max_f,locs_f] = max(c_f{iSpc});
    else
        [max_f,locs_f] = max(pks_f);
    end
    if isempty(pks_r)
        [max_r,locs_r] = max(c_r{iSpc});
    else
        [max_r,locs_r] = max(pks_r);
    end
    
    
    ymax_f = max(abs(dev{iSpc}(locs_f)))*100;
    ymax_f = eval(sprintf('%.2f',ymax_f));
    
    expon_f  = floor(log10(max_f));
    expon_r  = floor(log10(max_r));
    
    if expon_f > exp_thresh || expon_r > exp_thresh
        
        fprintf(id,'      %12s      |   %12.8e  |   %12.8e  |   %10.4f [%%]  | %-40s \n', ...
                 SpcList{iSpc}.phase, max_f, max_r, ymax_f, SpcList{iSpc}.name );
        
    end
end
fprintf(id,' \n');
end

% Plot the statistic for reaction numbers
function n = Analyse_ReactionNumers(id,full,red)

n_f = [ full.ReacNumbers.nreac,         ...
    full.ReacNumbers.nrgas,         ...
    full.ReacNumbers.nrgas_photo,   full.ReacNumbers.nrgas_const,       ...
    full.ReacNumbers.nrgas_temp,    full.ReacNumbers.nrgas_simp,        ...
    full.ReacNumbers.nrgas_lind,    full.ReacNumbers.nrgas_troe,        ...
    full.ReacNumbers.nrgas_spec,    full.ReacNumbers.nrgas_special,     ...
    full.ReacNumbers.nrhenry,       ...
    full.ReacNumbers.nrdiss,        ...
    full.ReacNumbers.nrdiss_dconst, full.ReacNumbers.nrdiss_dtemp,      ...
    full.ReacNumbers.nrdiss_special,...
    full.ReacNumbers.nraqua,        ...
    full.ReacNumbers.nraqua_photo,  full.ReacNumbers.nraqua_const,      ...
    full.ReacNumbers.nraqua_temp,   full.ReacNumbers.nraqua_spec,       ...
    full.ReacNumbers.nraqua_special, ...
    full.ReacNumbers.nrparti,       full.ReacNumbers.nrparti_special,   ...
    full.ReacNumbers.nrsolid,       ...
    full.ReacNumbers.nrsolid_dtemp3,full.ReacNumbers.nrsolid_equi,      ...
    full.ReacNumbers.nrsolid_spec,  full.ReacNumbers.nrsolid_special    ...
    full.ReacNumbers.nrmicro,       full.ReacNumbers.nrmicro_special ];
n_r = [ red.ReacNumbers.nreac,         ...
    red.ReacNumbers.nrgas,         ...
    red.ReacNumbers.nrgas_photo,   red.ReacNumbers.nrgas_const,       ...
    red.ReacNumbers.nrgas_temp,    red.ReacNumbers.nrgas_simp,        ...
    red.ReacNumbers.nrgas_lind,    red.ReacNumbers.nrgas_troe,        ...
    red.ReacNumbers.nrgas_spec,    red.ReacNumbers.nrgas_special,     ...
    red.ReacNumbers.nrhenry,       ...
    red.ReacNumbers.nrdiss,        ...
    red.ReacNumbers.nrdiss_dconst, red.ReacNumbers.nrdiss_dtemp,      ...
    red.ReacNumbers.nrdiss_special,...
    red.ReacNumbers.nraqua,        ...
    red.ReacNumbers.nraqua_photo,  red.ReacNumbers.nraqua_const,      ...
    red.ReacNumbers.nraqua_temp,   red.ReacNumbers.nraqua_spec,       ...
    red.ReacNumbers.nraqua_special, ...
    red.ReacNumbers.nrparti,       red.ReacNumbers.nrparti_special,   ...
    red.ReacNumbers.nrsolid,       ...
    red.ReacNumbers.nrsolid_dtemp3,red.ReacNumbers.nrsolid_equi,      ...
    red.ReacNumbers.nrsolid_spec,  red.ReacNumbers.nrsolid_special    ...
    red.ReacNumbers.nrmicro,       red.ReacNumbers.nrmicro_special ];

n_perc = zeros(1,30);
for i=1:30
    if ( n_f(i) > 0 )
        n_perc(i) = 100.0 - n_r(i)/n_f(i)*100;
    end
end
n = [ n_f ; n_r ; n_perc]';


fprintf(id,' \n');fprintf(id,' \n');fprintf(id,' \n');
fprintf(id,'           +----------------------------------------------------------------+ \n');
fprintf(id,'           |----------------> Analysis of reaction numbers <----------------| \n');
fprintf(id,'           |--+----------------------------------------------------------+--| \n');
fprintf(id,'           +--+                                                          +--+ \n');
fprintf(id,'                                                                              \n');
fprintf(id,'                                |      full     |    reduced   |  reduction-rate   \n');
fprintf(id,'    +---------------------------+---------------+--------------+------------------+ \n');
fprintf(id,'    |  number of all reactions  |    %6d     |    %6d    |    %8.4f [%%]  |\n',n(1,:));
fprintf(id,'    +---------------------------+---------------+--------------+------------------+ \n');

if ( full.ReacNumbers.nrgas > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'             gaseous reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(2,:));
    fprintf(id,'                      - photo   |    %6d     |    %6d    |    %8.4f [%%]\n',n(3,:));
    fprintf(id,'                      - const   |    %6d     |    %6d    |    %8.4f [%%]\n',n(4,:));
    fprintf(id,'                      - temp    |    %6d     |    %6d    |    %8.4f [%%]\n',n(5,:));
    fprintf(id,'                      - simp    |    %6d     |    %6d    |    %8.4f [%%]\n',n(6,:));
    fprintf(id,'                      - lind    |    %6d     |    %6d    |    %8.4f [%%]\n',n(7,:));
    fprintf(id,'                      - troe    |    %6d     |    %6d    |    %8.4f [%%]\n',n(8,:));
    fprintf(id,'                      - spec    |    %6d     |    %6d    |    %8.4f [%%]\n',n(9,:));
    fprintf(id,'                      - special |    %6d     |    %6d    |    %8.4f [%%]\n',n(10,:));
end

if ( full.ReacNumbers.nrhenry > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'               henry reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(11,:));
end

if ( full.ReacNumbers.nrdiss > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'        dissociation reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(12,:));
    fprintf(id,'                      - dconst  |    %6d     |    %6d    |    %8.4f [%%]\n',n(13,:));
    fprintf(id,'                      - dtemp   |    %6d     |    %6d    |    %8.4f [%%]\n',n(14,:));
    fprintf(id,'                      - special |    %6d     |    %6d    |    %8.4f [%%]\n',n(15,:));
end

if ( full.ReacNumbers.nraqua > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'             aqueous reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(16,:));
    fprintf(id,'                      - photo   |    %6d     |    %6d    |    %8.4f [%%]\n',n(17,:));
    fprintf(id,'                      - const   |    %6d     |    %6d    |    %8.4f [%%]\n',n(18,:));
    fprintf(id,'                      - temp    |    %6d     |    %6d    |    %8.4f [%%]\n',n(19,:));
    fprintf(id,'                      - spec    |    %6d     |    %6d    |    %8.4f [%%]\n',n(20,:));
    fprintf(id,'                      - special |    %6d     |    %6d    |    %8.4f [%%]\n',n(21,:));
end

if ( full.ReacNumbers.nrparti > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'           paricular reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(22,:));
    fprintf(id,'                      - special |    %6d     |    %6d    |    %8.4f [%%]\n',n(23,:));
end

if ( full.ReacNumbers.nrsolid > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'               solid reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(24,:));
    fprintf(id,'                      - dtemp3  |    %6d     |    %6d    |    %8.4f [%%]\n',n(25,:));
    fprintf(id,'                      - equi    |    %6d     |    %6d    |    %8.4f [%%]\n',n(26,:));
    fprintf(id,'                      - spec    |    %6d     |    %6d    |    %8.4f [%%]\n',n(27,:));
    fprintf(id,'                      - special |    %6d     |    %6d    |    %8.4f [%%]\n',n(28,:));
end

if ( full.ReacNumbers.nrmicro > 0 )
    fprintf(id,'           ---------------------+---------------+--------------+------------------ \n');
    fprintf(id,'       microphysical reactions  |    %6d     |    %6d    |    %8.4f [%%]\n',n(29,:));
    fprintf(id,'                      - special |    %6d     |    %6d    |    %8.4f [%%]\n',n(30,:));
end

fprintf(id,' \n');
end


% gather names of diagnose species
function list = GetDiagnoseSpecies(Path,SpcList)

fileID = fopen(Path,'r');
i = 0;
while ~feof(fileID)
    
    found = contains(strtrim(fgetl(fileID)),'BEGIN_DIAG');
    
    if found
        while true
            tline = strtrim(fgetl(fileID));
            ende  = contains(tline,'END_DIAG');
            if ende
                break;
            end
            if isempty(tline) || strncmp(tline,'#',1) || ~ischar(tline)
                continue
            end
            if any(strcmp(SpcList,tline))
                i = i + 1;
                
                
                k = contains(tline,{'a','p','m'});
                if ( k )
                    list{i}.ncdf = [tline,'_1_l'];
                    list{i}.name = tline;
                    list{i}.unit = 'mol/L]';
                    list{i}.phase = 'aqueous';
                else
                    list{i}.ncdf = tline;
                    list{i}.name = tline;
                    list{i}.unit = 'mol/m$^3$]';
                    list{i}.phase = 'gaseous';
                end
            end
        end
        break;
    end
end
end


% NETCDF routines

function [f,r] = GetNetCDF_Data(p_full,p_red,var)
f = ncread( p_full , var );
r = ncread( p_red  , var );
end

% for comparion full-reduced mechanism

function IDX=cutOff(T,b1,b2)
IDX(1) = 0;
IDX(2) = 0;
len = length(T);
for i=1:len
    if ( T(i) < b1 )
        IDX(1) = IDX(1) + 1;
    else
        break;
    end
end
for i=len:-1:1
    if ( T(i) > b2 )
        IDX(2) = IDX(2) + 1;
    else
        break;
    end
end
IDX(2) = len - IDX(2);
end

function DEV=max_deviation(C_f,C_r)

n1 = length(C_f);
n2 = length(C_r);

if n1~=n2
    disp('WARNING   len(full) /= len(red)')
end

DEV = zeros(n1,1);
for i=1:n1
    DEV(i) = (C_f(i) - C_r(i)) / max(C_f(i),C_r(i));
end
end

function DEV=mean_deviation(C_f,C_r)

n1 = length(C_f);
n2 = length(C_r);

if n1~=n2
    disp('WARNING   len(full) /= len(red)')
end

full = zeros(n1,1);
red  = zeros(n1,1);
DEV  = zeros(n1,1);
for i=1:n1
    full(i) = full(i) + C_f(i);
    red(i)  = red(i) + C_r(i);
end

    DEV(:)  = (full(:) - red(:))/n1 ;
end

% PLOT FUNCTIONS
function Plot_Concentrations(t_f,t_r,c_f,c_r,dev1,Diag)

global b1 b2 len_f len_r mech

fig = figure('Units', 'normalized', 'Position', [0.2, 0.1, 0.8, 0.725]);

% marker-, fontsize
ms = 15;
fs = 30;

% subplot concentrations
% CONCENTRATION FULL MECHANISM
subplot(2,1,1);
max_c  = max(c_f);
expon  = floor(log10(max_c));
strexp = num2str(expon);

plot( t_f, c_f/(10^expon) , '-' , ...
    'LineWidth',4   , ...
    'MarkerSize',ms , ...
    'MarkerIndices',1:1:len_f); hold on;
set(gca,'fontsize',fs);
xlim([b1,b2]);
title(Diag.name,'Interpreter','latex')

str_exp = ['[$10^{',strexp,'}$ ',Diag.unit];

ylabel( {'Concentration' ; str_exp } ,'Interpreter','latex');


% CONCENTRATION ISSA MECHANISM
plot( t_r, c_r/(10^expon) , 'o' , ...
    'LineWidth',4   , ...
    'MarkerSize',ms , ...
    'MarkerIndices',1:3:len_r); hold off;
set(gca,'fontsize',fs);
xlim([b1,b2]);


legend({[mech, ' full version'], [mech, ' reduced version']}, ...
    'Location', 'best',  'Interpreter','latex');
%            [mech, ' condensed vers.'], ...

% subplot deviations
subplot(2,1,2);
plot( t_f, dev1 , '.' , ...
    'LineWidth',4   , ...
    'MarkerSize',ms , ...
    'MarkerIndices',1:1:len_f); hold on;
set(gca,'fontsize',fs);


xlim([b1,b2]);
ylabel('Deviation [\%]','Interpreter','latex');

% Convert y-axis values to percentage values by multiplication
a = cellstr(num2str(get(gca,'ytick')'*100));

new_yticks = char(a);

% 'Reflect the changes on the plot
set(gca,'yticklabel',new_yticks);
ytickformat('%.2f');

xlabel('Time [h]','Interpreter','latex');

% add legend with maxdev
[pks,locs] = findpeaks(c_f);

if isempty(pks)
    [pks,locs] = max(c_f);
else
    [pks,locs] = max(pks);
end
ymax1 = max(abs(dev1(locs)))*100;
ymax1 = eval(sprintf('%.2f',ymax1));

strmax1 = {['$\delta_{\max}$ = ',num2str(ymax1),' [\%]']};


legend(strmax1,'Location', 'best','Interpreter','latex');

%     saveas(fig, [fig_path,'Full_VS_Red__',target_spc{iSpc}] ,'png');
%

end

%% END SUBROUTINES