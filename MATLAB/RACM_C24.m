%% print ncdf data for atmospheric systems:
clear all;
close all;

% paths
ncdf_path = '/Users/schimmel/Code/AtCSol/NetCDF/';
fig_path  = '/Users/schimmel/Code/AtCSol/OUTPUT/';

mech = 'RACM+C24';
% mech = 'C24';
% mech = 'MCM32+C40';

m_full = [mech,'_full.nc'];
% m_cond = [mech,'_cond.nc'];

p_full = [ncdf_path,m_full];
% p_red  = [ncdf_path,m_red];
% p_cond = [ncdf_path,m_cond];

% target species
iSpc = 1;
target_spc{1}  = 'HO';
target_spc{2}  = 'NO';
target_spc{3}  = 'NO2';
target_spc{4}  = 'NO3';
target_spc{5}  = 'O3';
target_spc{6}  = 'aHO_1_m3';
target_spc{7}  = 'aNO3_1_m3';
target_spc{8}  = 'aSO2_1_m3';


%
% target_spc{1}  = 'OH';
% target_spc{2}  = 'NO';
% target_spc{3}  = 'NO2';
% target_spc{4}  = 'NO3';
% target_spc{5}  = 'O3';
% target_spc{6}  = 'aHO_1_l';
% target_spc{7}  = 'aNO3_1_l';
% target_spc{8}  = 'aSULF_1_l';


%% Plot
ms = 15;
fs = 35;

% Time
T_f = ncread( p_full , 'time' );
len_f = length(T_f);

% begin time, end time
b1 = 1.0;   % 12 uhr mittags
b2 = 23.0;   % 12 uhr mittags nächster tag


fig = figure('Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.7]);

for iSpc=1:length(target_spc)
    k = strfind(target_spc{iSpc},'_');
    if ( ~isempty(k) )
        if ( k(1) > 0 )
            spc_label{iSpc} = target_spc{iSpc}(1:k(1)-1);
            unit{iSpc} = 'mol/L]';
        end
    else
        spc_label{iSpc} = target_spc{iSpc};
        unit{iSpc} = 'mol/m$^3$]';
    end
    
    % Concentration
    C_f{iSpc} = ncread( p_full , target_spc{iSpc});
end

C_f_out{1} = C_f{2}+C_f{3}+C_f{4};
s_out{1}   = '\ NOx';  u_out{1} = unit{2};
C_f_out{2} = C_f{1};        s_out{2} = ['\ ',spc_label{1}];    u_out{2} = unit{1};
C_f_out{3} = C_f{5};        s_out{3} = ['\ ',spc_label{5}];    u_out{3} = unit{5};
C_f_out{4} = C_f{6};        s_out{4} = ['\ ',spc_label{6}];    u_out{4} = unit{6};
C_f_out{5} = C_f{7}*6.0;    s_out{5} = ['\ ',spc_label{7}];    u_out{5} = unit{7};
fac(1) = 2;
fac(2) = 1;
fac(3) = 1;
fac(4) = 3;
fac(5) = 3;

%% sonne und wolke zuerst
sonne = ncread( p_full , 'Zenith');
sonne( sonne > pi/2.0 ) =  pi/2.0;
sonne = cos(sonne);
sonne(1) = 0;

patch(  'XData', T_f, ...
    'YData', sonne*4, ...
    'FaceColor', [255,215,0]/255, ...
    'FaceAlpha', 0.5);   hold on;

wolke = ncread( p_full , 'wetRadius_1');

[w_min,idx_min] = min(wolke);
[w_max,idx_max] = max(wolke);
t_min = T_f(idx_min);
t_max = T_f(idx_max);
% 
wolke_t     = [T_f(1), T_f(8), T_f(9), T_f(65), T_f(66), ...
    T_f(79), T_f(81), T_f(137), T_f(139),     ...
    T_f(151), T_f(152), T_f(1), T_f(1),];
wolke_patch = [w_max,  w_max,  w_min,  w_min,   w_max,  ...
    w_max,   w_min,   w_min,   w_max,       ...
    w_max,   0.0,    0.0, w_max];

patch(  'XData', wolke_t, ...
    'YData', wolke_patch*1.5e4, ...
    'FaceColor', [135,206,250]/255, ...
    'FaceAlpha', 0.5);    hold on;


xlim([b1,b2]);
ylim([0,7]);

ylabel('Concentration','Interpreter','latex','FontSize',fs);
xlabel('Time [h]','Interpreter','latex','FontSize',fs);


%% species danach

for iSpc=1:length(C_f_out)
    
    
    % subplot concentrations
    % CONCENTRATION FULL MECHANISM
    
    max_c = max(C_f_out{iSpc});
    expon = floor(log10(max_c)) ;
    str_exp = num2str(expon);
    
    s_out{iSpc} =[s_out{iSpc},' [',num2str(fac(iSpc)),'e',str_exp,' ',u_out{iSpc}];
    
    
    plot( T_f, fac(iSpc)*C_f_out{iSpc}/(10^expon) , '-' , ...
        'LineWidth',4   , ...
        'MarkerSize',ms , ...
        'MarkerIndices',1:5:len_f); hold on;
    set(gca,'fontsize',fs);
    %     title([ spc_label{iSpc}, ' , scenario URBAN' ])
    
    %     ylabel(['Concentration in [10^{',num2str(expon),unit]);
    
    
    
    
end

legend(['\ solar radiation','\ cloud phase',s_out], ...
    'Location', 'northwest',...
    'Interpreter','latex',  ...
    'FontSize',fs)



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

function DEV=deviation(C_f,C_r)

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
