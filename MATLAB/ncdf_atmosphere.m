    %% print ncdf data for atmospheric systems:
clear all;
close all;


% paths
ncdf_path = '/Users/schimmel/Code/AtCSol/NetCDF/';
fig_path  = '/Users/schimmel/Code/AtCSol/OUTPUT/';

% mech = 'RACM+C24';
mech = 'MCM32+C40';
% mech = 'RACM+C24_constLWC';
% mech = 'RACM+C24_noDUST';
% mech = 'C24';

m_full = [mech,'_full.nc'];
m_red  = [mech,'_red.nc'];
% m_cond = [mech,'_cond.nc'];

p_full = [ncdf_path,m_full];
p_red  = [ncdf_path,m_red];
% p_cond = [ncdf_path,m_cond];

% target species
iSpc = 1;
target_spc{1}  = 'OH';
target_spc{2}  = 'NO';
target_spc{3}  = 'NO2';
target_spc{4}  = 'NO3';
target_spc{5}  = 'O3';
target_spc{6}  = 'aH2O2_1_l';
target_spc{7}  = 'aHO_1_l';
target_spc{8}  = 'aNO3_1_l';
target_spc{9}  = 'aSO2_1_l';
target_spc{10} = 'PAN';
target_spc{11} = 'CO2';
target_spc{12} = 'O1D';
% target_spc{13} = 'OHm_1_l';
% target_spc{14} = 'Hp_1_l';
% target_spc{15} = 'SULF';

% Time
T_f = ncread( p_full , 'time' );
T_r = ncread( p_red  , 'time' );

 % begin time, end time
%  b1 = T_f(1)   + 1.0;   % 12 uhr mittags
%  b2 = T_f(end) - 1.0;   % 12 uhr mittags n�chster tag
 b1 = 12.0;   % 12 uhr mittags
 b2 = 36.0;   % 12 uhr mittags n�chster tag

idx_f = cutOff(T_f,b1,b2);
idx_r = cutOff(T_r,b1,b2);

T_f = T_f(idx_f(1):idx_f(2));
T_r = T_r(idx_f(1):idx_f(2));

for iSpc=1:length(target_spc)
    k = strfind(target_spc{iSpc},'_');
    if ( ~isempty(k) )
        if ( k(1) > 0 )
            spc_label{iSpc} = target_spc{iSpc}(1:k(1)-1);
            unit = 'mol/L]';
        end
    else
        spc_label{iSpc} = target_spc{iSpc};
        unit = 'mol/m$^3$]';
    end
    

    % Concentration
    C_f = ncread( p_full , target_spc{iSpc});
    C_r = ncread( p_red  , target_spc{iSpc});
    
    C_f = C_f(idx_f(1):idx_f(2));
    C_r = C_r(idx_f(1):idx_f(2));

    len_f = length(T_f);
    len_r = length(T_r);
    
    % calculate deviations
    dev1 = deviation(C_f,C_r);
   
    
    %% Plot
    ms = 15;
    fs = 30;
    fig = figure('Units', 'normalized', 'Position', [0.2, 0.1, 0.8, 0.725]);
    
    % subplot concentrations
    % CONCENTRATION FULL MECHANISM
    subplot(2,1,1);
    max_c  = max(C_f);
    expon  = floor(log10(max_c));
    strexp = num2str(expon);
        
    plot( T_f, C_f/(10^expon) , '-' , ...
        'LineWidth',4   , ...
        'MarkerSize',ms , ...
        'MarkerIndices',1:1:len_f); hold on;
    set(gca,'fontsize',fs);
    xlim([b1,b2]);
    title(spc_label{iSpc},'Interpreter','latex')
    
    str_exp = ['[$10^{',strexp,'}$ ',unit];
    
    ylabel( {'Concentration' ; str_exp } ,'Interpreter','latex');
%     ylabel(['Concentration', '[10^{-12}'],'Interpreter','latex');
    
%     % CONCENTRATION CONDENSED MECHANISM
%     plot( T_c, C_c , '--' , ...
%         'LineWidth',4   , ...
%         'MarkerSize',ms , ...
%         'MarkerIndices',1:2:len_r); hold on;
%     set(gca,'fontsize',fs);
%     xlim([b1,b2]);
    
    % CONCENTRATION ISSA MECHANISM
    plot( T_r, C_r/(10^expon) , 'o' , ...
        'LineWidth',4   , ...
        'MarkerSize',ms , ...
        'MarkerIndices',1:3:len_r); hold off;
    set(gca,'fontsize',fs);
    xlim([b1,b2]);
    
    
    legend([mech, ' full version'],    ...
           [mech, ' ISSA-reduced vers.'], ...
           'Location', 'best', ...
           'Interpreter','latex');    
%            [mech, ' condensed vers.'], ...
       
    % subplot deviations
    subplot(2,1,2);
    plot( T_f, dev1 , '.' , ...
        'LineWidth',4   , ...
        'MarkerSize',ms , ...
        'MarkerIndices',1:1:len_f); hold on;
    set(gca,'fontsize',fs);
    
%     plot( T_c, dev2 , '--' , ...
%         'LineWidth',4   , ...
%         'MarkerSize',ms , ...
%         'MarkerIndices',1:5:len_f); hold on;
%     set(gca,'fontsize',fs);
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
    [pks,locs] = findpeaks(C_f);
    
    if isempty(pks)
        [pks,locs] = max(C_f);
    end
    ymax1 = max(abs(dev1(locs)))*100;
    ymax1 = eval(sprintf('%.2f',ymax1));
    strmax1 = ['$\delta_{\max}$ = ',num2str(ymax1),' [\%]'];
    
    
    legend(strmax1,'Location', 'best','Interpreter','latex');
    
%     saveas(fig, [fig_path,'Full_VS_Red__',target_spc{iSpc}] ,'png');
%     
end


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
