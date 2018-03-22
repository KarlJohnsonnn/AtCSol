%% print ncdf data for atmospheric systems:
clear all;
close all;

% paths
ncdf_path = '/Users/schimmel/Code/AtCSol/NetCDF/';
fig_path  = '/Users/schimmel/Code/AtCSol/OUTPUT/';

mech = 'LLNL_MD';
% mech = 'ERC_nheptane';

m_full = [mech,'_full.nc'];

p_full = [ncdf_path,m_full];

% target species
iSpc = 1;
target_spc{1}  = 'h2o';
target_spc{2}  = 'co';
target_spc{3}  = 'co2';
target_spc{4}  = 'oh';
target_spc{5}  = 'md';
% target_spc{5}  = 'nc7h16';
target_spc{6}  = 'Temperature';


%% Plot
ms = 35;
fs = 35;

% Time
T_f = ncread( p_full , 'time' );
len_f = length(T_f);

% begin time, end time
b1 = 0.0;
b2 = 0.04;

% color map
map = [204, 0.0, 0.0
       204, 204, 0.0
       0.0, 204, 0.0
       0.0, 204, 204
       0.0, 0.0, 204
        50,  50,  50]/255;

fig = figure('Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.7]);
set(fig,'defaultAxesColorOrder',[map(5,:); map(6,:)]);    




for iSpc=1:length(target_spc)
    
    
    spc_label{iSpc} = upper(target_spc{iSpc});
    unit{iSpc} = '-]';
    
    
    % Concentration
    C_f{iSpc} = ncread( p_full , target_spc{iSpc});
end
spc_label{6} = 'Temperature';


s_out{1} = ['\ ',spc_label{1}];  u_out{1} = unit{1};
s_out{2} = ['\ ',spc_label{2}];  u_out{2} = unit{2};
s_out{3} = ['\ ',spc_label{3}];  u_out{3} = unit{3};
s_out{4} = ['\ ',spc_label{4}];  u_out{4} = unit{4};
s_out{5} = ['\ ',spc_label{5}];  u_out{5} = unit{5};
s_out{6} = ['\ ',spc_label{6}];  u_out{6} = 'K]';
fac(1) = 1; 
fac(2) = 1; 
fac(3) = 1; 
fac(4) = 1; 
fac(5) = 1; 
fac(6) = 1; 

% colorstring = 'kbgrym';

C_f_out = [C_f{1},C_f{2},C_f{3},C_f{4},C_f{5},C_f{6}];


%% species danach
yyaxis left

   
    
    
for iSpc=1:5
    
    plot( T_f, fac(iSpc)*C_f_out(:,iSpc), ...
        '-' , ...
        'Color',map(iSpc,:), ...
        'LineWidth',4   , ...
        'MarkerSize',ms , ...
        'MarkerIndices',1:5:len_f); hold on;
end

set(gca,'fontsize',fs);
set(gca, 'YScale', 'log')

ylabel('mole fraction [-]','Interpreter','latex','FontSize',fs);
xlabel('Time [s]','Interpreter','latex','FontSize',fs);
xticks([0 0.01 0.02 0.03 0.04]);

xlim([b1,b2]);
ylim([1.0e-12,10]);


yyaxis right
plot( T_f, fac(6)*C_f_out(:,6) , ...
        '.-' , ...
        'Color',map(6,:), ...
        'LineWidth',4   , ...
        'MarkerSize',ms , ...
        'MarkerIndices',1:100:len_f); hold off;
    
ylabel('Temperature [K]','Interpreter','latex','FontSize',fs);
%     
% ylim([1.0e-12,1]);

legend(s_out, 'Location', 'east', 'Interpreter','latex', 'FontSize',fs);
% legend(s_out{end}, 'Location', 'east', 'Interpreter','latex', 'FontSize',fs);
% 


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
