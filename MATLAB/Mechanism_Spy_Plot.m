%% paths and matrices
% close all

p_atcsol = '/Users/schimmel/Code/AtCSol/MATRICES/';

mechanism = cell(6,1);

mechanism(1) = {'SmallStratoKPP'};
mechanism(2) = {'RACM+C24'};
mechanism(3) = {'MCM32+CAPRAM40'};
mechanism(4) = {'ERC_nheptane'};
mechanism(5) = {'nheptane'};
mechanism(6) = {'LLNL_MD'};

linAlg = 'cl';
% linAlg = 'ex';


j = 1;

%% spy plots
fs = 20;
figure1 = figure('Units', 'normalized', 'Position', [0.0, 0.0, 1.0, 1.0],...
    'InvertHardcopy','off','Color',[1 1 1],...
    'Renderer','painters');

[~,sp_m]     = SparseInput([p_atcsol,'Miter_',mechanism{j},'_',linAlg,'.SparseMat']);

subplot(1,2,1);
spy(sp_m);      hold on;
%     xlabel([]);
xlabel(['nz = ',num2str(nnz(sp_m))], 'Interpreter','latex');
ylabel([]);     
set(gca,'xtick',[]);
set(gca,'ytick',[]);

set(gca,'color','none','FontSize',fs)

[~,sp_lum]     = SparseInput([p_atcsol,'LU_Miter_',mechanism{j},'_',linAlg,'.SparseMat']);

subplot(1,2,2);
spy(sp_lum);    hold on;
xlabel(['nz = ',num2str(nnz(sp_lum))], 'Interpreter','latex');
ylabel([]);     
set(gca,'xtick',[]);
set(gca,'ytick',[]);

set(gca,'color','none','FontSize',fs)

