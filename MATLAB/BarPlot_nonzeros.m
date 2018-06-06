%% nonzeros bar plot
clear all;
close all;


leg = {'standard','extended'};

solver{1} = 'sparse-LU';
solver{2} = 'MUMPS';

strat{1} = 'MD';
strat{2} = 'AMF';
strat{3} = 'AMD';
strat{4} = 'QAMD';
strat{5} = 'Metis';

nnz_std = [118356	577708	597638	651178	663958];
nnz_ext = [262250	921800	969018	1035790	1337674];

nnz = [nnz_std; nnz_ext]';


fs = 25;
fig = figure('Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.7]);


c = categorical(strat);
c = reordercats(c,strat);


blue1 = [65,105,225]/255;
blue2 = [176,196,222]/255;
clr = [blue1;blue2];

colormap(clr);


bar(c,nnz);
% xlabel('Ordering strategie','FontSize', fs, 'Interpreter','latex');
ylabel('Number of nozeros','FontSize', fs, 'Interpreter','latex');
set(gca,'fontsize',fs);

ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.TickLabelInterpreter='latex';

% xtickangle(30)

legend(leg, 'Location', 'northwest','FontSize', fs, 'Interpreter','latex');
