
close all

% %% RACM + CAPRAM 2.4
% Pfad1 = 'BA_RACM+C24.SparseMat';
% Pfad2 = 'Miter_RACM+C24.SparseMat';
% Pfad3 = 'LU_Miter_RACM+C24.SparseMat';
% Pfad4 = 'JAC_RACM+C24.SparseMat';
% Pfad5 = 'alpha_RACM+C24.SparseMat';
% Pfad6 = 'beta_RACM+C24.SparseMat';

%% Small Stato
Pfad1 = 'BA_SmallStratoKPP.SparseMat';
Pfad2 = 'Miter_SmallStratoKPP.SparseMat';
Pfad3 = 'LU_Miter_SmallStratoKPP.SparseMat';
Pfad4 = 'JAC_SmallStratoKPP.SparseMat';
Pfad5 = 'alpha_SmallStratoKPP.SparseMat';
Pfad6 = 'beta_SmallStratoKPP.SparseMat';
spc={'M  ';'NO ';'NO$_2$';'O  ';'O($^1$D)';'O$_2$ ';'O$_3$ '};

[ba,ba_sparse]       = SparseInput(Pfad1);
[m,m_sparse]         = SparseInput(Pfad2);
[lu,lu_sparse]       = SparseInput(Pfad3);
[jac,jac_sparse]     = SparseInput(Pfad4);
[alpha,alpha_sparse] = SparseInput(Pfad5);
[beta,beta_sparse]   = SparseInput(Pfad6);


spc_m  = alpha_sparse'*beta_sparse;
reac_m = alpha_sparse*beta_sparse';

big_m = [spc_m , beta_sparse'; alpha_sparse , reac_m];



fE = full(alpha_sparse);
fP = full(beta_sparse);

fEP = full(alpha_sparse'*beta_sparse);
fPE = full(beta_sparse*alpha_sparse');

fT1 = [fEP,                zeros(size(fE'));    zeros(size(fP)),    zeros(size(fPE))];
fT2 = [zeros(size(fEP)),       fP';             zeros(size(fP)),    zeros(size(fPE))];
fT3 = [zeros(size(fEP)),   zeros(size(fE'));     fE,                zeros(size(fPE))];
fT4 = [zeros(size(fEP)),   zeros(size(fE'));    zeros(size(fP)),    fPE];

figure; hold on;
spy(fT1,35);
spy(fT2,15);
spy(fT3,15);
spy(fT4,35);
reac = cell(size(fE,2),1);

for i=1:size(fE,1)
    reac(i) = {int2str(i)};
end

xtickangle(45)
set(gca, 'FontSize', 25)
set(gca,'YTick',1:1:size(fT1));
set(gca,'YTickLabel',[spc;reac]);

ytickangle(45)
set(gca,'XTick',1:1:size(fT1));
set(gca,'XTickLabel',[spc;reac]);
new = copyobj(gca,gcf);
set(new,'YAxisLocation','right','XAxisLocation','top');


