clc

%% Important species declared by user
%
%  - S_imp( 1) =     57     HO
%  - S_imp( 2) =     79     NO
%  - S_imp( 3) =     80     NO2
%  - S_imp( 4) =     83     O3
%  - S_imp( 5) =    211     aH2O2
%  - S_imp( 6) =    216     aHO
%  - S_imp( 7) =    228     aNO3
%  - S_imp( 8) =    242     aSO2
%  - S_imp( 9) =    246     Hp

s_imp = [57, 79, 80, 83, 211, 216, 228, 242, 246];

%% paths to matrices and cycles
%
mat_path = '/Users/schimmel/Code/BACKUP_241017_ISSA/chemie/MATRICES/';

ba_path  = [mat_path 'BA_RACM+C24_full.SparseMat'];
a_path   = [mat_path 'alpha_RACM+C24_full.SparseMat'];
b_path   = [mat_path 'beta_RACM+C24_full.SparseMat'];
m_path   = [mat_path 'Miter_RACM+C24_full_CL.SparseMat'];
lum_path = [mat_path 'LU_Miter_RACM+C24_full_CL.SparseMat'];

cyc_path = '/Users/schimmel/Code/BACKUP_241017_ISSA/matlab/RACM_C24_cycles.txt';

%% build sparse matricies
%
[a,sp_a]     = SparseInput(a_path);
[b,sp_b]     = SparseInput(b_path);
[ba,sp_ba]   = SparseInput(ba_path);
[m,sp_m]     = SparseInput(m_path);
[lum,sp_lum] = SparseInput(lum_path);


%% build cycles
cyclic_set = ReadCycles(cyc_path);


%% search reations and print them
iC = 1;

for iC = 1: 11
    for i = 1:cyclic_set(iC).len-1
        iE = cyclic_set(iC).Species(i);
        iP = cyclic_set(iC).Species(i+1);
        
        idxE = find(sp_a(:,iE)>0);
        idxP = find(sp_b(:,iP)>0);
        
        tidx = find(ismember(idxE,idxP));
        
        cyclic_set(iC).Reactions = idxE(tidx);
        g = sprintf('%d ', idxE(tidx));
        
        fprintf('iC = %4d :: E = %4d  -->  P = %4d  Reacs :: %s\n',iC,iE,iP,g);
    end
    fprintf('\n');
end















