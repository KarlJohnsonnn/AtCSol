
spc = {'M  ';'NO ';'NO2';'O  ';'O1D';'O2 ';'O3 '};
BA_dat = 'BA_SmallStratoKPP_ex.SparseMat';
JAC_dat = 'JAC_SmallStratoKPP_ex.SparseMat';

% BA_dat = 'BA_INORG_cl.SparseMat';
% JAC_dat = 'JAC_INORG_cl.SparseMat';

% BA_dat = 'BA_ERC_nheptane_cl.SparseMat';
% JAC_dat = 'JAC_ERC_nheptane_cl.SparseMat';

% BA_dat = 'BA_MCM+CAPRAM_cl.SparseMat';
% JAC_dat = 'JAC_MCM+CAPRAM_cl.SparseMat';



% BA = SparseInput(BA_dat);
% spBA = sparse(BA.ri,BA.ci,BA.val);

JAC     = SparseInput(JAC_dat);
JAC.val = 1.0d0;
spJAC   = sparse(JAC.ri,JAC.ci,JAC.val);
spJAC   = spJAC - sparse(diag(ones(1,JAC.ns)));

figure('units','normalized','outerposition',[0 0 1 1]);
gJAC = digraph( spJAC , spc , 'OmitSelfLoops' );
plot( gJAC );


