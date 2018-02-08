close all;

%% paths to matices alpha and beta, as well as path to .chem file
% reac_str{1} = 'O2 = 2.0 O';
% reac_str{2} = 'O + O2 = O3';
% reac_str{3} = 'O3 = O + O2';
% reac_str{4} = 'O + O3 = 2.0 O2';
% reac_str{5} = 'O3 = O1D + O2';
% reac_str{6} = 'O1D + M = O + M';
% reac_str{7} = 'O1D + O3 = 2.0 O2';
% reac_str{8} = 'NO + O3 = NO2 + O2';
% reac_str{9} = 'NO2 + O = NO + O2';
% reac_str{10} = 'NO2 = NO + O';


%% example: SmallStratoKPP
path_P = 'alpha_RACM+C24.SparseMat';
path_E = 'beta_RACM+C24.SparseMat';
path_chem = 'CHEM/RACM+C24.chem';
fs = 16; ds_big = 35; ds_small = 15; lw = 3;
plots = true;

%% example: INORG
% path_P = 'P_HG_INORG.SparseMat';
% path_E = 'E_HG_INORG.SparseMat';
% path_chem = 'CHEM/INORG.chem';
% fs = 5; ds_big = 20; ds_small = 10; lw = 1;

%% example: CAPRAM
% path_P = 'P_HG_CAPRAM.SparseMat';
% path_E = 'E_HG_CAPRAM.SparseMat';
% path_chem = 'CHEM/CAPRAM24.chem';
% fs = 4; ds_big = 10; ds_small = 10; lw = 1;

%% figure window positions
% figure predifine location and size:    [ x , y , width , hight]
% second screen top right
figPosM = [6.464285714285715e-01 , 1.658095238095238e+00,...
    4.999999999999999e-01 , 6.000000000000001e-01];

% second screen bottom right
figPosG = [6.464285714285715e-01 , 1.000000000000000e+00,...
    4.999999999999999e-01 , 6.000000000000001e-01];


%% generate sparse matrices of stoech. coefs of  educts and products

[E,sE] = SparseInput(path_E); 
[P,sP] = SparseInput(path_P); 
sE = sE';

spc = ReadSpeciesNames(path_chem);

%% building the Hypergraph matrix H
[ns,nr] = size(sE);
zeroNSNS = sparse(zeros(ns,ns));
zeroNSNR = sparse(zeros(ns,nr));
zeroNRNS = sparse(zeros(nr,ns));
zeroNRNR = sparse(zeros(nr,nr));

sEP = sE * sP;
sPE = sP * sE;

% change all values to 1
sEP(sEP~=0) = 1;
sPE(sPE~=0) = 1;

sT1 = [sEP,      zeroNSNR;  zeroNRNS, zeroNRNR];
sT2 = [zeroNSNS,       sE;  zeroNRNS, zeroNRNR];
sT3 = [zeroNSNS, zeroNSNR;     sP,    zeroNRNR];
sT4 = [zeroNSNS, zeroNSNR; zeroNRNS,   sPE    ];


% %% save species names in reac-array for plot
% reac = cell(nr,1);
% for i=1:nr
%     reac(i) = {int2str(i)};
% end


nLayers = 8;
widLine = 1;
xPos = 1 : 5 : 5*nLayers;
yPos = 1 : ns;

%% plotting a digraph
if plots
    figure('units','normalized','outerposition',[0 0 1 1]);
    % gJAC = digraph( sEP , spc  );
    gJAC = digraph( sEP );
%     g = plot( gJAC );
    h = plot(gJAC,'NodeLabel',[]);
%     for i=1:length(spc)
%         text(h.XData(i)+0.1,h.YData(i),spc(i),'fontsize',20);
%     end
    
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    % gJAC = digraph( sEP , spc  );
    gJAC1 = digraph( sPE );
%     g = plot( gJAC );
    h = plot(gJAC1,'NodeLabel',[]);
%     for i=1:length(reac_str)
%         text(h.XData(i)+0.1,h.YData(i),reac_str{i},'fontsize',20);
%     end
    
end



