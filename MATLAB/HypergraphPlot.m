close all;

%% paths to matices alpha and beta, as well as path to .chem file

%       O2 = 2.0 O
%   O + O2 = O3
%       O3 =   O + O2
%   O + O3 = 2.0 O2
%       O3 = O1D + O2
% O1D + M  =   O + M
% O1D + O3 = 2.0 O2
%  NO + O3 = NO2 + O2
% NO2 + O  =  NO + O2
%      NO2 =  NO + O

%% example: SmallStratoKPP
path_P = 'P_HG_KPP.SparseMat';
path_E = 'E_HG_KPP.SparseMat';
path_chem = 'CHEM/SmallStratoKPP.chem';
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
E = SparseInput(path_E); sE = sparse( E.ri, E.ci, E.val, E.m, E.n );
P = SparseInput(path_P); sP = sparse( P.ri, P.ci, P.val, P.m, P.n );

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


%% save species names in reac-array for plot
reac = cell(nr,1);
for i=1:nr
    reac(i) = {int2str(i)};
end


%% first figure predifine location and size:     x     y    wid   hig
if plots
    fG = figure('units','normalized', 'Position', figPosM);
    hold on;
    
    %  spy plot - Hypergraph H = [sEP,sE; sP,sPE]
    spy(sT1,ds_big);      spy(sT2,ds_small);
    spy(sT3,ds_small);    spy(sT4,ds_big);
    
    
    % figure settings
    set(gca, 'FontSize', fs)                % x,y label size
    set(gca, 'YTick',1:1:ns+nr);        % nr+ns y-ticks
    set(gca, 'YTickLabel',[spc;reac]);      % ispc + ireac
    set(gca, 'XTick',1:1:ns+nr);        % nr+ns x-ticks
    set(gca, 'XTickLabel',[spc;reac]);      % ispc + ireac
    new = copyobj(gca,gcf);                 % copy axis label to the top and right
    set(new,'YAxisLocation','right','XAxisLocation','top');
    hold off;
end
% %
% %% print species and reaction number in which it is produced
% for j=1:ns
%     for i=1:nr
%         if ( sP(i,j) * sE(j,i) > 0.0)
%             fprintf('Species:  %d (= %s)  wird hergestellt in Reaktion %d\n',j,spc{j},i);
%         end
%     end
% end



%% plotting bipartite graph using BiMat

% Species-Matrix graph
bp_EP = Bipartite(sEP);
bp_EP.row_labels = spc;
bp_EP.col_labels = spc;
bp_EP.row_class = 1:1:ns ;
bp_EP.col_class = 1:1:ns ;

% Reaction-Matrix graph
bp_PE = Bipartite(sPE);
bp_PE.row_labels = reac;
bp_PE.col_labels = reac;
bp_PE.row_class = 1:1:nr;
bp_PE.col_class = 1:1:nr;


% figure('units','normalized','Position',figPosG)
% subplot(1,2,1); bp_EP.plotter.PlotModularGraph();
% subplot(1,2,2); bp_PE.plotter.PlotModularGraph();

%% zum draufrumkitzeln
% figure('units','normalized','Position',figPosG)
% for i = 1 : 4
%     subplot(2,2,i); hold on;
%     spy(sT1,ds_big);      spy(sT2,ds_small);
%     spy(sT3,ds_small);    spy(sT4,ds_big);
%     set(gca, 'FontSize', fs)                % x,y label size
%     set(gca, 'YTick',1:1:ns+nr);        % nr+ns y-ticks
%     set(gca, 'YTickLabel',[spc;reac]);      % ispc + ireac
%     set(gca, 'XTick',1:1:ns+nr);        % nr+ns x-ticks
%     set(gca, 'XTickLabel',[spc;reac]);      % ispc + ireac
%     new = copyobj(gca,gcf);                 % copy axis label to the top and right
%     set(new,'YAxisLocation','right','XAxisLocation','top');
%
% end
% saveas(gcf,'test.png')
% export_fig test2.png

% figure;
% bp_sE = Bipartite(sT2+sT4);
% bp_sE.row_labels = spc2;
% bp_sE.col_labels = reac;
% bp_sE.row_class = [1 2 3 4 5 6 7] ;
% bp_sE.col_class = [1 2 3 4 5 6 7 8 9 10] ;
% bp_sE.plotter.font_size = 12.0;
% bp_sE.plotter.PlotModularGraph();



% % figure predifine location and size:     x     y    wid   hig
% figure('units','normalized','position',[ 0.3   0.8   0.75   0.6])
% bp_EP.plotter.use_type_interaction = true;
% bp_EP.plotter.color_interactions(1,:) = [1 0 0]; %Red color for clear lysis
% bp_EP.plotter.color_interactions(2,:) = [0 0 1]; %Blue color for turbid spots
% bp_EP.plotter.color_interactions(3,:) = [0 1 0]; %Green color for turbid spots
% bp_EP.plotter.color_interactions(4,:) = [0 1 1]; %Blue color for turbid spots
% bp_EP.plotter.color_interactions(5,:) = [1 1 0]; %yellow color for turbid spots
% bp_EP.plotter.back_color = 'white';
% bp_EP.plotter.use_isocline = true; %The isocline is used in the NTC algorithm.
% bp_EP.plotter.isocline_color = 'red';
% % After changing all the format we finally can call the plotting function.
% subplot(1,2,1); bp_EP.plotter.PlotNestedMatrix();
%
% bp.plotter.plot_iso_modules = true;
% subplot(1,2,2); bp_EP.plotter.PlotModularMatrix();
nLayers = 8;
widLine = 1;
xPos = 1 : 5 : 5*nLayers;
yPos = 1 : ns;


%% plotting a network of length nLayers using the matrix E*P
% figure('units','normalized','Position',figPosG); hold on;
% Smarker = 75;
% % plot(xPos,yPos);
% [row,col] = find(sEP);
%
% dp = -4;
%
% set(gca, 'YTick',1:1:ns);        % nr+ns y-ticks
% set(gca, 'YTickLabel',spc);      % ispc ##
%
% xlim([0 5*nLayers-1]);
%
%
% for layer = 1 : nLayers-1
%     subplot(2,1,1);
%     plot( [xPos(layer) xPos(layer+1)] , [row(:) col(:)] , '-k' , ...
%         'LineWidth' , widLine ); hold on;
%     xlim([0 5*nLayers-1]);
% end
% set(gca, 'YTick',1:1:ns);        % nr+ns y-ticks
% set(gca, 'YTickLabel',spc);
% for layer = 1 : nLayers
%     subplot(2,1,1);
%     plot(  xPos(layer) , yPos(:) , '.' , 'MarkerSize' , Smarker );hold on;
%     xlim([0 5*nLayers-1]);
% end
% set(gca, 'YTick',1:1:ns);        % nr+ns y-ticks
% set(gca, 'YTickLabel',spc);
% new = copyobj(gca,gcf);                 % copy axis label to the top and right
% set(new,'YAxisLocation','right');
%
%
%
% for layer = 1 : nLayers
%     subplot(2,1,2);
%     plot(  xPos(layer) , yPos(:) , '.' , 'MarkerSize' , Smarker );hold on
%     xlim([0 5*nLayers-1]);
% end
% set(gca, 'YTick',1:1:ns);        % nr+ns y-ticks
% set(gca, 'YTickLabel',spc);
% new = copyobj(gca,gcf);                 % copy axis label to the top and right
% set(new,'YAxisLocation','right');


%% plotting a digraph
if plots
    figure('units','normalized','outerposition',[0 0 1 1]);
    % gJAC = digraph( sEP , spc  );
    gJAC = digraph( sEP );
    plot( gJAC );
end
% for i=1:ns
%    sEP(i,i) = 0;
% end

A = sparse([2 3 4 5 5 6 6 7 8 4 9 5 10 6 9], ...
[1 2 2 3 4 3 5 6 4 8 8 9 9 10 6], ...
ones(1,15));

disp('================================================')
disp(' Calculating the elementary circuits in sEP:');
disp('================================================')
disp(' ')

t0 = cputime;
% [nCycles , cycles] = find_elem_circuits(A);
[nCycles , cycles] = find_elem_circuits(sEP);
te = cputime - t0;
if (nCycles>0)
    for i=1:size(cycles,2)
        fprintf(' %d ',cycles{1,i});
        fprintf(' \n ');
    end
end
fprintf(' Time needed for calculating the cirucits = %d \n',te);

