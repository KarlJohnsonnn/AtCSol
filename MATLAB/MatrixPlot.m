
close all

% Pfad1 = 'Miter_ERC_nheptane.SparseMat';
% Pfad2 = 'LU_Miter_ERC_nheptane.SparseMat';

% Pfad1 = 'Miter_LLNL_MD_cl.SparseMat';
% Pfad2 = 'LU_Miter_LLNL_MD_cl.SparseMat';
% Pfad3 = 'Miter_LLNL_MD_ex.SparseMat';
% Pfad4 = 'LU_Miter_LLNL_MD_ex.SparseMat';

Pfad1 = 'BA_SmallStratoKPP_ex.SparseMat'; spc=['M  ';'NO ';'NO2';'O  ';'O1D';'O2 ';'O3 '];
Pfad2 = 'Miter_SmallStratoKPP_ex.SparseMat';
% Pfad3 = 'LU_Miter_SmallStratoKPP_ex.SparseMat';
Pfad4 = 'JAC_SmallStratoKPP_ex.SparseMat';
Pfad3 = 'BA_SmallStratoKPP_ex.SparseMat';


Pfad6 = 'P_HG_KPP.SparseMat';P = SparseInput(Pfad6);spP = sparse(P.ri,P.ci,P.val);
Pfad7 = 'E_HG_KPP.SparseMat';E = SparseInput(Pfad7);spE = sparse(E.ri,E.ci,E.val);


fE = full(spE);
fP = full(spP);
fEP = full(spE*spP);
fPE = full(spP*spE);
fT1 = [fEP,                zeros(size(fE));    zeros(size(fP)),    zeros(size(fPE))];
fT2 = [zeros(size(fEP)),       fE;             zeros(size(fP)),    zeros(size(fPE))];
fT3 = [zeros(size(fEP)),   zeros(size(fE));     spP,                zeros(size(fPE))];
fT4 = [zeros(size(fEP)),   zeros(size(fE));    zeros(size(fP)),    fPE];

figure; hold on;
spy(fT1,35);
spy(fT2,15);
spy(fT3,15);
spy(fT4,35);
reac = cell(size(fE,2),1);

for i=1:size(fE,2)
    reac(i) = {int2str(i)};
end

set(gca, 'FontSize', 16)
set(gca,'YTick',1:1:size(fT1));
set(gca,'YTickLabel',[spc;reac]);

set(gca,'XTick',1:1:size(fT1));
set(gca,'XTickLabel',[spc;reac]);
new = copyobj(gca,gcf);
set(new,'YAxisLocation','right','XAxisLocation','top');


% BA_INORGdat = 'BA_INORG_cl.SparseMat';
% JAC_INORGdat = 'JAC_INORG_cl.SparseMat';
% BA_INORG = SparseInput(BA_INORGdat);
% JAC_INORG = SparseInput(JAC_INORGdat);


MS1 = SparseInput(Pfad1);
MS2 = SparseInput(Pfad2);
MS3 = SparseInput(Pfad3);
% JAC = SparseInput(Pfad4);
% % MS4 = SparseInput(Pfad4);
% JAC.val = 1.0d0;
% 
spMS1 = sparse(MS1.ri,MS1.ci,MS1.val);
spMS2 = sparse(MS2.ri,MS2.ci,MS2.val);
spMS3 = sparse(MS3.ri,MS3.ci,MS3.val);
% spJAC = sparse(JAC.ri,JAC.ci,JAC.val);
% spJAC = spJAC - sparse(diag(ones(1,JAC.ns)));
% fullJAC=full(spJAC);
% % spMS4 = sparse(MS4.ri,MS4.ci,MS4.val);

figure; spy(spMS3);
BA_st = full(spMS3);
summe = zeros(size(spMS3,2),1);
for i = 1 : size(BA_st,2)
    summe(i) = sum( BA_st(1:7,i) );
end
summe(abs(summe)<1.e-5) = 0;
round(summe) 
% 
% % finding values 
% consuming = spMS1 < 0.0e0;
% producing = spMS1 > 0.0e0;
% katalysator = (spMS1 == 4.440892098500626e-16);
% 
% 
% % spy(spMS1);hold on;
% figure;
% gJAC = digraph(spJAC);
% e = gJAC.Edges;
% plot(gJAC);
% 
% figure;
% spy(consuming,15);hold on;
% spy(producing,25);hold on;
% spy(katalysator,35);
% set(gca,'XTick',1:1:size(spMS1,2))
% set(gca,'XTickLabel',spc);
% set(gca,'YTick',1:1:size(spMS1,1));
% set(gca,'YTickLabel',1:1:size(spMS1,1))
% % subplot(2,2,4);spy(spMS4);
% 
% y=5; %dummyx
% cnt=0;
% 
% % zeige matrix einträge /= 1.0d0
% for i = 1: MS1.nnz
%     if (MS1.val(i)>1.00000001d0 || MS1.val(i)<0.99999999d0)
%         cnt = cnt+1;
%         rMSx(cnt) =  MS1.ri(i);
%         cMSx(cnt) =  MS1.ci(i);
%         vMSx(cnt) =  MS1.val(i);
%     end
% end
% spMSx = sparse(rMSx,cMSx,vMSx,MS1.m,MS1.n);
% spy(spMSx,30);
