function Positivity_neu5()
%02.11.2018 Dietz - aktuellste version Version(29.11.2018) 
% Programm positivity_neu, um die AtcSol Matrizen zu vergleichen
close all;
clear all; %#ok<*CLALL>
%%

global Reak Species

% load reactions,species,start and end time:
% [Reak,Species,tBegin,tEnd,dt]=SmallStratoKPP();  %SmallStratoKPP

%  [Reak,Species,tBegin,tEnd,dt]=Robertson();%Fehler?

  [Reak,Species,tBegin,tEnd,dt]=Brusselator();%Brusselator


%Positivitätsabfrage des Reagktionssystems
% positive()

%Anfangswert y0
y=zeros(size(Species));
for i=1:size(Species,1)
    y(Species(i).Pos)=Species(i).c;
end

% Einmalig ausführen für matfile der exakten Loesung:
% erzeugung_matfile_exakte_Loesung(tBegin,tEnd,y)
%  load('test_niedrige_tol_-3')
%load exakte Lösung für relTol=absTol=1e-8, tBegin=43200 und tEnd=36*3600

%%
%ODE Solver
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% options = odeset('RelTol',1e-3,'AbsTol',1e-3);

[tt_tempF,yy_tempF] = ode23s(@odefun,[tBegin,tEnd],y,options);


%Eigene Solver
% RK=RungeKutta('RK2a');
%RK=RungeKutta('RK2b');
RK=RungeKutta('SSP1');
%  RK=RungeKutta('meister3');
%  RK=RungeKutta('meister2');
%  RK=RungeKutta('meister1');

yRK=zeros(size(y,1),4);
yRKExp=zeros(size(y,1),4);
yRKPat=zeros(size(y,1),4);
errRK=zeros(4,1);
errRKExp=zeros(4,1);
errRKPat=zeros(4,1);
err1=zeros(4,1);
err2=zeros(4,1);
err3=zeros(4,1);
h=zeros(4,1);
for i=1:4
    h(i)=dt;
    err1(i)=h(i);
    err2(i)=h(i)^2;
    err3(i)=h(i)^3;
    [tRKExp,yyRKExp]=RungeKuttaExpMethod(tBegin,tEnd,y,dt,RK);
    [tRK,yyRK]=RungeKuttaMethod(tBegin,tEnd,y,dt,RK);
    %   [tRKPat,yyRKPat]=RungeKuttaPat2Method(tBegin,tEnd,y,dt,RK);
    [tRKPat,yyRKPat]=RungeKuttaPat3Method(tBegin,tEnd,y,dt,RK);
    yRK(:,i)=yyRK(end,:);
    yRKExp(:,i)=yyRKExp(end,:);
    yRKPat(:,i)=yyRKPat(end,:);
    errRK(i)=norm(yRK(:,i)-yy_tempF(end,:)');
    errRKExp(i)=norm(yRKExp(:,i)-yy_tempF(end,:)');
    errRKPat(i)=norm(yRKPat(:,i)-yy_tempF(end,:)');
    dt=0.5*dt;
end



figure;
p_1=loglog(h,errRKPat,'-o');hold on;
p_2=loglog(h,err2);hold on;
p_3=loglog(h,err3);hold off;
h=[p_1;p_2;p_3];
legend(h,'errRKPat','err2','err3', 'Location', 'NorthWest');

% Plots
fun_plot(tt_tempF,yy_tempF,Species,y,'odefun1')%ode
fun_plot(tt_temp2,yy_temp2,Species,y,'Verfahren2')

figure
fun_plot(tt_tempF,yy_tempF,Species,y,'odefun1')%ode
fun_plot(tt_temp3,yy_temp3,Species,y,'Verfahren3')
figure
fun_plot(tt_temp,yy_temp,Species,y,'odefun1')%ode
fun_plot(tt_temp4F,yy_temp4F,Species,y,'Verfahren4')


% fun_plot_brussel(tt_temp,yy_temp,tt_temp2,yy_temp2,Species,y,'odefun1','Verfahren2')
% figure
% fun_plot_brussel(tt_temp,yy_temp,tt_temp3,yy_temp3,Species,y,'odefun1','Verfahren3')
% figure
% fun_plot_brussel(tt_temp,yy_temp,tt_temp4,yy_temp4,Species,y,'odefun1','Verfahren4')

%Plot Schrittweite ode23s
figure
plot(0:length(tt_tempF)-1, diff([tBegin; tt_tempF]))
grid


end

%%
function dydt = odefun1(t,y)
global Reak Species
M=M_Berechnung(Reak,Species,t,y);%alt
dydt=M*y;
end

function dydt = odefun(t,y)%löschbar, test für odefun1
global Reak Species
dydt=zeros(size(y));
for i=1:size(Reak,1)
    r=Rate(Reak(i),t);
    for j=1:Reak(i).Left
        Pos=GetPos(Reak(i).NameL(j).Name,Species);
        r=r*y(Pos);
    end
    for j=1:Reak(i).Left
        Pos=GetPos(Reak(i).NameL(j).Name,Species);
        dydt(Pos)=dydt(Pos)-Reak(i).KoeffL(j)*r;
    end
    for j=1:Reak(i).Right
        Pos=GetPos(Reak(i).NameR(j).Name,Species);
        dydt(Pos)=dydt(Pos)+Reak(i).KoeffR(j)*r;
    end
end
%dydt1 = odefun1(t,y)
end

function fun_plot_brussel(tt,yy,tt2,yy2,Species,y,Verfahren,Verfahren2)
%Plot function
for i=1:size(y)
    for j=1:size(y)
        if i==Species(j).Pos
            PerPos=j;
        end
    end
    hold on
    % figure
    subplot(2,3,i);
    plot(tt,yy(:,i),'-o')
    hold on
    % figure
    plot(tt2,yy2(:,i),'--')
    legend( {Verfahren,Verfahren2},'FontSize',22);
    title( {Species(PerPos).Name} ,'FontSize',24);
end
end

function positive()
global Reak Species
%Positivitätsabfrage des Reagktionssystems
sum_left_side=0;
sum_right_side=0;
for i=1:size(Reak,1)
    for j=1:Reak(i).Left
        Pos=GetPosList(Reak(i).NameL(j).Name,Species);
        sum_left_side=sum_left_side+Species(Pos).MolMass*Reak(i).KoeffL(j);
    end
    for j=1:Reak(i).Right
        Pos=GetPosList(Reak(i).NameR(j).Name,Species);
        sum_right_side=sum_right_side+Species(Pos).MolMass*Reak(i).KoeffR(j);
    end
end
if sum_left_side>=sum_right_side && sum_right_side>=0
    disp('The PDS is conservative')% production-destruction systems (PDS)
else
    disp('The PDS is NOT conservative')% production-destruction systems (PDS)
    Error
end
end

function erzeugung_matfile_exakte_Loesung(tBegin,tEnd,y) %#ok<*DEFNU>
%Automatisches erzeugen und speichern der exakten Lösung in 1h Bereichen für
%24Stunden
%  save('StratoExakteLoesung_dt1h_24h')
%  save('StratoExakteLoesung_dt1h_24h_2')
stunde=3600;
nn = (tEnd-tBegin)/3600;
tNeu = tBegin;

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

for ii=1:nn
    tBegin=tNeu;
    tNeu=tNeu+stunde;
    [tt_temp,yy_temp] = ode23s(@odefun,[tBegin,tNeu],y,options);
    
    tm=length(tt_temp);
    exakte_Loesung(ii).tt_temp=tt_temp(tm); %#ok<*AGROW>
    exakte_Loesung(ii).yy_temp=yy_temp(tm,:); %#ok<*AGROW>
    save('test_niedrige_tol_-3','exakte_Loesung')
end
disp('>>>>>>>>>>>>>exakte Loesung erzeugt!')
end










