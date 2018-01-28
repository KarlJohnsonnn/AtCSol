function OutPut=ncdf_combustion()

global ls lw
% close all
% source= '/Users/schimmel/schimmel_chemie/NetCDF/ERC_nheptane.nc';OutPut.fuel='nc7h16';
% source= '/Users/schimmel/schimmel_chemie/NetCDF/nheptane.nc';OutPut.fuel='nc7h16';
source= '/Users/schimmel/Code/AtCSol_backup_ISSA_currentVersion/AtCSol/NetCDF/LLNL_MD.nc';OutPut.fuel='md';
% 

% source= '/Users/schimmel/schimmel_chemie/NetCDF/LLNL_MD_2504_hightol.nc';OutPut.fuel='md';


% sumSpcConc = 1.0d0;
sumSpcConc = 0.3207263001d-04;

if (sumSpcConc~=1.0d0)
    ystring = 'Conc in [mole fraction]';
else
    ystring = 'Conc in [mol/cm3]';
end

% ploty = false;
ploty = true;

OutPut.n2='n2';
OutPut.o2='o2';
OutPut.co='co';
OutPut.co2='co2';
OutPut.oh='oh';
OutPut.h2o='h2o';
OutPut.Temp='Temperature';
OutPut.time='time';
stpsze='Step_Size';
OutPut.fueldata = ncread(source,OutPut.fuel)/sumSpcConc;
OutPut.n2data   = ncread(source,OutPut.n2)/sumSpcConc;
OutPut.o2data   = ncread(source,OutPut.o2)/sumSpcConc;
OutPut.co2data  = ncread(source,OutPut.co2)/sumSpcConc;
OutPut.codata   = ncread(source,OutPut.co)/sumSpcConc;
OutPut.ohdata   = ncread(source,OutPut.oh)/sumSpcConc;
OutPut.h2odata  = ncread(source,OutPut.h2o)/sumSpcConc;
OutPut.Tempdata = ncread(source,OutPut.Temp);
OutPut.timedata = ncread(source,OutPut.time);
stpszedata = ncread(source,stpsze);

ls = 15;
lw = 2;

n = 3;
m = 3;
xstring = 'Time in [s]';

% figure;
if (ploty)
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot1(n,m,1,OutPut.timedata,OutPut.fueldata,['fuel: ',OutPut.fuel],xstring,ystring)
    subplot1(n,m,2,OutPut.timedata,OutPut.n2data,'n2',xstring,ystring)
    subplot1(n,m,3,OutPut.timedata,OutPut.o2data,'o2',xstring,ystring)
    subplot1(n,m,4,OutPut.timedata,OutPut.codata,'co',xstring,ystring)
    subplot1(n,m,5,OutPut.timedata,OutPut.ohdata,'oh',xstring,ystring)
    subplot1(n,m,6,OutPut.timedata,OutPut.co2data,'co2',xstring,ystring)
    subplot1(n,m,7,OutPut.timedata,OutPut.h2odata,'h2o',xstring,ystring)
    subplot1(n,m,8,OutPut.timedata,OutPut.Tempdata,'temperature',xstring,'Temp in [K]')
    
    
    subplot(n,m,9);
    semilogy(OutPut.timedata,OutPut.fueldata,'-.','LineWidth',lw);ylim([1.e-12,10]); xlim([OutPut.timedata(1),OutPut.timedata(end)]);hold on;
    semilogy(OutPut.timedata,OutPut.codata,'-.','LineWidth',lw);ylim([1.e-12,10]); xlim([OutPut.timedata(1),OutPut.timedata(end)]);hold on;
    semilogy(OutPut.timedata,OutPut.ohdata,'-.','LineWidth',lw);ylim([1.e-12,10]); xlim([OutPut.timedata(1),OutPut.timedata(end)]);hold on;
    semilogy(OutPut.timedata,OutPut.co2data,'-.','LineWidth',lw);ylim([1.e-12,10]); xlim([OutPut.timedata(1),OutPut.timedata(end)]);hold on;
    semilogy(OutPut.timedata,OutPut.h2odata,'-.','LineWidth',lw);ylim([1.e-12,10]); xlim([OutPut.timedata(1),OutPut.timedata(end)]);
    legend('fuel','co','oh','co2','h2o');
%     figure;
%     plot(OutPut.timedata,stpszedata)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function subplot1(n,m,i,x,y,tit,labx,laby)
%
% global ls lw
%
% subplot(n,m,i);
% plot(x,y,'-.','LineWidth',lw);
% title(tit);
% ylabel(laby, 'FontSize', ls);
% xlabel(labx, 'FontSize', ls);
%
%
% end

function subplot1(n,m,i,x,y,tit,labx,laby)

global ls lw

subplot(n,m,i);
if (i>7)
    plot(x,y,'-.','LineWidth',lw);
    xlim([x(1),x(end)]);
else
    semilogy(x,y,'-.','LineWidth',lw);
    ylim([1.e-12,10]);
    xlim([x(1),x(end)]);
end
title(tit);
ylabel(laby, 'FontSize', ls);
xlabel(labx, 'FontSize', ls);


end
