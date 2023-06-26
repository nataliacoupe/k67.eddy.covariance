% Reads satellite precipitation and temperature datasets
% Natalia Restrepo
%..........................................................................
flag_unix=0;                                            flag_MCWD=1;
fill_ET = 0; %if you want to fill the CWD ET
if flag_unix == 1
    slash    = '\';
    addpath \Users\ncoupe\Documents\OZ\
    addpath \Users\ncoupe\Documents\Amazon\
elseif flag_unix == 0
    FolderTRMMsites    = '/home/ncoupe/Documents/AmazonRS/TRMM/TRMMsites/';
    FolderTRMMdata     = '/media/ncoupe/Backup/Data/AmazonRS/TRMM/trmm_3hr_v7/K67/';
    FolderSoilMoisture = '/home/ncoupe/Documents/AmazonGEP/K67/SoilMoisture/';
    FolderBF     = '/home/ncoupe/Documents/Amazon/AllSites/';
    FolderTRMM   = '/home/ncoupe/Documents/AmazonRS/TRMM/';
    FolderMOD11  = '/home/ncoupe/Documents/AmazonRS/MOD11A2/';
    FolderCERES  = '/home/ncoupe/Documents/AmazonRS/CERES/CERESsites/';
    slash        = '/';
    addpath /home/ncoupe/Documents/OZ/;
    addpath /home/ncoupe/Documents/Amazon/;
    load ([FolderCERES 'siteAMceres.mat']);
    load([FolderTRMM 'siteAMtrmm.mat'])
    load([FolderTRMM 'siteAMgpm.mat'])
end

[K67.monthNetRad,~,~,~,~,K67.month_stgNetRad,K67.month] = AM_month2cycle(K67.netrad,K67.dateloc,tresholdCycle,3,2);
[K67.monthY,K67.monthM,K67.monthD] = datevec(K67.month);
K67.monthNetRad   = AM_rm_outlier(K67.monthNetRad,5);              	K67.monthNetRad = AM_fill(K67.monthNetRad);
[K67.month_avgNetRad,K67.month_avgstdNetRad] = AM_month_avg(K67.monthNetRad,K67.month);

K67.ET = K67.LE./K67.lambda;
K67.tasp_fill = AM_fill(K67.tasp);
K67.tasp_fill = AM_fill_day_cycle(K67.tasp,K67.dateloc,6,18);

K67.lambda = ((2.50*10^6)-((2.386*10^3).*(K67.tasp_fill)))./10^6;                %2.45 [MJ kg-1]     % Latent heat of vaporization, 2.45 [MJ kg-1],
K67.lambda = AM_fill_day_cycle(K67.lambda,K67.dateloc,6,18);
[K67.monthNetRad,~,~,~,~,K67.monthstg_NetRad,K67.month] = AM_month2cycle(K67.netrad,K67.dateloc,tresholdCycle,3,2);
K67.monthNetRad   = AM_rm_outlier(K67.monthNetRad,5);              	K67.monthNetRad = AM_fill(K67.monthNetRad);
K67.ET = K67.LE./K67.lambda;
K67.ET = K67.ET.*60.*60.*24/1000000;      %from s-1 to d-1 and from kg-1 to g-1 (mm)
[K67.weekET]       = AM_week2cycle_rs(K67.ET,K67.dateloc,tresholdCycle,3,2);

K67.tasp_fill = AM_fill(K67.tasp);
K67.tasp_fill = AM_fill_day_cycle(K67.tasp,K67.dateloc,6,18);
[K67.monthTasp] = AM_month2cycle(K67.tasp,K67.dateloc,tresholdCycle,3,2);
[K67.month_avgTasp,K67.month_avgstdTasp] = AM_month_avg(K67.monthTasp,K67.month);
[K67.monthTair4] = AM_month2cycle(K67.tair4,K67.dateloc,tresholdCycle,3,2);

K67.lambda = ((2.50*10^6)-((2.386*10^3).*(K67.tasp_fill)))./10^6;                %2.45 [MJ kg-1]     % Latent heat of vaporization, 2.45 [MJ kg-1],
K67.lambda = AM_fill_day_cycle(K67.lambda,K67.dateloc,6,18);
[K67.monthNetRad,~,~,~,~,K67.monthstg_NetRad,K67.month] = AM_month2cycle(K67.netrad,K67.dateloc,tresholdCycle,3,2);
K67.monthNetRad   = AM_rm_outlier(K67.monthNetRad,5);     K67.monthNetRad = AM_fill(K67.monthNetRad);
K67.ET = K67.LE./K67.lambda;
K67.ET = K67.ET.*60.*60.*24/1000000;      %from s-1 to d-1 and from kg to g (mm)
[K67.weekET]       = AM_week2cycle_rs(K67.ET,K67.dateloc,tresholdCycle,3,2);
[K67.eightET,~,~,~,~,~,K67.eight] = AM_eight2cycle_rs(K67.ET,K67.dateloc,tresholdCycle,3,2);
K67.monthET = AM_month2cycle(K67.ET,K67.dateloc,tresholdCycle,3,2);
K67.monthET = K67.monthET.*eomday(K67.monthY,K67.monthM);
K67.monthLE = AM_month2cycle(K67.LE,K67.dateloc,tresholdCycle,3,2);
[~,K67.monthPrecip,x] = AM_month(K67.rain,K67.dateloc);
K67.monthPrecip(x<(24*20)) = NaN;
[K67.month_avgPrecip] = AM_month_avg(K67.monthPrecip,K67.month);

K67.ETday = K67.ET;              K67.ETday(K67.So<tresholdSo) = NaN;
[K67.dayETday,~,x,~,K67.day,K67.dayETdayMIN,K67.dayETdayMAX]  = AM_day(K67.ETday,K67.dateloc,0);
K67.dayETday(x<6) = NaN;         K67.dayETdayMIN(x<6) = NaN;    K67.dayETdayMAX(x<6) = NaN;
[K67.weekETday,~,~,~,~,~,K67.week] = AM_week2cycle_rs(K67.ETday,K67.dateloc,tresholdCycleDay,3,2);
K67.weekETday     = AM_rm_outlier(K67.weekETday,6);             K67.weekETday  = AM_fill(K67.weekETday);
[K67.week_avgETday,~,K67.week_stdavgETday] = AM_week2_avg(K67.weekETday,K67.week);
[K67.eightETday]  = AM_eight2cycle_rs(K67.ETday,K67.dateloc,tresholdCycleDay,3,2);
K67.eightETday    = AM_rm_outlier(K67.eightETday,6);            K67.eightETday = AM_fill(K67.eightETday);
[K67.monthETday,~,~,~,K67.monthETdayMAX] = AM_month2cycle(K67.ETday,K67.dateloc,tresholdCycleDay,3,2);
K67.monthETday    = AM_rm_outlier(K67.monthETday,6);            K67.monthETday = AM_fill(K67.monthETday);
K67.monthETday    = K67.monthETday.*eomday(K67.monthY,K67.monthM);
K67.monthETdayMAX = AM_rm_outlier(K67.monthETdayMAX,6);         K67.monthETdayMAX = AM_fill(K67.monthETdayMAX);
K67.monthETdayMAX = K67.monthETdayMAX.*eomday(K67.monthY,K67.monthM);

[K67.monthY,K67.monthM,K67.monthD] = datevec(K67.month);
[K67.weekY,K67.weekM,K67.weekD]    = datevec(K67.week);
[K67.eightY,K67.eightM,K67.eightD] = datevec(K67.eight);

%..........................................................................
%% Call TRMM
%..........................................................................
ip = 5;
TRMM.dateloc = TRMMam(ip).date30;       TRMM.precip = TRMMam(ip).mean30;
GPM.dateloc  = GPMam(ip).date30;        GPM.precip  = GPMam(ip).mean30;

% Fill missing TRMM
ind = find(GPM.dateloc==TRMM.dateloc(end));
TRMM.dateloc = [TRMM.dateloc;GPM.dateloc(ind+1:end)];
TRMM.precip  = [TRMM.precip;GPM.precip(ind+1:end)];

% Calculations
[GPM.y,GPM.m,GPM.d] = datevec(GPM.dateloc);
GPM.year = unique(GPM.y);
[GPM.month_avgPrecip,GPM.month_avgPrecipstd] = AM_month_avg(GPM.precip,GPM.dateloc);
K67.monthGPM = interp1(GPM.dateloc,GPM.precip,K67.month);
K67.weekGPM = interp1(K67.month,K67.monthGPM,K67.week);

g = length(unique(GPM.y));
GPM.month_anomalyPrecip   = NaN((g*12),1);
GPM.REPmonth_avgPrecip    = NaN((g*12),1);     GPM.REPmonth_avgPrecipstd = NaN((g*12),1);
for ik = 1:g*12
    iq = GPM.m(ik);
    A  = GPM.month_avgPrecip(iq);
    B  = GPM.month_avgPrecipstd(iq);
    C  = GPM.precip(ik);
    D  = (C-A)./B;
    D(B==0) = 0;
    GPM.month_anomalyPrecip(ik)   = D;
    GPM.REPmonth_avgPrecip(ik)    = A;
    GPM.REPmonth_avgPrecipstd(ik) = B;
end

%% ........................................................................
[TRMM.y,TRMM.m,TRMM.d] = datevec(TRMM.dateloc);
TRMM.year = unique(TRMM.y);
[TRMM.month_avgPrecip,TRMM.month_avgPrecipstd] = AM_month_avg(TRMM.precip,TRMM.dateloc);
K67.monthTRMM = interp1(TRMM.dateloc,TRMM.precip,K67.month);
g = length(unique(TRMM.y));
TRMM.month_anomalyPrecip   = NaN((g*12),1);
TRMM.REPmonth_avgPrecip    = NaN((g*12),1);     TRMM.REPmonth_avgPrecipstd = NaN((g*12),1);
for ik = 1:g*12
    iq = TRMM.m(ik);
    A  = TRMM.month_avgPrecip(iq);
    B  = TRMM.month_avgPrecipstd(iq);
    C  = TRMM.precip(ik);
    D  = (C-A)./B;
    D(B==0) = 0;
    TRMM.month_anomalyPrecip(ik)   = D;
    TRMM.REPmonth_avgPrecip(ik)    = A;
    TRMM.REPmonth_avgPrecipstd(ik) = B;
end

if flag_MCWD == 1
    ip = 5;
    %ET from CERES when missing
    monthCERESnr   = CERESam(ip).mean30nr;        %monthCERESswd = [monthCERESswd(2:end);monthCERESswd(1)];
    monthCERESdate = CERESam(ip).date30;
    monthCERESswd  = CERESam(ip).mean30swd;
    %     monthCERESswd     = CERESam(ip).NetRadmean30;        %monthCERESswd = [monthCERESswd(2:end);monthCERESswd(1)];
    
    % %     [CERES.y,CERES.m,CERES.d] = datevec(monthCERESdate);
    % %     tableK67ceres = table(CERES.y,CERES.m,CERES.d,monthCERESnr, monthCERESswd);
    % %     tableK67.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'NetRadiationCERES_Wm2' 'SWdownCERES_Wm2'};
    % %     writetable (tableK67ceres,'out/K67.ceres.dat');
    
    [CERES.month_avgSWdown,CERES.month_avgSWdownstd] = AM_month_avg(monthCERESswd,monthCERESdate);
    [CERES.month_avgNR,CERES.month_avgNRstd]         = AM_month_avg(monthCERESnr,monthCERESdate);
    [CERES.y,CERES.m,~] = datevec(monthCERESdate);
    
    g = length(unique(CERES.y));
    CERES.month_anomalySWdown   = NaN.*CERES.m;
    CERES.REPmonth_avgSWdown    = NaN.*CERES.m;     CERES.REPmonth_avgSWdownstd = NaN.*CERES.m;
    for ik = 1:length(CERES.m)
        iq = CERES.m(ik);
        A  = CERES.month_avgSWdown(iq);
        B  = CERES.month_avgSWdownstd(iq);
        C  = monthCERESswd(ik);
        D  = (C-A)./B;
        D(B==0) = 0;
        CERES.month_anomalySWdown(ik)   = D;
        CERES.REPmonth_avgSWdown(ik)    = A;
        CERES.REPmonth_avgSWdownstd(ik) = B;
    end
    
    K67.monthCERESnr  = interp1(monthCERESdate,monthCERESnr,K67.month);
    K67.monthCERESswd = interp1(monthCERESdate,monthCERESswd,K67.month);
    K67.weekCERESnr  = interp1(monthCERESdate,monthCERESnr,K67.week);
    K67.weekCERESswd = interp1(monthCERESdate,monthCERESswd,K67.week);
    
    [K67.week_avgCERESswd,K67.week_avgstdCERESswd] = AM_week2_avg(K67.weekCERESswd,K67.week);
    [K67.week_avgCERESnr,K67.week_avgstdCERESnr] = AM_week2_avg(K67.weekCERESnr,K67.week);
    
    %                sm   =    standard deviation of the slope
    %%               sb   =    standard deviation of the y-intercept
    [p,r,~,p_val,sm,sb] = AM_compute_regression(K67.monthCERESnr,K67.monthNetRad,[.6 .6 .6],0);
    figure('color','white');
    subplot(2,2,1:2);  plot(K67.month,K67.monthNetRad);
    plot(K67.month,K67.monthNetRad);
    hold on;           plot(K67.month,K67.monthCERESnr);   plot(K67.month,K67.monthCERESswd);
    %     plot(K67.month,[NaN;K67.monthCERESnr(1:end-1)],'k');
    datetick('x');     ylabel(['NetRad (W m^-^2)']);
    legend('Rn_{EC}','NetRad_{CERES} K67','SW_{CERES}');
    subplot(2,2,3);
    [pRn,r2,~,p_val,sm,sb,K67.RNmodelET] = AM_compute_regression(K67.monthCERESnr,K67.monthET,[.6 .6 .6],1);
    ylabel([' ET (mm month^-^1)']);                           xlabel('NetRad_{CERES}(W m^-^2)');
    text(102,140,['ET =',num2str(pRn(1,1),3.2), 'NR_{CERES} ', num2str(pRn(1,2),2.1)])
    text(102,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(102,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    subplot(2,2,4);    [p,r2] = AM_compute_regression(K67.monthNetRad,K67.monthET,[.6 .6 .6],1);
    ylabel([' ET (mm month^-^1)']);                           xlabel('NetRad (W m^-^2)');
    text(102,140,['ET =',num2str(p(1,1),3.2), 'Rn ', num2str(p(1,2),2.1)])
    text(102,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(102,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    K67.monthETceres =  polyval(p,monthCERESnr);
    %%
    monthCERESnr      = CERESam(ip).mean30nr;          % monthCERESswd = [monthCERESswd(2:end);monthCERESswd(1)];
    K67.monthCERESnr  = interp1(monthCERESdate,monthCERESnr,K67.month);
    TRMM.monthCERESnr = interp1(monthCERESdate,monthCERESnr,TRMM.dateloc);
    monthCERESswd     = CERESam(ip).mean30swd;                % monthCERESswd = [monthCERESswd(2:end);monthCERESswd(1)];
    %     monthCERESswd           = CERESam(ip).NetRadmean30;        % fake as netrad
    K67.monthCERESswd  = interp1(monthCERESdate,monthCERESswd,K67.month);
    K67.monthCERESnr   = interp1(monthCERESdate,monthCERESnr,K67.month);
    TRMM.monthCERESswd = interp1(monthCERESdate,monthCERESswd,TRMM.dateloc);
    TRMM.monthCERESnr  = interp1(monthCERESdate,monthCERESnr,TRMM.dateloc);
    
    [K67.month_avgCERESnr,K67.month_avgstdCERESnr]   = AM_month_avg(K67.monthCERESnr,K67.month);
    [K67.month_avgCERESswd,K67.month_avgstdCERESswd] = AM_month_avg(K67.monthCERESswd,K67.month);
    [K67.month_avgCERESnr,K67.month_avgstdCERESnr]   = AM_month_avg(K67.monthCERESnr,K67.month);
    
    subplot(2,2,4);
    [pSW,r2,~,p_val,sm,sb,K67.SWmodelET] = AM_compute_regression(K67.monthCERESswd,K67.monthET,[.6 .6 .6],1);
    ylabel(['ET (mm month^-^1)']);                           xlabel('SW_{down_{CERES}}(W m^-^2)');
    text(150,140,['ET =',num2str(pSW(1,1),3.2), 'SW_{CERES} ', num2str(pSW(1,2),2.1)])
    text(150,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(150,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    %%
    figure('color','white');
    subplot(2,2,1);
    [p,r2] = AM_compute_regression(K67.monthCERESswd,K67.monthNetRad,[.6 .6 .6],1);
    xlabel(['SW{down}_{CERES}(W m^-^2)']);                    ylabel('NetRad (W m^-^2)');
    text(102,140,['ET =',num2str(pRn(1,1),3.2), 'SW_{down} ', num2str(pRn(1,2),2.1)])
    text(102,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(102,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    subplot(2,2,2);
    [p,r2] = AM_compute_regression(K67.monthCERESnr,K67.monthNetRad,[.6 .6 .6],1);
    xlabel([' NR_{CERES} (W m^-^2)']);                           ylabel('NetRad (W m^-^2)');
    text(150,140,['NetRad =',num2str(p(1,1),3.2), 'NR_{CERES} ', num2str(p(1,2),2.1)])
    text(150,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(150,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    subplot(2,2,3);
    [p,r2] = AM_compute_regression([NaN;K67.monthCERESnr(1:end-1)],K67.monthNetRad,[.6 .6 .6],1);
    xlabel([' NR_{CERES} (W m^-^2)']);                           ylabel('NetRad (W m^-^2)');
    text(150,140,['NetRad =',num2str(p(1,1),3.2), 'NR_{CERES} ', num2str(p(1,2),2.1)])
    text(150,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(150,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    subplot(2,2,4);
    [pSW,r2,~,p_val,sm,sb,K67.SWmodelET] = AM_compute_regression(K67.monthCERESswd,K67.monthET,[.6 .6 .6],1);
    ylabel(['ET (mm month^-^1)']);                           xlabel('SW_{down_{CERES}}(W m^-^2)');
    text(250,140,['ET =',num2str(pSW(1,1),3.2), 'SW_{CERES} ', num2str(pSW(1,2),2.1)])
    text(250,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
    text(250,120,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
    
    %% Fill ET K67 time series...............................................
    K67.monthETfill = K67.monthET;
    K67.monthETfill(isnan(K67.monthETfill)) = K67.SWmodelET(isnan(K67.monthETfill));
    [~,ix]          = AM_fill_month_cycle(K67.monthET,K67.month,6,18);
    K67.monthETfill(isnan(K67.monthETfill)) = ix(isnan(K67.monthETfill));
    %% CWD TRMM time series..................................................
    % EC ET on the TRMM time series
    ix = interp1(K67.month,K67.monthET,TRMM.dateloc);
    % EC ETceres on the TRMM time series
    % iy = interp1(monthCERESdate,K67.SWmodelET,TRMM.dateloc);
    TRMM.SWmodelETfill = polyval(pSW,TRMM.monthCERESswd);
    TRMM.RNmodelETfill = polyval(pRn,TRMM.monthCERESnr);
    % filled ET from EC ET monthly cycle
    iz  = AM_fill_month_cycle(ix,TRMM.dateloc,6,18);
    TRMM.monthET       = ix;
    TRMM.monthETfill   = ix;
    
    TRMM.monthETfill(isnan(ix))               = TRMM.SWmodelETfill(isnan(ix));
    TRMM.monthETfill(isnan(TRMM.monthETfill)) = iz(isnan(TRMM.monthETfill));
    
    figure('color','white');
    subplot(2,2,1:2);       plot(TRMM.dateloc,iz,'color',[0 0 .9],'LineWidth',2);
    hold on;                plot(TRMM.dateloc,TRMM.SWmodelETfill,'color',[.5 0 .9],'LineWidth',2);
    plot(TRMM.dateloc,ix,'color',[.9 0 0],'LineWidth',2);
    plot(TRMM.dateloc,TRMM.monthETfill,'color',[0 0 0],'LineWidth',1);
    datetick('x');  ylabel({' ET (mm month^-^1)'}); legend('Monthly EC','SW_{CERES}','EC','Filled');
    % .....................................................................
    %% CWD
    % .....................................................................
    TRMM.monthCWD = zeros(length(TRMM.dateloc),1);
    for im=2:length(TRMM.dateloc)
        if (TRMM.monthCWD(im-1)-TRMM.monthETfill(im)+TRMM.precip(im)) < 0
            TRMM.monthCWD(im) = TRMM.monthCWD(im-1)-TRMM.monthETfill(im)+TRMM.precip(im);
        else
            TRMM.monthCWD(im) = 0;
        end
    end
    
    TRMM.monthCWDceres = zeros(length(TRMM.dateloc),1);
    for im=2:length(TRMM.dateloc)
        if (TRMM.monthCWDceres(im-1)-TRMM.SWmodelETfill(im)+TRMM.precip(im)) < 0
            TRMM.monthCWDceres(im) = TRMM.monthCWDceres(im-1)-TRMM.SWmodelETfill(im)+TRMM.precip(im);
        else
            TRMM.monthCWDceres(im) = 0;
        end
    end
    
    TRMM.monthCWD100 = zeros(length(TRMM.dateloc),1);
    for im=2:length(TRMM.dateloc)
        if (TRMM.monthCWD100(im-1)-100+TRMM.precip(im)) < 0
            TRMM.monthCWD100(im) = TRMM.monthCWD100(im-1)-100+TRMM.precip(im);
        else
            TRMM.monthCWD100(im) = 0;
        end
    end
end

% .....................................................................
%% CWD as positive values
% .....................................................................
TRMM.monthCWDceres = -TRMM.monthCWDceres;
TRMM.monthCWD100   = -TRMM.monthCWD100;
TRMM.monthCWD      = -TRMM.monthCWD;

K67.monthCDW = interp1(TRMM.dateloc,TRMM.monthCWD,K67.month);
K67.month_avgCDW = AM_month_avg(K67.monthCDW,K67.month);
K67.monthTRMManomaly = interp1(TRMM.dateloc,TRMM.month_anomalyPrecip,K67.month);

%% ........................................................................
[TRMM.month_avgCWD,TRMM.month_avgCWDstd] = AM_month_avg(TRMM.monthCWD,TRMM.dateloc);

TRMM.precip_CWD = TRMM.precip./TRMM.monthCWD;
TRMM.monthCWDanomaly = NaN((g*12),1);
TRMM.REPmonth_avgCWD = NaN((g*12),1);   TRMM.REPmonth_avgCWDstd = NaN((g*12),1);
for ik = 1:length(TRMM.dateloc)
    iq = TRMM.m(ik);
    A  = TRMM.month_avgCWD(iq);
    B  = TRMM.month_avgCWDstd(iq);
    C  = TRMM.monthCWD(ik);
    D  = (C-A)./B;
    D(B==0) = 0;
    TRMM.monthCWDanomaly(ik)    = D;
    TRMM.REPmonth_avgCWD(ik)    = A;
    TRMM.REPmonth_avgCWDstd(ik) = B;
end

% .....................................................................
%% Open the text file.
% .....................................................................
FileName   = [FolderTRMM slash 'multivariateENSOindex.txt'];
delimiter  = '\t';
startRow   = 4;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(FileName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
MEI.year       = dataArray{1,1};
MEI.year_month = [];            MEI.value_month = [];      MEI.month = [];
MEI.trim       = [];            MEI.dateloc     = [];
for ik=1:length(MEI.year)
    for ip=2:13
        MEI.value_month = [MEI.value_month;dataArray{1,ip}(ik)];
        MEI.month       = [MEI.month;ip-1];
        MEI.trim        = [MEI.trim;ceil((ip-1)/3)];
        MEI.year_month  = [MEI.year_month;MEI.year(ik)];
        MEI.dateloc     = [MEI.dateloc;datenum(MEI.year(ik),(ip-1),1)];
    end
end
MEI.value_month(MEI.value_month==-999.)=NaN;
[MEI.y,MEI.m,MEI.d]=datevec(MEI.dateloc);


TRMM.mei = interp1(MEI.dateloc,MEI.value_month,TRMM.dateloc);
GPM.mei = interp1(MEI.dateloc,MEI.value_month,GPM.dateloc);

[TRMM.month_avgMEI,TRMM.month_avgstdMEI] = AM_month_avg(TRMM.mei,TRMM.dateloc);
% % %% Open the text file.
% % FileName   = [FolderTRMM slash 'multivariateENSOindexRANK.txt'];
% % delimiter  = '\t';
% % startRow   = 3;
% % formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
% %
% % fileID = fopen(FileName,'r');
% % dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% % fclose(fileID);
% % MEI.year       = dataArray{1,1};
% % MEI.year_month = [];            MEI.rank_month = [];      MEI.month = [];
% % MEI.trim       = [];            MEI.dateloc     = [];
% % for ik=1:length(MEI.year)
% %     for ip=2:13
% %         MEI.rank_month = [MEI.rank_month;dataArray{1,ip}(ik)];
% %         MEI.month       = [MEI.month;ip-1];
% %         MEI.trim        = [MEI.trim;ceil((ip-1)/3)];
% %         MEI.year_month  = [MEI.year_month;MEI.year(ik)];
% %         MEI.dateloc     = [MEI.dateloc;datenum(MEI.year(ik),(ip-1),1)];
% %     end
% % end
% % MEI.rank_month(MEI.rank_month==9999)=NaN;
% % [MEI.y,MEI.m,MEI.d]=datevec(MEI.dateloc);
% %
% .........................................................................
%% Check MOD11 quick
% Import TairMOD11.csv
% .........................................................................
ProductName   = {'LST_Day_1km';'LST_Night_1km';'Emis_31';'Emis_32'};
load([FolderMOD11 'siteMOD11A2lst.mat'],'siteMOD11A2lst');
ip = 5;          ia = 1;
% lst_timeseries = AM_rm_outlier(siteMOD11A2lst(ip).product(ia).mean08,3);
MOD11.LSTnight = siteMOD11A2lst(ip).product(2).mean08-273.15;
MOD11.LSTnight(MOD11.LSTnight<20) = NaN;
MOD11.LSTday   = siteMOD11A2lst(ip).product(1).mean08-273.15;
MOD11.LSTday((MOD11.LSTday<10)|(MOD11.LSTday>40))     = NaN;
MOD11.tair     = (AM_fill_plus(MOD11.LSTnight)+AM_fill_plus(MOD11.LSTday))/2;
MOD11.tair((MOD11.tair>40)|(MOD11.tair<10))= NaN;
MOD11.dateloc  = siteMOD11A2lst(ip).product(ia).date08;

%%
figure('color','white');                    plot(MOD11.dateloc,MOD11.tair);    hold on;
text(MOD11.dateloc(3),0.01,[site_name_BF{ip} '  ' ProductName{ia}])
% MOD11.tair = AM_rm_outlier(siteMOD11A2lst(ip).product(ia).mean08-273.15,3);
% [~,MOD11.tair] = AM_spike_detection(MOD11.tair,3);
MOD11.tair     = AM_fill_plus(MOD11.tair);
plot(MOD11.dateloc,MOD11.tair);        datetick('x');

K67.monthTairMODIS = interp1(MOD11.dateloc,MOD11.tair,K67.month);

ix = smooth(MOD11.tair,5);                  ix((ix>40)|(ix<10))= NaN;
ix(isnan(MOD11.tair)) = NaN;                MOD11.tair_sm = ix;
%%
figure('color','white');
hl2 = plot(MOD11.dateloc,MOD11.LSTnight+3,'g');       hold on;
% hl1 = plot(MOD11.dateloc,MOD11.LSTday,'c');       hold on;
hl3 = plot(MOD11.dateloc,MOD11.tair_sm,'LineWidth',2);
hl4 = plot(K67.month,K67.monthTasp,'k');
hl5 = plot(K67.month,K67.monthTair4,'color',[.6 .6 .6]);
hl1 = plot(MOD11.dateloc,MOD11.tair+1,'b');       hold on;
datetick('x');        grid on;      ylabel ('T_{air}(\circC)');
legend([hl1,hl2,hl3,hl4,hl5],'MOD11','MOD11night','MOD11_{mw}','Tair_{EC:50m}','T_{lv:29m}');

K67.weekLSTnight = interp1(MOD11.dateloc,MOD11.LSTnight+3,K67.week);
K67.weekLSTnight = AM_fill_plus2(K67.weekLSTnight);
K67.weekLSTnight = AM_fill_plus2(K67.weekLSTnight);

K67.eightLST = interp1(MOD11.dateloc,MOD11.tair,K67.eight);
K67.weekLST = interp1(MOD11.dateloc,MOD11.tair,K67.week);
K67.weekLST = AM_fill_plus2(K67.weekLST);
K67.weekLST = AM_fill_plus2(K67.weekLST);

%..........................................................................
%%   Linear regression models
%..........................................................................
ind = find (K67.monthY<=2022);
tbl = table(K67.monthCERESnr(ind),...
    K67.monthTairMODIS(ind),K67.monthLE(ind),'VariableNames',{'NetRad','Tair','LE'});
lm  = fitlm(tbl,'LE~NetRad+Tair')
AICvalue1 = lm.ModelCriterion.AIC

LE_model = lm.Coefficients.Estimate(1)+...
    (lm.Coefficients.Estimate(2).*K67.monthCERESswd)+...
    (lm.Coefficients.Estimate(3).*K67.monthTairMODIS);

ind = find((~isnan(K67.monthLE))&(~isnan(LE_model)));
display 'R^2 LE_{model} vs. LE_{observations} ' ; corrcoef(LE_model(ind),K67.monthLE(ind))
display 'diff LE_{model} - LE_{observations} ' ; nanmean(LE_model(ind)-K67.monthLE(ind))

% .........................................................................
ind = find (K67.monthY<2018);
tbl = table(K67.monthCERESnr(ind),K67.monthLE(ind),'VariableNames',{'NetRad','LE'});
lm  = fitlm(tbl,'LE~NetRad')
AICvalue1 = lm.ModelCriterion.AIC

LE_model = lm.Coefficients.Estimate(1)+...
    (lm.Coefficients.Estimate(2).*K67.monthCERESswd);

ind = find((~isnan(K67.monthLE))&(~isnan(LE_model)));
display 'R^2 LE_{model} vs. LE_{observations} ' ; corrcoef(LE_model(ind),K67.monthLE(ind))
display 'diff LE_{model} - LE_{observations} ' ;nanmean(LE_model(ind)-K67.monthLE(ind))

% .....................................................................
ind = find (K67.monthY<2018);
tbl = table(K67.monthTairMODIS(ind),K67.monthLE(ind),'VariableNames',{'Tair','LE'});
lm  = fitlm(tbl,'LE~Tair')
AICvalue1 = lm.ModelCriterion.AIC

LE_model = lm.Coefficients.Estimate(1)+...
    (lm.Coefficients.Estimate(2).*K67.monthTairMODIS);

ind = find((~isnan(K67.monthLE))&(~isnan(LE_model)));
display 'R^2 LE_{model} vs. LE_{observations} ' ; corrcoef(LE_model(ind),K67.monthLE(ind))
display 'diff LE_{model} - LE_{observations} ' ;nanmean(LE_model(ind)-K67.monthLE(ind))

%..........................................................................
%% Rainfall 3 hour
%..........................................................................
FileName = [FolderTRMMdata  'K67_3hrTRMM.dat.gz'];
gunzip(FileName)
[TRMM3h.y,TRMM3h.m,TRMM3h.d,TRMM3h.hr,TRMM3h.precipHQ,TRMM3h.precip,TRMM3h.precipError] =...
    textread([FolderTRMMdata  'K67_3hrTRMM.dat'],...
    '%f %f %f %f %f %f %f ','delimiter', ',','headerlines',1);
% gzip([FolderTRMMsites  'K67_3hrTRMM.dat.gzip'],[FolderTRMMsites  'K67_3hrTRMM.dat'])
delete ('K67_3hrTRMM.dat');
TRMM3h.dateloc = datenum(TRMM3h.y,TRMM3h.m,TRMM3h.d,TRMM3h.hr-utm,0,0);

TRMM1h.dateloc  = (datenum(TRMM3h.dateloc(1):datenum(0,0,0,1,0,0):K67.dateloc(end)))';
TRMM1h.precipHQ = interp1(TRMM3h.dateloc, TRMM3h.precipHQ, TRMM1h.dateloc);
TRMM1h.precip = interp1(TRMM3h.dateloc, TRMM3h.precip, TRMM1h.dateloc);
K67.TRMM1h = interp1(TRMM3h.dateloc, TRMM3h.precip, K67.dateloc);       %K67.TRMM1h = K67.TRMM1h./3;

TRMM1h.rain   = interp1(K67.dateloc(~isnan(K67.rain)), K67.rain(~isnan(K67.rain)), TRMM1h.dateloc);
TRMM1h.precip(isnan(TRMM1h.precip)) = TRMM1h.rain(isnan(TRMM1h.precip));
TRMM1h.precip(isnan(TRMM1h.precip)) = 0;

[~,TRMM1h.dayprecip,~,~,TRMM1h.day] = AM_day(TRMM3h.precip,TRMM3h.dateloc,0);
TRMM1h.dayprecip = TRMM1h.dayprecip.*3;

figure('color','white');
bar(K67.dateloc, K67.TRMM1h); hold on; bar(K67.dateloc,K67.rain);  datetick('x')

figure('color','white');
bar(K67.jdstartgmt, K67.TRMM1h); hold on; bar(K67.jdstartgmt,K67.rain);

[p,r2,~,p_val,sm,sb,K67.hourSWmodelET] = AM_compute_regression(K67.netrad,K67.ET,[.6 .6 .6],0);

K67.ETfill = K67.ET;
K67.ETfill(isnan(K67.ETfill)) = K67.hourSWmodelET(isnan(K67.ETfill));
[~,K67.ETfill] = AM_fill_month_cycle_force(K67.ET,K67.dateloc,6,18);

TRMM1h.ETfill = interp1(K67.dateloc(~isnan(K67.ETfill)), K67.ETfill(~isnan(K67.ETfill)), TRMM1h.dateloc);
TRMM1h.ETfill = AM_fill_month_cycle_force(TRMM1h.ETfill,TRMM1h.dateloc,6,18);
TRMM1h.dayETfill = AM_day(TRMM1h.ETfill,TRMM1h.dateloc,0);

%% CWD TRMM time series................................................
TRMM1h.CWD = zeros(length(TRMM1h.dateloc),1);
for im=2:length(TRMM1h.dateloc)
    if ((TRMM1h.CWD(im-1)-(TRMM1h.ETfill(im)/24)+TRMM1h.precip(im)) < 0)
        TRMM1h.CWD(im) = TRMM1h.CWD(im-1)-(TRMM1h.ETfill(im)/24)+TRMM1h.precip(im);
    else
        TRMM1h.CWD(im) = 0;
    end
end
TRMM1h.CWD       = -TRMM1h.CWD;
K67.TRMM1hCWD    = interp1(TRMM1h.dateloc,TRMM1h.CWD,K67.dateloc);
K67.TRMM1hPrecip = interp1(TRMM1h.dateloc,TRMM1h.precip,K67.dateloc);

%% CWD TRMM time series................................................
TRMM1h.dayCWD = zeros(length(TRMM1h.day),1);
for im=2:length(TRMM1h.day)
    if ((TRMM1h.dayCWD(im-1)-TRMM1h.dayETfill(im)+TRMM1h.dayprecip(im)) < 0)
        TRMM1h.dayCWD(im) = TRMM1h.dayCWD(im-1)-TRMM1h.dayETfill(im)+TRMM1h.dayprecip(im);
    else
        TRMM1h.dayCWD(im) = 0;
    end
end
TRMM1h.dayCWD    = -TRMM1h.dayCWD;
K67.dayTRMM1hCWD = interp1(TRMM1h.day,TRMM1h.dayCWD,K67.day);
K67.dayTRMM1h    = interp1(TRMM1h.day,TRMM1h.dayprecip,K67.day);

[~,~,TRMM1h.weekprecip,~,~,~,TRMM1h.week] = AM_week2day_rs(TRMM3h.precip,TRMM3h.dateloc);
TRMM1h.weekprecip = TRMM1h.weekprecip.*3;
[TRMM1h.weekCWD,~,~,~,~,~,TRMM1h.week] = AM_week2day_rs(TRMM1h.dayCWD,TRMM1h.day);
K67.TRMM1hCWD = interp1(TRMM1h.day,TRMM1h.dayCWD,K67.dateloc);
K67.weekTRMM1hCWD = interp1(TRMM1h.week,TRMM1h.weekCWD,K67.week);

K67.weekTRMM = interp1(TRMM.dateloc,TRMM.precip,K67.week);
K67.weekTRMM1h = interp1(TRMM1h.week,TRMM1h.weekprecip,K67.week);
K67.weekTRMM1h(isnan(K67.weekTRMM1h)) = K67.weekTRMM(isnan(K67.weekTRMM1h));
%..........................................................................
%% Rainfall
%..........................................................................
figure('color','white');
subplot(3,3,1:3);   hl1 = plot(K67.month,K67.monthTRMM);     hold on;
hl2 = plot(K67.month,K67.monthPrecip,'k');
[iy,~,~,ix] = AM_month(TRMM1h.precip,TRMM1h.dateloc);
hl3 = plot(ix,iy.*30.*24);
datetick('x');        grid on;      ylabel ('Precip(mm)');
legend([hl1,hl2,hl3],'TRMM','Precip_{mm}','TRMM_{3hour}');

subplot(3,3,4:6);   hl1 = plot(K67.month,K67.monthCDW);     hold on;
hl2 = plot(TRMM1h.dateloc,TRMM1h.CWD,'k');
[iy,~,~,ix] = AM_month(TRMM1h.CWD,TRMM1h.dateloc);
hl3 = plot(ix,iy);
datetick('x');        grid on;      ylabel ('CWD(mm)');
legend([hl1,hl2,hl3],'CWD_{TRMM month}','CWD_{TRMM hour}','CWD_{TRMM hour2month}');
text (K67.month(3),300,'K67');

subplot(3,3,7);   hl1 = plot(K67.monthPrecip(K67.monthPrecip>0),K67.monthTRMM(K67.monthPrecip>0),'.');
grid on;          ylabel ('TRMM(mm)');         xlabel ('Precip(mm)');
ind = find((K67.monthPrecip>0)&(~isnan(K67.monthPrecip))&(~isnan(K67.monthTRMM)));
disp('R^2 TRMM precipitation obs'); corrcoef(K67.monthPrecip(ind),K67.monthTRMM(ind))
%%
figure('color','white');    bar(K67.month,K67.monthTRMM);      grid on;
datetick('x');              ylabel({'Precipitation_{TRMM}';'(mm month^-^1)'});

% .........................................................................
%% Soil moisture
% .........................................................................
opts = delimitedTextImportOptions("NumVariables", 22);
opts.DataLines = [2, Inf];              opts.Delimiter = ",";
opts.VariableNames = ["hour", "day", "month", "year", "VWC05L", "VWC10L", "VWC20L", "VWC30L", "VWC40L", "VWC70L", "VWC100L", "VWC05R", "VWC10R", "VWC20R", "VWC30R", "VWC40R", "VWC50R", "VWC70R", "VWC100R", "VWCSrf_1", "VWCSrf_2", "VWCSrf_3"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";       opts.EmptyLineRule = "read";
pit1VWC = readtable([FolderSoilMoisture 'pit1VWC.c60.noNA.csv'], opts);
clear opts

pit1VWCdateloc = datenum(pit1VWC.year,pit1VWC.month,pit1VWC.day,pit1VWC.hour,0,0);
pit1VWC05L = pit1VWC.VWC05L;        pit1VWC10L = pit1VWC.VWC10L;    pit1VWC20L = pit1VWC.VWC20L;    pit1VWC30L = pit1VWC.VWC30L;
pit1VWC40L = pit1VWC.VWC40L;        pit1VWC70L = pit1VWC.VWC70L;    pit1VWC100L = pit1VWC.VWC100L;
pit1VWC05R = pit1VWC.VWC05R;        pit1VWC10R = pit1VWC.VWC10R;    pit1VWC20R = pit1VWC.VWC20R;
pit1VWC30R = pit1VWC.VWC30R;        pit1VWC40R = pit1VWC.VWC40R;    pit1VWC50R = pit1VWC.VWC50R;
pit1VWC70R = pit1VWC.VWC70R;        pit1VWC100R = pit1VWC.VWC100R;  pit1VWCSrf_1 = pit1VWC.VWCSrf_1;
pit1VWCSrf_2 = pit1VWC.VWCSrf_2;    pit1VWCSrf_3 = pit1VWC.VWCSrf_3;
pit1VWC05L(pit1VWC05L==-9999.) = NaN;   pit1VWC10L(pit1VWC10L==-9999.) = NaN;   pit1VWC20L(pit1VWC20L==-9999.) = NaN;
pit1VWC30L(pit1VWC30L==-9999.) = NaN;   pit1VWC40L(pit1VWC40L==-9999.) = NaN;   pit1VWC70L(pit1VWC70L==-9999.) = NaN;
pit1VWC100L(pit1VWC100L==-9999.) = NaN;
pit1VWC05R(pit1VWC05R==-9999.) = NaN;   pit1VWC10R(pit1VWC10R==-9999.) = NaN;   pit1VWC20R(pit1VWC20R==-9999.) = NaN;
pit1VWC30R(pit1VWC30R==-9999.) = NaN;   pit1VWC40R(pit1VWC40R==-9999.) = NaN;   pit1VWC70R(pit1VWC70R==-9999.) = NaN;
pit1VWC100R(pit1VWC100R==-9999.) = NaN;

pit1VWC10L = AM_fill_plus(pit1VWC10L);  pit1VWC10L = AM_fill_plus(pit1VWC10L);

% .........................................................................
%% Quick QAQC
% .........................................................................
ind = find(((pit1VWCdateloc>736400)&(pit1VWCdateloc<736500)));
pit1VWC05L(ind) = NaN;  pit1VWC10L(ind) = NaN;

ind = find((pit1VWCdateloc>=734585)&(pit1VWCdateloc<=734588));
pit1VWC20L(ind) = NaN;    pit1VWC30L(ind) = NaN;        pit1VWC40L(ind) = NaN;
pit1VWC70L(ind) = NaN;

ind = find((pit1VWCdateloc>=736400)&(pit1VWCdateloc<=736650));
pit1VWC20L(ind) = NaN;

ind = find((pit1VWCdateloc>=733640)&(pit1VWCdateloc<=733760));
pit1VWC70L(ind) = NaN;

ind = find((pit1VWCdateloc>=736500)&(pit1VWCdateloc<=736650));
pit1VWC100L(ind) = NaN;

ind = find((pit1VWCdateloc>=736430)&(pit1VWCdateloc<=736650));
pit1VWC05R(ind) = NaN;

ind = find((pit1VWCdateloc>=734654)&(pit1VWCdateloc<=734682));
pit1VWC10R(ind) = NaN;

ind = find(((pit1VWCdateloc>734655)&(pit1VWCdateloc<734685))|...
    ((pit1VWCdateloc>734780)&(pit1VWCdateloc<734815))|...
    ((pit1VWCdateloc>734540)&(pit1VWCdateloc<734550))|...
    ((pit1VWCdateloc>734699)&(pit1VWCdateloc<734710))|...
    ((pit1VWCdateloc>736400)&(pit1VWCdateloc<736700)));
pit1VWC10R(ind) = NaN;          pit1VWC20R(ind) = NaN;      pit1VWC30R(ind) = NaN;
pit1VWC40R(ind) = NaN;

ind = find(((pit1VWCdateloc>734585)&(pit1VWCdateloc<734587)));
pit1VWC20R(ind) = NaN;      pit1VWC30R(ind) = NaN;      pit1VWC70R(ind) = NaN;

ind = find(((pit1VWCdateloc>734363)&(pit1VWCdateloc<734480)));
pit1VWC30R(ind) = NaN;          pit1VWC40R(ind) = NaN;      pit1VWC70R(ind) = NaN;

ind = find(((pit1VWCdateloc>734492)&(pit1VWCdateloc<734494))|...
    ((pit1VWCdateloc>734800)&(pit1VWCdateloc<734860)));
pit1VWC40R(ind) = NaN;      pit1VWC70R(ind) = NaN;

ind = find(((pit1VWCdateloc>734660)&(pit1VWCdateloc<734685))|...
    ((pit1VWCdateloc>734480)&(pit1VWCdateloc<734495))|...
    ((pit1VWCdateloc>734586)&(pit1VWCdateloc<734589))|...
    ((pit1VWCdateloc>734340)&(pit1VWCdateloc<734550))|...
    ((pit1VWCdateloc>734699)&(pit1VWCdateloc<734730)));
pit1VWC70R(ind) = NaN;

pit1VWC70R = pit1VWC70R.*0.92;

ind = find(((pit1VWCdateloc>733786)&(pit1VWCdateloc<733796))|...
    ((pit1VWCdateloc>734650)&(pit1VWCdateloc<734720))|...
    ((pit1VWCdateloc>733697)&(pit1VWCdateloc<733710))|...
    ((pit1VWCdateloc>734340)&(pit1VWCdateloc<734550))|...
    ((pit1VWCdateloc>734541)&(pit1VWCdateloc<734544)));
pit1VWC100R(ind) = NaN;

ix = smooth(pit1VWC100R,12);    ix(isnan(pit1VWC100R)) = NaN;	pit1VWC100R = ix;
ix = smooth(pit1VWC70R,12);     ix(isnan(pit1VWC70R)) = NaN;    pit1VWC70R = ix;

% .........................................................................
% Quick filling Left with Rigth
% .........................................................................
% % [~,~,pit1VWC05L] = AM_compute_regression(pit1VWC05R,pit1VWC05L,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC05R] = AM_compute_regression(pit1VWC05L,pit1VWC05R,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC10L] = AM_compute_regression(pit1VWC10R,pit1VWC10L,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC10R] = AM_compute_regression(pit1VWC10L,pit1VWC10R,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC20L] = AM_compute_regression(pit1VWC20R,pit1VWC20L,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC20R] = AM_compute_regression(pit1VWC20L,pit1VWC20R,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC40L] = AM_compute_regression(pit1VWC40R,pit1VWC40L,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC40R] = AM_compute_regression(pit1VWC40L,pit1VWC40R,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC70L] = AM_compute_regression(pit1VWC70R,pit1VWC70L,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC70R] = AM_compute_regression(pit1VWC70L,pit1VWC70R,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC100L] = AM_compute_regression(pit1VWC100R,pit1VWC100L,[0.4,0.4,0.4],0);
% % [~,~,pit1VWC100R] = AM_compute_regression(pit1VWC100L,pit1VWC100R,[0.4,0.4,0.4],0);

% Depth scaling of soil moisture content from surface to profile: multistation testing of observation operators
%  Gao et al. 2017
%%
pit1VWCshallowL = ((2.*pit1VWC05L.*5)+((pit1VWC05L+pit1VWC10L).*10)+((pit1VWC10L+pit1VWC20L).*20)+...
    ((pit1VWC20L+pit1VWC40L).*40))./(2.*(5+10+20+40));

pit1VWCshallowR = ((2.*pit1VWC05R.*5)+((pit1VWC05R+pit1VWC10R).*10)+((pit1VWC10R+pit1VWC20R).*20)+...
    ((pit1VWC20R+pit1VWC40R).*40))./(2.*(5+10+20+40));

pit1VWCdeepR = ((2.*pit1VWC70R.*30)+((pit1VWC70R+pit1VWC100R).*60))./(2.*(30+60));

pit1VWCdeepL = ((2.*pit1VWC70L.*20)+((pit1VWC70L+pit1VWC100L).*60))./(2.*(30+60));

pit1VWCL = ((2.*pit1VWC05L.*5)+((pit1VWC05L+pit1VWC10L).*10)+((pit1VWC10L+pit1VWC20L).*20)+...
    ((pit1VWC20L+pit1VWC40L).*40)+((pit1VWC40L+pit1VWC70L).*70)+((pit1VWC70L+pit1VWC100L).*100))./...
    (2.*(5+10+20+40+70+100));

pit1VWCR = ((2.*pit1VWC05R.*5)+((pit1VWC05R+pit1VWC10R).*10)+((pit1VWC10R+pit1VWC20R).*20)+...
    ((pit1VWC20R+pit1VWC40R).*40)+((pit1VWC40R+pit1VWC70R).*70)+((pit1VWC70R+pit1VWC100R).*100))./...
    (2.*(5+10+20+40+70+100));

% .........................................................................
color_mtx = jet(8);
figure('color','white');
subplot(3,1,1:3);       hold on;
hl1=plot(pit1VWCdateloc,pit1VWC05L,'color',color_mtx(1,:)); hl2=plot(pit1VWCdateloc,pit1VWC10L,'color',color_mtx(3,:));
hl3=plot(pit1VWCdateloc,pit1VWC20L,'color',color_mtx(4,:)); hl4=plot(pit1VWCdateloc,pit1VWC30L,'color',color_mtx(5,:));
hl5=plot(pit1VWCdateloc,pit1VWC40L,'color',color_mtx(6,:)); hl6=plot(pit1VWCdateloc,pit1VWC70L,'color',color_mtx(7,:));
hl7=plot(pit1VWCdateloc,pit1VWC100L,'color',color_mtx(8,:));
hl8=plot(pit1VWCdateloc,pit1VWCL,'color',[.5 .5 .5],'LineWidth',4);

plot(pit1VWCdateloc,pit1VWC05R,'-','color',color_mtx(1,:)); plot(pit1VWCdateloc,pit1VWC10R,'-','color',color_mtx(3,:));
plot(pit1VWCdateloc,pit1VWC20R,'-','color',color_mtx(4,:)); plot(pit1VWCdateloc,pit1VWC30R,'-','color',color_mtx(5,:));
plot(pit1VWCdateloc,pit1VWC40R,'-','color',color_mtx(6,:)); plot(pit1VWCdateloc,pit1VWC70R,'-','color',color_mtx(7,:));
plot(pit1VWCdateloc,pit1VWC100R,'-','color',color_mtx(8,:));
plot(pit1VWCdateloc,pit1VWCR,'color',[.5 .5 .5],'LineWidth',4);
ylabel('Volumetric water content soil (m^3 m^-^3)','FontSize',12);                box on; grid on;
datetick('x');          axis([pit1VWCdateloc(1) pit1VWCdateloc(end) 0.14 0.34])
set(gca,'FontSize',12);
legend([hl1,hl2,hl3,hl4,hl5,hl6,hl7,hl8],'0.5 m','1.0 m','2.0 m','3.0 m',...
    '4.0 m','7.0 m','10.0 m','integrated','BagroundColor','white','FontSize',10);
legend box off;     set(legend,'color','white');
%%
pit1VWC = (pit1VWCR+pit1VWCL)./2;
pit1VWC(pit1VWC>0.306) = NaN;

[~,~,pit1VWC] = AM_compute_regression(pit1VWCR,pit1VWC,[0.4,0.4,0.4],0);
[~,~,pit1VWC] = AM_compute_regression(pit1VWCL,pit1VWC,[0.4,0.4,0.4],0);
plot(pit1VWCdateloc,pit1VWC,'color','k','LineWidth',4);

pit1VWCshallow = (pit1VWCshallowR+pit1VWCshallowL)./2;
[~,~,pit1VWCshallow] = AM_compute_regression(pit1VWCshallowR,pit1VWCshallow,[0.4,0.4,0.4],0);
[~,~,pit1VWCshallow] = AM_compute_regression(pit1VWCshallowL,pit1VWCshallow,[0.4,0.4,0.4],0);

pit1VWCdeep = (pit1VWCdeepR+pit1VWCdeepL)./2;
[~,~,pit1VWCdeep] = AM_compute_regression(pit1VWCdeepR,pit1VWCdeep,[0.4,0.4,0.4],0);
[~,~,pit1VWCdeep] = AM_compute_regression(pit1VWCdeepL,pit1VWCdeep,[0.4,0.4,0.4],0);

[pit1VWCweekALL,~,~,~,~,~,pit1VWCweek] = AM_week2day_rs(pit1VWCR,pit1VWCdateloc);
[pit1VWCweekDEEP] = AM_week2day_rs(pit1VWCdeep,pit1VWCdateloc);
[pit1VWCweekSHALLOW] = AM_week2day_rs(pit1VWCshallow,pit1VWCdateloc);

K67.pit1VWCweekALL = interp1(pit1VWCweek,pit1VWCweekALL,K67.week);
K67.pit1VWCshallow = interp1(pit1VWCdateloc,pit1VWCshallow,K67.dateloc);
K67.pit1VWCdeep    = interp1(pit1VWCdateloc,pit1VWCdeep,K67.dateloc);

% .........................................................................
%% Compare to CWD
% .........................................................................
% close all;
figure ('color','white');
subplot(2,3,1:3);       hold on;
hl1 = plot(pit1VWCdateloc,(nanmax(pit1VWC)-pit1VWC),'color','k','LineWidth',2);
ylabel({'\theta_{max}-\theta'},'FontSize',12,'Color','k');
set(gca,'FontSize',12);
yyaxis right
hl2 = plot(TRMM1h.dateloc,TRMM1h.CWD,'color',[.9 .3 .2]);  datetick('x');
ylabel('CWD (mm)','color',[.9 .3 .2],'FontSize',12);                        box on;
legend([hl1,hl2],'Volumetric water content soil (\theta_{max}-\theta)','CWD');     legend box off
set(gca,'FontSize',12);         box on;                 grid on;                xtickangle(0);

subplot(2,3,4);     hold on;
iq = interp1(TRMM1h.dateloc,TRMM1h.CWD,pit1VWCdateloc);
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,pit1VWC,[0.4,0.4,0.4],1);
RMSE = sqrt(nanmean((pit1VWC - yhat).^2));
[bin_x,bin_y] = AM_bin_noplot(iq,pit1VWC,20,10);
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(log(bin_x),bin_y,[0.3,0.3,0.7],0);
y_model = polyval(p,log(sort(iq)));
K67.SWintegrated = polyval(p,log(K67.TRMM1hCWD));
plot(sort(iq),y_model,'-k');
text(50,0.24,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
%     text(70,0.23,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
text(50,0.23,['\theta=' num2str(p(1),2.1) '+' num2str(p(2),2.1) 'ln(CWD)'],'color',[0 0 0]);
xlabel('CWD (mm)','FontSize',12);
ylabel({'\theta_{integrated}(m^3 m^-^3)'},'FontSize',12);               box on; grid on;
set(gca,'XTick',(0:200:600),'XTickLabel',(0:200:600),'YTick',(0.22:0.03:0.31),'YTickLabel',(0.22:0.03:0.31),...
    'Ylim',[0.22 0.31],'Xlim',[0 600],'FontSize',12);       xtickangle(0);

subplot(2,3,5);                 hold on;
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,pit1VWCshallow,[.9 .3 .2],1);
RMSE = sqrt(nanmean((pit1VWCshallow - yhat).^2));
[bin_x,bin_y] = AM_bin_noplot(iq,pit1VWCshallow,20,10);
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(log(bin_x),bin_y,[0.3,0.3,0.7],0);
y_model = polyval(p,log(sort(iq)));
plot(sort(iq),y_model,'-k');
K67.SWshallow = polyval(p,log(K67.TRMM1hCWD));
K67.SWshallow(isinf(K67.SWshallow)) = NaN;
text(50,0.23,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
%     text(70,0.22,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
text(50,0.22,['\theta=' num2str(p(1),2.1) '+' num2str(p(2),2.1) 'ln(CWD)'],'color',[0 0 0]);
xlabel('CWD (mm)','FontSize',12);
ylabel({'\theta_{shallow}(m^3 m^-^3)'},'FontSize',12);     box on;     grid on;
set(gca,'XTick',(0:200:600),'XTickLabel',(0:200:600),'YTick',(0.16:0.04:0.28),'YTickLabel',(0.16:0.04:0.28),...
    'Ylim',[0.16 0.28],'Xlim',[0 600],'FontSize',12);   xtickangle(0);

subplot(2,3,6);                 hold on;
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,pit1VWCdeep,[0.3,0.3,0.7],1);
RMSE = sqrt(nanmean((pit1VWCdeep - yhat).^2));
[bin_x,bin_y] = AM_bin_noplot(iq,pit1VWCdeep,20,10);
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(log(bin_x),bin_y,[0.3,0.3,0.7],0);
y_model = polyval(p,log(sort(iq)));
K67.SWdeep = polyval(p,log(K67.TRMM1hCWD));
plot(sort(iq),y_model,'-k');
text(50,0.27,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
%     text(70,0.26,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
text(50,0.26,['\theta=' num2str(p(1),2.1) '+' num2str(p(2),2.1) 'ln(CWD)'],'color',[0 0 0]);
xlabel('CWD (mm)','FontSize',12);
set(gca,'XTick',(0:200:600),'XTickLabel',(0:200:600),'YTick',(0.25:0.02:0.31),'YTickLabel',(0.25:0.02:0.31),...
    'Ylim',[0.25 0.31],'Xlim',[0 600],'FontSize',12);
ylabel({'\theta_{deep}(m^3 m^-^3)'},'FontSize',12);               box on; grid on;xtickangle(0);

%% ........................................................................
ix = (pit1VWCweekALL-nanmin(pit1VWCweekALL))./(nanmax(pit1VWCweekALL)-nanmin(pit1VWCweekALL));
iy = (TRMM1h.weekCWD-nanmin(TRMM1h.weekCWD))./(nanmax(TRMM1h.weekCWD)-nanmin(TRMM1h.weekCWD));
% ix = 1-ix;

iw = (pit1VWCweekSHALLOW-nanmin(pit1VWCweekSHALLOW))./(nanmax(pit1VWCweekSHALLOW)-nanmin(pit1VWCweekSHALLOW));
% iw = 1-iw;
iz = (pit1VWCweekDEEP-nanmin(pit1VWCweekDEEP))./(nanmax(pit1VWCweekDEEP)-nanmin(pit1VWCweekDEEP));
% iz = 1-iz;

figure ('color','white');
subplot(2,3,1:3);  hold on;
hl1 = plot(pit1VWCweek,pit1VWCweekALL,'color','k','LineWidth',2);
set(gca,'FontSize',12);         box on;                 grid on;
ylabel({'\theta_{normalized}'},'FontSize',12,'Color','k');
yyaxis right

iq = interp1(TRMM1h.week,TRMM1h.weekCWD,pit1VWCweek);
hl2 = plot(pit1VWCweek,iq,'color',[.7 .2 .2]);
datetick('x');       box on;                grid on;
legend([hl1,hl2],'Volumetric water content soil (\theta)','CWD','possition',...
    'southoutside','orientation','horizontal'); legend box off
ylabel('CWD (mm)','FontSize',12);                        box on;

iq = interp1(TRMM1h.week,iy,pit1VWCweek);
subplot(2,3,4);
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,ix,[0.4,0.4,0.4],1);
RMSE = sqrt(nanmean((ix - yhat).^2));
text(0.1,0.9,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
text(0.1,0.8,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
set(gca,'XTick',(0:.2:1),'XTickLabel',(0:.2:1),'YTick',(0.:0.2:1),'YTickLabel',(0.:0.2:1),...
    'Ylim',[0. 1.],'Xlim',[0 1],'FontSize',12);
xlabel('CWD (mm month^-^1)','FontSize',12);
ylabel({'\theta_{integrated normalized}'},'FontSize',12);               box on; grid on;
subplot(2,3,5);                 hold on;
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,iw,[0.7,0.3,0.3],1);
RMSE = sqrt(nanmean((iw - yhat).^2));
axis([0 1 0 1 ]);        grid on;
text(0.1,0.9,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
text(0.1,0.8,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
set(gca,'XTick',(0:.2:1),'XTickLabel',(0:.2:1),'YTick',(0.:0.2:1),'YTickLabel',(0.:0.2:1),...
    'Ylim',[0. 1.],'Xlim',[0 1],'FontSize',12);
xlabel('CWD (mm month^-^1)','FontSize',12);
ylabel({'\theta_{shallow normalized}'},'FontSize',12);               box on; grid on;
subplot(2,3,6);                 hold on;
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,iz,[0.3,0.3,0.7],1);
RMSE = sqrt(nanmean((iz - yhat).^2));
axis([0 1 0 1]);                    grid on;
text(0.1,0.9,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
text(0.1,0.8,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
set(gca,'XTick',(0:.2:1),'XTickLabel',(0:.2:1),'YTick',(0.:0.2:1),'YTickLabel',(0.:0.2:1),...
    'Ylim',[0. 1.],'Xlim',[0 1],'FontSize',12);
xlabel('CWD (mm month^-^1)','FontSize',12);
ylabel({'\theta_{deep normalized}'},'FontSize',12);               box on; grid on;

%% ........................................................................
iq = interp1(TRMM1h.week,TRMM1h.weekCWD,pit1VWCweek);

figure ('color','white');
subplot(2,3,1:3);  hold on;
hl1 = plot(pit1VWCweek,pit1VWCweekALL,'color','k','LineWidth',2);
set(gca,'Ylim',[0.22 0.32],'Xlim',[datenum(2001,12,31) datenum(2022,12,31)],'FontSize',12); 
box on;                 grid on;
ylabel({'\theta (m^3 m^-^3)'},'FontSize',12,'Color','k');
yyaxis right
hl2 = plot(pit1VWCweek,iq,'color',[.9 .3 .2]);
datetick('x');       box on;                grid on;        axis tight;
% set(gca,'Ylim',[0 1200],'Xlim',[datenum(2001,12,31) datenum(2020,12,31)],...
%    'FontSize',12);         box on;                 grid on;
legend([hl1,hl2],'Volumetric water content soil (\theta)','CWD'); legend box off
ylabel('CWD (mm)','FontSize',12);                        box on;     xtickangle(0);

pit1VWCweekALL(round(pit1VWCweekALL,4)==0.2643) = NaN;
subplot(2,3,4);
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,pit1VWCweekALL,[0.4,0.4,0.4],1);
RMSE = sqrt(nanmean((pit1VWCweekALL - yhat).^2));
text(10,0.295,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
text(10,0.29,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
set(gca,'XTick',(0:200:800),'XTickLabel',(0:200:800),'YTick',(0.22:0.02:.30),'YTickLabel',(0.22:0.02:.3),...
    'Ylim',[0.22 .30],'Xlim',[0 800],'FontSize',12);
xlabel('CWD (mm)','FontSize',12);
ylabel({'\theta_{integrated}(m^3 m^-^3)'},'FontSize',12);   box on; grid on; xtickangle(0);
subplot(2,3,5);                 hold on;
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,pit1VWCweekSHALLOW,[0.9,0.4,0.4],1);
RMSE = sqrt(nanmean((pit1VWCweekSHALLOW - yhat).^2));
axis([0 1 0 1 ]);        grid on;            xtickangle(0);
text(10,0.235,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
text(10,0.23,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
set(gca,'XTick',(0:200:800),'XTickLabel',(0:200:800),'YTick',(0.16:0.02:.24),'YTickLabel',(0.16:0.02:.24),...
    'Ylim',[0.16 .26],'Xlim',[0 800],'FontSize',12);
xlabel('CWD (mm)','FontSize',12);
ylabel({'\theta_{shallow}(m^3 m^-^3)'},'FontSize',12);  box on; grid on; xtickangle(0);
subplot(2,3,6);                 hold on;
[p,r2,~,p_value,~,~,yhat] = AM_compute_regression(iq,pit1VWCweekDEEP,[0.4,0.4,0.9],1);
RMSE = sqrt(nanmean((pit1VWCweekDEEP - yhat).^2));
axis([0 1 0 1]);        grid on;             xtickangle(0);
text(10,0.295,['R^2=',num2str(r2,2.1)],'color',[0 0 0]);
text(10,0.29,['RMSE=',num2str(RMSE,2.1)],'color',[0 0 0]);
set(gca,'XTick',(0:200:800),'XTickLabel',(0:200:800),'YTick',(0.25:0.01:.30),'YTickLabel',(0.25:0.01:.30),...
    'Ylim',[0.25 .30],'Xlim',[0 800],'FontSize',12);
xlabel('CWD (mm)','FontSize',12);
ylabel({'\theta_{deep}(m^3 m^-^3)'},'FontSize',12); box on; grid on; xtickangle(0);


K67.SWshallow((isnan(K67.SWshallow))|(isinf(K67.SWshallow)))          = NaN;
K67.SWdeep((isnan(K67.SWdeep))|(isinf(K67.SWdeep)))                   = NaN;
K67.SWintegrated((isnan(K67.SWintegrated))|(isinf(K67.SWintegrated))) = NaN;

K67.weekSWshallow = AM_week2day_rs(K67.SWshallow,K67.dateloc);
K67.weekSWdeep = AM_week2day_rs(K67.SWdeep,K67.dateloc);
K67.weekSWintegrated = AM_week2day_rs(K67.SWintegrated,K67.dateloc);
