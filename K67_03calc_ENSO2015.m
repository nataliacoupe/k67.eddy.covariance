% Calculate seasonal values for ENSO event or other analysis 
% Natalia Restrepo-Coupe
% Noisy PAR
% h2o mixing ratio improperly calculated
% u* improperly calculated
%..........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;

%..........................................................................
flag_timeseries   = 1;      flag_Rick = 1;
month_start_hydro = 7;
ix = (1:16:365);            ind = find(ix>=datenum(0,7,1));     week_start_hydro = ind(1);
tresholdSo    = 10;         tresholdCycleDay = 6;               tresholdCycle = 22;

eight_avg    = (1:8:365)';  week_avg = (1:16:365)';
[~,~,K67.eightTRMM] = AM_eightday_rs(K67.TRMM1h,K67.dateloc);                         
K67.eightTRMM = AM_fill(K67.eightTRMM);

K67.weekTRMM = K67.weekTRMM1h;

[K67.cycle_month_TRMM,K67.cycle_month_countTRMM,~,K67.cycle_month_stdTRMM] = ...
    AM_day_cycle_year(K67.TRMM1h,K67.dateloc);
[K67.week_avgTRMM,~,K67.week_avgstdTRMM] = AM_week2_avg(K67.weekTRMM,K67.week);

% .........................................................................
%% Missing precipitation
% ........................................................................
figure;     hold on;	plot(K67.dateloc,K67.TRMM1hPrecip);    plot(K67.dateloc,K67.rain);
legend('TRMM','K67');   grid on;    datetick('x');  ylabel('Precip (mm)');  box on;
% figure;     hold on;	plot(K67.rain,K67.TRMM1hPrecip);
[K67.CI,K67.PARtheory] = AM_CI(K67.par,K67.dateloc,lat,lon,utm);
K67.diffTair = [NaN;diff(K67.tasp_fill)];
%..........................................................................
%   Linear regression models
%..........................................................................
% % % CI and TRMM for daytime  R2=0.2626
% % ind = find ((K67.Y<2018)&(K67.TRMM1hPrecip>0)&(K67.rain>0)&(K67.So>500));
% % tbl = table(K67.CI(ind),...
% %     K67.TRMM1hPrecip(ind),K67.rain(ind),'VariableNames',{'CI','TRMM','rain'});
% % lm  = fitlm(tbl,'rain~CI+TRMM')
% % AICvalue1 = lm.ModelCriterion.AIC
% %
% % rain_model = lm.Coefficients.Estimate(1)+...
% %     (lm.Coefficients.Estimate(2).*K67.CI)+...
% %     (lm.Coefficients.Estimate(3).*K67.TRMM1hPrecip);
% % rain_model((rain_model<0)|(K67.TRMM1hPrecip==0)) = 0;
% %
% % ind = find((~isnan(K67.rain))&(~isnan(rain_model)));
% % corrcoef(rain_model(ind),K67.rain(ind))
% % nanmean(rain_model(ind)-K67.rain(ind))
% % plot(K67.dateloc,rain_model);
% %
% % ind = find((isnan(K67.rain))&(K67.So>500));
% % K67.rain_fill = K67.rain;      K67.rain_fill(ind) = rain_model(ind);

% tair and TRMM for nighttime
% ind = find ((K67.Y<2018)&(K67.TRMM1hPrecip>0)&(K67.rain>0)&(K67.So<=500));
ind = find ((K67.Y<2020)&(K67.TRMM1hPrecip>0)&(K67.rain>0)&(~isnan(K67.diffTair)));
tbl = table(K67.diffTair(ind),...
    K67.TRMM1hPrecip(ind),K67.rain(ind),'VariableNames',{'diffTair','TRMM','rain'});
lm  = fitlm(tbl,'rain~diffTair+TRMM')
AICvalue1 = lm.ModelCriterion.AIC

rain_model = lm.Coefficients.Estimate(1)+...
    (lm.Coefficients.Estimate(2).*K67.diffTair)+...
    (lm.Coefficients.Estimate(3).*K67.TRMM1hPrecip);
rain_model((rain_model<0)|(K67.TRMM1hPrecip==0)) = 0;

ind = find((~isnan(K67.rain))&(~isnan(rain_model)));
corrcoef(rain_model(ind),K67.rain(ind))
nanmean(rain_model(ind)-K67.rain(ind))
plot(K67.dateloc,rain_model);

K67.rain_fill = K67.rain;
ind = find(isnan(K67.rain_fill));
K67.rain_fill(ind) = rain_model(ind);

%% ........................................................................
figure('color','white');
subplot(2,1,1);     plot(K67.dateloc,K67.storagewt); hold on;  plot(K67am.dateloc,K67am.storage);
plot(K67.dateloc,K67.storagewt,'k');   	datetick('x');      ylabel('Sco_2');
legend('k67','Ameriflux');
subplot(2,1,2); plot(K67.dateloc,K67.fco2); hold on;    plot(K67am.dateloc,K67am.fco2);
datetick('x');      ylabel('Fc');           legend('k67','Ameriflux');

%..........................................................................
%% Fill meteorology
%..........................................................................
K67.ws_fill     = AM_rm_outlier(K67.ws,4);
K67.ws_fill     = AM_fill_month_cycle(K67.ws_fill,K67.dateloc,6,18);
K67.tc_fill     = AM_fill_day_cycle(K67.tasp,K67.dateloc,6,18);

%..........................................................................
%% Saturated and air water pressure
%..........................................................................
T_0C_K =  273.15;             %'Temperature in K when it is zero in C
R      =  8.3143e-3;          %'Universal gas constant [kPa m^3/(K mol)]
Rd     =  R/28.97;            %'Gas constant for dry air [kPa m^3/(K g)]
Cp     =  1.013*10^-3;        % Specific heat at constant pressure, 1.013 10-3 [MJ kg-1 �C-1];
k      =  0.41;               %'von Karman constant (Dyer & hicker 1970, Webb 1970)
A0 = 6.107799961;       A1 = 4.436518521e-1;        A2 = 1.428945805e-2;
A3 = 2.650648471e-4;    A4 = 3.031240396e-6;        A5 = 2.034080948e-8;
A6 = 6.136820929e-11;   RV = 0.0004615;             %Gas constant for water vapor [J/(mg K)]
K67.e_sat   = 0.1.*(A0+K67.tc_fill.*(A1+K67.tc_fill.*(A2+K67.tc_fill.*(A3+K67.tc_fill.*(A4+K67.tc_fill.*(A5+K67.tc_fill*A6))))));

K67.h2o_gm3 =  (0.01802.*K67.pambKpa)./((R.*K67.tc_fill+T_0C_K).*((1./K67.h2ommolm)+(0.01802/18)));

%'Calculate water vapor pressure using IRGASON measurements
K67.e     =  (K67.h2omix.*(K67.pambKpa./100)./(K67.h2omix+622)).*10^2; % vapor pressure (mb) then x10^6
K67.e     =  AM_rm_outlier(K67.e,3);
% K67.RH    =  K67.e./K67.e_sat;
ix        =  K67.e_sat.*(1-K67.RH./100);                                    % vapour pressure deficit;
K67.VPD   =  K67.e_sat.*(1-K67.RHec./100);
K67.VPD((K67.VPD<0)|(isnan(K67.VPD))) = ix((K67.VPD<0)|(isnan(K67.VPD)));
K67.VPD(K67.VPD<0) = 0;
K67.VPD_fill = AM_fill_month_cycle(K67.VPD,K67.dateloc,6,18);
K67.VPD_fill((K67.VPD_fill<0)|(isinf(K67.VPD_fill))) = 0;

K67.rho_d =  (K67.pambKpa-K67.e)./((K67.tc_fill+T_0C_K).*Rd);          % 'g/m^3
%K67.rho_a =  (K67.rho_d + K67.h2o_gm3)./1000;       %Natalia          % 'kg/m^3

% % %..........................................................................
K67.e_sat_fill = AM_fill_day_cycle(K67.e_sat,K67.dateloc,6,18);
K67.rho_a_fill = AM_fill_day_cycle(K67.rho_a,K67.dateloc,6,18);
K67.e_fill     = AM_fill_day_cycle(K67.e,K67.dateloc,6,18);

% .........................................................................
% lambda       =  2.45;                                                    %[MJ kg-1] FAO latent heat of vaporization of water
% K67.delta    =  (4098.*K67.e_sat_fill)./(power((K67.tc_fill+237.3),2));  %[kPa C-1]	 % Slope of the saturation vapour pressure temperature relationship
% http://agsys.cra-cin.it/tools/evapotranspiration/help/Slope_of_saturation_vapour_pressure_curve.html
% K67.delta   =  0.04145.*exp(0.06088.*K67.tc_fill);                       % 0.3 at 30degC [kPa C-1]	 % Slope of the saturation vapour pressure temperature relationship
K67.delta   = 4098*K67.e_sat_fill./((237.3+K67.tc_fill).^2);                % slope of vapor pressure curve, kPa/C

K67.sigma   = ((Cp.*K67.pambKpa)./(0.622.*K67.lambda));                     % sigma psychrometric constant 0.054 [kPa C-1]
%K67.sigma   = smooth(K67.sigma);
%where epsilon ratio molecular weight of water vapour/dry air  =  0.622
K67.epsilon =  K67.delta./K67.sigma;             %Unitless ~0.6219
% Rho_air: mean air density at constant pressure (Bruutsaet,p38)(kPa kg J-1)
% % if ((nanmax(K67.rho_air)>2)|(nanmin(K67.rho_air)<0.5));   '\rho air inacurate calculation please check';    end;
% % c = 1.01;                                                      % virtual absolute air temperature Tv = (c.*(K67.tc_fill+273.15)
% % K67.rho_air = K67.pambpa./((c.*(K67.tc_fill+273.15)).*R);      % multiply by 10 so we can get the 1.2value??? http://www.fao.org/docrep/x0490e/x0490e06.htm [kg m-3]
% % Qa  =  0.62197.*K67.e./((K67.pambpa./100)-(0.378030.*K67.e));  % specific humidity (kg/kg)
% %..........................................................................
% % 2002-2005
% %..........................................................................
% [K67.h_So]     = AM_solar_rad(K67.h_JD,K67.h_HR+1,lat,lon);
% K67.h_delta    =  0.04145.*exp(0.06088.*K67.h_tc_fill);                       %0.3 at 30degC [kPa C-1]	 % Slope of the saturation vapour pressure temperature relationship
% K67.h_lambda   =  (2.50*10^6)-((2.386*10^3).*(K67.h_tc_fill));                %2.45 [MJ kg-1]     % Latent heat of vaporization, 2.45 [MJ kg-1],
% K67.h_lambda   =  K67.h_lambda./10^6;
% K67.h_lambda(isnan(K67.h_lambda)) = 2.45;
% K67.h_sigma    =  ((Cp.*K67.h_pambKpa)./(0.622.*K67.h_lambda));                    %[kPa C-1]
% K67.h_sigma    =  smooth(K67.h_sigma);

%..........................................................................
%% ET ref
%..........................................................................
K67.LEref = AM_fill((K67.LE+K67.H).*(K67.delta./(K67.delta+(K67.sigma))));
% K67.LEref = AM_fill((K67.netrad).*(K67.delta./(K67.delta+(K67.sigma))));
K67.LEref((K67.So < 500)|(K67.rain_fill>0)) = NaN;
K67.ETref = AM_fill_plus(K67.LEref./K67.lambda);
K67.ETref = K67.ETref.*60.*60.*24/1000000;      %from s-1 to d-1 and from kg-1 to g-1 (mm)

%..........................................................................
K67.eightETref = AM_eight2cycle_rs(K67.ETref,K67.dateloc,tresholdCycleDay,4,2);    K67.eightETref = AM_fill(K67.eightETref);
[K67.weekETref,~,~,~,~,K67.week_stdETref]  = AM_week2cycle_rs(K67.ETref,K67.dateloc,tresholdCycleDay,4,2);     K67.weekETref  = AM_fill(K67.weekETref);
K67.monthETref = AM_month2cycle(K67.ETref,K67.dateloc,tresholdCycleDay,4,2);       K67.monthETref = AM_fill(K67.monthETref);
K67.monthETref = AM_month2cycle(K67.ETref,K67.dateloc,tresholdCycleDay,4,2);                      K67.monthETref  = AM_fill(K67.monthETref);
K67.monthETref = K67.monthETref.*eomday(K67.monthY,K67.monthM);

K67.ETrefday = K67.ETref;      K67.ETrefday(K67.So<=tresholdSo) = NaN;
[K67.weekETrefday,~,~,~,~,K67.week_stdETrefday]  = AM_week2cycle_rs(K67.ETrefday,K67.dateloc,tresholdCycleDay,4,2);     K67.weekETrefday  = AM_fill(K67.weekETrefday);
[K67.monthETrefday,~,~,~,~,K67.month_stdETrefday]= AM_month2cycle(K67.ETrefday,K67.dateloc,tresholdCycleDay,4,2);     K67.monthETrefday  = AM_fill(K67.monthETrefday);

K67.LErefday = K67.LEref;      K67.LErefday(K67.So<=tresholdSo) = NaN;
[K67.weekLErefday,~,~,~,~,K67.week_stdLErefday]  = AM_week2cycle_rs(K67.LErefday,K67.dateloc,tresholdCycleDay,4,2);     K67.weekLErefday  = AM_fill(K67.weekLErefday);

%..........................................................................
%% Rain
%..........................................................................
[~,~,K67.eightRain] = AM_eightday_rs(K67.rain,K67.dateloc);                         K67.eightRain = AM_fill(K67.eightRain);
[~,~,K67.weekRain,~,~,K67.week_stdRain]  = AM_week2day_rs(K67.rain,K67.dateloc);    K67.weekRain  = AM_fill(K67.weekRain);
[~,K67.monthRain] = AM_month(K67.rain,K67.dateloc);                                 K67.monthRain = AM_fill(K67.monthRain);
[K67.cycle_month_Rain,K67.cycle_month_countRain,~,K67.cycle_month_stdRain] = ...
    AM_day_cycle_year(K67.rain,K67.dateloc);
[K67.week_avgRain,~,K67.week_avgstdRain] = AM_week2_avg(K67.weekRain,K67.week);

% % [K67.eightCWD] = AM_eightday_rs(K67.dayTRMM1hCWD,K67.day);                            K67.eightCWD = AM_fill(K67.eightCWD);
% % [K67.weekCWD,~,~,~,~,K67.week_stdCWD]  = AM_week2day_rs(K67.dayTRMM1hCWD,K67.day);    K67.weekCWD  = AM_fill(K67.weekCWD);
% % [K67.week_avgCWD,~,K67.week_avgstdCWD] = AM_week2_avg(K67.weekCWD,K67.week);

[K67.eightCWD] = AM_eightday_rs(K67.dayTRMM1hCWD,K67.day);                            K67.eightCWD = AM_fill(K67.eightCWD);
[K67.weekCWD,~,~,~,~,K67.week_stdCWD]  = AM_week2day_rs(K67.dayTRMM1hCWD,K67.day);    K67.weekCWD  = AM_fill(K67.weekCWD);
[K67.week_avgCWD,~,K67.week_avgstdCWD] = AM_week2_avg(K67.weekCWD,K67.week);

%% ........................................................................
if flag_Rick == 1
    ind = find((K67.dateloc>datenum(2004,12,10))&(K67.dateloc<datenum(2004,12,20)));
    ind = find((K67.dateloc>datenum(2004,8,24))&(K67.dateloc<datenum(2004,9,7)));
    figure('color','white');
    y1 = -.2;     y2 = 1;       y3 = y1:((y2-y1)/4):y2;
    x1 = 0;       x2 = 2;       x3 = x1:((x2-x1)/4):x2;
    subplot(2,2,1:2);           hold on;
    hl1 = bar(K67.dateloc(ind),K67.rain(ind)./10,1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.8 0.8 0.8]);
    hl3 = plot(K67.dateloc(ind),K67.ET(ind)./24,'Color',[0. 0. .9],'LineWidth',1);
    hl4 = plot(K67.dateloc(ind),K67.ETref(ind)./24,'Color',[0. 0.9 .5],'LineWidth',1);
    legend ([hl3,hl4],'ET','ET_{ref}','EdgeColor','white')%,'position','southoutside');
    set(gca,'XTick',K67.day(1:1:end),'XTickLabel',[ '' '' ''],'YTick',y3,'YColor',[0. 0. .9]); hold on;
    xlim([K67.dateloc(ind(1)) K67.dateloc(ind(end))]);            ylim([y1 y2]);
    ylabel ({'ET (mm)'},'FontSize',12,'Color',[0. 0. .9]);
    ax1 = gca;          ax1_pos = ax1.Position;
    ax2 = axes('Position',ax1_pos,'YAxisLocation','right','Color','none'); hold on;
    %     hl7 = plot(K67.dateloc(ind),K67.tasp(ind),'Parent',ax2,'Color',[0.4 0.4 0.4],'LineWidth',1);    hold on;
    hl7 = plot(K67.dateloc(ind),K67.ustar(ind),'Parent',ax2,'Color',[0.4 0.4 0.4],'LineWidth',1);    hold on;
    set(ax2,'XColor','k','YColor',[ 0.0 0.0 0.0]);
    %     ylabel ({'NetRad (W m^-^2)'},'FontSize',12,'Color',[0.2 0.2 .2]);
    %     ylabel ({'T_{air} (\circC)'},'FontSize',12,'Color',[0.4 0.4 .4]);
    ylabel ({'u_* (m s^-^1)'},'FontSize',12,'Color',[0.4 0.4 .4]);
    set(ax2,'XTick',K67.dateloc(ind(1):48:ind(end)),'XTickLabel',datestr(K67.dateloc(ind(1):48:ind(end)),'yyyy/mm/dd'),'YTick',x3,'YColor',[0.4 0.4 .4]);
    xlim([K67.dateloc(ind(1)) K67.dateloc(ind(end))]);            ylim([x1 x2]);
    box on; grid on;	text(K67.dateloc(ind(20)),x3(5),site_name,'FontSize',8)
    
    ind = find(K67.So>100);
    subplot(2,2,3);     hold on;
    AM_bin(K67.netrad(ind),K67.LE(ind),15,1,1,[0 0 0.8]);
    xlabel({'Rn (W m^-^2)'},'FontSize',10);ylabel ({'LE (W m^-^2)'},'FontSize',10);
    [p,r2,~,~,sm,sb] = AM_compute_regression(K67.netrad(ind),K67.LE(ind),[0.4,0.4,0.4],0);
    y_model=polyval(p,K67.netrad(ind));   hold on;
    plot(K67.netrad(ind),y_model,'-','color',[0.4,0.4,0.4]);
    xlim([-200 1200]);      ylim([-100 600]);
    text (-50,350,site_name{:},'FontSize',10);
    text (-50,300,{['LE=' num2str(p(1),'%6.4f') 'Rn+' num2str(p(2),'%6.4f')];...
        ['std slope:' num2str(sm(1),'%3.4f')];['std int.:' num2str(sb(1),'%3.4f')];...
        ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10);         box on;
end

%..........................................................................
%% Canopy conductance
%..........................................................................
z   = height;
d   = 2.*canopy./3;
zo  = 0.123*canopy;
zoh = 0.0123.*canopy;

K67.A = (K67.LE+K67.H);
% weekA = AM_week2cycle_rs(K67.A,K67.dateloc,20);  weekH = AM_week2cycle_rs(K67.H,K67.dateloc,20);  weekLE = AM_week2cycle_rs(K67.LE,K67.dateloc,3);
%%
%         ACHTUNG!!!!!!!!!!!!!
% ik Yebra it apears as K67.ga but looks like K67.ra!!!!!!!!!!!
% K67.ra =
% (log((height-d)./zo)).*(log((height-d)./zoh))./(k^2.*K67.ws_fill);
% Natalia's original I use both zo and zoh
Ra        = (log((z-d)./zo)).*(log((z-d)./zoh))./(k.^2.*K67.ws_fill);
Ra(K67.So<tresholdSo) = NaN;        Ra = AM_rm_outlier(Ra,5);
K67.Ga        = 1./Ra;              K67.Ga = AM_rm_outlier(K67.Ga,5);
Ra_ustar  = K67.ws_fill./(K67.ustar.^2);
Ra_ustar(K67.So<tresholdSo) = NaN;	Ra_ustar = AM_rm_outlier(Ra_ustar,5);
Ga_ustar  = 1./Ra_ustar;            Ga_ustar = AM_rm_outlier(Ga_ustar,5);

%..........................................................................
%% PET from PM equation
%..........................................................................
% rho_a 1.2 kg m-3  Cp = 1.013 10-3 [MJ kg-1 C-1] > [J kg-1 k-1]   >>> delta [kPa m-2]
% lambda 2.45 [MJ kg-1]
% K67.LErefPM = ((AM_fill((K67.LE+K67.H).*K67.delta)+(K67.rho_a_fill.*Cp.*10^6.*K67.VPD_fill.*K67.Ga)))./(K67.delta+K67.sigma);
K67.LErefPM = ((AM_fill((K67.netrad).*K67.delta)+(K67.rho_a_fill.*Cp.*10^6.*K67.VPD_fill.*K67.Ga)))./(K67.delta+K67.sigma);
K67.ETrefPM = K67.LErefPM./K67.lambda;
K67.ETrefPM = K67.ETrefPM.*60.*60.*24/1000000;      %from s-1 to d-1 and from kg-1 to g-1 (mm)
[K67.weekETrefPM,~,~,~,~,K67.week_stdETrefPM]  = AM_week2cycle_rs(K67.ETrefPM,K67.dateloc,tresholdCycleDay,4,2); 
K67.weekETrefPM  = AM_fill(K67.weekETrefPM);
[K67.week_avgETrefPM,~,K67.week_avgstdETrefPM] = AM_week2_avg(K67.weekETrefPM,K67.week);

K67.ETrefdayPM = K67.ETrefPM;      K67.ETrefdayPM(K67.So<=tresholdSo) = NaN;
[K67.weekETrefdayPM,~,~,~,~,K67.week_stdETrefdayPM]  = AM_week2cycle_rs(K67.ETrefdayPM,K67.dateloc,tresholdCycleDay,4,2);
K67.weekETrefdayPM  = AM_fill(K67.weekETrefdayPM);
[K67.week_avgETrefdayPM,~,K67.week_avgstdETrefdayPM] = AM_week2_avg(K67.weekETrefdayPM,K67.week);

%..........................................................................
%% ET dry conditions
%..........................................................................
ik = find(K67.rain_fill>0);
ij = zeros(length(K67.rain_fill),1);
ij(ik) = 1;                         indDry24 = ij;           im = 1;
while im < length(ij)-25
    if ij(im) == 1
        for in = 1:24
            indDry24(im+in) = 1;
        end
    end
    im = im + 1;
end
K67.LEdry24 = K67.LE;             K67.ETdry24 = K67.ET;
K67.LEdry24(indDry24==1) = NaN;     K67.ETdry24(indDry24==1) = NaN;
K67.ETdryday24 = K67.ETdry24;     K67.ETdryday24(K67.So<tresholdSo) = NaN;

K67.ETrefdry24 = K67.ETref;           K67.ETrefdry24(indDry24==1) = NaN;
K67.ETrefdryday24 = K67.ETrefdry24;   K67.ETrefdryday24(K67.So<tresholdSo) = NaN;

K67.LEwet24 = K67.LE;             K67.ETwet24 = K67.ET;
K67.LEwet24(indDry24==0) = NaN;     K67.ETwet24(indDry24==0) = NaN;
K67.ETwetday24 = K67.ETwet24;     K67.ETwetday24(K67.So<tresholdSo) = NaN;

K67.ETrefwet24 = K67.ETref;           K67.ETrefwet24(indDry24==0) = NaN;
K67.ETrefwetday24 = K67.ETrefwet24;   K67.ETrefwetday24(K67.So<tresholdSo) = NaN;

%% ........................................................................
ik = find(K67.rain_fill>0);
ij = zeros(length(K67.rain_fill),1);
ij(ik) = 1;                         indDry = ij;           im = 1;
while im < length(ij)-13
    if ij(im) == 1
        for in = 1:12
            indDry(im+in) = 1;
        end
    end
    im = im + 1;
end
K67.LEdry = K67.LE;             K67.ETdry = K67.ET;
K67.LEdry(indDry==1) = NaN;     K67.ETdry(indDry==1) = NaN;
K67.ETdryday = K67.ETdry;       K67.ETdryday(K67.So<tresholdSo) = NaN;

K67.ETrefdry = K67.ETref;       K67.ETrefdry(indDry==1) = NaN;
K67.ETrefdryday = K67.ETrefdry; K67.ETrefdryday(K67.So<tresholdSo) = NaN;

K67.LEwet = K67.LE;             K67.ETwet = K67.ET;
K67.LEwet(indDry==0) = NaN;     K67.ETwet(indDry==0) = NaN;
K67.ETwetday = K67.ETwet;       K67.ETwetday(K67.So<tresholdSo) = NaN;

K67.ETrefwet = K67.ETref;         K67.ETrefwet(indDry==0) = NaN;
K67.ETrefwetday = K67.ETrefwet;   K67.ETrefwetday(K67.So<tresholdSo) = NaN;

ix = K67.HR.*indDry; ix(indDry==0) = NaN;
figure ('color','white');       histogram(ix,24);
ylabel('Frequency time rainfall');
xlabel('Local time rainfall K67');

%..................................................................
%%         As in Yebra
%..................................................................
A   = K67.LE+K67.H;
% A   = K67.netrad;       ix = K67.LE+K67.H;              A(isnan(A)) = ix(isnan(A));
one = K67.epsilon.*A;   two = (K67.epsilon+1).*K67.LE;	three = K67.rho_a.*10^6.*Cp.*(K67.e_sat-K67.e).*K67.Ga./K67.sigma;
Gs_OptionY = K67.LE.*K67.Ga./(one-two+three);
Gs_OptionY = Gs_OptionY.*1000;              %m s-1 to mm s-1
Gs_OptionY((Gs_OptionY<=0)|(Gs_OptionY>200)) = NaN;
Gs_OptionY((K67.So<tresholdSo)|(K67.rain>0)) = NaN;         Gs_OptionY = AM_rm_outlier(Gs_OptionY,3);
%..................................................................
one = K67.epsilon.*A;   two = (K67.epsilon+1).*K67.LE;
three = K67.rho_a.*10^6.*Cp.*(K67.e_sat-K67.e).*Ga_ustar./K67.sigma;    three = AM_rm_outlier(three,3);
Gs_OptionYustar = K67.LE.*Ga_ustar./(one-two+three);
Gs_OptionYustar = Gs_OptionYustar.*1000;    %m s-1 to mm s-1
Gs_OptionYustar((Gs_OptionYustar<=0)|(Gs_OptionYustar>1000)) = NaN;
Gs_OptionYustar((K67.So<tresholdSo)|(K67.rain>0)) = NaN;    Gs_OptionYustar = AM_rm_outlier(Gs_OptionYustar,3);

[weekthree,~,~,~,~,~,week] = AM_week2day_rs(three,K67.dateloc);
weekone = AM_week2day_rs(one,K67.dateloc);                  weektwo   = AM_week2day_rs(two,K67.dateloc);
weekRa  = AM_week2day_rs(Ra,K67.dateloc);                   weekSigma = AM_week2day_rs(K67.sigma,K67.dateloc);

%..................................................................
%%         As in Blanken & Black, 2004
%..................................................................
K67.Bowen   = K67.H./K67.LE;
Rs_OptionBB = (((K67.epsilon.*K67.Bowen)-1).*Ra)+((K67.rho_a.*Cp.*10^6.*K67.VPD_fill)./(K67.sigma.*K67.LE));
Gs_OptionBB = 1000./Rs_OptionBB;
Gs_OptionBB((Gs_OptionBB<=0)|(Gs_OptionBB>400)) = NaN;
Gs_OptionBB((K67.So<tresholdSo)|(K67.rain>0)) = NaN;        Gs_OptionBB = AM_rm_outlier(Gs_OptionBB,3);

Rs_OptionBB_ustar = (((K67.epsilon.*K67.Bowen)-1).*Ra_ustar)+((K67.rho_a.*Cp.*10^6.*K67.VPD_fill)./(K67.sigma.*K67.LE));
Gs_OptionBB_ustar = 1000./Rs_OptionBB_ustar;
Gs_OptionBB_ustar((Gs_OptionBB_ustar<=0)|(Gs_OptionBB_ustar>200)) = NaN;
Gs_OptionBB_ustar((K67.So<tresholdSo)|(K67.rain>0)) = NaN;	Gs_OptionBB_ustar = AM_rm_outlier(Gs_OptionBB_ustar,3);

% .........................................................................
K67.BowenDry = K67.H./K67.LEdry;
Rs_OptionBBdry = (((K67.epsilon.*K67.BowenDry)-1).*Ra)+((K67.rho_a.*Cp.*10^6.*K67.VPD_fill)./(K67.sigma.*K67.LEdry));
Gs_OptionBBdry = 1000./Rs_OptionBBdry;
Gs_OptionBBdry((Gs_OptionBBdry<=0)|(Gs_OptionBBdry>400)) = NaN;
Gs_OptionBBdry((K67.So<tresholdSo)|(K67.rain>0)) = NaN;     Gs_OptionBBdry = AM_rm_outlier(Gs_OptionBBdry,3);

Rs_OptionBBdry_ustar = (((K67.epsilon.*K67.BowenDry)-1).*Ra_ustar)+((K67.rho_a.*Cp.*10^6.*K67.VPD_fill)./(K67.sigma.*K67.LEdry));
Gs_OptionBBdry_ustar = 1000./Rs_OptionBBdry_ustar;
Gs_OptionBBdry_ustar((Gs_OptionBBdry_ustar<=0)|(Gs_OptionBBdry_ustar>200)) = NaN;
Gs_OptionBBdry_ustar((K67.So<tresholdSo)|(K67.rain>0)) = NaN;
Gs_OptionBBdry_ustar = AM_rm_outlier(Gs_OptionBBdry_ustar,3);

%..................................................................
% Climatological resistance
%..................................................................
K67.Ri = (K67.rho_a.*Cp.*10^6.*K67.VPD_fill)./(K67.sigma.*(K67.LE+K67.H));
%K67.Ri = (K67.rho_a.*Cp.*10^6.*K67.VPD_fill)./(K67.sigma.*(K67.netrad));
K67.Gi = 1000./K67.Ri;
K67.Gi(K67.So<tresholdSo) = NaN;        K67.Gi(K67.Gi<0) = NaN;
K67.Gi = AM_rm_outlier(K67.Gi,3);

% Equilibrium LE
K67.PTalpha = K67.LE./K67.LEref;
% .........................................................................
% Gs units mm s-1
% .........................................................................
K67.GsPMall   = Gs_OptionBB;
K67.GsPM      = Gs_OptionBBdry;
K67.GsPMustar = Gs_OptionBB_ustar;

% From Pearcy et all. Measurement of transpiration and leaf conductance
% P in kPa, tair in degC, Gs from cm/s to mol m-2s-1
% K67.GsPMallmol = K67.GsPMall.*0.446.*(273./(273+K67.tasp_fill)).*((K67.pambpa_fill./1000)./101.3);

% Gas ideal gas law Gs [m s-1] x Pa[kPa]/R [m3 Pa K-1 mol] x Ta [k]   % average sea-level pressure is 101.325 kPa
% R = 8.3144621 m3 Pa K−1 mol−1   8.3144621 m3 /1000 kPa-1 K-1 mol
% Gs mol m-2 s-1
K67.GsPMallmmol = (K67.GsPMall./1000).*(K67.pambpa_fill./1000)./(R.*(273+K67.tasp_fill));
K67.GsPMallmmol = K67.GsPMallmmol.*1000;

K67.GsPMmmol = (K67.GsPM./1000).*(K67.pambpa_fill./1000)./(R.*(273+K67.tasp_fill));
K67.GsPMmmol = K67.GsPMmmol.*1000;                        
% Typical values for irrigated trials are: 300-700 mmol m-2s-1;
% and for mildly water stressed trials are: 80-300 mmol m-2s-1

%% ........................................................................
K67.eightGsPMall = AM_eight2cycle_rs(K67.GsPMall,K67.dateloc,6,2,3);    K67.eightGsPMall = AM_fill(K67.eightGsPMall);
[K67.weekGsPMall,~,~,~,~,K67.week_stdGsPMall] = AM_week2cycle_rs(K67.GsPMall,K67.dateloc,6,2,3);  K67.weekGsPMall  = AM_fill(K67.weekGsPMall);
[K67.monthGsPMall,~,~,K67.monthGsPMallmin,K67.monthGsPMallmax] = AM_month2cycle(K67.GsPMall,K67.dateloc,6,2,2);
K67.monthGsPMall = AM_fill(K67.monthGsPMall);                           K67.monthGsPMallmax = AM_fill(K67.monthGsPMallmax);
[K67.dayGsPMall,~,x,~,~,K67.dayGsPMallmin,K67.dayGsPMallmax]   = AM_day(K67.GsPMall,K67.dateloc,0);
K67.dayGsPMall(x<6) = NaN;         K67.dayGsPMallmin(x<6) = NaN;        K67.dayGsPMallmax(x<6) = NaN;

K67.eightGsPMallmmol = AM_eight2cycle_rs(K67.GsPMallmmol,K67.dateloc,6,2,3);    K67.eightGsPMallmmol = AM_fill(K67.eightGsPMallmmol);
[K67.weekGsPMallmmol,~,~,~,~,K67.week_stdGsPMallmmol] = AM_week2cycle_rs(K67.GsPMallmmol,K67.dateloc,6,2,3);  K67.weekGsPMallmmol  = AM_fill(K67.weekGsPMallmmol);
[K67.monthGsPMallmmol,~,~,K67.monthGsPMallmmolmin,K67.monthGsPMallmmolmax] = AM_month2cycle(K67.GsPMallmmol,K67.dateloc,6,2,2);
K67.monthGsPMallmmol = AM_fill(K67.monthGsPMallmmol);                           K67.monthGsPMallmmolmax = AM_fill(K67.monthGsPMallmmolmax);
[K67.dayGsPMallmmol,~,x,~,~,K67.dayGsPMallmmolmin,K67.dayGsPMallmmolmax]   = AM_day(K67.GsPMallmmol,K67.dateloc,0);
K67.dayGsPMallmmol(x<6) = NaN;         K67.dayGsPMallmmolmin(x<6) = NaN;        K67.dayGsPMallmmolmax(x<6) = NaN;

K67.eightGsPM = AM_eight2cycle_rs(K67.GsPM,K67.dateloc,6,2,3);    K67.eightGsPM = AM_fill(K67.eightGsPM);
[K67.weekGsPM,~,~,~,~,K67.week_stdGsPM] = AM_week2cycle_rs(K67.GsPM,K67.dateloc,6,2,3);  K67.weekGsPM  = AM_fill(K67.weekGsPM);
[K67.monthGsPM,~,~,K67.monthGsPMmin,K67.monthGsPMmax] = AM_month2cycle(K67.GsPM,K67.dateloc,6,2,2);
K67.monthGsPM = AM_fill(K67.monthGsPM);                           K67.monthGsPMmax = AM_fill(K67.monthGsPMmax);
[K67.dayGsPM,~,x,~,~,K67.dayGsPMmin,K67.dayGsPMmax]   = AM_day(K67.GsPM,K67.dateloc,0);
K67.dayGsPM(x<6) = NaN;         K67.dayGsPMmin(x<6) = NaN;        K67.dayGsPMmax(x<6) = NaN;

K67.eightGsPMmmol  = AM_eight2cycle_rs(K67.GsPMmmol,K67.dateloc,6,2,2); K67.eightGsPMmmol = AM_fill(K67.eightGsPMmmol);
[K67.weekGsPMmmol,~,~,~,~,K67.week_stdGsPMmmol] = AM_week2cycle_rs(K67.GsPMmmol,K67.dateloc,tresholdCycleDay,4,1);  K67.weekGsPMmmol  = AM_fill(K67.weekGsPMmmol);
K67.monthGsPMmmol  = AM_month2cycle(K67.GsPMmmol,K67.dateloc,6,2,2);    K67.monthGsPMmmol = AM_fill(K67.monthGsPMmmol);
[K67.dayGsPMmmol,~,x] = AM_day(K67.GsPMmmol,K67.dateloc,0);             K67.dayGsPMmmol(x<6) = NaN;

[K67.cycleGsPMmmol,K67.cycle_stdGsPMmmol,~,K67.acycle_stdGsPMmmol,K67.acycleGsPMmmol] =....
    AM_cycle_composite(K67.GsPMmmol,K67.dateloc);
[K67.cycle_monthGsPMmmol,~,~,K67.cycle_month_stdGsPMmmol] = AM_day_cycle_year(K67.GsPMmmol,K67.dateloc);

[K67.cycleGsPMallmmol,K67.cycle_stdGsPMallmmol,~,K67.acycle_stdGsPMallmmol,K67.acycleGsPMallmmol] =....
    AM_cycle_composite(K67.GsPMallmmol,K67.dateloc);

% % % % K67.eightGsPMustar  = AM_eight2cycle_rs(K67.GsPMustar,K67.dateloc,6,2,2);   K67.eightGsPMustar = AM_fill(K67.eightGsPMustar);
% % % % [K67.weekGsPMustar,~,~,~,~,K67.week_stdGsPMustar] = AM_week2cycle_rs(K67.GsPMustar,K67.dateloc,6,2,2);  K67.weekGsPMustar  = AM_fill(K67.weekGsPMustar);
% % % % [K67.monthGsPMustar,~,~,K67.monthGsPMustarMIN,K67.monthGsPMustarMAX] = AM_month2cycle(K67.GsPMustar,K67.dateloc,6,2,2);
% % % % K67.monthGsPMustar = AM_fill(K67.monthGsPMustar);                           K67.monthGsPMustarMAX = AM_fill(K67.monthGsPMustarMAX);
% % % % [K67.dayGsPMallPMustar,~,x,~,~,K67.dayGsPMallPMustarMIN,K67.dayGsPMallPMustarMAX]   = AM_day(K67.GsPMustar,K67.dateloc,0);
% % % % K67.dayGsPMallPMustar(x<6) = NaN;  K67.dayGsPMallPMustarMIN(x<6) = NaN;     K67.dayGsPMallPMustarMAX(x<6) = NaN;
% % % % 

%% ........................................................................
[K67.dayETdry,~,x] = AM_day(K67.ETdry,K67.dateloc,0);               K67.dayETdry(x<6) = NaN;
K67.eightETdry = AM_eight2cycle_rs(K67.ETdry,K67.dateloc,6,1,2);    K67.eightETdry  = AM_fill(K67.eightETdry);
[K67.weekETdry,~,~,~,~,K67.week_stdETdry]  = AM_week2cycle_rs(K67.ETdry,K67.dateloc,6,2,2);     K67.weekETdry   = AM_fill(K67.weekETdry);
K67.monthETdry = AM_month2cycle(K67.ETdry,K67.dateloc,6,4,2);       K67.monthETdry  = AM_fill(K67.monthETdry);
K67.monthETdry = K67.monthETdry.*eomday(K67.monthY,K67.monthM);

[K67.dayETdryday,~,x] = AM_day(K67.ETdryday,K67.dateloc,0);              K67.dayETdryday(x<6) = NaN;
K67.eightETdryday = AM_eight2cycle_rs(K67.ETdryday,K67.dateloc,6,1,2);   K67.eightETdryday = AM_fill(K67.eightETdryday);
[K67.weekETdryday,~,~,~,~,K67.week_stdETdryday] = AM_week2cycle_rs(K67.ETdryday,K67.dateloc,6,2,2);    K67.weekETdryday  = AM_fill(K67.weekETdryday);
K67.monthETdryday = AM_month2cycle(K67.ETdryday,K67.dateloc,6,4,2);      K67.monthETdryday = AM_fill(K67.monthETdryday);
K67.monthETdryday = K67.monthETdryday.*eomday(K67.monthY,K67.monthM);

[K67.dayETrefdryday,~,x] = AM_day(K67.ETrefdryday,K67.dateloc,0);              K67.dayETrefdryday(x<6) = NaN;
K67.eightETrefdryday = AM_eight2cycle_rs(K67.ETrefdryday,K67.dateloc,6,1,2);   K67.eightETrefdryday = AM_fill(K67.eightETrefdryday);
[K67.weekETrefdryday,~,~,~,~,K67.week_stdETrefdryday] = AM_week2cycle_rs(K67.ETrefdryday,K67.dateloc,6,2,2);    K67.weekETrefdryday  = AM_fill(K67.weekETrefdryday);
K67.monthETrefdryday = AM_month2cycle(K67.ETrefdryday,K67.dateloc,6,4,2);      K67.monthETrefdryday = AM_fill(K67.monthETrefdryday);
K67.monthETrefdryday = K67.monthETrefdryday.*eomday(K67.monthY,K67.monthM);

[K67.dayETdryday24,~,x] = AM_day(K67.ETdryday24,K67.dateloc,0);              K67.dayETdryday24(x<6) = NaN;
K67.eightETdryday24 = AM_eight2cycle_rs(K67.ETdryday24,K67.dateloc,6,1,2);   K67.eightETdryday24 = AM_fill(K67.eightETdryday24);
[K67.weekETdryday24,~,~,~,~,K67.week_stdETdryday24] = AM_week2cycle_rs(K67.ETdryday24,K67.dateloc,6,2,2);    K67.weekETdryday24  = AM_fill(K67.weekETdryday24);
K67.monthETdryday24 = AM_month2cycle(K67.ETdryday24,K67.dateloc,6,4,2);      K67.monthETdryday24 = AM_fill(K67.monthETdryday24);
K67.monthETdryday24 = K67.monthETdryday24.*eomday(K67.monthY,K67.monthM);

[K67.dayETrefdryday24,~,x] = AM_day(K67.ETrefdryday24,K67.dateloc,0);              K67.dayETrefdryday24(x<6) = NaN;
K67.eightETrefdryday24 = AM_eight2cycle_rs(K67.ETrefdryday24,K67.dateloc,6,1,2);   K67.eightETrefdryday24 = AM_fill(K67.eightETrefdryday24);
[K67.weekETrefdryday24,~,~,~,~,K67.week_stdETrefdryday24] = AM_week2cycle_rs(K67.ETrefdryday24,K67.dateloc,6,2,2);    K67.weekETrefdryday24  = AM_fill(K67.weekETrefdryday24);
K67.monthETrefdryday24 = AM_month2cycle(K67.ETrefdryday24,K67.dateloc,6,4,2);      K67.monthETrefdryday24 = AM_fill(K67.monthETrefdryday24);
K67.monthETrefdryday24 = K67.monthETrefdryday24.*eomday(K67.monthY,K67.monthM);

[K67.dayETwetday24,~,x] = AM_day(K67.ETwetday24,K67.dateloc,0);              K67.dayETwetday24(x<6) = NaN;
K67.eightETwetday24 = AM_eight2cycle_rs(K67.ETwetday24,K67.dateloc,6,1,2);   K67.eightETwetday24 = AM_fill(K67.eightETwetday24);
[K67.weekETwetday24,~,~,~,~,K67.week_stdETwetday24] = AM_week2cycle_rs(K67.ETwetday24,K67.dateloc,6,2,2);    K67.weekETwetday24  = AM_fill(K67.weekETwetday24);
K67.monthETwetday24 = AM_month2cycle(K67.ETwetday24,K67.dateloc,6,4,2);      K67.monthETwetday24 = AM_fill(K67.monthETwetday24);
K67.monthETwetday24 = K67.monthETwetday24.*eomday(K67.monthY,K67.monthM);

[K67.dayETrefwetday24,~,x] = AM_day(K67.ETrefwetday24,K67.dateloc,0);              K67.dayETrefwetday24(x<6) = NaN;
K67.eightETrefwetday24 = AM_eight2cycle_rs(K67.ETrefwetday24,K67.dateloc,6,1,2);   K67.eightETrefwetday24 = AM_fill(K67.eightETrefwetday24);
[K67.weekETrefwetday24,~,~,~,~,K67.week_stdETrefwetday24] = AM_week2cycle_rs(K67.ETrefwetday24,K67.dateloc,6,2,2);    K67.weekETrefwetday24  = AM_fill(K67.weekETrefwetday24);
K67.monthETrefwetday24 = AM_month2cycle(K67.ETrefwetday24,K67.dateloc,6,4,2);      K67.monthETrefwetday24 = AM_fill(K67.monthETrefwetday24);
K67.monthETrefwetday24 = K67.monthETrefwetday24.*eomday(K67.monthY,K67.monthM);

[K67.dayETwetday,~,x] = AM_day(K67.ETwetday,K67.dateloc,0);              K67.dayETwetday(x<6) = NaN;
K67.eightETwetday = AM_eight2cycle_rs(K67.ETwetday,K67.dateloc,6,1,2);   K67.eightETwetday = AM_fill(K67.eightETwetday);
[K67.weekETwetday,~,~,~,~,K67.week_stdETwetday] = AM_week2cycle_rs(K67.ETwetday,K67.dateloc,6,2,2);    K67.weekETwetday  = AM_fill(K67.weekETwetday);
K67.monthETwetday = AM_month2cycle(K67.ETwetday,K67.dateloc,6,4,2);      K67.monthETwetday = AM_fill(K67.monthETwetday);
K67.monthETwetday = K67.monthETwetday.*eomday(K67.monthY,K67.monthM);

[K67.dayETrefwetday,~,x] = AM_day(K67.ETrefwetday,K67.dateloc,0);              K67.dayETrefwetday(x<6) = NaN;
K67.eightETrefwetday = AM_eight2cycle_rs(K67.ETrefwetday,K67.dateloc,6,1,2);   K67.eightETrefwetday = AM_fill(K67.eightETrefwetday);
[K67.weekETrefwetday,~,~,~,~,K67.week_stdETrefwetday] = AM_week2cycle_rs(K67.ETrefwetday,K67.dateloc,6,2,2);    K67.weekETrefwetday  = AM_fill(K67.weekETrefwetday);
K67.monthETrefwetday = AM_month2cycle(K67.ETrefwetday,K67.dateloc,6,4,2);      K67.monthETrefwetday = AM_fill(K67.monthETrefwetday);
K67.monthETrefwetday = K67.monthETrefwetday.*eomday(K67.monthY,K67.monthM);

%% .........................................................................
x = K67.Gi;       x(K67.So<tresholdSo) = NaN;
[K67.weekGiday,~,~,~,~,K67.week_stdGiday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
% K67.weekGiday     = AM_rm_outlier(K67.weekGiday,3);
[~,K67.weekGiday] = AM_spike_detection(K67.weekGiday,5);
K67.weekGiday     = AM_fill(K67.weekGiday);
[K67.eightGiday]  = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
% K67.eightGiday    = AM_rm_outlier(K67.eightGiday,3);
[~,K67.eightGiday]= AM_spike_detection(K67.eightGiday,5);
K67.eightGiday    = AM_fill(K67.eightGiday);
[K67.monthGiday]  = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
% K67.monthGiday    = AM_rm_outlier(K67.monthGiday,3);
[~,K67.eightGiday] = AM_spike_detection(K67.eightGiday,3);
K67.monthGiday = AM_fill(K67.monthGiday);

x = K67.LEdry;      x(K67.So<tresholdSo) = NaN;
[K67.weekLEdry,~,~,~,~,K67.week_stdLEdry] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekLEdry] = AM_spike_detection(K67.weekLEdry,5);
% K67.weekLEdry   = AM_rm_outlier(K67.weekLEdry,3);
K67.weekLEdry     = AM_fill(K67.weekLEdry);
[K67.eightLEdry]  = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,2,2);
[~,K67.eightLEdry]= AM_spike_detection(K67.eightLEdry,5);
% K67.eightLEdry  = AM_rm_outlier(K67.eightLEdry,3);
K67.eightLEdry    = AM_fill(K67.eightLEdry);
[K67.monthLEdry] = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,2,2);
[~,K67.monthLEdry] = AM_spike_detection(K67.monthLEdry,5);
% K67.monthLEdry   = AM_rm_outlier(K67.monthLEdry,3);
K67.monthLEdry   = AM_fill(K67.monthLEdry);

x = K67.PTalpha;      x(K67.So<tresholdSo) = NaN;
[K67.weekPTalphaday,~,~,~,~,K67.week_stdPTalphaday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekPTalphaday] = AM_spike_detection(K67.weekPTalphaday,5);
% K67.weekPTalphaday    = AM_rm_outlier(K67.weekPTalphaday,3);
K67.weekPTalphaday    = AM_fill(K67.weekPTalphaday);
[K67.eightPTalphaday] = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightPTalphaday] = AM_spike_detection(K67.eightPTalphaday,5);
% K67.eightPTalphaday   = AM_rm_outlier(K67.eightPTalphaday,3);
K67.eightPTalphaday   = AM_fill(K67.eightPTalphaday);
[K67.monthPTalphaday] = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthPTalphaday] = AM_spike_detection(K67.monthPTalphaday,5);
% K67.monthPTalphaday   = AM_rm_outlier(K67.monthPTalphaday,3);
K67.monthPTalphaday   = AM_fill(K67.monthPTalphaday);

%%
x = K67.RHec;
[K67.dayRH,~,y] = AM_day(x,K67.dateloc,0);        K67.dayRH(y<6) = NaN;
[K67.weekRH,~,~,~,~,K67.week_stdRH] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekRH]  = AM_spike_detection(K67.weekRH,5);
% K67.weekRH    = AM_rm_outlier(K67.weekRH,3);
K67.weekRH      = AM_fill(K67.weekRH);
[K67.eightRH]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycle,4,2);
[~,K67.eightRH] = AM_spike_detection(K67.eightRH,5);
% K67.eightRH   = AM_rm_outlier(K67.eightRH,3);
K67.eightRH     = AM_fill(K67.eightRH);
[K67.monthRH]   = AM_month2cycle(x,K67.dateloc,tresholdCycle,4,2);
[~,K67.monthRH] = AM_spike_detection(K67.monthRH,5);
% K67.monthRH   = AM_rm_outlier(K67.monthRH,3);
K67.monthRH  = AM_fill(K67.monthRH);

x = K67.RHec;      x(K67.So<tresholdSo) = NaN;
[K67.dayRHday,~,y] = AM_day(x,K67.dateloc,0);        K67.dayRHday(y<6) = NaN;
[K67.weekRHday,~,~,~,~,K67.week_stdRHday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.weekRHday]  = AM_spike_detection(K67.weekRHday,5);
% K67.weekRHday    = AM_rm_outlier(K67.weekRHday,3);
K67.weekRHday      = AM_fill(K67.weekRHday);
[K67.eightRHday]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.eightRHday] = AM_spike_detection(K67.eightRHday,5);
% K67.eightRHday   = AM_rm_outlier(K67.eightRHday,3);
K67.eightRHday   = AM_fill(K67.eightRHday);
[K67.monthRHday] = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.monthRHday] = AM_spike_detection(K67.monthRHday,5);
% K67.monthRHday   = AM_rm_outlier(K67.monthRHday,3);
K67.monthRHday     = AM_fill(K67.monthRHday);
[K67.month_avgRHday,K67.month_avgstdRHday] = AM_month_avg(K67.monthRHday,K67.month);

%% ........................................................................
x = K67.VPD;      x(K67.So<tresholdSo) = NaN;
[K67.dayVPDday,~,y] = AM_day(x,K67.dateloc,0);        K67.dayVPDday(y<6) = NaN;
[K67.weekVPDday,~,~,~,~,K67.week_stdVPDday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekVPDday] = AM_spike_detection(K67.weekVPDday,5);
% K67.weekVPDday   = AM_rm_outlier(K67.weekVPDday,3);
K67.weekVPDday     = AM_fill(K67.weekVPDday);
K67.week_stdVPDday(isnan(K67.weekVPDday)) = NaN;
[K67.eightVPDday]  = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightVPDday] = AM_spike_detection(K67.eightVPDday,5);
% K67.eightVPDday   = AM_rm_outlier(K67.eightVPDday,3);
K67.eightVPDday     = AM_fill(K67.eightVPDday);
[K67.monthVPDday]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthVPDday] = AM_spike_detection(K67.monthVPDday,5);
% K67.monthVPDday   = AM_rm_outlier(K67.monthVPDday,3);
K67.monthVPDday  = AM_fill(K67.monthVPDday);

figure; subplot(2,2,1);plot(week,weekone);  datetick('x'); ylabel('one');  subplot(2,2,2);plot(week,weektwo);    datetick('x');   ylabel('two');
subplot(2,2,3);plot(week,K67.weekGsPMall); hold on; plot(week,K67.weekGsPM); datetick('x'); ylabel('Gs');   
subplot(2,2,4);plot(week,weekthree);    datetick('x');   ylabel('three');

%%
figure('color','white');
color_mtx=jet(length(K67.year));
for ik = 1:length(K67.year)-1
    ind = find((K67.week>=datenum(K67.year(ik),1,1))&(K67.week<datenum(K67.year(ik)+1,1,1)));
    if length(ind)>1
        plot(week_avg,K67.weekGsPMall(ind),'color',color_mtx(ik,:));   hold on;
        plot(week_avg(2),10+ik,'.','color',color_mtx(ik,:));
        text (week_avg(3),10+ik,num2str(K67.year(ik)),'FontSize',9);
    end
end

% [K67.h_weekGs] = AM_week2cycle_rs(K67.h_gs_OptionXustar,K67.h_dateloc,6,2,2);
% K67.h_week_avgGs = AM_week2_avg(K67.h_weekGs,K67.h_week);
% K67.h_weekGs_smooth  =  smooth(K67.h_week_avgGs,c);     K67.h_weekGs_smooth(isnan(K67.h_week_avgGs)) = NaN;
% hold on; plot(K67.h_weekGs_smooth);

%..........................................................................
%% Bowen and Gs
%..........................................................................
[K67.eightLE,~,~,K67.eightLE_min,K67.eightLE_max,~,K67.eight] = AM_eight2cycle_rs(K67.LE,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightLE] = AM_spike_detection(K67.eightLE,5);
% K67.eightLE   = AM_rm_outlier(K67.eightLE,5);
K67.eightLE = AM_fill(K67.eightLE);
[K67.weekLE,~,~,~,~,K67.week_stdLE] = AM_week2cycle_rs(K67.LE,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekLE] = AM_spike_detection(K67.weekLE,5);
% K67.weekLE   = AM_rm_outlier(K67.weekLE,3);
K67.weekLE     = AM_fill(K67.weekLE);
[K67.week_avgLE,~,K67.week_avgstdLE] = AM_week2_avg(K67.weekLE,K67.week);
[K67.monthLE,~,~,~,~,K67.monthstg_LE] = AM_month2cycle(K67.LE,K67.dateloc,tresholdCycle,4,2);
[~,K67.monthLE] = AM_spike_detection(K67.monthLE,5);
% K67.monthLE   = AM_rm_outlier(K67.monthLE,5);
K67.monthLE     = AM_fill(K67.monthLE);

K67.LEday = K67.LE;            K67.LEday(K67.So<tresholdSo)  = NaN;
[K67.dayLEday,~,x,~,K67.day,K67.dayLEdayMIN,K67.dayLEdayMAX] = AM_day(K67.LEday,K67.dateloc,0);
K67.dayLEday(x<6) = NaN;       K67.dayLEdayMIN(x<6) = NaN;   K67.dayLEdayMAX(x<6) = NaN;
[K67.weekLEday,~,~,~,~,K67.week_stdLEday] = AM_week2cycle_rs(K67.LEday,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.weekLEday] = AM_spike_detection(K67.weekLEday,5);
% K67.weekLEday   = AM_rm_outlier(K67.weekLEday,6);
K67.weekLEday     = AM_fill_plus(K67.weekLEday);
[K67.week_avgLEday,~,K67.week_avgstdLEday] = AM_week2_avg(K67.weekLEday,K67.week);
[K67.eightLEday]   = AM_eight2cycle_rs(K67.LEday,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.eightLEday] = AM_spike_detection(K67.eightLEday,5);
% K67.eightLEday   = AM_rm_outlier(K67.eightLEday,6);
K67.eightLEday     = AM_fill(K67.eightLEday);
[K67.monthLEday,~,~,~,K67.monthLEdayMAX] = AM_month2cycle(K67.LEday,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.monthLEday] = AM_spike_detection(K67.monthLEday,5);
% K67.monthLEday    = AM_rm_outlier(K67.monthLEday,6);
K67.monthLEday    = AM_fill(K67.monthLEday);
% K67.monthLEday    = K67.monthLEday.*eomday(K67.monthY,K67.monthM);

%%
[K67.monthH,~,~,~,~,K67.monthstg_H] = AM_month2cycle(K67.H,K67.dateloc,tresholdCycle,4,2);
[~,K67.monthH] = AM_spike_detection(K67.monthH,5);
% K67.monthH   = AM_rm_outlier(K67.monthH,5);
K67.monthH     = AM_fill(K67.monthH);
[K67.weekH,~,~,~,~,K67.week_stdH] = AM_week2cycle_rs(K67.H,K67.dateloc,tresholdCycle,4,2);
[K67.week_avgH,~,K67.week_avgstdH] = AM_week2_avg(K67.weekH,K67.week);
[~,K67.weekH]  = AM_spike_detection(K67.weekH,5);
% K67.weekH = AM_rm_outlier(K67.weekH,3);
K67.weekH      = AM_fill(K67.weekH);
[K67.eightH,~,~,~,~,K67.eightH_std] = AM_eight2cycle_rs(K67.H,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightH] = AM_spike_detection(K67.eightH,5);
% K67.eightH   = AM_rm_outlier(K67.eightH,5);
K67.eightH     = AM_fill(K67.eightH);

K67.weekBowen  = AM_fill(K67.weekH./K67.weekLE);     K67.weekEF  = AM_fill(K67.weekLE./(K67.weekLE+K67.weekH));
K67.eightBowen = AM_fill(K67.eightH./K67.eightLE);   K67.eightEF = AM_fill(K67.eightLE./(K67.eightLE+K67.eightH));
K67.monthBowen = AM_fill(K67.monthH./K67.monthLE);   K67.monthEF = AM_fill(K67.monthLE./(K67.monthLE+K67.monthH));

K67.Bowenday = AM_fill(K67.H./K67.LE);                  K67.Bowenday(K67.So<tresholdSo) = NaN;
K67.Bowenday((K67.Bowenday<0)|(K67.Bowenday>1)) = NaN;
K67.Bowenday = AM_rm_outlier(K67.Bowenday,3);
K67.eightBowenday = AM_eight2cycle_rs(K67.Bowenday,K67.dateloc,tresholdCycleDay,4,2);	K67.eightBowenday = AM_fill(K67.eightBowenday);
[K67.weekBowenday,~,~,~,~,K67.week_stdBowenday] = AM_week2cycle_rs(K67.Bowenday,K67.dateloc,tresholdCycleDay,4,2);    K67.weekBowenday  = AM_fill(K67.weekBowenday);
K67.monthBowenday = AM_month2cycle(K67.Bowenday,K67.dateloc,tresholdCycleDay,4,2);      K67.monthBowenday = AM_fill(K67.monthBowenday);

K67.EFdaydry = AM_fill(K67.LEdry./(K67.LE+K67.H));      K67.EFdaydry(K67.So<tresholdSo) = NaN;
K67.EFdaydry((K67.EFdaydry<0)|(K67.EFdaydry>1)) = NaN;
K67.EFdaydry = AM_rm_outlier(K67.EFdaydry,3);
K67.eightEFdaydry = AM_eight2cycle_rs(K67.EFdaydry,K67.dateloc,tresholdCycleDay,2,2);	K67.eightEFdaydry = AM_fill(K67.eightEFdaydry);
[K67.weekEFdaydry,~,~,~,~,K67.week_stdEFdaydry] = AM_week2cycle_rs(K67.EFdaydry,K67.dateloc,tresholdCycleDay,2,2);    K67.weekEFdaydry  = AM_fill(K67.weekEFdaydry);
K67.monthEFdaydry = AM_month2cycle(K67.EFdaydry,K67.dateloc,tresholdCycleDay,2,2);      K67.monthEFdaydry = AM_fill(K67.monthEFdaydry);

%..........................................................................
%% Storage heat
%..........................................................................
[Sensor,HeightSensorTable] = ...
    textread('/home/ncoupe/tapajos/R.codes/heights.txt','%s %f',...
    'delimiter', '\t','headerlines',1);
SensorHeight=NaN(8,1);
for ik=1:8
    for ij = 1:length(Sensor)
        ix = Sensor{ij};    iy = strcat('Tair_', num2str(ik));
        if strcmp(ix,iy)
            SensorHeight(ik)=HeightSensorTable(ij);
        end
    end
end
K67.tair8clr = K67.tair8;             	K67.tair8clr(K67.dateloc>datenum(2016,9,15)) = NaN;
K67.tair5    = AM_rm_outlier(K67.tair5,4);
K67.tair6(88340:89130)  = NaN;          K67.tair6(K67.tair6<20) = NaN;              K67.tair6 = AM_rm_outlier(K67.tair6,4);
K67.tair7(77490:100000) = NaN;          K67.tair7 = AM_rm_outlier(K67.tair7,4);
DataTable = [K67.tair1,K67.tair2,K67.tair3,K67.tair4,K67.tair5,K67.tair6,K67.tair7,K67.tair8clr];
%K67.Sheat = AM_Sco2_calc(K67.dateloc,K67.rho_a_fill,DataTable,SensorHeight,'H');

%..........................................................................
%% GEP at saturated light
%..........................................................................
K67.radiation   = K67.par;                    
K67.netrad_plus = K67.netrad;
K67.netrad_plus(isnan(K67.netrad_plus)) = K67.H(isnan(K67.netrad_plus))+K67.LE(isnan(K67.netrad_plus));
K67.radiation_plus  = K67.radiation;
K67.radiation_plus  = AM_fill_plus(K67.radiation_plus);
[~,~,K67.radiation_plus] = AM_compute_regression(K67.netrad_plus,K67.radiation,[.4 .4 .4],0);
K67.radiation_plus((K67.radiation_plus<0)&(K67.So<10)) = 0;
K67.radiation_plus(K67.radiation_plus<0) = NaN;

ix = K67.radiation_plus;    ix((K67.So<10)|(ix<5)) = NaN;
[iy] = AM_hour_year(K67.dateloc,ix);
PARdayAVG      = 500;                % nanmean(iy);

K67.weekGEPsat = K67.week.*NaN;     K67.weekTairsat = K67.week.*NaN;    K67.weekVPDsat = K67.week.*NaN;
for im=1:length(K67.week)
    ik = find((K67.dateloc>=K67.week(im))&(K67.dateloc<(K67.week(im)+16))&...
        (K67.So>100)&...
        (~isnan(K67.gee))&(-K67.gee>=0)&(~isnan(K67.tasp)));
    if ~isempty(ik)
        [~,LUE,GEPsat,R2,~,Pc,re_LUE,PARsat]=AM_LUE_toolbox(K67.radiation_plus(ik),-K67.gee(ik),0,PARdayAVG);
        PARsat(isempty(PARsat)) = PARdayAVG;            PARsat(PARsat>=1000) = PARdayAVG;
        display(PARsat)
        PARsat = PARdayAVG;
        ik = find((K67.dateloc>=K67.week(im))&(K67.dateloc<(K67.week(im)+16))&...
            (K67.So>100)&...
            (~isnan(K67.gee))&(-K67.gee>=0)&(~isnan(K67.tasp))&...
            (K67.So>100)&(K67.radiation_plus>=PARsat));
        K67.weekGEPsat(im)  = nanmean(-K67.gee(ik));
        K67.weekTairsat(im) = nanmean(K67.tasp(ik));
        K67.weekVPDsat(im)  = nanmean(K67.VPD(ik));
    end
end
K67.monthGEPsat = interp1(K67.week,K67.weekGEPsat,K67.month);


.........................................................................
%% VPD vs GEPsat
.........................................................................
    figure('color','white');  hold on;
K67.weekGEPsatVPDslope = K67.week.*NaN;     K67.weekGEPsatVPDintercept = K67.week.*NaN;
for im=1:length(K67.week)
    ik = find((K67.dateloc>=K67.week(im))&(K67.dateloc<(K67.week(im)+16))&...
        (K67.So>100)&(K67.radiation_plus>=1000)&...
        (~isnan(K67.gee))&(-K67.gee>=0)&(~isnan(K67.VPD)));
    if ~isempty(ik)
        [p,r2,~,~,sm,sb] = AM_compute_regression(K67.VPD(ik),-K67.gee(ik).*1000./K67.radiation_plus(ik),[0.4,0.4,0.4],1);
        y_model=polyval(p,K67.VPD(ik));   hold on;
        plot(K67.VPD(ik),y_model,'-','color',[0.4,0.4,0.4]);
        xlim([0 3]);      ylim([0 50]);
        text (1,40,{['std slope:' num2str(sm(1),'%3.4f')];['std int.:' num2str(sb(1),'%3.4f')];...
            ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10);         box on;
        % close
        K67.weekGEPsatVPDslope(im)     = p(1);
        K67.weekGEPsatVPDintercept(im) = p(2);
    end
end

K67.weekGEPsatVPDslope((K67.weekGEPsatVPDslope<-100)|(K67.weekGEPsatVPDslope>0)) = NaN; 
K67.weekGEPsatVPDintercept((K67.weekGEPsatVPDintercept<0)|(K67.weekGEPsatVPDintercept>100)) = NaN;

K67.weekGEPsatVPDslope     = AM_rm_outlier(K67.weekGEPsatVPDslope,3);       
K67.weekGEPsatVPDintercept = AM_rm_outlier(K67.weekGEPsatVPDintercept,3);
K67.weekGEPsatVPDslope     = AM_fill(K67.weekGEPsatVPDslope);
K67.weekGEPsatVPDintercept = AM_fill(K67.weekGEPsatVPDintercept);

ix = smooth(K67.weekGEPsatVPDslope);    ix(isnan(K67.weekGEPsatVPDslope)) = NaN;    K67.weekGEPsatVPDslope = ix;
ix = smooth(K67.weekGEPsatVPDintercept);    ix(isnan(K67.weekGEPsatVPDintercept)) = NaN;    K67.weekGEPsatVPDintercept = ix;

[K67.week_avgGEPsatVPDslope,~,K67.week_avgstdGEPsatVPDslope]     = AM_week2_avg(K67.weekGEPsatVPDslope,K67.week); 
[K67.week_avgGEPsatVPDintercept,~,K67.week_avgGEPsatVPDintercept] = AM_week2_avg(K67.weekGEPsatVPDintercept,K67.week); 

.........................................................................
%% Eight
% .........................................................................
[K67.eightNetRad,~,~,~,~,K67.eightNetRad_std] = AM_eight2cycle_rs(K67.netrad,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightNetRad] = AM_spike_detection(K67.eightNetRad,5);
% K67.eightNetRad     = AM_rm_outlier(K67.eightNetRad,5);
K67.eightNetRad = AM_fill(K67.eightNetRad);

[K67.weekNetRad,~,~,~,~,K67.week_stdNetRad] = AM_week2cycle_rs(K67.netrad,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekNetRad] = AM_spike_detection(K67.weekNetRad,5);
% K67.weekNetRad     = AM_rm_outlier(K67.weekNetRad,5);
K67.weekNetRad = AM_fill(K67.weekNetRad);
[K67.week_avgNetRad,~,K67.week_avgstdNetRad] = AM_week2_avg(K67.weekNetRad,K67.week);

[K67.eightTotRad]   = AM_eight2cycle_rs(K67.totsun,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightTotRad] = AM_spike_detection(K67.eightTotRad,5);
% K67.eightTotRad     = AM_rm_outlier(K67.eightTotRad,5);
K67.eightTotRad = AM_fill(K67.eightTotRad);

[K67.weekTotRad,~,~,~,~,K67.week_stdTotRad] = AM_week2cycle_rs(K67.totsun,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekTotRad] = AM_spike_detection(K67.weekTotRad,5);
% K67.weekTotRad     = AM_rm_outlier(K67.weekTotRad,5);
K67.weekTotRad = AM_fill(K67.weekTotRad);

[K67.eightDifRad,~,~,~,~,K67.eightDifRad_min] = AM_eight2cycle_rs(K67.difsun,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightDifRad] = AM_spike_detection(K67.eightDifRad,5);
% K67.eightDifRad   = AM_rm_outlier(K67.eightDifRad,5);
K67.eightDifRad  = AM_fill(K67.eightDifRad);

[K67.weekDifRad,~,~,~,~,K67.week_stdDifRad] = AM_week2cycle_rs(K67.difsun,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekDifRad] = AM_spike_detection(K67.weekDifRad,5);
% K67.weekDifRad     = AM_rm_outlier(K67.weekDifRad,5);
K67.weekDifRad = AM_fill(K67.weekDifRad);

[K67.eightPAR,~,~,~,~,~,K67.eight] = AM_eight2cycle_rs(K67.par,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightPAR] = AM_spike_detection(K67.eightPAR,5);
[K67.weekPAR,~,~,~,~,K67.week_stdPAR] = AM_week2cycle_rs(K67.par,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekPAR]  = AM_spike_detection(K67.weekPAR,5);
[K67.week_avgPAR,~,K67.week_avgstdPAR] = AM_week2_avg(K67.weekPAR,K67.week);
[K67.monthPAR,~,~,~,~,K67.monthPAR_std,K67.month] = AM_month2cycle(K67.par,K67.dateloc,tresholdCycle,4,2);
[~,K67.monthPAR] = AM_spike_detection(K67.monthPAR,5);
% K67.monthPAR   = AM_rm_outlier(K67.monthPAR,3);
K67.monthPAR     = AM_fill(K67.monthPAR);

[K67.eightTasp]   = AM_eight2cycle_rs(K67.tasp,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightTasp] = AM_spike_detection(K67.eightTasp,5);
% K67.eightTasp   = AM_rm_outlier(K67.eightTasp,5);
K67.eightTasp     = AM_fill(K67.eightTasp);

[K67.weekTasp,~,~,~,~,K67.week_stdTasp] = AM_week2cycle_rs(K67.tasp,K67.dateloc,tresholdCycle,2,2);
[~,K67.weekTasp] = AM_spike_detection(K67.weekTasp,5);
% K67.weekTasp   = AM_rm_outlier(K67.weekTasp,5);
K67.weekTasp     = AM_fill(K67.weekTasp);

[K67.eightH2O]   = AM_eight2cycle_rs(K67.h2ommolm,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightH2O] = AM_spike_detection(K67.eightH2O,5);
K67.eightH2O     = AM_fill(K67.eightH2O);

[K67.weekH2O,~,~,~,~,K67.week_stdH2O] = AM_week2cycle_rs(K67.h2ommolm,K67.dateloc,tresholdCycle,2,2);
[~,K67.weekH2O] = AM_spike_detection(K67.weekH2O,5);
K67.weekH2O     = AM_fill(K67.weekH2O);

[K67.eightPambpa]   = AM_eight2cycle_rs(K67.pambpa,K67.dateloc,tresholdCycle,2,2);
[~,K67.eightPambpa] = AM_spike_detection(K67.eightPambpa,5);
% K67.eightPambpa   = AM_rm_outlier(K67.eightPambpa,5);
K67.eightPambpa     = AM_fill(K67.eightPambpa);

% .........................................................................
%% Daytime
% .........................................................................
x = K67.ustar;      x(K67.So<tresholdSo) = NaN;
[K67.weekUstarday,~,~,~,~,K67.week_stdUstarday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekUstarday] = AM_spike_detection(K67.weekUstarday,5);
% K67.weekUstarday   = AM_rm_outlier(K67.weekUstarday,5);
K67.weekUstarday     = AM_fill(K67.weekUstarday);

x = K67.netrad;         x(K67.So<tresholdSo) = NaN;
[K67.weekNetRadday,~,~,~,~,K67.week_stdNetRadday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
% [~,K67.weekNetRadday] = AM_spike_detection(K67.weekNetRadday,5);
% K67.weekNetRadday = AM_rm_outlier(K67.weekNetRadday,5);
K67.weekNetRadday   = AM_fill_plus(K67.weekNetRadday);
% ix = smooth(K67.weekNetRadday,c);       ix(isnan(K67.weekNetRadday)) = NaN; K67.weekNetRadday_smooth = ix;
disp('Check this data time step');
% K67.weekNetRadday(K67.weekNetRadday<150) = NaN;
[K67.week_avgNetRadday,~,K67.week_avgstdNetRadday] = AM_week2_avg(K67.weekNetRadday,K67.week);

[K67.eightNetRadday]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightNetRadday] = AM_spike_detection(K67.eightNetRadday,5);
% K67.eightNetRadday   = AM_rm_outlier(K67.eightNetRadday,5);
K67.eightNetRadday     = AM_fill(K67.eightNetRadday);
% ix = smooth(K67.eightNetRadday,c);       ix(isnan(K67.eightNetRadday)) = NaN; K67.eightNetRadday_smooth = ix;
[K67.eight_avgNetRadday,~,K67.eight_stdNetRadday] = AM_eight_avg(K67.eightNetRadday,K67.eight);

[K67.monthNetRadday]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthNetRadday] = AM_spike_detection(K67.monthNetRadday,5);
% K67.monthNetRadday   = AM_rm_outlier(K67.monthNetRadday,5);
K67.monthNetRadday     = AM_fill(K67.monthNetRadday);
% ix = smooth(K67.monthNetRadday,c);       ix(isnan(K67.monthNetRadday)) = NaN; K67.monthNetRadday_smooth = ix;
K67.month_avgNetRadday = AM_month_avg(K67.monthNetRadday,K67.month);

K67.ETday = K67.ET;              K67.ETday(K67.So<tresholdSo) = NaN;
[K67.dayETday,~,x,~,~,K67.dayETdayMIN,K67.dayETdayMAX]  = AM_day(K67.ETday,K67.dateloc,0);
K67.dayETday(x<6) = NaN;         K67.dayETdayMIN(x<6) = NaN;    K67.dayETdayMAX(x<6) = NaN;
[K67.weekETday,~,~,~,~,K67.week_stdETday] = AM_week2cycle_rs(K67.ETday,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekETday] = AM_spike_detection(K67.weekETday,5);
% K67.weekETday     = AM_rm_outlier(K67.weekETday,6);
K67.weekETday     = AM_fill(K67.weekETday);
[K67.week_avgETday,~,K67.week_avgstdETday] = AM_week2_avg(K67.weekETday,K67.week);

[K67.eightETday]   = AM_eight2cycle_rs(K67.ETday,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightETday] = AM_spike_detection(K67.eightETday,5);
% K67.eightETday    = AM_rm_outlier(K67.eightETday,6);
K67.eightETday = AM_fill(K67.eightETday);

[K67.monthETday,~,~,~,K67.monthETdayMAX] = AM_month2cycle(K67.ETday,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthETday] = AM_spike_detection(K67.monthETday,5);
% K67.monthETday    = AM_rm_outlier(K67.monthETday,6);
K67.monthETday = AM_fill(K67.monthETday);
K67.monthETday    = K67.monthETday.*eomday(K67.monthY,K67.monthM);
[~,K67.monthETday] = AM_spike_detection(K67.monthETday,5);
% K67.monthETdayMAX = AM_rm_outlier(K67.monthETdayMAX,6);
K67.monthETdayMAX = AM_fill(K67.monthETdayMAX);
K67.monthETdayMAX = K67.monthETdayMAX.*eomday(K67.monthY,K67.monthM);

%%
K67.Hday = K67.H;              K67.Hday(K67.So<tresholdSo) = NaN;
[K67.dayHday,~,x,~,~,K67.dayHdayMIN,K67.dayHdayMAX]  = AM_day(K67.Hday,K67.dateloc,0);
K67.dayHday(x<6) = NaN;        K67.dayHdayMIN(x<6) = NaN;    K67.dayHdayMAX(x<6) = NaN;
[K67.weekHday,~,~,~,~,K67.week_stdHday] = AM_week2cycle_rs(K67.Hday,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekHday] = AM_spike_detection(K67.weekHday,5);
% K67.weekHday   = AM_rm_outlier(K67.weekHday,6);
K67.weekHday     = AM_fill(K67.weekHday);
[K67.week_avgHday,~,K67.week_avgstdHday] = AM_week2_avg(K67.weekHday,K67.week);
[K67.eightHday]   = AM_eight2cycle_rs(K67.Hday,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightHday] = AM_spike_detection(K67.eightHday,5);
% K67.eightHday   = AM_rm_outlier(K67.eightHday,6);
K67.eightHday     = AM_fill(K67.eightHday);
[K67.monthHday,~,~,~,K67.monthHdayMAX] = AM_month2cycle(K67.Hday,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightGiday] = AM_spike_detection(K67.eightGiday,5);
% K67.monthHday    = AM_rm_outlier(K67.monthHday,6);
K67.monthHday = AM_fill(K67.monthHday);
% K67.monthHday    = K67.monthHday.*eomday(K67.monthY,K67.monthM);

K67.EFday = AM_fill(K67.LE./(K67.LE+K67.H));         K67.EFday(K67.So<tresholdSo) = NaN;
K67.EFday((K67.EFday<0)|(K67.EFday>1)) = NaN;
K67.EFday = AM_rm_outlier(K67.EFday,3);
K67.eightEFday = AM_eight2cycle_rs(K67.EFday,K67.dateloc,tresholdCycleDay,4,2);	K67.eightEFday = AM_fill(K67.eightEFday);
[K67.weekEFday,~,~,~,~,K67.week_stdEFday] = AM_week2cycle_rs(K67.EFday,K67.dateloc,tresholdCycleDay,4,2);     K67.weekEFday  = AM_fill(K67.weekEFday);
K67.monthEFday = AM_month2cycle(K67.EFday,K67.dateloc,tresholdCycleDay,4,2);       K67.monthEFday = AM_fill(K67.monthEFday);

%..........................................................................
%% Seasonality
%..........................................................................
[K67.week_avgTasp,~,K67.week_avgstdTasp]   = AM_week2_avg(K67.weekTasp,K67.week);
[K67.week_avgUstarday,~,K67.week_avgstdUstarday] = AM_week2_avg(K67.weekUstarday,K67.week);
[K67.week_avgBowen,~,K67.week_avgstdBowen] = AM_week2_avg(K67.weekBowen,K67.week);
[K67.week_avgEF,~,K67.week_avgstdEF] = AM_week2_avg(K67.weekEF,K67.week);
[K67.week_avgET,~,K67.week_avgstdET] = AM_week2_avg(K67.weekET,K67.week);
[K67.week_avgETref,~,K67.week_avgstdETref] = AM_week2_avg(K67.weekETref,K67.week);
[K67.week_avgGsPM,~,K67.week_avgstdGsPM]       = AM_week2_avg(K67.weekGsPM,K67.week);
[K67.week_avgGsPMall,~,K67.week_avgstdGsPMall] = AM_week2_avg(K67.weekGsPMall,K67.week);
[K67.week_avgETdry,~,K67.week_avgstdETdry]       = AM_week2_avg(K67.weekETdry,K67.week);
[K67.week_avgETdryday,~,K67.week_avgstdETdryday] = AM_week2_avg(K67.weekETdryday,K67.week);
[K67.week_avgETref,~,K67.week_avgstdETref]       = AM_week2_avg(K67.weekETref,K67.week);
[K67.week_avgETrefday,~,K67.week_avgstdETrefday] = AM_week2_avg(K67.weekETrefday,K67.week);
[K67.week_avgGsPMmmol,~,K67.week_avgstdGsPMmmol] = AM_week2_avg(K67.weekGsPMmmol,K67.week);
[K67.week_avgGsPMallmmol,~,K67.week_avgstdGsPMallmmol] = AM_week2_avg(K67.weekGsPMallmmol,K67.week);
[K67.week_avgGiday,~,K67.week_avgstdGiday] = AM_week2_avg(K67.weekGiday,K67.week);
[K67.week_avgLEdry,~,K67.week_avgstdLEdry] = AM_week2_avg(K67.weekLEdry,K67.week);
[K67.week_avgPTalphaday,~,K67.week_avgstdPTalphaday] = AM_week2_avg(K67.weekPTalphaday,K67.week);
[K67.week_avgVPDday,~,K67.week_avgstdVPDday] = AM_week2_avg(K67.weekVPDday,K67.week);
[K67.week_avgEFday,~,K67.week_avgstdEFday]   = AM_week2_avg(K67.weekEFday,K67.week);
[K67.week_avgEFdaydry,~,K67.week_avgstdEFdaydry] = AM_week2_avg(K67.weekEFdaydry,K67.week);
[K67.week_avgBowenday,~,K67.week_avgstdBowenday] = AM_week2_avg(K67.weekBowenday,K67.week);
[K67.week_avgH2O,~,K67.week_avgstdH2O]   = AM_week2_avg(K67.weekH2O,K67.week);

% .........................................................................
%% Eight seasonality
% .........................................................................
[K67.eight_avgBowen,~,K67.eight_stdBowen] = AM_eight_avg(K67.eightBowen,K67.eight);
[K67.eight_avgEF,~,K67.eight_stdEF]       = AM_eight_avg(K67.eightEF,K67.eight);
[K67.eight_avgET,~,K67.eight_stdET]       = AM_eight_avg(K67.eightET,K67.eight);
[K67.eight_avgGs,~,K67.eight_stdGs]       = AM_eight_avg(K67.eightGsPMall,K67.eight);
[K67.eight_avgGsPM,~,K67.eight_stdGsPM]   = AM_eight_avg(K67.eightGsPM,K67.eight);
[K67.month_avgGsPMall,K67.month_avgstdGsPMall] = AM_month_avg(K67.monthGsPMall,K67.month);
[K67.eight_avgETdry,~,K67.eight_stdETdry] = AM_eight_avg(K67.eightETdry,K67.eight);
[K67.eight_avgETdryday,~,K67.eight_stdETdryday] = AM_eight_avg(K67.eightETdryday,K67.eight);
[K67.eight_avgETref,~,K67.eight_stdETref]       = AM_eight_avg(K67.eightETref,K67.eight);
[K67.eight_avgGsPMmmol,~,K67.eight_stdGsPMmmol] = AM_eight_avg(K67.eightGsPMmmol,K67.eight);
[K67.eight_avgGsPMallmmol,~,K67.eight_stdGsPMallmmol] = AM_eight_avg(K67.eightGsPMallmmol,K67.eight);
[K67.eight_avgPTalphaday,~,K67.eight_stdPTalphaday] = AM_eight_avg(K67.eightPTalphaday,K67.eight);
[K67.eight_avgVPDday,~,K67.eight_stdVPDday] = AM_eight_avg(K67.eightVPDday,K67.eight);
[K67.eight_avgLEdry,~,K67.eight_stdLEdry]   = AM_eight_avg(K67.eightLEdry,K67.eight);
[K67.eight_avgH,~,K67.eight_stdH]   = AM_eight_avg(K67.eightH,K67.eight);
[K67.eight_avgLE,~,K67.eight_stdLE] = AM_eight_avg(K67.eightLE,K67.eight);
[K67.eight_avgGiday,~,K67.eight_stdGiday] = AM_eight_avg(K67.eightGiday,K67.eight);
[K67.eight_avgEFday,~,K67.eight_stdEFday] = AM_eight_avg(K67.eightEFday,K67.eight);
[K67.eight_avgEFdaydry,~,K67.eight_stdEFdaydry] = AM_eight_avg(K67.eightEFdaydry,K67.eight);
[K67.eight_avgBowenday,~,K67.eight_stdBowenday] = AM_eight_avg(K67.eightBowenday,K67.eight);

[K67.month_avgBowen,K67.month_avgstdBowen] = AM_month_avg(K67.monthBowen,K67.month);
[K67.month_avgEF,K67.month_avgstdEF] = AM_month_avg(K67.monthEF,K67.month);
[K67.month_avgET,K67.month_avgstdET] = AM_month_avg(K67.monthET,K67.month);
[K67.month_avgGs,K67.month_avgstdGs] = AM_month_avg(K67.monthGsPMall,K67.month);
[K67.month_avgGsPMmmol,K67.month_avgstdGsPMmmol] = AM_month_avg(K67.monthGsPMmmol,K67.month);
[K67.month_avgGsPMallmmol,K67.month_avgstdGsPMallmmol] = AM_month_avg(K67.monthGsPMallmmol,K67.month);
[K67.month_avgGsPM,K67.month_avgstdGsPM]   = AM_month_avg(K67.monthGsPM,K67.month);
[K67.month_avgGsPMall,K67.month_avgstdGsPMall]   = AM_month_avg(K67.monthGsPMall,K67.month);
[K67.month_avgETdry,K67.month_avgstdETdry]   = AM_month_avg(K67.monthETdry,K67.month);
[K67.month_avgETdryday,K67.month_avgstdETdryday] = AM_month_avg(K67.monthETdryday,K67.month);
[K67.month_avgETref,K67.month_avgstdETref]   = AM_month_avg(K67.monthETref,K67.month);
[K67.month_avgETday,K67.month_avgstdETday]   = AM_month_avg(K67.monthETday,K67.month);
[K67.month_avgGiday,K67.month_avgstdGiday]   = AM_month_avg(K67.monthGiday,K67.month);
[K67.month_avgLEdry,K67.month_avgstdLEdry]   = AM_month_avg(K67.monthLEdry,K67.month);
[K67.month_avgPTalphaday,K67.month_avgstdPTalphaday] = AM_month_avg(K67.monthPTalphaday,K67.month);
[K67.month_avgVPDday,K67.month_avgstdVPDday] = AM_month_avg(K67.monthVPDday,K67.month);
[K67.month_avgEFday,K67.month_avgstdEFday]   = AM_month_avg(K67.monthEFday,K67.month);
[K67.month_avgEFdaydry,K67.month_avgstdEFdaydry] = AM_month_avg(K67.monthEFdaydry,K67.month);
[K67.month_avgBowenday,K67.month_avgstdBowenday] = AM_month_avg(K67.monthBowenday,K67.month);

%..........................................................................
%% Hydro year
%..........................................................................
K67.HY_Y = K67.Y;
K67.HY_Y(K67.M < month_start_hydro)   = K67.HY_Y(K67.M < month_start_hydro)-1;

K67.HY_dayY = K67.dayY;
K67.HY_dayY(K67.dayM < month_start_hydro)   = K67.HY_dayY(K67.dayM < month_start_hydro)-1;

K67.HY_weekY = K67.weekY;
K67.HY_weekY(K67.weekM < month_start_hydro)   = K67.HY_weekY(K67.weekM < month_start_hydro)-1;

K67.HY_eightY = K67.eightY;
K67.HY_eightY(K67.eightM < month_start_hydro) = K67.HY_eightY(K67.eightM < month_start_hydro)-1;

K67.HY_monthY = K67.monthY;
K67.HY_monthY(K67.monthM < month_start_hydro) = K67.HY_monthY(K67.monthM < month_start_hydro)-1;

% K67.HY_week = K67.Yweek - datenum(0,7,1);
% [K67.HY_weekY,K67.HY_Mweek,K67.HY_Dweek] = datevec(K67.HY_week);

all_yearENSO_tick = (datenum(2015,1,1):1:datenum(2016,12,31))';
[~,~,x1] = datevec(all_yearENSO_tick);          all_yearENSO_tick(x1~= 1) = [];
[~,MdateENSO,~]   = datevec(all_yearENSO_tick);
all_year_tick     = (datenum(2008,1,1):1:datenum(2016,12,31))';
[~,~,x1]          = datevec(all_year_tick);     all_year_tick(x1~= 1) = [];
[~,Mdate,~]       = datevec(all_year_tick);
allHA_year_tick   = (datenum(2002,1,1):1:datenum(2005,12,31))';
[~,~,x1]          = datevec(allHA_year_tick); 	allHA_year_tick(x1~= 1) = [];

TRMM.hy_dateloc   = TRMM.dateloc-datenum(0,7,1);
[TRMM.hy_y,TRMM.hy_m,TRMM.hy_d] = datevec(TRMM.hy_dateloc);

MEI.hy_dateloc    = MEI.dateloc-datenum(0,7,1);
[MEI.hy_y,MEI.hy_m,MEI.hy_d] = datevec(MEI.hy_dateloc);

% K67.weekCWD   = interp1(TRMM.dateloc,TRMM.monthCWD,K67.week);    % Better from 3 hr data
% K67.weekTRMM  = interp1(TRMM.dateloc,TRMM.precip,K67.week);
K67.monthCWD  = interp1(TRMM.dateloc,TRMM.monthCWD,K67.month);
K67.monthCWD100  = interp1(TRMM.dateloc,TRMM.monthCWD100,K67.month);
K67.monthTRMM = interp1(TRMM.dateloc,TRMM.precip,K67.month);

% [K67.week_avgCWD,~,K67.week_avgstdCWD]   = AM_week2_avg(K67.weekCWD,K67.week);
% [K67.week_avgTRMM,~,K67.week_avgstdTRMM] = AM_week2_avg(K67.weekTRMM,K67.week);
[K67.month_avgCWD,K67.month_avgstdCWD]   = AM_month_avg(K67.monthCWD,K67.month);
[K67.month_avgTRMM,K67.month_avgstdTRMM] = AM_month_avg(K67.monthTRMM,K67.month);

%..........................................................................
%% Meteorological values
%..........................................................................
x = K67.VPD;      x(K67.So<tresholdSo) = NaN;
[K67.weekVPDday,~,~,~,~,K67.week_stdVPDday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekVPDday] = AM_spike_detection(K67.weekVPDday,5);
% K67.weekVPDday   = AM_rm_outlier(K67.weekVPDday,3);
K67.weekVPDday     = AM_fill(K67.weekVPDday);
[K67.eightVPDday]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightVPDday] = AM_spike_detection(K67.eightVPDday,5);
% K67.eightVPDday   = AM_rm_outlier(K67.eightVPDday,3);
K67.eightVPDday     = AM_fill(K67.eightVPDday);
[K67.monthVPDday]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthVPDday] = AM_spike_detection(K67.monthVPDday,5);
% K67.monthVPDday   = AM_rm_outlier(K67.monthVPDday,3);
K67.monthVPDday     = AM_fill(K67.monthVPDday);
K67.monthVPD        = AM_month2cycle(K67.VPD,K67.dateloc,6,2,2); K67.monthVPD = AM_fill(K67.monthVPD);
[K67.month_avgVPD,K67.month_avgstdVPD] = AM_month_avg(K67.monthVPD,K67.month);

x = K67.tasp;      x(K67.So<tresholdSo) = NaN;
[K67.weekTaspday,~,~,~,~,K67.week_stdTaspday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekTaspday] = AM_spike_detection(K67.weekTaspday,5);
% K67.weekTaspday   = AM_rm_outlier(K67.weekTaspday,3);
K67.weekTaspday     = AM_fill(K67.weekTaspday);
[K67.week_avgTaspday,~,K67.week_avgstdTaspday] = AM_week2_avg(K67.weekTaspday,K67.week);
K67.week_stdTaspday(isnan(K67.week_avgTaspday)) = NaN;
[K67.eightTaspday]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightTaspday] = AM_spike_detection(K67.eightTaspday,5);
% K67.eightTaspday   = AM_rm_outlier(K67.eightTaspday,3);
K67.eightTaspday     = AM_fill(K67.eightTaspday);
[K67.monthTaspday]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthTaspday] = AM_spike_detection(K67.monthTaspday,5);
% K67.monthTaspday   = AM_rm_outlier(K67.monthTaspday,3);
K67.monthTaspday     = AM_fill(K67.monthTaspday);

x = K67.tair3;      x(K67.So<tresholdSo) = NaN;
[K67.weektair3day,~,~,~,~,K67.week_stdtair3day] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weektair3day] = AM_spike_detection(K67.weektair3day,5);
K67.weektair3day     = AM_fill(K67.weektair3day);
[K67.week_avgtair3day,~,K67.week_avgstdtair3day] = AM_week2_avg(K67.weektair3day,K67.week);
[K67.eighttair3day]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eighttair3day] = AM_spike_detection(K67.eighttair3day,5);
K67.eighttair3day     = AM_fill(K67.eighttair3day);
[K67.monthtair3day]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthtair3day] = AM_spike_detection(K67.monthtair3day,5);
K67.monthtair3day     = AM_fill(K67.monthtair3day);

x = K67.ustar;      x(K67.So<tresholdSo) = NaN;
[K67.weekustarday,~,~,~,~,K67.week_stdustarday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.weekustarday] = AM_spike_detection(K67.weekustarday,5);
% K67.weekustarday   = AM_rm_outlier(K67.weekustarday,3);
K67.weekustarday     = AM_fill(K67.weekustarday);
[K67.week_avgustarday,K67.week_avgstdustarday] = AM_week2_avg(K67.weekustarday,K67.week);
[K67.eightustarday]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightustarday] = AM_spike_detection(K67.eightustarday,5);
% K67.eightustarday   = AM_rm_outlier(K67.eightustarday,3);
K67.eightustarday     = AM_fill(K67.eightustarday);
[K67.monthustarday]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthustarday] = AM_spike_detection(K67.monthustarday,5);
% K67.monthustarday   = AM_rm_outlier(K67.monthustarday,3);
K67.monthustarday     = AM_fill(K67.monthustarday);

x = K67.ws;      x(K67.So<tresholdSo) = NaN;
[K67.weekWSday,~,~,~,~,K67.week_stdWSday] = AM_week2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightGiday] = AM_spike_detection(K67.eightGiday,5);
% K67.weekWSday    = AM_rm_outlier(K67.weekWSday,3);
K67.weekWSday  = AM_fill(K67.weekWSday);
[K67.week_avgWSday,K67.week_avgstdWSday] = AM_week2_avg(K67.weekWSday,K67.week);
[K67.eightWSday]   = AM_eight2cycle_rs(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.eightWSday] = AM_spike_detection(K67.eightWSday,5);
% K67.eightWSday   = AM_rm_outlier(K67.eightWSday,3);
K67.eightWSday     = AM_fill(K67.eightWSday);
[K67.monthWSday]   = AM_month2cycle(x,K67.dateloc,tresholdCycleDay,4,2);
[~,K67.monthWSday] = AM_spike_detection(K67.monthWSday,5);
% K67.monthWSday   = AM_rm_outlier(K67.monthWSday,3);
K67.monthWSday = AM_fill(K67.monthWSday);

[K67.weekH2O,~,~,~,~,K67.week_stdH2O] = AM_week2cycle_rs(K67.h2ommolm,K67.dateloc,22,2,2);

[K67.weekSo,~,~,~,~,K67.week_stdSo]   = AM_week2cycle_rs(K67.So,K67.dateloc,tresholdCycle,4,2);

% .........................................................................
% .........................................................................
[K67.weekStorageHeat,~,~,~,K67.weekStorageHeatMax,K67.week_stdStorageHeat] = AM_week2cycle_rs(K67.StorageHeat,K67.dateloc,tresholdCycle,4,2);
K67.weekStorageHeat       = AM_fill(K67.weekStorageHeat);
K67.weekStorageHeatMax    = AM_fill(K67.weekStorageHeatMax);
[K67.week_avgStorageHeat,~,K67.week_avgstdStorageHeat] = AM_week2_avg(K67.weekStorageHeat,K67.week);
[K67.week_avgStorageHeatMax,~,K67.week_avgstdStorageHeatMax] = AM_week2_avg(K67.weekStorageHeatMax,K67.week);

[K67.week_avgStorageHeat,~,K67.week_avgstdStorageHeat] = AM_week2_avg(K67.weekStorageHeat,K67.week);
[K67.weekStorageHeatCanopy,~,~,~,~,K67.week_stdStorageHeatCanopy]  = AM_week2cycle_rs(K67.StorageHeatCanopy,K67.dateloc,tresholdCycle,4,2);
K67.weekStorageHeatCanopy = AM_fill(K67.weekStorageHeatCanopy);
[K67.week_avgStorageHeatCanopy,~,K67.week_avgstdStorageHeatCanopy] = AM_week2_avg(K67.weekStorageHeatCanopy,K67.week);

[K67.weekStorageLE,~,~,~,K67.weekStorageLEMax,K67.week_stdStorageLE] = AM_week2cycle_rs(K67.StorageLE,K67.dateloc,tresholdCycle,4,2);
K67.weekStorageLE       = AM_fill(K67.weekStorageLE);
K67.weekStorageLEMax    = AM_fill(K67.weekStorageLEMax);
[K67.week_avgStorageLE,~,K67.week_avgstdStorageLE] = AM_week2_avg(K67.weekStorageLE,K67.week);
[K67.week_avgStorageLEMax,~,K67.week_avgstdStorageLEMax] = AM_week2_avg(K67.weekStorageLEMax,K67.week);

%..........................................................................
%% Select years
%..........................................................................
[K67.Yweek,K67.Mweek,K67.Dweek] = datevec(K67.week);
for ik = 1:length(K67.year)
    weekK67(ik).date = K67.year(ik);
    weekK67(ik).LE   = K67.weekLE(K67.Yweek == K67.year(ik));
    weekK67(ik).NetRad = K67.weekNetRad(K67.Yweek == K67.year(ik));
    weekK67(ik).DifRad = K67.weekDifRad(K67.Yweek == K67.year(ik));
    weekK67(ik).TotRad = K67.weekTotRad(K67.Yweek == K67.year(ik));
    weekK67(ik).H     = K67.weekH(K67.Yweek == K67.year(ik));
    weekK67(ik).Bowen = K67.weekH(K67.Yweek == K67.year(ik))./K67.weekLE(K67.Yweek == K67.year(ik));
    weekK67(ik).EF    = K67.weekLE(K67.Yweek == K67.year(ik))./(K67.weekH(K67.Yweek == K67.year(ik))+K67.weekLE(K67.Yweek == K67.year(ik)));
    weekK67(ik).EFday = K67.weekEFday(K67.Yweek == K67.year(ik));
    weekK67(ik).EFdaydry  = K67.weekEFdaydry(K67.Yweek == K67.year(ik));
    weekK67(ik).NetRadday = K67.weekNetRadday(K67.Yweek == K67.year(ik));
    weekK67(ik).LEday = K67.weekLEday(K67.Yweek == K67.year(ik));
    weekK67(ik).ETday = K67.weekETday(K67.Yweek == K67.year(ik));
    weekK67(ik).Hday  = K67.weekHday(K67.Yweek == K67.year(ik));
    weekK67(ik).Bowenday = K67.weekBowenday(K67.Yweek == K67.year(ik));
    weekK67(ik).Gs       = K67.weekGsPMall(K67.Yweek == K67.year(ik));
    weekK67(ik).GsPMmmol = K67.weekGsPMmmol(K67.Yweek == K67.year(ik));
    weekK67(ik).GsPM     = K67.weekGsPM(K67.Yweek == K67.year(ik));
    weekK67(ik).ETdry    = K67.weekETdry(K67.Yweek == K67.year(ik));
    weekK67(ik).ETdryday = K67.weekETdryday(K67.Yweek == K67.year(ik));
    weekK67(ik).ETref    = K67.weekETref(K67.Yweek == K67.year(ik));
    weekK67(ik).PTalphaday = K67.weekPTalphaday(K67.Yweek == K67.year(ik));
    weekK67(ik).VPDday= K67.weekVPDday(K67.Yweek == K67.year(ik));
    weekK67(ik).LEdry = K67.weekLEdry(K67.Yweek == K67.year(ik));
    weekK67(ik).Giday = K67.weekGiday(K67.Yweek == K67.year(ik));
    weekK67(ik).Tasp  = K67.weekTasp(K67.Yweek == K67.year(ik));
    weekK67(ik).H2O   = K67.weekH2O(K67.Yweek == K67.year(ik));
    weekK67(ik).ustarday = K67.weekustarday(K67.Yweek == K67.year(ik));
    weekK67(ik).VPDday   = K67.weekVPDday(K67.Yweek == K67.year(ik));
    weekK67(ik).WSday    = K67.weekWSday(K67.Yweek == K67.year(ik));
    weekK67(ik).Taspday  = K67.weekTaspday(K67.Yweek == K67.year(ik));
    weekK67(ik).tair3day = K67.weektair3day(K67.Yweek == K67.year(ik));
    weekK67(ik).CWD      = K67.weekCWD(K67.Yweek == K67.year(ik));
    weekK67(ik).Rain     = K67.weekRain(K67.Yweek == K67.year(ik));
    
    weekK67hy(ik).date   = K67.year(ik);
    weekK67hy(ik).Y      = K67.HY_weekY(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).LE     = K67.weekLE(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).NetRad = K67.weekNetRad(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).DifRad = K67.weekDifRad(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).TotRad = K67.weekTotRad(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).H      = K67.weekH(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Bowen  = K67.weekH(K67.HY_weekY == K67.year(ik))./K67.weekLE(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).EF     = K67.weekLE(K67.HY_weekY == K67.year(ik))./(K67.weekH(K67.HY_weekY == K67.year(ik))+K67.weekLE(K67.HY_weekY == K67.year(ik)));
    weekK67hy(ik).EFday  = K67.weekEFday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).EFdaydry  = K67.weekEFdaydry(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).NetRadday = K67.weekNetRadday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Hday   = K67.weekHday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).LEday  = K67.weekLEday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).ETday  = K67.weekETday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Bowenday = K67.weekBowenday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).GsPMall  = K67.weekGsPMall(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).GsPMmmol = K67.weekGsPMmmol(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).GsPM    = K67.weekGsPM(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).ETdry   = K67.weekETdry(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).ETdryday= K67.weekETdryday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).ETref   = K67.weekETref(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Giday   = K67.weekGiday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).PTalphaday = K67.weekPTalphaday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).VPDday  = K67.weekVPDday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).LEdry   = K67.weekLEdry(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Tasp    = K67.weekTasp(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).H2O     = K67.weekH2O(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).ustarday = K67.weekustarday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).VPDday   = K67.weekVPDday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).WSday    = K67.weekWSday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Taspday  = K67.weekTaspday(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).tair3day = K67.weektair3day(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).CWD      = K67.weekCWD(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Rain     = K67.weekRain(K67.HY_weekY == K67.year(ik));
    
    monthK67(ik).date = K67.HY_weekY(ik);
    monthK67(ik).CWD  = TRMM.monthCWD(TRMM.y == K67.year(ik));
    monthK67(ik).TRMM = TRMM.precip(TRMM.y == K67.year(ik));
    monthK67(ik).TRMManomaly = TRMM.month_anomalyPrecip(TRMM.y == K67.year(ik));
    monthK67(ik).MEI  = MEI.value_month(MEI.y == K67.year(ik));
    monthK67(ik).CERESnr  = TRMM.monthCERESnr(TRMM.y == K67.year(ik));
    monthK67(ik).CERESswd = TRMM.monthCERESswd(TRMM.y == K67.year(ik));
    monthK67(ik).Rain = K67.monthRain(K67.monthY == K67.year(ik));
    
    monthK67hy(ik).date = K67.HY_weekY(ik);
    monthK67hy(ik).CWD  = TRMM.monthCWD(TRMM.hy_y == K67.year(ik));
    monthK67hy(ik).TRMM = TRMM.precip(TRMM.hy_y == K67.year(ik));
    monthK67hy(ik).TRMManomaly = TRMM.month_anomalyPrecip(TRMM.hy_y == K67.year(ik));
    monthK67hy(ik).MEI  = MEI.value_month(MEI.hy_y == K67.year(ik));
    monthK67hy(ik).CERESnr  = TRMM.monthCERESnr(TRMM.hy_y == K67.year(ik));
    monthK67hy(ik).CERESswd = TRMM.monthCERESswd(TRMM.hy_y == K67.year(ik));
    monthK67hy(ik).Rain = K67.monthRain(K67.HY_monthY == K67.year(ik));
    
    eightK67(ik).date = K67.year(ik);
    eightK67(ik).LE   = K67.eightLE(K67.eightY == K67.year(ik));
    eightK67(ik).NetRad = K67.eightNetRad(K67.eightY == K67.year(ik));
    eightK67(ik).DifRad = K67.eightDifRad(K67.eightY == K67.year(ik));
    eightK67(ik).TotRad = K67.eightTotRad(K67.eightY == K67.year(ik));
    eightK67(ik).H     = K67.eightH(K67.eightY == K67.year(ik));
    eightK67(ik).Bowen = K67.eightH(K67.eightY == K67.year(ik))./K67.eightLE(K67.eightY == K67.year(ik));
    eightK67(ik).EF    = K67.eightLE(K67.eightY == K67.year(ik))./(K67.eightH(K67.eightY == K67.year(ik))+K67.eightLE(K67.eightY == K67.year(ik)));
    eightK67(ik).EFday = K67.eightEFday(K67.eightY == K67.year(ik));
    eightK67(ik).EFdaydry = K67.eightEFdaydry(K67.eightY == K67.year(ik));
    eightK67(ik).LEday = K67.eightLEday(K67.eightY == K67.year(ik));
    eightK67(ik).ETday = K67.eightETday(K67.eightY == K67.year(ik));
    eightK67(ik).Hday  = K67.eightHday(K67.eightY == K67.year(ik));
    eightK67(ik).NetRadday = K67.eightNetRadday(K67.eightY == K67.year(ik));
    eightK67(ik).Bowenday  = K67.eightBowenday(K67.eightY == K67.year(ik));
    eightK67(ik).Gs      = K67.eightGsPMall(K67.eightY == K67.year(ik));
    eightK67(ik).GsPM   = K67.eightGsPM(K67.eightY == K67.year(ik));
    eightK67(ik).ETdry   = K67.eightETdry(K67.eightY == K67.year(ik));
    eightK67(ik).ETdryday= K67.eightETdryday(K67.eightY == K67.year(ik));
    eightK67(ik).ETref   = K67.eightETref(K67.eightY == K67.year(ik));
    eightK67(ik).LEdry = K67.eightLEdry(K67.eightY == K67.year(ik));
    eightK67(ik).Giday   = K67.eightGiday(K67.eightY == K67.year(ik));
    eightK67(ik).PTalphaday = K67.eightPTalphaday(K67.eightY == K67.year(ik));
    eightK67(ik).VPDday  = K67.eightVPDday(K67.eightY == K67.year(ik));
    eightK67(ik).Tasp = K67.eightTasp(K67.eightY == K67.year(ik));
    eightK67(ik).H2O  = K67.eightH2O(K67.eightY == K67.year(ik));
    eightK67(ik).Taspday  = K67.eightTaspday(K67.eightY == K67.year(ik));
    eightK67(ik).tair3day  = K67.eighttair3day(K67.eightY == K67.year(ik));
    eightK67(ik).ustarday = K67.eightustarday(K67.eightY == K67.year(ik));
    eightK67(ik).VPDday   = K67.eightVPDday(K67.eightY == K67.year(ik));
    eightK67(ik).WSday    = K67.eightWSday(K67.eightY == K67.year(ik));
    eightK67(ik).Rain    = K67.eightRain(K67.eightY == K67.year(ik));
    
    eightK67hy(ik).date   = K67.HY_eightY(ik);
    eightK67hy(ik).LE     = K67.eightLE(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).NetRad = K67.eightNetRad(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).DifRad = K67.eightDifRad(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).TotRad = K67.eightTotRad(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).H      = K67.eightH(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Bowen  = K67.eightH(K67.HY_eightY == K67.year(ik))./K67.eightLE(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).EF     = K67.eightLE(K67.HY_eightY == K67.year(ik))./(K67.eightH(K67.HY_eightY == K67.year(ik))+K67.eightLE(K67.HY_eightY == K67.year(ik)));
    eightK67hy(ik).EFday  = K67.eightEFday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).EFdaydry  = K67.eightEFdaydry(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).NetRadday = K67.eightNetRadday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Hday   = K67.eightHday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).LEday  = K67.eightLEday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).ETday  = K67.eightETday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Bowenday = K67.eightBowenday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Gs       = K67.eightGsPMall(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).GsPM    = K67.eightGsPM(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).ETdry    = K67.eightETdry(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).ETdryday = K67.eightETdryday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).ETref    = K67.eightETref(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).PTalphaday = K67.eightPTalphaday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).VPDday = K67.eightVPDday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).LEdry  = K67.eightLEdry(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Giday  = K67.eightGiday(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Tasp   = K67.eightTasp(K67.HY_eightY == K67.year(ik));
    eightK67hy(ik).Rain   = K67.eightRain(K67.HY_eightY == K67.year(ik));
end

%..........................................................................
%% Cycle
%..........................................................................
[K67.cycle_avgLE]     = AM_day_cycle(K67.LE,K67.dateloc,6,18);
[K67.cycle_avgH]      = AM_day_cycle(K67.H,K67.dateloc,6,18);
[K67.cycle_avgNetRad] = AM_day_cycle(K67.netrad,K67.dateloc,6,18);
[K67.cycle_avgNEE]    = AM_day_cycle(K67.neefill,K67.dateloc,6,18);
[K67.cycleLE,~,K67.cycle] = AM_day_cycle_year(K67.LE,K67.dateloc);
[K67.cycleH]          = AM_day_cycle_year(K67.H,K67.dateloc);
[K67.cycleNetRad]     = AM_day_cycle_year(K67.netrad,K67.dateloc);
[K67.cycleNEE]        = AM_day_cycle_year(K67.neefill,K67.dateloc);

%..........................................................................
%% Repeating annual cycle
%..........................................................................
site_LE_mean_week = [];     site_H_mean_week = [];
site_LE_std_week  = [];     site_H_std_week  = [];
for ik = 1:length(unique(K67.weekY))
    site_LE_mean_week = [site_LE_mean_week;K67.week_avgLE];
    site_LE_std_week  = [site_LE_std_week;K67.week_avgstdLE];
    site_H_mean_week  = [site_H_mean_week;K67.week_avgH];
    site_H_std_week   = [site_H_std_week;K67.week_avgstdH];
end

%..........................................................................
%% Time vectors seasonality
%..........................................................................
week_avg  = (1:16:365)';
month_avg = (1:1:365)';       [~,~,x] = datevec(month_avg);     month_avg(x~= 1) = [];
year_trmm = (1998:1:2016)';

month_avg_2year   = [month_avg;month_avg+366];
all_yearENSO_week = [(datenum(2015,1,week_avg));(datenum(2016,1,week_avg))];

HYmonth_avg  = [month_avg(month_start_hydro:12);month_avg(1:month_start_hydro-1)];

save ('K67data.mat','K67*','monthK67*','weekK67*','TRMM*')

%........................................................................
%% Save daytime Gs and ET
%..........................................................................
[K67.dayY,K67.dayM,K67.dayD]= datevec(K67.day);
ind = [];       % ind = find(K67.dayY==2013);
K67.dayGsPMall(ind)    = NaN;    K67.dayGsPMallmax(ind) = NaN;        K67.dayGsPMallPM(ind) = NaN;
K67.dayETday(ind) = NaN;    K67.dayETdayMAX(ind) = NaN;     K67.dayETdry(ind) = NaN; K67.dayVPDday(ind);
K67.dayETdryday(ind) = NaN;
ind = 1:length(K67.day);
tableK67day = table(K67.dayY(ind),K67.dayM(ind),K67.dayD(ind),K67.dayGsPMall(ind),K67.dayGsPMallmax(ind),K67.dayGsPM(ind),K67.dayETday(ind),K67.dayETdayMAX(ind),K67.dayETdry(ind),K67.dayETdryday(ind),K67.dayVPDday(ind));
tableK67day.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'Gs_mm_s1'  'GsMAX_mm_s1'   'Gs_fnT_mm_s1' 'ETdaytime_mm_d1'  'ETdaytime_max_mm_d1' 'ETfnT_mm_d1' 'ETdayfnT_mm_d1' 'VPDday_kPa'};
writetable (tableK67day,'out/dayK67.ET.Gs.dat');

ind = [];% ind = find(K67.HY_monthY==2013);
K67.monthGsPMall(ind) = NaN;    K67.monthGsPMallmin(ind) = NaN;     K67.monthGsPM(ind) = NaN;
K67.monthETday(ind) = NaN;      K67.monthETdayMAX(ind) = NaN;       K67.monthETdry(ind) = NaN;
K67.monthVPDday(ind);
K67.monthETdryday(ind) = NaN;
ind = 1:length(K67.month);    %
tableK67month = table(K67.monthY(ind),K67.monthM(ind),K67.monthD(ind),K67.monthGsPMall(ind),K67.monthGsPMallmax(ind),K67.monthGsPM(ind),K67.monthETday(ind),...
    K67.monthETdayMAX(ind),K67.monthETdry(ind),K67.monthETdryday(ind),K67.monthVPDday(ind),K67.monthTRMM(ind),K67.monthCWD(ind));
tableK67month.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'Gs_mm_s1'  'GsMAX_mm_s1'   'Gs_fnT_mm_s1' 'ETdaytime_mm_d1'  'ETdaytime_max_mm_d1' 'ETfnT_mm_d1' 'ETdayfnT_mm_d1' 'VPDday_kPa' 'PrecipTRMM_mm_m1' 'CWD_mm_m1'};
writetable (tableK67month,'out/monthK67.ET.Gs.dat');

ind = 1:length(K67.week);    %
tableK67week = table(K67.weekY(ind),K67.weekM(ind),K67.weekD(ind),K67.weekTRMM(ind));
tableK67week.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'PrecipTRMM_mm_m1'};
writetable (tableK67week,'out/weekK67.TRMM.dat');

%..........................................................................
%% Save 
%..........................................................................
tableK67 = table(K67.Y,K67.M,K67.D,K67.HR,K67.H,K67.LE,K67.fco2,K67.ws,...
    K67.tson,K67.h2omix,K67.co2,K67.wsdir,K67.pambpa,K67.fmom,K67.tasp,K67.RH,K67.totsun,...
    K67.difsun,K67.Y.*NaN,K67.Y.*NaN,K67.tair3,K67.tair2,K67.tair3,K67.tair4,K67.rain,...
    K67.ustar,K67.par,K67.netrad);
tableK67.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'HR' 'H_Wm2','LE_Wm2','Fc_umolm2s1','ws_ms1','tson_C','h2omix_gkg1','CO2_pmm','wdir_azimut','pamb_pa','Tau_Nm2','tair_C','RH','totsun_Wm2','difsun_Wm2','rgl_Wm2','rglout_Wm2','tair1_C','tair2_C','tair3_C','tair4_C','rain_mm','ustar_ms1','par_umolm2s1','netrad_Wm2'};
writetable (tableK67,'out/K67.data.dat');

tableK67 = [];
tableK67 = table(K67.Y,K67.M,K67.D,K67.HR,K67.H,K67.LE,K67.GsPM,K67.GsPMmmol,K67.netrad,K67.VPD,K67.TRMM1hPrecip,K67.TRMM1hCWD,K67.tasp,K67.par);
tableK67.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'HR' 'H_Wm2','LE_Wm2','Gs_mms1','Gs_mmolm2s1','netrad_Wm2','VPD_kPa','TRMM_mm','CWD_mm','Tair_C','PAR_umolm2s1'};
writetable (tableK67,'out/K67.data.hour.dat');

%% ........................................................................
K67.dateloc30minGMT = (K67.dateloc(1):datenum(0,0,0,0,30,0):K67.dateloc(end))';
K67.dateloc30minGMT = K67.dateloc30minGMT + datenum(0,0,0,utm,0,0);
[K67.Y30minGMT,K67.M30minGMT,K67.D30minGMT,K67.HR30minGMT,K67.MI30minGMT,K67.SS30minGMT] = datevec(K67.dateloc30minGMT);
K67.MI30minGMT((K67.MI30minGMT>27)&(K67.MI30minGMT<32)) = 30;
K67.TRMM1hCWD30minGMT = interp1(K67.dateloc,K67.TRMM1hCWD,K67.dateloc30minGMT);

K67.TRMM1hPrecip30minGMT = interp1(K67.dateloc,K67.TRMM1hPrecip,K67.dateloc30minGMT);
K67.TRMM1hPrecip30minGMT = K67.TRMM1hPrecip30minGMT./2;
% K67.TRMM1hPrecip30minGMT = interp(K67.TRMM1hPrecip,2);

K67.ETfill30minGMT = interp1(K67.dateloc,K67.ETfill,K67.dateloc30minGMT);
K67.ET30minGMT     = interp1(K67.dateloc,K67.ET,K67.dateloc30minGMT);
K67.GsPM30minGMT  = interp1(K67.dateloc,K67.GsPM,K67.dateloc30minGMT);
K67.GsPMmmol30minGMT = interp1(K67.dateloc,K67.GsPMmmol,K67.dateloc30minGMT);
K67.netrad30minGMT = interp1(K67.dateloc,K67.netrad,K67.dateloc30minGMT);
K67.tasp30minGMT = interp1(K67.dateloc,K67.tasp,K67.dateloc30minGMT);
K67.VPD30minGMT = interp1(K67.dateloc,K67.VPD,K67.dateloc30minGMT);
K67.par30minGMT = interp1(K67.dateloc,K67.par,K67.dateloc30minGMT);

TableGsCWD = table(K67.Y30minGMT, K67.M30minGMT, K67.D30minGMT, K67.HR30minGMT, K67.MI30minGMT,....
    K67.ET30minGMT,K67.GsPM30minGMT,K67.GsPMmmol30minGMT,K67.netrad30minGMT,K67.VPD30minGMT,...
    K67.TRMM1hPrecip30minGMT,K67.TRMM1hCWD30minGMT,K67.tasp30minGMT,K67.par30minGMT,...
    'VariableNames',{'YYYY' 'MM' 'DD' 'HR' 'MI' 'ET_mm','Gs_mms1','Gs_mmolm2s1','netrad_Wm2','VPD_kPa','TRMM_mm','CWD_mm','Tair_C','PAR_umolm2s1'}) ;
writetable(TableGsCWD,'K6730min.Gs.CWD.dat');

%% ........................................................................
tableK67 =[];
tableK67 = table(K67.weekY,K67.weekM,K67.weekD,K67.weekD.*0,K67.weekH,K67.weekLE,K67.weekGsPM,K67.weekGsPMmmol,K67.weekNetRad,K67.weekTasp,K67.weekPAR);
tableK67.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'HR' 'H_Wm2','LE_Wm2','Gs_mms1','Gs_mmolm2s1','netrad_Wm2','Tair_C','PAR_umolm2s1'};
writetable (tableK67,'out/K67.data.16day.dat');

tableK67 =[];
tableK67 = table(K67.HY_monthY,K67.monthM,K67.monthD,K67.monthD.*0,K67.monthH,K67.monthLE,K67.monthGsPM,K67.monthGsPMmmol,K67.monthNetRad,K67.monthVPD,K67.monthTRMM,K67.monthCWD,K67.monthTasp,K67.monthPAR);
tableK67.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'HR' 'H_Wm2','LE_Wm2','Gs_mms1','Gs_mmolm2s1','netrad_Wm2','VPD_kPa','TRMM_mm','CWD_mm','Tair_C','PAR_umolm2s1'};
writetable (tableK67,'out/K67.data.month.dat');

%%  .......................................................................
% K67.monthGsPMnofill = K67.monthGsPM;
figure('color','white');    subplot(2,1,1);
plot(K67.month,K67.monthET,'LineWidth',2);  hold on;
plot(K67.month,K67.monthETdry,'LineWidth',2,'color',[.3,.8,.2])
legend ('ET','ET_{dry}');     datetick('x');              ylabel('ET K67 (mm d^-^1)');
subplot(2,1,2);  plot(K67.month,K67.monthGsPMall,'LineWidth',2);  hold on;
plot(K67.month,K67.monthGsPM,'color',[.3,.8,.2],'LineWidth',2)
datetick('x');              ylabel('Gs (mm s^-^1)');
plot(K67.month,K67.monthGsPMmax,'color',[.6,.2,.4],'LineWidth',2)
% plot(K67.month,K67.monthGsPMnofill,'LineWidth',2,'color',[.8,.3,.2]);
legend ('Gs','Gs_{dry}','Gs_{max}','Gs_{dry TRMM rain filled K67}');

%% ........................................................................
[K67.weeknee,~,~,~,~,K67.week_stdnee] = AM_week2cycle_rs(K67.neefill,K67.dateloc,tresholdCycle,4,2);
% [~,K67.weeknee] = AM_spike_detection(K67.weeknee,5);
K67.weeknee = AM_rm_outlier(K67.weeknee,2.5);
K67.weeknee = AM_fill(K67.weeknee);

[K67.weekneelight,~,~,~,~,K67.week_stdneelight] = AM_week2cycle_rs(K67.neefilllight,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekneelight] = AM_spike_detection(K67.weekneelight,5);
% K67.weeknee   = AM_rm_outlier(K67.weeknee,3);
K67.weekneelight = AM_fill(K67.weekneelight);

ix = K67.weeknee.*gC_m2_d1;
ix(isnan(K67.weeknee)) = 0;
iz = [diff(K67.week);13];
iy = cumsum(ix.*iz);        iy = iy.*10000./1000;
figure('color','white');    plot(K67.week,iy);      grid on;
datetick('x');              ylabel({'Cumulative NEE';'(kgC ha^-^1)'});

ind = find((K67.week>=datenum(2002,1,1))&(K67.week<=datenum(2002,9,1)));
figure('color','white');    plot(K67.week(ind),iy(ind));      grid on;
datetick('x','m');              ylabel({'Cumulative NEE';'(kgC ha^-^1)'});

%% ........................................................................
[K67.daynee,~,ix] = AM_day(K67.neefill,K67.dateloc,0);          K67.daynee(ix<22) = NaN;
K67.daynee        = AM_rm_outlier(K67.daynee,3);                K67.daynee = AM_fill_plus(K67.daynee);
[K67.dayneelight,~,ix] = AM_day(K67.neefilllight,K67.dateloc,0);K67.dayneelight(ix<22) = NaN;
K67.dayneelight        = AM_rm_outlier(K67.dayneelight,3);      K67.dayneelight = AM_fill_plus(K67.dayneelight);

ix = K67.daynee.*gC_m2_d1;
ix(isnan(K67.daynee)) = 0;
iz = ones(length(K67.daynee),1);             %[diff(K67.day);13];
iy = cumsum(ix.*iz);       iy = iy.*10000./1000;
figure('color','white');    plot(K67.day,iy,'LineWidth',2);      grid on;
set(gca,'XTick',K67.month(1:12:end),'XTickLabel',datestr(K67.month(1:12:end),'yyyy'),...
    'XLim',[K67.month(1) datenum(2020,8,1)],'FontSize',12);
ylabel({'Cumulative NEE';'(kgC ha^-^1)'},'FontSize',12);

% ........................................................................
ind = find((K67.day>=datenum(2002,1,1))&(K67.day<datenum(2003,8,1)));
figure('color','white');    plot(K67.day(ind),iy(ind));      grid on;
datetick('x','m');              ylabel({'Cumulative NEE';'(kgC ha^-^1)'});

% ........................................................................
ix = K67.daynee.*gC_m2_d1;
ix(isnan(K67.daynee)) = 0;
iz = ones(length(K67.daynee),1);             %[diff(K67.day);13];
iy = cumsum(ix.*iz);       iy = iy.*10000./1000;
figure('color','white');    plot(K67.day,iy);      grid on;
datetick('x');              ylabel({'Cumulative NEE';'(kgC ha^-^1)'});

ind = find((K67.day>=datenum(2002,1,1))&(K67.day<datenum(2003,8,1)));
figure('color','white');    plot(K67.day(ind),iy(ind));      grid on;
datetick('x','m');              ylabel({'Cumulative NEE';'(kgC ha^-^1)'});

% .........................................................................
%% Cycles
% .........................................................................
[K67.cycleRn,K67.cycle_stdRn,~,K67.acycle_stdRn,K67.acycleRn] = AM_cycle_composite(K67.netrad,K67.dateloc);
[K67.cycle_monthRn,~,~,K67.cycle_month_stdRn] = AM_day_cycle_year(K67.netrad,K67.dateloc);
[K67.cycleSo,K67.cycle_stdSo,~,K67.acycle_stdSo,K67.acycleSo] = AM_cycle_composite(K67.So,K67.dateloc);
[K67.cycle_monthSo,~,~,K67.cycle_month_stdSo] = AM_day_cycle_year(K67.So,K67.dateloc);

% .........................................................................
%% Trend
% .........................................................................
g = length(unique(TRMM.y));
TRMM.month_trendPrecip   = NaN((g*12),1);
for ik = 1:g*12
    iq = TRMM.m(ik);
    A  = TRMM.month_avgPrecip(iq);
    C  = TRMM.precip(ik);
    D  = (C-A);
    TRMM.month_trendPrecip(ik) = D;
end
p = polyfit(TRMM.dateloc(~isnan(TRMM.month_trendPrecip)),TRMM.month_trendPrecip(~isnan(TRMM.month_trendPrecip)),1);
y = polyval(p,TRMM.dateloc);

figure('color','white');
plot(TRMM.dateloc,TRMM.month_trendPrecip);                    hold on;
plot(TRMM.dateloc,y,'.','color',[0.9 .5 .5]);
K67.monthTRMMtrend = interp1(TRMM.dateloc,TRMM.month_trendPrecip,K67.month);

% % figure('color','white');
% % plot(monthCERESdate,month_trendCERESswd);                    hold on;
% % wts = [1/24;repmat(1/12,11,1);1/24];
% % yS = conv(y,monthCERESswd,'valid');
% % h = plot(7:monthCERESdate-6,yS,'r','LineWidth',2);
% % legend(h,'13-Term Moving Average')
% % hold off

% % g = length(unique(TRMM.y));
% % TRMM.month_trendPrecip   = NaN((g*12),1);
% % for ik = 1:g*12
% %     iq = TRMM.m(ik);
% %     A  = TRMM.month_avgPrecip(iq);
% %     C  = TRMM.precip(ik);
% %     D  = (C-A);
% %     TRMM.month_trendPrecip(ik) = D;
% % end
% % p  = polyfit(TRMM.dateloc(~isnan(TRMM.month_trendPrecip)),TRMM.month_trendPrecip(~isnan(TRMM.month_trendPrecip)),1);
% % y   = polyval(p,TRMM.dateloc);
% % figure('color','white');
% % plot(TRMM.dateloc,TRMM.month_trendPrecip);                    hold on;
% % plot(TRMM.dateloc,y,'.','color',[0.9 .5 .5]);
% % K67.monthTRMMtrend = interp1(TRMM.dateloc,TRMM.month_trendPrecip,K67.month);


%% ........................................................................
K67.weekBelterraCI = AM_week2day_rs(K67b.CI,K67.dateloc);
[K67.week_avgBelterraCI,~,K67.week_avgstdBelterraCI] = AM_week2_avg(K67.weekBelterraCI,K67.week);
[K67.weekBelterraCIrecurrent] = AM_week2dayAVG_recurrent(K67.week_avgBelterraCI,K67.week);

K67.weekBelterraCInoon = AM_week2day_rs(K67b.CInoon,K67.dateloc);
[K67.week_avgBelterraCInoon,~,K67.week_avgstdBelterraCInoon] = AM_week2_avg(K67.weekBelterraCInoon,K67.week);
% [K67.weekBelterraCInoonrecurrent] = AM_week2dayAVG_recurrent(K67.week_avgBelterraCI,K67.week);

K67.monthBelterraCInoon = AM_month(K67b.CInoon,K67.dateloc);
[K67.month_avgBelterraCInoon,~,K67.month_avgstdBelterraCInoon] = AM_month_avg(K67.monthBelterraCInoon,K67.month);


K67.weekCI = AM_week2day_rs(K67.CI,K67.dateloc);
[K67.week_avgCI,~,K67.week_avgstdCI] = AM_week2_avg(K67.weekCI,K67.week);
[K67.weekCIrecurrent] = AM_week2dayAVG_recurrent(K67.week_avgCI,K67.week);

