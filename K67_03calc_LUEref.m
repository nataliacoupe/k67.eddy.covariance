% Calculates LUEref as in Wu et al  science manuscript
% Natalia Restrepo-Coupe
% Toronto
% .........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;

month_name = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];

color_mtx = [0.1020 0.5882 0.2549; 215/255  25/255 28/255; 43/255 131/255 186/255; ...
    253/255 174/255 97/255; 0.2549 0.1020 0.5882; 81/255 161/255 74/255;0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

umol2gC  = 60.*60.*24.*12.*(10.^-6);
cycle_hour = (0:23)';

K67.gepDry    = -K67.gee;           K67.gepDry(indDry==1)            = NaN;
K67.gepDryDay = K67.gepDry;         K67.gepDryDay(K67.So<tresholdSo) = NaN;
K67.gepday    = -K67.gee;           K67.gepday(K67.So<tresholdSo)    = NaN;

K67.rDry   = K67.r;                 K67.rDry(indDry==1)   = NaN;
K67.neeDry = K67.neefill;           K67.neeDry(indDry==1) = NaN;

%% ..........................................................................
% K67.NEE       = K67.fco2 + K67.storagewt;
% K67.NEE(K67.ustar<ust_tresh) = NaN;                 K67.NEE = AM_rm_outlier(K67.NEE,3);
% [~,K67.NEE]   = AM_spike_detection(K67.NEE,5);
%
% [~,K67.storage] = AM_spike_detection(K67.storagewt,5);
% K67.storage = AM_rm_outlier(K67.storagewt,3);
% K67.storage(95470:99120) = NaN;                 K67.storage(95470:99120) = NaN;
% K67.storage(130480:131700) = NaN;               K67.fco2(130480:131700) = NaN;
% K67.nee   = K67.fco2 + K67.storage;
% [K67.GEPmodel_ust_min,K67.GEPmodel_ust_max,K67.GEP_ec,K67.Remodel_ust_min,K67.Remodel_ust_max,...
%     K67.Re_ec,K67.NEE_ec_ust_min,K67.NEE_ec_ust_max,K67.NEE_ec] = ...
%     AM_gep_ust(K67.nee,K67.radiation_plus,K67.ustar,ust_min,ust_tresh,ust_max,...
%     lat,lon,K67.dateloc+datenum(0,0,0,1,0,0),0,-3,-3,30);    %utm

%..........................................................................
%% LUEref
%..........................................................................
figure('color','white');
subplot(2,2,1);           hold on;
iy = interp1(K67am.dateloc,K67am.VPD./1000,K67.dateloc);
[pVPDam,r2,~,~,sm] = AM_compute_regression(K67.VPD_fill,iy,[0.4,0.4,0.4],1);
xlabel('VPD'); ylabel('VPD_{ameriflux}');

subplot(2,2,3:4);           hold on;
plot(K67.dateloc,K67.VPD);  plot(K67am.dateloc,K67am.VPD./1000)
plot(K67.dateloc,(polyval(p,K67.VPD_fill)));  legend('UofA','Ameriflux');

[K67.weekLUEref,weekLUEref_count,weekLUEref_std] = ...
    AM_week2day_LUEref(K67.gepday,K67.radiation_plus,polyval(pVPDam,K67.VPD_fill),K67.tc_fill,K67.dateloc,lat,lon,-1);
K67.weekLUEref(K67.weekLUEref>0.026) = NaN;
K67.weekLUEref = AM_rm_outlier(K67.weekLUEref,3);
ix = smooth(K67.weekLUEref,3);    ix(isnan(K67.weekLUEref)) = NaN;    K67.weekLUEref = ix;
K67.weekLUEref  = AM_fill_plus(K67.weekLUEref);

[K67.monthLUEref,monthLUEref_count,monthLUEref_std] = ...
    AM_month_LUEref(K67.gepday,K67.radiation_plus,polyval(pVPDam,K67.VPD_fill),K67.tc_fill,K67.dateloc,lat,lon,-1);
K67.monthLUEref(K67.monthLUEref>0.027) = NaN;
K67.monthLUEref = AM_rm_outlier(K67.monthLUEref,3);
ix = smooth(K67.monthLUEref,3);   ix(isnan(K67.monthLUEref)) = NaN;   K67.monthLUEref = ix;
K67.monthLUEref = AM_fill_plus(K67.monthLUEref);

[K67.weekLUErefRELAX,weekLUErefRELAX_count,weekLUErefRELAX_std] = ...
    AM_week2day_LUErefRELAX(K67.gepday,K67.radiation_plus,polyval(pVPDam,K67.VPD_fill),K67.tc_fill,K67.dateloc,lat,lon,-1);
K67.weekLUErefRELAX(K67.weekLUErefRELAX>0.03) = NaN;  
K67.weekLUErefRELAX = AM_rm_outlier(K67.weekLUErefRELAX,3);
K67.weekLUErefRELAX  = AM_fill_plus(K67.weekLUErefRELAX);
ix = smooth(K67.weekLUErefRELAX,3);    ix(isnan(K67.weekLUErefRELAX)) = NaN;    K67.weekLUErefRELAX = ix;

[K67.monthLUErefRELAX,monthLUErefRELAX_count,monthLUErefRELAX_std] = ...
    AM_month_LUErefRELAX(K67.gepday,K67.radiation_plus,polyval(pVPDam,K67.VPD_fill),K67.tc_fill,K67.dateloc,lat,lon,-1);
K67.monthLUErefRELAX(K67.monthLUErefRELAX>0.028) = NaN;
K67.monthLUErefRELAX = AM_rm_outlier(K67.monthLUErefRELAX,3);
ix = smooth(K67.monthLUErefRELAX,3);    ix(isnan(K67.monthLUErefRELAX)) = NaN;    K67.monthLUErefRELAX = ix;
K67.monthLUErefRELAX  = AM_fill_plus(K67.monthLUErefRELAX);

% K67.weekLUEref(105) = NaN;
% K67.weekLUEref  = AM_rm_outlier(K67.weekLUEref,3);
% K67.monthLUEref = AM_rm_outlier(K67.monthLUEref,3);
% K67.monthLUEref(K67.monthLUEref>0.028) = NaN;
% K67.monthLUEref((K67.monthY<=2010)&(K67.monthLUEref<0.011)) = NaN;

%% ........................................................................
[K67.week_avgLUEref,~,K67.week_avgstdLUEref]   = AM_week2_avg(K67.weekLUEref,K67.week);
[K67.month_avgLUEref,~,K67.month_avgstdLUEref] = AM_month_avg(K67.monthLUEref,K67.month);
[K67.week_recurrentLUEref]  = AM_week2dayAVG_recurrent(K67.week_avgLUEref,K67.week);
[K67.month_recurrentLUEref] = AM_monthAVG_recurrent(K67.month_avgLUEref,K67.month);

% ........................................................................
%% ........................................................................
K67.weekLUErefRELAX = AM_rm_outlier(K67.weekLUErefRELAX,3);
ix = smooth(K67.weekLUErefRELAX,3);    ix(isnan(K67.weekLUErefRELAX)) = NaN;    K67.weekLUErefRELAX = ix;
K67.weekLUErefRELAX  = AM_fill_plus(K67.weekLUErefRELAX);

K67.monthLUErefRELAX = AM_rm_outlier(K67.monthLUErefRELAX,3);
ix = smooth(K67.monthLUErefRELAX,3);   ix(isnan(K67.monthLUErefRELAX)) = NaN;   K67.monthLUErefRELAX = ix;
K67.monthLUErefRELAX = AM_fill_plus(K67.monthLUErefRELAX);

[K67.week_avgLUErefRELAX,~,K67.week_avgstdLUErefRELAX]   = AM_week2_avg(K67.weekLUErefRELAX,K67.week);
[K67.month_avgLUErefRELAX,~,K67.month_avgstdLUErefRELAX] = AM_month_avg(K67.monthLUErefRELAX,K67.month);
[K67.week_recurrentLUErefRELAX]  = AM_week2dayAVG_recurrent(K67.week_avgLUErefRELAX,K67.week);
[K67.month_recurrentLUErefRELAX] = AM_monthAVG_recurrent(K67.month_avgLUErefRELAX,K67.month);

% ........................................................................
%% Model
% ........................................................................
%Calculates a teoretical PAR
PAR = K67.radiation_plus;
PAR((K67.So<10)|(K67.radiation_plus<0))=0;
[PARmaxtest_cycle,PARmax_cycle] = AM_day_cycle_max(PAR,K67.dateloc);
[~,So_cycle] = AM_day_cycle_max(K67.So,K67.dateloc);
PARmaxtest_cycle = reshape(PARmaxtest_cycle,[12*24,1]);   PARmax_cycle = reshape(PARmax_cycle,[12*24,1]);
So_cycle         = reshape(So_cycle,[12*24,1]);
p_test = polyfit(So_cycle(~isnan(So_cycle)&~isnan(PARmaxtest_cycle)),PARmaxtest_cycle(~isnan(So_cycle)&~isnan(PARmaxtest_cycle)),1);
p      = polyfit(So_cycle(~isnan(So_cycle)&~isnan(PARmax_cycle)),PARmax_cycle(~isnan(So_cycle)&~isnan(PARmax_cycle)),1);
PARtheory = polyval(p,K67.So);          PARtheorytest = polyval(p_test,K67.So);
% % figure('color','white');            plot(PARtheory);   hold on;
% % plot(PARtheorytest,'m');            plot(PAR,'g');
% Calculates cloudiness index
K67.CI = 1-(PAR./PARtheory);            K67.CI(K67.CI<0) = 0;    K67.CI(K67.CI>1) = 1;
K67.CI(K67.So<10) = NaN;
K67.monthCI = AM_month2cycle(K67.CI,K67.dateloc,tresholdCycleDay,4,2);
K67.monthCI = AM_fill(K67.monthCI);


%% VPD
VPDfit = polyval(pVPDam,K67.VPD).*1000;
VPDfit(VPDfit<0) = 0;

VPDam = interp1(K67am.dateloc,K67am.VPD,K67.dateloc);
% % mdl = fitlm(VPD,VPDam,'y ~ x1');
% % m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
% % VPDfit = (m.*VPD)+b;                        VPDfit(VPDfit<0) = 0;

figure('color','white');
subplot(2,2,1:2);                           hold on;
plot(K67.dateloc,K67.VPD.*1000,'color',[0.0,0.4,0.4]);
plot(K67.dateloc,VPDfit,'color',[0.0,0.4,0.9]);
plot(K67.dateloc,VPDam,'color',[0.9,0.4,0.4]);
subplot(2,2,3);                           hold on;
plot(K67.VPD.*1000,VPDam,'.','color',[0.4,0.4,0.4]);  plot(K67.VPD.*1000,VPDfit,'.','color',[0.9,0.4,0.4]);
plot((0:500:4000),(0:500:4000),'--g');
ylabel('VPD_{ameriflux}');      xlabel('VPD_{k67}');

%% ........................................................................
% % K67.LUEref = interp1(K67.week+8,K67.weekLUEref,K67.dateloc,'nearest');
% % K67.LUErefRELAX = interp1(K67.week+8,K67.weekLUErefRELAX,K67.dateloc,'nearest');
% % K67.LUEref_recurrent = interp1(K67.week+8,K67.week_recurrentLUEref,K67.dateloc,'nearest');
% % K67.LUEref_constant  = (K67.dateloc.*0)+nanmean(K67.weekLUEref);
K67.LUEref = interp1(K67.month+15,K67.monthLUEref,K67.dateloc,'nearest');
K67.LUErefRELAX = interp1(K67.month+15,K67.monthLUErefRELAX,K67.dateloc,'nearest');
K67.LUEref_recurrent = interp1(K67.month+15,K67.month_recurrentLUEref,K67.dateloc,'nearest');
K67.LUEref_constant  = (K67.dateloc.*0)+nanmean(K67.monthLUEref);

kCI  = 2.056;                   fENVref = 1.225;
kVPD = 191.1.*10^-6;            PARo    = 6216.4;

K67.GEPmodelRELAX = (K67.LUErefRELAX./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
    (1./(1+(K67.radiation_plus./PARo)));
K67.GEPmodel = (K67.LUEref./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
    (1./(1+(K67.radiation_plus./PARo)));
% % K67.GEPmodelVPD = (K67.LUEref./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*K67.VPD.*1000)).*....
% %     (1./(1+(K67.radiation_plus./PARo)));
K67.GEPmodel_LUErecurrent = (K67.LUEref_recurrent./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
    (1./(1+(K67.radiation_plus./PARo)));
K67.GEPmodel_LUEavg = (K67.LUEref_constant./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
    (1./(1+(K67.radiation_plus./PARo)));

%%
K67.GEPmodelRELAX(K67.So==0)=0;              K67.GEPmodelRELAX(K67.GEPmodelRELAX<-1) = NaN;                 K67.GEPmodelRELAX = AM_fill_plus(K67.GEPmodelRELAX);
K67.GEPmodel(K67.So==0)=0;                   K67.GEPmodel(K67.GEPmodel<-1) = NaN;                           K67.GEPmodel = AM_fill_plus(K67.GEPmodel);
K67.GEPmodel_LUErecurrent(K67.So==0)=0;      K67.GEPmodel_LUErecurrent(K67.GEPmodel_LUErecurrent<-1) = NaN; K67.GEPmodel_LUErecurrent = AM_fill_plus(K67.GEPmodel_LUErecurrent);
K67.GEPmodel_LUEavg(K67.So==0)=0;            K67.GEPmodel_LUEavg(K67.GEPmodel_LUEavg<-1) = NaN;             K67.GEPmodel_LUEavg = AM_fill_plus(K67.GEPmodel_LUEavg);


[K67.weekGEPmodelRELAX,~,~,~,~,K67.week_stdGEPmodelRELAX] = AM_week2cycle_rs(K67.GEPmodelRELAX,K67.dateloc,tresholdCycle,4,2);
K67.weekGEPmodelRELAX   = AM_rm_outlier(K67.weekGEPmodelRELAX,3);               K67.weekGEPmodelRELAX    = AM_fill(K67.weekGEPmodelRELAX);

[K67.monthGEPmodelRELAX,~,~,~,~,K67.month_stdGEPmodelRELAX] = AM_month2cycle(K67.GEPmodelRELAX,K67.dateloc,tresholdCycle,4,2);
K67.monthGEPmodelRELAX((K67.monthGEPmodelRELAX<3)|(K67.monthGEPmodelRELAX>12)) =  NaN;
K67.monthGEPmodelRELAX    = AM_fill(K67.monthGEPmodelRELAX);

[K67.weekGEPmodel,~,~,~,~,K67.week_stdGEPmodel] = AM_week2cycle_rs(K67.GEPmodel,K67.dateloc,tresholdCycle,4,2);
K67.weekGEPmodel   = AM_rm_outlier(K67.weekGEPmodel,3);               K67.weekGEPmodel    = AM_fill(K67.weekGEPmodel);

[K67.monthGEPmodel,~,~,~,~,K67.month_stdGEPmodel] = AM_month2cycle(K67.GEPmodel,K67.dateloc,tresholdCycle,4,2);
K67.monthGEPmodel    = AM_fill(K67.monthGEPmodel);

[K67.monthGEPmodel_LUErecurrent,~,~,~,~,K67.month_stdGEPmodel_LUErecurrent] =...
    AM_month2cycle(K67.GEPmodel_LUErecurrent,K67.dateloc,tresholdCycle,4,2);
% K67.monthGEPmodel_LUErecurrent = AM_rm_outlier(K67.monthGEPmodel_LUErecurrent,3);
K67.monthGEPmodel_LUErecurrent = AM_fill(K67.monthGEPmodel_LUErecurrent);

[K67.monthGEPmodel_LUEavg,~,~,~,~,K67.month_stdGEPmodel_LUEavg] =...
    AM_month2cycle(K67.GEPmodel_LUEavg,K67.dateloc,tresholdCycle,4,2);

% K67.monthGEPmodel_LUEavg = AM_rm_outlier(K67.monthGEPmodel_LUEavg,3);
K67.monthGEPmodel_LUEavg = AM_fill(K67.monthGEPmodel_LUEavg);
[K67.week_recurrentgee]  = AM_week2dayAVG_recurrent(K67.week_avggee,K67.week);
K67.month_avggee         = AM_month_avg(K67.monthgee,K67.month);
[K67.month_recurrentgee] = AM_monthAVG_recurrent(K67.month_avggee,K67.month);

%% ........................................................................
figure('color','white');
v1 = 0.006;        v2 = 0.028;  v3 = v1:((v2-v1)/3):v2;
x1 = datenum(K67.Y(1),1,1);      x2 = datenum(K67.Y(end)+1,1,1);    x3 = datenum(K67.year,1,1)';
subplot(2,1,1);    hold on;
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);           hold on;
plot(K67.week,K67.week_recurrentLUEref,'color',[0 0 0.7],'LineWidth',2)
plot(K67.week,K67.week.*0+nanmean(K67.week_avgLUEref),'color',[0.9 0.5 0.],'LineWidth',2)
plot(K67.week,K67.weekLUEref,'color',[0.0 0.7 0.3],'LineWidth',2)
set(gca,'XTick',x3,'XTickLabel',datestr(x3,'yyyy'),'YTick',v3,'YTickLabel',v3,...
    'Ylim',[v1 v2],'Xlim',[x1 x2],'FontSize',12);
ylabel ({'LUE_{ref}';'(\mumolCO_2/\mumol photons)'});

v1 = 0.008;         v2 = 0.026;  v3 = v1:((v2-v1)/3):v2;
x1 = datenum(K67.Y(1),1,1);      x2 = datenum(K67.Y(end)+1,1,1);    x3 = datenum(K67.year,1,1)';
subplot(2,1,2);     hold on;
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
plot(K67.month,K67.month_recurrentLUEref,'color',[0 0 0.7],'LineWidth',2)
plot(K67.month,K67.month.*0+nanmean(K67.month_avgLUEref),'color',[0.9 0.5 0.],'LineWidth',2)
plot(K67.month,K67.monthLUEref,'color',[0.0 0.7 0.3],'LineWidth',2)
set(gca,'XTick',x3,'XTickLabel',datestr(x3,'yyyy'),'YTick',v3,'YTickLabel',v3,...
    'Ylim',[v1 v2],'Xlim',[x1 x2],'FontSize',12);
ylabel ({'LUE_{ref}';'(\mumolCO_2/mol photons)'});

%% ........................................................................
figure('color','white');
v1 = 3;         v2 = 12;    v3 = v1:((v2-v1)/4):v2;
x1 = datenum(K67.Y(1),1,1);      x2 = datenum(K67.Y(end)+1,1,1);    x3 = datenum(K67.year,1,1)';
subplot(2,1,1);    hold on;
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);           hold on;
plot(K67.week,-K67.week_recurrentgee,'color',[0 0 0.7],'LineWidth',2)
plot(K67.week,-K67.weekgee,'color','k','LineWidth',2)
% plot(K67.week,-K67.weekgee,'color',[0.9 0.5 0.],'LineWidth',2)
plot(K67.week,K67.weekGEPmodel,'color',[0. 0.7 0.0],'LineWidth',2)
set(gca,'XTick',x3,'XTickLabel',datestr(x3,'yyyy'),'YTick',v3,'YTickLabel',v3,...
    'Ylim',[v1 v2],'Xlim',[x1 x2],'FontSize',12);
ylabel ({'GEP';'(\mumolCO_2 m^-^2 s^-^1)'});

v1 = 3;         v2 = 11;    v3 = v1:((v2-v1)/4):v2;
x1 = datenum(K67.Y(1),1,1);         x2 = datenum(K67.Y(end)+1,1,1);    x3 = datenum(K67.year,1,1)';
subplot(2,1,2);     hold on;
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;

% % %     upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
% % %     y_area = [(-K67.monthgee+K67.month_stdgee./1000),K67.month_stdgee./1000.*2];
% % %     h  =  area([K67.month,K67.month],y_area);       set(gca,'Layer','top');
% % %     set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
% % %     set(h(1),'FaceColor','none','EdgeColor','none');
% % %     set(h,'BaseValue',0);                       alpha(h,.5);

hl1 = plot(K67.month,-K67.month_recurrentgee,'color',[0 0.5 0.7],'LineWidth',1);
hl7 = plot(K67.month,(K67.month.*0)-nanmean(K67.monthgee),':','color',[0.5 0.0 0.5],'LineWidth',1);
hl2 = plot(K67.month,-K67.monthgee,'color','k','LineWidth',2);
hl3 = plot(K67.month,K67.monthGEPmodel,'color',[0. 0.7 0.0],'LineWidth',2);
hl4 = plot(K67.month,K67.monthGEPmodel_LUErecurrent,'color',[0.0 0.0 0.7],'LineWidth',2);
hl5 = plot(K67.month,K67.monthGEPmodel_LUEavg,'color',[0.9 0.5 0.0],'LineWidth',2);
hl6 = plot(K67.month,K67.monthGEPmodelRELAX,'color',[0.5 0. 0.0],'LineWidth',2);
set(gca,'XTick',x3,'XTickLabel',datestr(x3,'yyyy'),'YTick',v3,'YTickLabel',v3,...
    'Ylim',[v1 v2],'Xlim',[x1 x2],'FontSize',12);
ylabel ({'GEP';'(\mumolCO_2 m^-^2 s^-^1)'});
legend([hl1,hl2,hl3,hl4,hl5,hl6],'Mean seasonal value','Observations','f_{env}+Monthly LUE_{ref}',...
    'f_{env}+Mean seasonal cycle of LUE_{ref}','f_{env}+Constant LUE_{ref}','f_{env}+Monthly LUE_{refRELAX}','orientation','horizontal');
legend box off;
box on;


%%
iy = K67.GEPmodel_LUEavg;
figure('color','white');
x1 = 0;           x2 = 40;         x3 = x1:((x2-x1)/4):x2;
y1 = 0;           y2 = 40;            y3 = y1:((y2-y1)/4):y2;
subplot(2,2,1);                         hold on;
ind = find((K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
plot(iy(ind),-K67.gee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.gee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.gee)&~isnan(iy)&(K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
r2 = corrcoef(iy(ind),-K67.gee(ind));
text (x1,y2,{['R^2:' num2str(r2(2),'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;                                 grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_Y == ik(ij));
    hl2  =  plot(iy(ind),-K67.gee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.gee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{model} (gC m^-^2 d^-^1)';'f_{env}+Constant LUE_{ref}';},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


iy = K67.GEPmodel_LUErecurrent;
% y1 = 4;           y2 = 12;            y3 = y1:((y2-y1)/4):y2;
subplot(2,2,2);                         hold on;
ind = find((K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
plot(iy(ind),-K67.gee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.gee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.gee)&~isnan(iy)&(K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
r2 = corrcoef(iy(ind),-K67.gee(ind));
text (x1,y2,{['R^2:' num2str(r2(2),'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_Y == ik(ij));
    hl2  =  plot(iy(ind),-K67.gee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.gee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{model} (gC m^-^2 d^-^1)';'f_{env}+Mean seasonal cycle of LUE_{ref}';},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


iy = K67.GEPmodel;
% y1 = 4;           y2 = 12;            y3 = y1:((y2-y1)/4):y2;
subplot(2,2,3);                         hold on;
ind = find((K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
plot(iy(ind),-K67.gee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.gee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.gee)&~isnan(iy)&(K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
r2 = corrcoef(iy(ind),-K67.gee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_Y == ik(ij));
    hl2  =  plot(iy(ind),-K67.gee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.gee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{model} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);

iy = K67.GEPmodel;
subplot(2,2,4);                         hold on;
ind = find((K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
plot(iy(ind),-K67.gee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.gee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.gee)&~isnan(iy)&(K67.HY_Y ~= 2008)&(K67.HY_Y ~= 2015));
r2 = corrcoef(iy(ind),-K67.gee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_Y == ik(ij));
    hl2  =  plot(iy(ind),-K67.gee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.gee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{modelVPD} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);

%%
iy = K67.monthGEPmodel_LUEavg;
figure('color','white');
x1 = 5;           x2 = 11;         x3 = x1:((x2-x1)/3):x2;
y1 = 5;           y2 = 11;            y3 = y1:((y2-y1)/3):y2;
subplot(2,2,1);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2),'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(iy(ind),-K67.monthgee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{model} (gC m^-^2 d^-^1)';'f_{env}+Constant LUE_{ref}';},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


iy = K67.monthGEPmodel_LUErecurrent;
% y1 = 4;           y2 = 12;            y3 = y1:((y2-y1)/4):y2;
subplot(2,2,2);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2),'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(iy(ind),-K67.monthgee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{model} (gC m^-^2 d^-^1)';'f_{env}+Mean seasonal cycle of LUE_{ref}';},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


iy = K67.monthGEPmodel;
% y1 = 4;           y2 = 12;            y3 = y1:((y2-y1)/4):y2;
subplot(2,2,3);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(iy(ind),-K67.monthgee(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(iy(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,iy(ind));                     hold on;
    plot(iy(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
xlabel ({'GEP_{model} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


% ix = -K67.monthgee./K67.monthGEPmodel;
ix = -K67.monthgee-K67.monthGEPmodel;
iy =  K67.monthGsFGmmol./1000;
y1 = 0;           y2 = 1.2;            y3 = y1:((y2-y1)/4):y2;
% x1 = 0.6;         x2 = 1.8;         x3 = x1:((x2-x1)/4):x2;
x1 = -2;         x2 = 2;         x3 = x1:((x2-x1)/4):x2;
subplot(2,2,4);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(ix(ind),iy(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(ix(ind),iy(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(ix)):.5:ceil(nanmax(ix)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
ind = find(~isnan(iy)&~isnan(ix));
r2 = corrcoef(ix(ind),iy(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(ix(ind),iy(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(ix(ind),iy(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,ix(ind));                     hold on;
    plot(ix(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
plot((y3.*0)+1,y3,'-','color',[.5 .5 .5]);
xlabel ('GEP-GEP_{model}','FontSize',12);
ylabel ('Gs (mmol m^-^2 s^-^1)','FontSize',12);

%% ........................................................................
figure('color','white');
% ix = -K67.monthgee./K67.monthGEPmodel;
iy =  K67.monthVPDday;
y1 = 0;           y2 = 2.4;            y3 = y1:((y2-y1)/4):y2;
% x1 = 0.6;         x2 = 1.8;         x3 = x1:((x2-x1)/4):x2;
subplot(2,2,1);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(ix(ind),iy(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(ix(ind),iy(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(ix)):.5:ceil(nanmax(ix)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
ind = find(~isnan(iy)&~isnan(ix));
r2 = corrcoef(ix(ind),iy(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(ix(ind),iy(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(ix(ind),iy(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,ix(ind));                     hold on;
    plot(ix(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
plot((y3.*0)+1,y3,'-','color',[.5 .5 .5]);
xlabel ('GEP/GEP_{model}','FontSize',12);
ylabel ('VPD (kPa)','FontSize',12);


% ix = -K67.monthgee./K67.monthGEPmodel;
iy =  K67.monthPAR;
y1 = 0;           y2 = 1000;            y3 = y1:((y2-y1)/4):y2;
% x1 = 0.6;         x2 = 1.8;         x3 = x1:((x2-x1)/4):x2;
subplot(2,2,2);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(ix(ind),iy(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(ix(ind),iy(ind),'y ~ x1');      ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(ix)):.5:ceil(nanmax(ix)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
ind = find(~isnan(iy)&~isnan(ix));          r2 = corrcoef(ix(ind),iy(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(ix(ind),iy(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(ix(ind),iy(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,ix(ind));                     hold on;
    plot(ix(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
plot((y3.*0)+1,y3,'-','color',[.5 .5 .5]);
xlabel ('GEP/GEP_{model}','FontSize',12);
ylabel ('PAR (\mumol m^-^2 s^-^1)','FontSize',12);

% ix = -K67.monthgee./K67.monthGEPmodel;
iy =  K67.monthLUEref;
y1 = 0;           y2 = .04;            y3 = y1:((y2-y1)/4):y2;
% x1 = 0.6;         x2 = 1.8;         x3 = x1:((x2-x1)/4):x2;
subplot(2,2,3);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(ix(ind),iy(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(ix(ind),iy(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(ix)):.5:ceil(nanmax(ix)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
ind = find(~isnan(iy)&~isnan(ix));          r2 = corrcoef(ix(ind),iy(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(ix(ind),iy(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(ix(ind),iy(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,ix(ind));                     hold on;
    plot(ix(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
plot((y3.*0)+1,y3,'-','color',[.5 .5 .5]);
xlabel ('GEP/GEP_{model}','FontSize',12);
ylabel ({'LUE_{ref}';'(\mumolCO_2 mmol^-^1)'},'FontSize',12);

% ix = -K67.monthgee./K67.monthGEPmodel;
iy =  K67.monthCI;
y1 = 0.2;           y2 = 0.8;            y3 = y1:((y2-y1)/4):y2;
% x1 = 0.6;         x2 = 1.8;         x3 = x1:((x2-x1)/4):x2;
subplot(2,2,4);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(ix(ind),iy(ind),'.','color',[.5 .5 .5],'MarkerSize',10);
mdl = fitlm(ix(ind),iy(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(ix)):.5:ceil(nanmax(ix)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
ind = find(~isnan(iy)&~isnan(ix));          r2 = corrcoef(ix(ind),iy(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ik = [2015;2008];                       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(ix(ind),iy(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(ix(ind),iy(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,ix(ind));                     hold on;
    plot(ix(ind),y_model,'-','Color',color_mtx(count,:));
    %         text (x3(2),y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    %             ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
    text (x1,y2,{['R^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
plot((y3.*0)+1,y3,'-','color',[.5 .5 .5]);
xlabel ('GEP/GEP_{model}','FontSize',12);
ylabel ('CI','FontSize',12);


%% ........................................................................
[K67.week_deseasonLUEref,K67.week_trendLUEref,K67.week_mtrendLUEref] = AM_week2trend(K67.weekLUEref.*1000,K67.week,1,'LUE_{ref}'); hold on;

v1 = 11;                         v2 = 26;       v3 = v1:((v2-v1)/5):v2;
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl1=plot(K67.week,K67.weekLUEref.*1000,'LineWidth',1,'color',[0.4 0.4 0.6]); hold on;
hl2=plot(K67.week,K67.week_deseasonLUEref,'LineWidth',2,'color',[0.6 0.4 0.4]); hold on;
hl3=plot(K67.week,(K67.week.*0)+nanmean(K67.weekLUEref.*1000),':','LineWidth',1,'color',[0.4 0.4 0.6]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',14);
ylabel ('LUE_{ref} (\mumolC mmol)','FontSize',12);  xtickangle(0);

%%
figure('color','white','Position',[50 50 600 800],'units','inches');
v1 = 3;         v2 = 11;    v3 = v1:((v2-v1)/4):v2;
x1 = datenum(K67.Y(1),1,1);         x2 = datenum(K67.Y(end)+1,1,1);    x3 = datenum(K67.year,1,1)';
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');    hold on;
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);    set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;

% % %     upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
% % %     y_area = [(-K67.monthgee+K67.month_stdgee./1000),K67.month_stdgee./1000.*2];
% % %     h  =  area([K67.month,K67.month],y_area);       set(gca,'Layer','top');
% % %     set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
% % %     set(h(1),'FaceColor','none','EdgeColor','none');
% % %     set(h,'BaseValue',0);                       alpha(h,.5);

hl1 = plot(K67.month,-K67.month_recurrentgee,'color',[0 0.5 0.7],'LineWidth',1);
hl7 = plot(K67.month,(K67.month.*0)-nanmean(K67.monthgee),':','color',[0.5 0.0 0.5],'LineWidth',1);
hl2 = plot(K67.month,-K67.monthgee,'color','k','LineWidth',2);
hl3 = plot(K67.month,K67.monthGEPmodel,'color',[0. 0.7 0.0],'LineWidth',2);
hl4 = plot(K67.month,K67.monthGEPmodel_LUErecurrent,'color',[0.0 0.0 0.7],'LineWidth',2);
hl5 = plot(K67.month,K67.monthGEPmodel_LUEavg,'color',[0.9 0.5 0.0],'LineWidth',2);
hl6 = plot(K67.month,K67.monthGEPmodelRELAX,'color',[0.6 0.0 0.6],'LineWidth',1.5);
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',14);
ylabel ({'GEP';'(\mumolCO_2 m^-^2 s^-^1)'});
legend([hl1,hl2,hl3,hl4,hl5,hl6],'Mean seasonal value','Observations','f_{env}+Monthly LUE_{ref}',...
    'f_{env}+Mean seasonal cycle of LUE_{ref}','f_{env}+Constant LUE_{ref}',...
    'f_{env}+Monthly LUE_{refRELAX}','orientation','vertical');
legend box off;     box on;        xtickangle(0);

% ........................................................................
%% sensitivity check
% ........................................................................
percentChange = (0.75:0.05:1.25)';
x1 = 5;           x2 = 11;         x3 = x1:((x2-x1)/3):x2;
y1 = 5;           y2 = 11;            y3 = y1:((y2-y1)/3):y2;
figure('color','white');
iy = K67.monthGEPmodel;

% Biology ...............................................................
subplot(2,2,1);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',5);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
count = 1;
for ix = 1:length(percentChange)
    count = count*0.95;
    multiplyChange = percentChange(ix);
    iw = K67.LUEref.*multiplyChange;  iw(iw<0) = 0;
    GEPmodelChange = (iw./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
        (1./(1+(K67.radiation_plus./PARo)));
    GEPmodelChange(K67.So==0)=0;            GEPmodelChange = AM_fill(GEPmodelChange);
    monthGEPmodelChange = AM_month2cycle(GEPmodelChange,K67.dateloc,tresholdCycle,4,2);
    ik = [2015;2008];
    for ij = 1:length(ik)
        if ik(ij)==2015
            colorChange = [count 0 0];
        elseif ik(ij)==2008
            colorChange = [0 0 count];
        end
        ind = find(K67.HY_monthY == ik(ij));
        hl2  =  plot(monthGEPmodelChange(ind),-K67.monthgee(ind),'.','Color',colorChange,'MarkerSize',10); hold on;
        [p,r2,~,~,sm] = AM_compute_regression(monthGEPmodelChange(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
        y_model = polyval(p,monthGEPmodelChange(ind));                     hold on;
        plot(monthGEPmodelChange(ind),y_model,'-','Color',colorChange);
        %         text (x1,y2,{['R^2_{multiplier:}' num2str(multiplyChange) ': ' num2str(r2(1),'%3.2f')]},...
        %             'FontSize',10,'color',colorChange);
    end
end

xlabel ({'GEP_{model biology} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


% radiation_plus ...............................................................
subplot(2,2,2);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',5);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
count = 1;
for ix = 1:length(percentChange)
    count = count*0.95;
    multiplyChange = percentChange(ix);
    iw = K67.radiation_plus.*multiplyChange; iw(iw<0) = 0;
    GEPmodelChange = (K67.LUEref./fENVref).*iw.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
        (1./(1+(iw./PARo)));
    GEPmodelChange(K67.So==0)=0;            GEPmodelChange = AM_fill(GEPmodelChange);
    monthGEPmodelChange = AM_month2cycle(GEPmodelChange,K67.dateloc,tresholdCycle,4,2);
    ik = [2015;2008];
    for ij = 1:length(ik)
        if ik(ij)==2015
            colorChange = [count 0 0];
        elseif ik(ij)==2008
            colorChange = [0 0 count];
        end
        ind = find(K67.HY_monthY == ik(ij));
        hl2  =  plot(monthGEPmodelChange(ind),-K67.monthgee(ind),'.','Color',colorChange,'MarkerSize',10); hold on;
        [p,r2,~,~,sm] = AM_compute_regression(monthGEPmodelChange(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
        y_model = polyval(p,monthGEPmodelChange(ind));                     hold on;
        plot(monthGEPmodelChange(ind),y_model,'-','Color',colorChange);
        %         text (x1,y2,{['R^2_{multiplier:}' num2str(multiplyChange) ': ' num2str(r2(1),'%3.2f')]},...
        %             'FontSize',10,'color',colorChange);
    end
end

xlabel ({'GEP_{model radiation_plus} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


% radiation_plus ...............................................................
subplot(2,2,3);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',5);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
count = 1;
for ix = 1:length(percentChange)
    count = count*0.95;
    multiplyChange = percentChange(ix);
    iw = VPDfit.*multiplyChange;        iw(iw<0) = 0;
    GEPmodelChange = (K67.LUEref./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*(iw))).*....
        (1./(1+(K67.radiation_plus./PARo)));
    GEPmodelChange(K67.So==0)=0;            GEPmodelChange = AM_fill(GEPmodelChange);
    monthGEPmodelChange = AM_month2cycle(GEPmodelChange,K67.dateloc,tresholdCycle,4,2);
    ik = [2015;2008];
    for ij = 1:length(ik)
        if ik(ij)==2015
            colorChange = [count 0 0];
        elseif ik(ij)==2008
            colorChange = [0 0 count];
        end
        ind = find(K67.HY_monthY == ik(ij));
        hl2  =  plot(monthGEPmodelChange(ind),-K67.monthgee(ind),'.','Color',colorChange,'MarkerSize',10); hold on;
        [p,r2,~,~,sm] = AM_compute_regression(monthGEPmodelChange(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
        y_model = polyval(p,monthGEPmodelChange(ind));                     hold on;
        plot(monthGEPmodelChange(ind),y_model,'-','Color',colorChange);
        %         text (x1,y2,{['R^2_{multiplier:}' num2str(multiplyChange) ': ' num2str(r2(1),'%3.2f')]},...
        %             'FontSize',10,'color',colorChange);
    end
end

xlabel ({'GEP_{model VPD} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


% CI ...............................................................
subplot(2,2,4);                         hold on;
ind = find((K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
plot(iy(ind),-K67.monthgee(ind),'.','color',[.5 .5 .5],'MarkerSize',5);
mdl = fitlm(iy(ind),-K67.monthgee(ind),'y ~ x1');
ci  = coefCI(mdl,0.05);
m = mdl.Coefficients.Estimate(2);           b = mdl.Coefficients.Estimate(1);
xValues = (floor(nanmin(iy)):.5:ceil(nanmax(iy)))';
upper_line= (ci(2,1).*xValues)+ci(1,2);     lower_line = (ci(2,2).*xValues)+ci(1,1);
y_area = [lower_line,upper_line-lower_line];
h  =  area([xValues,xValues],y_area);       set(gca,'Layer','top');
set(h(2),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);                       alpha(h,.5);
y_model = (m.*xValues)+b;
plot(xValues,y_model,'-','color',[0.4,0.4,0.4]);
plot(xValues,xValues,'--','color',[0.4,0.4,0.4]);
ind = find(~isnan(-K67.monthgee)&~isnan(iy)&(K67.HY_monthY ~= 2008)&(K67.HY_monthY ~= 2015));
r2 = corrcoef(iy(ind),-K67.monthgee(ind));
text (x1,y2,{['R^2:' num2str(r2(2)^2,'%3.2f')]},'FontSize',10);
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
count = 1;
for ix = 1:length(percentChange)
    count = count*0.95;
    multiplyChange = percentChange(ix);
    iw = K67.CI.*multiplyChange;        iw(iw<0) = 0; iw(iw>1) = 1;
    GEPmodelChange = (K67.LUEref./fENVref).*K67.radiation_plus.*(1+(kCI.*(iw))).*(1-(kVPD.*VPDfit)).*....
        (1./(1+(K67.radiation_plus./PARo)));
    GEPmodelChange(K67.So==0)=0;            GEPmodelChange = AM_fill(GEPmodelChange);
    monthGEPmodelChange = AM_month2cycle(GEPmodelChange,K67.dateloc,tresholdCycle,4,2);
    ik = [2015;2008];
    for ij = 1:length(ik)
        if ik(ij)==2015
            colorChange = [count 0 0];
        elseif ik(ij)==2008
            colorChange = [0 0 count];
        end
        ind = find(K67.HY_monthY == ik(ij));
        hl2  =  plot(monthGEPmodelChange(ind),-K67.monthgee(ind),'.','Color',colorChange,'MarkerSize',10); hold on;
        [p,r2,~,~,sm] = AM_compute_regression(monthGEPmodelChange(ind),-K67.monthgee(ind),[0.4,0.4,0.4],0);
        y_model = polyval(p,monthGEPmodelChange(ind));                     hold on;
        plot(monthGEPmodelChange(ind),y_model,'-','Color',colorChange);
        %         text (x1,y2,{['R^2_{multiplier:}' num2str(multiplyChange) ': ' num2str(r2(1),'%3.2f')]},...
        %             'FontSize',10,'color',colorChange);
    end
end

xlabel ({'GEP_{model cloudiness} (gC m^-^2 d^-^1)';'f_{env}+Monthly LUE_{ref}'},'FontSize',12);
ylabel ({'GEP (gC m^-^2 d^-^1)'},'FontSize',12);


%% Testing the model using the 16-day average
% % [K67.monthradiation_plus,~,~,~,~,~,~,K67.month_cycleradiation_plus,K67.month_cycledateloc] = AM_month2cycle(K67.radiation_plus,K67.dateloc,tresholdCycleDay,4,2);
% % [K67.monthVPD,~,~,~,~,~,~,K67.month_cycleVPD] = AM_month2cycle(K67.VPD,K67.dateloc,tresholdCycleDay,4,2);
% % [K67.monthCI,~,~,~,~,~,~,K67.month_cycleCI]   = AM_month2cycle(K67.CI,K67.dateloc,tresholdCycleDay,4,2);
% % [~,~,~,~,~,~,~,K67.month_cycleLUEref]         = AM_month2cycle(K67.LUEref,K67.dateloc,tresholdCycleDay,4,2);
% % [~,~,~,~,~,~,~,K67.month_cycleGEP]            = AM_month2cycle(-K67.gee,K67.dateloc,tresholdCycleDay,4,2);

ind = find((K67.HY_Y~=2008)&(K67.HY_Y~=2015));
[K67.month_cycleradiation_plus] = AM_day_cycle(K67.radiation_plus(ind),K67.dateloc(ind),6,18);
[K67.month_cycleCI]        = AM_day_cycle(K67.CI(ind),K67.dateloc(ind),6,18);
[K67.month_cycleVPD]       = AM_day_cycle(VPDfit(ind),K67.dateloc(ind),6,18);
[K67.month_cycleLUEref]    = AM_day_cycle(K67.LUEref(ind),K67.dateloc(ind),6,18);
[K67.month_cycleGEP]       = AM_day_cycle(-K67.gee(ind),K67.dateloc(ind),6,18);

hour_avg = (1/24:1/24:1)';
figure('color','white');
newcolors = jet(12);
colororder(newcolors)
subplot(2,3,1); plot(cycle_hour,K67.month_cycleradiation_plus','LineWidth',2);ylabel('PAR (\mumol m^-^2 s^-^1)','FontSize',12);       
set(gca,'XTick',cycle_hour(1:6:24),'YTick',0:400:1600,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});	grid on;
axis([0 24 0 1600]);      xlabel('Hour','FontSize',12);      xtickangle(0);
subplot(2,3,2); plot(cycle_hour,K67.month_cycleCI','LineWidth',2);            ylabel('CI','FontSize',12);                       
set(gca,'XTick',cycle_hour(1:6:24),'YTick',0:0.2:1,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});	grid on;
axis([0 24 0 1]);      xlabel('Hour','FontSize',12);      xtickangle(0);
subplot(2,3,3); plot(cycle_hour,K67.month_cycleVPD','LineWidth',2);           ylabel('VPD (kPa)','FontSize',12);                
set(gca,'XTick',cycle_hour(1:6:24),'YTick',0:500:2000,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});	grid on;
axis([0 24 0 2000]);      xlabel('Hour','FontSize',12);      xtickangle(0);
subplot(2,3,4); plot(cycle_hour,K67.month_cycleLUEref'.*1000,'LineWidth',2);        ylabel('LUE_{ref} (/mumol_{CO_{2}} /mol^-^1)','FontSize',12); 
set(gca,'XTick',cycle_hour(1:6:24),'YTick',16:4:24,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});	grid on;
axis([0 24 16 24]);      xlabel('Hour','FontSize',12);      xtickangle(0);
subplot(2,3,5); plot(cycle_hour,K67.month_cycleGEP','LineWidth',2);           ylabel('GEP (\mumol_{CO_{2}} m^-^2 s^-^1)','FontSize',12);  
set(gca,'XTick',cycle_hour(1:6:24),'YTick',0:10:30,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});	grid on;
axis([0 24 0 30]);      xlabel('Hour','FontSize',12);      xtickangle(0);

K67.month_cycleradiation_plus = reshape(K67.month_cycleradiation_plus',[12*24,1]);
K67.month_cycleCI        = reshape(K67.month_cycleCI',[12*24,1]);
K67.month_cycleVPD       = reshape(K67.month_cycleVPD',[12*24,1]);
K67.month_cycleLUEref    = reshape(K67.month_cycleLUEref',[12*24,1]);
K67.month_cycleGEP       = reshape(K67.month_cycleGEP',[12*24,1]);

K67.month_cycleGEPmodel = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel((isnan(K67.month_cycleGEPmodel))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
ix = reshape(K67.month_cycleGEPmodel,[24,12]);

subplot(2,3,6); plot(cycle_hour,ix,'LineWidth',2); ylabel('GEP_{model} (\mumol_{CO_{2}} m^-^2 s^-^1)','FontSize',12)
set(gca,'XTick',cycle_hour(1:6:24),'YTick',0:10:30,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});	grid on;
axis([0 24 0 30]);      xlabel('Hour','FontSize',12);      xtickangle(0);

%% ........................................................................
cycle_hour = [];
for ik = 1:12
    cycle_hour = [cycle_hour;datenum(0,ik,1,hour_avg,0,0)];
end

K67.month_avgGEPmodelcycle = AM_month_avg(K67.month_cycleGEPmodel,cycle_hour);
K67.month_avgGEPcycle      = AM_month_avg(K67.month_cycleGEP,cycle_hour);
ind = find((K67.HY_monthY~=2008)&(K67.HY_monthY~=2015));
K67.month_avgGEPmodel      = AM_month_avg(K67.monthGEPmodel(ind),K67.month(ind));
K67.monthradiation_plus = AM_month2cycle(K67.radiation_plus,K67.dateloc,tresholdCycle,4,2);

%% ........................................................................
ind = find(K67.HY_Y==2008);
K67.month_cycleradiation_plus2008 = AM_day_cycle(K67.radiation_plus(ind),K67.dateloc(ind),6,18);
K67.month_cycleCI2008        = AM_day_cycle(K67.CI(ind),K67.dateloc(ind),6,18);
K67.month_cycleVPD2008       = AM_day_cycle(VPDfit(ind),K67.dateloc(ind),6,18);
K67.month_cycleLUEref2008    = AM_day_cycle(K67.LUErefRELAX(ind),K67.dateloc(ind),6,18);
K67.month_cycleGEP2008       = AM_day_cycle(-K67.gee(ind),K67.dateloc(ind),6,18);

K67.month_cycleradiation_plus2008 = reshape(K67.month_cycleradiation_plus2008',[12*24,1]);
K67.month_cycleCI2008 = reshape(K67.month_cycleCI2008',[12*24,1]);
K67.month_cycleVPD2008 = reshape(K67.month_cycleVPD2008',[12*24,1]);
K67.month_cycleLUEref2008 = reshape(K67.month_cycleLUEref2008',[12*24,1]);
K67.month_cycleGEP2008 = reshape(K67.month_cycleGEP2008',[12*24,1]);

% % K67.GEPmodelRELAX = (K67.LUErefRELAX./fENVref).*K67.radiation_plus.*(1+(kCI.*K67.CI)).*(1-(kVPD.*VPDfit)).*....
% %     (1./(1+(K67.radiation_plus./PARo)));
% % K67.month_cycleGEPmodel2008radiation_plus = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus2008.*...
% %     (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD)).*....
% %     (1./(1+(K67.month_cycleradiation_plus2008./PARo)));

K67.month_cycleGEPmodel2008radiation_plus = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus2008.*...
    (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD)).*....
    (1./(1+(K67.month_cycleradiation_plus2008./PARo)));
K67.month_cycleGEPmodel2008radiation_plus((isnan(K67.month_cycleGEPmodel2008radiation_plus))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2008cyclePAR = AM_month_avg(K67.month_cycleGEPmodel2008radiation_plus,cycle_hour);

K67.month_cycleGEPmodel2008vpd = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD2008)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2008vpd((isnan(K67.month_cycleGEPmodel2008vpd))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2008cycleVPD = AM_month_avg(K67.month_cycleGEPmodel2008vpd,cycle_hour);

K67.month_cycleGEPmodel2008ci = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI2008)).*(1-(kVPD.*K67.month_cycleVPD)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2008ci((isnan(K67.month_cycleGEPmodel2008ci))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2008cycleCI = AM_month_avg(K67.month_cycleGEPmodel2008ci,cycle_hour);

K67.month_cycleGEPmodel2008lue = (K67.month_cycleLUEref2008./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI2008)).*(1-(kVPD.*K67.month_cycleVPD)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2008lue((isnan(K67.month_cycleGEPmodel2008lue))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2008cycleLUEref = AM_month_avg(K67.month_cycleGEPmodel2008lue,cycle_hour);
K67.month_avgGEPmodel2008cycleLUEref(K67.month_avgGEPmodel2008cycleLUEref<3) = NaN;
K67.month_avgGEPmodel2008cycleLUEref = AM_fill(K67.month_avgGEPmodel2008cycleLUEref);

%% ........................................................................
v1 = 2;         v2 = 11;    v3 = v1:((v2-v1)/3):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avggee-K67.week_avgstdgee,2.*K67.week_avgstdgee];
y = -y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y.*gC_m2_d1);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl1  =  plot(week_avg,-K67.weekgee(K67.HY_weekY == ik).*gC_m2_d1,'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl3= plot(month_avg,-K67.month_avggee([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',2);
hl4  =  plot(week_avg,-K67.weekgee(K67.HY_weekY == 2008).*gC_m2_d1,'Color',[.0 .0 .0],'LineWidth',2); hold on;
hl5 = plot(month_avg,K67.month_avgGEPmodel([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0 0 0],'LineStyle',':','LineWidth',2);
hl6 = plot(month_avg,K67.month_avgGEPmodelcycle([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.2 0.2 0.2],'LineStyle','-.','LineWidth',2);
hl7 = plot(month_avg,K67.month_avgGEPmodel2008cyclePAR([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.5 0.2 0.9],'LineStyle','--','LineWidth',2);
hl8 = plot(month_avg,K67.month_avgGEPmodel2008cycleVPD([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.5 0.2 0.9],'LineStyle','-.','LineWidth',2);
hl9 = plot(month_avg,K67.month_avgGEPmodel2008cycleCI([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.5 0.2 0.9],'LineStyle',':','LineWidth',2);
hl10 = plot(month_avg,K67.month_avgGEPmodel2008cycleLUEref([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.5 0.2 0.9],'LineStyle','-','LineWidth',2);
hl11 =  plot(week_avg,K67.weekGEPmodelRELAX(K67.HY_weekY == 2008).*gC_m2_d1,'Color',[.0 .7 .0],'LineWidth',2); hold on;
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',12);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'GEP';'(gC m^-^2 d^-^1)'},'FontSize',12);        box on;
legend([hl3,hl4,hl5,hl6,hl7,hl8,hl9,hl10,hl11],'average all observations','observations HY2008','model','model mean annual cycle',...
    'PAR_{HY2008}','VPD_{HY2008}','CI_{HY2008}','LUE_{ref HY2008}','model_{HY2008}','FontSize',10);
legend box off; xtickangle(0);

%% ........................................................................
ind = find(K67.HY_Y==2015);
K67.month_cycleradiation_plus2015 = AM_day_cycle(K67.radiation_plus(ind),K67.dateloc(ind),6,18);
K67.month_cycleCI2015        = AM_day_cycle(K67.CI(ind),K67.dateloc(ind),6,18);
K67.month_cycleVPD2015       = AM_day_cycle(VPDfit(ind),K67.dateloc(ind),6,18);
K67.month_cycleLUEref2015    = AM_day_cycle(K67.LUErefRELAX(ind),K67.dateloc(ind),6,18);
K67.month_cycleGEP2015       = AM_day_cycle(-K67.gee(ind),K67.dateloc(ind),6,18);

K67.month_cycleradiation_plus2015 = reshape(K67.month_cycleradiation_plus2015',[12*24,1]);
K67.month_cycleCI2015 = reshape(K67.month_cycleCI2015',[12*24,1]);
K67.month_cycleVPD2015 = reshape(K67.month_cycleVPD2015',[12*24,1]);
K67.month_cycleLUEref2015 = reshape(K67.month_cycleLUEref2015',[12*24,1]);
K67.month_cycleGEP2015 = reshape(K67.month_cycleGEP2015',[12*24,1]);

K67.month_cycleGEPmodel2015radiation_plus = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2015radiation_plus((isnan(K67.month_cycleGEPmodel2015radiation_plus))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2015cyclePAR = AM_month_avg(K67.month_cycleGEPmodel2015radiation_plus,cycle_hour);

K67.month_cycleGEPmodel2015vpd = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD2015)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2015vpd((isnan(K67.month_cycleGEPmodel2015vpd))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2015cycleVPD = AM_month_avg(K67.month_cycleGEPmodel2015vpd,cycle_hour);

K67.month_cycleGEPmodel2015ci = (K67.month_cycleLUEref./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI2015)).*(1-(kVPD.*K67.month_cycleVPD)).*....
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2015ci((isnan(K67.month_cycleGEPmodel2015ci))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2015cycleCI = AM_month_avg(K67.month_cycleGEPmodel2015ci,cycle_hour);

K67.month_cycleGEPmodel2015lue = (K67.month_cycleLUEref2015./fENVref).*K67.month_cycleradiation_plus.*...
    (1+(kCI.*K67.month_cycleCI)).*(1-(kVPD.*K67.month_cycleVPD)).*....￼
    (1./(1+(K67.month_cycleradiation_plus./PARo)));
K67.month_cycleGEPmodel2015lue((isnan(K67.month_cycleGEPmodel2015lue))|(K67.month_cycleradiation_plus<=tresholdSo)) = 0;
K67.month_avgGEPmodel2015cycleLUEref = AM_month_avg(K67.month_cycleGEPmodel2015lue,cycle_hour);


%% ........................................................................
v1 = 2;         v2 = 11;    v3 = v1:((v2-v1)/3):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avggee-K67.week_avgstdgee,2.*K67.week_avgstdgee];
y = -y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y.*gC_m2_d1);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl1  =  plot(week_avg,-K67.weekgee(K67.HY_weekY == ik).*gC_m2_d1,'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl3 = plot(month_avg,-K67.month_avggee([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',2);
hl4  =  plot(week_avg,-K67.weekgee(K67.HY_weekY == 2015).*gC_m2_d1,'Color',[.0 .0 .0],'LineWidth',2); hold on;
hl5 = plot(month_avg,K67.month_avgGEPmodel([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0 0 0],'LineStyle',':','LineWidth',2);
hl6 = plot(month_avg,K67.month_avgGEPmodelcycle([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.2 0.2 0.2],'LineStyle','-.','LineWidth',2);
hl7 = plot(month_avg,K67.month_avgGEPmodel2015cyclePAR([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.9 0.2 0.7],'LineStyle','--','LineWidth',2);
hl8 = plot(month_avg,K67.month_avgGEPmodel2015cycleVPD([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.9 0.2 0.7],'LineStyle','-.','LineWidth',2);
hl9 = plot(month_avg,K67.month_avgGEPmodel2015cycleCI([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.9 0.2 0.7],'LineStyle',':','LineWidth',2);
hl10 = plot(month_avg,K67.month_avgGEPmodel2015cycleLUEref([month_start_hydro:12,1:(month_start_hydro-1)]).*gC_m2_d1,'Color',[0.9 0.2 0.7],'LineStyle','-','LineWidth',2);
hl11 =  plot(week_avg,K67.weekGEPmodelRELAX(K67.HY_weekY == 2015).*gC_m2_d1,'Color',[.0 .7 .0],'LineWidth',2); hold on;
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',12);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'GEP';'(gC m^-^2 d^-^1)'},'FontSize',12);        box on;
%     text (week_avg(2),v2-10,site_name,'FontSize',12);
legend([hl3,hl4,hl5,hl6,hl7,hl8,hl9,hl10,hl11],'average all observations','observations HY2015','model','model mean annual cycle',...
    'PAR_{HY2015}','VPD_{HY2015}','CI_{HY2015}','LUE_{ref HY2015}','model_{HY2015}','Location','south','FontSize',10);
legend box off;     xtickangle(0);