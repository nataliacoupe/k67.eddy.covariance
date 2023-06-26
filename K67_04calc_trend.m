% Linear trend -- long time scales -- removes seasonality
% Natalia Restrepo-Coupe  
% Toronto
% .........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;

month_name = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];
% % 
% % % .........................................................................
% % %% Just removing the average  
% % % .........................................................................
% % g = length(unique(CERES.y));
% % month_trendCERESnr   = NaN(length(CERES.m),1);       month_trendCERESswd   = NaN(length(CERES.m),1);
% %     for ik = 1:length(CERES.m)
% %         iq = CERES.m(ik);
% %     A  = CERES.month_avgNR(iq);
% %     C  = monthCERESnr(ik);
% %     D  = (C-A);
% %     month_trendCERESnr(ik) = D;
% % 
% %     A  = CERES.month_avgSWdown(iq);
% %     C  = monthCERESswd(ik);
% %     D  = (C-A);
% %     month_trendCERESswd(ik) = D;
% % end
% % p  = polyfit(monthCERESdate(~isnan(month_trendCERESswd)),month_trendCERESswd(~isnan(month_trendCERESswd)),1);
% % y   = polyval(p,monthCERESdate);
% % figure('color','white');
% % plot(monthCERESdate,month_trendCERESswd);                    hold on;
% % plot(monthCERESdate,y,'.','color',[0.9 .5 .5]);
% % p  = polyfit(monthCERESdate(~isnan(month_trendCERESnr)),month_trendCERESnr(~isnan(month_trendCERESswd)),1);
% % y   = polyval(p,monthCERESdate);
% % plot(monthCERESdate,month_trendCERESnr);                    hold on;
% % plot(monthCERESdate,y,'.','color',[0.5 .9 .5]);
% % 
% % K67.month_trendCERESswd = interp1(monthCERESdate,month_trendCERESswd,K67.month);
% % K67.month_trendCERESnr  = interp1(monthCERESdate,month_trendCERESnr,K67.month);

% .........................................................................
% Seasonal Adjustment Using a Stable Seasonal Filter
% .........................................................................

% .........................................................................
%% Moving Average Trend Estimation
% .........................................................................
y = monthCERESswd;                      T = length(y);

figure('color','white');
subplot(3,1,1);     plot(y);            h1 = gca;           hold on
h1.XLim = [0,T];    h1.XTick = 1:12:T;
h1.XTickLabel = datestr(monthCERESdate(1:12:T),10);
title 'Shortwave radiation CERES';      ylabel 'SW_{down} (W m^-^2)';

sW13 = [1/24;repmat(1/12,11,1);1/24];   yS = conv(y,sW13,'same');
yS(1:6) = yS(7); yS(T-5:T) = yS(T-6);   xt = y-yS;

p  = polyfit(monthCERESdate(~isnan(yS)),yS(~isnan(yS)),1);   display(p);
ylinear   = polyval(p,monthCERESdate);

hl1 = plot(yS,'r','LineWidth',2);         hl2 = plot(ylinear,'--','LineWidth',2);
legend([hl1,hl2],'Annual Moving Average','Linear trend');     hold off
K67.month_trendCERESswd = interp1(monthCERESdate,ylinear,K67.month);
K67.month_deseasonCERESswd = interp1(monthCERESdate,yS,K67.month);

% .........................................................................
y = monthCERESnr;                      T = length(y);

subplot(3,1,2);     plot(y);            h1 = gca;           hold on
h1.XLim = [0,T];    h1.XTick = 1:12:T;
h1.XTickLabel = datestr(monthCERESdate(1:12:T),10);
title 'Net radiation CERES';      ylabel 'Rn (W m^-^2)';

sW13 = [1/24;repmat(1/12,11,1);1/24];   yS = conv(y,sW13,'same');
yS(1:6) = yS(7); yS(T-5:T) = yS(T-6);   xt = y-yS;

p  = polyfit(monthCERESdate(~isnan(yS)),yS(~isnan(yS)),1);   display(p);
ylinear   = polyval(p,monthCERESdate);

hl1 = plot(yS,'r','LineWidth',2);         hl2 = plot(ylinear,'--','LineWidth',2);
legend([hl1,hl2],'Annual Moving Average','Linear trend');     hold off
K67.month_trendCERESnr    = interp1(monthCERESdate,ylinear,K67.month);
K67.month_deseasonCERESnr = interp1(monthCERESdate,yS,K67.month);

%% .........................................................................
y = TRMM1h.weekprecip;                  T = length(y);

subplot(3,1,3);     plot(y);            h1 = gca;           hold on
h1.XLim = [0,T];    h1.XTick = 1:23:T;
h1.XTickLabel = datestr(TRMM1h.week(1:23:T),10);
title 'Precipitation TRMM';      ylabel 'Precip (mm)';

sW13 = [1/46;repmat(1/23,22,1);1/46];   yS = conv(y,sW13,'same');
yS(1:6) = yS(7); yS(T-5:T) = yS(T-6);   xt = y-yS;

p  = polyfit(TRMM1h.week(~isnan(yS)),yS(~isnan(yS)),1);   display(p);
ylinear   = polyval(p,TRMM1h.week);

hl1 = plot(yS,'r','LineWidth',2);         hl2 = plot(ylinear,'--','LineWidth',2);
legend([hl1,hl2],'Annual Moving Average','Linear trend');     hold off
K67.week_trendTRMM = interp1(TRMM1h.week,ylinear,K67.week);
K67.week_deseasonTRMM = interp1(TRMM1h.week,yS,K67.week);

% .........................................................................
% .........................................................................
%% Other test
%..........................................................................
% .........................................................................
y = AM_fill_plus(K67.weekPc);              T = length(y);

figure('color','white');
hl1 = plot(y);      h1 = gca;       hold on;
h1.XLim = [0,T];    h1.XTick = 1:23:T;
h1.XTickLabel = datestr(K67.week(1:23:T),10);

% Seasonal indices ........................................................
sW13 = [1/46;repmat(1/23,22,1);1/46];
yS = conv(y,sW13,'same');           
yS(1:23) = yS(22);              yS(T-22:T) = yS(T-23);
xt = y-yS;
hl2 = plot(yS,'r','LineWidth',2);     %hold off

% Moving average ..........................................................
s = 23;
sidx = cell(s,1);
for i = 1:s
 sidx{i,1} = i:s:T;
end
% sidx{1:2}

%Apply a stable seasonal filter............................................
sst = cellfun(@(x) nanmean(xt(x)),sidx);

% Put smoothed values back into a vector of length N
nc = floor(T/s); % no. complete years
rm = mod(T,s);   % no. extra months
sst = [repmat(sst,nc,1);sst(1:rm)];

% Center the seasonal estimate (additive)
sBar = nanmean(sst); % for centering
sst = sst-sBar;

hl3 = plot(sst);

% Deseasonalize the series.................................................
dt = y - sst;

hl4 = plot(dt,'LineWidth',1);
hl5 = plot(y*0+(nanmean(y)),':k','LineWidth',1);
legend([hl1,hl2,hl3,hl4,hl5],'Pc','Annual Moving Average',...
    'Stable Seasonal Component','Deseasonalized Series','Mean 16-day precipitation');
ylabel 'Pc (gCO_2 m^-^2 d^-^1)';

% .........................................................................
%% Other
% .........................................................................
z  = AM_fill_plus(-K67.weekgee);    T = length(y);
y = AM_fill_month_cycle(-K67.weekgee,K67.week,6,18);     T = length(y);
yS = conv(y,sW13,'same');           
yS(1:23) = yS(22);                  yS(T-22:T) = yS(T-23);
yS(isnan(z)) = NaN;                 xt = y-yS;
sst = cellfun(@(x) nanmean(xt(x)),sidx);
sst = [repmat(sst,nc,1);sst(1:rm)];
sst(isnan(z)) = NaN;                sBar = nanmean(sst); % for centering
sst = sst-sBar;                     dt = y - sst;
p  = polyfit(K67.week(~isnan(yS)),yS(~isnan(yS)),1);   display(p);
ylinear   = polyval(p,K67.week);
K67.week_trendGPP = ylinear;        K67.week_deseasonGPP = yS;
%%
figure('color','white');
hl1 = plot(z);      h1 = gca;       hold on;
h1.XLim = [0,T];    h1.XTick = 1:23:T;
h1.XTickLabel = datestr(K67.week(1:23:T),10);
hl2 = plot(yS,'r','LineWidth',2);     %hold off
% hl6 = plot(ySsmooth,'c','LineWidth',2);     %hold off
hl3 = plot(sst);
hl4 = plot(dt,'LineWidth',1);
hl5 = plot(y*0+(nanmean(y)),':k','LineWidth',1);
legend([hl1,hl2,hl3,hl4,hl5],'Pc','Annual Moving Average',...
    'Stable Seasonal Component','Deseasonalized Series','Mean 16-day precipitation');
ylabel 'GEP (gCO_2 m^-^2 d^-^1)';

[K67.week_deseasongee,K67.week_trendgee,K67.week_mtrendgee] = AM_week2trend(K67.weekgee,K67.week,1,'GEE');
[K67.week_deseasonPc,K67.week_trendPc,K67.week_mtrendPc]    = AM_week2trend(K67.weekPc,K67.week,1,'Pc');
[K67.week_deseasonr,K67.week_trendr,K67.week_mtrendr]       = AM_week2trend(K67.weekr,K67.week,1,'Reco');
[K67.week_deseasonnee,K67.week_trendnee,K67.week_mtrendnee] = AM_week2trend(K67.weeknee,K67.week,1,'NEE');
[K67.week_deseasonET,K67.week_trendET,K67.week_mtrendET]    = AM_week2trend(K67.weekET,K67.week,1,'ET');
[K67.week_deseasonETday,K67.week_trendETday,K67.week_mtrendETday]    = AM_week2trend(K67.weekETday,K67.week,1,'ETday');
[K67.week_deseasonETrefday,K67.week_trendETrefday,K67.week_mtrendETrefday]    = AM_week2trend(K67.weekETrefday,K67.week,1,'ETday&ref');
[K67.week_deseasonH,K67.week_trendH,K67.week_mtrendH]       = AM_week2trend(K67.weekH,K67.week,1,'H');
[K67.week_deseasonNetRad,K67.week_trendNetRad,K67.week_mtrendNetRad]    = AM_week2trend(K67.weekNetRad,K67.week,1,'Rn');
[K67.week_deseasonTaspday,K67.week_trendTaspday,K67.week_mtrendTaspday] = AM_week2trend(K67.weekTaspday,K67.week,1,'T_{air day}');
[K67.week_deseasonVPDday,K67.week_trendVPDday,K67.week_mtrendVPDday]    = AM_week2trend(K67.weekVPDday,K67.week,1,'VPD_{day}');
[K67.week_deseasonLEday,K67.week_trendLEday,K67.week_mtrendLEday] = AM_week2trend(K67.weekLEday,K67.week,1,'LE_{day}');
[K67.week_deseasonHday,K67.week_trendHday,K67.week_mtrendHday] = AM_week2trend(K67.weekHday,K67.week,1,'H_{day}');
[K67.week_deseasonBowenday,K67.week_trendBowenday,K67.week_mtrendBowenday] = AM_week2trend(K67.weekBowenday,K67.week,1,'Bowen');
[K67.week_deseasonGsFGmmol,K67.week_trendGsFGmmol,K67.week_mtrendGsFGmmol] = AM_week2trend(K67.weekGsFGmmol,K67.week,1,'Gs');
[K67.week_deseasonWUEdryday,K67.week_trendWUEdryday,K67.week_mtrendWUEdryday] = AM_week2trend(K67.weekWUEdryday,K67.week,1,'WUE_{dry&day}');
[K67.week_deseasonLAInormal,K67.week_trendLAInormal,K67.week_mtrendLAInormal] = AM_week2trend(K67.weekLAInormal,K67.week,1,'LAI');
[K67.week_deseasonFlush,K67.week_trendFlush,K67.week_mtrendFlush] = AM_week2trend(K67.weekFlush,K67.week,1,'Flush'); hold on;
[K67.week_recurrentFlush] = AM_week2dayAVG_recurrent(K67.week_avgFlush,K67.week); plot(K67.week_recurrentFlush,'k')
[K67.week_recurrentLAI] = AM_week2dayAVG_recurrent(K67.week_avgLAI,K67.week); plot(K67.week_recurrentLAI,'k')

[K67.week_deseasonLeafLitter,K67.week_trendLeafLitter,K67.week_mtrendLeafLitter] = AM_week2trend(K67.weekLeafLitterFill,K67.week,1,'NPP_{leaf litter}'); hold on;
[K67.week_deseasonAllLitter,K67.week_trendAllLitter,K67.week_mtrendAllLitter] = AM_week2trend(K67.weekAllLitter,K67.week,1,'NPP_{litter}'); hold on;
[K67.week_deseasonLAI,K67.week_trendLAI,K67.week_mtrendLAI] = AM_week2trend(K67.weekLAI,K67.week,1,'LAI'); hold on;
[K67.week_deseasonCWD,K67.week_trendCWD,K67.week_mtrendCWD] = AM_week2trend(K67.weekCWD,K67.week,1,'CWD'); hold on;


[K67.week_deseasonTRMM,K67.week_trendTRMM,K67.week_mtrendTRMM] = AM_week2trend(K67.weekTRMM,K67.week,1,'TRMM'); hold on;
[K67.week_deseasonLST,K67.week_trendLST,K67.week_mtrendLST] = AM_week2trend(K67.weekLST,K67.week,1,'LST'); hold on;
[K67.week_deseasonLSTnight,K67.week_trendLSTnight,K67.week_mtrendLSTnight] = AM_week2trend(K67.weekLSTnight,K67.week,1,'LSTnight'); hold on;

[K67.week_deseasonLAIfpar,K67.week_trendLAIfpar,K67.week_mtrendLAIfpar] = AM_week2trend(K67.weekLAIfpar,K67.week,1,'LAIfpar'); hold on;
[K67.week_deseasonLAIfparDG,K67.week_trendLAIfparDG,K67.week_mtrendLAIfparDG] = AM_week2trend(K67.weekLAIfparDG,K67.week,1,'LAIfpar'); hold on;
[K67.week_deseasonrPAR,K67.week_trendrPAR,K67.week_mtrendrPAR] = AM_week2trend(K67.weekrPAR,K67.week,1,'rPAR'); hold on;
[K67.week_deseasonPAR,K67.week_trendPAR,K67.week_mtrendPAR] = AM_week2trend(K67.weekPAR,K67.week,1,'PAR'); hold on;

[K67.week_deseasonGmax,K67.week_trendGmax,K67.week_mtrendGmax] = AM_week2trend(K67.weekGmax,K67.week,1,'Gmax'); hold on;
[K67.week_deseasonALPHA,K67.week_trendALPHA,K67.week_mtrendALPHA] = AM_week2trend(K67.weekALPHA,K67.week,1,'ALPHA'); hold on;
[K67.week_deseasonGmaxAM,K67.week_trendGmaxAM,K67.week_mtrendGmaxAM] = AM_week2trend(K67.weekGmaxAM,K67.week,1,'GmaxAM'); hold on;
[K67.week_deseasonALPHAam,K67.week_trendALPHAam,K67.week_mtrendALPHAam] = AM_week2trend(K67.weekALPHAam,K67.week,1,'ALPHAam'); hold on;
% [K67.week_deseasonExtCoeff,K67.week_trendExtCoeff,K67.week_mtrendExtCoeff] = AM_week2trend(K67.weekExtCoeff,K67.week,1,'ExtCoeff'); hold on;

[K67.week_deseasonUstarday,K67.week_trendUstarday,K67.week_mtrendUstarday] = AM_week2trend(K67.weekUstarday,K67.week,1,'u_{*day}'); hold on;
[K67.week_deseasonBowenday,K67.week_trendBowenday,K67.week_mtrendBowenday] = AM_week2trend(K67.weekBowenday,K67.week,1,'Bowen_{daytime}'); hold on;

[K67.week_deseasonEvaporation,K67.week_trendEvaporation,K67.week_mtrendEvaporation] = AM_week2trend(K67.weekEvaporation,K67.week,1,'Evaporation'); hold on;
[K67.week_deseasonTranspiration,K67.week_trendTranspiration,K67.week_mtrendTranspiration] = AM_week2trend(K67.weekTranspiration,K67.week,1,'Transpiration'); hold on;

[K67.week_deseasonBelterraCInoon,K67.week_trendBelterraCInoon,K67.week_mtrendBelterraCInoon] = AM_week2trend(K67.weekBelterraCInoon,K67.week,1,'BelterraCInoon'); hold on;
[K67.week_deseasonBelterraCI,K67.week_trendBelterraCI,K67.week_mtrendBelterraCI] = AM_week2trend(K67.weekBelterraCI,K67.week,1,'BelterraCI'); hold on;