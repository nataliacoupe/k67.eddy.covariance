% Canopy conductance ecosystem following Weir and Saleska
% Natalia Restrepo-Coupe
%..........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;
% FolderK67litter = '/home/ncoupe/Documents/AmazonGEP/K67/Litter/';
FolderK67lai    = '/home/ncoupe/Documents/AmazonGEP/K67/SecaForesta/';
FolderK67laiLIDAR    = '/home/ncoupe/Documents/AmazonGEP/K67/LAI/';
% FolderK67dendrometry = '/home/ncoupe/Documents/AmazonGEP/K67/Dendrometry/';

month_name = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];
cycle_hour = (0:23)';


color_mtx = [0.1020 0.5882 0.2549; 215/255  25/255 28/255; 43/255 131/255 186/255; ...
    253/255 174/255 97/255; 0.2549 0.1020 0.5882; 171/255 221/255 164/255;...
    0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

[K67.weekSo]   = AM_week2cycle_rs(K67.So,K67.dateloc,tresholdCycle,3,2);

flag_MonthlyCycleGs = 1;

% .........................................................................
%% Leaf temperature
% .........................................................................
canopy   = canopy_BF(ip);
extCoeff = 4.39-(3.97.*exp(-0.258.*K67.LAIfparDGfill));
extCoeff = AM_fill_month_cycle(extCoeff,K67.dateloc,6,18);
zeta     = canopy./height;
Uzeta    = K67.ws;

Uh = Uzeta./exp((extCoeff.*zeta)-1);

        % Malhado et al 2009 Spatial trends in leaf size of Amazonian rainforest trees : mesophyll size class (20.25–182.25 cm2)
        % Leaf = sqrt((182.25-20.25)./(2*100));     %[m]
Leaf = sqrt((182.25-20.25)./2);   % cm
Leaf = Leaf./100;                 % cm to m

one = (Leaf./Uh);           two = extCoeff.*(1-zeta)./2;
rbH = (150./K67.LAIfparDGfill).*sqrt(one.*exp(two).*1);
rbH((isinf(rbH))|(rbH<0)|(K67.So<tresholdSo)) = NaN;  %|(rbH>50)

Sc  = 0.67;
Pr  = 0.71;
f   = 0.5;
rbV = (1./f).*rbH.*((Sc./Pr).^(2/3));            % nanmean == 31
CpJ = 1.013*10^-3.*10^6;  %J kg-1 C-1
K67.TL  = ((K67.H.*rbH)./(K67.rho_a.*CpJ))+K67.tair_fill;
% K67.TL = K67.TL+1;   %%%this makes PM and FG equal
% [K67.cycleTL,K67.cycle_stdTL,~,K67.acycle_stdTL,K67.acycleTL,cycle_hour] = AM_cycle_composite(K67.TL,K67.dateloc);

% Bolton formula for Saturation Vapor Pressure in Pa
%es = 611.2 * exp( (17.67*Tc)./(Tc + 243.5) );
% es_TL = 611.2.*exp((17.62.*TL)./(243.12+TL));
K67.es_TL = 611.2.*exp((17.67.*K67.TL)./(243.5+K67.TL));

% rsv = canopy resistance to water vapor (s m-1)
% la       leaf boundary layer resistance to water vapor (mol m-2 s-1)
% E = Fh2oDry [mmol m-2 s-1]= flux of transpired water vapor (mol m-2 s-1)
% Fh2o = 0-10 eddy flux of H2O (mol/m2/sec) [relative to dry air]
% R = molar gas constant (8.314472 J mol-1 K-1)
% ea == e_fill [~3 kPa]= vapor pressure in the canopy air (Pa)
% x= TairC;         x(So<tresholdSo) = NaN;
Rmol  = 8.314472;         % J mol-1 K-1
K67.E = K67.fh2o./1000;   % mmol m-2 s-1 to mol m-2 s-1
K67.E((K67.So<tresholdSo)|(isnan(K67.LE))|(K67.LE<0)) = NaN;
K67.E(indDry==1) = NaN;
K67.E((K67.E>0.02)|(K67.E<-0.001)) = NaN;
% [K67.cycleE,K67.cycle_stdE,~,K67.acycle_stdE,K67.acycleE,cycle_hour] = AM_cycle_composite(K67.E,K67.dateloc);

K67.ea = K67.e_fill.*1000;    % kPa to Pa
K67.ea(K67.So<tresholdSo) = NaN;
one    = ((K67.es_TL-K67.ea)./(Rmol.*(K67.tair_fill+273.15).*K67.E));   %one((one<0)|(one>1000)) = NaN;
one(indDry==1) = NaN;
K67.rsV  = one-rbV;
K67.rsV((K67.rsV<0)|(K67.So<tresholdSo)|(K67.rsV>5E4)) = NaN;                 %|(rsV>1050)
K67.rsV = AM_fill(K67.rsV);
%        Pa        /J mol-1 K-1 * C * mol m-2 s-1      -    s m-1
%        Pa        /Pa m3                 m-2  s-1     -    s m-1
%        s         /m                                  -    s m-1

% GsFGmmol: gsv: canopy stomatal conductance to water vapor -flux method
% P == Psurf = atmospheric pressure [Pa] = 101,325 Pa = 101.3 kPa
K67.P = K67.pambpa_fill;
K67.GsFGmmol = K67.P./(Rmol.*(K67.TL+273.15).*K67.rsV);
K67.GsFGmmol(indDry==1) = NaN;
K67.GsFGmmol((K67.GsFGmmol>20)|(K67.So<tresholdSo)|(K67.GsFGmmol<0))  = NaN;
K67.GsFGmmol = AM_rm_outlier(K67.GsFGmmol,3);
K67.GsFGmmol = K67.GsFGmmol.*1000;        % mol to mmol m-2 s-1

K67.GsFGmmol24 = K67.GsFGmmol;
K67.GsFGmmol24(indDry24==1) = NaN;
%      Pa/J mol-1 K-1  K         1 /s m-1
%      Pa/Pa m3 mol-1            m s-1
%      mol m-2 s-1

%% From Pearcy et all. Measurement of transpiration and leaf conductance
% P in kPa, tair in degC, Gs from cm/s to mol m-2 s-1
% GsFGmmol = GsFGmmol.*1000; %from mol m-2 s-1 to mmol m-2 s-1
% one = 44.6.*(273./(273+K67.tair_fill)).*(K67.P./(1000.*101.3));
% K67.GsFG = K67.GsFGmmol.*1000./one;                            % mol to mmol to m s-1

% Gas ideal gas law Gs [m s-1] x Pa[kPa]/R [m3 kPa K-1 mol] x Ta [k]
% K67.GsPMallmmol = K67.GsPMall.*(K67.pambpa_fill./1000)./(R.*(273+K67.tasp_fill));
K67.GsFG = (K67.GsFGmmol./1000).*(R.*(273+K67.tasp_fill))./(K67.pambpa_fill./1000);
K67.GsFG = K67.GsFG.*1000;               % from m s-1 to mm s-1

K67.GsFG((isinf(K67.GsFG))|(K67.So<tresholdSo)|((K67.es_TL-K67.ea)<1)) = NaN;
K67.GsFG(indDry==1) = NaN;
[~,K67.GsFG] = AM_spike_detection(K67.GsFG,5);

K67.GsFG24 = K67.GsFG;
K67.GsFG24(indDry24==1) = NaN;
% Typical values for irrigated trials are: 300-700 mmol m-2s-1;
% and for mildly water stressed trials are: 80-300 mmol m-2s-1

%..........................................................................
%% week
%..........................................................................
% % % K67.GsV_week     = AM_week2cycle_rs(K67.GsV,K67.dateloc,tresholdCycleDay,4,1);
% % % K67.GsV_week     = AM_rm_outlier(K67.GsV_week,3);
% % % % K67.GsV_week    = AM_fill_composite(K67.GsV_week,K67.week,23);
% % % [K67.weekGsFG,~,K67.week_stdGsFG]       = AM_week2_avg(K67.GsV_week,K67.week);
% % % K67.weekGsFG = AM_fill_plus(K67.weekGsFG);
% % %
% % % K67.monthGsFGmmol  = AM_month2cycle(K67.GsV,K67.dateloc,tresholdCycleDay,4,1);
% % % K67.monthGsFGmmol  = AM_rm_outlier(K67.monthGsFGmmol,3);
% % % % K67.monthGsFGmmol = AM_fill_composite(K67.monthGsFGmmol,K67.month,12);
% % % [K67.GsV_month_avg,~,K67.GsV_month_avg_std]       = AM_month_avg(K67.monthGsFGmmol,K67.month);
% % % K67.GsV_month_avg = AM_fill_plus(K67.GsV_month_avg);
% % %

[K67.cycleGsFGmmol,K67.cycle_stdGsFGmmol,~,K67.acycle_stdGsFGmmol,K67.acycleGsFGmmol] =....
    AM_cycle_composite(K67.GsFGmmol,K67.dateloc);
[K67.cycle_monthGsFGmmol,~,K67.cycle_month,K67.cycle_month_stdGsFGmmol] = ...
    AM_day_cycle_year(K67.GsFGmmol,K67.dateloc);
[K67.weekGsFGmmol,~,~,~,~,K67.week_stdGsFGmmol] =...
    AM_week2cycle_rs(K67.GsFGmmol,K67.dateloc,tresholdCycleDay,4,1);
K67.weekGsFGmmol(K67.weekGsFGmmol>1700) = NaN;
% K67.weekGsFGmmol     = AM_rm_outlier(K67.weekGsFGmmol,3);
iw = smooth(K67.weekGsFGmmol);   iw(isnan(K67.weekGsFGmmol)) = NaN;   K67.weekGsFGmmol = iw;
% K67.weekGsFGmmol   = AM_fill_composite(K67.weekGsFGmmol,K67.week,23);

[K67.cycleGsFGmmol24,K67.cycle_stdGsFGmmol24,~,K67.acycle_stdGsFGmmol24,K67.acycleGsFGmmol24] =....
    AM_cycle_composite(K67.GsFGmmol24,K67.dateloc);
[K67.cycle_monthGsFGmmol24,~,K67.cycle_month,K67.cycle_month_stdGsFGmmol24] = ...
    AM_day_cycle_year(K67.GsFGmmol24,K67.dateloc);
[K67.weekGsFGmmol24,~,~,~,~,K67.week_stdGsFGmmol24] =...
    AM_week2cycle_rs(K67.GsFGmmol24,K67.dateloc,tresholdCycleDay,4,1);
% K67.weekGsFGmmol24     = AM_rm_outlier(K67.weekGsFGmmol24,3);
iw = smooth(K67.weekGsFGmmol24);   iw(isnan(K67.weekGsFGmmol24)) = NaN;   K67.weekGsFGmmol24 = iw;

[K67.week_avgGsFGmmol,~,K67.week_avgstdGsFGmmol] = AM_week2_avg(K67.weekGsFGmmol,K67.week);
K67.week_avgGsFGmmol = AM_fill_plus(K67.week_avgGsFGmmol);
K67.monthGsFGmmol    = AM_month2cycle(K67.GsFGmmol,K67.dateloc,tresholdCycleDay,4,1);
% K67.monthGsFGmmol    = AM_rm_outlier(K67.monthGsFGmmol,3);
% K67.monthGsFGmmol     = AM_fill_composite(K67.monthGsFGmmol,K67.month,12);
[K67.month_avgGsFGmmol,~,K67.month_avgstdGsFGmmol] = AM_month_avg(K67.monthGsFGmmol,K67.month);
K67.month_avgGsFGmmol = AM_fill_plus(K67.month_avgGsFGmmol);

[K67.week_avgGsFGmmol24,~,K67.week_avgstdGsFGmmol24] = AM_week2_avg(K67.weekGsFGmmol24,K67.week);
K67.week_avgGsFGmmol24 = AM_fill_plus(K67.week_avgGsFGmmol24);
K67.monthGsFGmmol24    = AM_month2cycle(K67.GsFGmmol24,K67.dateloc,tresholdCycleDay,4,1);
K67.monthGsFGmmol24    = AM_rm_outlier(K67.monthGsFGmmol24,3);
% K67.monthGsFGmmol24     = AM_fill_composite(K67.monthGsFGmmol24,K67.month,12);
[K67.month_avgGsFGmmol24,~,K67.month_avgstdGsFGmmol24] = AM_month_avg(K67.monthGsFGmmol24,K67.month);
K67.month_avgGsFGmmol24 = AM_fill_plus(K67.month_avgGsFGmmol24);

%% ........................................................................
[K67.cycleGsFG,K67.cycle_stdGsFG,~,K67.acycle_stdGsFG,K67.acycleGsFG] =....
    AM_cycle_composite(K67.GsFG,K67.dateloc);
[K67.cycle_monthGsFG,~,K67.cycle_month,K67.cycle_month_stdGsFG] = ...
    AM_day_cycle_year(K67.GsFG,K67.dateloc);
[K67.weekGsFG,~,~,~,~,K67.week_stdGsFG] =...
    AM_week2cycle_rs(K67.GsFG,K67.dateloc,tresholdCycleDay,4,1);
K67.weekGsFG     = AM_rm_outlier(K67.weekGsFG,3);
iw = smooth(K67.weekGsFG,3);   iw(isnan(K67.weekGsFG)) = NaN;   K67.weekGsFG = iw;
% K67.weekGsFG     = AM_fill_composite(K67.weekGsFG,K67.week,23);
[K67.week_avgGsFG,~,K67.week_avgstdGsFG] = AM_week2_avg(K67.weekGsFG,K67.week);
K67.week_avgGsFG = AM_fill_plus(K67.week_avgGsFG);
K67.monthGsFG    = AM_month2cycle(K67.GsFG,K67.dateloc,tresholdCycleDay,4,1);
K67.monthGsFG    = AM_rm_outlier(K67.monthGsFG,3);
% K67.monthGsFG   = AM_fill_composite(K67.monthGsFG,K67.month,12);
[K67.month_avgGsFG,~,K67.month_avgstdGsFG]       = AM_month_avg(K67.monthGsFG,K67.month);
K67.month_avgGsFG = AM_fill_plus(K67.month_avgGsFG);

%% ........................................................................
figure('color','white');
subplot(2,2,1);   hold on;
y = [K67.week_avgGsFGmmol-K67.week_avgstdGsFGmmol,2.*K67.week_avgstdGsFGmmol];
h  =  area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0.0);
hold on;
y = [K67.week_avgGsPMmmol-K67.week_avgstdGsPMmmol,2.*K67.week_avgstdGsPMmmol];
h  =  area(week_avg,y);         set(gca,'Layer','top');
set(h(2),'FaceColor',[.7 .9 .7],'EdgeColor',[.7 .9 .7]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0.0);         alpha(.5);          hold on;
hl1 = plot(week_avg,K67.week_avgGsFGmmol,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
hl2 = plot(week_avg,K67.week_avgGsPMmmol,'Color',[0. 0.7 0.],'LineStyle','-','LineWidth',2);
hl3 = plot(week_avg,K67.week_avgGsFGmmol24,'Color',[.0 0 0.7],'LineStyle','-','LineWidth',2);
% hl3 = plot(week_avg,K67.week_avgGsPMall,'Color',[0. 0.5 0.],'LineStyle',':','LineWidth',2);
datetick('x','mmm');    grid on; box on;                legend box off;
legend([hl1,hl2],'Flux-gradient','PM','FontSize',8);    ylabel('Gs (m s^-^1)','FontSize',11);
ylim([200 1200]);

subplot(2,2,3:4);  hold on
plot(K67.week,K67.weekGsFGmmol,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
plot(K67.week,K67.weekGsPMmmol,'Color',[0. 0.7 0.],'LineStyle','-','LineWidth',2);
plot(K67.week,K67.weekGsFGmmol24,'Color',[.0 0 0.7],'LineStyle','-','LineWidth',2);
datetick('x','yyyy');        xlim([K67.week(1) K67.week(end)]);       grid on; box on;
ylabel('Gs (m s^-^1)','FontSize',11);   legend('Flux','Penman');

% subplot(2,2,3);  hold on
% plot(K67.week,K67.weekGsFGmmol,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
% plot(K67.week,K67.weekGsPMmmol,'Color',[0. 0.7 0.],'LineStyle','-','LineWidth',2);
% datetick('x','yyyy');    xlim([K67.week(1) K67.week(end)]); grid on; box on;
% ylabel('Gs (mol m^-^2 s^-^1)','FontSize',11);

subplot(2,2,2);  hold on
y = [K67.acycleGsFGmmol-K67.acycle_stdGsFGmmol,2.*K67.acycle_stdGsFGmmol];
h  =  area(cycle_hour,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);
set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0.33);
hold on;
y = [K67.acycleGsPMmmol-K67.acycle_stdGsPMmmol,2.*K67.acycle_stdGsPMmmol];
h  =  area(cycle_hour,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.7 .9 .7],'EdgeColor',[.7 .9 .7]);    set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0.33);        alpha(.5);                  hold on;
% hl1 = plot(cycle_hour,K67.acycle_Rn./300,'Color',[0 0.5 0],'LineStyle','-','LineWidth',1);
% hl2 = plot(cycle_hour,K67.acycle_So./600,'Color',[0 0.5 0.5],'LineStyle','-','LineWidth',2);
hl1 = plot(cycle_hour,K67.acycleGsFGmmol,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
hl2 = plot(cycle_hour,K67.acycleGsPMmmol,'Color',[0. 0.7 0],'LineStyle','-','LineWidth',2);
hl3 = plot(cycle_hour,K67.acycleGsFGmmol24,'Color',[.0 0 0.0],'LineStyle',':','LineWidth',2);
% hl3 = plot(cycle_hour,K67.acycleGsPMmmol12,'Color',[.0 0 0.7],'LineStyle','-','LineWidth',2);
grid on;                box on;            ylabel('Gs (mmol m^-^2 s^-^1)','FontSize',11);
% text (10,0.1,site_name_BF(ip));
box on;
set(gca,'XTick',(0:6:24),'XTickLabel',{'6:00','12:00','18:00'},...
    'YTick',(0:400:1600),'YAxisLocation','left','FontSize',10,'FontName','Arial',...
    'XLim',[0,24],'YLim',[0,1600]);
% legend([hl1,hl3,hl4],'Rn','R&S','PM'); legend box off; xlabel('UTC');    grid on;
legend([hl1,hl2,hl3],'Gs_{FG}',"Gs_{PM}'",'Gs_{PM}','FontSize',8); legend box off; xlabel('UTC','FontSize',11);
grid on;                xtickangle(0)

%% ........................................................................
figure('color','white');
subplot(2,2,1);   hold on;
y = [K67.week_avgGsFG-K67.week_avgstdGsFG,2.*K67.week_avgstdGsFG];
h  =  area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);
set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0.0);
hold on;
y = [K67.week_avgGsPM-K67.week_avgstdGsPM,2.*K67.week_avgstdGsPM];
h  =  area(week_avg,y);         set(gca,'Layer','top');
set(h(2),'FaceColor',[.7 .9 .7],'EdgeColor',[.7 .9 .7]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0.0);         alpha(.5);          hold on;
hl1 = plot(week_avg,K67.week_avgGsFG,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
hl2 = plot(week_avg,K67.week_avgGsPM,'Color',[0. 0.7 0.],'LineStyle','-','LineWidth',2);
hl3 = plot(week_avg,K67.week_avgGsPMall,'Color',[0. 0.5 0.],'LineStyle',':','LineWidth',2);
datetick('x','mmm');    grid on; box on;                legend box off;
% legend([hl1,hl2],'Flux-gradient','PM','FontSize',8);    
ylabel('Gs (m s^-^1)','FontSize',11);                   ylim([5 25]);

v1 = 0;                         v2 = 1500;       v3 = v1:((v2-v1)/3):v2;
subplot(2,2,3:4);                   hold on
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
plot(K67.week,K67.weekGsFGmmol,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
plot(K67.week,K67.weekGsPMmmol,'Color',[0. 0.7 0.],'LineStyle','-','LineWidth',2);
plot(K67.week,K67.weekGsPMallmmol,'Color',[0. 0.5 0],'LineStyle',':','LineWidth',2);
datetick('x','yyyy');        xlim([K67.week(1) K67.week(end)]);       grid on; box on;
ylabel('Gs (mmol m^-^2 s^-^1)','FontSize',11);
% legend('Flux','Penman');

subplot(2,2,2);  hold on
y = [K67.acycleGsFGmmol-K67.acycle_stdGsFGmmol,2.*K67.acycle_stdGsFGmmol];
% h  =  area(cycle_hour,[y(4:end,:);y(1:3,:)]);       set(gca,'Layer','top');
h  =  area(cycle_hour,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);
set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0.33);
hold on;
y = [K67.acycleGsPMmmol-K67.acycle_stdGsPMmmol,2.*K67.acycle_stdGsPMmmol];
% h  =  area(cycle_hour,[y(4:end,:);y(1:3,:)]);       set(gca,'Layer','top');
h  =  area(cycle_hour,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.7 .9 .7],'EdgeColor',[.7 .9 .7]);    set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0.33);        alpha(.5);                  hold on;
% hl3 = plot(cycle_hour,[K67.acycleGsFGmmol(4:end,:);K67.acycleGsFGmmol(1:3,:)],'Color',[0 0 0],'LineStyle','-','LineWidth',2);
% hl4 = plot(cycle_hour,[K67.acycleGsPMmmol(4:end,:);K67.acycleGsPMmmol(1:3,:)],'Color',[0. 0.7 0],'LineStyle','-','LineWidth',2);
% hl5 = plot(cycle_hour,[K67.acycleGsPMallmmol(4:end,:);K67.acycleGsPMallmmol(1:3,:)],'Color',[0. 0.5 0],'LineStyle',':','LineWidth',2);
hl3 = plot(cycle_hour,[K67.acycleGsFGmmol],'Color',[0 0 0],'LineStyle','-','LineWidth',2);
hl4 = plot(cycle_hour,[K67.acycleGsPMmmol],'Color',[0. 0.7 0],'LineStyle','-','LineWidth',2);
hl5 = plot(cycle_hour,[K67.acycleGsPMallmmol],'Color',[0. 0.5 0],'LineStyle',':','LineWidth',2);
legend([hl3,hl4,hl5],'Gs_{FG}',"Gs_{PM}'",'Gs_{PM}'); legend box off;
grid on;                box on;            ylabel('Gs (mmol m^-^2 s^-^1)','FontSize',11);
box on;
set(gca,'XTick',(0:6:18),'XTickLabel',{'0:00','6:00','12:00','18:00'},...
    'YTick',(0:400:1600),'YAxisLocation','left','FontSize',10,'FontName','Arial',...
    'XLim',[0,24],'YLim',[0,1600]);
grid on;                xtickangle(0)

% .........................................................................
%% Monthly cycles
% .........................................................................
if flag_MonthlyCycleGs == 1
    figure('color','white');
    for ik = 1:12
        subplot(4,3,ik);  hold on
        y = [K67.cycleGsFGmmol(ik,:)'-K67.cycle_stdGsFGmmol(ik,:)',2.*K67.cycle_stdGsFGmmol(ik,:)'];
        h  =  area(cycle_hour,y./1000);       set(gca,'Layer','top');
        set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
        set(h,'BaseValue',0.33);
        hold on;
        y = [K67.cycleGsPMmmol(ik,:)'-K67.cycle_stdGsPMmmol(ik,:)',2.*K67.cycle_stdGsPMmmol(ik,:)'];
        h  =  area(cycle_hour,y./1000);       set(gca,'Layer','top');
        set(h(2),'FaceColor',[.8 .9 .8],'EdgeColor',[.8 .9 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
        set(h,'BaseValue',0.33);
        alpha(.5)
        hold on;
        %         y = [K67.cycleGsPMmmol(ik,:)'-K67.cycle_stdGsPMmmol(ik,:)',2.*K67.cycle_stdGsPMmmol(ik,:)'];
        %         h  =  area(cycle_hour,y);       set(gca,'Layer','top');
        %         set(h(2),'FaceColor',[.9 .7 .7],'EdgeColor',[.9 .7 .7]);       set(h(1),'FaceColor','none','EdgeColor','none');
        %         set(h,'BaseValue',0.33);
        %         hold on;
        %         hl1 = plot(cycle_hour,K67.cycle_Rn(ik,:)'./1000,'Color',[0 0.5 0],'LineStyle','-','LineWidth',1);
        %         hl2 = plot(cycle_hour,K67.cycle_So(ik,:)'./2000,'Color',[0 0.5 0.5],'LineStyle','-','LineWidth',1);
        hl3 = plot(cycle_hour,K67.cycleGsFGmmol(ik,:)'./1000,'Color',[0 0 0],'LineStyle','-','LineWidth',2);
        hl4 = plot(cycle_hour,K67.cycleGsPMmmol(ik,:)'./1000,'Color',[0 .8 0],'LineStyle','-','LineWidth',2);
        %         hl4 = plot(cycle_hour,K67.cycleGsPMmmol(ik,:)','Color',[0.9 0 0],'LineStyle','-','LineWidth',2);
grid on;                xtickangle(0);        box on;
        set(gca,'XTick',(0:6:18),'XTickLabel',{'0:00','6:00','12:00','18:00'},...
            'YTick',(0:.500:1.500),'YAxisLocation','left','FontSize',11,'FontName','Arial',...
            'XLim',[0,24],'YLim',[0,1.500]);    grid on;
        title(month_name(ik,:),'FontSize',10,'FontName','Arial','FontWeight','normal');
        %         if ik==12;      legend([hl1,hl3,hl4],'Rn','flux-grad','PM');  legend box off; end
        if ((ik==1)|(ik==4)|(ik==7)|(ik==10));  ylabel('Gs (mol m^-^2 s^-^1)','FontSize',11); end
        if (ik>9);      xlabel('local time');                          end
        %         if (ik==1);     text (7,2.4,site_name_BF(ip));          end
    end
    
    %% ....................................................................
    [K67.yearcycle,K67.monthcycle,~,K67.HRcycle,K67.MMcycle,~] = datevec(K67.cycle_month);
    Ycycle_year = unique(K67.yearcycle);        HRcycle_year = sort(unique(K67.HRcycle));
    for im=1:length(Ycycle_year)
        figure('color','white');
        for in=1:12
            subplot(3,4,in);
            ik = find((K67.monthcycle==in)&(Ycycle_year(im)==K67.yearcycle));
            if ~isempty(ik)
                y = [K67.cycle_monthGsFGmmol(ik,:)-K67.cycle_month_stdGsFGmmol(ik,:),2.*K67.cycle_month_stdGsFGmmol(ik,:)];
                h  =  area(HRcycle_year',y);       set(gca,'Layer','top');
                set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
                set(h,'BaseValue',0.33);
                hold on;
                %                 y = [K67.cycle_monthGsPMmmol(ik,:)-K67.cycle_month_stdGsPMmmol(ik,:),2.*K67.cycle_month_stdGsPMmmol(ik,:)];
                %                 h  =  area(HRcycle_year,y);       set(gca,'Layer','top');
                %                 set(h(2),'FaceColor',[.9 .7 .7],'EdgeColor',[.9 .7 .7]);       set(h(1),'FaceColor','none','EdgeColor','none');
                %                 set(h,'BaseValue',0.33);
                %                 hold on;
                %                 hl4 = plot(HRcycle_year,K67.cycle_month_Rn(ik,:)'./1000,'Color',[0 0.5 0],'LineStyle','-','LineWidth',1);
                %                 hl2 = plot(HRcycle_year,K67.cycle_month_So(ik,:)'./2000,'Color',[0 0.5 0.5],'LineStyle','-','LineWidth',1);
                hl3 = plot(cycle_hour,K67.cycleGsFGmmol(in,:)','Color',[0 .9 0],'LineStyle','-','LineWidth',2);
                hl1=plot(HRcycle_year,K67.cycle_monthGsFGmmol(ik,:),'LineWidth',2,'color',[0.4 0.4 0.4]); hold on;
                %                 hl3=plot(HRcycle_year,K67.cycle_monthGsPMmmol(ik,:),'LineWidth',2,'color',[0.9 0.4 0.4]);
                title (month_name(in,:));
                axis ([4 18 0 1500]);
                set(gca,'XTick',cycle_hour(4:4:18),'YTick',0:500:1500);
                set(gca,'XTickLabel',{'4:00','8:00','12:00','16:00','20;00'});         grid on;
                if ((in==10)||(in==12));           xlabel ('UTM');      end
                if in==1;        title(strcat(num2str(Ycycle_year(im)),'    ',month_name(in,:)));  end
                if (in==1)|(in==5)|(in==9);        ylabel ('Gs (mol m^-^2 s^-^1)');      end
                %                 if in==12;       legend([hl1,hl3,hl4],'flux-grad','PM','Rn'); legend box off;   end
            end
        end
    end
    
    %% ....................................................................
    [K67.cycleTL,K67.cycle_stdTL,~,K67.acycle_stdTL,K67.acycleTL] =....
        AM_cycle_composite(K67.TL,K67.dateloc);
    [K67.cycle_monthTL,~,~,K67.cycle_month_stdTL] = ...
        AM_day_cycle_year(K67.TL,K67.dateloc);
    
    [K67.cycleE,K67.cycle_stdE,~,K67.acycle_stdE,K67.acycleE] =....
        AM_cycle_composite(K67.E,K67.dateloc);
    [K67.cycle_monthE,~,~,K67.cycle_month_stdE] = ...
        AM_day_cycle_year(K67.E,K67.dateloc);
    
    [K67.cycle_rsV,K67.cycle_stdrsV,~,K67.acycle_stdrsV,K67.acycle_rsV] =....
        AM_cycle_composite(K67.rsV,K67.dateloc);
    [K67.cycle_month_rsV,~,~,K67.cycle_month_stdrsV] = ...
        AM_day_cycle_year(K67.rsV,K67.dateloc);
    
    %% .........................................................................
    ix = K67.cycle_month_rsV;       iy = K67.cycle_month_stdrsV;
    iz = K67.cycle_rsV;
    for im=1:length(Ycycle_year)
        figure('color','white');
        for in=1:12
            subplot(3,4,in);
            ik = find((K67.monthcycle==in)&(Ycycle_year(im)==K67.yearcycle));
            if ~isempty(ik)
                y = [ix(ik,:)-iy(ik,:),2.*iy(ik,:)];
                h  =  area(HRcycle_year',y);       set(gca,'Layer','top');
                set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
                set(h,'BaseValue',0.33);
                hold on;
                hl4 = plot(HRcycle_year,K67.cycle_monthRn(ik,:)'./2,'Color',[0 0.5 0],'LineStyle','-','LineWidth',1);
                hl3 = plot(cycle_hour,iz(in,:)','Color',[0 .9 0],'LineStyle','-','LineWidth',2);
                hl1=plot(HRcycle_year,ix(ik,:),'LineWidth',2,'color',[0.4 0.4 0.4]); hold on;
                %                 hl3=plot(HRcycle_year,K67.cycle_monthGsPMmmol(ik,:),'LineWidth',2,'color',[0.9 0.4 0.4]);
                title (month_name(in,:));
                axis ([0 24 0 200]);
                set(gca,'XTick',cycle_hour(1:6:24),'YTick',0:40:200);
                set(gca,'XTickLabel',{'0:00','6:00','12:00','18:00','0:00'});         grid on;
                if ((in==10)||(in==12));           xlabel ('UTM');      end
                if in==1;        title(strcat(num2str(Ycycle_year(im)),'    ',month_name(in,:)));  end
                if (in==1)|(in==5)|(in==9);        ylabel ('rsV (m^2 s mol^-^1)');      end
                %                 if in==12;       legend([hl1,hl3,hl4],'flux-grad','PM','Rn'); legend box off;   end
            end
        end
    end
end

%% ....................................................................
v1 = 0;         v2 = 1200;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgGsFGmmol-K67.week_avgstdGsFGmmol,2.*K67.week_avgstdGsFGmmol];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekGsFGmmol(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgGsFGmmol([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekGsFGmmol(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekGsFGmmol(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekGsFGmmol(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekGsFGmmol(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'Gs';'(mmol m^-^2 s^-^1)'},'FontSize',11);     box on;         xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');
% % %
% % % %% ....................................................................
% % % v1 = 0;         v2 = 3.0;  v3 = v1:((v2-v1)/2):v2;
% % % figure('color','white');
% % % subplot(3,1,1);             hold on;
% % % x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
% % % hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
% % % hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
% % % y = [K67.week_avgLeafLitter-K67.week_avgstdLeafLitter,2.*K67.week_avgstdLeafLitter];
% % % y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
% % % h = area(week_avg,y);       set(gca,'Layer','top');
% % % set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
% % % set(h,'BaseValue',0);
% % % ik = 2002;
% % % while ik<2019
% % %     hl2  =  plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
% % %     % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
% % %     ik = ik+1;
% % % end
% % % hl1 = plot(week_avg,K67.week_avgLeafLitter([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
% % % %     hl6 =K67.weekLeafLitter plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
% % % hl7 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
% % % %     hl3 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
% % % hl4 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
% % % % hl5 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2014),'Color',color_mtx(3,:),'LineWidth',2);
% % % xlim([week_avg(1) 365]);    ylim([v1 v2]);
% % % set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
% % % ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
% % % ylabel ({'NPP_{litter_{leaf}}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
% % % text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
% % % %     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
% % % legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');
% % %
% % % %% ....................................................................
% % % v1 = 0.6;         v2 = 3.4;  v3 = v1:((v2-v1)/2):v2;
% % % figure('color','white');
% % % subplot(3,1,1);             hold on;
% % % x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
% % % hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
% % % hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
% % % y = [K67.week_avgAllLitter-K67.week_avgstdAllLitter,2.*K67.week_avgstdAllLitter];
% % % y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
% % % h = area(week_avg,y);       set(gca,'Layer','top');
% % % set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
% % % set(h,'BaseValue',0);
% % % ik = 2002;
% % % while ik<2019
% % %     hl2  =  plot(week_avg,K67.weekAllLitter(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
% % %     % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
% % %     ik = ik+1;
% % % end
% % % hl1 = plot(week_avg,K67.week_avgAllLitter([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
% % % %     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
% % % hl7 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
% % % %     hl3 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
% % % hl4 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
% % % %     hl5 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2014),'Color',color_mtx(3,:),'LineWidth',2);
% % % xlim([week_avg(1) 365]);    ylim([v1 v2]);
% % % set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
% % % ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
% % % ylabel ({'NPP_{litter_{all}}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
% % % text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
% % % %     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
% % % legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');
% % %
% % % %%
% % % v1 = 0;         v2 = 3.0;  v3 = v1:((v2-v1)/2):v2;
% % % figure('color','white');
% % % subplot(3,1,1);             hold on;
% % % x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
% % % hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
% % % hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
% % % y = [K67.week_avgAllLitter-K67.week_avgstdAllLitter,2.*K67.week_avgstdAllLitter];
% % % y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
% % % h = area(week_avg,y);       set(gca,'Layer','top');
% % % set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
% % % set(h,'BaseValue',0);
% % % ik = 2002;
% % % while ik<2019
% % %     hl2  =  plot(week_avg,K67.weekAllLitter(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
% % %     ik = ik+1;
% % % end
% % % hl1 = plot(week_avg,K67.week_avgAllLitter([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
% % % %     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
% % % hl7 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
% % % %     hl3 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
% % % hl4 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
% % % %     hl5 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
% % % xlim([week_avg(1) 365]);    ylim([v1 v2]);
% % % set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
% % % ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
% % % ylabel ({'NPP_{litter}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
% % % % text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
% % % %     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
% % % % legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

%%
v1 = 5;         v2 = 7;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgLAI-K67.week_avgstdLAI,2.*K67.week_avgstdLAI];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekLAI(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgLAI([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekLAI(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekLAI(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekLAI(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekLAI(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'LAI'},'FontSize',11);       box on;           xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

%%
v1 = 0;         v2 = 1;  v3 = v1:((v2-v1)/2):v2;        v2 = 1.1;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgLAInormal-K67.week_avgstdLAInormal,2.*K67.week_avgstdLAInormal];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekLAInormal(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgLAInormal([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekLAInormal(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekLAInormal(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekLAInormal(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekLAInormal(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'LAI_{normalized}'},'FontSize',11);            box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;     xtickangle(0);
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

%% ....................................................................
x1 = 100;       x2 = 500;          	x3 = x1:((x2-x1)/3):x2;
y1 = 0;         y2 = 1200;           y3 = y1:((y2-y1)/3):y2;
figure('color','white','position',[100,100,400,1200]);
subplot(4,1,1);                     hold on;
plot(K67.weekNetRadday,K67.weekGsFGmmol,'.','color',[.5 .5 .5],'MarkerSize',10);
[p,r2,~,~,sm] = AM_compute_regression(K67.weekNetRadday,K67.weekGsFGmmol,[0.4,0.4,0.4],0);
y_model = polyval(p,K67.weekNetRadday);                     hold on;
plot(K67.weekNetRadday,y_model,'-','color',[0.4,0.4,0.4]);
text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10);                            box on;
ik = [2015;2008];       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_weekY == ik(ij));
    hl2  =  plot(K67.weekNetRadday(ind),K67.weekGsFGmmol(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(K67.weekNetRadday(ind),K67.weekGsFGmmol(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,K67.weekNetRadday(ind));                     hold on;
    plot(K67.weekNetRadday(ind),y_model,'-','Color',color_mtx(count,:));
    text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
        ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'XTickLabel',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ylabel ({'Gs (mmol m^-^2 s^-^1)'},'FontSize',12);   xlabel ({'Rn_{day} (W m^-^2)'},'FontSize',12);
text (x3(2),y3(3),site_name,'FontSize',12);        box on;

x1 = 0;             x2 = 600;                    x3 = x1:((x2-x1)/3):x2;
subplot(4,1,2);     hold on;
plot(K67.monthCWD,K67.monthGsFGmmol,'.','color',[.5 .5 .5],'MarkerSize',10);
ind = (K67.monthCWD>0);
[p,r2,~,~,sm] = AM_compute_regression(K67.monthCWD(ind),K67.monthGsFGmmol(ind),[0.4,0.4,0.4],0);
y_model = polyval(p,K67.monthCWD);                     hold on;
plot(K67.monthCWD,y_model,'-','color',[0.4,0.4,0.4]);
text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10);                            box on;
ik = [2015;2008];       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(K67.monthCWD(ind),K67.monthGsFGmmol(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(K67.monthCWD(ind),K67.monthGsFGmmol(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,K67.monthCWD(ind));                     hold on;
    plot(K67.monthCWD(ind),y_model,'-','Color',color_mtx(count,:));
    text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
        ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ylabel ({'Gs (mmol m^-^2 s^-^1)'},'FontSize',12);   xlabel ({'CWD (mm)'},'FontSize',12);

% ........................................................................
x1 = 160;           x2 = 280;           x3 = x1:((x2-x1)/4):x2;
subplot(4,1,3);     hold on;
plot(K67.monthCERESswd,K67.monthGsFGmmol,'.','color',[.5 .5 .5],'MarkerSize',10);
[p,r2,~,~,sm] = AM_compute_regression(K67.monthCERESswd,K67.monthGsFGmmol,[0.4,0.4,0.4],0);
y_model = polyval(p,K67.monthCERESswd);                     hold on;
plot(K67.monthCERESswd,y_model,'-','color',[0.4,0.4,0.4]);
text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10);                            box on;
ik = [2015;2008];       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_monthY == ik(ij));
    hl2  =  plot(K67.monthCERESswd(ind),K67.monthGsFGmmol(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(K67.monthCERESswd(ind),K67.monthGsFGmmol(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,K67.monthCERESswd(ind));                     hold on;
    plot(K67.monthCERESswd(ind),y_model,'-','Color',color_mtx(count,:));
    text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
        ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ylabel ({'Gs (mmol m^-^2 s^-^1)'},'FontSize',12);   xlabel ({'SW_{down CERES} (W m^-^2)'},'FontSize',12);

x1 = 0;             x2 = 2;           x3 = x1:((x2-x1)/4):x2;
% [K67.cycleTL,K67.cycle_stdTL,~,K67.acycle_stdTL,K67.acycleTL,cycle_hour] = AM_cycle_composite(K67.TL,K67.dateloc);
subplot(4,1,4);     hold on;
plot(K67.weekVPDday,K67.weekGsFGmmol,'.','color',[.5 .5 .5],'MarkerSize',10);
[p,r2,~,~,sm] = AM_compute_regression(K67.weekVPDday,K67.weekGsFGmmol,[0.4,0.4,0.4],0);
y_model = polyval(p,K67.weekVPDday);                     hold on;
plot(K67.weekVPDday,y_model,'-','color',[0.4,0.4,0.4]);
text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10);                            box on;
ik = [2015;2008];       count = 0;
for ij = 1:length(ik)
    count(ik(ij)==2015) = 2;            count(ik(ij)==2008) = 5;
    ind = find(K67.HY_weekY == ik(ij));
    hl2  =  plot(K67.weekVPDday(ind),K67.weekGsFGmmol(ind),'.','Color',color_mtx(count,:),'MarkerSize',10); hold on;
    [p,r2,~,~,sm] = AM_compute_regression(K67.weekVPDday(ind),K67.weekGsFGmmol(ind),[0.4,0.4,0.4],0);
    y_model = polyval(p,K67.weekVPDday(ind));                     hold on;
    plot(K67.weekVPDday(ind),y_model,'-','Color',color_mtx(count,:));
    text (x1+10,y1+20,{['Gs=' num2str(p(1),'%4.2f') 'Rn_{day}+' num2str(p(2),'%4.2f')];...
        ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color',color_mtx(count,:));
end
set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],'YTick',y3,'Ylim',[y1 y2],'FontSize',12);
box on;             grid on;
ylabel ({'Gs (mmol m^-^2 s^-^1)'},'FontSize',12);   xlabel ({'VPD_{day} (kPa)'},'FontSize',12);

%..........................................................................
%% save 
%..........................................................................
tableK67 = table(K67.Y,K67.M,K67.D,K67.HR,K67.GsFGmmol,K67.GsPMmmol);
tableK67.Properties.VariableNames = {'YYYY' 'MM' 'DD' 'HR'  'GsFG_molm2s1' 'GsPM_molm2s1'};
writetable (tableK67,'out/K67.Gs.dat');
