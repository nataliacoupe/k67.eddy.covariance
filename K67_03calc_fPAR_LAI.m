% Calculates FPAR and LAI as a functions of PAR observations
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
%% Leaf area
% .........................................................................
% Seca Floresta
LAIsf = load([FolderK67lai 'LeafAreaIndex.dat']);
dayLAIsf = datenum((LAIsf(:,1)),(LAIsf(:,2)),(LAIsf(:,3)));
dataLAIsfINI = LAIsf(:,5);
% dataLAIsfINI((dataLAIsfINI<4.4)|(dataLAIsfINI>6.9)) = NaN;

ix = dataLAIsfINI(dayLAIsf>=datenum(2002,1,1));     % &(~isnan(dataLAIsf)))
iz = nanmean(dataLAIsfINI(dayLAIsf<datenum(2002,1,1)));

iy =  detrend(ix)+iz;

dataLAIsf = [dataLAIsfINI(dayLAIsf<datenum(2002,1,1));iy];
dataLAIsf((dataLAIsf>6.6)|(dataLAIsf<5.4)) = NaN;

dailyLAIsf = (dayLAIsf(1):1:dayLAIsf(end))';
dataLAIsf  = interp1(dayLAIsf,dataLAIsf,dailyLAIsf,'linear');   %x(86:end) = NaN;

[K67.weekLAIsf,~,~,~,~,~,K67.weeksf] = AM_week2day_rs(dataLAIsf,dailyLAIsf);
[K67.weekYsf,K67.weekMsf,~] = datevec(K67.weeksf);
K67.weekHY_Ysf = K67.weekYsf;
K67.weekHY_Ysf(K67.weekMsf < month_start_hydro) = K67.weekHY_Ysf(K67.weekMsf < month_start_hydro)-1;
[K67.week_avgLAIsf,~,K67.week_avgstdLAIsf] = AM_week2_avg(K67.weekLAIsf,K67.weeksf);

%% .........................................................................
opts = delimitedTextImportOptions("NumVariables", 42);
opts.DataLines = [5, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["ht", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule   = "read";
LAIabshightSmith2019 = readtable([FolderK67laiLIDAR 'LAI_abs_hight_Smith_2019.csv'], opts);
LAIabshightSmith2019 = table2array(LAIabshightSmith2019);

% row 2 = year, 3: month, 4 = day
opts.DataLines = [2, 4];
LAIabshightSmith2019date = readtable([FolderK67laiLIDAR 'LAI_abs_hight_Smith_2019.csv'], opts);
LAIabshightSmith2019date = table2array(LAIabshightSmith2019date);

dayLAIlidar = datenum(LAIabshightSmith2019date(1,2:end)',LAIabshightSmith2019date(2,2:end)',LAIabshightSmith2019date(3,2:end)');
heightLAIlidar = LAIabshightSmith2019(:,1);
dataLAIlidar = NaN(length(dayLAIlidar),1);
for ik = 1:length(dayLAIlidar)
    ix = LAIabshightSmith2019(:,ik+1);
    dataLAIlidar(ik) = sum(ix(~isnan(ix)));
end
dailyLAIlidar = (dayLAIlidar(1):1:dayLAIlidar(end))';
x = interp1(dayLAIlidar,dataLAIlidar,dailyLAIlidar,'linear');   %x(86:end) = NaN;
% Available data: 2009, 22 Aug 2012-3 Dec 2012 and 3 Nov 2015-21 Jun 2017
x((dailyLAIlidar>datenum(2011,1,1))&(dailyLAIlidar<datenum(2012,8,21))) = NaN;
x((dailyLAIlidar>datenum(2012,12,4))&(dailyLAIlidar<datenum(2015,11,2))) = NaN;
x((dailyLAIlidar>datenum(2017,6,22))&(dailyLAIlidar<datenum(2020,1,1))) = NaN;
x2010_2013 = x;

x((dailyLAIlidar>datenum(2017,6,22))&(dailyLAIlidar<datenum(2020,1,1))) = NaN;
x2010_2013((dailyLAIlidar>datenum(2015,11,1))&(dailyLAIlidar<datenum(2018,1,1))) = NaN;

dataLAIlidar((dayLAIlidar>datenum(2011,1,1))&(dayLAIlidar<datenum(2012,8,21))) = NaN;
dataLAIlidar((dayLAIlidar>datenum(2012,12,4))&(dayLAIlidar<datenum(2015,11,2))) = NaN;
dataLAIlidar((dayLAIlidar>datenum(2017,6,22))&(dayLAIlidar<datenum(2020,1,1))) = NaN;

a = nanmin(dataLAIsf);          b = nanmax(dataLAIsf);
normalLAIsf = (dataLAIsf-nanmin(dataLAIsf))./(nanmax(dataLAIsf)-nanmin(dataLAIsf));
% normalLAIlidar = (x-nanmin(x))./(nanmax(x)-nanmin(x));
normalLAIlidar = (x-nanmin(x2010_2013))./(nanmax(x2010_2013)-nanmin(x2010_2013));

LAIlidar = a+((b-a).*normalLAIlidar);
LAIlidarSFmean = dataLAIlidar.*nanmean(dataLAIsf)./nanmean(dataLAIlidar);
% % a = nanmin(x);          b = nanmax(x);
% % LAIlidar = a+(((b-a).*(x-nanmin(dataLAIsf)))./(nanmax(dataLAIsf)-nanmin(dataLAIsf)));

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["VarName1", "VarName2"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
lbakm67LAIMCosta = readtable("/home/ncoupe/Documents/AmazonGEP/K67/LAI/lba_km67_LAI_MCosta.dat", opts);
dataLAIMCosta = lbakm67LAIMCosta.VarName2(1:end-1);
dayLAIMCosta = datenum(2009,lbakm67LAIMCosta.VarName1(1:end-1),1);

clear opts
%%
figure('color','white');
plot(dailyLAIlidar,LAIlidar,'-');         hold on;
plot(dailyLAIsf,dataLAIsf,'-');            plot(dayLAIMCosta,dataLAIMCosta,'color',[0.3 .8 .0]);  
plot(dayLAIlidar,dataLAIlidar,'.','color',[0.0 .5 .8],'MarkerSize',10);
plot(dayLAIlidar,LAIlidarSFmean,'.','color',[0.5 .8 0],'MarkerSize',10);
legend('Lidar scaled','SecaFloresta','Marcos Costa','Lidar original','Lidar x avg(SF)/avg(LiDAR)'); legend box off
datetick('x');        ylabel('LAI')
%%
figure;
plot(dailyLAIlidar,normalLAIlidar,'.');         hold on;
datetick('x');   ylabel('LAI');
% plot(dailyLAIlidar,LAIlidar,'.');         

%%
dailyLAI  = [dailyLAIsf;dailyLAIlidar];
LAI       = [dataLAIsf;LAIlidar];
normalLAI = [normalLAIsf;normalLAIlidar];

[K67.LAI_week,~,~,~,~,~,ix] = AM_week2day_rs(LAI,dailyLAI);
K67.weekLAI      = interp1(ix,K67.LAI_week,K67.week,'linear');
% K67.weekLAI = AM_fill_composite(K67.weekLAI,K67.week,23);
% x = smooth(K67.weekLAI);    x(isnan(K67.weekLAI)) = NaN;
% K67.weekLAI = x;
[K67.week_avgLAI,~,K67.week_avgstdLAI] = AM_week2_avg(K67.weekLAI,K67.week);

[K67.normalLAI_week,~,~,~,~,~,ix] = AM_week2day_rs(normalLAI,dailyLAI);
K67.weekLAInormal      = interp1(ix,K67.normalLAI_week,K67.week,'linear');
[K67.week_avgLAInormal,~,K67.week_avgstdLAInormal] = AM_week2_avg(K67.weekLAInormal,K67.week);

% .........................................................................
[K67.week_recurrentLAIsf]  = AM_week2dayAVG_recurrent(K67.week_avgLAIsf,K67.week);


% .........................................................................
v1 = 5.;         v2 = 6.6;  v3 = v1:((v2-v1)/4):v2;
figure('color','white');        hold on
plot(K67.week,K67.week_recurrentLAIsf,'color',[0.6 0 0]);
plot(K67.week,K67.weekLAI,'color',[0. 0 0.8], 'LineWidth',2);
datetick('x');        ylabel('LAI'); box on; grid on;
set(gca,'XTick',K67.week(1:23:end),'XTickLabel',datestr(K67.week(1:24:end),'yyyy'),'XLim',[K67.week(1) K67.week(end)],...
    'YTick',v3,'YLim',[v1 v2],'FontSize',11);


%% ........................................................................
K67.par1dn_sample = K67.par1dn;
ind = find((K67.jdstartgmt>=993.)&(K67.jdstartgmt<= 995.646));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1081.)&(K67.jdstartgmt<= 1084));        K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1412.438)&(K67.jdstartgmt<= 1416));     K67.par1dn_sample(ind) = NaN;
% ind = find((K67.jdstartgmt>=2100.)&(K67.jdstartgmt<= 2300.0));      K67.par1dn_sample(ind) = NaN; % no good quality

ind = find((K67.jdstartgmt>=3244.0)&(K67.jdstartgmt<= 3392.0));     K67.par1dn_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=3610.)&(K67.jdstartgmt<=3640));         K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3652.)&(K67.jdstartgmt<=3658.583));     K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3973.0)&(K67.jdstartgmt<=3985.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4185.0)&(K67.jdstartgmt<=4200.688));    K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4476)&(K67.jdstartgmt<=4480));          K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4268.688)&(K67.jdstartgmt<=4275.812));  K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4550.)&(K67.jdstartgmt<=4564.604));     K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4910.0)&(K67.jdstartgmt<=4980.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=5990.229)&(K67.jdstartgmt<=5998.271));  K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6108.750)&(K67.jdstartgmt<=6114.646));  K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6219.5)&(K67.jdstartgmt<=6220.5));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6261.0)&(K67.jdstartgmt<=6264.0));      K67.par1dn_sample(ind) = NaN;  %check
ind = find((K67.jdstartgmt>=6278.0)&(K67.jdstartgmt<=6285.0));      K67.par1dn_sample(ind) = NaN;  %check
ind = find((K67.jdstartgmt>=6354.0)&(K67.jdstartgmt<=6357.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6555.0)&(K67.jdstartgmt<=6564.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6612.0)&(K67.jdstartgmt<=6630.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6790.5)&(K67.jdstartgmt<=6794.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7082.0)&(K67.jdstartgmt<=7097.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7171.0)&(K67.jdstartgmt<=7175.5));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7487.0)&(K67.jdstartgmt<=7503.0));      K67.par1dn_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=1405.0)&(K67.jdstartgmt<=1407.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1405.0)&(K67.jdstartgmt<=1407.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6838.0)&(K67.jdstartgmt<=6860.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1458.0)&(K67.jdstartgmt<=1459.0));      K67.par1dn_sample(ind) = NaN;

K67.par1dn_sample(K67.par1dn_sample>150) = NaN;
K67.par1dn_sample((K67.So==tresholdSo)&(~isnan(K67.par1dn_sample))) = 0;
K67.par1dn_sample(K67.par1dn_sample<=0)  = NaN;
K67.par1dn_sample = AM_fill_plus(K67.par1dn_sample);
% ind = find((K67.jdstartgmt>=6000));         K67.par1dn_sample(ind) = K67.par1dn_sample(ind).*.9;

K67.par2up_sample = K67.par2up;

ind = find((K67.jdstartgmt>=2130.0)&(K67.jdstartgmt<=2215.0));      K67.par2up_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=3272.0)&(K67.jdstartgmt<=3273.0));      K67.par2up_sample(ind) = NaN;
% % ind = find((K67.jdstartgmt>=5994.0)&(K67.jdstartgmt<=6002.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7227.0)&(K67.jdstartgmt<=7233.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7541.0)&(K67.jdstartgmt<=7560.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7583.0)&(K67.jdstartgmt<=7586.0));      K67.par2up_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=1273.0)&(K67.jdstartgmt<=1275.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1384.0)&(K67.jdstartgmt<=1385.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1458.0)&(K67.jdstartgmt<=1459.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1875.0)&(K67.jdstartgmt<=1877.0));      K67.par2up_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=3615.0)&(K67.jdstartgmt<=3618.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3626.0)&(K67.jdstartgmt<=3635.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3652.0)&(K67.jdstartgmt<=3657.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3973.0)&(K67.jdstartgmt<=3985.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=5589.0)&(K67.jdstartgmt<=5590.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6169.0)&(K67.jdstartgmt<=6175.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6180.0)&(K67.jdstartgmt<=6181.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6193.0)&(K67.jdstartgmt<=6196.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6257.0)&(K67.jdstartgmt<=6260.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=6612.0)&(K67.jdstartgmt<=6615.0));      K67.par2up_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=6624.0)&(K67.jdstartgmt<=6629.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7164.0)&(K67.jdstartgmt<=7166.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7208.0)&(K67.jdstartgmt<=7213.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7216.0)&(K67.jdstartgmt<=7221.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7228.0)&(K67.jdstartgmt<=7238.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7270.0)&(K67.jdstartgmt<=7271.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7272.0)&(K67.jdstartgmt<=7273.0));      K67.par2up_sample(ind) = NaN;

% ind = find((K67.jdstartgmt>=3260.0)&(K67.jdstartgmt<= 3280.0));     K67.par2up_sample(ind) = NaN;


ind = find((K67.jdstartgmt>=6125.0)&(K67.jdstartgmt<= 6140.0));         %2006 bad PAR
ind = find((K67.jdstartgmt>=2100.0)&(K67.jdstartgmt<= 2250.0));         %2006 bad PAR
ind = find((K67.jdstartgmt>=5790.0)&(K67.jdstartgmt<= 5580.0));         %2006 bad PAR
K67.par1dn_sample(ind) = NaN;       K67.par2up_sample(ind) = NaN;

K67.par2up_sample((~isnan(K67.par2up))&(K67.So==0)) = 0;
K67.par2up_sample(K67.par2up>1400) = NaN;
K67.par2up_sample(K67.par2up_sample<0) = 0;
K67.par2up_sample = AM_fill_plus(K67.par2up_sample);
% ind = find((K67.jdstartgmt>=6000));         K67.par2up_sample(ind) = K67.par2up_sample(ind).*.9;

% % % K67.radiation = K67.par;
% % % [~,~,K67.radiation] = AM_compute_regression(K67.netrad,K67.radiation,[.4 .4 .4],0);
% % % K67.radiation((K67.radiation<=0)&(K67.So==0)) = 0;
% % % K67.radiation(K67.radiation<0) = NaN;


[K67.PAR2up_recurrent,K67.PAR2up_arecurrent] = AM_cycle_repeat(K67.par2up_sample,K67.dateloc);
[K67.PAR1dn_recurrent,K67.PAR1dn_arecurrent] = AM_cycle_repeat(K67.par1dn_sample,K67.dateloc);
[K67.radiation_recurrent,K67.radiation_arecurrent] = AM_cycle_repeat(K67.radiation,K67.dateloc);

%% ........................................................................
figure('color','white');
subplot(3,1,1);                     hold on;
plot(K67.jdstartgmt,K67.radiation./3,'color','r');
plot(K67.jdstartgmt,K67.PAR2up_arecurrent,'color','k');
plot(K67.jdstartgmt,K67.PAR2up_recurrent,'color',[0.0 0.3 0.7]);
plot(K67.jdstartgmt,K67.par2up,'linewidth',2,'color',[0.3 0.7 0]);
plot(K67.jdstartgmt,K67.par2up_sample,'color',[.9 .7 0.0]);
legend('radiation','par_{2up recurrent}','par_{year recurrent}','par_{2cpy}','par_{2cal clean}');

subplot(3,1,2);                     hold on;
plot(K67.jdstartgmt,K67.radiation./12,'color','m');
plot(K67.jdstartgmt,K67.PAR1dn_arecurrent,'color','k');
plot(K67.jdstartgmt,K67.PAR1dn_recurrent,'color',[0.0 0.3 0.7]);
plot(K67.jdstartgmt,K67.par1dn,'linewidth',2,'color',[0.3 0.7 0]);
plot(K67.jdstartgmt,K67.par1dn_sample,'color',[.9 .7 0.0]);
legend('radiation','par_{1dn recurrent}','par_{year recurrent}','par_{1dn}','par_{1dn clean}');

subplot(3,1,3);                     hold on;
plot(K67.jdstartgmt,K67.radiation_arecurrent,'color','k');
plot(K67.jdstartgmt,K67.radiation_recurrent,'color',[0.0 0.3 0.7]);
plot(K67.jdstartgmt,K67.radiation,'color',[.9 .7 0.0]);
plot(K67.jdstartgmt,K67.par,'color',[0.9 0.0 0]);
legend('par_{recurrent}','par_{year recurrent}','par_{radiation}','par');

%% ........................................................................
figure('color','white');
subplot(2,1,1);       hold on;
plot(K67.jdstartgmt,K67.par,'color',[0.0 0.3 0.7]);     plot(K67.jdstartgmt,K67.par2up_sample,'color',[0.3 0.7 0]);
plot(K67.jdstartgmt,K67.par1dn,'color',[.9 .7 0.0]);    plot(K67.jdstartgmt,K67.par1dn_sample.*20,'color',[.9 0 0]);
plot(K67.jdstartgmt,K67.par1dn,'color',[.0 .0 0.0]);
legend('par','par_{cpy}','par_{dwn incl cal}','par_{dwn}');

subplot(2,1,2);       hold on;
plot(K67.dateloc,K67.So./1.5,'color',[0.6 0.6 0.6]);
plot(K67.dateloc,K67.par,'color',[0.0 0.3 0.7]);         plot(K67.dateloc,K67.netrad.*2,'color',[0.8 0. 0.8]);
plot(K67.dateloc,K67.par2up_sample,'color',[0.3 0.7 0],'linewidth',2);
plot(K67.dateloc,K67.par1dn_sample.*20,'color',[.9 0 0]);
ind = find((K67.dateloc>=datenum(2010,5,2))&(K67.dateloc<=datenum(2011,7,1)));   %dates set to clean and calibrate
plot(K67.dateloc(ind),K67.par1dn_sample(ind).*20,'color',[.8 0.8 0.0]);
legend('TOA','par','netrad','par_{2up}','par_{dwn}');
datetick('x');

%% ........................................................................
% K67.fPAR = (K67.parfill - K67.par1dn_sample)./K67.par;
[K67.weekPAR1dn,~,~,~,~,K67.week_stdPAR1dn,K67.week] = AM_week2cycle_rs(K67.par1dn_sample,K67.dateloc,tresholdCycle,4,2);
[K67.weekPAR2up,~,~,~,~,K67.week_stdPAR2up] = AM_week2cycle_rs(K67.par2up_sample,K67.dateloc,tresholdCycle,4,2);
[K67.weekPAR,~,~,~,~,K67.week_stdPAR] = AM_week2cycle_rs(K67.radiation,K67.dateloc,tresholdCycle,4,2);
K67.weekjdstartgmt = AM_week2cycle_rs(K67.jdstartgmt,K67.dateloc,tresholdCycle,4,2);


K67.par1dn_fill = K67.par1dn_sample;
[p,r2,ix,~,sm,iz] = AM_compute_regression_zero(K67.par2up_sample,K67.par1dn_sample,[0.4,0.4,0.4],0);
K67.par1dn_fill(isnan(K67.par1dn_fill)) = iz(isnan(K67.par1dn_fill));
K67.par1dn_fill((K67.par1dn_fill<0)|(K67.So==0)) = 0;

K67.par1dn_sample(isnan(K67.par1dn_sample(ind))) = K67.par1dn_fill(isnan(K67.par1dn_sample(ind)));

K67.fPAR = (K67.radiation - K67.par1dn_fill - K67.par2up_sample)./K67.radiation;
ind = find((K67.So<=tresholdSo)|(K67.radiation<100)|(K67.par1dn_sample==0)|(K67.par2up_sample==0)|....
    (isnan(K67.par2up_sample))|(isnan(K67.radiation))|(isnan(K67.par1dn_sample)));
K67.fPAR(ind) = NaN;
K67.fPAR((K67.fPAR<0.)|(K67.fPAR>1)) = NaN;

K67.rPAR = K67.par1dn_fill./K67.radiation;
K67.rPAR(ind) = NaN;
K67.rPAR((K67.rPAR<0.)|(K67.rPAR>1)) = NaN;


K67.par1dn_day = K67.par1dn_sample;     K67.par1dn_day(K67.So<=tresholdSo) = NaN;
K67.par2up_day = K67.par2up_sample;     K67.par2up_day(K67.So<=tresholdSo) = NaN;
K67.par_day = K67.par;                  K67.par_day(K67.So<=tresholdSo) = NaN;

color_mtx = jet(12);            iy = ([2002:2005,2008:2012,2015:2020])';
figure('color','white');        count = 0;
for ik = 1:length(iy)
    count = 1 +count;
    ind = find(K67.Y==iy(ik));
    subplot(3,5,count);             hold on;
    [cycle_month,~,~,~,cycle] = AM_day_cycle(K67.par2up_sample(ind), K67.dateloc(ind),6,18);
    for ij = 1:12
        plot((1/24:1/24:1),cycle_month(ij,:),'color',color_mtx(ij,:));
        plot((1/24:1/24:1),cycle,'k','linewidth',2);
    end
    datetick('x','HH:MM');  xticks(0:1/4:1);
    xlim([0 1]);   ylim([0 500]);      %datetick('x','HR');
    title([num2str(iy(ik))]);
    if ((ik==1)|(ik==6)|(ik==12));  ylabel('PAR_{2up}(mumol m^-^2 s^-^1)'); end
end

%
color_mtx = jet(12);            iy = ([2002:2005,2008:2012,2015:2020])';
figure('color','white');        count = 0;
for ik = 1:length(iy)
    count = 1 +count;
    ind = find(K67.Y==iy(ik));
    subplot(3,5,count);             hold on;
    [cycle_month,~,~,~,cycle] = AM_day_cycle(K67.par1dn_sample(ind), K67.dateloc(ind),6,18);
    for ij = 1:12
        plot((1/24:1/24:1),cycle_month(ij,:),'color',color_mtx(ij,:));
        plot((1/24:1/24:1),cycle,'r','linewidth',2);
    end
    datetick('x','HH:MM');  xticks(0:1/4:1);
    xlim([0 1]);   ylim([0 80]);      %datetick('x','HR');
    title([num2str(iy(ik))]);
    if ((ik==1)|(ik==6)|(ik==12));  ylabel('PAR_{1dn}(mumol m^-^2 s^-^1)'); end
end

%% ........................................................................
figure('color','white');
subplot(3,1,1);       hold on;
plot(K67.jdstartgmt,K67.fPAR,'color',[0.0 0.3 0.7]);

subplot(3,1,2);       hold on;
plot(K67.week,K67.weekPAR./4,'color',[0.7 0.3 0.]);    plot(K67.week,K67.weekPAR1dn,'color',[0.0 0.3 0.7]);
plot(K67.week,K67.weekPAR2up,'color',[0.0 0.7 0.]);    datetick('x');
legend('PAR','PAR1dn','PAR2up')

subplot(3,1,3);       hold on;
plot(K67.weekjdstartgmt,K67.weekPAR./4,'color',[0.7 0.3 0.]);    plot(K67.weekjdstartgmt,K67.weekPAR1dn,'color',[0.0 0.3 0.7]);
plot(K67.weekjdstartgmt,K67.weekPAR2up,'color',[0.0 0.7 0.]);
legend('PAR','PAR1dn','PAR2up')

%%
figure('color','white');
v1 = .7;                v2 = .95;       v3 = v1:((v2-v1)/4):v2;
subplot(3,1,1);         hold on;
K67.weekfPAR = (K67.weekPAR - K67.weekPAR1dn - K67.weekPAR2up)./K67.weekPAR;
plot(K67.week,K67.weekfPAR)
[K67.weekfPAR,~,~,~,K67.weekfPARmax,K67.week_stdfPAR] = AM_week2cycle_rs(K67.fPAR,K67.dateloc,tresholdCycleDay,4,2);
plot(K67.week,K67.weekfPAR); hold on;   plot(K67.week,K67.weekfPARmax)
% [~,K67.weekfPAR] = AM_spike_detection(K67.weekfPAR,3);
% K67.weekfPAR(K67.weekfPAR<.7) = NaN;        K67.weekfPAR    = AM_fill_plus2(K67.weekfPAR);
iq = smooth(K67.weekfPAR,5);      iq(isnan(K67.weekfPAR)) = NaN;  K67.weekfPAR = iq;
plot(K67.week,K67.weekfPAR)
% K67.weekfPAR     = AM_rm_outlier(K67.weekfPAR,3);
% K67.weekfPAR  = OZ_out_bad_neighbor(K67.weekfPAR);   %K67.weekfPAR = OZ_out_bad_neighbor(K67.weekfPAR);
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);  grid on;
ylabel('fPAR');
legend('week','hourly','hourly maximum','moving window');

%%
v1 = .03;         v2 = .07;       v3 = v1:((v2-v1)/4):v2;
subplot(3,1,2)
[K67.weekrPAR,~,~,~,~,K67.week_stdrPAR] = AM_week2cycle_rs(K67.rPAR,K67.dateloc,tresholdCycleDay,4,2);
K67.weekrPAR(K67.weekrPAR<0.034) = NaN;
iq = smooth(K67.weekrPAR,5);      iq(isnan(K67.weekrPAR)) = NaN;  K67.weekrPAR = iq;
K67.weekrPAR = AM_fill_plus(K67.weekrPAR);
plot(K67.week,K67.weekrPAR);                hold on;
K67.weekrPAR = K67.weekPAR1dn./K67.weekRadiation;
plot(K67.week,K67.weekrPAR);                hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);  grid on;
ylabel('reflected PAR')

v1 = 5;         v2 = 7;       v3 = v1:((v2-v1)/4):v2;
subplot(3,1,3)
plot(K67.week,K67.weekLAI)
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);  grid on;
ylabel('LAI')

%% ........................................................................
% % K67.weekfPAR = 0.95.*(1-exp(-k.*K67.LAI);
% exp(-k.*LAI) = (1-(fPAR/0.95));   1/exp(k.*LAI) = (1-(fPAR/0.95))
ix = (1./(1-(K67.weekfPARmax./0.95)));
ix((ix<=0)|(ix>100)) = NaN;
iy = K67.weekLAI;                       %iy(K67.week>datenum(2008,1,1)) = NaN;
K67.weekExtCoeff = log(ix)./iy;         %K67.weekExtCoeff(K67.weekExtCoeff<0.22) = NaN;
K67.weekExtCoeff = AM_fill_plus(K67.weekExtCoeff);
% iq = smooth(K67.weekExtCoeff,5);        iq(isnan(K67.weekExtCoeff)) = NaN;  K67.weekrPAR = iq;

[K67.week_deseasonLAI,K67.week_trendLAI,K67.week_mtrendLAI] = AM_week2trend(K67.weekLAI,K67.week,1,'LAI'); hold on;
[K67.week_deseasonExtCoeff,K67.week_trendExtCoeff,K67.week_mtrendExtCoeff] = AM_week2trend(K67.weekExtCoeff,K67.week,1,'ExtCoeff'); hold on;
close; close;

% K67.extCoeff_recurrent = (K67.weekExtCoeff_recurrent.*0)+0.25;
K67.week_avgExtCoeff = AM_week2_avg(K67.weekExtCoeff,K67.week);
K67.weekExtCoeff_recurrent = AM_week2dayAVG_recurrent(K67.week_avgExtCoeff,K67.week);
% K67.weekExtCoeff_recurrent = AM_fill_month_cycle_force(K67.weekExtCoeff,K67.week,6,18);
% K67.weekExtCoeff_recurrent = K67.weekExtCoeff_recurrent+.05;
% K67.weekExtCoeff_recurrent = (K67.week.*0)+nanmean(K67.week_avgExtCoeff);

K67.week_timeseriesExtCoeff = K67.week_trendExtCoeff+K67.weekExtCoeff_recurrent-nanmean(K67.weekExtCoeff_recurrent);

% iq = (0.8:0.2/(length(K67.week)-1):1.0)';     %iq = 1-iq(1:end-1);
% % iq((K67.week>datenum(2008,1,1))&(K67.week<datenum(2015,1,1))) = 1.1;
% % iq(K67.week<datenum(2011,1,1)) = 1;
% K67.weekExtCoeff_recurrent = K67.weekExtCoeff_recurrent.*iq;

K67.weekLAIfpar = log(ix)./K67.week_timeseriesExtCoeff;
% [~,K67.weekLAIfpar] = AM_spike_detection(K67.weekLAIfpar,5);
iq = smooth(K67.weekLAIfpar,5);      iq(isnan(K67.weekLAIfpar)) = NaN;  K67.weekLAIfpar = iq;
K67.weekLAIfpar = AM_fill_plus2(K67.weekLAIfpar);
[K67.week_avgLAIfpar,~,K67.week_avgstdLAIfpar] = AM_week2_avg(K67.weekLAIfpar,K67.week);

% .........................................................................
%% Doughty and Goulden 2008
% .........................................................................
% % % % % ind = ((K67.sza>30)|(K67.par<1400)|(K67.So==0));       
% % % % % ix = (K67.par2up-200)./K67.par;         ix(ind) = NaN;     ix(ix<=0) = NaN;    scaleLAI = 1;
% % % % % ix = (K67.par1dn.*2)./K67.par;          ix(ind) = NaN;  ix(ix<=0) = NaN;    ix(ix>0.14) = NaN; scaleLAI = 0;
% % % % % ix = .01;                               ix(ind) = NaN;  ix(ix<=0) = NaN;    scaleLAI = 0;
% % % % % iy = 0.41./cosd(K67.sza);               iy(ind) = NaN;
% % % % % K67.LAIfparDG = -log(ix)./iy;
% % % % % K67.LAIfparDG((K67.LAIfparDG>7)|(K67.LAIfparDG<=4.5)) = NaN;
% % % % % K67.LAIfparDG = K67.LAIfparDG+scaleLAI;
% % % % % [~,~,~,~,~,~,K67.dayLAIfparDG] = 2(K67.LAIfparDG,K67.dateloc,0);


ind = ((K67.sza>30)|(K67.par<1400)|(K67.So==0));       
ix = K67.par1dn.*2./K67.par;             ix(ind) = NaN;  %ix(ix<=0) = NaN;    ix(ix>0.14) = NaN; scaleLAI = 0;
iy = 0.41./cosd(K67.sza);                %iy(ind) = NaN;
K67.LAIfparDG = -log(ix)./iy;
K67.LAIfparDG((K67.LAIfparDG>7)|(K67.LAIfparDG<=4)) = NaN;
K67.LAIfparDG = K67.LAIfparDG+scaleLAI;
[~,~,~,~,~,~,K67.dayLAIfparDG] = AM_day(K67.LAIfparDG,K67.dateloc,0);

K67.monthLAIfparDGfill = AM_month(K67.dayLAIfparDG,K67.day);

K67.monthLAIfparDGfill = AM_fill_plus2(K67.monthLAIfparDGfill);
K67.monthLAIfparDGfill = AM_fill_month_cycle(K67.monthLAIfparDGfill,K67.month,6,18);
K67.LAIfparDGfill = interp1(K67.month,K67.monthLAIfparDGfill,K67.dateloc);

K67.weekLAIfparDG = AM_week2day_rs(K67.dayLAIfparDG,K67.day);
% K67.weekLAIfparDG = K67.weekLAIfparDG+3;
% ix = smooth(K67.weekLAIfparDG,3);    ix(isnan(K67.weekLAIfparDG )) = NaN;     K67.weekLAIfparDG  = ix;
K67.weekLAIfparDG = AM_fill_plus2(K67.weekLAIfparDG);
[K67.week_avgLAIfparDG,~,K67.week_avgstdLAIfparDG] = AM_week2_avg(K67.weekLAIfparDG,K67.week);

display(nanmean(K67.weekLAIsf));
ind = find(K67.week<datenum(2008,1,1));     display(nanmean(K67.weekLAIfparDG(ind)));

K67.dayLAIMCosta = interp1(dayLAIMCosta,dataLAIMCosta,K67.day);
K67.weekLAIMCosta = AM_week2day_rs(K67.dayLAIMCosta,K67.day);

%%
close all
    color_mtx = turbo(5);
figure('color','white');    hold on;    
hl6 = plot(K67.week,K67.week_recurrentLAIsf,'color',[0.6 0.6 0.6]);
hl1 = plot(K67.dateloc,K67.LAIfparDG,'.','color',color_mtx(1,:));
hl2 = plot(K67.day,K67.dayLAIfparDG,'.','color',color_mtx(2,:));
hl3 = plot(K67.week,K67.weekLAIfparDG,'color',color_mtx(3,:),'LineWidth',2);
hold on
% plot(K67.week,K67.weekLAIfpar);
hl4 = plot(K67.week,K67.weekLAI,'color',color_mtx(4,:),'LineWidth',2);
ind = find(K67.week>datenum(2008,1,1));
hl5 = plot(K67.week(ind),K67.weekLAI(ind),'color','r','LineWidth',2);
% hl5 = plot(K67.week,K67.weekLAIMCosta,'yellow','LineWidth',2);
datetick('x');              grid on;                ylabel('LAI_{fPAR} (m^2 m^-^2)')
xtickangle(0)
legend([hl1,hl2,hl3,hl4,hl5,hl6],'LAI_{fnPAR hourly}','LAI_{fnPAR daily max}','LAI_{fnPAR}','LAI','LiDAR','seasonal average');
legend box off;             xlim([K67.week(1) K67.week(end)]);ylim([4.8 7]);
box on;
%% ........................................................................
[K67.week_deseasonLAIfparDG,K67.week_trendLAIfparDG,K67.week_mtrendLAIfparDG] = AM_week2trend(K67.weekLAIfparDG,K67.week,1,'LAIfpar'); hold on;
[K67.week_deseasonLAIfpar,K67.week_trendLAIfpar,K67.week_mtrendLAIfpar] = AM_week2trend(K67.weekLAIfpar,K67.week,1,'LAIfpar'); hold on;
[K67.week_deseasonPAR1dn,K67.week_trendPAR1dn,K67.week_mtrendPAR1dn]    = AM_week2trend(K67.weekPAR1dn,K67.week,1,'PAR1dn'); hold on;
[K67.week_deseasonPAR2up,K67.week_trendPAR2up,K67.week_mtrendPAR2up]    = AM_week2trend(K67.weekPAR2up,K67.week,1,'PAR2up'); hold on;
[K67.week_deseasonPAR,K67.week_trendPAR,K67.week_mtrendPAR] = AM_week2trend(K67.weekPAR,K67.week,1,'PAR'); hold on;
[K67.week_deseasonNetRadday,K67.week_trendNetRadday,K67.week_mtrendNetRadday] = AM_week2trend(K67.weekNetRadday,K67.week,1,'Rn');
[K67.week_deseasonCERESnr,K67.week_trendCERESnr,K67.week_mtrendCERESnr]       = AM_week2trend(K67.weekCERESnr,K67.week,1,'Rn');
[K67.week_deseasonfPAR,K67.week_trendfPAR,K67.week_mtrendfPAR]          = AM_week2trend(K67.weekfPAR,K67.week,1,'fPAR'); hold on;
K67.weekLAIfpar_recurrent = AM_week2dayAVG_recurrent(K67.weekLAIfpar,K67.week);
close; close; close; close; close; close; close; close;

%% ........................................................................
figure('color','white');
plot(K67.week,K67.weekLAIfpar,'k');         hold on;
plot(K67.week,K67.weekLAI,'b');             datetick('x');

%% ........................................................................
figure('color','white','Position',[50 50 600 400],'units','inches');
subplot(3,1,1)
v1 = 12;                        v2 = 22;                            v3 = v1:((v2-v1)/4):v2;
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);     iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);     set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);     set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl3=plot(K67.week,AM_fill_plus(K67.weekPAR1dn),'LineWidth',1,'color',[0.4 0.4 0.6 ]); hold on;
hl6=plot(K67.week,K67.week_deseasonPAR1dn,'LineWidth',2,'color',[0.6 0. 0.]); hold on;
hl4=plot(K67.week,(K67.week.*0)+nanmean(K67.weekPAR1dn),':','LineWidth',1,'color',[0.6 0.4 0.4 ]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
ylabel ({'PAR_{1dn}';'(\mumol m^-^2 s^-^1)'},'FontSize',12);

subplot(3,1,2)
v1 = 40;                        v2 = 100;                           v3 = v1:((v2-v1)/4):v2;
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);     iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);     set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);     set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl3=plot(K67.week,AM_fill_plus(K67.weekPAR2up),'LineWidth',1,'color',[0.4 0.4 0.6 ]); hold on;
hl6=plot(K67.week,K67.week_deseasonPAR2up,'LineWidth',2,'color',[0.6 0. 0.]); hold on;
hl4=plot(K67.week,(K67.week.*0)+nanmean(K67.weekPAR2up),':','LineWidth',1,'color',[0.6 0.4 0.4 ]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
ylabel ({'PAR_{2up}';'(\mumol m^-^2 s^-^1)'},'FontSize',12);

subplot(3,1,3)
v1 = 200;                       v2 = 600;                           v3 = v1:((v2-v1)/4):v2;
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);     iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);     set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);     set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl3=plot(K67.week,AM_fill_plus(K67.weekPAR),'LineWidth',1,'color',[0.4 0.4 0.6 ]); hold on;
hl6=plot(K67.week,K67.week_deseasonPAR,'LineWidth',2,'color',[0.6 0. 0.]); hold on;
hl4=plot(K67.week,(K67.week.*0)+nanmean(K67.weekPAR),':','LineWidth',1,'color',[0.6 0.4 0.4 ]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
ylabel ({'PAR';'(\mumol m^-^2 s^-^1)'},'FontSize',12);

%% ....................................................................
figure('color','white','Position',[50 50 600 400],'units','inches');
subplot(3,1,1)
v1 = 150;                         v2 = 350;                           v3 = v1:((v2-v1)/4):v2;
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);     iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);     set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);     set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl3=plot(K67.week,AM_fill_plus(K67.weekNetRadday),'LineWidth',1,'color',[0.4 0.4 0.6 ]); hold on;
hl6=plot(K67.week,K67.week_deseasonNetRadday,'LineWidth',2,'color',[0.6 0. 0.]); hold on;
hl4=plot(K67.week,(K67.week.*0)+nanmean(K67.weekNetRadday),':','LineWidth',1,'color',[0.6 0.4 0.4 ]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
ylabel ({'NetRad_{day}';'(W m^-^2)'},'FontSize',12);

subplot(3,1,2)
v1 = 125;                         v2 = 175;                           v3 = v1:((v2-v1)/4):v2;
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);     iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);     set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);     set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl3=plot(K67.week,AM_fill_plus(K67.weekCERESnr),'LineWidth',1,'color',[0.4 0.4 0.6 ]); hold on;
hl6=plot(K67.week,K67.week_deseasonCERESnr,'LineWidth',2,'color',[0.6 0. 0.]); hold on;
hl4=plot(K67.week,(K67.week.*0)+nanmean(K67.weekCERESnr),':','LineWidth',1,'color',[0.6 0.4 0.4 ]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
ylabel ({'NetRad_{CERES}';'(W m^-^2)'},'FontSize',12);

subplot(3,1,3)
v1 = 4;                         v2 = 7;                           v3 = v1:((v2-v1)/4):v2;
ix = find(K67.monthTRMM<100);   iy = ones(length(K67.month),1);     iy = iy.*v1;
iz = iy;                        iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);     set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);     set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl3=plot(K67.week,AM_fill_plus(K67.weekLAIfpar),'LineWidth',1,'color',[0.4 0.4 0.6 ]); hold on;
hl6=plot(K67.week,K67.week_deseasonLAIfpar,'LineWidth',2,'color',[0.6 0. 0.]); hold on;
hl5=plot(K67.week,K67.weekLAIfpar_recurrent,'LineWidth',1,'color',[0.6 0. 0.9]); hold on;
hl4=plot(K67.week,(K67.week.*0)+nanmean(K67.weekLAIfpar),':','LineWidth',1,'color',[0.6 0.4 0.4 ]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
ylabel ({'LAI_{fPAR}'},'FontSize',12);


%%
color_mtx = [0.1020 0.5882 0.2549; 215/255  25/255 28/255; 43/255 131/255 186/255; ...
    253/255 174/255 97/255; 0.2549 0.1020 0.5882; 81/255 161/255 74/255;0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

v1 = 4.8;         v2 = 6.8;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgLAIfparDG-K67.week_avgstdLAIfparDG,2.*K67.week_avgstdLAIfparDG];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekLAIfparDG(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgLAIfparDG([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekLAIfparDG(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
    hl7 = plot(week_avg,K67.weekLAIfparDG(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
    %     hl3 = plot(week_avg,K67.weekLAIfparDG(K67.HY_weekY == 2009),'Color',color_mtx(1,:),'LineWidth',2);
    hl4 = plot(week_avg,K67.weekLAIfparDG(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'LAI_{fPAR}'},'FontSize',11);       box on;        xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','09-10','Jul2015-Jun2016','2016-2017');