%% Calculates NEE and GEP quickly
% Natalia Restrepo-Coupe
% Toronto
% .........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;

month_name = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];

umol2gC  = 60.*60.*24.*12.*(10.^-6);
cycle_hour = (0:23)';

K67.gepDry    = -K67.gee;           K67.gepDry(indDry==1)            = NaN;
K67.gepDryDay = K67.gepDry;         K67.gepDryDay(K67.So<tresholdSo) = NaN;
K67.gepday    = -K67.gee;           K67.gepday(K67.So<tresholdSo)    = NaN;

K67.rDry   = K67.r;                 K67.rDry(indDry==1)   = NaN;
K67.neeDry = K67.neefill;           K67.neeDry(indDry==1) = NaN;

%% ........................................................................
color_mtx = jet(8);
figure('color','white');
subplot(3,1,1:3);       hold on;	plot(K67.dateloc,K67.co2pf1,'color',color_mtx(1,:));
plot(K67.dateloc,K67.co2pf3,'color',color_mtx(3,:));
plot(K67.dateloc,K67.co2pf4,'color',color_mtx(4,:));
plot(K67.dateloc,K67.co2pf5,'color',color_mtx(5,:));
plot(K67.dateloc,K67.co2pf6,'color',color_mtx(6,:));
plot(K67.dateloc,K67.co2pf7,'color',color_mtx(7,:));
plot(K67.dateloc,K67.co2pf8,'color',color_mtx(8,:));
plot(K67.dateloc,K67.co2pf2,'color',[.0 .0 .9]);
plot(K67.dateloc,K67.co2,'color','k');%[.7 .7 .7]);
datetick('x');                      ylabel('CO_2 (ppm)');

figure('color','white');
subplot(3,1,1:3);       hold on;    plot(K67.dateloc,K67.h2opf1);
plot(K67.dateloc,K67.h2opf3);       plot(K67.dateloc,K67.h2opf4);
plot(K67.dateloc,K67.h2opf5);       plot(K67.dateloc,K67.h2opf6);
plot(K67.dateloc,K67.h2opf7);       plot(K67.dateloc,K67.h2opf8);
plot(K67.dateloc,K67.h2opf2,'color',[.0 .0 .9]);
plot(K67.dateloc,K67.h2ommolm,'color','k');%[.7 .7 .7]);
datetick('x');                      ylabel('H_2O (ppt)');

%% ........................................................................
[K67.daynee,~,ix] = AM_day(K67.neefill,K67.dateloc,0);      K67.daynee(ix<22) = NaN;
K67.daynee        = AM_rm_outlier(K67.daynee,3);            K67.daynee = AM_fill_plus(K67.daynee);
[K67.dayneelight,~,ix] = AM_day(K67.neefilllight,K67.dateloc,0);K67.dayneelight(ix<22) = NaN;
K67.dayneelight        = AM_rm_outlier(K67.dayneelight,3);  K67.dayneelight = AM_fill_plus(K67.dayneelight);

[K67.weeknee,~,~,~,~,K67.week_stdnee] = AM_week2cycle_rs(K67.neefill,K67.dateloc,tresholdCycle,4,2);
% K67.weeknee   = AM_rm_outlier(K67.weeknee,3); 
K67.weeknee      = AM_fill(K67.weeknee);

[K67.dayr,~,ix] = AM_day(K67.r,K67.dateloc,0);      K67.dayr(ix<20) = NaN;
% K67.dayr        = AM_rm_outlier(K67.dayr,3);            
K67.dayr = AM_fill_plus(K67.dayr);

[K67.weekr,~,~,~,~,K67.week_stdr]     = AM_week2cycle_rs(K67.r,K67.dateloc,tresholdCycle,4,2);
K67.weekr   = AM_rm_outlier(K67.weekr,3);                   K67.weekr      = AM_fill(K67.weekr);
[K67.monthr,~,~,~,~,K67.month_stdr]   = AM_month2cycle(K67.r,K67.dateloc,tresholdCycle,4,2);
K67.monthr   = AM_rm_outlier(K67.monthr,3);                 K67.monthr      = AM_fill(K67.monthr);
[K67.weekrlight,~,~,~,~,K67.week_stdrlight] = AM_week2cycle_rs(K67.rlight,K67.dateloc,tresholdCycle,4,2);
K67.weekrlight   = AM_rm_outlier(K67.weekrlight,3);         K67.weekrlight = AM_fill(K67.weekrlight);

[K67.daygee,~,~,~,~,K67.day_stdgee] = AM_day(K67.gee,K67.dateloc,0);
K67.daygee   = AM_rm_outlier(K67.daygee,3);               K67.daygee    = AM_fill(K67.daygee);
[K67.daygee,~,~,~,~,K67.day_stdgee] = AM_day(K67.gee,K67.dateloc,0);
K67.daygee   = AM_rm_outlier(K67.daygee,3);               K67.daygee    = AM_fill(K67.daygee);

[K67.eightgee,~,~,~,~,K67.eight_stdgee] = AM_eight2cycle_rs(K67.gee,K67.dateloc,tresholdCycle,4,2);
K67.eightgee   = AM_rm_outlier(K67.eightgee,3);               K67.eightgee    = AM_fill(K67.eightgee);
[K67.weekgee,~,~,~,~,K67.week_stdgee] = AM_week2cycle_rs(K67.gee,K67.dateloc,tresholdCycle,4,2);
K67.weekgee   = AM_rm_outlier(K67.weekgee,3);               K67.weekgee    = AM_fill(K67.weekgee);
[K67.monthgee,~,~,~,~,K67.month_stdgee] = AM_month2cycle(K67.gee,K67.dateloc,tresholdCycle,4,2);
K67.monthgee   = AM_rm_outlier(K67.monthgee,3);             K67.monthgee    = AM_fill(K67.monthgee);
[K67.weekgeelight,~,~,~,~,K67.week_stdgeelight] = AM_week2cycle_rs(K67.geefilllight,K67.dateloc,tresholdCycle,4,2);
K67.weekgeelight   = AM_rm_outlier(K67.weekgeelight,3);     K67.weekgeelight = AM_fill(K67.weekgeelight);
[K67.monthnee,~,~,~,~,K67.month_stdnee] = AM_month2cycle(K67.neefill,K67.dateloc,tresholdCycle,4,2);
K67.monthnee = AM_rm_outlier(K67.monthnee,3);           K67.monthnee = AM_fill(K67.monthnee);

K67.weeknee = umol2gC.*(K67.weeknee);                   K67.week_stdnee = umol2gC.*(K67.week_stdnee);
K67.weekgee = umol2gC.*(K67.weekgee);                   K67.week_stdgee = umol2gC.*(K67.week_stdgee);
K67.monthgee = umol2gC.*(K67.monthgee);                 K67.month_stdgee = umol2gC.*(K67.month_stdgee);
K67.weekr   = umol2gC.*(K67.weekr);                     K67.week_stdr = umol2gC.*(K67.week_stdr);
K67.monthr   = umol2gC.*(K67.monthr);                   K67.month_stdr = umol2gC.*(K67.month_stdr);
K67.monthnee = umol2gC.*(K67.monthnee);                 K67.month_stdnee = umol2gC.*(K67.month_stdnee);

[K67.weekneelight,~,~,~,~,K67.week_stdneelight] = AM_week2cycle_rs(K67.neefilllight,K67.dateloc,tresholdCycle,4,2);
% K67.weekneelight   = AM_rm_outlier(K67.weekneelight,3); 
K67.weekneelight = AM_fill(K67.weekneelight);

K67.weekneelight = umol2gC.*(K67.weekneelight);         K67.week_stdneelight = umol2gC.*(K67.week_stdneelight);
K67.weekgeelight = umol2gC.*(K67.weekgeelight);         K67.week_stdgeelight = umol2gC.*(K67.week_stdgeelight);
K67.weekrlight   = umol2gC.*(K67.weekrlight);           K67.week_stdrlight = umol2gC.*(K67.week_stdrlight);

[K67.week_avggee,~,K67.week_avgstdgee] = AM_week2_avg(K67.weekgee,K67.week);
[K67.week_avggeelight,~,K67.week_avgstdneelight] = AM_week2_avg(K67.weekgeelight,K67.week);
[K67.week_avgnee,~,K67.week_avgstdnee] = AM_week2_avg(K67.weeknee,K67.week);
[K67.week_avgneelight,~,K67.week_avgstdneelight] = AM_week2_avg(K67.weekneelight,K67.week);
[K67.week_avgr,~,K67.week_avgstdr]     = AM_week2_avg(K67.weekr,K67.week);
[K67.week_avgrlight,~,K67.week_avgstdrlight]     = AM_week2_avg(K67.weekrlight,K67.week);
[K67.month_avgnee,K67.month_avgstdnee] = AM_month_avg(K67.monthnee,K67.month);

% .........................................................................
[K67.weekgepDry,~,~,~,~,K67.week_stdgeeDry]  = AM_week2cycle_rs(K67.gepDry,K67.dateloc,tresholdCycle,4,2);
K67.weekgepDry   = AM_rm_outlier(K67.weekgepDry,3);           K67.weekgepDry = AM_fill(K67.weekgepDry);
[K67.weekgepDryDay,~,~,~,~,K67.week_stdgepDryDay] = AM_week2cycle_rs(K67.gepDryDay,K67.dateloc,tresholdCycleDay,3,2);
K67.weekgepDryDay   = AM_rm_outlier(K67.weekgepDryDay,3);     K67.weekgepDryDay = AM_fill(K67.weekgepDryDay);
[K67.weekrDry,~,~,~,~,K67.week_stdrDry]      = AM_week2cycle_rs(K67.rDry,K67.dateloc,tresholdCycle,4,2);
K67.weekrDry   = AM_rm_outlier(K67.weekrDry,3);               K67.weekrDry   = AM_fill(K67.weekrDry);
[K67.weekneeDry,~,~,~,~,K67.week_stdneeDry]  = AM_week2cycle_rs(K67.neeDry,K67.dateloc,tresholdCycle,4,2);
K67.weekneeDry   = AM_rm_outlier(K67.weekneeDry,3);           K67.weekneeDry = AM_fill(K67.weekneeDry);

[K67.week_avggepDry,~,K67.week_avgstdgepDry] = AM_week2_avg(K67.weekgepDry,K67.week);
[K67.week_avggepDryDay,~,K67.week_avgstdneelight] = AM_week2_avg(K67.weekgepDryDay,K67.week);
[K67.week_avgneeDry,~,K67.week_avgstdneeDry] = AM_week2_avg(K67.weekneeDry,K67.week);
[K67.week_avgrDry,~,K67.week_avgstdrDry] = AM_week2_avg(K67.weekrDry,K67.week);

K67.weekgepDry = umol2gC.*(K67.weekgepDry);       	K67.week_stdgeeDry = umol2gC.*(K67.week_stdgeeDry);
K67.weekgepDryDay = umol2gC.*(K67.weekgepDryDay);   K67.week_stdgepDryDay = umol2gC.*(K67.week_stdgepDryDay);
K67.weekneeDry = umol2gC.*(K67.weekneeDry);       	K67.week_stdneeDry = umol2gC.*(K67.week_stdneeDry);
K67.weekrDry = umol2gC.*(K67.weekrDry);             K67.week_stdrDry = umol2gC.*(K67.week_stdrDry);

[K67.weekgepday,~,~,~,~,K67.week_stdgepday] = AM_week2cycle_rs(K67.gepday,K67.dateloc,tresholdCycleDay,3,2);
K67.weekgepday   = AM_rm_outlier(K67.weekgepday,3);               K67.weekgepday    = AM_fill(K67.weekgepday);
[K67.monthgepday,~,~,~,~,K67.month_stdgepday] = AM_month2cycle(K67.gepday,K67.dateloc,tresholdCycleDay,3,2);
K67.monthgepday   = AM_rm_outlier(K67.monthgepday,3);             K67.monthgepday    = AM_fill(K67.monthgepday);
[K67.week_avggepday,~,K67.week_avgstdgepday] = AM_week2_avg(K67.weekgepday,K67.week);

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
K67.radiation = K67.par;
K67.radiation((K67.dateloc>datenum(2012,10,31))&(K67.dateloc<datenum(2015,1,25))) = NaN;
[p,r2,K67.radiation,~,sm,~,iz] = AM_compute_regression(K67.netrad,K67.radiation,[0.4,0.4,0.4],0);
K67.radiation(isnan(K67.radiation)) = iz(isnan(K67.radiation));
K67.radiation((K67.radiation<0)|(K67.So==0)) = 0;
K67.weekRadiation = AM_week2cycle_rs(K67.radiation,K67.dateloc,tresholdCycle,3,2);
% ix = K67.LE+K67.H;
% K67.radiation(isnan(K67.radiation)) = ix(isnan(K67.radiation));

figure('color','white');    plot(K67.dateloc,K67.neefill);  datetick('x');  hold on;

[K67.GEPmodel_ust_min,K67.GEPmodel_ust_max,K67.GEP_ec,K67.Remodel_ust_min,K67.Remodel_ust_max,...
    K67.Re_ec,K67.NEE_ec_ust_min,K67.NEE_ec_ust_max,K67.NEE_ec] = ...
    AM_gep_ust(K67.nee,K67.radiation,K67.ustar,ust_min,ust_tresh,ust_max,...
    lat,lon,K67.dateloc+datenum(0,0,0,1,0,0),0,-3,-3,30);    %utm

% % K67.GEP_ec = umol2gC.*(K67.GEP_ec);
% % K67.GEPmodel_ust_min = umol2gC.*(K67.GEPmodel_ust_min);	 K67.GEPmodel_ust_max = umol2gC.*(K67.GEPmodel_ust_max);
% % K67.Re_ec  = umol2gC.*(K67.Re_ec);
% % K67.Remodel_ust_min  = umol2gC.*(K67.Remodel_ust_min);   K67.Remodel_ust_max  = umol2gC.*(K67.Remodel_ust_max);
% % K67.NEE_ec = umol2gC.*(K67.NEE_ec);
% % K67.NEE_ec_ust_min   = umol2gC.*(K67.NEE_ec_ust_min);    K67.NEE_ec_ust_max   = umol2gC.*(K67.NEE_ec_ust_max);
%% ........................................................................
K67.GEPmodel_ust_min(K67.GEPmodel_ust_min>60) = NaN;
K67.GEPmodel_ust_max(K67.GEPmodel_ust_max>60) = NaN;
K67.GEP_ec(K67.GEP_ec>60) = NaN;

[K67.weekGEP,~,~,~,~,K67.week_stdGEP] = AM_week2cycle_rs(K67.GEP_ec,K67.dateloc,tresholdCycle,4,2);
K67.weekGEP = AM_fill_composite(K67.weekGEP,K67.week,23);
[K67.weekRe,~,~,~,~,K67.week_stdRe]   = AM_week2cycle_rs(K67.Re_ec,K67.dateloc,tresholdCycle,4,2);
K67.weekRe  = AM_fill_composite(K67.weekRe,K67.week,23);

[K67.weekNEE,~,~,~,~,K67.week_stdNEE] = AM_week2cycle_rs(K67.NEE_ec,K67.dateloc,tresholdCycle,4,2);
[~,K67.weekNEE] = AM_spike_detection(K67.weekNEE,5);
K67.weekNEE = AM_fill_composite(K67.weekNEE,K67.week,23);

K67.weekNEE = umol2gC.*(K67.weekNEE);
K67.weekGEP = umol2gC.*(K67.weekGEP);
K67.weekRe  = umol2gC.*(K67.weekRe);

% % % QA ......................................................................
% % ind = find((K67.week>datenum(2015,1,1))&(K67.week<datenum(2015,5,16)));
% % K67.weekNEE(ind) = NaN;       K67.weekRe(ind) = NaN;    K67.weekGEP(ind) = NaN;
% % % .........................................................................

[K67.week_avgNEE,~,K67.week_avgstdNEE] = AM_week2_avg(K67.weekNEE,K67.week);
[K67.week_avgGEP,~,K67.week_avgstdGEP] = AM_week2_avg(K67.weekGEP,K67.week);
[K67.week_avgRe,~,K67.week_avgstdRe]   = AM_week2_avg(K67.weekRe,K67.week);

K67.monthGEP = AM_month2cycle(K67.GEP_ec,K67.dateloc,tresholdCycle,4,2);    K67.monthGEP = AM_fill(K67.monthGEP);
K67.monthNEE = AM_month2cycle(K67.NEE_ec,K67.dateloc,tresholdCycle,4,2);    K67.monthNEE = AM_rm_outlier(K67.monthNEE,3); 
K67.monthNEE = AM_fill(K67.monthNEE);
K67.monthRe  = AM_month2cycle(K67.Re_ec,K67.dateloc,tresholdCycle,4,2);     K67.monthRe  = AM_fill(K67.monthRe);
% [~,K67.mont/hNEE] = AM_spike_detection(K67.monthNEE,5);                      K67.monthNEE = AM_fill(K67.monthNEE);
%%
K67.weekWUEdryday = K67.weekgepDryDay.*umol2gC./K67.weekETdryday;
K67.weekWUEdryday(K67.weekWUEdryday>6) = NaN;
K67.weekWUEdryday = AM_fill(K67.weekWUEdryday);
K67.week_stdWUEdryday = K67.week_stdgepDryDay.*umol2gC./K67.week_stdETdryday;
[K67.week_avgWUE,~,K67.week_avgstdWUE] = AM_week2_avg(K67.weekWUEdryday,K67.week);

%..........................................................................
%% Pc
%..........................................................................
ind = find((K67.dateloc>datenum(2011,1,1))&(K67.dateloc<datenum(2012,1,1)));
figure('color','white');
plot(K67.dateloc(ind),K67.netrad(ind),'k');        hold on;     plot(K67.dateloc(ind),K67.netrad(ind),'.k');
plot(K67.dateloc(ind),K67.netradcal(ind),'g');     plot(K67.dateloc(ind),K67.radiation(ind)./4,'y');
plot(K67.dateloc(ind),K67.So(ind)./4,'m');
plot(K67.dateloc(ind),K67.H(ind),'r');             plot(K67.dateloc(ind),K67.LE(ind),'b'); hold on;
datetick('x')
%%
[K67.cycle_month_So] = AM_day_cycle_year(K67.So,K67.dateloc);
[iq,~,K67.cycle_month,it] = AM_day_cycle_year(K67.VPD,K67.dateloc);
[K67.yearcycle,K67.monthcycle,~,K67.HRcycle,K67.MMcycle,~] = datevec(K67.cycle_month);
Ycycle_year = unique(K67.yearcycle);        HRcycle_year = sort(unique(K67.HRcycle));
y1= 0;                     y2= 3;             y3 = y1:((y2-y1)/4):y2;
for im=1:length(Ycycle_year)
    figure('color','white');
    for in=1:12
        subplot(3,4,in);
        ik = find((K67.monthcycle==in)&(Ycycle_year(im)==K67.yearcycle));
        if ~isempty(ik)
            y = [iq(ik,:)-it(ik,:),2.*it(ik,:)];
            h  =  area(HRcycle_year',y);       set(gca,'Layer','top');
            set(h(2),'FaceColor',[.8 .8 .9],'EdgeColor',[.8 .8 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
            set(h,'BaseValue',0.33);
            hold on;
            hl2 = plot(HRcycle_year,K67.cycle_month_So(ik,:)'./500,'Color',[0 0.5 0.5],'LineStyle','-','LineWidth',1);
            hl1=plot(HRcycle_year,iq(ik,:),'LineWidth',2,'color',[0.4 0.4 0.4]); hold on;
            title (month_name(in,:));                           axis ([6 18 y1 y2]);
            set(gca,'XTick',cycle_hour(7:6:19),'YTick',y3);
            set(gca,'XTickLabel',{'6:00','12:00','18:00'});     grid on;
            if ((in==10)||(in==12));           xlabel ('UTM');  end
            if in==1;        title(strcat(num2str(Ycycle_year(im)),'    ',month_name(in,:)));  end
            if (in==1)|(in==5)|(in==9);        ylabel ('radiation');      end
            %             if in==12;       legend([hl1,hl3,hl4],'flux-grad','PM','Rn'); legend box off;   end
        end
    end
end
%% ........................................................................
% ind = find((K67.So>tresholdSo)&(K67.radiation>5));
[K67.weekPc,weekPc_count,weekPc_std,PARdayAVG,VPDdayAVG,TAdayAVG,CIdayAVG,...
    PARdaySTD,VPDdaySTD,TAdaySTD,CIdaySTD] = ...
    AM_week2day_Pc_TaPARvpdCI(K67.gepday,K67.radiation,K67.VPD_fill,...
    K67.tc_fill,K67.dateloc,lat,lon,-1);

K67.weekPc = umol2gC.*(K67.weekPc);
K67.weekPc = AM_fill_plus(K67.weekPc);

K67.monthPc = interp1(K67.week,K67.weekPc,K67.month);

[K67.week_avgPc,~,K67.week_avgstdPc] = AM_week2_avg(K67.weekPc,K67.week);

%..........................................................................
%% Select years
%..........................................................................
for ik = 1:length(K67.year)
    weekK67(ik).NEE = K67.weekNEE(K67.Yweek == K67.year(ik));
    weekK67(ik).GEP = K67.weekGEP(K67.Yweek == K67.year(ik));
    weekK67(ik).Re  = K67.weekRe(K67.Yweek  == K67.year(ik));
    
    weekK67hy(ik).NEE = K67.weekNEE(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).GEP = K67.weekGEP(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).Re  = K67.weekRe(K67.HY_weekY  == K67.year(ik));
    
    weekK67(ik).nee = K67.weeknee(K67.Yweek == K67.year(ik));
    weekK67(ik).gee = K67.weekgee(K67.Yweek == K67.year(ik));
    weekK67(ik).r  = K67.weekr(K67.Yweek  == K67.year(ik));
    weekK67(ik).Pc  = K67.weekPc(K67.Yweek  == K67.year(ik));
    
    weekK67hy(ik).nee = K67.weeknee(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).gee = K67.weekgee(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).r  = K67.weekr(K67.HY_weekY  == K67.year(ik));
    weekK67hy(ik).Pc  = K67.weekPc(K67.HY_weekY  == K67.year(ik));
    
    weekK67(ik).neeDry = K67.weekneeDry(K67.Yweek == K67.year(ik));
    weekK67(ik).gepDry = K67.weekgepDry(K67.Yweek == K67.year(ik));
    weekK67(ik).gepDryDay = K67.weekgepDryDay(K67.Yweek == K67.year(ik));
    weekK67(ik).rDry  = K67.weekrDry(K67.Yweek  == K67.year(ik));
    
    weekK67hy(ik).neeDry = K67.weekneeDry(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).gepDryDay = K67.weekgepDryDay(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).gepDry = K67.weekgepDry(K67.HY_weekY == K67.year(ik));
    weekK67hy(ik).rDry  = K67.weekrDry(K67.HY_weekY  == K67.year(ik));
    
end

%..........................................................................
%% Repeating annual cycle
%..........................................................................
site_LE_mean_month = [];     site_H_mean_month = [];
site_LE_std_month  = [];     site_H_std_month  = [];
for ik = 1:length(unique(K67.Y))
    site_LE_mean_month = [site_LE_mean_month;K67.week_avgLE];
    site_LE_std_month  = [site_LE_std_month;K67.week_avgstdLE];
    site_H_mean_month  = [site_H_mean_month;K67.week_avgH];
    site_H_std_month   = [site_H_std_month;K67.week_avgstdH];
end

%..........................................................................
%% Repeating annual cycle
%..........................................................................
figure('color','white');    subplot(2,1,1);
plot(K67.week,-umol2gC.*K67.weekGEP); datetick('x');  hold on;
plot(K67.week,K67.weekgee);
% plot(K67.week,-umol2gC.*K67.weekgee);
ylabel('GEP'); legend ('umolC0_2 m^-^2 s^-^1','gC0_2 m^-^2 d^-^1')

subplot(2,1,2);             plot(K67.week,K67.weeknee);      datetick('x');  hold on;
plot(K67.week,K67.weekNEE);
ylabel('NEE (umolC0_2 m^-^2 s^-^1)');legend ('','gC0_2 m^-^2 d^-^1')

GEPtable = [K67.monthY,K67.monthM,K67.monthD,K67.monthGEP,K67.monthRe,K67.monthNEE,K67.monthPc];
save ('K67.cflux.170913.dat','GEPtable','-ascii');

%% ........................................................................
figure('color','white');
subplot(2,1,1);                             plot(K67.dateloc,K67.neefill);
hold on;                                    plot(K67am.dateloc,K67am.NEEturb);
plot(K67am.dateloc,K67am.NEEfill,'k');      plot(K67.dateloc,K67.NEE_ec);
datetick('x'); legend('k67_{fill}','Ameriflux','Ameriflux_{fill}','k67 mine');
ylabel('NEE (\mumol m^-^2 s^-^1)');

subplot(2,1,2); plot(K67.week,K67.weeknee); hold on;    plot(K67am.week,K67am.weekNEEturb);
plot(K67am.week,K67am.weekNEEfill,'k');                 plot(K67.week,K67.weekNEE);
plot(K67h.week,K67h.weekNEE,'m');                       plot(K67.week,K67.weekneelight,'g');
datetick('x');	legend('k67','Ameriflux','Ameriflux_{fill}','K67_{mine}','Harvard','K67 light');
ylabel('NEE (\mumol m^-^2 s^-^1)');

%% ........................................................................
% [K67.GEPam_hyperbola,K67.GEPpm_hyperbola,K67.GEP_hyperbola,K67.GmaxAM,K67.GmaxPM,K67.Gmax,...
%     K67.ALPHAam,K67.ALPHApm,K67.ALPHA]  = ...
%     AM_gep_fourier_Gmax(K67.nee,K67.radiation,K67.ustd,ust_tresh,lat,lon,K67.dateloc,-0.5,-1,-1);
% 
%%
[~,~,~,K67.GmaxAM,K67.GmaxPM,K67.Gmax,K67.ALPHAam,K67.ALPHApm,K67.ALPHA]  = ...
    AM_Gmax_Quantum(K67.gee,K67.radiation,K67.ustd,ust_tresh,lat,lon,K67.dateloc,-1.,-3,-3);

%% .........................................................................
ind = find((K67.dateloc>datenum(2011,6,1))&(K67.dateloc<datenum(2013,12,31)));
K67.Gmax(ind) = NaN;       K67.GmaxAM(ind) = NaN;    K67.ALPHAam(ind) = NaN;    K67.ALPHA(ind) = NaN;

ind = find((K67.dateloc>datenum(2005,12,1))&(K67.dateloc<datenum(2008,1,1)));
K67.Gmax(ind) = NaN;       K67.GmaxAM(ind) = NaN;    K67.ALPHAam(ind) = NaN;    K67.ALPHA(ind) = NaN;

K67.ALPHA(K67.ALPHA>0.12) = NaN;        K67.ALPHAam(K67.ALPHAam>0.12) = NaN;

[K67.weekGmaxIni,~,~,~,~,K67.week_stdGmax]  = AM_week2cycle_rs(K67.Gmax,K67.dateloc,tresholdCycleDay,3,2);
K67.weekGmax = K67.weekGmaxIni;
% K67.weekGmax = AM_rm_outlier(K67.weekGmax,3);
ix = smooth(K67.weekGmax);      ix(isnan(K67.weekGmax)) = NaN;    K67.weekGmax = ix;   
K67.weekGmax = AM_fill(K67.weekGmax);
[K67.week_avgGmax,~,K67.week_avgstdGmax] = AM_week2_avg(K67.weekGmax,K67.week);

[K67.weekALPHAini,~,~,~,~,K67.week_stdALPHA] = AM_week2cycle_rs(K67.ALPHA,K67.dateloc,tresholdCycleDay,3,2);
K67.weekALPHA = K67.weekALPHAini;
K67.weekALPHA = K67.weekALPHA.*10^6./(60.*60.*24.*2);   K67.week_stdALPHA = K67.week_stdALPHA.*10^6./(60.*60.*24.*2);
% K67.weekALPHA(K67.weekALPHA>0.12) = NaN;
% K67.weekALPHA = AM_rm_outlier(K67.weekALPHA,3);
ix = smooth(K67.weekALPHA);      ix(isnan(K67.weekALPHA)) = NaN;      K67.weekALPHA = ix;     
K67.weekALPHA = AM_fill(K67.weekALPHA);

% .........................................................................
[K67.weekGmaxAMini,~,~,~,~,K67.week_stdGmaxAM]  = AM_week2cycle_rs(K67.GmaxAM,K67.dateloc,tresholdCycleDay,3,2);
K67.weekGmaxAM = K67.weekGmaxAMini;
% K67.weekGmaxAM = AM_rm_outlier(K67.weekGmaxAM,3);
ix = smooth(K67.weekGmaxAM);      ix(isnan(K67.weekGmaxAM)) = NaN;    K67.weekGmaxAM = ix;   
K67.weekGmaxAM = AM_fill(K67.weekGmaxAM);
[K67.week_avgGmaxAM,~,K67.week_avgstdGmaxAM] = AM_week2_avg(K67.weekGmaxAM,K67.week);

[K67.weekALPHAamini,~,~,~,~,K67.week_stdALPHAam] = AM_week2cycle_rs(K67.ALPHAam,K67.dateloc,tresholdCycleDay,3,2);
K67.weekALPHAam = K67.weekALPHAamini;
K67.weekALPHAam = K67.weekALPHAam.*10^6./(60.*60.*24.*2);   K67.week_stdALPHAam = K67.week_stdALPHAam.*10^6./(60.*60.*24.*2);
% K67.weekALPHAam(K67.weekALPHAam>0.12) = NaN;
% K67.weekALPHAam = AM_rm_outlier(K67.weekALPHAam,3);
ix = smooth(K67.weekALPHAam);      ix(isnan(K67.weekALPHAam)) = NaN;      K67.weekALPHAam = ix;     
K67.weekALPHAam = AM_fill(K67.weekALPHAam);


[K67.monthALPHAamini,~,~,~,~,K67.month_stdALPHAam] = AM_month2cycle(K67.ALPHAam,K67.dateloc,tresholdCycleDay,3,2);
K67.monthALPHAam = K67.monthALPHAamini;
K67.monthALPHAam = K67.monthALPHAam.*10^6./(60.*60.*24.*2);   K67.month_stdALPHAam = K67.month_stdALPHAam.*10^6./(60.*60.*24.*2);

% .........................................................................
[K67.week_avgALPHAam,~,K67.week_avgstdALPHAam] = AM_week2_avg(K67.weekALPHAam,K67.week);
[K67.week_avgGmaxAM,~,K67.week_avgstdGmaxAM]   = AM_week2_avg(K67.weekGmaxAM,K67.week);
[K67.week_avgALPHA,~,K67.week_avgstdALPHA] = AM_week2_avg(K67.weekALPHA,K67.week);
[K67.week_avgGmax,~,K67.week_avgstdGmax]   = AM_week2_avg(K67.weekGmax,K67.week);
