% % Read ed1.finfile:///home/ncoupe/Documents/AmazonRS/MCD43C1/AM_rs3read_sites_MCD43C1.m
% % Natalia Restrepo-Coupe
% % Tucson AZ lovely summer of 2009 (June 21) Modified December 2016 Toronto
% %........................................................................
% %........................................................................
flag_unix=0;        flag_read=1;        flag_figures=10;     %1 go for figures
flag_qa=1;
if flag_unix==1
    slash='\';
    addpath \Users\ncoupe\Documents\OZ\
    addpath \Users\ncoupe\Documents\Amazon\
elseif flag_unix==0
    K67path='/media/ncoupe/seagate/Data/tapajos/km67data/';
    FolderK67='/home/ncoupe/tapajos/R.codes/out.files/';
    FolderFlux='/home/ncoupe/Documents/Amazon/AllSites/';
    FolderBF='/home/ncoupe/Documents/Amazon/AllSites/';
    slash='/';
    addpath /home/ncoupe/Documents/OZ/;
    addpath /home/ncoupe/Documents/Amazon/;
end
%..........................................................................
ip = 5;
[site_name_BF, lat_BF, lon_BF, site_long_BF, year_ini, year_end, siteelevation_BF,...
    height_BF, canopy_BF, ust_tresh_BF, sunrise_BF, sunset_BF, utm_BF, ust_min_BF, ust_max_BF,...
    LAI_BF]  =  textread([FolderFlux 'sites_loc.txt'],...
    '%s	%f	%f	%s	%f	%f	%f	%f %f	%f	%f	%f	%f	%f	%f	%f',...
    'delimiter',',','headerlines', 1);
site_name = site_name_BF(ip);         lat = lat_BF(ip);         lon = lon_BF(ip);     site_long = site_long_BF(ip);
siteelevation = siteelevation_BF(ip); canopy = canopy_BF(ip);   ust_tresh = ust_tresh_BF(ip);
sunrise = sunrise_BF(ip);             sunset = sunset_BF(ip);   utm = utm_BF(ip);     ust_min = ust_min_BF(ip);
ust_max = ust_max_BF(ip);             LAI = LAI_BF(ip);         height = height_BF(ip);

tresholdSo    = 10;
tresholdCycle = 22;                   tresholdCycleDay = 6;
%..........................................................................
% gC_m2_d1=60.*60.*24.*12.*(10.^-6);
% kgC_ha_m1=1.2e-05 * 10000. * 86400 * 365.25/1000/1000;
gC_m2_d1=60.*60.*24.*12.*(10.^-6);
% kgC_ha_m1=60.*60.*24.*12.*(10.^-6).*220;
%u*22 replaced by u* 0.25 and consequently u*22 by 0.19

%..........................................................................
month_id = ['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'];

% K67FileName = [FolderK67 'km67.eddyflux.30min.txt'];
% K67FileName = [FolderK67 'km67.eddyfluxEA.txt'];
% K67FileName = ['/home/ncoupe/tapajos/R.codes/out.files.h2o/km67.eddyflux.EBco2cal.2005.txt'];
K67FileName = [FolderK67 'km67.eddyflux.2021.txt'];

%..........................................................................
if flag_read==1
    %% Close the text file.
    opts = delimitedTextImportOptions("NumVariables", 91);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["hours", "jdstartgmt", "YYYY", "JD", "HRMI", "ws", "wsdir", "tc", "RH", "tson", "tasp", "fheat", "fmom", "co2", "fco2", "h2ommolm", "fh2o", "h2omix", "pambpa", "H", "LHdry", "ustar", "ppm2umol", "netrad", "snetrad", "netradcal", "par", "spar1up", "par2up", "spar2up", "par1dn", "spar1dn", "parcal", "tair1", "tair2", "tair3", "tair4", "tair5", "tair6", "tair7", "tair8", "ws1", "sws1", "ws2", "sws2", "ws3", "sws3", "ws4", "sws4", "rain", "totsun", "difsun", "e_sat", "RHec", "VPDec", "u", "v", "w", "ustd", "vstd", "wstd", "CO2colwt", "storagewt", "co2pf1", "co2pf2", "co2pf3", "co2pf4", "co2pf5", "co2pf6", "co2pf7", "co2pf8", "H2Ocolwt", "H2Ostoragewt", "h2opf1", "h2opf2", "h2opf3", "h2opf4", "h2opf5", "h2opf6", "h2opf7", "h2opf8", "nee", "neefill", "r", "gee", "parfill", "TOA", "neefilllight", "rlight", "geefilllight", "storagefill"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    k67 = readtable("/home/ncoupe/tapajos/R.codes/out.files/km67.eddyflux.txt", opts);
    k67 = standardizeMissing(k67,[999.99 -999.99]);
    
    %     K67(K67==9.999900000000000e+02) = NaN;
    K67 = table2struct(k67,'ToScalar',true);
    
    %% Clear temporary variables
    clear opts k67
    
    ind = find((K67.jdstartgmt>5530)&(K67.jdstartgmt<5538));    K67.H(ind) = NaN;
    %%
    %..........................................................................
    %Set time vector to local time and start of the time series
    %..........................................................................
    %     K67.dateloc = (datenum(2000,0,K67.jdstartgmt,-utm,0,0));
    K67.dateloc = (datenum(2000,0,K67.jdstartgmt(1),-utm,0,0):1/24:datenum(2000,0,K67.jdstartgmt(end),-utm,0,0))';
    [K67.Y,K67.M,K67.D,K67.HR,K67.MI,K67.SS] = datevec(K67.dateloc);
    % %     ind = find(K67.MI>=55); K67.MI(ind)=0; K67.HR(ind)=1+K67.HR(ind);
    % %     ind = find(K67.MI<=5);  K67.MI(ind)=0;
    
    K67.dateloc_day = datenum(K67.Y,K67.M,K67.D);
    K67.JD   = round(K67.dateloc-(datenum(K67.Y,1,0,K67.HR,K67.MI,K67.SS)));
    
    [K67.So] = AM_solar_rad(K67.JD,K67.HR+1,lat,lon);
    %     [K67.So,~,~,~,~,K67.sza] = AM_solar_rad(K67.JD,K67.HR+1,lat,lon);
    K67.sza = AM_solar_radNOAA(K67.JD,K67.HR+1,K67.Y,lat,lon);
    %run K67._01read_ed1par_cal.m
    
    %......................................................................
    %% Quality
    %......................................................................
    K67.LE = K67.LHdry;
    K67record.fco2 = K67.fco2;    K67record.co2 = K67.co2;  K67record.storagewt = K67.storagewt;
    K67record.nee  = K67.neefill;
    K67record.LE   = K67.LHdry;   K67record.H   = K67.H;    K67record.netrad = K67.netrad;
    
    K67.LE((K67.dateloc>datenum(2008,8,4,12,0,0))&(K67.dateloc<datenum(2008,8,6,12,0,0))) = NaN;
    save('km67.eddyflux.mat','K67');
end
%..........................................................................
%..........................................................................
%Time vector
%..........................................................................
%..........................................................................
K67.a = length(K67.dateloc);
K67.year = unique(K67.Y);

K67.day=(datenum(K67.Y(1),K67.M(1),K67.D(1)):1:datenum(K67.Y(end),K67.M(end),K67.D(end)))';
[K67.dayY,K67.dayM,K67.dayD]=datevec(K67.day);
K67.dayJD=K67.day-(datenum(K67.dayY,1,0));

K67.year=(K67.Y(1):1:K67.Y(end))';     K67.g=length(K67.year);

K67.month=(datenum(K67.Y(1),K67.M(1),1):1:(datenum(K67.Y(end),K67.M(end),1)))';
[~,~,z]=datevec(K67.month);            K67.month(z~=1)=[];

ALLmonth_avg=(datenum(1,1,1):1:datenum(1,12,31))';
[~,~,z]=datevec(ALLmonth_avg);         ALLmonth_avg(z~=1)=[];

ALLmonth_hydro_avg=(datenum(1,8,1,0,0,0):1:datenum(1,7,31,0,0,0))';
[~,~,z]=datevec(ALLmonth_hydro_avg);   ALLmonth_hydro_avg(z~=1)=[];

% % %% ........................................................................
% % K67.netrad((K67.dateloc>datenum(2002,2,15))&(K67.dateloc<datenum(2002,2,27))) = NaN;
% % K67.netrad((K67.dateloc>datenum(2002,1,2))&(K67.dateloc<datenum(2002,1,21))) = NaN;
% % ind = (K67.dateloc>datenum(2002,1,21))&(K67.dateloc<datenum(2002,2,16));
% % K67.netrad(ind) = K67.netrad(ind)+20;
% % K67.par(ind) = K67.par(ind)+40;
% % ind = (K67.dateloc>datenum(2002,1,21))&(K67.dateloc<datenum(2002,2,16)&(K67.TOA==0));
% % K67.par(ind) = 0;
% %
% % K67.netrad((K67.dateloc>datenum(2020,1,1))&(K67.dateloc<datenum(2020,1,10))) = NaN;
% % K67.netrad((K67.dateloc>datenum(2005,2,13))&(K67.dateloc<datenum(2005,3,1))) = NaN;
% % K67.netrad((K67.dateloc>datenum(2012,7,17))&(K67.dateloc<datenum(2012,12,1))&(K67.netrad<0)) = NaN;
% %
% % K67.totsun((K67.dateloc>datenum(2004,8,31.5))&(K67.dateloc<datenum(2004,9,1.5))) = NaN;


K67.tair_fill = K67.tasp;
[~,~,K67.tair_fill] = AM_compute_regression(K67.tson,K67.tair_fill,[.4 .4 .4],0);
K67.tair_fill = AM_fill_day_cycle(K67.tair_fill,K67.dateloc,6,18);
pamb = AM_fill_day_cycle(K67.pambpa,K67.dateloc,6,18);
A0 = 6.107799961;       A1 = 4.436518521e-1;        A2 = 1.428945805e-2;
A3 = 2.650648471e-4;    A4 = 3.031240396e-6;        A5 = 2.034080948e-8;
A6 = 6.136820929e-11;   RV = 0.0004615;             %Gas constant for water vapor [J/(mg K)]
K67.e_sat   = 0.1.*(A0+K67.tair_fill.*(A1+K67.tair_fill.*(A2+K67.tair_fill.*(A3+K67.tair_fill.*(A4+K67.tair_fill.*(A5+K67.tair_fill*A6))))));
%     e_sat2=0.611.*10.^((7.5.*K67.tair_fill)./(237.7+K67.tair_fill));                       %saturation vapor pressure (kPa)
K67.h2o_sat = ((0.622.*K67.e_sat)./(pamb-K67.e_sat)).*1000000;      %%g/kg
% h2o_sat_mmolm = h2o_sat.* 1e-3.*18.016./0.02895; % (g h2o/ kg dry air)

%%
figure('color','white');
subplot(3,1,1);    hold on;                     plot(K67.dateloc,K67.h2o_sat);
plot(K67.dateloc,K67.h2opf2.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2omix);
% plot(K67.dateloc,K67.h2ommolm,'k');
axis([K67.dateloc(1) K67.dateloc(end) 0 50]);   grid on; legend('H2O_{sat}','H2O_{PF}','H2O_{EC}');
datetick('x');                                  ylabel('H_2O (g kg^-^1)');
subplot(3,1,2);    plot(K67.dateloc,K67.RH);	hold on;    plot(K67.dateloc,K67.dateloc.*0+100);
% plot(K67.dateloc,K67.RHpf2.*100);
ylabel('RH (%)');     axis([K67.dateloc(1) K67.dateloc(end) 30 101]);   grid on; datetick('x');
legend('H2O_{EC}','RH==100%','H2O_{PF}');
subplot(3,1,3);    plot(K67.dateloc,K67.tasp);	hold on;    plot(K67.dateloc,K67.tair1);
ylabel('T_{air} (C)');     axis([K67.dateloc(1) K67.dateloc(end) 20 40]);   grid on; datetick('x');
legend('T_{air}','T_{layer 1}');

%%
figure('color','white');
subplot(3,1,1);     hold on;                     plot(K67.jdstartgmt,K67.h2o_sat);
plot(K67.jdstartgmt,K67.h2opf2.* 1e-3.*18.016./0.02895);
plot(K67.jdstartgmt,K67.h2omix);
axis([K67.jdstartgmt(1) K67.jdstartgmt(end) 0 50]); grid on;
legend('H2O_{sat}','H2O_{PF}','H2O_{EC}');          ylabel('H_2O (g kg^-^1)');
subplot(3,1,2);     plot(K67.jdstartgmt,K67.RH);     hold on;
plot(K67.jdstartgmt,K67.dateloc.*0+1);
ylabel('RH (%)');   axis([K67.jdstartgmt(1) K67.jdstartgmt(end) .3 1.1]);
grid on;            legend('H2O_{EC}','RH==100%');
subplot(3,1,3);     plot(K67.jdstartgmt,K67.tasp);	hold on;
plot(K67.jdstartgmt,K67.tair1);
ylabel('T_{air} (C)');     axis([K67.jdstartgmt(1) K67.jdstartgmt(end) 20 40]);
grid on;            legend('T_{air}','T_{layer 1}');

figure('color','white');
subplot(3,1,1);    hold on;     plot(K67.dateloc,K67.So./2);
plot(K67.dateloc,K67.LE);
axis([K67.dateloc(1) K67.dateloc(end) 0 1000]);   grid on; legend('LE','TOA');
datetick('x');                                  ylabel('W m^-^2)');

%% ........................................................................
figure('color','white');
subplot(2,1,1);    hold on;   plot(K67.dateloc,K67.h2opf1.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2opf3.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2opf4.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2opf5.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2opf6.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2opf7.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2opf8.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2o_sat,'color',[.9 .5 .5]);       hold on;
plot(K67.dateloc,K67.h2omix,'color','k');%[.7 .7 .7]);
plot(K67.dateloc,K67.h2opf2.* 1e-3.*18.016./0.02895,'color',[.0 .0 .9]);
% plot(K67.dateloc,K67.h2opf2,'color',[.0 .0 .5]);
% plot(K67.dateloc,K67.h2ommolm,'color',[.5 .0 .5]);
datetick('x'); ylabel('H_2O (g kg^-^1)');
%%
figure('color','white');
ind= find(K67.Y==2002);
subplot(2,1,1);    hold on;   plot(K67.dateloc(ind),K67.h2opf1(ind).* 1e-3.*18.016./0.02895);
plot(K67.dateloc(ind),K67.h2opf3(ind).* 1e-3.*18.016./0.02895);
plot(K67.dateloc(ind),K67.h2opf4(ind).* 1e-3.*18.016./0.02895);
plot(K67.dateloc(ind),K67.h2opf5(ind).* 1e-3.*18.016./0.02895);
plot(K67.dateloc(ind),K67.h2opf6(ind).* 1e-3.*18.016./0.02895);
plot(K67.dateloc(ind),K67.h2opf7(ind).* 1e-3.*18.016./0.02895);
plot(K67.dateloc(ind),K67.h2opf8(ind).* 1e-3.*18.016./0.02895);
hl1 = plot(K67.dateloc(ind),K67.h2o_sat(ind),'color',[.9 .5 .5]);       hold on;
hl2 = plot(K67.dateloc(ind),K67.h2omix(ind),'color','k');%[.7 .7 .7]);
hl3 = plot(K67.dateloc(ind),K67.h2opf2(ind).* 1e-3.*18.016./0.02895,'color',[.0 .0 .9]);
% plot(K67.dateloc(ind),K67.h2opf2(ind),'color',[.0 .0 .5]);
% plot(K67.dateloc(ind),K67.h2ommolm(ind),'color',[.5 .0 .5]);
datetick('x'); ylabel('H_2O (g kg^-^1)'); legend([hl1,hl2,hl3],'sat','EB','PF2');

subplot(2,1,2);    hold on;   plot(K67.jdstartgmt(ind),K67.h2opf1(ind));
% plot(K67.jdstartgmt,K67.h2o_sat,'color',[.9 .5 .5]);       hold on;
% plot(K67.jdstartgmt,K67.h2omix,'color','k');%[.7 .7 .7]);
plot(K67.jdstartgmt(ind),K67.h2opf2(ind),'color',[.0 .0 .9]);
% plot(K67.dateloc,K67.h2opf2,'color',[.0 .0 .5]);
plot(K67.jdstartgmt(ind),K67.h2ommolm(ind),'color',[.5 .0 .5]);
ylabel('H_2O (mmol mol^-^1)');  legend('PF1','PF2','EB');

%% ........................................................................
figure('color','white');
subplot(3,1,1);      plot(K67.dateloc,K67.co2);         hold on;	grid on;    ylabel('CO_2');
axis([K67.dateloc(1) K67.dateloc(end) 300 500]);        datetick('x');
subplot(3,1,2);      plot(K67.dateloc,K67.storagewt);  hold on;    grid on;    ylabel('storage');
axis([K67.dateloc(1) K67.dateloc(end) -40 40]);         datetick('x');
subplot(3,1,3);      plot(K67.dateloc,K67.fco2);        hold on;    grid on;    ylabel('fCO_2');
axis([K67.dateloc(1) K67.dateloc(end) -50 50]);         datetick('x');

%% ........................................................................
figure('color','white');        hold on;
plot(K67.jdstartgmt,400+K67.So./30,'g');
plot(K67.jdstartgmt,K67.netrad,'k');    plot(K67.jdstartgmt,K67.H,'r');
plot(K67.jdstartgmt,K67.LE,'b');        plot(K67.jdstartgmt,K67.H+K67.LE,'m');
grid on;                                ylabel('energy');
axis([K67.jdstartgmt(1) K67.jdstartgmt(end) 300 500]);    %datetick('x');

figure('color','white');
subplot(3,1,1);     hold on;	plot(K67.jdstartgmt,400+K67.So./30,'g');
plot(K67.jdstartgmt,K67.co2);   plot(K67.jdstartgmt,K67record.co2);
grid on;                        ylabel('CO_2');
axis([K67.jdstartgmt(1) K67.jdstartgmt(end) 300 500]);    %datetick('x');
subplot(3,1,2);     hold on;        plot(K67.jdstartgmt,K67record.storagewt);
plot(K67.jdstartgmt,K67.storagewt); grid on;    ylabel('storage');
plot(K67.jdstartgmt,-K67.So./100);
axis([K67.jdstartgmt(1) K67.jdstartgmt(end) -40 40]);     %datetick('x');
subplot(3,1,3);    hold on;     plot(K67.jdstartgmt,K67record.fco2);
plot(K67.jdstartgmt,K67.fco2,'k');  plot(K67.jdstartgmt,-K67.So./100,'m');
bar(K67.jdstartgmt,K67.rain.*10,'r');
%plot(K67.jdstartgmt,K67.h2ommolm.*2-60,'g');
grid on;    ylabel('fCO_2');
axis([K67.jdstartgmt(1) K67.jdstartgmt(end) -50 50]);     %datetick('x');

%% ........................................................................
[K67.storagefill]       = AM_sco2fill_LUT(K67.storagewt,K67.ustar,K67.dateloc);

figure('color','white');            hold on;     plot(K67.jdstartgmt,K67record.nee);
plot(K67.jdstartgmt,K67.nee);
ylabel('NEE (\mumol m^-^2 s^-^1)');

%%  Cal precipitation
run /home/ncoupe/tapajos/process.u.az/K67_02read_BelterraGOamazon.m
% % % 18 December 2003,  100 mL into raingauge for calibration (14:30 LT)
% % % 27 July 2017, total 1000 ml, start: local time 10:50 UTM 13:50, end: local time 11:01 UTM 14:01
% % % 31 July 2017, total 2000 ml, start: local time 14:32 UTM 17:32, end: local time 14:55 UTM 17:55
Belterra.pressure = Belterra.pressure.*100;
Belterra.e_sat = 0.1.*(A0+Belterra.temp.*(A1+Belterra.temp.*(A2+Belterra.temp.*(A3+Belterra.temp.*(A4+Belterra.temp.*(A5+Belterra.temp*A6))))));
% Belterra.e  = Belterra.RH.*Belterra.e_sat./100;
% Belterra.h2omix = 621.9907.*Belterra.e./((Belterra.pressure./1000)-Belterra.e);

%..........................................................................
RHBelterra = interp1(Belterra.dateloc,Belterra.RH,K67.dateloc);
hrBelterra = interp1(Belterra.dateloc,Belterra.hour,K67.dateloc);
% h2omixBelterra = interp1(Belterra.dateloc,Belterra.h2omix,K67w.dateloc);
eBelterra  = RHBelterra.*K67.e_sat./100;
h2omixBelterra = 621.9907.*eBelterra./((pamb./1000)-eBelterra);
dayh2omixBelterra = AM_day(h2omixBelterra,K67.dateloc,1);
%..........................................................................
RHGoAmazon    = interp1(DateTime,RH72,K67.dateloc);
tairGoAmazon  = interp1(DateTime,Tair72,K67.dateloc);
e_satGoAmazon = 0.1.*(A0+tairGoAmazon.*(A1+tairGoAmazon.*(A2+tairGoAmazon.*(A3+tairGoAmazon.*(A4+tairGoAmazon.*(A5+tairGoAmazon*A6))))));
% eGoAmazon      = RHGoAmazon.*e_satGoAmazon./100;
eGoAmazon      = RHGoAmazon.*K67.e_sat./100;
h2omixGoAmazon = 621.9907.*eGoAmazon./(pamb./1000-eGoAmazon);
% QA GoAmazon
ind = find((K67.jdstartgmt>=5845)&(K67.jdstartgmt<=5848));
RHGoAmazon(ind) = NaN;  eGoAmazon(ind) = NaN;   h2omixGoAmazon(ind) = NaN;

%%.........................................................................
figure('color','white');
subplot(2,2,1:2); hold on;  plot(K67.dateloc,K67.rain,'LineWidth',2);
% plot(K67.dateloc,K67.gee_filled./100);
ylabel('rain (mm)');        datetick('x');

DateTime   = datenum(DateTimeGMT0700)+datenum(0,0,0,1,0,0);     %+datenum(0,0,3,0,0);
[DateTimeY,DateTimeM,DateTimeD,DateTimeHR,~,~] = datevec(DateTime);
DateTimeHour = ((datenum(DateTimeY(1),DateTimeM(1),DateTimeD(1),DateTimeHR(1),0,0)):...
    datenum(0,0,0,1,0,0):...
    (datenum(DateTimeY(end),DateTimeM(end),DateTimeD(end),DateTimeHR(end)+1,0,0)))';

Rain71Hour = DateTimeHour.*NaN;     Rain72Hour = DateTimeHour.*NaN;
for ik = 1:length(DateTimeHour)
    ind = find((DateTime>=DateTimeHour(ik))&(DateTime<DateTimeHour(ik)+datenum(0,0,0,0,30,0)));
    Rain71Hour(ik) = nansum(Rain71(ind));
    Rain72Hour(ik) = nansum(Rain72(ind));
end
plot(DateTimeHour,Rain71Hour,'k');
plot(DateTimeHour,Rain72Hour,'r');

K67.Rain71 = interp1(DateTimeHour,Rain71Hour,K67.dateloc);
K67.Rain72 = interp1(DateTimeHour,Rain72Hour,K67.dateloc);
[~,K67.monthRain71] = AM_month(K67.Rain71,K67.dateloc);
[~,K67.monthRain]   = AM_month(K67.rain,K67.dateloc);
[~,K67.monthRain72] = AM_month(K67.Rain72,K67.dateloc);

subplot(2,2,3:4);  hold on;
y = [K67.monthRain,K67.monthRain71,K67.monthRain72];
xlim([datenum(DateTimeY(1),DateTimeM(1),1) datenum(DateTimeY(end),DateTimeM(end)+1,1)]);
for ik=1:3;     bar(K67.month,y(ik));  end
datetick('x'); ylabel('rain (mm)');
legend('k67','GoAmazon 71','GoAmazon 72');     datetick('x');
%%
figure('color','white');
subplot(2,2,1); hold on;  plot(K67.Rain71,K67.rain,'.');
subplot(2,2,2); hold on;  plot(K67.Rain72,K67.rain,'.');
subplot(2,2,3:4); hold on;  plot(K67.dateloc,K67.tasp,'LineWidth',2);
% plot(K67.dateloc,K67.gee_filled./100);
ylabel('T_{air}');     datetick('x')
plot(DateTime,Tair72,'k');

%% ........................................................................
figure('color','white');
ind = find((K67.dateloc>datenum(2017,7,31))&(K67.dateloc<datenum(2017,8,1)));   display(sum(K67.rain(ind)))
% display(K67.jdstartgmt(ind(1)))
subplot(2,2,1); plot(K67.dateloc(ind),K67.rain(ind));   title('31 July 2017');  datetick('x','HHPM')
ind = find((K67.dateloc>datenum(2017,7,30))&(K67.dateloc<datenum(2017,8,2)));
subplot(2,2,2); plot(K67.dateloc(ind),K67.rain(ind));	title('31 July 2017');  datetick('x','HHPM')
K67.rain(ind) = 0;
ind = find((K67.dateloc>datenum(2017,7,27))&(K67.dateloc<datenum(2017,7,28)));  display(sum(K67.rain(ind)))
subplot(2,2,3); plot(K67.dateloc(ind),K67.rain(ind));   title('28 July 2017');  datetick('x','HHPM')
K67.rain(ind) = 0;
ind = find((K67.dateloc>datenum(2003,12,18))&(K67.dateloc<datenum(2003,12,19)));     display(sum(K67.rain(ind)))
subplot(2,2,4); plot(K67.dateloc(ind),K67.rain(ind));   title('18 December 2003');  datetick('x','HHPM')

%% ........................................................................
figure;                             hold on;
plot(K67.dateloc,K67.h2opf2.* 1e-3.*18.016./0.02895);
plot(K67.dateloc,K67.h2omix);
plot(K67.dateloc,K67.h2o_sat);    plot(K67.dateloc,h2omixGoAmazon);
plot(K67.dateloc,h2omixBelterra,'k');
legend('PF','EB','sat','GoAmazon','Belterra');
ylim([ 0 40]);datetick('x'); ylabel('H_2O (g/kg)');

%% ........................................................................
color_mtx = jet(12);
figure('color','white');
subplot(2,2,1);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.tasp, K67.dateloc,6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');  xticks(0:1/4:1);
ylabel('T_{asp}(C)');   xlabel('HR');

subplot(2,2,2);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.tair1, K67.dateloc,6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');  xticks(0:1/4:1);
ylabel('T_{air1}(C)');   xlabel('HR');
subplot(2,2,3);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.tair2, K67.dateloc,6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');  xticks(0:1/4:1);
ylabel('T_{air2}(C)');   xlabel('HR');
subplot(2,2,4);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.tson, K67.dateloc,6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');  xticks(0:1/4:1);
ylabel('T_{csat}(C)');  xlabel('HR');

%% ........................................................................
ind = find((K67.Y==2016)&(K67.M>=7));%|(K67.Y==2018)
color_mtx = jet(12);
figure('color','white');
subplot(2,2,1);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.co2(ind), K67.dateloc(ind),6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');      xticks(0:1/4:1);
ylabel('CO2');              xlabel('HR');

%% ........................................................................
color_mtx = jet(length(K67.year));      K67.YcycleCO2 = NaN(length(K67.year),1)
K67.YcycleCO2min = K67.YcycleCO2;       K67.YcycleCO2max = K67.YcycleCO2;
figure('color','white');
subplot(2,2,1);             hold on;
for ik = 1:length(K67.year)
    ind = find((K67.Y==K67.year(ik)));
    [month_cycle,~,~,~,cycle] = AM_day_cycle(K67.co2(ind), K67.dateloc(ind),6,18);
    %     plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'color',color_mtx(ik,:),'linewidth',2);
    K67.YcycleCO2(ik) = nanmean(cycle);
    K67.YcycleCO2min(ik) = nanmin(month_cycle(:));
    K67.YcycleCO2max(ik) = nanmax(month_cycle(:));
end
K67.YcycleCO2max(7) = NaN;
datetick('x','HH:MM');      xticks(0:1/4:1);
ylabel('CO2');              xlabel('HR');

%% ........................................................................
subplot(2,2,2);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.fco2, K67.dateloc,6,18);
plot((1/24:1/24:1),cycle,'k','linewidth',2);
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.storagewt, K67.dateloc,6,18);
plot((1/24:1/24:1),cycle,'b','linewidth',2);
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.nee, K67.dateloc,6,18);
plot((1/24:1/24:1),cycle,'color',[0 .7 .7],'linewidth',2);
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.fco2+K67.storagewt, K67.dateloc,6,18);
plot((1/24:1/24:1),cycle,'color',[.7 0 .7],'linewidth',1);
datetick('x','HH:MM');  xticks(0:1/4:1);
ylabel('CO_{2}flux');   xlabel('HR');   legend('fco_{2}','Sco2','nee','fco2+Sco2')
%%
color_mtx = jet(12);
ind = find((K67.Y==2002)&(K67.M>0));
subplot(2,2,3);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.co2(ind), K67.dateloc(ind),6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');      xticks(0:1/4:1);
ylabel('CO2');              xlabel('HR');
%%
ind = find((K67.Y==2016)&(K67.M<7));
subplot(2,2,4);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.co2(ind), K67.dateloc(ind),6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');      xticks(0:1/4:1);
ylabel('CO2');              xlabel('HR');
%%
figure('color','white');
%%
ind = find((K67.Y==2003)&(K67.M<7));
subplot(2,2,1);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.h2ommolm(ind), K67.dateloc(ind),6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');      xticks(0:1/4:1);
ylabel('H_2O (mmol m)');              xlabel('HR');
%%
ind = find((K67.Y==2015)&(K67.M<7));
subplot(2,2,4);             hold on;
[cycle_month,~,~,~,cycle] = AM_day_cycle(K67.h2ommolm(ind), K67.dateloc(ind),6,18);
for ik = 1:12
    plot((1/24:1/24:1),cycle_month(ik,:),'color',color_mtx(ik,:));
    plot((1/24:1/24:1),cycle,'k','linewidth',2);
end
datetick('x','HH:MM');      xticks(0:1/4:1);
ylabel('H_2O (mmol m)');              xlabel('HR');
%%
close all;
%..........................................................................
%%  Check figures
%..........................................................................
if flag_figures==1
    %AM_plot_bimonth(K67.ws,K67.So,K67.Y,K67.M,K67.HR,'K67.','wsp (m s-1)');
    %close all;
    AM_plot_bimonth(K67.tamb,K67.So,K67.Y,K67.M,K67.HR,'K67.','tair (degC)');
    close all;
    AM_plot_bimonth(K67.tson,K67.So,K67.Y,K67.M,K67.HR,'K67.','tson (degC)');
    close all;
    AM_plot_bimonth(K67.tasp,K67.So,K67.Y,K67.M,K67.HR,'K67.','tasp (degC)');
    close all;
    AM_plot_bimonth(K67.fheat,K67.So,K67.Y,K67.M,K67.HR,'K67.','fheat (m s-1)');
    close all;
    AM_plot_bimonth(K67.fmom,K67.So,K67.Y,K67.M,K67.HR,'K67.','fmomentum (m s-1)');
    close all;
    AM_plot_bimonth(K67.co2,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2 (ppm)');
    close all;
    AM_plot_bimonth(K67.fco2,K67.So,K67.Y,K67.M,K67.HR,'K67.','fco2 (\mumolCO_2m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.h2ommolm,K67.So,K67.Y,K67.M,K67.HR,'K67.','H2O (mmol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.fh2o,K67.So,K67.Y,K67.M,K67.HR,'K67.','fH2O (mmol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.pambpa,K67.So,K67.Y,K67.M,K67.HR,'K67.','Pamb (Pa)');
    close all;
    AM_plot_bimonth(K67.H,K67.So,K67.Y,K67.M,K67.HR,'K67.','H (W m^{-2})');
    close all;
    AM_plot_bimonth(K67.LE,K67.dateloc,'K67.','LE (W m^{-2})');
    close all;
    AM_plot_bimonth(K67.netrad,K67.So,K67.Y,K67.M,K67.HR,'K67.','NetRad (W m^{-2})');
    close all;
    AM_plot_bimonth(K67.par,K67.So,K67.Y,K67.M,K67.HR,'K67.','PAR (mmol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.par2up,K67.So,K67.Y,K67.M,K67.HR,'K67.','PAR2up (mmol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.par1dn,K67.So,K67.Y,K67.M,K67.HR,'K67.','PAR1dn (mmol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.tair1,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair1 (degC)');
    close all;
    AM_plot_bimonth(K67.tair2,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair2 (degC)');
    close all;
    AM_plot_bimonth(K67.tair3,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair3 (degC)');
    close all;
    AM_plot_bimonth(K67.tair4,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair4 (degC)');
    close all;
    AM_plot_bimonth(K67.tair5,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair5 (degC)');
    close all;
    AM_plot_bimonth(K67.tair6,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair6 (degC)');
    close all;
    AM_plot_bimonth(K67.tair7,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair7 (degC)');
    close all;
    AM_plot_bimonth(K67.tair8,K67.So,K67.Y,K67.M,K67.HR,'K67.','Tair8 (degC)');
    close all;
    AM_plot_bimonth(K67.ws1,K67.So,K67.Y,K67.M,K67.HR,'K67.','ws1 (m s^{-2})');
    close all;
    AM_plot_bimonth(K67.ws2,K67.So,K67.Y,K67.M,K67.HR,'K67.','ws2 (m s^{-2})');
    close all;
    AM_plot_bimonth(K67.ws3,K67.So,K67.Y,K67.M,K67.HR,'K67.','ws3 (m s^{-2})');
    close all;
    AM_plot_bimonth(K67.ws4,K67.So,K67.Y,K67.M,K67.HR,'K67.','ws4 (m s^{-2})');
    close all;
    AM_plot_bimonth(K67.totsun,K67.So,K67.Y,K67.M,K67.HR,'K67.','TotalRad (W m^{-2})');
    close all;
    AM_plot_bimonth(K67.difsun,K67.So,K67.Y,K67.M,K67.HR,'K67.','DiffRad (W m^{-2})');
    close all;
    AM_plot_bimonth(K67.storagewt,K67.So,K67.Y,K67.M,K67.HR,'K67.','Sco2 (/mumol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.co2pf1,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv1 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf2,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv2 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf3,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv3 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf4,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv4 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf5,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv5 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf6,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv6 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf7,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv7 (pm)');
    close all;
    AM_plot_bimonth(K67.co2pf8,K67.So,K67.Y,K67.M,K67.HR,'K67.','CO2.lv8 (pm)');
    close all;
    %    AM_plot_bimonth(K67.storagewt,K67.ee,K67.Y,K67.M,K67.HR,'K67.','NEE (/mumol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.storagewt,K67.r,K67.Y,K67.M,K67.HR,'K67.','R (/mumol m^{-2}s^{-1})');
    close all;
    AM_plot_bimonth(K67.storagewt,K67.gee,K67.Y,K67.M,K67.HR,'K67.','GEE (/mumol m^{-2}s^{-1})');
    close all;
    
    ix =   AM_cycle_composite(K67.H,4);
    %     AM_cycle_check(K67.ws,K67.dateloc, 7, 18,'  ws');         close all;
    % %     AM_cycle_check_day(K67.ws,K67.dateloc,'ws','  K67.');
    %     AM_cycle_check(K67.tamb,K67.dateloc, 7, 18,'  Tamb');     close all;
    %     AM_cycle_check(K67.tdew,K67.dateloc,'  Tdew');            close all;
    %     AM_cycle_check(K67.tson,K67.dateloc,7, 18,'  Tson');      close all;
    %     AM_cycle_check(K67.tasp,K67.dateloc,7, 18,'  Tasp');      close all;
    %     AM_cycle_check(K67.fheat,K67.dateloc, 7, 18,'  fheat');   close all;
    %     AM_cycle_check(K67.fmom,K67.dateloc, 7, 18,'  fmom');     close all;
    AM_cycle_check(K67.co2,K67.dateloc, 7, 18,'  co2');       close all;
    AM_cycle_check(K67.fco2,K67.dateloc, 7, 18,'  fco2');     close all;
    %     AM_cycle_check(K67.h2ommolm,K67.dateloc,7, 18,'  h2o');   close all;
    %     AM_cycle_check(K67.fh2o,K67.dateloc,7, 18,'  fh2o');      close all;
    AM_cycle_check(K67.h2omix,K67.dateloc,7, 18,'  h2o.mix'); close all;
    %     AM_cycle_check(K67.pambpa,K67.dateloc, 7, 18,'  Pamb');   close all;
    %     AM_cycle_check(K67.H,K67.dateloc, 7, 18,'  H');           close all;
    AM_cycle_check(K67.LE,K67.dateloc, 7, 18,'  LE');         close all;
    AM_cycle_check(K67.H,K67.dateloc, 7, 18,'  H');         close all;
    %     AM_cycle_check(K67.ustar,K67.dateloc, 7, 18,'  ustar');   close all;
    %     AM_cycle_check(K67.ppm2umol,K67.dateloc,7, 18,'  ppm2mol'); close all;
    AM_cycle_check(K67.netrad,K67.dateloc, 7, 18,'  NetRad'); close all;
    %     AM_cycle_check(K67.par,K67.dateloc, 7, 18,'  PAR');       close all;
    %     AM_cycle_check(K67.par2up,K67.dateloc, 7, 18,'  PAR2up'); close all;
    %     AM_cycle_check(K67.par1dn,K67.dateloc, 7, 18,'  PAR1dn'); close all;
    %     AM_cycle_check(K67.tair1,K67.dateloc, 7, 18,'  Tair1');   close all;
    %     AM_cycle_check(K67.tair2,K67.dateloc,7, 18,'  Tair2');    close all;
    %     AM_cycle_check(K67.tair3,K67.dateloc, 7, 18,'  Tair3');   close all;
    %     AM_cycle_check(K67.tair4,K67.dateloc,7, 18,'  Tair4');    close all;
    %     AM_cycle_check(K67.tair5,K67.dateloc, 7, 18,'  Tair5');   close all;
    %     AM_cycle_check(K67.tair6,K67.dateloc,7, 18,'  Tair6');    close all;
    %     AM_cycle_check(K67.tair7,K67.dateloc, 7, 18,'  Tair7');   close all;
    %     AM_cycle_check(K67.tair8,K67.dateloc,7, 18,'  Tair8');    close all;
    %     AM_cycle_check(K67.ws1,K67.dateloc,7, 18,'  WS1');        close all;
    %     AM_cycle_check(K67.ws2,K67.dateloc,7, 18,'  WS2');        close all;
    %     AM_cycle_check(K67.ws3,K67.dateloc,7, 18,'  WS3');        close all;
    %     AM_cycle_check(K67.ws4,K67.dateloc,7, 18,'  WS4');        close all;
    %     AM_cycle_check(K67.rain,K67.dateloc, 7, 18,'  Rain');     close all;
    %     AM_cycle_check(K67.totsun,K67.dateloc, 7, 18,'  TotRad'); close all;
    %     AM_cycle_check(K67.difsun,K67.dateloc, 7, 18,'  DifRad'); close all;
    %
    AM_cycle_check(K67.co2colwt,K67.dateloc, 7, 18,'  col2col.wt');   close all;
    AM_cycle_check(K67.storagewt,K67.dateloc, 7, 18,'  Sco2.wt');     close all;
    %     AM_cycle_check(K67.co2pf1,K67.dateloc, 7, 18,'  co2.lv1'); close all;
    %     AM_cycle_check(K67.co2pf2,K67.dateloc, 7, 18,'  co2.lv2'); close all;
    %     AM_cycle_check(K67.co2pf3,K67.dateloc, 7, 18,'  co2.lv3'); close all;
    %     AM_cycle_check(K67.co2pf4,K67.dateloc, 7, 18,'  co2.lv4'); close all;
    %     AM_cycle_check(K67.co2pf5,K67.dateloc, 7, 18,'  co2.lv5'); close all;
    %     AM_cycle_check(K67.co2pf6,K67.dateloc, 7, 18,'  co2.lv6'); close all;
    %     AM_cycle_check(K67.co2pf7,K67.dateloc, 7, 18,'  co2.lv7'); close all;
    %     AM_cycle_check(K67.co2pf8,K67.dateloc, 7, 18,'  co2.lv8'); close all;
    %
    %     AM_cycle_check(K67.nee,K67.dateloc,7, 18,'  NEE');        close all;
    %    AM_cycle_check(K67.r,K67.dateloc,7, 18,'  R'); close all;
    %    AM_cycle_check(K67.gee,K67.dateloc,7, 18,'  GEE'); close all;
    %    AM_cycle_check(K67.stg_filled,K67.dateloc,7, 18,'  Sco2.filled'); close all;
    %    AM_cycle_check(K67.par_filled,K67.dateloc,7, 18,'  PAR.filled'); close all;
    %    AM_cycle_check(K67.netrad_filled,K67.dateloc,7, 18,'  NetRad.filled'); close all;
    %     AM_cycle_check(hr_2,K67.dateloc,7, 18,'  hr2'); close all;
    
    %     AM_cycle_check_two(difsun, DifRad_noQAQC, K67.dateloc, 7, 18,' DifRad');
    %     AM_cycle_check_two(totsun, TotRad_noQAQC, K67.dateloc, 7, 18,' TotRad');
    %     AM_cycle_check_two(h2ommolm, h2ommolm_noQAQC, K67.dateloc, 7, 18,' DifRad');
    %     AM_cycle_check_two(par, par_noQAQC,K67.dateloc, 7, 18,' PAR');
    %     AM_cycle_check_two(fco2, fco2_noQAQC,K67.dateloc, 7, 18,' fco2');
    %     AM_cycle_check_two(co2, co2_noQAQC,K67.dateloc, 7, 18,' co2');
    %     AM_cycle_check_two(LE, LE_noQAQC,K67.dateloc, 7, 18,' LE');
    %     AM_cycle_check_two(K67.netrad, NetRad_noQAQC,K67.dateloc, 7, 18,' NetRad');
    %     AM_cycle_check_two(K67.netrad, NetRad_noQAQC,K67.dateloc, 7, 18,' NetRad');
end

%% ........................................................................
figure('color','white');
ind = find((K67.dateloc>datenum(2012,11,1))&(K67.dateloc<datenum(2013,1,1)));
subplot(2,2,1:2);               hold on;
hl1 = plot(K67.dateloc(ind),K67.netrad(ind),'color',[0 0 .9],'Marker','.','LineWidth',1);
hl2 = plot(K67.dateloc(ind),K67.par(ind),'color',[0.9 0 .5],'Marker','.','LineWidth',1);
datetick('x'); ylabel('radiation'); grid on;    legend('Rn (W m^-^2)','PAR (\mumol m^-^2 s^-^1)');

ind = find((K67.dateloc>datenum(2015,11,1))&(K67.dateloc<datenum(2016,1,1)));
subplot(2,2,3:4);               hold on;
hl1 = plot(K67.dateloc(ind),K67.netrad(ind),'color',[0 0 .9],'Marker','.','LineWidth',1);
hl2 = plot(K67.dateloc(ind),K67.par(ind),'color',[0.9 0 .5],'Marker','.','LineWidth',1);
datetick('x'); ylabel('radiation'); grid on;    legend('Rn (W m^-^2)','PAR (\mumol m^-^2 s^-^1)');

%% ........................................................................
figure('color','white');
ind = find((K67.dateloc>datenum(2015,1,1))&(K67.dateloc<datenum(2015,3,1)));
subplot(2,2,1:2);               hold on;
hl1 = plot(K67.dateloc(ind),K67.netrad(ind),'color',[0 0 .9],'Marker','.','LineWidth',1);
hl2 = plot(K67.dateloc(ind),K67.H(ind),'color',[0.9 0 .5],'Marker','.','LineWidth',1);
hl2 = plot(K67.dateloc(ind),K67.LE(ind),'color',[0. 0.9 .5],'Marker','.','LineWidth',1);
datetick('x'); grid on;         ylabel('flux (W m^-^2)');       legend('Rn','H','LE');

%% ........................................................................
figure('color','white');
ind = find((K67.dateloc>datenum(2019,1,1))&(K67.dateloc<datenum(2019,12,31)));
subplot(2,2,1:2);               hold on;
hl4 = plot(K67.jdstartgmt(ind),K67.h2ommolm(ind),'color',[0 0 .0],'Marker','.','LineWidth',1);
hl10 = plot(K67.jdstartgmt(ind),K67.h2omix(ind),'color',[0.9 0 .0],'Marker','.','LineWidth',1);
hl1 = plot(K67.jdstartgmt(ind),K67.h2opf1(ind),'color',[0 0 .9],'Marker','.','LineWidth',1);
hl2 = plot(K67.jdstartgmt(ind),K67.h2opf2(ind),'color',[0.9 0 .5],'Marker','.','LineWidth',1);
hl3 = plot(K67.jdstartgmt(ind),K67.h2opf3(ind),'color',[0. 0.9 .5],'Marker','.','LineWidth',1);
hl6 = plot(K67.jdstartgmt(ind),K67.h2opf4(ind),'color',[0. 0.5 .5],'Marker','.','LineWidth',1);
hl7 = plot(K67.jdstartgmt(ind),K67.h2opf5(ind),'color',[0.2 0.6 .5],'Marker','.','LineWidth',1);
hl8 = plot(K67.jdstartgmt(ind),K67.h2opf6(ind),'color',[0.6 0. .5],'Marker','.','LineWidth',1);
grid on;         ylabel('h2o_{pf} (ppt)');      %datetick('x');
legend('EB','h2opf1','h2opf2','h2opf3');

subplot(2,2,3);             hold on;
hl4 = plot(K67.HR(ind),K67.h2ommolm(ind),'.','color',[0 0 .0]);
hl1 = plot(K67.HR(ind),K67.h2opf1(ind),'.','color',[0 0 .9]);
hl2 = plot(K67.HR(ind),K67.h2opf2(ind),'.','color',[0.9 0 .5]);
hl3 = plot(K67.HR(ind),K67.h2opf3(ind),'.','color',[0. 0.9 .5]);
xticks([0:6:24]);           grid on;    ylabel('h2o_{pf} (ppt)');

%% ........................................................................
ws1Cycle =  AM_cycle(K67.ws1,K67.dateloc);
ws2Cycle =  AM_cycle(K67.ws2,K67.dateloc);
ws3Cycle =  AM_cycle(K67.ws3,K67.dateloc);
ws4Cycle =  AM_cycle(K67.ws4,K67.dateloc);
WindHeight = [64.1;52;38.22;30.66];
color_matrix = jet(25);

figure('color','white');
subplot(2,2,1);     hold on;
plot(1:24,ws1Cycle,'color',color_matrix(1,:));
plot(1:24,ws2Cycle,'color',color_matrix(5,:));
plot(1:24,ws3Cycle,'color',color_matrix(9,:));
plot(1:24,ws4Cycle,'color',color_matrix(13,:));

subplot(2,2,2);     hold on;
wsCycle=[ws1Cycle,ws2Cycle,ws3Cycle,ws4Cycle];
for ik=1:24
    plot(WindHeight,wsCycle(ik,:),'o','color',color_matrix(ik,:));
end
plot(WindHeight,wsCycle(10,:),'ko');

y = (wsCycle(10,:))';
x = WindHeight;

% % f = fit(x,y,'poly1');   % 'exp1'
% % f = fit(x,y,'exp1');
[f,r] = AM_compute_regression_zero(x,y,[0.7 0. 0.],0);
z = f.*x;

subplot(2,2,3);
plot(x,y,'o');  hold on;
% plot(x, f(B,x), '-r') %
% plot(f,x,y);
plot(x,z,'-r');     plot(x,z,'.r');
xlim([0,70]);
grid on;            xlabel('height');        ylabel('wind speed (m s^-^1)')

xLAI = (1.5:1.5:6)';         x = (10:10:40)';
% y = (f.a.*x.^f.b);
% y = (f.a.*x);
z = f.*x;
subplot(2,2,4);
plot(xLAI,z,'o');  hold on;
grid on;           xlabel('LAI');        ylabel('wind speed (m s^-^1)')
%%
color_mtx = jet(4);
figure('color','white');
plot(K67.jdstartgmt,K67.ws,'k'); hold on;
plot(K67.jdstartgmt,K67.ws1,'color',color_mtx(1,:));    plot(K67.jdstartgmt,K67.ws2,'color',color_mtx(2,:))
plot(K67.jdstartgmt,K67.ws3,'color',color_mtx(3,:));    plot(K67.jdstartgmt,K67.ws4,'color',color_mtx(4,:))
%%
figure('color','white');
plot(K67.jdstartgmt,K67.netrad,'k');        hold on;
plot(K67.jdstartgmt,K67.netradcal,'g');     plot(K67.jdstartgmt,K67.par./2,'m');
plot(K67.jdstartgmt,K67.H+K67.LE,'r');      plot(K67.jdstartgmt,K67.LE,'b'); hold on;
%%
figure('color','white');
plot(K67.dateloc,K67.netrad,'k');        hold on;
plot(K67.dateloc,K67.netradcal,'g');     plot(K67.dateloc,K67.par./2,'m');
plot(K67.dateloc,K67.H+K67.LE,'r');      plot(K67.dateloc,K67.LE,'b'); hold on;
plot(K67.dateloc,K67.H,'c');
% datetick('x');
% plot(K67.week,AM_fill(K67.weekNetRad),'k','LineWidth',2);
%%
figure
iy = K67.H+K67.LE;  iy(K67.So>0)=NaN;
plot(iy,K67.netrad,'.');

%%
figure('color','white');            hold on;
plot(K67.jdstartgmt,K67.fco2,'c');
plot(K67.jdstartgmt,K67.gee,'k');
plot(K67.jdstartgmt,K67.geefilllight,'g');     plot(K67.jdstartgmt,K67.nee,'m');
plot(K67.jdstartgmt,K67.rlight,'r');           plot(K67.jdstartgmt,K67.r,'b'); hold on;
datetick('x');      legend('fco2','gee','gee_{fill light}','nee','r_{light}','r_{eco}');
%%
figure('color','white');            hold on;
plot(K67.dateloc,K67.fco2,'c');
plot(K67.dateloc,K67.gee,'k');
plot(K67.dateloc,K67.geefilllight,'g');     plot(K67.dateloc,K67.nee,'m');
plot(K67.dateloc,K67.rlight,'r');           plot(K67.dateloc,K67.r,'b'); hold on;
datetick('x');

%% ........................................................................
figure('color','white');
plot(K67.jdstartgmt,K67.pambpa,'k');        hold on;

%% ........................................................................
color_mtx = jet(9);
ix = K67.co2;     ix((K67.dateloc<datenum(2003,12,31))|(K67.dateloc>datenum(2005,1,1))) = NaN;
[K67.co2_recurrentYR,K67.co2_arecurrent] = AM_cycle_repeat(ix,K67.dateloc);

figure('color','white');
subplot(3,1,1:3);       hold on;
plot(K67.jdstartgmt,K67.co2pf1,'color',color_mtx(1,:)); plot(K67.jdstartgmt,K67.co2pf3,'color',color_mtx(3,:));
plot(K67.jdstartgmt,K67.co2pf4,'color',color_mtx(4,:)); plot(K67.jdstartgmt,K67.co2pf5,'color',color_mtx(5,:));
plot(K67.jdstartgmt,K67.co2pf6,'color',color_mtx(6,:)); plot(K67.jdstartgmt,K67.co2pf7,'color',color_mtx(7,:));
plot(K67.jdstartgmt,K67.co2pf8,'color',color_mtx(8,:)); plot(K67.jdstartgmt,K67.co2pf2,'color',[.0 .0 .9]);
plot(K67.jdstartgmt,K67.co2,'color','k');%[.7 .7 .7]);
plot(K67.jdstartgmt,K67.co2_recurrentYR,'color','m');%[.7 .7 .7]);
ylabel('CO_2 (ppm)');
legend('top','3','4','5','6','7','ground','2','eddy');  legend box off;

%% ........................................................................
close all
figure('color','white');
subplot(3,1,1:3);       hold on;
plot(K67.jdstartgmt,K67.par,'color',[0.0 0.3 0.7]);     plot(K67.jdstartgmt,K67.par2up,'color',[0.3 0.7 0]);
plot(K67.jdstartgmt,K67.par1dn,'color',[.9 0 0]);


K67.radiation = K67.par;
[~,~,K67.radiation] = AM_compute_regression(K67.radiation,K67.netrad,[.4 .4 .4],0);
K67.radiation((K67.radiation<0)|(K67.So==0)) = 0;


%% ........................................................................
K67.par1dn_sample = K67.par1dn;
ind = find((K67.jdstartgmt>=993.)&(K67.jdstartgmt<= 995.646));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1081.)&(K67.jdstartgmt<= 1084));        K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1412.438)&(K67.jdstartgmt<= 1416));     K67.par1dn_sample(ind) = NaN;
% ind = find((K67.jdstartgmt>=2100.)&(K67.jdstartgmt<= 2300.0));      K67.par1dn_sample(ind) = NaN; % no good quality

ind = find((K67.jdstartgmt>=3244.0)&(K67.jdstartgmt<= 3392.0));     K67.par1dn_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=3610.)&(K67.jdstartgmt<=3640));         K67.par1dn_sample(ind) = NaN;
% ind = find((K67.jdstartgmt>=3652.)&(K67.jdstartgmt<=3658.583));     K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3973.0)&(K67.jdstartgmt<=3985.0));      K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4185.0)&(K67.jdstartgmt<=4200.688));    K67.par1dn_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=4476)&(K67.jdstartgmt<=4486));          K67.par1dn_sample(ind) = NaN;
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

ind = find((K67.jdstartgmt>=3272.0)&(K67.jdstartgmt<=3273.0));      K67.par2up_sample(ind) = NaN;
% % ind = find((K67.jdstartgmt>=5994.0)&(K67.jdstartgmt<=6002.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7227.0)&(K67.jdstartgmt<=7233.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7541.0)&(K67.jdstartgmt<=7560.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=7583.0)&(K67.jdstartgmt<=7586.0));      K67.par2up_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=1273.0)&(K67.jdstartgmt<=1275.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1384.0)&(K67.jdstartgmt<=1385.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1458.0)&(K67.jdstartgmt<=1459.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=1875.0)&(K67.jdstartgmt<=1877.0));      K67.par2up_sample(ind) = NaN;

ind = find((K67.jdstartgmt>=3610.0)&(K67.jdstartgmt<=3630.0));      K67.par2up_sample(ind) = NaN;
ind = find((K67.jdstartgmt>=3626.0)&(K67.jdstartgmt<=3633.0));      K67.par2up_sample(ind) = NaN;
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
plot(K67.jdstartgmt,K67.par2up,'color',[0.3 0.7 0]);
plot(K67.jdstartgmt,K67.par2up_sample,'color',[.9 .7 0.0]);
legend('radiation','par_{2up recurrent}','par_{year recurrent}','par_{2cpy}','par_{2cal clean}');

subplot(3,1,2);                     hold on;
plot(K67.jdstartgmt,K67.radiation./12,'color','m');
plot(K67.jdstartgmt,K67.PAR1dn_arecurrent,'color','k');
plot(K67.jdstartgmt,K67.PAR1dn_recurrent,'color',[0.0 0.3 0.7]);
plot(K67.jdstartgmt,K67.par1dn,'color',[0.3 0.7 0]);
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
ind = find((K67.dateloc>=datenum(2010,6,1))&(K67.dateloc<=datenum(2011,7,1)));   %dates set to clean and calibrate
plot(K67.jdstartgmt(ind),K67.par1dn_sample(ind).*20,'color',[.8 0.8 0.0]);
legend('par','par_{cpy}','par_{dwn incl cal}','par_{dwn}');
% K67.jdstartgmt(ind(1))

subplot(2,1,2);       hold on;
plot(K67.dateloc,K67.So./1.5,'color',[0.6 0.6 0.6]);
plot(K67.dateloc,K67.netrad.*3,'color',[0.8 0. 0.8]);       plot(K67.dateloc,K67.par,'color',[0.0 0.3 0.7]);
plot(K67.dateloc,K67.par2up_sample,'color',[0.3 0.7 0],'linewidth',2);
plot(K67.dateloc,K67.par1dn_sample.*20,'color',[.9 0 0]);
ind = find((K67.dateloc>=datenum(2010,6,1))&(K67.dateloc<=datenum(2011,7,1)));   %dates set to clean and calibrate
plot(K67.dateloc(ind),K67.par1dn_sample(ind).*20,'color',[.8 0.8 0.0]);
% plot(K67am.dateloc,K67am.PAR,'color','k','LineWidth',1);
legend('TOA','netrad','par','par_{2up}','par_{dwn}');
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
plot(K67.week,K67.weekPAR2up,'color',[0.0 0.7 0.]);     datetick('x');
legend('PAR','PAR1dn','PAR2up')

subplot(3,1,3);       hold on;
plot(K67.weekjdstartgmt,K67.weekPAR./4,'color',[0.7 0.3 0.]);    plot(K67.weekjdstartgmt,K67.weekPAR1dn,'color',[0.0 0.3 0.7]);
plot(K67.weekjdstartgmt,K67.weekPAR2up,'color',[0.0 0.7 0.]);
legend('PAR','PAR1dn','PAR2up')

%%
figure('color','white');
v1 = .0;         v2 = .95;       v3 = v1:((v2-v1)/4):v2;
subplot(3,1,1);       hold on;
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

%% ........................................................................
v1 = .03;         v2 = .07;       v3 = v1:((v2-v1)/4):v2;
subplot(3,1,2)
[K67.weekRadiation,~,~,~,~,K67.week_stdrRadiation] = AM_week2cycle_rs(K67.radiation,K67.dateloc,tresholdCycleDay,4,2);
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


%mmol mol-1 >> 
% .* 1e-3.*18.016./0.02895