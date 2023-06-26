% % Energy balance components calculatins
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

load('/home/ncoupe/Documents/AndesAmazon/matlabBFallocation.mat')

%..........................................................................
ip=5;
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
HeightEC = 57.79;
%..........................................................................
% gC_m2_d1=60.*60.*24.*12.*(10.^-6);
% kgC_ha_m1=1.2e-05 * 10000. * 86400 * 365.25/1000/1000;
gC_m2_d1=60.*60.*24.*12.*(10.^-6);
% kgC_ha_m1=60.*60.*24.*12.*(10.^-6).*220;
%u*22 replaced by u* 0.25 and consequently u*22 by 0.19

%% ........................................................................
% [~,K67.tair1] = AM_spike_detection(K67.tair1,5);
% [~,K67.tair6] = AM_spike_detection(K67.tair6,5);   plot(K67.jdstartgmt,K67.tair6);
% [~,K67.tair3] = AM_spike_detection(K67.tair3,5);

%% ........................................................................
color_jet = jet(8);
figure('color','white');    hold on;
plot(K67.dateloc,K67.tair1,'color',color_jet(1,:));  plot(K67.dateloc,K67.tair2,'color','k');%color_jet(2,:)
plot(K67.dateloc,K67.tair3,'color',color_jet(3,:));  plot(K67.dateloc,K67.tair4,'color',color_jet(4,:));
plot(K67.dateloc,K67.tair5,'color',color_jet(5,:));  plot(K67.dateloc,K67.tair6,'color',color_jet(6,:));
plot(K67.dateloc,K67.tair7,'color',color_jet(7,:));  plot(K67.dateloc,K67.tair8,'color',color_jet(8,:));
plot(K67.dateloc,K67.tson,'color','m');     ylabel('Temperature (C)'); 
legend('1','2','3','4','5','6','7','8','tson');   datetick('x')
%%
figure('color','white');    hold on;
plot(K67.jdstartgmt,K67.tair1,'color',color_jet(1,:));  plot(K67.jdstartgmt,K67.tair2,'color','k');%color_jet(2,:)
plot(K67.jdstartgmt,K67.tair3,'color',color_jet(3,:));  plot(K67.jdstartgmt,K67.tair4,'color',color_jet(4,:));
plot(K67.jdstartgmt,K67.tair5,'color',color_jet(5,:));  plot(K67.jdstartgmt,K67.tair6,'color',color_jet(6,:));
plot(K67.jdstartgmt,K67.tair7,'color',color_jet(7,:));  plot(K67.jdstartgmt,K67.tair8,'color',color_jet(8,:));
plot(K67.jdstartgmt,K67.tson,'color','m');  ylabel('Temperature (C)'); 
legend('1','2','3','4','5','6','7','8','tson');
%%
[~,~,K67.tair1_fill] = AM_compute_regression(K67.tair_fill,K67.tair1,[.4 .4 .4],0);
[~,~,K67.tair2_fill] = AM_compute_regression(K67.tair_fill,K67.tair2,[.4 .4 .4],0);
[~,~,K67.tair3_fill] = AM_compute_regression(K67.tair_fill,K67.tair3,[.4 .4 .4],0);
[~,~,K67.tair4_fill] = AM_compute_regression(K67.tair_fill,K67.tair4,[.4 .4 .4],0);
[~,~,K67.tair5_fill] = AM_compute_regression(K67.tair_fill,K67.tair5,[.4 .4 .4],0);
[~,~,K67.tair6_fill] = AM_compute_regression(K67.tair_fill,K67.tair6,[.4 .4 .4],0);
[~,~,K67.tair7_fill] = AM_compute_regression(K67.tair_fill,K67.tair7,[.4 .4 .4],0);
[~,~,K67.tair8_fill] = AM_compute_regression(K67.tair_fill,K67.tair8,[.4 .4 .4],0);
%..........................................................................

%..........................................................................
%% QA by energy balance closure
%..........................................................................
ind = find(((K67.par>=-1)&(K67.par<=1))&((K67.netrad>=100)));
figure('color','white');    plot(K67.jdstartgmt,K67.netrad,'-r');  hold on;
plot(K67.jdstartgmt,K67.par,'-m');  plot(K67.jdstartgmt(ind),K67.netrad(ind),'.g');
%%
% ix = (K67.netrad-K67.LE-K67.H);       ind = find(abs(ix)>=100);
ix = (K67.netrad-K67.LE-K67.H);         ind = find((abs(K67.netrad)<10)&(K67.H>=100));
figure('color','white');                plot(K67.jdstartgmt,K67.netrad,'-r');
hold on;    plot(K67.jdstartgmt,K67.H,'-m');
plot(K67.jdstartgmt,K67.LE,'-b');
plot(K67.jdstartgmt(ind),K67.netrad(ind),'.g','MarkerSize',20);
legend ('Net Rad','H','LE');
%% ........................................................................
K67.year = unique(K67.Y);       color_mtx = jet(length(K67.year));
figure('color','white');
subplot(2,2,1);
y1 = K67.LE+K67.H;
% [p,r,x]=AM_compute_regression(K67.netrad,y1,[.8 .8 .8],1);
hold on;
record = [];
for ik = 1:length(K67.year)
    ind = find((~isnan(K67.netrad))&(~isnan(y1))&(K67.Y==K67.year(ik)));
    %     [p,r,x]=AM_compute_regression(K67.netrad(ind),y1(ind),color_mtx(ik,:),1);
    p  = polyfit(K67.netrad(ind),y1(ind),1);
    x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
    y   = polyval(p,x);    plot (x,y,'-','LineWidth',2,'Color','k');
    plot(K67.netrad(ind),y1(ind),'.','color',[0.5 .5 .5]);
    E   = 3.*std(y)*ones(size(x)); %error bars that are two standard deviation units
    errorbar(x,y,E)
    i   = p(1,1);   j=-1;   k=p(1,2);   m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
    ind =[];        ind= find((abs(m)>(mean(E))));%&(rgs>500));
    plot (K67.netrad(ind),y1(ind),'.r');
    
    record = [record;ind];
end
% % % ind = find((~isnan(K67.netrad))&(~isnan(y1)));
% % % y2  = polyfit(K67.netrad(ind),y1(ind),1);
% % % x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
% % % y   = polyval(y2,x);    plot (x,y,'-','LineWidth',2,'Color','k');
% % % xlabel('Rn');           ylabel('LE');
% % % E   = std(y).*3.*ones(size(x)); %error bars that are two standard deviation units
% % % errorbar(x,y,E)
% % % i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
% % % ind = find((abs(m)>(mean(E))));%&(rgs>500));
% % % plot (K67.netrad(ind),y1(ind),'.m');

subplot(2,2,2);
y1 = K67.H;
[p,r,x]=AM_compute_regression(K67.netrad,y1,[.8 .8 .8],1);         hold on;
ind = find((~isnan(K67.netrad))&(~isnan(y1)));
y2  = polyfit(K67.netrad(ind),y1(ind),1);
x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
y   = polyval(y2,x);    plot (x,y,'-','LineWidth',2,'Color','k');
ylabel('Rn');           xlabel('H');
E   = std(y).*3.*ones(size(x)); %eRnor bars that are two standard deviation units
errorbar(x,y,E)
i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
ind = find((abs(m)>(mean(E))));%&(rgs>500));
plot (K67.netrad(ind),y1(ind),'.m');
% K67.H(ind)  = NaN;

y1 = K67.H + K67.LE;
subplot(2,2,3);
[p,r,x]=AM_compute_regression(K67.netrad,y1,[.8 .8 .8],1);         hold on;
ind = find((~isnan(K67.netrad))&(~isnan(y1)));
y2  = polyfit(K67.netrad(ind),y1(ind),1);
x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
y   = polyval(y2,x);    plot (x,y,'-','LineWidth',2,'Color','k');
ylabel('Rn');           xlabel('H+LE');
E   = std(y).*2.*ones(size(x)); %eRnor bars that are two standard deviation units
errorbar(x,y,E)
i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
ind = find((abs(m)>(mean(E))));
plot (K67.netrad(ind),y1(ind),'.r');
ind = find(abs(K67.netrad-K67.LE-K67.H)>400);%&(rgs>500));
plot (K67.netrad(ind),y1(ind),'.r');
K67.H(ind)  = NaN;          K67.netrad(ind) = NaN;
K67.LE(ind) = NaN;          fco2(ind)=NaN;

subplot(2,2,4);
[p,r,x]=AM_compute_regression(K67.par,y1,[.8 .8 .8],1);         hold on;
ind = find((K67.dateloc>datenum(2007,1,1))&(~isnan(K67.par))&(~isnan(y1)));
y2  = polyfit(K67.par(ind),y1(ind),2);
x   = nanmin(K67.par):((nanmax(K67.par)-nanmin(K67.par))/20):nanmax(K67.par);
y   = polyval(y2,x);     plot (x,y,'-','LineWidth',2,'Color','k');
ylabel('par');           xlabel('H+LE');
E   = 2.*std(y).*ones(size(x)); %eRnor bars that are two standard deviation units
errorbar(x,y,E)
i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.par.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
ind = find((abs(m)>(mean(E))));
plot (K67.par(ind),y1(ind),'.m');
% K67.H(ind)=NaN;          K67.par(ind)=NaN;     K67.LE(ind)=NaN;

figure('color','white');                plot(K67.jdstartgmt,K67.netrad,'-r');
hold on;    plot(K67.jdstartgmt,K67.H,'-m');
plot(K67.jdstartgmt,K67.LE,'-b');
plot(K67.jdstartgmt(record),K67.netrad(record),'.g','MarkerSize',20);

% ..........................................................................
%% Aditional variables
% ..........................................................................
%  CO2 storage
% ..........................................................................
K67.pambpa_fill = AM_fill_month_cycle(K67.pambpa,K67.dateloc,6,18);
K67.pambKpa     = K67.pambpa_fill./1000;
K67.tc_fill     = AM_fill_day_cycle(K67.tasp,K67.dateloc,6,18);
K67.RH_fill     = K67.RHec;                 K67.RH_fill(isnan(K67.RHec)) = K67.RH(isnan(K67.RHec))./100;
K67.RH_fill(K67.RH_fill>1.05) = 1.01;       K67.RH_fill = AM_fill_day_cycle(K67.RH_fill,K67.dateloc,6,18);
K67.rho_a = (K67.pambKpa./(0.28705*(K67.tc_fill+273.15))).*(1-(0.378*K67.e_sat.*K67.RH_fill./K67.pambKpa));
%Phillip % density of air, kg/m3
%Density of air 1.2041 kg/m-3 at 20degC

[level_CO2, height_CO2]   =   textread(['/home/ncoupe/tapajos/process.u.az/heights.txt'],...
    '%s	%f','headerlines', 1);

DataTable  = [K67.co2pf1,K67.co2pf2,K67.co2pf3,K67.co2pf4,K67.co2pf5,K67.co2pf6,K67.co2pf7,K67.co2pf8];
K67.storageCO2 = AM_Sco2_calc(K67.dateloc,K67.rho_a,DataTable,height_CO2,HeightEC,'CO2');

DataTableTair  = [K67.tair1_fill,K67.tair2_fill,K67.tair3_fill,K67.tair4_fill,...
    K67.tair5_fill,K67.tair6_fill,K67.tair7_fill,K67.tair8_fill];
K67.StorageHeat = AM_Sco2_calc(K67.dateloc,K67.rho_a,DataTableTair,height_CO2,HeightEC,'H');
K67.StorageHeat((K67.dateloc>datenum(2011,12,18))&(K67.dateloc<datenum(2013,12,31))) = NaN;   %few tair profile working

K67.StorageHeatCanopy = AM_Sco2_calc(K67.dateloc,K67.rho_a,DataTableTair(:,3:end),height_CO2(3:end),HeightEC,'H');
K67.StorageHeatCanopy((K67.dateloc>datenum(2011,12,18))&(K67.dateloc<datenum(2013,12,31))) = NaN;   %few tair profile working

DataTable = [K67.h2opf1,K67.h2opf2,K67.h2opf3,K67.h2opf4,K67.h2opf5,K67.h2opf6,K67.h2opf7,K67.h2opf8];
DataTable = DataTable.* 1e-3.*18.016./0.02895;              % mixing ratio g kg-1
DataTable = (DataTable.*(K67.pambpa)./(DataTable+622));     % vapour pressure in Pa
DataTable = DataTable.*2.16679./(K67.tc_fill+273.15);     % absolute humidity g m-3  AH = C × (Pw/TK), where C = 2.16679 gk/J

K67.tasp_fill = AM_fill(K67.tasp);
K67.tasp_fill = AM_fill_day_cycle(K67.tasp,K67.dateloc,6,18);

K67.lambda = ((2.50*10^6)-((2.386*10^3).*(K67.tasp_fill)))./10^6;                %2.45 [MJ kg-1]     % Latent heat of vaporization, 2.45 [MJ kg-1],
K67.lambda = AM_fill_day_cycle(K67.lambda,K67.dateloc,6,18);

K67.StorageLE = AM_Sco2_calc(K67.dateloc,K67.lambda,DataTable,height_CO2,HeightEC,'LE');
K67.StorageLE = AM_rm_outlier(K67.StorageLE,4);
K67.StorageLE = AM_fill_plus(K67.StorageLE);

dummy = AM_Sco2_calc(K67.dateloc,K67.lambda,DataTable(:,2),35,HeightEC,'LE');

Cveg = 2958;   % J kg-1 K-1 Medorow
%         K67 325 Mg ha−1  Hunter
%         mveg = 325.*1000./10000;    % kg m-2
mveg = 325.*1000.*2./10000;    % kg m-2   x2 from dry to wet
if ip==9;  LevelCpyTemp = 2;  end
iy = [NaN(1,1);diff(K67.tair3_fill)];      %iy(end) = [];
iy(iy==0) = NaN;
iy = AM_rm_outlier(iy,3);                             iy = AM_fill(iy);
K67.BiomassHeat = Cveg.*mveg.*iy./(60.*60);

%% .........................................................................
figure('color','white');                    hold on;
plot(K67.dateloc,K67.netrad,'k');
plot(K67.dateloc,K67.H,'r');
plot(K67.dateloc,K67.LE,'b');
plot(K67.dateloc,K67.StorageLE,'c');        %plot(K67.dateloc,dummy,'g');
plot(K67.dateloc,K67.BiomassHeat,'color',[0 .4 0]);
plot(K67.dateloc,K67.StorageHeat,'color',[.7 .7 0]);      ylabel('\Delta Storage heat');
datetick('x');

%% .........................................................................
% Energy for phothosyntesis
% 2.5mgCO2 m-2 s-2 ~ 28 W m-2 Masseroni, Corbari and Mancini, 2014
K67.Ephotosyntesis = -K67.gee.*1000.*28./(2.5.*60.*60.*24);
K67.Ephotosyntesis_light = -K67.geefilllight.*1000.*28./(2.5.*60.*60.*24);

%% .........................................................................
color_mtx = [0 0 0; 171/255 221/255 164/255;43/255 131/255 186/255;...
    215/255  25/255 28/255; 253/255 174/255 97/255; ...
    0.020 0.0549 0.6882; 0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

figure('color','white');
subplot(2,2,4);              hold on;
% x = AM_cycle(K67.So./2,K67.dateloc);   plot((0:1/23:1),x,'color',color_mtx(1,:),'LineWidth',2);
y = AM_cycle(K67.netrad,K67.dateloc);   plot((0:1/23:1),y./y,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE,K67.dateloc);       plot((0:1/23:1),x./y,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H,K67.dateloc);        plot((0:1/23:1),x./y,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat,K67.dateloc); 	plot((0:1/23:1),x./y,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE,K67.dateloc);    plot((0:1/23:1),x./y,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis,K67.dateloc);    plot((0:1/23:1),x./y,'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -.400 1]);     grid on;
legend('Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}','FontSize',12);
legend box off;
datetick('x','HH:MM');      box on;          ylabel({'Flux and Storage';'Normalized by Rn'},'FontSize',12);

subplot(2,2,2);                         hold on;
x = AM_cycle(K67.So./2,K67.dateloc);    plot((0:1/23:1),x,'color',color_mtx(1,:),'LineWidth',2);
x = AM_cycle(K67.netrad,K67.dateloc);   plot((0:1/23:1),x,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE,K67.dateloc);       plot((0:1/23:1),x,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H,K67.dateloc);        plot((0:1/23:1),x,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat,K67.dateloc); 	plot((0:1/23:1),x,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE,K67.dateloc);    plot((0:1/23:1),x,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis,K67.dateloc);    plot((0:1/23:1),x,'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -100 700]);     grid on;
legend('0.5TOA','Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}');
legend box off;
datetick('x','HH:MM');      box on;          ylabel({'Energy Flux and';'Storage (W m^-^2)'},'FontSize',12);

%%
figure('color','white');
subplot(2,2,1);              hold on;
ind = find((K67.H~=0));
y = histogram(K67.H(ind));          box on;    xlabel ('H (W m^-^2)');
axis([-100 250 0 12000]);     grid on;    
subplot(2,2,2);              hold on;
ind = find((K67.So==0)&(K67.H~=0));
y = histogram(K67.H(ind));          box on; xlabel ('H (W m^-^2)');
axis([-100 20 0 3000]);     grid on;    title('Night-time frequency distribution');
% datetick('x','HH:MM');      box on;          ylabel('Flux and Storage Normalized by Rn','FontSize',12);
subplot(2,2,4);              hold on;
ind = find((K67.So>=10)&(K67.H~=0));
y = histogram(K67.H(ind));   hold on;       
ind = find((K67.So>=400)&(K67.H~=0));
y = histogram(K67.H(ind));          box on; xlabel ('H (W m^-^2)');
axis([-100 250 0 3000]);     grid on;  title('Day-time frequency distribution');
legend ('TOA>10 W m^-^2','TOA>400 W m^-^2');  legend box off;


%% .........................................................................
figure('color','white');
subplot(2,2,1);              hold on;
x = AM_cycle(K67.So./2,K67.dateloc);   plot((0:1/23:1),x,'color',color_mtx(1,:),'LineWidth',2);
x = AM_cycle(K67.netrad,K67.dateloc);   plot((0:1/23:1),x,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE,K67.dateloc);       plot((0:1/23:1),x,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H,K67.dateloc);        plot((0:1/23:1),x,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat,K67.dateloc); 	plot((0:1/23:1),x,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE,K67.dateloc);    plot((0:1/23:1),x,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis,K67.dateloc);    plot((0:1/23:1),x,'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -100 700]);     grid on;
legend('0.5TOA','netrad','LE','H','stgH','stgLE','E_{photosynthesis}');
legend box off;
datetick('x','HH:MM');      box on;          ylabel('(W m^-^2)');

color_mtx = jet(8);
subplot(2,2,3);              hold on;
ind = find((K67.dateloc>=datenum(2019,12,11))&(K67.dateloc<datenum(2019,12,12)));
plot((0:1/23:1),K67.So(ind)./2,'color',color_mtx(1,:),'LineWidth',2);
plot((0:1/23:1),K67.netrad(ind),'color',color_mtx(2,:),'LineWidth',2);
plot((0:1/23:1),K67.LE(ind),'color',color_mtx(3,:),'LineWidth',2);
plot((0:1/23:1),K67.H(ind),'color',color_mtx(4,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageHeat(ind),'color',color_mtx(5,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageLE(ind),'color',color_mtx(6,:),'LineWidth',2);
plot((0:1/23:1),K67.Ephotosyntesis(ind),'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -100 700]);     grid on;
datetick('x','HH:MM');      box on;          ylabel('(W m^-^2)');
title(datestr(datenum(2019,12,11)));

subplot(2,2,4);              hold on;
ind = find((K67.dateloc>=datenum(2017,5,10))&(K67.dateloc<datenum(2017,5,11)));
plot((0:1/23:1),K67.So(ind)./2,'color',color_mtx(1,:),'LineWidth',2);
plot((0:1/23:1),K67.netrad(ind),'color',color_mtx(2,:),'LineWidth',2);
plot((0:1/23:1),K67.LE(ind),'color',color_mtx(3,:),'LineWidth',2);
plot((0:1/23:1),K67.H(ind),'color',color_mtx(4,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageHeat(ind),'color',color_mtx(5,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageLE(ind),'color',color_mtx(6,:),'LineWidth',2);
plot((0:1/23:1),K67.Ephotosyntesis(ind),'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -100 700]);     grid on;
datetick('x','HH:MM');      box on;          ylabel('(W m^-^2)');
title(datestr(datenum(2017,5,10)));

%% .........................................................................
figure('color','white');
subplot(2,2,1);              hold on;
% x = AM_cycle(K67.So./2,K67.dateloc);   plot((0:1/23:1),x,'color',color_mtx(1,:),'LineWidth',2);
y = AM_cycle(K67.netrad,K67.dateloc);   plot((0:1/23:1),y./y,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE,K67.dateloc);       plot((0:1/23:1),x./y,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H,K67.dateloc);        plot((0:1/23:1),x./y,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat,K67.dateloc); 	plot((0:1/23:1),x./y,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE,K67.dateloc);    plot((0:1/23:1),x./y,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis,K67.dateloc);    plot((0:1/23:1),x./y,'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -.400 1]);     grid on;
legend('Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}');
legend box off;
datetick('x','HH:MM');      box on;          ylabel('Flux and Storage Normalized by Rn');

subplot(2,2,3);              hold on;
ind = find((K67.dateloc>=datenum(2019,12,11))&(K67.dateloc<datenum(2019,12,12)));
% plot((0:1/23:1),K67.So(ind)./2,'color',color_mtx(1,:),'LineWidth',2);
plot((0:1/23:1),K67.netrad(ind)./K67.netrad(ind),'color',color_mtx(2,:),'LineWidth',2);
plot((0:1/23:1),K67.LE(ind)./K67.netrad(ind),'color',color_mtx(3,:),'LineWidth',2);
plot((0:1/23:1),K67.H(ind)./K67.netrad(ind),'color',color_mtx(4,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageHeat(ind)./K67.netrad(ind),'color',color_mtx(5,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageLE(ind)./K67.netrad(ind),'color',color_mtx(6,:),'LineWidth',2);
plot((0:1/23:1),K67.Ephotosyntesis(ind)./K67.netrad(ind),'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -.600 1.2]);     grid on;
datetick('x','HH:MM');      box on;          ylabel('normalized by netrad');
title(datestr(datenum(2019,12,11)));

subplot(2,2,4);              hold on;
ind = find((K67.dateloc>=datenum(2017,5,10))&(K67.dateloc<datenum(2017,5,11)));
plot((0:1/23:1),K67.netrad(ind)./K67.netrad(ind),'color',color_mtx(2,:),'LineWidth',2);
plot((0:1/23:1),K67.LE(ind)./K67.netrad(ind),'color',color_mtx(3,:),'LineWidth',2);
plot((0:1/23:1),K67.H(ind)./K67.netrad(ind),'color',color_mtx(4,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageHeat(ind)./K67.netrad(ind),'color',color_mtx(5,:),'LineWidth',2);
plot((0:1/23:1),K67.StorageLE(ind)./K67.netrad(ind),'color',color_mtx(6,:),'LineWidth',2);
plot((0:1/23:1),K67.Ephotosyntesis(ind)./K67.netrad(ind),'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -.600 1.2]);     grid on;
datetick('x','HH:MM');      box on;         ylabel('normalized by netrad');
title(datestr(datenum(2017,5,10)));

%% ........................................................................
color_mtx = jet(8);
figure ('color','white');
subplot(2,2,1);              hold on;
ind = find((K67.dateloc>=datenum(2019,12,11))&(K67.dateloc<datenum(2019,12,12)));
plot((0:1/23:1),K67.tair1(ind),'color',color_mtx(1,:),'LineWidth',2);
plot((0:1/23:1),K67.tair2(ind),'color',color_mtx(2,:),'LineWidth',2);
plot((0:1/23:1),K67.tair3(ind),'color',color_mtx(3,:),'LineWidth',2);
plot((0:1/23:1),K67.tair4(ind),'color',color_mtx(4,:),'LineWidth',2);
plot((0:1/23:1),K67.tair5(ind),'color',color_mtx(5,:),'LineWidth',2);
plot((0:1/23:1),K67.tair6(ind),'color',color_mtx(6,:),'LineWidth',2);
plot((0:1/23:1),K67.tair7(ind),'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 22 32]);     grid on;
datetick('x','HH:MM');  box on;         ylabel('T_{air}(\circC)');
legend('62m','50m','39m','29m','20m','10m','3m');
title(datestr(datenum(2019,12,11)));     legend box off;

[K67.weekLE,~,~,~,~,~,K67.week] = AM_week2cycle_rs(K67.LE,K67.dateloc,tresholdCycleDay,3,2);
LAI_week = weekNPPsite(ip).LAI;
LAI_week_avg = AM_week2_avg(LAI_week,K67.week(1:length(LAI_week)));
[LAI_week] = AM_week2dayAVG_recurrent(LAI_week_avg,K67.week);
LAI_dateloc = interp1(K67.week,LAI_week,K67.dateloc);

% .........................................................................
%% Leaf temperature
% .........................................................................
canopy   = canopy_BF(ip);
extCoeff = 4.39-(3.97.*exp(-0.258.*LAI_dateloc));
zeta     = canopy./height;
Uzeta    = K67.ws;

Uh = Uzeta./exp((extCoeff.*zeta)-1);

Leaf = sqrt((182.25-20.25)./2);   % cm

Leaf = Leaf./100;                 % cm to m
Leaf = 0.1;
one = (Leaf./Uh);           two = extCoeff.*(1-zeta)./2;
rbH = (150./LAI_dateloc).*sqrt(one.*exp(two).*1);
Sc  = 0.67;
Pr  = 0.71;
f   = 0.5;
rbV = (1./f).*rbH.*((Sc./Pr).^(2/3));            % nanmean == 31
CpJ = 1.013*10^-3.*10^6;  %J kg-1 C-1
K67.TL  = ((K67.H.*rbH)./(K67.rho_a.*CpJ))+K67.tasp;
plot((0:1/23:1),AM_fill(K67.TL(ind)),'color',[0 0 0],'LineWidth',2);

legend('62m','50m','39m','29m','20m','10m','3m','LeafTemp_{W&S}');

%%
x1 = (-200:10:1000);
figure('color','white');
subplot(2,2,1);                       hold on;
[p,r2,~,~,sm] = AM_compute_regression_zero(K67.netrad,K67.H+K67.LE,[0.4,0.4,0.4],1);    p = [p,0];
y_model = polyval(p,K67.netrad);                     hold on;
plot(x1,x1,':','Color',[.5 .5 .5]);
plot(K67.netrad,y_model,'-','Color',[.7 0 0]);
text (-50,200,{['LE+H=' num2str(p(1),'%4.2f') 'Rn'];...% num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color','k');
grid on;    box on;     axis([-200 1000 -200 1000]);
ylabel('LE+H (W m^-^2)','FontSize',12);     xlabel('Rn (W m^-^2)','FontSize',12);

subplot(2,2,3);                       hold on;
[p,r2,~,~,sm] = AM_compute_regression_zero(K67.netrad-K67.StorageHeat-K67.StorageLE-K67.Ephotosyntesis,K67.H+K67.LE,...
[0.4,0.4,0.4],1);       p = [p,0];
y_model = polyval(p,K67.netrad);                     hold on;
plot(K67.netrad,y_model,'-','Color',[.7 0 0]);
plot(x1,x1,':','Color',[.5 .5 .5]);
text (-50,200,{['LE+H=' num2str(p(1),'%4.2f') '(Rn-\Delta S_{H}-\Delta S_{LE}-E_{p})'];...% num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color','k');
grid on;    box on;     axis([-200 1000 -200 1000]);
ylabel('LE+H (W m^-^2)','FontSize',12);     xlabel('Rn-\Delta S_{H}-\Delta S_{LE}-E_{p}(W m^-^2)','FontSize',12);

% % ind = (K67.So>200);
% % subplot(2,2,4);                       hold on;
% % [p,r2,~,~,sm] = AM_compute_regression_zero(K67.netrad(ind)-K67.StorageHeat(ind)-K67.StorageLE(ind)-K67.Ephotosyntesis(ind),...
% %     K67.H(ind)+K67.LE(ind),[0.4,0.4,0.4],1);    p = [p,0];
% % % [p,r2,~,~,sm] = AM_compute_regression(K67.netrad,K67.H+K67.LE+K67.StorageLE+K67.StorageHeat,[0.4,0.4,0.4],1);
% % y_model = polyval(p,K67.netrad(ind)-K67.StorageHeat(ind)-K67.StorageLE(ind)-K67.Ephotosyntesis(ind));                     hold on;
% % plot(K67.netrad(ind)-K67.StorageHeat(ind)-K67.StorageLE(ind)-K67.Ephotosyntesis(ind),y_model,'-','Color',[.7 0 0]);
% % text (-50,200,{['Daytime LE+H=' num2str(p(1),'%4.2f') 'Rn-\Delta S_{H}-\Delta S_{LE}-E_{p}'];... % num2str(p(2),'%4.2f')];...
% %     ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color','k');
% % grid on;    box on;
% % ylabel('LE+H (W m^-^2)','FontSize',12);     xlabel('Rn-\Delta S_{H}-\Delta S_{LE}-E_{p}(W m^-^2)','FontSize',12);

% ind = (K67.So>200);
[K67.dayNetRad,~,x]      = AM_day(K67.netrad,K67.dateloc,0);        K67.dayNetRad(x<21) = NaN; 
[K67.dayStorageHeat,~,x] = AM_day(K67.StorageHeat,K67.dateloc,0);   K67.dayStorageHeat(x<21) = NaN; 
[K67.dayStorageLE,~,x]   = AM_day(K67.StorageLE,K67.dateloc,0);     K67.dayStorageLE(x<21) = NaN; 
[K67.dayEphotosyntesis,~,x] = AM_day(K67.Ephotosyntesis,K67.dateloc,0);     K67.dayEphotosyntesis(x<21) = NaN; 
[K67.dayH,~,x] = AM_day(K67.H,K67.dateloc,0);       K67.dayH(x<21) = NaN; 
[K67.dayLE,~,x] = AM_day(K67.LE,K67.dateloc,0);     K67.dayLE(x<21) = NaN; 
subplot(2,2,4);                       hold on;
[p,r2,~,~,sm] = AM_compute_regression_zero(K67.dayNetRad-K67.dayStorageHeat-K67.dayStorageLE-K67.dayEphotosyntesis,...
    K67.dayH+K67.dayLE,[0.4,0.4,0.4],1);    p = [p,0];
y_model = polyval(p,K67.dayNetRad-K67.dayStorageHeat-K67.dayStorageLE-K67.dayEphotosyntesis);   hold on;
plot(K67.dayNetRad-K67.dayStorageHeat-K67.dayStorageLE-K67.dayEphotosyntesis,y_model,'-','Color',[.7 0 0]);
plot(x1,x1,':','Color',[.5 .5 .5]);
text (-50,200,{['Daily LE+H=' num2str(p(1),'%4.2f') '(Rn-\Delta S_{H}-\Delta S_{LE}-E_{p})'];... % num2str(p(2),'%4.2f')];...
    ['r^2:' num2str(r2(1),'%3.2f')]},'FontSize',10,'color','k');
grid on;    box on;     axis([0 250 0 250]);
ylabel('LE+H (W m^-^2)','FontSize',12);     xlabel('Rn-\Delta S_{H}-\Delta S_{LE}-E_{p}(W m^-^2)','FontSize',12);

%% ........................................................................
tresholdEBclosue = 100;
ind = find(((K67.netrad-K67.LE-K67.H-K67.StorageHeat-K67.StorageLE-K67.Ephotosyntesis)./K67.netrad)>tresholdEBclosue);

figure('color','white');
subplot(2,2,1);
y1 = K67.H+K67.LE;
[p,r,x]=AM_compute_regression(K67.netrad,y1,[.8 .8 .8],1);         hold on;
plot (K67.netrad(ind),y1(ind),'.r');
ind = find((~isnan(K67.netrad))&(~isnan(y1)));
% y2  = polyfit(K67.netrad(ind),y1(ind),1);
x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
y   = polyval(p,x);         plot (x,y,'-','LineWidth',2,'Color','k');
ylabel('Rn');               xlabel('H+LE');
E   = 2.*std(y).*ones(size(x));
errorbar(x,y,E)
i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
ind = find((abs(m)>(mean(E))));
plot (K67.netrad(ind),y1(ind),'.m');
%K67.H(ind)  = NaN;

subplot(2,2,2);
y1 = K67.LE;
[p,r,x]=AM_compute_regression(K67.netrad,y1,[.8 .8 .8],1);         hold on;
ind = find((~isnan(K67.netrad))&(~isnan(y1)));
% y2  = polyfit(K67.netrad(ind),y1(ind),1);
x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
y   = polyval(p,x);        plot (x,y,'-','LineWidth',2,'Color','k');
ylabel('Rn');               xlabel('LE');
E   = 2.*std(y).*ones(size(x));
errorbar(x,y,E)
i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
ind = find((abs(m)>(mean(E))));
plot (K67.netrad(ind),y1(ind),'.m');

subplot(2,2,3);
y1 = K67.H;
[p,r,x]=AM_compute_regression(K67.netrad,y1,[.8 .8 .8],1);         hold on;
ind = find((~isnan(K67.netrad))&(~isnan(y1)));
% y2  = polyfit(K67.netrad(ind),y1(ind),1);
x   = nanmin(K67.netrad):((nanmax(K67.netrad)-nanmin(K67.netrad))/20):nanmax(K67.netrad);
y   = polyval(p,x);            plot (x,y,'-','LineWidth',2,'Color','k');
ylabel('Rn');                   xlabel('H');
E   = 2.*std(y).*ones(size(x));
errorbar(x,y,E)
i   = p(1,1); j=-1; k=p(1,2); m=abs(((K67.netrad.*i)+(y1.*j)+k)./(sqrt((i.^2)+1))); %distance point to line
ind = find((abs(m)>(mean(E))));
plot (K67.netrad(ind),y1(ind),'.m');

% .........................................................................
%% Net Rad check
% .........................................................................
run('/home/ncoupe/Documents/AmazonRSsite/ICI/AM_02ici_logger_noresample.m');
close all

% % % %% ........................................................................
% % % [K67.weekNetRad,~,~,~,K67.week_maxNetRad] = AM_week2cycle_rs(K67.netrad,K67.dateloc,tresholdCycle,3,2);
% % % figure('color','white');
% % % subplot(2,2,1:2);       plot(K67.dateloc,K67.netrad,'LineWidth',2);
% % % hold on;                plot(K67cameraMin.minutes,K67cameraMin.Rn_Avg);
% % % legend ('K67','camera'); datetick('x');
% % % ix = find(K67.dateloc>=K67camera.dateloc(1));           iy = find(K67camera.dateloc<=K67.dateloc(end));
% % % subplot(2,2,3);
% % % [ip,r2] = AM_compute_regression(K67camera.Rn_Avg(iy),K67.netrad(ix),[.4 .4 .4],1);
% % % 
% % % %%
% % % figure('color','white');
% % % subplot(2,2,1:2);       plot(K67.dateloc,K67.tc,'LineWidth',2);
% % % hold on;                plot(K67cameraMin.minutes,K67cameraMin.AirTC_Avg);
% % % legend ('K67','camera');datetick('x');
% % % ix = find(K67.dateloc>=K67camera.dateloc(1));           iy = find(K67camera.dateloc<=K67.dateloc(end));
% % % subplot(2,2,3);
% % % [ip,r2] = AM_compute_regression(K67camera.AirTC_Avg(iy),K67.tc(ix),[.4 .4 .4],1);
% % % %%
% % % 
% % % figure('color','white');
% % % subplot(2,2,1:2);      hold on;         plot(K67.dateloc(ix),K67.tc(ix),'r');
% % % % plot(K67.dateloc(ix),(K67.tc(ix).*-78.854)-1517.127,'k');
% % % plot(K67camera.dateloc(iy),K67camera.AirTC_Avg(iy),'b');
% % % plot(K67camera.dateloc(iy),K67camera.TempTC_C_Avg2(iy),'g');
% % % plot(K67camera.dateloc(iy),K67camera.TempTC_C_Avg1(iy),'c');
% % % plot(K67.dateloc(ix),K67.tair1(ix),'m');
% % % ylabel('t_{air}');                      datetick('x');
% % % legend ('K67','camera HMP','termocouple ICI 1','ICI2','profile 1');
% % % 
% % % %%
% % % % % % % iy = K67cameraMin.AirTC_Avg;        K67camera.AirTC_Avg = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.RH;               K67camera.RH = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.RH_Max;           K67camera.RH_Max = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.RH_Min;           K67camera.RH_Min = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.PAR_Den_Avg;      K67camera.PAR_Den_Avg = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.PAR_Tot_Tot;      K67camera.PAR_Tot_Tot = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.TempTC_C_Avg1;    K67camera.TempTC_C_Avg1 = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.TempTC_C_Avg2;    K67camera.TempTC_C_Avg2 = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.AirT_therm_Avg;   K67camera.AirT_therm_Avg = AM_hour(iy,ix);       
% % % % % % % iy = K67cameraMin.CNR4TC_Avg;       K67camera.CNR4TC_Avg = AM_hour(iy,ix);       


% .........................................................................
%%    Calculate VPD                                                                      #
% ........................................................................
% % tasp.fill should be the temperature at the level in degrees
% % pamb.pa: is ambient pressure in Pascals
% % the h20.mix needs to be in gh20/kg dry air.... if it is in mmol/mol use this to run the conversion
% % h2o.mix      <- ed1.halfhr[,"h2o.mmol.m"] * 1e-3*18.016/0.02895 # (g h2o/ kg dry air)
% % 
A0 = 6.107799961;       A1 = 4.436518521e-1;        A2 = 1.428945805e-2;
A3 = 2.650648471e-4;    A4 = 3.031240396e-6;        A5 = 2.034080948e-8;
A6 = 6.136820929e-11;  
R  = 0.28704;            %R specific gas constant = 0.287 kJ kg-1 K-1

K67.h2omixpf1 = K67.h2opf1.* 1e-3.*18.016./0.02895; K67.h2omixpf2 = K67.h2opf2.* 1e-3.*18.016./0.02895;
K67.h2omixpf3 = K67.h2opf3.* 1e-3.*18.016./0.02895; K67.h2omixpf4 = K67.h2opf4.* 1e-3.*18.016./0.02895;
K67.h2omixpf5 = K67.h2opf5.* 1e-3.*18.016./0.02895; K67.h2omixpf6 = K67.h2opf6.* 1e-3.*18.016./0.02895;
K67.h2omixpf7 = K67.h2opf7.* 1e-3.*18.016./0.02895; K67.h2omixpf8 = K67.h2opf8.* 1e-3.*18.016./0.02895;

e_sat = 0.1.*(A0+K67.tair1_fill.*(A1+K67.tair1_fill.*(A2+K67.tair1_fill.*(A3+...
K67.tair1_fill.*(A4+K67.tair1_fill.*(A5+K67.tair1_fill.*A6))))));
e_air = (K67.h2omixpf1.*((K67.pambpa./1000)./100)./(K67.h2omixpf1+622)).*10^2;
K67.RHpf1    = e_air./e_sat;
K67.VPDpf1   = e_sat.*(1-K67.RHpf1);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair2_fill.*(A1+K67.tair2_fill.*(A2+K67.tair2_fill.*(A3+...
    K67.tair2_fill.*(A4+K67.tair2_fill.*(A5+K67.tair2_fill.*A6))))));
e_air = (K67.h2omixpf2.*((K67.pambpa./1000)./100)./(K67.h2omixpf2+622)).*10^2;
K67.RHpf2    = e_air./e_sat;
K67.VPDpf2   = e_sat.*(1-K67.RHpf2);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair3_fill.*(A1+K67.tair3_fill.*(A2+K67.tair3_fill.*(A3+...
    K67.tair3_fill.*(A4+K67.tair3_fill.*(A5+K67.tair3_fill.*A6))))));
e_air = (K67.h2omixpf3.*((K67.pambpa./1000)./100)./(K67.h2omixpf3+622)).*10^2;
K67.RHpf3    = e_air./e_sat;
K67.VPDpf3   = e_sat.*(1-K67.RHpf3);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair4_fill.*(A1+K67.tair4_fill.*(A2+K67.tair4_fill.*(A3+...
    K67.tair4_fill.*(A4+K67.tair4_fill.*(A5+K67.tair4_fill.*A6))))));
e_air = (K67.h2omixpf4.*((K67.pambpa./1000)./100)./(K67.h2omixpf4+622)).*10^2;
K67.RHpf4    = e_air./e_sat;
K67.VPDpf4   = e_sat.*(1-K67.RHpf4);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair5_fill.*(A1+K67.tair5_fill.*(A2+K67.tair5_fill.*(A3+...
    K67.tair5_fill.*(A4+K67.tair5_fill.*(A5+K67.tair5_fill.*A6))))));
e_air = (K67.h2omixpf5.*((K67.pambpa./1000)./100)./(K67.h2omixpf5+622)).*10^2;
K67.RHpf5    = e_air./e_sat;
K67.VPDpf5   = e_sat.*(1-K67.RHpf5);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair6_fill.*(A1+K67.tair6_fill.*(A2+K67.tair6_fill.*(A3+...
    K67.tair6_fill.*(A4+K67.tair6_fill.*(A5+K67.tair6_fill.*A6))))));
e_air = (K67.h2omixpf6.*((K67.pambpa./1000)./100)./(K67.h2omixpf6+622)).*10^2;
K67.RHpf6    = e_air./e_sat;
K67.VPDpf6   = e_sat.*(1-K67.RHpf6);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair7_fill.*(A1+K67.tair7_fill.*(A2+K67.tair7_fill.*(A3+...
    K67.tair7_fill.*(A4+K67.tair7_fill.*(A5+K67.tair7_fill.*A6))))));
e_air = (K67.h2omixpf7.*((K67.pambpa./1000)./100)./(K67.h2omixpf7+622)).*10^2;
K67.RHpf7    = e_air./e_sat;
K67.VPDpf7   = e_sat.*(1-K67.RHpf7);                                       % vapour pressure deficit;

e_sat = 0.1.*(A0+K67.tair8_fill.*(A1+K67.tair8_fill.*(A2+K67.tair8_fill.*(A3+...
    K67.tair8_fill.*(A4+K67.tair8_fill.*(A5+K67.tair8_fill.*A6))))));
e_air = (K67.h2omixpf8.*((K67.pambpa./1000)./100)./(K67.h2omixpf8+622)).*10^2;
K67.RHpf8    = e_air./e_sat;
K67.VPDpf8   = e_sat.*(1-K67.RHpf8);                                       % vapour pressure deficit;

%%
figure('color','white');    hold on;
plot(K67.jdstartgmt,K67.VPDpf1,'color',color_jet(1,:));  plot(K67.jdstartgmt,K67.VPDpf2,'color','k');%color_jet(2,:)
plot(K67.jdstartgmt,K67.VPDpf3,'color',color_jet(3,:));  plot(K67.jdstartgmt,K67.VPDpf4,'color',color_jet(4,:));
plot(K67.jdstartgmt,K67.VPDpf5,'color',color_jet(5,:));  plot(K67.jdstartgmt,K67.VPDpf6,'color',color_jet(6,:));
plot(K67.jdstartgmt,K67.VPDpf7,'color',color_jet(7,:));  plot(K67.jdstartgmt,K67.VPDpf8,'color',color_jet(8,:));
plot(K67.jdstartgmt,K67.VPDec,'color','k');
legend('1','2','3','4','5','6','7','8','VPD_{ec}');


%% .........................................................................
color_mtx = [0 0 0; 171/255 221/255 164/255;43/255 131/255 186/255;...
    215/255  25/255 28/255; 253/255 174/255 97/255; ...
    0.020 0.0549 0.6882; 0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

figure('color','white');
subplot(2,2,3);              hold on;
ind = find((K67.dateloc>=datenum(2020,6,15))&(K67.dateloc<datenum(2020,6,16)));
y = AM_cycle(K67.netrad(ind),K67.dateloc(ind));   plot(K67.dateloc(ind),y./y,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE(ind),K67.dateloc(ind));       plot(K67.dateloc(ind),x./y,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H(ind),K67.dateloc(ind));        plot(K67.dateloc(ind),x./y,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat(ind),K67.dateloc(ind)); 	plot(K67.dateloc(ind),x./y,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE(ind),K67.dateloc(ind));      plot(K67.dateloc(ind),x./y,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis(ind),K67.dateloc(ind)); plot(K67.dateloc(ind),x./y,'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -.400 1]);     grid on;
legend('Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}','FontSize',12);
legend box off;             xtickangle(0);
datetick('x','HH:MM');      box on;          ylabel({'Flux and Storage';'Normalized by Rn'},'FontSize',12);
title ('2020/06/15'); 
subplot(2,2,1);                         hold on;
x = AM_cycle(K67.So(ind)./2,K67.dateloc(ind));    plot(K67.dateloc(ind),x,'color',color_mtx(1,:),'LineWidth',2);
x = AM_cycle(K67.netrad(ind),K67.dateloc(ind));   plot(K67.dateloc(ind),x,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE(ind),K67.dateloc(ind));       plot(K67.dateloc(ind),x,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H(ind),K67.dateloc(ind));        plot(K67.dateloc(ind),x,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat(ind),K67.dateloc(ind)); 	plot(K67.dateloc(ind),x,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE(ind),K67.dateloc(ind));      plot(K67.dateloc(ind),x,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis(ind),K67.dateloc(ind)); plot(K67.dateloc(ind),x,'color',color_mtx(7,:),'LineWidth',2);
 plot(datenum(2020,6,15,09-3.5,46,05),187.494,'diamond','color',[0.7 0.3 0],'MarkerSize',7);
axis([0 24 -100 700]);      grid on;
legend('0.5TOA','Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}','ECOSTRESS');
legend box off;             xtickangle(0);
datetick('x','HH:MM');      box on;          ylabel({'Energy Flux and';'Storage (W m^-^2)'},'FontSize',12);

%%
figure('color','white');
subplot(2,2,3);              hold on;
ind = find((K67.dateloc>=datenum(2019,8,12))&(K67.dateloc<datenum(2019,8,13)));
y = AM_cycle(K67.netrad(ind),K67.dateloc(ind));   plot(K67.dateloc(ind),y./y,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE(ind),K67.dateloc(ind));       plot(K67.dateloc(ind),x./y,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H(ind),K67.dateloc(ind));        plot(K67.dateloc(ind),x./y,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat(ind),K67.dateloc(ind)); 	plot(K67.dateloc(ind),x./y,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE(ind),K67.dateloc(ind));      plot(K67.dateloc(ind),x./y,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis(ind),K67.dateloc(ind)); plot(K67.dateloc(ind),x./y,'color',color_mtx(7,:),'LineWidth',2);
axis([0 24 -.400 1]);     grid on;
legend('Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}','FontSize',12);
legend box off;             xtickangle(0);
datetick('x','HH:MM');      box on;          ylabel({'Flux and Storage';'Normalized by Rn'},'FontSize',12);
title ('2020/06/15'); 
subplot(2,2,1);                         hold on;
x = AM_cycle(K67.So(ind)./2,K67.dateloc(ind));    plot(K67.dateloc(ind),x,'color',color_mtx(1,:),'LineWidth',2);
x = AM_cycle(K67.netrad(ind),K67.dateloc(ind));   plot(K67.dateloc(ind),x,'color',color_mtx(2,:),'LineWidth',2);
x = AM_cycle(K67.LE(ind),K67.dateloc(ind));       plot(K67.dateloc(ind),x,'color',color_mtx(3,:),'LineWidth',2);
x = AM_cycle(K67.H(ind),K67.dateloc(ind));        plot(K67.dateloc(ind),x,'color',color_mtx(4,:),'LineWidth',2);
x = AM_cycle(K67.StorageHeat(ind),K67.dateloc(ind)); 	plot(K67.dateloc(ind),x,'color',color_mtx(5,:),'LineWidth',2);
x = AM_cycle(K67.StorageLE(ind),K67.dateloc(ind));      plot(K67.dateloc(ind),x,'color',color_mtx(6,:),'LineWidth',2);
x = AM_cycle(K67.Ephotosyntesis(ind),K67.dateloc(ind)); plot(K67.dateloc(ind),x,'color',color_mtx(7,:),'LineWidth',2);
 plot(datenum(2019,8,12,11-3.5,30,24),187.494,'diamond','Color',[0.7 0.3 0],'MarkerSize',7);
axis([0 24 -100 700]);      grid on;         xtickangle(0);
% legend('0.5TOA','Rn','LE','H','\DeltaS_{H}','\DeltaS_{LE}','E_{p}','ECOSTRESS');
% legend box off;            
datetick('x','HH:MM');      box on;          ylabel({'Energy Flux and';'Storage (W m^-^2)'},'FontSize',12);

%%
color_jet = jet(8);
figure('color','white');   subplot(2,2,1); hold on;
plot(K67.dateloc(ind),K67.tson(ind),'color','k','LineWidth',2);     
plot(K67.dateloc(ind),K67.tair1(ind),'color',color_jet(1,:),'LineWidth',2); 
plot(K67.dateloc(ind),K67.tair2(ind),'color',color_jet(2,:),'LineWidth',2); 
plot(K67.dateloc(ind),K67.tair3(ind),'color',color_jet(3,:),'LineWidth',2);
plot(K67.dateloc(ind),K67.tair4(ind),'color',color_jet(4,:),'LineWidth',2);
plot(K67.dateloc(ind),K67.tair5(ind),'color',color_jet(5,:),'LineWidth',2);
plot(K67.dateloc(ind),K67.tair6(ind),'color',color_jet(6,:),'LineWidth',2);
plot(K67.dateloc(ind),K67.tair7(ind),'color',color_jet(7,:),'LineWidth',2);
plot(K67.dateloc(ind),K67.tair8(ind),'color',color_jet(8,:),'LineWidth',2);
legend('Sonic anemometer','62 m','53 m','39 m','28 m','20 m','10 m','3 m','1 m'); 
legend box off;
axis([0 24 22 34]);      grid on;         xtickangle(0);
datetick('x','HH:MM');      box on;          ylabel({'Energy Flux and';'Storage (W m^-^2)'},'FontSize',12);
ylabel('Temperature (C)','FontSize',12);  