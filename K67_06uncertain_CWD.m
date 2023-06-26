% NGARS meeting 2016 December figures
% Natalia Restrepo-Coupe
% Noisy PAR
% h2o mixing ratio improperly calculated
% u* improperly calculated
%..........................................................................
flag_timeseries   = 1;
month_start_hydro = 7;

ip = 5;
[site_name_BF, lat_BF, lon_BF, site_long_BF, year_ini, year_end, siteelevation_BF,...
    height_BF, canopy_BF, ust_tresh_BF, sunrise_BF, sunset_BF, utm_BF, ust_min_BF, ust_max_BF,...
    LAI_BF]   =   textread([FolderFlux 'sites_loc.txt'],...
    '%s	%f	%f	%s	%f	%f	%f	%f %f	%f	%f	%f	%f	%f	%f	%f',...
    'delimiter',',','headerlines', 1);

site_name  =  site_name_BF(ip);         lat     =  lat_BF(ip);      lon        =  lon_BF(ip);         site_long  =  site_long_BF(ip);
siteelevation  =  siteelevation_BF(ip); canopy  =  canopy_BF(ip);   ust_tresh  =  ust_tresh_BF(ip);   height = height_BF(ip);
sunrise    =  sunrise_BF(ip);           sunset  =  sunset_BF(ip);   utm        =  utm_BF(ip);         ust_min    =  ust_min_BF(ip);
ust_max    =  ust_max_BF(ip);           LAI     =  LAI_BF(ip);

eight_avg    = (1:8:365)';              week_avg = (1:16:365)';
K67.weekTRMM = interp1(TRMM.dateloc,TRMM.precip,K67.week);

K67.year = unique(K67.Y);

%% ........................................................................
figure('color','white');
subplot(2,2,1);   [p,r2] = AM_compute_regression(K67.SWmodelET,K67.monthET,[.6 .6 .6],1);
ylabel([' ET (mm month^-^1)']);                           xlabel('ET_{SW_{downCERES}}(mm month^-^1)');
text(50,140,['ET =',num2str(p(1,1),3.2), 'ET_{SW_{CERES}} ', num2str(p(1,2),2.1)])
text(50,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
text(50,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);

DiffETmodel    = K67.SWmodelET-K67.monthET;
[DiffETmodelAVG,DiffETmodelSTD] = AM_month_avg(DiffETmodel,K67.month);

v1 = -25;         v2 = 25;  v3 = v1:((v2-v1)/2):v2;
subplot(2,2,2);                         hold on;
x4  = zeros(12,1);                      x4(K67.month_avgTRMM<101) = (v2);
hl1 = bar(month_avg,x4,1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y   = [DiffETmodelAVG-DiffETmodelSTD,2.*DiffETmodelSTD];
h   = area(month_avg,y);                set(gca,'Layer','top');
set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);
hl2 = plot(month_avg,DiffETmodelAVG,'Color',[0 0 0],'LineWidth',2);
hl3 = plot(month_avg,zeros(12,1)+(nanmean(DiffETmodelAVG)),':','Color',[0 0 0],'LineWidth',1);
xlim([month_avg(1) month_avg(end)]);    ylim([v1 v2]);
set(gca,'XTick',[month_avg],'XTickLabel',datestr(month_avg,'m'),'YTick',v3,'FontSize',10);
set(gca,'XColor','k','YColor',[ 0.0 0.0 0.0]);%,'Position',[0.47 0.78 0.22 0.18]);
ylabel ({'\Delta ET_{SWCERES}''(mm month^-^1)'},'FontSize',10);
text (month_avg(2),15,site_name,'FontSize',10);         box on;

DiffCWDmodel    = TRMM.monthCWDceres-TRMM.monthCWD;
[DiffCWDmodelAVG,DiffCWDmodelSTD] = AM_month_avg(DiffCWDmodel,TRMM.dateloc);

subplot(2,2,3);   [p,r2] = AM_compute_regression(TRMM.monthCWDceres,TRMM.monthCWD,[.6 .6 .6],1);
ylabel([' CWD (mm)']);                           xlabel('CWD_{SW_{downCERES}}(mm)');
text(50,140,['CWD =',num2str(p(1,1),3.2), 'CWD_{SW_{CERES}} ', num2str(p(1,2),2.1)])
text(50,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
text(50,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);

v1 = -25;         v2 = 25;  v3 = v1:((v2-v1)/2):v2;
subplot(2,2,4);                         hold on;
x4  = zeros(12,1);                      x4(K67.month_avgTRMM<101) = (v2);
hl1 = bar(month_avg,x4,1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y   = [DiffCWDmodelAVG-DiffCWDmodelSTD,2.*DiffCWDmodelSTD];
h   = area(month_avg,y);                set(gca,'Layer','top');
set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);
hl2 = plot(month_avg,DiffCWDmodelAVG,'Color',[0 0 0],'LineWidth',2);
hl3 = plot(month_avg,zeros(12,1)+(nanmean(DiffCWDmodelAVG)),':','Color',[0 0 0],'LineWidth',1);
xlim([month_avg(1) month_avg(end)]);    ylim([v1 v2]);
set(gca,'XTick',[month_avg],'XTickLabel',datestr(month_avg,'m'),'YTick',v3,'FontSize',10);
set(gca,'XColor','k','YColor',[ 0.0 0.0 0.0]);%,'Position',[0.47 0.78 0.22 0.18]);
ylabel ({'\Delta CWD_{SWCERES}''(mm)'},'FontSize',10);
text (month_avg(2),15,site_name,'FontSize',10);         box on;

figure('color','white');
subplot(2,2,1);   [p,r2] = AM_compute_regression(K67.monthNetRad,K67.monthTRMM,[.6 .6 .6],1);
ylabel(['Precip_{TRMM} (mm month^-^1)']);                           xlabel('Rn (W m^-^2)');
text(150,140,['Precip =',num2str(p(1,1),3.2), 'Rn ', num2str(p(1,2),2.1)])
text(150,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
text(150,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]); axis([50 250 0 600]);

%..........................................................................
%% Identify the drought periods creat the regression
%..........................................................................
ix = nanmean(K67.monthEFday);     iy = nanstd(K67.monthEFday);
K67.monthETdry = K67.monthET;     K67.monthETdry(K67.monthEFday>(ix-iy)) = NaN;

figure('color','white');          % LE/H.....................................
subplot(3,2,1:2);                 hold on;
hl1 = plot(K67.month,K67.monthEFday,'LineWidth',2,'color',[0 0 0]); datetick('x');
ylabel('EF_{daytime}');         grid on;        box on;
hl2 = plot(K67.month,K67.month.*0+ix+iy,'color',[0 0 .9],'LineWidth',2);
hl3 = plot(K67.month,K67.month.*0+ix,':','color',[0 0 0],'LineWidth',1);
hl4 = plot(K67.month,K67.month.*0+ix-iy,'color',[0.9 0 0 ],'LineWidth',2);
% legend([hl1,hl2,hl3,hl4],'EF_{daytime}','ET_{noQA}','ET','ET_{EBCfilter}');

%% Get the regression for ET during drought... identified as EF<EFmu-EFstd
subplot(3,2,3);                     hold on;
[p,r2] = AM_compute_regression(K67.monthCERESswd,K67.monthET,[.6 .6 .6],1);
text(160,140,['ET =',num2str(p(1,1),3.2), 'SW_{CERES} ', num2str(p(1,2),2.1)])
text(160,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
text(160,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);  grid on; box on;
hold on;        grid on;        box on;          ylabel('ET (mm month^-^1)');     xlabel('SW_{down{CERES}}(W m^-^2)');
[p,r2,~,p_val,sm,sb,K67.SWmodelETdry] = AM_compute_regression(K67.monthCERESswd,K67.monthETdry,[.9 .2 .2],1);
text(160,140,['ET_{dry} =',num2str(p(1,1),3.2), 'SW_{CERES} ', num2str(p(1,2),2.1)],'color',[.9 .2 .2])
grid on;        box on;

subplot(3,2,4);                 hold on;
[p,r2] = AM_compute_regression(K67.SWmodelET,K67.monthET,[.4 .4 .4],1);
text(60,140,['ET =',num2str(p(1,1),3.2), 'ET_{SW_{CERES_{dry}}}', num2str(p(1,2),2.1)])
text(60,130,['std slp :',num2str(sm,3.2), ' std int:', num2str(sb,2.1)])
text(60,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);  grid on; box on;
[p,r2] = AM_compute_regression(K67.SWmodelETdry,K67.monthET,[.9 .2 .2],1);
ylabel('ET (mm month^-^1)');                           xlabel('ET_{SW_{downCERES}}(mm month^-^1)');
text(60,140,['ET =',num2str(p(1,1),3.2), 'ET_{SW_{CERES_{dry}}', num2str(p(1,2),2.1)],'color',[.9 .2 .2])

%..........................................................................
%% Regression CERES for dry
%..........................................................................
TRMM.monthETfillDRY = interp1(K67.month,K67.SWmodelETdry,TRMM.dateloc);
TRMM.monthETfillDRY = AM_fill_month_cycle(TRMM.monthETfillDRY,TRMM.dateloc,6,18);

TRMM.monthCWDdry      = zeros(length(TRMM.dateloc),1);
for im=2:length(TRMM.dateloc)
    if (TRMM.monthCWDdry(im-1)-TRMM.monthETfillDRY(im)+TRMM.precip(im)) < 0
        TRMM.monthCWDdry(im) = -(TRMM.monthCWDdry(im-1)-TRMM.monthETfillDRY(im)+TRMM.precip(im));
    else
        TRMM.monthCWDdry(im) = 0;
    end
end

K67.monthCWDdry  = interp1(TRMM.dateloc,TRMM.monthCWDdry,K67.month);

subplot(3,2,5:6);       plot(K67.month,K67.monthCWD,'color',[0. 0. 0.9],'LineWidth',2);
ylabel('CWD (mm)');     grid on;        hold on;
plot(K67.month,K67.monthCWDdry,'color',[0.9 0. 0.],'LineWidth',2);      datetick('x');
text (K67.month(3),-20,site_name);
legend ('CWD_{fill(CERES)}','CWD_{fill(CERESdry)}');

%% Uncertainity analysis
CWDmtx = NaN(length(TRMM.precip),100);    ETmtx = CWDmtx;
ik=1;
while ik<=100
    [monthCWD,ETsw] = AM_uncertainityCWD(TRMM.monthET,TRMM.monthCERESswd,TRMM.precip,TRMM.dateloc,40);
    CWDmtx(:,ik) = monthCWD;
    ETmtx(:,ik)  = ETsw;
    ik = ik+1;
end
%%
color_mtx = jet(100);
figure ('color','white');
subplot(2,1,1);          hold on;
for ik =1:100; plot(TRMM.dateloc,CWDmtx(:,ik),'color',color_mtx(ik,:),'LineWidth',1);   end;
ylabel('CWD (mm)');     grid on;
hl1 = plot(TRMM.dateloc,nanmean(CWDmtx,2),'color',[0 0. 0],'LineWidth',2);      datetick('x');
hl2 = plot(TRMM.dateloc,TRMM.monthCWD,'color',[0 0. 0.9],'LineWidth',2);      datetick('x');
hl3 = plot(TRMM.dateloc,TRMM.monthCWD100,'color',[0.0 0.7 0.7],'LineWidth',2);      datetick('x');
hl4 = plot(TRMM.dateloc,TRMM.monthCWDdry,'color',[0.3 0.3 0.3],'LineWidth',2);      datetick('x');
text (K67.month(3),-20,site_name);      box on;
legend ([hl1,hl2,hl3,hl4],'CWD_{\muMonteCarlo}','CWD','CWD_{ETfix100}','CWD_{fill(CERESdry)}')
subplot(2,1,2);          hold on;
for ik =1:100; plot(TRMM.dateloc,ETmtx(:,ik),'color',color_mtx(ik,:),'LineWidth',1);   end;
ylabel('ET (mm)');     grid on;        box on;
plot(TRMM.dateloc,nanmean(ETmtx,2),'color',[0 0. 0],'LineWidth',2);      datetick('x');
plot(TRMM.dateloc,TRMM.monthETfill,'color',[0 0. 0.9],'LineWidth',2);      datetick('x');
plot(TRMM.dateloc,TRMM.dateloc.*0+100,'color',[0.0 0.7 0.7],'LineWidth',2);      datetick('x');
plot(TRMM.dateloc,TRMM.monthETfillDRY,'color',[0.3 0.3 0.3],'LineWidth',2);      datetick('x');
%%
monthCWD = nanmean(CWDmtx,2);
ETsw = nanmean(ETmtx,2);

figure('color','white');
subplot(2,2,1);
[p,r2] = AM_compute_regression(ETsw,TRMM.monthET,[.9 .2 .2],1);
ylabel('ET (mm)');                           xlabel('ET_{monte carlo}(mm)');
text(50,140,['ET =',num2str(p(1,1),3.2), 'ET_{model}+', num2str(p(1,2),2.1)])
text(50,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);  grid on; box on;

subplot(2,2,2);
[p,r2] = AM_compute_regression(monthCWD,TRMM.monthCWD,[.9 .2 .2],1);
ylabel('CWD (mm)');                           xlabel('CWD_{monte carlo}(mm)');
text(50,140,['CWD =',num2str(p(1,1),3.2), 'CWD_{model} ', num2str(p(1,2),2.1)])
text(50,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);  grid on; box on;

subplot(2,2,3);
[p,r2] = AM_compute_regression(TRMM.monthCWD,TRMM.monthCWD100,[.9 .2 .2],1);
ylabel('CWD_{100mm} (mm)');                           xlabel('CWD(mm)');
text(50,140,['CWD_{100mm} =',num2str(p(1,1),3.2), 'CWD+', num2str(p(1,2),2.1)])
text(50,120,['R2=',num2str(r2,2.1)],'color',[0 0 0]);  grid on; box on;
