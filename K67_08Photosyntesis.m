% Compare radiation from Harvard 

figure;
subplot(2,2,1);
plot(K67h.totsun(K67h.totsun>0),K67h.netrad(K67h.totsun>0),'.')
ylabel('netrad (W m^-^2)'); xlabel('totsun (\mumol m^-^2 s^-^1)');
axis([0 2000 0 1000]);
subplot(2,2,2);
plot(K67.totsun(K67.totsun>0),K67.netrad(K67.totsun>0),'.')
ylabel('netrad (W m^-^2)'); xlabel('totsun (\mumol m^-^2 s^-^1)');
axis([0 1000 0 1000]);
subplot(2,2,3);
plot(K67h.totsun(K67h.totsun>0),K67h.netrad(K67h.totsun>0)./0.794,'.')
ylabel('netrad (W m^-^2)'); xlabel('totsun (\mumol m^-^2 s^-^1)');
axis([0 2000 0 1000]);

K67h.difsun_Wm2=K67h.difsun.*0.4431;                        K67h.totsun_Wm2=K67h.totsun.*0.4431;

K67.dateloc_all=[K67h.dateloc;K67.dateloc];
[K67.Y_all,K67.M_all,K67.D_all,K67.HR_all,~,~]=datevec(K67.dateloc_all);
K67.difsun_all=[K67h.difsun_Wm2;K67.difsun];
K67.totsun_all=[K67h.totsun_Wm2;K67.totsun];

subplot(2,2,4);
plot(K67h.totsun_Wm2(K67h.totsun_Wm2>0),K67h.netrad(K67h.totsun_Wm2>0),'.')
ylabel('netrad (W m^-^2)'); xlabel('totsun (\mumol m^-^2 s^-^1)');
axis([0 1000 0 1000]);

%QAQC Dif and totsun
% AM_plot_bimonth_qaqc(K67.totsun_all,K67.difsun_all,K67.difsun_all.*NaN,K67.difsun_all.*NaN,K67.Y_all,K67.M_all,K67.HR_all,'K67',' Sunshine (W m^-^2)');

K67.DifFraction=K67.difsun./K67.totsun;        K67.DifFraction(K67.totsun<10)=NaN;

% ........................................................................................
%Daily averages
% ........................................................................................
[K67.day_difsun,K67.day_difsun_sum,K67.day_difsun_count]=AM_day(K67.difsun,K67.dateloc,0); 
K67.day_difsun(K67.day_difsun_count<21)= NaN;
K67.day_difsun_smooth = filter(b,a,K67.day_difsun);
K67.day_difsun_smooth(1)=K67.day_difsun(1);                   K67.day_difsun_smooth(end)=K67.day_difsun(end);

[K67.day_totsun,K67.day_totsun_sum,K67.day_totsun_count]=AM_day(K67.totsun,K67.dateloc,0); 
K67.day_totsun(K67.day_totsun_count<21)= NaN;
K67.day_totsun_smooth = filter(b,a,K67.day_totsun);
K67.day_totsun_smooth(1)=K67.day_totsun(1);      K67.day_totsun_smooth(end)=K67.day_totsun(end);

[K67.day_DifFraction,K67.day_DifFraction_sum,K67.day_DifFraction_count]=AM_day(K67.DifFraction,K67.dateloc,0); 
K67.day_DifFraction(K67.day_DifFraction_count<8)= NaN;
K67.day_DifFraction_smooth = filter(b,a,K67.day_DifFraction);
K67.day_DifFraction_smooth(1)=K67.day_DifFraction(1);         K67.day_DifFraction_smooth(end)=K67.day_DifFraction(end);
% ........................................................................................
% Eight day averages
% ........................................................................................
[K67eight_difsun,K67eight_count_difsun,~,~,~,K67eight_difsun_std]=AM_eightday_rs(K67.day_difsun,K67.day); 
K67eight_difsun((K67eight_count_difsun<2))= NaN;
K67eight_difsun=AM_fill(K67eight_difsun);
K67eight_difsun_smooth = filter(b,a,K67eight_difsun);
K67eight_difsun_smooth(1)=K67eight_difsun(1);               K67eight_difsun_smooth(end)=K67eight_difsun(end);
% K67eight_avg_difsun_smooth=AM_eight_avg(K67eight_difsun_smooth);

[K67eight_totsun,K67eight_count_totsun,~,~,~,K67eight_totsun_std]=AM_eightday_rs(K67.day_totsun,K67.day); 
K67eight_totsun((K67eight_count_totsun<2))= NaN;
K67eight_totsun=AM_fill(K67eight_totsun);
K67eight_totsun_smooth = filter(b,a,K67eight_totsun);
K67eight_totsun_smooth(1)=K67eight_totsun(1);               K67eight_totsun_smooth(end)=K67eight_totsun(end);
% K67eight_avg_totsun_smooth=AM_eight_avg(K67eight_totsun_smooth);

[K67eight_DifFraction,K67eight_count_DifFraction,~,~,~,K67eight_DifFraction_std]=AM_eightday_rs(K67.day_DifFraction,K67.day); 
K67eight_DifFraction((K67eight_count_DifFraction<2))= NaN;
K67eight_DifFraction=AM_fill(K67eight_DifFraction);
K67eight_DifFraction_smooth = filter(b,a,K67eight_DifFraction);
K67eight_DifFraction_smooth(1)=K67eight_DifFraction(1);     K67eight_DifFraction_smooth(end)=K67eight_DifFraction(end);
% K67eight_avg_DifFraction_smooth=AM_eight_avg(K67eight_DifFraction_smooth);

% ........................................................................................
% Previous years
% ........................................................................................
[K67h.day_difsun,~,x]=AM_day(K67h.difsun_Wm2,K67h.dateloc);         K67h.day_difsun(x<20)=NaN;
[K67h.day_totsun,~,x]=AM_day(K67h.totsun_Wm2,K67h.dateloc);         K67h.day_totsun(x<20)=NaN;
[K67h.day_DifFraction,~,x]=AM_day(K67h.DifFraction,K67h.dateloc);   K67h.day_DifFraction(x<20)=NaN;

[K67h.day_DifFraction,~,x]=AM_day(K67h.DifFraction,K67h.dateloc);	K67h.day_DifFraction(x<10)=NaN;
[K67.day_DifFraction,~,x]=AM_day(K67.DifFraction,K67.dateloc);      K67.day_DifFraction(x<10)=NaN;
[K67h.day_totsun,~,x]=AM_day(K67h.totsun_Wm2,K67h.dateloc);         K67h.day_totsun(x<10)=NaN;
[K67.day_totsun,~,x]=AM_day(K67.totsun,K67.dateloc);                K67.day_totsun(x<10)=NaN;
[K67h.day_difsun,~,x]=AM_day(K67h.difsun_Wm2,K67h.dateloc);         K67h.day_difsun(x<10)=NaN;
[K67.day_difsun,~,x]=AM_day(K67.difsun,K67.dateloc);                K67.day_difsun(x<10)=NaN;

K67h.week2_difsun=AM_week2day_rs(K67h.day_difsun,K67h.day);
K67h.week2_totsun=AM_week2day_rs(K67h.day_totsun,K67h.day);
K67h.week2_DifFraction=AM_week2day_rs(K67h.day_DifFraction,K67h.day);

K67h.week2_avg_difsun=AM_week2_avg(K67h.week2_difsun);
K67h.week2_avg_totsun=AM_week2_avg(K67h.week2_totsun);
K67h.week2_avg_DifFraction=AM_week2_avg(K67h.week2_DifFraction);

[K67eight_difsun,K67eight_count_difsun,~,~,~,K67eight_difsun_std] = AM_eightday_rs(K67.day_difsun,K67.day); 
K67eight_difsun((K67eight_count_difsun<2))= NaN;
K67eight_difsun=AM_fill(K67eight_difsun);
K67eight_difsun_smooth = filter(b,a,K67eight_difsun);
K67eight_difsun_smooth(1)=K67eight_difsun(1);            K67eight_difsun_smooth(end)=K67eight_difsun(end);
[K67eight_avg_difsun_smooth,~,K67eight_avg_difsun_std] = AM_eight_avg(K67eight_difsun_smooth);

[K67eight_totsun,K67eight_count_totsun,~,~,~,K67eight_totsun_std]=AM_eightday_rs(K67.day_totsun,K67.day); 
K67eight_totsun((K67eight_count_totsun<2))= NaN;
K67eight_totsun=AM_fill(K67eight_totsun);
K67eight_totsun_smooth = filter(b,a,K67eight_totsun);
K67eight_totsun_smooth(1)=K67eight_totsun(1);            K67eight_totsun_smooth(end)=K67eight_totsun(end);
[K67eight_avg_totsun_smooth,~,K67eight_avg_totsun_std] = AM_eight_avg(K67eight_totsun_smooth);

[K67eight_DifFraction,K67eight_count_DifFraction,~,~,~,K67eight_DifFraction_std]=AM_eightday_rs(K67.day_DifFraction,K67.day); 
K67eight_DifFraction((K67eight_count_DifFraction<2)) = NaN;
K67eight_DifFraction=AM_fill(K67eight_DifFraction);
K67eight_DifFraction_smooth = filter(b,a,K67eight_DifFraction);
K67eight_DifFraction_smooth(1)=K67eight_DifFraction(1);      K67eight_DifFraction_smooth(end)=K67eight_DifFraction(end);
[K67eight_avg_DifFraction_smooth,~,K67eight_avg_DifFraction_std]=AM_eight_avg(K67eight_DifFraction_smooth);

[K67h.eight_DifFraction,x,~,~,~,K67h.eight_DifFraction_std] = AM_eightday_rs(K67h.day_DifFraction,K67h.day); 
K67h.eight_DifFraction((x<2))=NaN;
K67h.eight_DifFraction=AM_fill(K67h.eight_DifFraction);
K67h.eight_DifFraction_smooth = filter(b,a,K67h.eight_DifFraction);
K67h.eight_DifFraction_smooth(1)=K67h.eight_DifFraction(1);     K67h.eight_DifFraction_smooth(end)=K67h.eight_DifFraction(end);
[K67h.eight_avg_DifFraction_smooth,~,K67h.eight_avg_DifFraction_std]=AM_eight_avg(K67h.eight_DifFraction_smooth);

[K67h.eight_totsun,x,~,~,~,K67h.eight_totsun_std] = AM_eightday_rs(K67h.day_totsun,K67h.day); 
K67h.eight_totsun((x<2))= NaN;
K67h.eight_totsun=AM_fill(K67h.eight_totsun);
K67h.eight_totsun_smooth = filter(b,a,K67h.eight_totsun);
K67h.eight_totsun_smooth(1)=K67h.eight_totsun(1);            K67h.eight_totsun_smooth(end)=K67h.eight_totsun(end);
[K67h.eight_avg_totsun_smooth,~,K67h.eight_avg_totsun_std] = AM_eight_avg(K67h.eight_totsun_smooth);

[K67h.eight_difsun,x,~,~,~,K67h.eight_difsun_std] = AM_eightday_rs(K67h.day_difsun,K67h.day); 
K67h.eight_difsun((x<2))= NaN;
K67h.eight_difsun = AM_fill(K67h.eight_difsun);
K67h.eight_difsun_smooth = filter(b,a,K67h.eight_difsun);
K67h.eight_difsun_smooth(1) = K67h.eight_difsun(1);          K67h.eight_difsun_smooth(end)=K67h.eight_difsun(end);
[K67h.eight_avg_difsun_smooth,~,K67h.eight_avg_difsun_std] = AM_eight_avg(K67h.eight_difsun_smooth);

% ........................................................................................
% Diffsun
% ........................................................................................
K67h.eight_difsun_smooth_2002=K67h.eight_difsun_smooth;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67h.eight_difsun_smooth_2002(ind)=[]; 

K67h.eight_difsun_smooth_2003=K67h.eight_difsun_smooth;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67h.eight_difsun_smooth_2003(ind)=[]; 

K67h.eight_difsun_smooth_2004=K67h.eight_difsun_smooth;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67h.eight_difsun_smooth_2004(ind)=[]; 

K67h.eight_difsun_smooth_2005=K67h.eight_difsun_smooth;   
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67h.eight_difsun_smooth_2005(ind)=[]; 

K67eight_difsun_smooth_2008=K67eight_difsun_smooth;   
ind=find((K67eight_days<datenum(2008,1,1))|(K67eight_days>datenum(2008,12,31))); 
K67eight_difsun_smooth_2008(ind)=[]; 

K67eight_difsun_smooth_2009=K67eight_difsun_smooth;   
ind=find((K67eight_days<datenum(2009,1,1))|(K67eight_days>datenum(2009,12,31))); 
K67eight_difsun_smooth_2009(ind)=[]; 

K67eight_difsun_smooth_2010=K67eight_difsun_smooth;   
ind=find((K67eight_days<datenum(2010,1,1))|(K67eight_days>datenum(2010,12,31))); 
K67eight_difsun_smooth_2010(ind)=[]; 

% ........................................................................................
%totsun
% ........................................................................................
K67h.eight_totsun_smooth_2002=K67h.eight_totsun_smooth;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67h.eight_totsun_smooth_2002(ind)=[]; 

K67h.eight_totsun_smooth_2003=K67h.eight_totsun_smooth;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67h.eight_totsun_smooth_2003(ind)=[]; 

K67h.eight_totsun_smooth_2004=K67h.eight_totsun_smooth;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67h.eight_totsun_smooth_2004(ind)=[]; 

K67h.eight_totsun_smooth_2005=K67h.eight_totsun_smooth;   
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67h.eight_totsun_smooth_2005(ind)=[]; 

K67eight_totsun_smooth_2008=K67eight_totsun_smooth;   
ind=find((K67eight_days<datenum(2008,1,1))|(K67eight_days>datenum(2008,12,31))); 
K67eight_totsun_smooth_2008(ind)=[]; 

K67eight_totsun_smooth_2009=K67eight_totsun_smooth;   
ind=find((K67eight_days<datenum(2009,1,1))|(K67eight_days>datenum(2009,12,31))); 
K67eight_totsun_smooth_2009(ind)=[]; 

K67eight_totsun_smooth_2010=K67eight_totsun_smooth;   
ind=find((K67eight_days<datenum(2010,1,1))|(K67eight_days>datenum(2010,12,31))); 
K67eight_totsun_smooth_2010(ind)=[]; 

% ........................................................................................
%DifFraction
% ........................................................................................
K67h.eight_DifFraction_smooth_2002=K67h.eight_DifFraction_smooth;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67h.eight_DifFraction_smooth_2002(ind)=[]; 

K67h.eight_DifFraction_smooth_2003=K67h.eight_DifFraction_smooth;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67h.eight_DifFraction_smooth_2003(ind)=[]; 

K67h.eight_DifFraction_smooth_2004=K67h.eight_DifFraction_smooth;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67h.eight_DifFraction_smooth_2004(ind)=[]; 

K67h.eight_DifFraction_smooth_2005=K67h.eight_DifFraction_smooth;   
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67h.eight_DifFraction_smooth_2005(ind)=[]; 

K67eight_DifFraction_smooth_2008=K67eight_DifFraction_smooth;   
ind=find((K67eight_days<datenum(2008,1,1))|(K67eight_days>datenum(2008,12,31))); 
K67eight_DifFraction_smooth_2008(ind)=[]; 

K67eight_DifFraction_smooth_2009=K67eight_DifFraction_smooth;   
ind=find((K67eight_days<datenum(2009,1,1))|(K67eight_days>datenum(2009,12,31))); 
K67eight_DifFraction_smooth_2009(ind)=[]; 

K67eight_DifFraction_smooth_2010=K67eight_DifFraction_smooth;   
ind=find((K67eight_days<datenum(2010,1,1))|(K67eight_days>datenum(2010,12,31))); 
K67eight_DifFraction_smooth_2010(ind)=[]; 

% ........................................................................................
% ........................................................................................
% Figures
% ........................................................................................
% ........................................................................................
%%
v1=-2;              v2=4;      	v3=v1:((v2-v1)/3):v2; %NEE
x1=0;               x2=600;     x3=x1:((x2-x1)/3):x2; %prec
w1=20;              w2=320;     w3=w1:((w2-w1)/3):w2; %totsun
y1=4;               y2=12;      y3=y1:((y2-y1)/4):y2; %GEE, R
m=(1:12);

figure('color','white','PaperOrientation','portrait');
subplot(3,1,1);
hl3 = plot(K67eight_days,(K67eight_NEE_smooth.*gC_m2_d1),...
    'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
hl4 = plot(K67eight_days,(K67eight_NEE_smooth.*gC_m2_d1),...
    'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% box off;                                         %datetick ('x',4);                      
set(gca,'XTick',[K67_month],'Xlim',[K67_month(1) K67_month(end)],'Ylim',[v1 v2]);   
set(gca,'Layer','top')
set(gca,'XTickLabel',[ '' '' ''])
set(gca,'YTick',v3,'FontSize',10); 
ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'NEE_u_*_2_2';'(gC m^-^2 day^-^1)'},'FontSize',10);
text (K67_month(2),9,'K67','FontSize',10);
xlim=([K67_month(1) K67_month(end)]);                   ylim=([v1 v2]);
set(gca,'XTickLabel',[ '' '' ''])

subplot(3,1,2);
hl3 = plot(K67eight_days,(-K67eight_GEE_smooth.*gC_m2_d1),...
    'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% box off;                                         %datetick ('x',4);                      
set(gca,'XTick',[K67_month],'Xlim',[K67_month(1) K67_month(end)],'Ylim',[y1 y2]);   
set(gca,'Layer','top')
set(gca,'XTickLabel',[ '' '' ''])
set(gca,'YTick',y3,'FontSize',10); 
ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.46 0.28 0.20]);
ylabel ({'GEE';'(gC m^-^2 day^-^1)'},'FontSize',10);

ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
hold on;
hl2 = plot(K67eight_days,(K67eight_R_smooth.*gC_m2_d1),...
    'Color',[0.5 0.0 0.0],'LineWidth',2); hold on;
hold on;
set(gca,'XTick',[K67_month],'Xlim',[K67_month(1) K67_month(end)],'Ylim',[y1 y2]);   
set(gca,'XTickLabel','','FontSize',10);    
set(gca,'YTick',y3,'FontSize',10); 
set(ax2,'XColor','k','YColor',[ 0.5 0.0 0.0],'Position',[0.17 0.46 0.28 0.2]);
ylabel ({'R';'(gC m^-^2 day^-^1)'},'FontSize',10);

subplot(3,1,3);
% hl1 = bar(K67eight_days,K67eight_rain,1,'FaceColor',[ 0.6 0.6 0.6],'EdgeColor',[ 0.5 0.5 0.5]); hold on;
% hl2 = plot(K67eight_days,(K67eight_rain.*0)+100,'Color','k','LineWidth',1','LineStyle',':'); hold on;
% hl1 = bar(K67eight_days,K67eight_Belterra_rain,1,'FaceColor',[ 0.6 0.6 0.6],'EdgeColor',[ 0.5 0.5 0.5]); hold on;
% hl2 = plot(K67eight_days,(K67eight_Belterra_rain.*0)+100,'Color','k','LineWidth',1','LineStyle',':'); hold on;
hl1 = bar(K67eight_days,K67eight_rain_corr,1,'FaceColor',[ 0.6 0.6 0.6],'EdgeColor',[ 0.5 0.5 0.5]); hold on;
hl2 = plot(K67eight_days,(K67eight_rain_corr.*0)+100,'Color','k','LineWidth',1','LineStyle',':'); hold on;
set(gca,'XTick',[K67_month],'Xlim',[K67_month(1) K67_month(end)],'Ylim',[x1 x2]);   
set(gca,'Layer','top')
set(gca,'XTickLabel',[ '' '' ''])
set(gca,'YTick',x3,'FontSize',10); 
ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.20 0.28 0.2]);
ylabel ({'prec';'(mm mo^-^1)'},'FontSize',10);
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
hold on;
hl5 = plot(K67eight_days,K67eight_totsun+K67eight_totsun_std,...
    'Color',[0.9 0.4 0.4],'LineWidth',1,'LineStyle',':'); hold on;
hold on;
hl6 = plot(K67eight_days,K67eight_totsun-K67eight_totsun_std,...
    'Color',[0.9 0.4 0.4],'LineWidth',1,'LineStyle',':'); hold on;
hl7 = plot(K67eight_days,K67eight_totsun,...
    'Color',[0.5 0.0 0.0],'LineWidth',2); hold on;

hl8 = plot(K67eight_days,K67eight_difsun+K67eight_difsun_std,...
    'Color',[0.9 0.4 0.4],'LineWidth',1,'LineStyle',':'); hold on;
hold on;
hl9 = plot(K67eight_days,K67eight_difsun-K67eight_difsun_std,...
    'Color',[0.9 0.4 0.4],'LineWidth',1,'LineStyle',':'); hold on;
hl0 = plot(K67eight_days,K67eight_difsun,...
    'Color',[0.5 0.0 0.0],'LineWidth',2,'LineStyle','--'); hold on;

set(gca,'XTick',[K67_month],'Xlim',[K67_month(1) K67_month(end)],'Ylim',[w1 w2]);   
set(gca,'XTickLabel',{datestr(K67_month,4)},'FontSize',10);    
set(gca,'YTick',w3,'FontSize',10); 
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.20 0.28 0.2]);
ylabel ({'totsun';'(W m^-^2)'},'FontSize',10);
%
% ........................................................................................
%Environment
% ........................................................................................
x1=0.3;           x2=.8;        x3=x1:((x2-x1)/2):x2; %DifFrac
y1=60;            y2=260;       y3=y1:((y2-y1)/2):y2; %difsun
w1=4;             w2=8;      	w3=w1:((w2-w1)/2):w2; %LAI
v1=4.5;           v2=6.5;      	v3=v1:((v2-v1)/2):v2; %ta
z1=ALL_eight_avg(1)-8;          z2=ALL_eight_avg(end)+8;
p1=-1;            p2=3;      	p3=p1:((p2-p1)/2):p2; %Wood increment
q1=-1;            q2=3;       q3=q1:((q2-q1)/2):q2; %prec

% ........................................................................................
%difsun
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(y2);
x4(23:end)=y2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl2 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2003,...
    'Color',[0.6 0.6 0.6],'LineStyle',':','LineWidth',2.5); hold on;
hl3 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2004,...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2005,...
    'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_difsun_smooth_2008,...
    'Color',[0.9 0.1 0.1],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_difsun_smooth_2009))), K67eight_difsun_smooth_2009,...
    'Color',[0.1 0.1 0.9],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_difsun_smooth_2010))), K67eight_difsun_smooth_2010,...
    'Color',[0.1 0.9 0.9],'LineWidth',2.5); hold on;
% xlim([ALL_month_avg(1) 349]);                       ylim([v1 v2]);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[y1 y2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',y3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'difsun (W m^-^2)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
% ........................................................................................
%totsun
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(y2);
x4(23:end)=y2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67h.eight_totsun_smooth_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl2 = plot(ALL_eight_avg, K67h.eight_totsun_smooth_2003,...
    'Color',[0.6 0.6 0.6],'LineStyle',':','LineWidth',2.5); hold on;
hl3 = plot(ALL_eight_avg, K67h.eight_totsun_smooth_2004,...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67h.eight_totsun_smooth_2005,...
    'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_totsun_smooth_2008,...
    'Color',[0.9 0.1 0.1],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_totsun_smooth_2009))), K67eight_totsun_smooth_2009,...
    'Color',[0.1 0.1 0.9],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_totsun_smooth_2010))), K67eight_totsun_smooth_2010,...
    'Color',[0.1 0.9 0.9],'LineWidth',2.5); hold on;
% xlim([ALL_month_avg(1) 349]);                       ylim([v1 v2]);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[y1 y2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',y3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'totsun (W m^-^2)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
%DifFraction
% ........................................................................................
%%
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(x2);
x4(23:end)=x2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67h.eight_DifFraction_smooth_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl2 = plot(ALL_eight_avg, K67h.eight_DifFraction_smooth_2003,...
    'Color',[0.6 0.6 0.6],'LineStyle',':','LineWidth',2.5); hold on;
hl3 = plot(ALL_eight_avg, K67h.eight_DifFraction_smooth_2004,...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67h.eight_DifFraction_smooth_2005,...
    'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_DifFraction_smooth_2008,...
    'Color',[0.9 0.1 0.1],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_DifFraction_smooth_2009))), K67eight_DifFraction_smooth_2009,...
    'Color',[0.1 0.1 0.9],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_DifFraction_smooth_2010))), K67eight_DifFraction_smooth_2010,...
    'Color',[0.1 0.9 0.9],'LineWidth',2.5); hold on;
% xlim([ALL_month_avg(1) 349]);                       ylim([v1 v2]);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[x1 x2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',x3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Diffuse.Total Radiation^-^1'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

figure('color','white');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(x2);
x4(23:end)=x2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2002./K67h.eight_totsun_smooth_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl2 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2003./K67h.eight_totsun_smooth_2003,...
    'Color',[0.6 0.6 0.6],'LineStyle',':','LineWidth',2.5); hold on;
hl3 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2004./K67h.eight_totsun_smooth_2004,...
    'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67h.eight_difsun_smooth_2005./K67h.eight_totsun_smooth_2005,...
    'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_difsun_smooth_2008./K67eight_totsun_smooth_2008,...
    'Color',[0.9 0.1 0.1],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg, K67eight_difsun_smooth_2009./K67eight_totsun_smooth_2009,...
    'Color',[0.1 0.1 0.9],'LineWidth',2.5); hold on;
hl6 = plot(ALL_eight_avg(1:max(size(K67eight_difsun_smooth_2010))),...
    K67eight_difsun_smooth_2010./K67eight_totsun_smooth_2010,...
    'Color',[0.1 0.9 0.9],'LineWidth',2.5); hold on;
% xlim([ALL_month_avg(1) 349]);                       ylim([v1 v2]);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[x1 x2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',x3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.38 0.44 0.36]);
ylabel ({'DifFraction'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
%%
% ........................................................................................
% ........................................................................................
%LUE
% ........................................................................................
% ........................................................................................
%%

monthstr = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';...
    'Nov';'Dec'];
hour=0:1:23;
x=-100:10:900;
% 12 months and 7 years: 2002 2003 2004 2005 2008 2009
K67Pc=NaN(12,6);    K67GEPhypR2=NaN(12,6);      K67GEPmax=NaN(12,6);
K67alpha=NaN(12,6);

K67.DifFraction_all=[K67h.DifFraction;K67.DifFraction];
K67.difsun_all=[K67h.difsun;K67.difsun];
K67.totsun_all=[K67h.totsun;K67.totsun];
K67GEE_all=[K67h.GEE;K67GEE];
K67par_all=[K67h.par_fill;K67par];
K67.M_all=[K67h.M;K67.M];

% 
% for k=1:12;       %month
% figure ('color','white','name',['K67 LUE' k]);
%     subplot(2,2,1);
%     % DifFraction >.85 and A.M.
%     ind=find((K67.DifFraction_all<0.85)|(isnan(K67.DifFraction_all))|(K67.M_all~=k)|...
%         (K67par_all<40)|(isnan(K67par_all))|(isnan(K67GEE_all))&(-K67GEE_all>=0)|...
%         (K67.M_all>=12));
%     CKvar1=K67par_all;              CKvar1(ind)=[];
%     CKvar2=-K67GEE_all;             CKvar2(ind)=[];
%     if max(size(CKvar1))>1;
%         [GEPhypPARnum,GEPavg,GEPstd,GEPserror]=AM_LUE_a_GEPsat_clean(CKvar1,CKvar2,[0.2 0.2 0.2],[0.3 0.3 0.3],'DifFrac >0.85 AM');
%         K67alpha(k,1)=GEPhypPARnum(1);  K67GEPmax(k,1)=GEPhypPARnum(2);     K67GEPhypR2(k,1)=GEPhypPARnum(3);      
%         K67Pc(k,1)=GEPavg(1);
%     end;
%     hold on;
%     % DifFraction >.85 and P.M.
%     ind=find((K67.DifFraction_all<0.85)|(isnan(K67.DifFraction_all))|(K67.M_all~=k)|...
%         (K67par_all<40)|(isnan(K67par_all))|(isnan(K67GEE_all))&(-K67GEE_all>=0)|...
%         (K67.M_all<12));
%     CKvar1=K67par_all;              CKvar1(ind)=[];
%     CKvar2=-K67GEE_all;             CKvar2(ind)=[];
%     if max(size(CKvar1))>1;
%         [GEPhypPARnum,GEPavg,GEPstd,GEPserror]=AM_LUE_a_GEPsat_clean(CKvar1,CKvar2,[0.4 0.4 0.4],[0.5 0.5 0.5],'DifFrac >0.85 PM');
%         K67alpha(k,1)=GEPhypPARnum(1);  K67GEPmax(k,1)=GEPhypPARnum(2);     K67GEPhypR2(k,1)=GEPhypPARnum(3);      
%         K67Pc(k,1)=GEPavg(1);
%     end;
%     % DifFraction <.25 and A.M.
%     ind=find((K67.DifFraction_all>0.25)|(isnan(K67.DifFraction_all))|(K67.M_all~=k)|...
%         (K67par_all<40)|(isnan(K67par_all))|(isnan(K67GEE_all))&(-K67GEE_all>=0)|...
%         (K67.M_all>=12));
%     CKvar1=K67par_all;            	CKvar1(ind)=[];
%     CKvar2=-K67GEE_all;             CKvar2(ind)=[];
%     if max(size(CKvar1))>1;
%         [GEPhypPARnum,GEPavg,GEPstd,GEPserror]=AM_LUE_a_GEPsat_clean(CKvar1,CKvar2,[0.9 0.1 0.1],[0.9 0.5 0.5],'DifFrac <0.25 AM');
%         if (max(size(GEPhypPARnum)))>1;
%             K67alpha(k,5)=GEPhypPARnum(1);  
%             K67GEPmax(k,5)=GEPhypPARnum(2);
%             K67GEPhypR2(k,5)=GEPhypPARnum(3);
%             K67Pc(k,5)=GEPavg(1);
%         end;
%     end;
%     % DifFraction <.25 and P.M.
%     ind=find((K67.DifFraction_all>0.25)|(isnan(K67.DifFraction_all))|(K67.M_all~=k)|...
%         (K67par_all<40)|(isnan(K67par_all))|(isnan(K67GEE_all))&(-K67GEE_all>=0)|...
%         (K67.M_all<12));
%     CKvar1=K67par_all;            	CKvar1(ind)=[];
%     CKvar2=-K67GEE_all;             CKvar2(ind)=[];
%     if max(size(CKvar1))>1;
%         [GEPhypPARnum,GEPavg,GEPstd,GEPserror]=AM_LUE_a_GEPsat_clean(CKvar1,CKvar2,[0.9 0.5 0.5],[0.9 0.7 0.7],'DifFrac <0.25 PM');
%         if (max(size(GEPhypPARnum)))>1;
%             K67alpha(k,5)=GEPhypPARnum(1);  
%             K67GEPmax(k,5)=GEPhypPARnum(2);
%             K67GEPhypR2(k,5)=GEPhypPARnum(3);
%             K67Pc(k,5)=GEPavg(1);
%         end;
%     end;
% 
%     FileNameTitle = [ monthstr(k,:) ];
%     title ([FileNameTitle],'FontName','Arial','FontWeight','Bold');
%     xlabel ('PAR (\mumol m^-^2 s^-^1)','FontName','Arial','FontWeight','Bold');
%     ylabel ('GEP (\mumol m^-^2 s^-^1)','FontName','Arial','FontWeight','Bold');
%     axis square;    axis([ 0 1600 0 50]);
%     set(gca,'XTick',0:800:1600,'YTick',0:25:50,'YAxisLocation','left');
%     set(gca,'XTickLabel',{'0','800','1600'})
%     set(gca,'FontName','Arial','FontWeight','Bold');
% end;
%%
% ........................................................................................
%LAI
% ........................................................................................
%Seca Foresta Nepstad et al.
%(ftp://lba.cptec.inpe.br/lba_archives/LC/LC-14/Nepstad/LeafAreaIndex.csv)
%1:Y 2:M 3:day 4:no samples 5:mLAI4 average	6:sLAI4	standard error of the mean 
%7:mGAP4 maximum (?)	8:sGAP4 minimum(?)
K67_LAI=load('C:\Users\ncoupe\Documents\AmazonGEP\K67\SecaForesta\LeafAreaIndex.dat');
K67.dayLAI=datenum((K67_LAI(:,1)),(K67_LAI(:,2)),(K67_LAI(:,3)));
K67LAI=(K67_LAI(:,5));
% K67LAI(K67LAI<4.5)=NaN;
K67eight_LAI=interp1(K67.dayLAI,K67LAI,K67h.eight_days,'linear');
K67eight_LAI_smooth=smooth(K67eight_LAI,.1,'loess');
K67eight_LAI_smooth(isnan(K67eight_LAI))=NaN;
K67eight_avg_LAI=AM_eight_avg(K67eight_LAI);

% ........................................................................................
% ........................................................................................
figure;
subplot(2,2,1);
plot(K67h.eight_days,K67eight_LAI); hold on;
plot(K67.dayLAI,K67LAI,'.r'); 

K67eight_LAI_2002=K67eight_LAI;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67eight_LAI_2002(ind)=[]; 

K67eight_LAI_2003=K67eight_LAI;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67eight_LAI_2003(ind)=[]; 

K67eight_LAI_2004=K67eight_LAI;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67eight_LAI_2004(ind)=[]; 

K67eight_LAI_2005=K67eight_LAI;   
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67eight_LAI_2005(ind)=[]; 

% %Marcos Costa
% K67_LAI=load('C:\Users\ncoupe\Documents\AmazonGEP\K67\Litter\lba_km67_LAI_MCosta.dat');
% K67.dayLAI_Costa=datenum(2003,K67_LAI(:,1),1);
% K67LAI_Costa=K67_LAI(:,2);
% K67eight_LAI_Costa=interp1(K67.dayLAI_Costa,K67LAI_Costa,K67h.eight_days,'linear');
% subplot(2,2,2);
% plot(K67h.eight_days,K67eight_LAI_Costa); hold on;
% plot(K67.dayLAI_Costa,K67LAI_Costa,'.r'); 

% ........................................................................................
%LAI
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(w2);
x4(23:end)=w2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl3 = plot(ALL_eight_avg, K67eight_LAI_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67eight_LAI_2003,...
    'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_LAI_2004,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2.5); hold on;
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[w1 w2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',w3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'LAI (m^2 m^-^2'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
legend ([hl3,hl4,hl5],'2002','2003','2004');
%
% ........................................................................................
%LAI
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(v2);
x4(23:end)=v2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_LAI,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[v1 v2],'Xlim',[z1 z2])
set(gca,'XTickLabel',['',''],'FontSize',10,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',v3,'FontSize',10); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58  0.28 0.2]);
ylabel ({'LAI (m^2 m^-^2)'},'FontSize',10);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
% ........................................................................................

% [K67_sci5x5TERRAday,K67_sci5x5TERRAEVI_day_avgQA] = ...
%     textread('C:\Users\ncoupe\Documents\AmazonRS\EVI\K67_sci5x5TERRAEVI_day_avgQA.dat',...
%     '%f %f','delimiter', ',','headerlines',0);

% ........................................................................................
%Wood
% ........................................................................................
[K67.dayDBH,K67survey_wood_increment] = ...
    textread('C:\Users\ncoupe\Documents\AmazonGEP\AllSites091231\K67survey_wood_increment_litter.dat',...
    '%f %f','delimiter', ',','headerlines',0);

K67eight_Wood=interp1(K67.dayDBH,K67survey_wood_increment,K67h.eight_days,'linear');
K67eight_avg_Wood=AM_eight_avg(K67eight_Wood);

% ........................................................................................
% ........................................................................................
figure;
subplot(2,2,1);
plot(K67h.eight_days,K67eight_Wood); hold on;
plot(K67.dayDBH,K67survey_wood_increment,'.r'); 

K67eight_Wood_2002=K67eight_Wood;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67eight_Wood_2002(ind)=[]; 

K67eight_Wood_2003=K67eight_Wood;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67eight_Wood_2003(ind)=[]; 

K67eight_Wood_2004=K67eight_Wood;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67eight_Wood_2004(ind)=[]; 

K67eight_Wood_2005=K67eight_Wood;   
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67eight_Wood_2005(ind)=[]; 
%%
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl3 = plot(ALL_eight_avg, K67eight_Wood_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67eight_Wood_2003,...
    'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_Wood_2004,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2.5); hold on;
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Tree Wood Increment';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
legend ([hl3,hl4,hl5],'2002','2003','2004');
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_Wood,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Tree Wood Increment';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
%Leaffall
% ........................................................................................
[K67.daylitter,K67leaf] = ...
    textread('C:\Users\ncoupe\Documents\AmazonGEP\AllSites091231\K67survey_litter.dat',...
    '%f %f','delimiter', ',','headerlines',0);

K67eight_Leaffall=interp1(K67.daylitter,K67leaf,K67h.eight_days,'linear');
K67eight_avg_Leaffall=AM_eight_avg(K67eight_Leaffall);

% ........................................................................................
% ........................................................................................
figure;
subplot(2,2,1);
plot(K67h.eight_days,K67eight_Leaffall); hold on;
plot(K67.daylitter,K67leaf,'.r'); 

K67eight_Leaffall_2002=K67eight_Leaffall;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67eight_Leaffall_2002(ind)=[]; 

K67eight_Leaffall_2003=K67eight_Leaffall;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67eight_Leaffall_2003(ind)=[]; 

K67eight_Leaffall_2004=K67eight_Leaffall;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67eight_Leaffall_2004(ind)=[]; 

K67eight_Leaffall_2005=K67eight_Leaffall;   
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67eight_Leaffall_2005(ind)=[]; 
% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl3 = plot(ALL_eight_avg, K67eight_Leaffall_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67eight_Leaffall_2003,...
    'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_Leaffall_2004,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2.5); hold on;
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Leaf-fall';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
legend ([hl3,hl4,hl5],'2002','2003','2004');

% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_Leaffall,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Leaf-fall';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');


% ........................................................................................
% ........................................................................................

K67eight_LAI_der=NaN((K67h.g*46),1);
for m=2:((K67g*46)-1);
    z=K67eight_LAI_smooth(m-1:m+1);  
    w=K67h.eight_days(m-1:m+1)-K67h.eight_days(m-1);
    p = polyfit(w,z,2);
%     K67_week2_Pc_dt_smooth(m)=(-2*p(1,1)/p(1,2));
    K67eight_LAI_der(m)=(2*p(1,1)*(w(2)))+(p(1,2));
end;
K67eight_LAI_der=AM_rm_outlier(K67eight_LAI_der,2);
% %Tomas Domingues 
% %SLA km_67	Dry	0.012062088     Wet	0.012488142 m2 g-1
% K67eight_avg_SLA=NaN(46,1);     K67eight_SLA=NaN(46.*K67h.g,1);
% K67eight_avg_SLA(1:23)=0.012488142.*2; %I am assuming this is dry matter as well
% K67eight_avg_SLA(24:46)=0.012062088.*2; 
% K67eight_SLA(1:46)=K67eight_avg_SLA;    
% ind=find(isnan(K67eight_SLA));  K67eight_SLA(ind)=K67eight_SLA(ind-46);
% LAM g m-2 148.3

K67eight_LAI_der=K67eight_LAI_der.*74.15;
K67eight_avg_LAI_der=AM_eight_avg(K67eight_LAI_der);

K67eight_flush=K67eight_LAI_der+K67eight_Leaffall;
K67eight_avg_flush=AM_eight_avg(K67eight_flush);

% ........................................................................................
% ........................................................................................
K67eight_flush_2002=K67eight_flush;     K67eight_LAI_der_2002=K67eight_LAI_der;
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67eight_flush_2002(ind)=[];            K67eight_LAI_der_2002(ind)=[];

K67eight_flush_2003=K67eight_flush;     K67eight_LAI_der_2003=K67eight_LAI_der;
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67eight_flush_2003(ind)=[];             K67eight_LAI_der_2003(ind)=[];

K67eight_flush_2004=K67eight_flush;     K67eight_LAI_der_2004=K67eight_LAI_der;
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67eight_flush_2004(ind)=[];             K67eight_LAI_der_2004(ind)=[];

K67eight_flush_2005=K67eight_flush;     K67eight_LAI_der_2005=K67eight_LAI_der;
ind=find((K67h.eight_days<datenum(2005,1,1))|(K67h.eight_days>datenum(2005,12,31))); 
K67eight_flush_2005(ind)=[];             K67eight_LAI_der_2005(ind)=[];

% ........................................................................................
p1=-1;           p2=3;      	p3=p1:((p2-p1)/2):p2; %Wood increment
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_flush_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl2 = plot(ALL_eight_avg, K67eight_flush_2003,...
    'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',2.5); hold on;
hl3 = plot(ALL_eight_avg, K67eight_flush_2004,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2.5); hold on;
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Leaf-flush';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_flush,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'Leaf-flush';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_LAI_der_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl2 = plot(ALL_eight_avg, K67eight_LAI_der_2003,...
    'Color',[0.2 0.2 0.2],'LineStyle',':','LineWidth',2.5); hold on;
hl3 = plot(ALL_eight_avg, K67eight_LAI_der_2004,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2.5); hold on;
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[q1 q2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',q3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
% ylabel ({'dLAI dt^-^1';'(day^-^1)'},'FontSize',14);
ylabel ({'LMA * dLAI dt^-^1';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(q2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_LAI_der,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[q1 q2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',q3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
% ylabel ({'dLAI dt^-^1';'(day^-^1)'},'FontSize',14);
ylabel ({'LMA * dLAI dt^-^1';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

%%
% % ........................................................................................
% % ........................................................................................
% 
% % ........................................................................................
% % ........................................................................................
%%
p1=0;           p2=4;      	p3=p1:((p2-p1)/2):p2; %Wood increment
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=0;   x4(1:23)=p2;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_flush,...
    'Color',[0.4 0.4 0.9],'LineWidth',2);
hl3 = plot(ALL_eight_avg, K67eight_avg_Wood,...
    'Color',[0.9 0.4 0.4],'LineWidth',2);
hl4 = plot(ALL_eight_avg, K67eight_avg_flush+K67eight_avg_Wood,...
    'Color',[0.2 0.2 0.2],'LineWidth',2);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',['',''],'FontSize',10,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',10); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.28 0.2]);
ylabel ({'Productivity';'(gC m^-^2 day^-^1)'},'FontSize',10);
set(gca,'FontName','Arial','FontWeight','Normal');
legend('prec','LAI flush','Wood Increment','AGNPP');

figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(q2);
x4(23:end)=q2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_LAI_der,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[q1 q2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',q3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
% ylabel ({'dLAI dt^-^1';'(day^-^1)'},'FontSize',14);
ylabel ({'';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

p1=-1;           p2=3;      	p3=p1:((p2-p1)/2):p2; %Wood increment
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_LAI_der,...
    'Color',[0.4 0.4 0.9],'LineWidth',2.5);
hl3 = plot(ALL_eight_avg, K67eight_avg_Leaffall,...
    'Color',[0.9 0.4 0.4],'LineWidth',2.5);
hl4 = plot(ALL_eight_avg, K67eight_avg_flush,...
    'Color',[0.2 0.2 0.2],'LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
legend('prec','prec','LMA * dLAI dt^-^1','Leaf-fall','LAI flush');
%%

q1=7;           q2=10;      	q3=q1:((q2-q1)/2):q2; %Wood increment
x=AM_rm_outlier(K67h.eight_GEE_smooth_2002,3);
x(2)=(x(1)+x(3))./2;
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(q2);
x4(23:end)=q2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg,-1.*(x+K67h.eight_GEE_smooth_2003+...
    K67h.eight_GEE_smooth_2004)./3,...
    'Color','k','LineWidth',2.5);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[q1 q2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',q3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
% ylabel ({'dLAI dt^-^1';'(day^-^1)'},'FontSize',14);
ylabel ({'GEP';'(gC m^-^2 day^-^1)'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');

% ........................................................................................
% ........................................................................................
%CUE = AGNPP/ GPP
% ........................................................................................
% ........................................................................................
K67h.eight_GEE_smooth(2)=(K67h.eight_GEE_smooth(1)+K67h.eight_GEE_smooth(3))./2;
K67eight_CUE=(K67eight_flush+K67eight_Wood)./...
    (-1.*K67h.eight_GEE_smooth);
K67eight_avg_CUE=AM_eight_avg(K67eight_CUE);

K67eight_CUE_2002=K67eight_CUE;   
ind=find((K67h.eight_days<datenum(2002,1,1))|(K67h.eight_days>datenum(2002,12,31))); 
K67eight_CUE_2002(ind)=[]; 

K67eight_CUE_2003=K67eight_CUE;   
ind=find((K67h.eight_days<datenum(2003,1,1))|(K67h.eight_days>datenum(2003,12,31))); 
K67eight_CUE_2003(ind)=[]; 

K67eight_CUE_2004=K67eight_CUE;   
ind=find((K67h.eight_days<datenum(2004,1,1))|(K67h.eight_days>datenum(2004,12,31))); 
K67eight_CUE_2004(ind)=[]; 

% % ........................................................................................
p1=0;           p2=1;      	p3=p1:((p2-p1)/2):p2; %Wood increment
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=p2;   x4(1:23)=0;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl3 = plot(ALL_eight_avg, K67eight_CUE_2002,...
    'Color','k','LineWidth',2.5); hold on;
hl4 = plot(ALL_eight_avg, K67eight_CUE_2003,...
    'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',2.5); hold on;
hl5 = plot(ALL_eight_avg, K67eight_CUE_2004,...
    'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2.5); hold on;
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',[{datestr(ALL_month_avg(m),4)}],'FontSize',12,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',12); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58 0.44 0.36]);
ylabel ({'CUE';'(Leaf+Wood Production) GEP^-^1'},'FontSize',14);
set(gca,'FontName','Arial','FontWeight','Normal');
legend ([hl3,hl4,hl5],'2002','2003','2004');
%%
% % ........................................................................................
figure('color','white','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(p2);
x4(23:end)=0;   x4(1:23)=p2;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67eight_avg_CUE,...
    'Color','k','LineWidth',2);
m=1:12;
set(gca,'XTickLabel',[ '' '' '']);
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[p1 p2],'Xlim',[z1 z2])
set(gca,'XTickLabel',['',''],'FontSize',10,...
    'FontName','Arial','FontWeight','Normal')
set(gca,'YTick',p3,'FontSize',10); ax1 = gca;  
set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58  0.28 0.20]);
ylabel ({'CUE'},'FontSize',10);
set(gca,'FontName','Arial','FontWeight','Normal');

%%
y1=280;            y2=420;       y3=y1:((y2-y1)/2):y2; %difsun
x1=0.4;              x2=0.8;         x3=x1:((x2-x1)/2):x2; %difsun
figure('color','white');
subplot(3,1,2);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(y2);
x4(23:end)=0;   x4(1:23)=y2;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]); hold on;
hl3 = plot(ALL_eight_avg,K67h.eight_avg_par_smooth,...
    'Color',[0.9 0.4 0.2],'LineWidth',2); hold on;
% box off;                                         %datetick ('x',4);                      
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[y1 y2],'Xlim',[z1 z2])
set(gca,'Layer','top')
set(gca,'XTickLabel',[ '' '' ''])
set(gca,'YTick',y3,'FontSize',10,'YAxisLocation','right'); 
ax1 = gca;  
set(ax1,'XColor','k','YColor',[0.9 0.4 0.2],'Position',[0.17 0.46 0.28 0.20]);
ylabel ({'PAR';'(\mumol m^-^2 s^-^1)'},'FontSize',10);

ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k');
hold on;
hl2 = plot(ALL_eight_avg,K67h.eight_avg_DifFraction_smooth,...
    'Color',[0.2 0.2 0.2],'LineWidth',2,'LineStyle',':'); hold on;
hold on;
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[x1 x2],'Xlim',[z1 z2])
set(gca,'XTickLabel','','FontSize',10);    
set(gca,'YTick',x3,'FontSize',10); 
set(ax2,'XColor','k','YColor',[0.2 0.2 0.2],'Position',[0.17 0.46 0.28 0.2]);
set(gca,'XTickLabel',['',''],'FontSize',10,...
    'FontName','Arial','FontWeight','Normal')
ylabel ({'Diffuse Fraction'},'FontSize',10);

%%
% % ........................................................................................
% %Dif Fraction
% % ........................................................................................
v1=0.4;              v2=0.8;         v3=v1:((v2-v1)/2):v2; %difsun
figure('color','none','PaperOrientation','portrait');
subplot(4,1,1);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(v2);
x4(23:end)=0;   x4(1:23)=v2;
% hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
% hl2 = bar(ALL_eight_avg,-x4,1,'FaceColor',[ 0.8 0.8 0.8],'EdgeColor',[ 0.8 0.8 0.8]); hold on;
hl1 = plot(ALL_eight_avg, K67h.eight_avg_DifFraction_smooth,...
    'Color','k','LineWidth',2.5,'LineStyle','--'); hold on;
hl2 = plot([ALL_eight_avg(1)-6;ALL_eight_avg(1)-6],...
    [v1;v2],'Color','k','LineWidth',2,'LineStyle','-');
hl3 = plot([ALL_eight_avg(1)-6;ALL_eight_avg(1)+4],...
    [v3(1);v3(1)],'Color','k','LineWidth',2,'LineStyle','-');
hl4 = plot([ALL_eight_avg(1)-6;ALL_eight_avg(1)+4],...
    [v3(2);v3(2)],'Color','k','LineWidth',2,'LineStyle','-');
hl5 = plot([ALL_eight_avg(1)-6;ALL_eight_avg(1)+4],...
    [v3(3);v3(3)],'Color','k','LineWidth',2,'LineStyle','-');
m=1:12;
set(gca,'XTickLabel',[ '' '' ''],'XTick',[],'Xlim',[z1 z2],...
    'XColor','white');
set(gca,'Ylim',[v1 v2])
set(gca,'YTick',v3,'FontSize',10,'YTickLabel',[ '' '' '']); ax1 = gca;  
set(ax1,'YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.58  0.28 0.2]);
% ylabel ({'Diffuse Fraction'},'FontSize',10);
set(gca,'FontName','Arial','FontWeight','Normal','Color','none');
axis off


%%
y1=5;            y2=6.2;       y3=y1:((y2-y1)/2):y2; 
x1=0;            x2=4;         x3=x1:((x2-x1)/2):x2;
figure('color','white');
subplot(3,1,2);
x4=zeros(46,1); x4((K67eight_avg_Belterra_rain_smooth.*8)<10000)=(y2);
x4(23:end)=0;   x4(1:23)=y2;
hl2 = bar(ALL_eight_avg,x4,1,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]); hold on;
hl3 = plot(ALL_eight_avg,K67eight_avg_LAI,...
    'Color',[0.2 0.4 0.9],'LineWidth',2); hold on;
% box off;                                         %datetick ('x',4);                      
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[y1 y2],'Xlim',[z1 z2])
set(gca,'Layer','top')
set(gca,'XTickLabel',[ '' '' ''])
set(gca,'YTick',y3,'FontSize',10,'YAxisLocation','right'); 
ax1 = gca;  
set(ax1,'XColor','k','YColor',[0.2 0.4 0.9],'Position',[0.17 0.46 0.28 0.20]);
ylabel ({'LAI';'(m^2 m^-^2)'},'FontSize',10);

ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k');
hold on;
hl2 = plot(ALL_eight_avg,K67eight_avg_Leaffall,...
    'Color',[0.2 0.2 0.2],'LineWidth',2,'LineStyle',':'); hold on;
hold on;
set(gca,'XTick',[ALL_month_avg(m)],'Ylim',[x1 x2],'Xlim',[z1 z2])
set(gca,'XTickLabel','','FontSize',10);    
set(gca,'YTick',x3,'FontSize',10); 
set(ax2,'XColor','k','YColor',[0.2 0.2 0.2],'Position',[0.17 0.46 0.28 0.2]);
set(gca,'XTickLabel',['',''],'FontSize',10,...
    'FontName','Arial','FontWeight','Normal')
ylabel ({'Litterfall';'(gC m^-^2 day^-^1)'},'FontSize',10);
