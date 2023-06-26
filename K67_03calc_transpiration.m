% Partition of ET between transpiration and transpiration --model
% Natalia Restrepo-Coupe
%Toronto, June 2021
%..........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;
FolderK67lai       = '/home/ncoupe/Documents/AmazonGEP/K67/SecaForesta/';
FolderK67laiLIDAR  = '/home/ncoupe/Documents/AmazonGEP/K67/LAI/';

month_name = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];
cycle_hour = (0:23)';

color_mtx = [0.1020 0.5882 0.2549; 215/255  25/255 28/255; 43/255 131/255 186/255; ...
    253/255 174/255 97/255; 0.2549 0.1020 0.5882; 171/255 221/255 164/255;...
    0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

tresholdSo = 100;

% .........................................................................
%% Leaf level VPD
% .........................................................................
% rho_a, density of air:  1.225 kg/m3     % Cp, specific heat capacity of air: 0.001013 J/kg K
% sigma, psychrometric constant: 0.054 [kPa C-1]
% Gs: from mm s-1 to m s-1
K67.VPDleaf = K67.sigma.*K67.LE./(K67.rho_a.*Cp.*K67.GsPMall.*1000);
K67.VPDleaf(K67.VPDleaf>10) = NaN;    %kPa
K67.VPDleaf((K67.So==0)&(isnan(K67.VPDleaf))) = 0.01;     K67.VPDleaf = AM_fill_plus(K67.VPDleaf);
K67.VPDleaf(K67.So<=tresholdSo) = NaN;

figure('color','white'); plot(K67.VPD,K67.VPDleaf,'.'); hold on; plot(K67.VPD,K67.VPD,'-k');
axis([0 4 0 4]);    xlabel('VPD'); ylabel('VPD_{leaf}');

%% ........................................................................
% Gs = Go + G1.*(GPP./VPDleaf.^mCoef);
% m, G0:Gsoil and G1:Gveg are fitted coefficients
% modelspec = 'Gs ~ 1 + ix';
modelGs = NaN;
mCoef_range = (0.5:0.02:1.0)';            color_mtx = jet(length(mCoef_range));
mdlRMSE = NaN(length(mCoef_range),1);     mdlCoefficients = [mdlRMSE,mdlRMSE];
mdlR2   = mdlRMSE;

K67.weekGsoil = NaN(length(K67.week),1);    K67.weekGveg = K67.weekGsoil;
K67.weekmCoef = K67.weekGsoil;              K67.weekGvegGsoilRMSE = K67.weekGsoil;
K67.weekGvegGsoilR2 = K67.weekGsoil;        K67.weekGvegTMP = K67.weekGsoil;

K67.weekVPDleaf = AM_week2cycle_rs(K67.VPDleaf,K67.dateloc,tresholdCycleDay,3,2);

weekMultiplier = NaN(length(K67.week),length(mCoef_range));

for im = 1:length(mCoef_range)
    gepVPD = -K67.gee./(K67.VPDleaf.^mCoef_range(im));      % only 16-day
    gepVPD((isinf(gepVPD))|(K67.So<tresholdSo)|(-K67.gee<=0)|(K67.VPDleaf<=0)) = NaN;
    weekMultiplier(:,im) = AM_week2cycle_rs(gepVPD,K67.dateloc,tresholdCycleDay,3,2);
end

%% ........................................................................
figure ('color','white');
plot(K67.week,weekMultiplier(:,1),'g');  hold on;
plot(K67.week,weekMultiplier(:,2)','y')
plot(K67.week,weekMultiplier(:,3)','b');  plot(K67.week,weekMultiplier(:,4));
% plot(weekMultiplier(:,5));  plot(weekMultiplier(:,6));    plot(weekMultiplier(:,7));
plot(K67.week,weekMultiplier(:,length(mCoef_range))','r');    plot(K67.week,K67.weekGsPMallmmol./10,'k')
legend('GEP/(VPD_{leaf}^{m}) m=0.4','m=0.45','m=0.5','m=0.55','m=0.6',... %%     ,'m=0.65','m=0.7', 'm=1',...
    'Gs_{PM}x100 (mol m^-^1 s^-^2)'); %,'m=0.85','m=0.9','m=0.95'
ylabel('(mol m^-^1 s^-^2)'); xlim([K67.week(1) K67.week(end)]);
datetick('x'); grid on;
%% ........................................................................
figure('color','white');
ix = find((-K67.gee>0)&(-K67.gee>0)&(K67.VPDleaf>0.01));
plot(-K67.gee(ix)./(K67.VPDleaf(ix).^0.5),K67.GsPMallmmol(ix)./1000,'.');
xlabel('GEP/VPD_{l}^m');      ylabel('Gs (mol m^-^2 s^-^1)');
axis([0 150 0 6])
K67.Gveg  = NaN(length(K67.dateloc),1);     K67.Gsoil = K67.Gveg;

%% ........................................................................
GsPMallTMP = K67.GsPMallmmol;       GsPMallTMP(K67.So==0) = 0;      GsPMallTMP = AM_fill(GsPMallTMP);
for in = 1:length(K67.week)-1
    ind = find((K67.dateloc>=K67.week(in))&(K67.dateloc<K67.week(in+1)));
    
    color_mtx = jet(length(mCoef_range));
    figure('color','white');    hold on;
    for im = 1:length(mCoef_range)
        mCoef = mCoef_range(im);
        So = K67.So(ind);       gep = -K67.gee(ind);    VPDleaf = K67.VPDleaf(ind);
        GsPMall = K67.GsPMallmmol(ind)./1000;           % mmol m-2 s-1 to mol m-2 s-1
        gepVPD = gep./(VPDleaf.^mCoef);
        iw = find ((isinf(gepVPD))|(So<tresholdSo)|(gep<=0)|(VPDleaf<=0.1)|(isnan(gepVPD))|(isnan(GsPMall)));
        gepVPD(iw) = [];   GsPMall(iw) = [];
        %         gepVPD = AM_rm_outlier(gepVPD,3);     GsPMall = AM_rm_outlier(GsPMall,3);
        display(length(gepVPD))
        if ((~isempty(gepVPD))&(length(gepVPD)>40))
            %             mdl = fitlm(gepVPD,GsPMall);
            %             CoefficientsEstimate = mdl.Coefficients.Estimate;
            %             iy  = polyval(CoefficientsEstimate,gepVPD);
            % Orthogonal linear regression method == AM_compute_regression
            %             [~,I] = fit_2D_data(gepVPD,GsPMall,'yes');
            
            [CoefficientsEstimate,Rsquared] = AM_compute_regression(gepVPD,GsPMall,[.4 .4 .4],0);
            iy  = polyval(CoefficientsEstimate,gepVPD);
            
            
            if (CoefficientsEstimate(1)~=0)
                %% first iteration ............................................
                % % %             figure; plot(iy,GsPMall,'s');       hold on
                iw  = (min(gepVPD):((max(gepVPD)-min(gepVPD))/20):max(gepVPD))';
                iz  = polyval(CoefficientsEstimate,iw);
                % % %             plot (iw,iz,'-','LineWidth',1,'Color','k');     ylabel('Gs_{PM}');   xlabel('model');
                E   = std(iz)*ones(size(iw));
                % % %             errorbar(iw,iz,E)
                slope = abs(((gepVPD.*CoefficientsEstimate(1))+(GsPMall.*-1)+CoefficientsEstimate(2))./(sqrt((CoefficientsEstimate(1).^2)+1))); % distance point to line
                out = find((abs(slope)>(mean(E))));  %for figure
                plot(gepVPD,iy);              hold on;        plot(gepVPD,GsPMall,'s');
                plot (gepVPD(out),GsPMall(out),'.r');
                
                % %                 mdl = fitlm(gepVPD(out),GsPMall(out));      CoefficientsEstimate = mdl.Coefficients.Estimate;
                gepVPD(out) = NaN;             GsPMall(out) = NaN;
                [CoefficientsEstimate,Rsquared] = AM_compute_regression(gepVPD,GsPMall,[.4 .4 .4],0);
                iz  = polyval(CoefficientsEstimate,iw);
                plot(iw,iz,'-k');
                xlabel('GEP VPD_{leaf}^-^1');   ylabel('Gs_{PM}'); text(5,0.8,['m:' num2str(mCoef)],'color',color_mtx(im,:));
                % %                 % second iteration ............................................
                % %                 iz  = polyval(CoefficientsEstimate,iw);
                % %                 E   = std(iz)*ones(size(iw));
                % %                 slope = abs(((gepVPD.*CoefficientsEstimate(1))+(GsPMall.*-1)+CoefficientsEstimate(2))./(sqrt((CoefficientsEstimate(1).^2)+1))); %distance point to line
                % %                 out = find((abs(slope)>(mean(E))));
                % %                 plot (gepVPD(out),GsPMall(out),'.g');
                % %                 % %                 mdl = fitlm(gepVPD(out),GsPMall(out));
                % %                 gepVPD(~out) = NaN;             GsPMall(~out) = NaN;
                % %                 [CoefficientsEstimate,Rsquared] = AM_compute_regression(gepVPD,GsPMall,[.4 .4 .4],0);
                
                if ((CoefficientsEstimate(1)*CoefficientsEstimate(2))<=0)
                    [CoefficientsEstimate,Rsquared] = AM_compute_regression_zero(gepVPD,GsPMall,[.4 .4 .4],0);
                    CoefficientsEstimate = [CoefficientsEstimate,0];
                    iz  = polyval(CoefficientsEstimate,iw);
                    plot(iw,iz,'-','color',color_mtx(im,:));
                end
            end
            mdlRMSE(im) = sqrt(nanmean((GsPMall - iy).^2));
            mdlCoefficients(im,:) = CoefficientsEstimate;
            mdlR2(im)   = Rsquared;
            %         plot(mdl,'color',color_mtx(im,:));            hold on;
            %             mdlCoefficients(mdlCoefficients<0) = NaN;
        end
    end
    close all
    
    % select the m that minimized the RMSE
    if ~isempty(mdlRMSE(~isnan(mdlRMSE)))
        
        iw   = find(mdlR2 == nanmax(mdlR2));                selectmCoef = iw;
        if (~isempty(iw))
            K67.weekGsoil(in) = mdlCoefficients(iw(1),2);
            K67.weekmCoef(in) = mCoef_range(iw(1));
            K67.weekGvegGsoilRMSE(in) = mdlRMSE(iw(1));
            K67.weekGvegGsoilR2(in)   = mdlR2(iw(1));
            
            So = K67.So(ind);       gep = -K67.gee(ind);    VPDleaf = K67.VPDleaf(ind);
            GsPMall = K67.GsPMallmmol(ind)./1000;           % mmol m-2 s-1 to mol m-2 s-1
            gepVPD = gep./(VPDleaf.^mCoef_range(iw(1)));
            gepVPD((isinf(gepVPD))|(So<tresholdSo)|(gep<=0)|(VPDleaf<=0)) = NaN;
            
            K67.Gveg(ind)  = mdlCoefficients(iw(1),1).*gepVPD;
            K67.Gsoil(ind) = mdlCoefficients(iw(1),2)+(0.*gepVPD);
            
            figure;            hold on;
            plot(mdlCoefficients(iw(1),1).*gepVPD);
            plot(mdlCoefficients(iw(1),2)+(0.*gepVPD))
            plot((mdlCoefficients(iw(1),2)+(0.*gepVPD))+(mdlCoefficients(iw(1),1).*gepVPD))
            plot(GsPMall,'g'); legend('Gveg','Gsoil','Gveg+Gsoil','Gs_{PM}')
            
            figure;             hold on;
            plot((mdlCoefficients(iw(1),2)+(0.*gepVPD))+(mdlCoefficients(iw(1),1).*gepVPD),GsPMall,'.g')
            plot(GsPMall,GsPMall,'--k')
            
            K67.weekGvegTMP(in) = weekMultiplier(in,iw(1)).*mdlCoefficients(iw(1),1);
            
        end
    else
        K67.Gsoil(ind) = NaN;       K67.Gveg(ind) = NaN;
        K67.weekmCoef(in) = NaN;    K67.weekGvegGsoilRMSE(in) = NaN;
        K67.weekGvegGsoilR2(in) = NaN;
    end
end

K67.Gveg(K67.Gveg>8) = NaN;

figure;                         hold on;
plot(K67.GsPMallmmol./1000);    plot(K67.Gveg);     plot(K67.Gsoil);

K67.weekGveg = AM_week2cycle_rs(K67.Gveg,K67.dateloc,tresholdCycleDay,3,2);
K67.weekGsoil = AM_week2cycle_rs(K67.Gsoil,K67.dateloc,tresholdCycleDay,3,2);

%%. ........................................................................
figure('color','white');
subplot(2,1,1);
plot(K67.weekGveg,K67.weekGsFGmmol./1000,'.');      hold on;
plot(K67.weekGveg,K67.weekGveg,'--k');
xlabel('Gs_{vegetation}');     ylabel('Flux gradient');

subplot(2,1,2);
plot(K67.Gveg,K67.GsFGmmol./1000,'.');  hold on;
plot(K67.Gveg,K67.Gveg,'--k');          axis([0 5 0 5]);
xlabel('Gs_{vegetation}');              ylabel('Flux gradient');

%%
figure('color','white');            hold on;
plot(K67.weekGsPMmmol./1000,'r')
plot(K67.weekGsFGmmol./1000,'g');     plot(K67.weekGsPMallmmol./1000,'k')
plot(K67.weekGvegTMP,'c');
plot(K67.weekGveg,'m');
plot(K67.weekGsoil,'y');
plot((K67.weekGsoil+K67.weekGveg),'b')
legend('Gs_{veg}','Gs_{vegTMP}','Gs_{PM_{dry}}','Gs_{FluxGradient}','Gs_{PM_{all}}','Gs_{soil}','Gs_{veg+soil}');

ix = K67.weekGsoil+K67.weekGveg-(K67.weekGsPMall);

% ........................................................................
%% Seasonal calculations
% .........................................................................
K67.Transpiration = K67.Gveg.*1000.*K67.ETday./K67.GsPMallmmol;
K67.Evaporation   = K67.Gsoil.*1000.*K67.ETday./K67.GsPMallmmol;

K67.dayTranspiration = AM_day(K67.Transpiration,K67.dateloc,0);   K67.dayTranspiration(K67.dayTranspiration>16)=NaN;
K67.dayEvaporation = AM_day(K67.Evaporation,K67.dateloc,0);

K67.weekTranspiration = AM_week2cycle_rs(K67.Transpiration,K67.dateloc,tresholdCycleDay,3,2);
K67.weekEvaporation = AM_week2cycle_rs(K67.Evaporation,K67.dateloc,tresholdCycleDay,3,2);

K67.weekTranspiration((K67.week==datenum(2009,2,18))|(K67.week>=datenum(2020,11,25))|...
    (K67.week==datenum(2012,8,28))|(K67.week==datenum(2018,06,10))|...
    (K67.week==datenum(2019,3,06))|(K67.week==datenum(2020,11,16))) = NaN;
K67.weekEvaporation((K67.week==datenum(2009,2,18))|(K67.week>=datenum(2020,11,25))) = NaN;
% % K67.weekTranspiration((K67.week==datenum(2002,9,30))|(K67.week==datenum(2006,1,17))|(K67.week==datenum(2008,7,27))|(K67.week>=datenum(2020,11,15))) = NaN;
% % K67.weekEvaporation((K67.week==datenum(2002,9,30))|(K67.week==datenum(2006,1,17))|(K67.week==datenum(2008,7,27))|(K67.week>=datenum(2020,11,15))) = NaN;
K67.weekETrefdryday((K67.week==datenum(2002,1,1))) = NaN;
K67.weekETday((K67.week==datenum(2002,1,1))) = NaN;

ix = smooth(K67.weekTranspiration);     ix(isnan(K67.weekTranspiration)) = NaN;     K67.weekTranspiration = ix;
ix = smooth(K67.weekEvaporation);       ix(isnan(K67.weekEvaporation)) = NaN;       K67.weekEvaporation = ix;
%%
K67.weekTranspiration = AM_fill(K67.weekTranspiration);
K67.weekEvaporation   = AM_fill(K67.weekEvaporation);
K67.weekETrefdryday   = AM_fill(K67.weekETrefdryday);
K67.weekETday   = AM_fill(K67.weekETday);

[K67.week_avgEvaporation,~,K67.week_avgstdEvaporation] = AM_week2_avg(K67.weekEvaporation,K67.week);
[K67.week_avgTranspiration,~,K67.week_avgstdTranspiration] = AM_week2_avg(K67.weekTranspiration,K67.week);

% .........................................................................
figure('color','white');            hold on;
plot(K67.weekETday,'r')
plot(K67.weekTranspiration,'m');     plot(K67.weekEvaporation,'k')

% ........................................................................
K67.monthTranspiration = AM_month2cycle(K67.Transpiration,K67.dateloc,tresholdCycleDay,3,2);
K67.monthEvaporation = AM_month2cycle(K67.Evaporation,K67.dateloc,tresholdCycleDay,3,2);

K67.monthTranspiration((K67.month==datenum(2008,8,1))|(K67.month>=datenum(2020,11,1))) = NaN;
K67.monthEvaporation((K67.month==datenum(2008,8,1))|(K67.month>=datenum(2020,11,1))) = NaN;


%..........................................................................
%% Gs model
%..........................................................................
X = [K67.CI,K67.VPD_fill,K67.radiation,K67.LUEref];
mdl = fitlm(X,K67.GsFGmmol);
CoefficientsEstimateGs = mdl.Coefficients.Estimate;
K67.modelGs  = CoefficientsEstimateGs(1)+(CoefficientsEstimateGs(2).*K67.CI)+...
    (CoefficientsEstimateGs(3).*K67.VPD_fill)+(CoefficientsEstimateGs(4).*K67.radiation)+...
    (CoefficientsEstimateGs(5).*K67.LUEref);
K67.modelGs((isinf(K67.modelGs))|(K67.modelGs<=0)|(K67.So<=tresholdSo)) = NaN;
[K67.weekGsmodel,~,~,~,~,K67.week_stdGsmodel] =...
    AM_week2cycle_rs(K67.modelGs,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsVPDavg  = CoefficientsEstimateGs(1)+(CoefficientsEstimateGs(2).*K67.CI)+...
    (CoefficientsEstimateGs(3).*mean(K67.VPD_fill,'omitnan'))+(CoefficientsEstimateGs(4).*K67.radiation)+...
    (CoefficientsEstimateGs(5).*K67.LUEref);
[K67.weekGsmodelVPDavg,~,~,~,~,K67.week_stdGsmodelVPDavg] =...
    AM_week2cycle_rs(K67.modelGsVPDavg,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsVPDlow  = CoefficientsEstimateGs(1)+(CoefficientsEstimateGs(2).*K67.CI)+...
    (CoefficientsEstimateGs(4).*K67.radiation)+...
    (CoefficientsEstimateGs(5).*K67.LUEref);
[K67.weekGsmodelVPDlow,~,~,~,~,K67.week_stdGsmodelVPDlow] =...
    AM_week2cycle_rs(K67.modelGsVPDlow,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsVPDhigh  = CoefficientsEstimateGs(1)+(CoefficientsEstimateGs(2).*K67.CI)+...
    (CoefficientsEstimateGs(3).*(K67.VPD_fill.*0+2))+(CoefficientsEstimateGs(4).*K67.radiation)+...
    (CoefficientsEstimateGs(5).*K67.LUEref);
[K67.weekGsmodelVPDhigh,~,~,~,~,K67.week_stdGsmodelVPDhigh] =...
    AM_week2cycle_rs(K67.modelGsVPDhigh,K67.dateloc,tresholdCycleDay,4,1);

% ........................................................................
%% Gs model
% ........................................................................
X = [K67.CI,K67.VPD_fill,K67.radiation,K67.TRMM1hCWD];
mdl = fitlm(X,K67.GsFGmmol);
CoefficientsEstimateGsCWD = mdl.Coefficients.Estimate;
K67.modelGsCWD  = CoefficientsEstimateGsCWD(1)+(CoefficientsEstimateGsCWD(2).*K67.CI)+...
    (CoefficientsEstimateGsCWD(3).*K67.VPD_fill)+(CoefficientsEstimateGsCWD(4).*K67.radiation)+...
    (CoefficientsEstimateGsCWD(5).*K67.TRMM1hCWD);
K67.modelGsCWD((isinf(K67.modelGsCWD))|(K67.modelGsCWD<=0)|(K67.So<=tresholdSo)) = NaN;
[K67.weekGsmodelCWD,~,~,~,~,K67.week_stdGsmodelCWD] =...
    AM_week2cycle_rs(K67.modelGsCWD,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsCWDVPDavg  = CoefficientsEstimateGsCWD(1)+(CoefficientsEstimateGsCWD(2).*K67.CI)+...
    (CoefficientsEstimateGsCWD(3).*mean(K67.VPD_fill,'omitnan'))+(CoefficientsEstimateGsCWD(4).*K67.radiation)+...
    (CoefficientsEstimateGsCWD(5).*K67.TRMM1hCWD);
[K67.weekGsmodelCWDVPDavg,~,~,~,~,K67.week_stdGsmodelCWDVPDavg] =...
    AM_week2cycle_rs(K67.modelGsCWDVPDavg,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsCWDVPDlow  = CoefficientsEstimateGsCWD(1)+(CoefficientsEstimateGsCWD(2).*K67.CI)+...
    (CoefficientsEstimateGsCWD(4).*K67.radiation)+...
    (CoefficientsEstimateGsCWD(5).*K67.TRMM1hCWD);
[K67.weekGsmodelCWDVPDlow,~,~,~,~,K67.week_stdGsmodelCWDVPDlow] =...
    AM_week2cycle_rs(K67.modelGsCWDVPDlow,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsCWDVPDhigh  = CoefficientsEstimateGsCWD(1)+(CoefficientsEstimateGsCWD(2).*K67.CI)+...
    (CoefficientsEstimateGsCWD(3).*(K67.VPD_fill.*0+2))+(CoefficientsEstimateGsCWD(4).*K67.radiation)+...
    (CoefficientsEstimateGsCWD(5).*K67.TRMM1hCWD);
[K67.weekGsmodelCWDVPDhigh,~,~,~,~,K67.week_stdGsmodelCWDVPDhigh] =...
    AM_week2cycle_rs(K67.modelGsCWDVPDhigh,K67.dateloc,tresholdCycleDay,4,1);

%..........................................................................
%% Gs model
%..........................................................................
K67.SWshallow(isinf(K67.SWshallow)) = NaN;
X = [K67.CI,K67.VPD_fill,K67.radiation,K67.SWshallow.*100];
mdl = fitlm(X,K67.GsFGmmol);
CoefficientsEstimateGsSM = mdl.Coefficients.Estimate;
K67.modelGsSM  = CoefficientsEstimateGsSM(1)+(CoefficientsEstimateGsSM(2).*K67.CI)+...
    (CoefficientsEstimateGsSM(3).*K67.VPD_fill)+(CoefficientsEstimateGsSM(4).*K67.radiation)+...
    (CoefficientsEstimateGsSM(5).*K67.SWshallow.*100);
K67.modelGsSM((isinf(K67.modelGsSM))|(K67.modelGsSM<=0)|(K67.So<=tresholdSo)) = NaN;
[K67.weekGsmodelSM,~,~,~,~,K67.week_stdGsmodelSM] =...
    AM_week2cycle_rs(K67.modelGsSM,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsSMVPDavg  = CoefficientsEstimateGsSM(1)+(CoefficientsEstimateGsSM(2).*K67.CI)+...
    (CoefficientsEstimateGsSM(3).*mean(K67.VPD_fill,'omitnan'))+(CoefficientsEstimateGsSM(4).*K67.radiation)+...
    (CoefficientsEstimateGsSM(5).*K67.SWshallow.*100);
[K67.weekGsmodelSMVPDavg,~,~,~,~,K67.week_stdGsmodelSMVPDavg] =...
    AM_week2cycle_rs(K67.modelGsSMVPDavg,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsSMVPDlow  = CoefficientsEstimateGsSM(1)+(CoefficientsEstimateGsSM(2).*K67.CI)+...
    (CoefficientsEstimateGsSM(4).*K67.radiation)+...
    (CoefficientsEstimateGsSM(5).*K67.SWshallow.*100);
[K67.weekGsmodelSMVPDlow,~,~,~,~,K67.week_stdGsmodelSMVPDlow] =...
    AM_week2cycle_rs(K67.modelGsSMVPDlow,K67.dateloc,tresholdCycleDay,4,1);

K67.modelGsSMVPDhigh  = CoefficientsEstimateGsSM(1)+(CoefficientsEstimateGsSM(2).*K67.CI)+...
    (CoefficientsEstimateGsSM(3).*(K67.VPD_fill.*0+2))+(CoefficientsEstimateGsSM(4).*K67.radiation)+...
    (CoefficientsEstimateGsSM(5).*K67.SWshallow.*100);
[K67.weekGsmodelSMVPDhigh,~,~,~,~,K67.week_stdGsmodelSMVPDhigh] =...
    AM_week2cycle_rs(K67.modelGsSMVPDhigh,K67.dateloc,tresholdCycleDay,4,1);

%% ........................................................................
figure('color','white');
plot(K67.week,K67.weekGsmodel,'k');         hold on;    plot(K67.week,K67.weekGsFGmmol,'color',[0,.5,0])
plot(K67.week,K67.weekGsmodelVPDavg,'color',[.5,.9,0]);   plot(K67.week,K67.weekGsmodelVPDlow,'color',[0,.5,.9])
plot(K67.week,K67.weekGsmodelVPDhigh,'color',[.7,0,0]);
datetick('x');
legend('Gs_{model}','Gs_{FG}','Gs_{avgVPD}','Gs_{VPD=0}','Gs_{VPD=2}');
ylabel('Gs (mmol m^-^2 s^-^1)');

%% ........................................................................
figure('color','white');
plot(K67.week,K67.weekGsmodelCWD,'k');         hold on;    plot(K67.week,K67.weekGsFGmmol,'color',[0,.5,0])
plot(K67.week,K67.weekGsmodelCWDVPDavg,'color',[.5,.9,0]);   plot(K67.week,K67.weekGsmodelCWDVPDlow,'color',[0,.5,.9])
plot(K67.week,K67.weekGsmodelCWDVPDhigh,'color',[.7,0,0]);
datetick('x');
legend('Gs_{model}','Gs_{FG}','Gs_{avgVPD}','Gs_{VPD=0}','Gs_{VPD=2}');
ylabel('Gs (mmol m^-^2 s^-^1)');

%..........................................................................
%% Gs as PC 
%..........................................................................
ind = find((K67.HY_Y~=2008)&(K67.HY_Y~=2015));
[K67.weekGsasPc,weekGsasPc_count,weekGsasPc_std,PARdayAVG,VPDdayAVG,TAdayAVG,CIdayAVG,...
    PARdaySTD,VPDdaySTD,TAdaySTD,CIdaySTD] = ...
    AM_week2day_Pc_TaPARvpdCI(K67.GsFGmmol(ind),K67.radiation(ind),K67.VPD_fill(ind),...
    K67.tc_fill(ind),K67.dateloc(ind),lat,lon,-1);
K67.weekGsasPc = AM_rm_outlier(K67.weekGsasPc,3);
K67.weekGsasPc = AM_fill(K67.weekGsasPc);
ix = smooth(K67.weekGsasPc);    ix(isnan(K67.weekGsasPc)) = NaN;    K67.weekGsasPc = ix;
[K67.week_avgGsasPc,~,K67.week_avgstdGsasPc] = AM_week2_avg(K67.weekGsasPc,K67.week);

[K67.week_deseasonGsasPc,K67.week_trendGsasPc,K67.week_mtrendGsasPc] =...
    AM_week2trend(K67.weekGsasPc,K67.week,1,'Gs_{FG asPc}(mmol m^-^2 s^-^1)');
PARdayAVG
VPDdayAVG
TAdayAVG
CIdayAVG
PARdaySTD
VPDdaySTD
TAdaySTD
CIdaySTD

%% ........................................................................
ind = find((K67.HY_Y~=2008)&(K67.HY_Y~=2015));
[K67.monthGsasPc,monthGsasPc_count,monthGsasPc_std,PARdayAVG,VPDdayAVG,TAdayAVG,CIdayAVG,...
    PARdaySTD,VPDdaySTD,TAdaySTD,CIdaySTD] = ...
    AM_month_Pc_TaPARvpdCI(K67.GsFGmmol(ind),K67.radiation(ind),K67.VPD_fill(ind),...
    K67.tc_fill(ind),K67.dateloc(ind),lat,lon,-1);
K67.monthGsasPc = AM_rm_outlier(K67.monthGsasPc,3);
K67.monthGsasPc = AM_fill(K67.monthGsasPc);
ix = smooth(K67.monthGsasPc);    ix(isnan(K67.monthGsasPc)) = NaN;    K67.monthGsasPc = ix;
[K67.month_avgGsasPc,~,K67.month_avgstdGsasPc] = AM_month_avg(K67.monthGsasPc,K67.month);

%% ........................................................................
[K67.weekGsasPc24,weekGsasPc24_count,weekGsasPc24_std] = ...
    AM_week2day_Pc_TaPARvpdCI(K67.GsFGmmol24(ind),K67.radiation(ind),K67.VPD_fill(ind),...
    K67.tc_fill(ind),K67.dateloc(ind),lat,lon,-1);
K67.weekGsasPc24 = AM_rm_outlier(K67.weekGsasPc24,3);
K67.weekGsasPc24 = AM_fill(K67.weekGsasPc24);
ix = smooth(K67.weekGsasPc24);    ix(isnan(K67.weekGsasPc24)) = NaN;    K67.weekGsasPc24 = ix;
[K67.week_avgGsasPc24,~,K67.week_avgstdGsasPc24] = AM_week2_avg(K67.weekGsasPc24,K67.week);

[K67.week_deseasonGsasPc24,K67.week_trendGsasPc24,K67.week_mtrendGsasPc24] = AM_week2trend(K67.weekGsasPc24,K67.week,1,'Gs_{FG asPc24}(mmol m^-^2 s^-^1)');

% ........................................................................
%% dry season only
% ........................................................................
ind = find((K67.HY_Y~=2008)&(K67.HY_Y~=2015)&(K67.M>6));
[K67.weekGsasPcDrySeason,weekGsasPcDrySeason_count,weekGsasPcDrySeason_std,PARdayAVG,VPDdayAVG,TAdayAVG,CIdayAVG,...
    PARdaySTD,VPDdaySTD,TAdaySTD,CIdaySTD] = ...
    AM_week2day_Pc_TaPARvpdCI(K67.GsFGmmol(ind),K67.radiation(ind),K67.VPD_fill(ind),...
    K67.tc_fill(ind),K67.dateloc(ind),lat,lon,-1);
K67.weekGsasPcDrySeason = AM_rm_outlier(K67.weekGsasPcDrySeason,3);
K67.weekGsasPcDrySeason = AM_fill(K67.weekGsasPcDrySeason);
ix = smooth(K67.weekGsasPcDrySeason);    ix(isnan(K67.weekGsasPcDrySeason)) = NaN;    K67.weekGsasPcDrySeason = ix;
[K67.week_avgGsasPcDrySeason,~,K67.week_avgstdGsasPcDrySeason] = AM_week2_avg(K67.weekGsasPcDrySeason,K67.week);

[K67.week_deseasonGsasPcDrySeason,K67.week_trendGsasPcDrySeason,K67.week_mtrendGsasPcDrySeason] = AM_week2trend(K67.weekGsasPcDrySeason,K67.week,1,'Gs_{FG asPc dry season}(mmol m^-^2 s^-^1)');
PARdayAVG
VPDdayAVG
TAdayAVG
CIdayAVG
PARdaySTD
VPDdaySTD
TAdaySTD
CIdaySTD

%% ....................................................................
v1 = 200;                         v2 = 1400;       v3 = v1:((v2-v1)/4):v2;
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
hl1=plot(K67.week,K67.weekGsasPc,'LineWidth',1,'color',[0.4 0.4 0.8]); hold on;
hl2=plot(K67.week,K67.week_deseasonGsasPc,'LineWidth',3,'color',[0.6 0.3 0.3]); hold on;
hl4=plot(K67.week,K67.week_trendGsasPc,'LineWidth',2,'color',[0.8 0.6 0.4]); hold on;
hl3=plot(K67.week,(K67.week.*0)+nanmean(K67.weekGsasPc),':','LineWidth',1,'color',[0.4 0.4 0.6]); hold on;
set(gca,'XTick',K67.month(1:36:end),'XTickLabel',datestr(K67.month(1:36:end),'yyyy'),'YTick',v3,...
    'XLim',[K67.month(1) datenum(2021,1,1)],'YLim',[v1 v2],'FontSize',12);
text(K67.month(2),v3(4),['trend: ' num2str((nanmin(K67.week_trendGsasPc)-nanmax(K67.week_trendGsasPc)),'%4.2f')]);
ylabel ('Gs_{FG asPc}(mmol m^-^2 s^-^1)','FontSize',12);

