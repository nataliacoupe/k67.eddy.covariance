% Biometry calculations, mortality, changes in density, recruitment
% Natalia Restrepo-Coupe
%..........................................................................
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;
FolderK67litter = '/home/ncoupe/Documents/AmazonGEP/K67/Litter/';
FolderK67dendrometry = '/home/ncoupe/Documents/AmazonGEP/K67/Dendrometry/';
FolderK67mortality = '/home/ncoupe/Documents/AmazonGEP/K67/TreeMortality/';

month_name = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];
cycle_hour = (0:23)';

color_mtx = [0.1020 0.5882 0.2549; 215/255  25/255 28/255; 43/255 131/255 186/255; ...
    253/255 174/255 97/255; 0.2549 0.1020 0.5882; 171/255 221/255 164/255;0.1020 0.5882 0.2549; 0.1020 0.5882 0.2549];

% ........................................................................
% ........................................................................
%% Litter
% ........................................................................
% ........................................................................
opts = delimitedTextImportOptions("NumVariables", 15);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["YYYY", "MM", "DD", "contDOY", "dT", "leaf", "leafse", "fruit", "fruitse", "wood", "woodse", "Misc", "Miscse", "total", "totalse"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
Litter = readtable([FolderK67litter 'litter.km67.join.csv'], opts);
clear opts

dayLitter = datenum(Litter{:,1},Litter{:,2},Litter{:,3});

LeafLitter = Litter{:,6};           LeafLitter = LeafLitter.*0.5;
AllLitter = Litter{:,14};           AllLitter  = AllLitter.*0.5;        %dry mass CO2 to C
FruitLitter = Litter{:,8};          WoodLitter = Litter{:,10};
FruitLitter = FruitLitter.*0.5;     WoodLitter  = WoodLitter.*0.5;        %dry mass CO2 to C
MiscLitter = Litter{:,12};          MiscLitter = MiscLitter.*0.5;

LeafLitterStdIni = Litter{:,7};             AllLitterStdIni = Litter{:,15};
LeafLitterStdIni = LeafLitterStdIni.*0.5;   AllLitterStdIni = AllLitterStdIni.*0.5;        %dry mass CO to C
FruitLitterStdIni = Litter{:,9};            WoodLitterStdIni = Litter{:,11};
FruitLitterStdIni = FruitLitterStdIni.*0.5; WoodLitterStdIni = WoodLitterStdIni.*0.5;        %dry mass CO to C
MiscLitterStdIni = Litter{:,13};            MiscLitterStdIni = MiscLitterStdIni.*0.5;
dT = Litter{:,5};

% AllLitter(AllLitter>3500) = NaN;

figure('color','white');        hold on;
plot(LeafLitter(dayLitter<datenum(2007,1,1)));
plot(LeafLitter((dayLitter>datenum(2011,1,1))&(dayLitter<datenum(2014,1,1))));
plot(LeafLitter((dayLitter>datenum(2014,1,1))&(dayLitter<datenum(2017,1,1))));
plot(LeafLitter((dayLitter>datenum(2017,1,1))&(dayLitter<datenum(2018,11,1))));
plot(LeafLitter((dayLitter>datenum(2018,11,1))&(dayLitter<datenum(2021,1,1))));
legend ('Harvard','Mauricio','GoAmazon','Cosme','Isotraits');      grid on;

figure('color','white');
plot(dayLitter,LeafLitter);      hold on
plot(dayLitter,LeafLitter,'.');  datetick('x')
grid on;                         ylabel('Litter_{leaf}')

% ........................................................................
%% from basket to forest
% ........................................................................
ix = LeafLitter;                iy = dayLitter;
[ah,~,ch] = unique(iy);         out_size = (histcounts(iy,ah)); out_size = [1,out_size]';
out_hour = [ah,accumarray(ch,ix,[],@nanmean)];
out_hourstd = [ah,accumarray(ch,ix,[],@nanstd)];
ix = out_hour(:,2);                 iy = out_hour(:,1);
in = out_hourstd(:,2)./sqrt(out_size);
LeafLitter = ix;                    LeafLitterStd = in;

iz = AllLitter;
out_hour = [ah,accumarray(ch,iz,[],@nanmean)];
out_hourstd = [ah,accumarray(ch,iz,[],@nanstd)];
iz = out_hour(:,2);
im = out_hourstd(:,2)./sqrt(out_size);
AllLitter = iz;                     AllLitterStd = im;

ir = FruitLitter;                   iy = dayLitter;
[ah,~,ch] = unique(iy);             out_size = (histcounts(iy,ah)); out_size = [1,out_size]';
out_hour = [ah,accumarray(ch,ir,[],@nanmean)];
out_hourstd = [ah,accumarray(ch,ir,[],@nanstd)];
ir = out_hour(:,2);                 iy = out_hour(:,1);
it = out_hourstd(:,2)./sqrt(out_size);
FruitLitter = ir;     FruitLitterStd = it;

iq = WoodLitter;
out_hour = [ah,accumarray(ch,iq,[],@nanmean)];
out_hourstd = [ah,accumarray(ch,iq,[],@nanstd)];
iq = out_hour(:,2);
iu = out_hourstd(:,2)./sqrt(out_size);
WoodLitter = iq;        WoodLitterStd = iu;

iq = MiscLitter;                   iy = dayLitter;
[ah,~,ch] = unique(iy);             out_size = (histcounts(iy,ah)); out_size = [1,out_size]';
out_hour = [ah,accumarray(ch,iq,[],@nanmean)];
out_hourstd = [ah,accumarray(ch,iq,[],@nanstd)];
iq = out_hour(:,2);                 iy = out_hour(:,1);
it = out_hourstd(:,2)./sqrt(out_size);
MiscLitter = iq;     MiscLitterStd = it;

iw = dT;
out_hour = [ah,accumarray(ch,iw,[],@nanmean)];
iw = out_hour(:,2);
dT = iw;        dayLitter  = iy;


% Harvard change units
ind = find(dayLitter<datenum(2007,1,1));            %% leaf litter mass (Mg/ha/yr) 10^6 g / 10000 m2 365 day
LeafLitter(ind) = LeafLitter(ind).*10^6./(10000.*365);
AllLitter(ind)  = AllLitter(ind).*10^6./(10000.*365);
LeafLitterStd(ind) = LeafLitterStdIni(ind).*10^6./(10000.*365);
AllLitterStd(ind)  = AllLitterStdIni(ind).*10^6./(10000.*365);
FruitLitter(ind) = FruitLitter(ind).*10^6./(10000.*365);
WoodLitter(ind)  = WoodLitter(ind).*10^6./(10000.*365);
FruitLitterStd(ind) = FruitLitterStdIni(ind).*10^6./(10000.*365);
WoodLitterStd(ind)  = WoodLitterStdIni(ind).*10^6./(10000.*365);
MiscLitter(ind) = MiscLitter(ind).*10^6./(10000.*365);
MiscLitterStd(ind) = MiscLitterStdIni(ind).*10^6./(10000.*365);

% % % Cosme change units  kg/ha   traps are (0.25 m2) .*1000./(10000)
ind = find((dayLitter>datenum(2017,10,1))&(dayLitter<datenum(2018,9,1))); %% leaf litter mass (kg/ha) 10^3 g / 10000 m2
LeafLitter(ind) = LeafLitter(ind)./(0.25.*dT(ind));
AllLitter(ind)  = AllLitter(ind)./(0.25.*dT(ind));
LeafLitterStd(ind) = LeafLitterStd(ind)./(0.25.*dT(ind));
AllLitterStd(ind)  = AllLitterStd(ind)./(0.25.*dT(ind));
FruitLitter(ind) = FruitLitter(ind)./(0.25.*dT(ind));
WoodLitter(ind)  = WoodLitter(ind)./(0.25.*dT(ind));
FruitLitterStd(ind) = FruitLitterStd(ind)./(0.25.*dT(ind));
WoodLitterStd(ind)  = WoodLitterStd(ind)./(0.25.*dT(ind));
MiscLitter(ind) = MiscLitter(ind)./(0.25.*dT(ind));
MiscLitterStd(ind) = MiscLitterStd(ind)./(0.25.*dT(ind));

% % % % % Mauricio change units  kg/ha    traps are (0.22 m2) and values
% are in grams
% % ind = find((dayLitter>datenum(2011,1,1))&(dayLitter<datenum(2014,1,1))&(LeafLitter>30)); %% leaf litter mass (kg/ha) 10^3 g / 10000 m2
% % LeafLitter(ind) = NaN;
ind = find((dayLitter>datenum(2011,9,30))&(dayLitter<datenum(2013,7,2))); %% leaf litter mass (kg/ha) 10^3 g / 10000 m2
offset = 2.8;   offsetAll = 4.;   offsetMisc = 0.3;   offsetWood = .6;    offsetFruit = 0.3;
% offset = 0;   offsetAll = 0;   offsetMisc = 0;   offsetWood = 0.;    offsetFruit = 0.;
multiplier = 1;
LeafLitter(ind)    = multiplier.*(LeafLitter(ind)+offset);
LeafLitterStd(ind) = multiplier.*(LeafLitterStd(ind)+offset);
AllLitter(ind)     = multiplier.*(AllLitter(ind)+offsetAll);
AllLitterStd(ind)  = multiplier.*(AllLitterStd(ind)+offsetAll);
FruitLitter(ind)   = multiplier.*(FruitLitter(ind)+offsetFruit);
FruitLitterStd(ind)= multiplier.*(FruitLitterStd(ind)+offsetFruit);
WoodLitterStd(ind) = multiplier.*(WoodLitterStd(ind)+offsetWood);
WoodLitter(ind)    = multiplier.*(WoodLitter(ind)+offsetWood);
MiscLitter(ind)    = multiplier.*(MiscLitter(ind)+offsetMisc);
MiscLitterStd(ind) = multiplier.*(MiscLitterStd(ind)+offsetMisc);


multiplier = 1;%0.5;
ind = find((dayLitter>datenum(2011,1,1))&(dayLitter<datenum(2014,1,1))); %% leaf litter mass (kg/ha) 10^3 g / 10000 m2
LeafLitter(ind)    = (LeafLitter(ind))./(multiplier.*0.220619.*dT(ind));
LeafLitterStd(ind) = (LeafLitterStd(ind))./(multiplier.*0.220619.*dT(ind));
AllLitter(ind)     = (AllLitter(ind))./(multiplier.*0.220619.*dT(ind));
AllLitterStd(ind)  = (AllLitterStd(ind))./(multiplier.*0.220619.*dT(ind));
FruitLitter(ind)   = (FruitLitter(ind))./(multiplier.*0.220619.*dT(ind));
FruitLitterStd(ind)= (FruitLitterStd(ind))./(multiplier.*0.220619.*dT(ind));
WoodLitterStd(ind) = (WoodLitterStd(ind))./(multiplier.*0.220619.*dT(ind));
WoodLitter(ind)    = (WoodLitter(ind))./(multiplier.*0.220619.*dT(ind));
MiscLitter(ind)    = (MiscLitter(ind))./(multiplier.*0.220619.*dT(ind));
MiscLitterStd(ind) = (MiscLitterStd(ind))./(multiplier.*0.220619.*dT(ind));

% .........................................................................
% kg/ha   10^3 g / 10000 m2 >>> for the day use DeltaT   traps are (0.25 m2)
ind = find((dayLitter>datenum(2014,1,1))&(dayLitter<datenum(2017,10,1)));
LeafLitter(ind) = LeafLitter(ind).*10^3./(10000.*dT(ind));
AllLitter(ind)  = AllLitter(ind).*10^3./(10000.*dT(ind));
LeafLitterStd(ind) = LeafLitterStd(ind).*10^3./(10000.*dT(ind));
AllLitterStd(ind)  = AllLitterStd(ind).*10^3./(10000.*dT(ind));
FruitLitter(ind) = FruitLitter(ind).*10^3./(10000.*dT(ind));
WoodLitter(ind)  = WoodLitter(ind).*10^3./(10000.*dT(ind));
FruitLitterStd(ind) = FruitLitterStd(ind).*10^3./(10000.*dT(ind));
WoodLitterStd(ind)  = WoodLitterStd(ind).*10^3./(10000.*dT(ind));
MiscLitter(ind) = MiscLitter(ind).*10^3./(10000.*dT(ind));
MiscLitterStd(ind) = MiscLitterStd(ind).*10^3./(10000.*dT(ind));

ind = find((dayLitter>datenum(2018,09,1))&(dayLitter<datenum(2021,1,1)));
LeafLitter(ind) = LeafLitter(ind).*10^3./(10000.*dT(ind));
AllLitter(ind)  = AllLitter(ind).*10^3./(10000.*dT(ind));
LeafLitterStd(ind) = LeafLitterStd(ind).*10^3./(10000.*dT(ind));
AllLitterStd(ind)  = AllLitterStd(ind).*10^3./(10000.*dT(ind));
FruitLitter(ind) = FruitLitter(ind).*10^3./(10000.*dT(ind));
WoodLitter(ind)  = WoodLitter(ind).*10^3./(10000.*dT(ind));
FruitLitterStd(ind) = FruitLitterStd(ind).*10^3./(10000.*dT(ind));
WoodLitterStd(ind)  = WoodLitterStd(ind).*10^3./(10000.*dT(ind));
MiscLitter(ind) = MiscLitter(ind).*10^3./(10000.*dT(ind));
MiscLitterStd(ind) = MiscLitterStd(ind).*10^3./(10000.*dT(ind));

% LeafLitter(LeafLitter>4) = NaN;
% AllLitter(AllLitter>6) = NaN;

figure('color','white');        hold on;
plot(LeafLitter(dayLitter<datenum(2007,1,1)));
plot(LeafLitter((dayLitter>datenum(2011,1,1))&(dayLitter<datenum(2014,1,1))));
plot(LeafLitter((dayLitter>datenum(2014,1,1))&(dayLitter<datenum(2017,1,1))));
plot(LeafLitter((dayLitter>datenum(2017,1,1))&(dayLitter<datenum(2018,11,1))));
plot(LeafLitter((dayLitter>datenum(2018,11,1))&(dayLitter<datenum(2021,1,1))));
legend ('Harvard','Mauricio','GoAmazon','Cosme','Isotraits');      grid on;
ylabel('Litter_{leaf} (gC m^-^2 d^-^1)');

%% ........................................................................
% K67.dayLeafLitter  = interp1(dayLitter,dataLeafLitter,K67.day,'linear');
K67.dayLeafLitter  = K67.day.*NaN;       K67.dayAllLitter = K67.day.*NaN;
K67.dayLeafLitterStd  = K67.day.*NaN;    K67.dayAllLitterStd = K67.day.*NaN;
K67.dayFruitLitter  = K67.day.*NaN;      K67.dayWoodLitter = K67.day.*NaN;
K67.dayFruitLitterStd  = K67.day.*NaN;   K67.dayWoodLitterStd = K67.day.*NaN;
K67.dayMiscLitter  = K67.day.*NaN;       K67.dayMiscLitterStd = K67.day.*NaN;
for ik = 1:length(K67.day)
    ind = find(dayLitter==K67.day(ik));
    if ~isempty(ind)
        K67.dayLeafLitter(ik) = LeafLitter(ind);
        K67.dayLeafLitterStd(ik) = LeafLitterStd(ind);
        K67.dayAllLitter(ik)  = AllLitter(ind);
        K67.dayAllLitterStd(ik)  = AllLitterStd(ind);
        K67.dayFruitLitter(ik) = FruitLitter(ind);
        K67.dayFruitLitterStd(ik) = FruitLitterStd(ind);
        K67.dayWoodLitter(ik)  = WoodLitter(ind);
        K67.dayWoodLitterStd(ik)  = WoodLitterStd(ind);
        K67.dayMiscLitter(ik) = MiscLitter(ind);
        K67.dayMiscLitterStd(ik) = MiscLitterStd(ind);
    end
end

K67.weekLeafLitter     = AM_week2day_rs(K67.dayLeafLitter,K67.day);
K67.week_stdLeafLitter = AM_week2day_rs(K67.dayLeafLitterStd,K67.day);
K67.weekLeafLitter     = AM_fill_plus2(K67.weekLeafLitter);       K67.weekLeafLitter = AM_fill_plus2(K67.weekLeafLitter);
K67.week_stdLeafLitter = AM_fill_plus2(K67.week_stdLeafLitter);   K67.week_stdLeafLitter = AM_fill_plus2(K67.week_stdLeafLitter);

K67.monthLeafLitter    = AM_month(K67.dayLeafLitter,K67.day);
K67.monthLeafLitter    = AM_fill_plus2(K67.monthLeafLitter);      K67.monthLeafLitter = AM_fill_plus2(K67.monthLeafLitter);

[K67.week_avgLeafLitter,~,K67.week_avgstdLeafLitter] = AM_week2_avg(K67.weekLeafLitter,K67.week);

ind = find((K67.day>datenum(2011,1,1))&(K67.day<datenum(2014,1,1))); %% leaf litter mass (kg/ha) 10^3 g / 10000 m2
ix = K67.dayLeafLitter;     ix(ind) = NaN;
iy = AM_week2day_rs(ix,K67.day);        [ix] = AM_week2_avg(iy,K67.week);
K67.weekLeafLitter_recurrent = AM_week2dayAVG_recurrent(ix,K67.week);

figure('color','white');
plot(K67.week,K67.weekLeafLitter_recurrent,'color',[.6 .6 0.6]);    hold on
plot(K67.week,K67.weekLeafLitter,'b');      hold on
plot(K67.week,K67.weekLeafLitter,'.r');  datetick('x')
grid on;                                ylabel('Litter_{leaf} (gC m^-^2 d^-^1)')

figure('color','white');
plot(K67.month,K67.monthLeafLitter.*365.*10000./10^6);      hold on
grid on;                                ylabel('Litter_{leaf} (MgC ha^-^2 year^-^1)')
ind = ((K67.month>datenum(2011,1,1))&(K67.month<datenum(2014,1,1)));
display(K67.monthLeafLitter(ind).*365.*10000./10^6)

% .........................................................................
% iy = AM_fill_composite(K67.weekLeafLitter,K67.week,23);
% K67.weekLeafLitter(~isnan(ix)) = iy(~isnan(ix));
% x = smooth(K67.weekLeafLitter);    x(isnan(K67.weekLeafLitter)) = NaN;
% K67.weekLeafLitter = x;

K67.weekAllLitter     = AM_week2day_rs(K67.dayAllLitter,K67.day);
K67.week_stdAllLitter = AM_week2day_rs(K67.dayAllLitterStd,K67.day);
K67.weekAllLitter     = AM_fill_plus2(K67.weekAllLitter);       K67.weekAllLitter = AM_fill_plus2(K67.weekAllLitter);
K67.week_stdAllLitter = AM_fill_plus2(K67.week_stdAllLitter);   K67.week_stdAllLitter = AM_fill_plus2(K67.week_stdAllLitter);
[K67.week_avgAllLitter,~,K67.week_avgstdAllLitter] = AM_week2_avg(K67.weekAllLitter,K67.week);

ind = find((K67.day>datenum(2011,1,1))&(K67.day<datenum(2014,1,1))); %% leaf litter mass (kg/ha) 10^3 g / 10000 m2
ix = K67.dayAllLitter;     ix(ind) = NaN;
iy = AM_week2day_rs(ix,K67.day);        [ix] = AM_week2_avg(iy,K67.week);
K67.weekAllLitter_recurrent = AM_week2dayAVG_recurrent(ix,K67.week);

figure('color','white');
plot(K67.week,K67.weekAllLitter_recurrent,'color',[.6 .6 0.6]);    hold on
plot(K67.week,K67.weekAllLitter,'b');      hold on
plot(K67.week,K67.weekAllLitter,'.r');  datetick('x')
grid on;                               ylabel('Litter_{all} (gC m^-^2 d^-^1)')

% .........................................................................
K67.weekFruitLitter     = AM_week2day_rs(K67.dayFruitLitter,K67.day);
K67.week_stdFruitLitter = AM_week2day_rs(K67.dayFruitLitterStd,K67.day);
K67.weekFruitLitter     = AM_fill_plus2(K67.weekFruitLitter);       K67.weekFruitLitter = AM_fill_plus2(K67.weekFruitLitter);
K67.week_stdFruitLitter = AM_fill_plus2(K67.week_stdFruitLitter);   K67.week_stdFruitLitter = AM_fill_plus2(K67.week_stdFruitLitter);
[K67.week_avgFruitLitter,~,K67.week_avgstdFruitLitter] = AM_week2_avg(K67.weekFruitLitter,K67.week);

ind = find((K67.day>datenum(2011,1,1))&(K67.day<datenum(2014,1,1))); %% Fruit litter mass (kg/ha) 10^3 g / 10000 m2
ix = K67.dayFruitLitter;     ix(ind) = NaN;
iy = AM_week2day_rs(ix,K67.day);        [ix] = AM_week2_avg(iy,K67.week);
K67.weekFruitLitter_recurrent = AM_week2dayAVG_recurrent(ix,K67.week);

figure('color','white');
plot(K67.week,K67.weekFruitLitter_recurrent,'color',[.6 .6 0.6]);    hold on
plot(K67.week,K67.weekFruitLitter,'b');      hold on
plot(K67.week,K67.weekFruitLitter,'.r');  datetick('x')
grid on;                                ylabel('Litter_{Fruit} (gC m^-^2 d^-^1)')

% .........................................................................
%% Litter Oxford 
% .........................................................................
opts = delimitedTextImportOptions("NumVariables", 9);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["TransAF8Cod", "PreAF8ENAF8disturbance", "ENAF8impact", "Trap", "TrapAF8Code", "Date", "Month", "Year", "LeafAF8monthAF8MgAF8ha"];
opts.VariableTypes = ["categorical", "categorical", "categorical", "double", "categorical", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["TransAF8Cod", "PreAF8ENAF8disturbance", "ENAF8impact", "TrapAF8Code"], "EmptyFieldRule", "auto");
LitterOxford = readtable([FolderK67litter '210611 Leaf litter Oxford undisturbed.csv'], opts);
clear opts
% % Date", "Month", "Year"
dayLitterOxford = datenum(LitterOxford{:,"Year"},LitterOxford{:,"Month"},LitterOxford{:,"Date"});
dayLitterOxford = dayLitterOxford+15;

LeafLitterOxford = LitterOxford{:,"LeafAF8monthAF8MgAF8ha"};          LeafLitterOxford = LeafLitterOxford.*0.5;
LeafLitterOxford = LeafLitterOxford.*10^6./(10000);      % Mg/ha/month
dT = [31;diff(unique(dayLitterOxford))];

% ........................................................................
% from basket to forest
% ........................................................................
ix = LeafLitterOxford;                  iy = dayLitterOxford;
[ah,~,ch] = unique(iy);                 out_size = (histcounts(iy,ah)); out_size = [1,out_size]';
out_hour = [ah,accumarray(ch,ix,[],@nanmean)];
out_hourstd = [ah,accumarray(ch,ix,[],@nanstd)];
ix = out_hour(:,2);                     iy = out_hour(:,1);
in = out_hourstd(:,2)./sqrt(out_size);
LeafLitterOxfordDay = ix./dT;               DayLeafLitterOxfordStd = in./dT;

% % iw = dT;
% % out_hour = [ah,accumarray(ch,iw,[],@nanmean)];
% % iw = out_hour(:,2);
% % dT = iw;        dayLitterOxford  = iy;
dayLitterOxfordDate = unique(dayLitterOxford);

K67.dayLeafLitterOxford  = interp1(dayLitterOxfordDate,LeafLitterOxfordDay,K67.day);
K67.weekLeafLitterOxford = AM_week2day_rs(K67.dayLeafLitterOxford,K67.day);

figure('color','white');
hl3 = plot(K67.week,K67.weekLeafLitter_recurrent,'color',[.6 .6 0.6]);    hold on
hl1 = plot(K67.week,K67.weekLeafLitter,'b');              plot(K67.week,K67.weekLeafLitter,'.r');  
plot(dayLitterOxfordDate,LeafLitterOxfordDay,'.','color',[0 0.6 0]); 
hl2 = plot(K67.week,K67.weekLeafLitterOxford,'-g');
datetick('x');              grid on;                ylabel('Litter_{leaf} (gC m^-^2 d^-^1)')
xtickangle(0)
legend([hl1,hl3,hl2],'UofA','seasonal average','undisturbed Madalena et al. 2020');
legend box off;             xlim([K67.week(1) K67.week(end)]);

% .........................................................................
%% Wood 
% .........................................................................
K67.weekWoodLitter     = AM_week2day_rs(K67.dayWoodLitter,K67.day);
K67.week_stdWoodLitter = AM_week2day_rs(K67.dayWoodLitterStd,K67.day);
K67.weekWoodLitter     = AM_fill_plus2(K67.weekWoodLitter);       K67.weekWoodLitter = AM_fill_plus2(K67.weekWoodLitter);
K67.week_stdWoodLitter = AM_fill_plus2(K67.week_stdWoodLitter);   K67.week_stdWoodLitter = AM_fill_plus2(K67.week_stdWoodLitter);
[K67.week_avgWoodLitter,~,K67.week_avgstdWoodLitter] = AM_week2_avg(K67.weekWoodLitter,K67.week);

ind = find((K67.day>datenum(2011,1,1))&(K67.day<datenum(2014,1,1))); %% Wood litter mass (kg/ha) 10^3 g / 10000 m2
ix = K67.dayWoodLitter;     ix(ind) = NaN;
iy = AM_week2day_rs(ix,K67.day);        [ix] = AM_week2_avg(iy,K67.week);
K67.weekWoodLitter_recurrent = AM_week2dayAVG_recurrent(ix,K67.week);

figure('color','white');
plot(K67.week,K67.weekWoodLitter_recurrent,'color',[.6 .6 0.6]);    hold on
plot(K67.week,K67.weekWoodLitter,'b');      hold on
plot(K67.week,K67.weekWoodLitter,'.r');  datetick('x')
grid on;                                ylabel('Litter_{Wood} (gC m^-^2 d^-^1)')

% .........................................................................
K67.weekMiscLitter     = AM_week2day_rs(K67.dayMiscLitter,K67.day);
K67.week_stdMiscLitter = AM_week2day_rs(K67.dayMiscLitterStd,K67.day);
K67.weekMiscLitter     = AM_fill_plus2(K67.weekMiscLitter);       K67.weekMiscLitter = AM_fill_plus2(K67.weekMiscLitter);
K67.week_stdMiscLitter = AM_fill_plus2(K67.week_stdMiscLitter);   K67.week_stdMiscLitter = AM_fill_plus2(K67.week_stdMiscLitter);
[K67.week_avgMiscLitter,~,K67.week_avgstdMiscLitter] = AM_week2_avg(K67.weekMiscLitter,K67.week);

ind = find((K67.day>datenum(2011,1,1))&(K67.day<datenum(2014,1,1))); %% Misc litter mass (kg/ha) 10^3 g / 10000 m2
ix = K67.dayMiscLitter;     ix(ind) = NaN;
iy = AM_week2day_rs(ix,K67.day);        [ix] = AM_week2_avg(iy,K67.week);
K67.weekMiscLitter_recurrent = AM_week2dayAVG_recurrent(ix,K67.week);

figure('color','white');
plot(K67.week,K67.weekMiscLitter_recurrent,'color',[.6 .6 0.6]);    hold on
plot(K67.week,K67.weekMiscLitter,'b');      hold on
plot(K67.week,K67.weekMiscLitter,'.r');  datetick('x')
grid on;                                ylabel('Litter_{Misc} (gC m^-^2 d^-^1)')

%% ....................................................................
v1 = 0;         v2 = 3.4;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgLeafLitter-K67.week_avgstdLeafLitter,2.*K67.week_avgstdLeafLitter];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgLeafLitter([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 =K67.weekLeafLitter plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl3 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2014),'Color',[0.8 .8 0],'LineWidth',2);
hl4 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
hl5 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2016),'Color',[0.5 0 .5],'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'NPP_{litter_{leaf}}';'(gC m^-^2 d^-^1)'},'FontSize',11);        xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl3,hl4,hl5],'Jul2008-Jun2009','Jul2014-Jun2015','Jul2015-Jun2016','Jul2016-Jun2017');

%% ....................................................................
v1 = 0;         v2 = 36;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgLeafLitter-K67.week_avgstdLeafLitter,2.*K67.week_avgstdLeafLitter];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,cumsum(y));       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,cumsum(K67.weekLeafLitter(K67.HY_weekY == ik)),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,cumsum(K67.week_avgLeafLitter([week_start_hydro:23,1:(week_start_hydro-1)])),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 =K67.weekLeafLitter plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,cumsum(K67.weekLeafLitter(K67.HY_weekY == 2008)),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl3 = plot(week_avg,cumsum(K67.weekLeafLitter(K67.HY_weekY == 2014)),'Color',[0.8 .8 0],'LineWidth',2);
hl4 = plot(week_avg,cumsum(K67.weekLeafLitter(K67.HY_weekY == 2015)),'Color',color_mtx(2,:),'LineWidth',2);
hl5 = plot(week_avg,cumsum(K67.weekLeafLitter(K67.HY_weekY == 2016)),'Color',[0.5 0 .5],'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'Litter_{leaf}';'(gC m^-^2 d^-^1)'},'FontSize',11);        xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl3,hl4,hl5],'Jul2008-Jun2009','Jul2014-Jun2015','Jul2015-Jun2016','Jul2016-Jun2017');

%% ....................................................................
v1 = 0.6;         v2 = 3.4;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgLeafLitter-K67.week_avgstdLeafLitter,2.*K67.week_avgstdLeafLitter];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgLeafLitter([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekLeafLitter(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);        xtickangle(0);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'Litter_{leaf}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

% .........................................................................
v1 = 0.6;         v2 = 3.4;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgAllLitter-K67.week_avgstdAllLitter,2.*K67.week_avgstdAllLitter];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekAllLitter(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgAllLitter([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekAllLitter(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);        xtickangle(0);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'Litter_{all}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;

% ........................................................................
% .........................................................................
%%  DBH
% .........................................................................
% ........................................................................
dayDBH = [196	345	414	445	489	551	592	634	678	702	775	800	838	873	914	950	992	1029	1075	1103	1129	1159	1195	1230	1271	1299	1327	1356	1385	1412	1446	1480	1506	1547	1580	1632	1657	1709	1733	1762	1799	1825	1855	1888	1982	2063	2099	2127	2168	2197	2280];
dayDBH = (datenum(1999,1,0)+dayDBH)';
DBH = load([FolderK67dendrometry 'lba_km67_dend.dat']);
DBHid = DBH(:,1);         DBH(:,1:3) = [];    DBH = DBH';     DBH(DBH<11) = NaN;            %only DBH >10 cm
x = (-1.06+(0.333.*(log(DBH)))+(0.933.*((log(DBH)).^2))+(-0.122.*((log(DBH)).^3)));
TreeMass = exp(x);
TreeMass = TreeMass.*0.5./1000;             %MgC tree-1
disp('standing biomass')
x = 1:length(dayDBH); 	x = find(dayDBH>(datenum(2002,1,1)));     x = nanmin(x);
y = TreeMass(x,:);        y(isnan(y)) = [];                       disp((nanmean(y))*470)   %tree ha-1
% at Rice et al. 2001 Live Biomass is 147.4 MgC ha-1
TreeMass = TreeMass.*(10^6);                        %gC tree-1
TreeMass_delta = diff(TreeMass);                    dayDBH_delta = diff(dayDBH);
TreeMass_delta(TreeMass_delta<-5000) = NaN;
Wood = dayDBH_delta.*NaN;
for ik = 1:length(dayDBH_delta)
    x = TreeMass_delta(ik,:); x(isnan(x)) = [];
    Wood(ik) = nanmean(x);
end
TreeDensity       = 470/10000;                             %tree m-2
K67.WoodIncrement = TreeDensity.*Wood./dayDBH_delta;       %gC m-2 d-1
K67.WoodIncrement = [NaN;K67.WoodIncrement];
K67.weekWood      = interp1(dayDBH,K67.WoodIncrement,K67.week,'linear');

[K67.week_avgWood,~,K67.week_avgstdWood]  = AM_week2_avg(K67.weekWood,K67.week);

[K67.yearWood,~,~,K67.yearHY] = AM_year(K67.weekWood,K67.week+datenum(0,0,365/2));       %gC m-2 d-1
K67.yearWood      = K67.yearWood.*365./1000;

%
v1 = 0.2;         v2 = 0.5;  v3 = v1:((v2-v1)/3):v2;
figure('color','white');
subplot(3,1,1);             hold on;
color_yyplot = [.2,.2,.2;.2 ,.7,.2];
colororder(color_yyplot)
yyaxis left
bar(K67.yearHY,K67.yearWood,'FaceColor',[.7 .7 .7]);
iy = [NaN;diff(K67.yearWood)];      hold on;
set(gca,'XLim',[2002 2020],'XTickLabel',[K67.yearHY],'XLim',[K67.yearHY(1) 2020],...
    'YTick',v3,'YLim',[v1 v2],'FontSize',11);
ylabel('Wood (kg m^-^2 y^-^1)','FontSize',12,'Color',[.2 .2 .2]);
v1 = -0.2;         v2 = 0.2;  v3 = v1:((v2-v1)/2):v2;
yyaxis right
plot(K67.yearHY,iy,'color',[0 .7 .0],'LineWidth',2);
set(gca,'XLim',[2002 2020],'XTickLabel',[K67.year],'XLim',[K67.year(1) 2020],...
    'YTick',v3,'YLim',[v1 v2],'FontSize',11);
box on;
ylabel({'Wood growth';'(kg m^-^2 y^-^1)'},'FontSize',12,'Color',[.2 .7 .2]);
grid on;

%% ........................................................................
v1 = 0.;         v2 = 3.;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgWood-K67.week_avgstdWood,2.*K67.week_avgstdWood];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekWood(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    % hl3  =  plot(month_avg,TRMM.anomaly_precip2009,'Color',[.0 0 .9],'LineWidth',2); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgWood([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
hl7 = plot(week_avg,K67.weekWood(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekWood(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'NPP_{wood}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

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
ylabel ({'LAI'},'FontSize',11);       box on;
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
ylabel ({'LAI_{normalized}'},'FontSize',11);       box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');


% ........................................................................
% .........................................................................
%% Flush
% .........................................................................
% ........................................................................
K67.weekDt = [16;diff(K67.week-K67.week(1))];

% K67.weekPc_smooth = smooth(K67.weekPc,3);         K67.weekPc_smooth(isnan(K67.weekPc)) = NaN;
K67.weekPc_smooth = K67.weekPc;
K67.weekPc_smooth = AM_fill_plus(K67.weekPc_smooth);

K67.weekDtPc = diff(K67.weekPc_smooth);         K67.weekDtPc = [NaN;K67.weekDtPc];
% K67.weekDtPc = diff(K67.weekPc_smooth);         K67.weekDtPc = [NaN;NaN;NaN;NaN;K67.weekDtPc(1:end-4)];
% iw = smooth(K67.weekDtPc);                      iw(isnan(K67.weekDtPc)) = NaN;
iw = K67.weekDtPc;                              iw = AM_rm_outlier(iw,3);       iw = AM_fill(iw);
K67.weekDtPc = iw./K67.weekDt;                  K67.weekDtPc = AM_rm_outlier(K67.weekDtPc,3);
ix = smooth(K67.weekDtPc);      ix(isnan(K67.weekDtPc)) = NaN;      K67.weekDtPc = ix;
%%
% % Old SLA = 1/(53.3 + (0.95*100));                  SLA = SLA./.48;     %SLA = SLA.*1.5;
% % Old SLA = 1/(40 + (0.95*100));                    SLA = SLA./.48;     %SLA = SLA.*1.5;
% % Old K67.week_avgSLA = ones(23,1).*SLA;

% % LMA (or its inverse, specific leaf area, SLA)
% % Wu et al. 2017
% % Our 759 leaves from 11 trees in Brazil encompassed variation inleaf mass per area (LMA; 35–270 g m-2)
% % Convergence in relationships between leaf traits, spectra and age across diverse canopy environments and two contrasting tropical forests
% % From figure 2a  logaritmic regression leaf age g m-2 and age in days
xLeafAge_yLMA = [21.393,57.237; 31.841,70.066;...
    44.776,81.908; 61.692,93.75; 78.607,103.618; 103.483,110.526; 124.378,117.434;...
    147.264,125.329; 170.149,128.289; 195.025,133.224; 219.9,137.171; 251.741,140.132;...
    283.582,145.066; 309.453,148.026; 336.318,151.974; 366.169,153.947; 397.512,157.895];
x =xLeafAge_yLMA(:,1);      y = xLeafAge_yLMA(:,2);
% y = x*loga(b);
modelfun = @(b,x)(b(1).*log(x) + b(2));
beta = nlinfit(x,xLeafAge_yLMA(:,2),modelfun,[0;0]);
figure('color','white');            plot(x,y,'.','MarkerSize',10);  hold on;
w = beta(1).*log(x)+beta(2);        plot(x,w,'-','LineWidth',2);
ylabel('LMA (g m^-^2)');            xlabel('Age (d)');
set(gca,'XTick',[0:100:400],'XTickLabel',[0:100:400],'Xlim',[0 400],...
    'YTick',[50:50:200],'YTickLabel',[50:50:200],'Ylim',[50 200],'FontSize',12);

% The GCC from the SOC indicates young leaves are dominant in July - mature
% in January and old in June
% where x is date and y is the age of leaf
x = (datenum(0,7,1):datenum(0,7,1)+365)';     y = (50:50+365)';
K67.week_avgLeafAge = AM_week2_avg(y,x);
K67.week_avgLMA = beta(1).*log(K67.week_avgLeafAge)+beta(2);
K67.weekLMA     = AM_week2dayAVG_recurrent(K67.week_avgLMA,K67.week);
% K67.weekLMA = K67.weekLMA.*0.8;

umol2gC = 60.*60.*24.*12.*(10.^-6);
Amax = (2.97.*((0.014.*(100))+1.4)) + 1.18;      %at height 100 Amax = 15;%  umol CO2 m-2 s-1 (max)
K67.week_avgAmax = ones(23,1).*umol2gC.*Amax;

Amax80_decrease = [.95, .9, .85, .8];
Amax60_decrease = [.93, .86, .8, .73, .66, .6];
Amax40_decrease = [.925, .85, .775, .7, .625, .55, .475, .4];

K67.weekAmax = [];                 
for ik = 1:length(K67.year)
    %     K67.weekSLA  = [K67.weekSLA;K67.week_avgSLA];
    K67.weekAmax = [K67.weekAmax;K67.week_avgAmax];
end

K67.weekLeafLitterFill = K67.weekLeafLitter;
ind = find((K67.week>datenum(2016,5,20))&(K67.week<datenum(2019,5,20)));
K67.weekLeafLitterFill(ind) = K67.weekLeafLitterOxford(ind);
plot(K67.week,K67.weekLeafLitterFill,'k');             
[K67.week_avgLeafLitterFill,~,K67.week_avgstdLeafLitterFill] = AM_week2_avg(K67.weekLeafLitterFill,K67.week);

K67.weekLeafLitterFill = AM_fill_month_cycle(K67.weekLeafLitterFill,K67.week,6,18);
K67.weekLeafLitterFill_smooth = smooth(K67.weekLeafLitterFill,3);

Amax40_week = K67.weekAmax;         Amax80_week = K67.weekAmax;      Amax60_week = K67.weekAmax;

% ix = smooth(K67.weekLeafLitter);         ix(isnan(K67.weekLeafLitter)) = NaN;

for ik = 1:length(K67.year)
    ind = find(K67.Yweek == K67.year(ik));
    ix = K67.weekLeafLitterFill(ind);   iy = K67.week(ind);     %ix(isnan(ix)) = K67.week_avgLeafLitter(isnan(ix));
    iw = find(ix == nanmax(ix));
    if ~isempty(iy(iw))
        iz = iy(iw(1));        ind = find(K67.week==iz);
        Amax80_week(ind-4:ind-1) = K67.weekAmax(ind-4:ind-1).*Amax80_decrease';
        Amax60_week(ind-6:ind-1) = K67.weekAmax(ind-6:ind-1).*Amax60_decrease';
        Amax40_week(ind-8:ind-1) = K67.weekAmax(ind-8:ind-1).*Amax40_decrease';
    end
end
close;
Amax80_week = smooth(Amax80_week,3);      Amax60_week = smooth(Amax60_week,3);
Amax40_week = smooth(Amax40_week,3);

K67.weekFlush = (K67.weekDtPc.*K67.weekLMA.*1./Amax80_week)+K67.weekLeafLitterFill;
K67.weekFlush = AM_rm_outlier(K67.weekFlush,3);

K67.weekFlush(K67.weekFlush<0) = 0;
ix = smooth(K67.weekFlush,5);     ix(isnan(K67.weekFlush)) = NaN;       K67.weekFlush = ix;

[K67.week_avgFlush,~,K67.week_avgstdFlush] = AM_week2_avg(K67.weekFlush,K67.week);
K67.monthFlush = interp1(K67.week,K67.weekFlush,K67.month);

[K67.week_avgDtPc,~,K67.week_avgstdDtPc] = AM_week2_avg(K67.weekDtPc,K67.week);

%             0.0032 for K67 and K83 (Domingues et al., 2005) and to 0.0038 m2 gC-1 for K34 (Carswell et al., 2000)
%             Delta LAIxSLA
% LAI change as leaf flux driver ..........................................
% K67.weekDtLAI = diff(K67.weekLAI);      K67.weekDtLAI = [NaN;NaN;NaN;NaN;K67.weekDtLAI(1:end-4)];
ix = K67.weekLAIfparDG;                   ix(isempty(K67.weekLAIfparDG)) = K67.weekLAI(isempty(K67.weekLAIfparDG));
K67.weekDtLAI = [NaN;diff(ix)];
K67.weekDtLAI = K67.weekDtLAI./K67.weekDt;
ix = smooth(K67.weekDtLAI);           ix(isnan(K67.weekDtLAI)) = NaN;        K67.weekDtLAI = ix;
K67.week_avgDtLAI = AM_week2_avg(K67.weekDtLAI,K67.week);
K67.weekDtLAIleafFlush = (K67.weekDtLAI.*K67.weekLMA.*1.)+K67.weekLeafLitterFill;
K67.weekDtLAIleafFlush(K67.weekDtLAIleafFlush<0) = 0;

[K67.week_avgDtLAIleafFlush,~,K67.week_avgstdDtLAIleafFlush] = AM_week2_avg(K67.weekDtLAIleafFlush,K67.week);

%%
ix = smooth(K67.weekFlush,5);     ix(isnan(K67.weekFlush)) = NaN;
figure('color','white');
plot(K67.week,K67.weekFlush); hold on;      plot(K67.week,ix);
plot(K67.week,K67.weekDtLAIleafFlush);
plot(K67.week,(K67.weekDtPc.*1.0.*K67.weekLMA./Amax40_week)); plot(K67.week,K67.weekLeafLitter);
legend('Flush_{fn Pc}','Flush_{fn Pc smooth}','Flush_{dLAI/dt}','dPc/dt LMA/A_{max}','Litter');
%%
% % SLA=1/(53.3 + (0.95*100));      SLA=SLA./.48;   SLA_week_avg=ones(23,1).*SLA;
% % umol2gC=60.*60.*24.*12.*(10.^-6);
% % Amax=(2.97.*((0.014.*(80))+1.4)) + 1.18;    %Amax=15;%umol CO2 m-2 s-1 (max)
% % Amax_week_avg=ones(23,1).*umol2gC.*Amax;
% % %find the avg annual peak in litterfall and then do the linear decrease to reach 80% 60% and 40% at that time
% % Leaf_week_avg=AM_week2_avg(Leaf_week);
% % Amax80_decrease=[.95, .9, .85, .8];
% % Amax60_decrease=[.9, .80, .78, .72, .66, .6];
% % Amax40_decrease=[.96, .88, .8, .72, .64, .56, .48, .4];
% % [Yweek,~,~]=datevec(week);                  year=unique(Yweek);
% % Amax_week=[];                               SLA_week=[];
% % for ik=1:length(year)
% %     SLA_week=[SLA_week;SLA_week_avg];
% %     Amax_week=[Amax_week;Amax_week_avg];
% % end
% % Amax40_week=Amax_week;         Amax80_week=Amax_week;      Amax60_week=Amax_week;
% % for ik=1:length(year)
% %     ind=find(Yweek==year(ik));
% %     ix=Leaf_week(ind);          iy=week(ind);       ix(isnan(ix))=Leaf_week_avg(isnan(ix));
% %     iw=find(ix==nanmax(ix));    iz=iy(iw);          ind=find(week==iz);
% %     Amax80_week(ind-4:ind-1)=Amax_week(ind-4:ind-1).*Amax80_decrease';
% %     Amax60_week(ind-6:ind-1)=Amax_week(ind-6:ind-1).*Amax60_decrease';
% %     Amax40_week(ind-8:ind-1)=Amax_week(ind-8:ind-1).*Amax40_decrease';
% % end
% % Amax80_week=smooth(Amax80_week);Amax60_week=smooth(Amax60_week);Amax40_week=smooth(Amax40_week);
% % LAI_week_dt=diff(LAI_week);                    LAI_week_dt=[NaN;LAI_week_dt];
% % LAI_week_dt=LAI_week_dt./week_dt;              x=smooth(LAI_week_dt); x(isnan(LAI_week_dt))=NaN;    LAI_week_dt=x;%./SLA_week;
% % LAI_week_dt_avg=AM_week2_avg(LAI_week_dt);
% % Pc_week_dt=Pc_week_dt;%./(SLA_week.*Amax80_week);
% % Flush_week=(Pc_week_dt./(Amax60_week.*SLA_week))+Leaf_week;
% % x=smooth(Flush_week);     x(isnan(Flush_week))=NaN;   Flush_week=x;

% %             %             0.0032 for K67 and K83 (Domingues et al., 2005) and to 0.0038 m2 gC-1 for K34 (Carswell et al., 2000)
% %             %             Delta LAIxSLA
% %             hl4 = line(week_avg,LAI_week_dt_avg./SLA_week_avg,'Color',[0.35 0.35 0.35],'LineStyle','--','LineWidth',2,'Parent',ax2);

%% ........................................................................
v1 = 0;         v2 = 4.8;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgFlush-K67.week_avgstdFlush,2.*K67.week_avgstdFlush];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekFlush(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgFlush([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekFlush(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekFlush(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekFlush(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekFlush(K67.HY_weekY == 2014),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'NPP_{leaf-flush}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');  legend box off;


%% ........................................................................
v1 = 0;         v2 = 4.0;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgDtLAIleafFlush-K67.week_avgstdDtLAIleafFlush,2.*K67.week_avgstdDtLAIleafFlush];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekDtLAIleafFlush(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgDtLAIleafFlush([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekDtLAIleafFlush(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekFlush(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekDtLAIleafFlush(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekDtLAIleafFlush(K67.HY_weekY == 2014),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'NPP_{leaf-flush LAI}';'(gC m^-^2 d^-^1)'},'FontSize',11);       box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);        xtickangle(0);
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');  legend box off;

%%
v1 = -.4;         v2 = .4;  v3 = v1:((v2-v1)/3):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgDtPc-K67.week_avgstdDtPc,2.*K67.week_avgstdDtPc];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekDtPc(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgDtPc([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekDtPc(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekDtPc(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekDtPc(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekDtPc(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'DtPc (gC m^-^2)'},'FontSize',11);             box on;
text (week_avg(2),v2-1,site_name,'FontSize',11);        xtickangle(0);
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');  legend box off;

%% ........................................................................
v1 = 6;         v2 = 28;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgPc-K67.week_avgstdPc,2.*K67.week_avgstdPc];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2019
    hl2  =  plot(week_avg,K67.weekPc(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgPc([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekPc(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekPc(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekPc(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekPc(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'Pc';'(gC m^-^2 d^-^1)'},'FontSize',11);       xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

%% ........................................................................
v1 = 0;         v2 = 5;  v3 = v1:((v2-v1)/2):v2;
figure('color','white');
subplot(3,1,1);             hold on;
x4=0.*month_avg;            x4((K67.month_avgTRMM<105))=(v2);
hl1 = bar(month_avg+15,x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
hl2 = bar(month_avg+15,-x4([month_start_hydro:12,1:month_start_hydro-1]),1,'FaceColor',[ 0.9 0.9 0.9],'EdgeColor',[ 0.9 0.9 0.9]);
y = [K67.week_avgWUE-K67.week_avgstdWUE,2.*K67.week_avgstdWUE];
y = y([week_start_hydro:23,1:(week_start_hydro-1)],:);
h = area(week_avg,y);       set(gca,'Layer','top');
set(h(2),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',0);
ik = 2002;
while ik<2020
    hl2  =  plot(week_avg,K67.weekWUEdryday(K67.HY_weekY == ik),'Color',[.7 .7 .7],'LineWidth',1); hold on;
    ik = ik+1;
end
hl1 = plot(week_avg,K67.week_avgWUE([week_start_hydro:23,1:(week_start_hydro-1)]),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
%     hl6 = plot(week_avg,K67.weekGsmmol(K67.HY_weekY == 2005),'Color',color_mtx(4,:),'LineWidth',2);
hl7 = plot(week_avg,K67.weekWUEdryday(K67.HY_weekY == 2008),'Color',color_mtx(5,:),'LineWidth',2);
%     hl3 = plot(week_avg,K67.weekWUEdryday(K67.HY_weekY == 2008),'Color',color_mtx(1,:),'LineWidth',2);
hl4 = plot(week_avg,K67.weekWUEdryday(K67.HY_weekY == 2015),'Color',color_mtx(2,:),'LineWidth',2);
%     hl5 = plot(week_avg,K67.weekWUEdryday(K67.HY_weekY == 2016),'Color',color_mtx(3,:),'LineWidth',2);
xlim([week_avg(1) 365]);    ylim([v1 v2]);
set(gca,'XTick',month_avg,'XTickLabel',[ '' '' ''],'YTick',v3,'XTickLabel',datestr(HYmonth_avg,4),'FontSize',11);
ax1 = gca;      set(ax1,'XColor','k','YColor',[ 0.0 0.0 0.0],'Position',[0.17 0.7 0.28 0.20]);
ylabel ({'WUE';'(gC mm^-^1)'},'FontSize',11);           xtickangle(0);
text (week_avg(2),v2-1,site_name,'FontSize',11);        box on;
%     legend([hl3,hl4,hl5],'Jul2009-Jun2008','Jul2015-Jun2016','Jul2016-Jun2017');
legend([hl7,hl4],'Jul2008-Jun2009','Jul2015-Jun2016');

% ........................................................................
% ........................................................................
%% Mortality as biomass loss
% ........................................................................
% ........................................................................
opts = delimitedTextImportOptions("NumVariables", 15);
opts.DataLines = [3, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["YYYY" "avg", "max", "min"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
mortality = readtable([FolderK67mortality 'K67mortality.biomass.csv'], opts);
clear opts

mortalityY = mortality{:,1};
mortalityAvg = mortality{:,2};     mortalityMax = mortality{:,3};    mortalityMin = mortality{:,4};

K67.year_avgMortality = interp1(mortalityY,mortalityAvg,K67.year);
K67.year_maxMortality = interp1(mortalityY,mortalityMax,K67.year);
K67.year_minMortality = interp1(mortalityY,mortalityMin,K67.year);
K67.year_dateloc = datenum(K67.year,1,1);
%%
xaxis = [];                 yaxis = [];
for ix = 1:length(K67.year)
    xaxis = [xaxis;(K67.year(ix).*ones(15,1))];
    yaxis = [yaxis;K67.year_minMortality(ix);K67.year_maxMortality(ix);...
        K67.year_minMortality(ix);K67.year_maxMortality(ix);...
        K67.year_minMortality(ix);K67.year_maxMortality(ix);...
        0.1+K67.year_avgMortality(ix);-0.1+K67.year_avgMortality(ix);K67.year_avgMortality(ix);...
        0.1+K67.year_avgMortality(ix);-0.1+K67.year_avgMortality(ix);K67.year_avgMortality(ix);...
        0.1+K67.year_avgMortality(ix);-0.1+K67.year_avgMortality(ix);K67.year_avgMortality(ix)];
end

%
figure('color','white','Position',[50 50 600 400],'units','inches');
boxplot(yaxis,xaxis);
hold on;
plot(K67.year_dateloc,(K67.year_dateloc.*0)+nanmean(K67.year_avgMortality),':','LineWidth',1,'color',[0.4 0.4 0.6]);
hold off;
ylabel ('Mortality (%AGB year^-^1)','FontSize',14);
%%
figure('color','white');
for ix = 1:length(K67.year)
    low = K67.year_minMortality(ix);
    high = (K67.year_maxMortality(ix)-K67.year_minMortality(ix));
    y = [low,high;low,high];
    x = [K67.year(ix)-.2,K67.year(ix)+0.2];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.2 .2 .2]);
    set(h(1),'FaceColor','none','EdgeColor','none');
    set(h,'BaseValue',0);   hold on;
    
    plot([K67.year(ix)-.2;K67.year(ix)+.2],[K67.year_avgMortality(ix);K67.year_avgMortality(ix)],'-k','MarkerSize',30)
    
    count = count+1;
end
Xlim=([K67.year(1)-1 K67.year(end)+1]);
set(gca,'XTick',K67.year(1)-1:K67.year(end)+1,'XTickLabel',datenum(K67.year,1,1),...
    'Ylim',[v1 v2],'YTick',v3,'FontSize',11);
ylabel ({'Mortality (gC m^-^2 y^-^1)'},'FontSize',12);     xlabel('Sites in this proposal','FontSize',12);
grid on;                box on;     datetick('x','YY')

%% ........................................................................
v1 = 0;                         v2 = 4;       v3 = v1:((v2-v1)/4):v2;
mortalityYdate = datenum(mortalityY,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
for ix = 1:length(mortalityYdate)
    low = mortalityMin(ix);
    high = (mortalityMax(ix)-mortalityMin(ix));
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',[.7 .7 .9],'EdgeColor',[.0 .0 .2]);
    set(h(1),'FaceColor','none','EdgeColor','none');
    set(h,'BaseValue',0);   hold on;
    
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[mortalityAvg(ix);mortalityAvg(ix)],'-k','MarkerSize',30)
    
    count = count+1;
end
%plot(K67.year_dateloc,(K67.year_dateloc.*0)+nanmean(K67.year_avgMortality),':','LineWidth',1,'color',[0.4 0.4 0.6]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
%plot(K67.year_dateloc+140,K67.year_avgMortality,':','LineWidth',1,'color',[0.4 0.4 0.6]);
grid on;                box on;
ylabel ({'Mortality (%AGB year^-^1)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;   xtickangle(0);


% ........................................................................
% ........................................................................
%% Demography
% ........................................................................
% ........................................................................
%
opts = delimitedTextImportOptions("NumVariables", 13);
opts.DataLines = [2, Inf];              opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["Var1", "yr1999", "yr2001", "yr2005", "yr2008", "yr2009", "yr2010", "yr2011", "yr2012", "yr2013", "yr2014", "yr2015", "yr2017"];
opts.SelectedVariableNames = ["yr1999", "yr2001", "yr2005", "yr2008", "yr2009", "yr2010", "yr2011", "yr2012", "yr2013", "yr2014", "yr2015", "yr2017"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";       opts.EmptyLineRule = "read";
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
% rowname transect 1 to 4
AGBtotal = readtable([FolderK67dendrometry 'K67.AGB.total.2021.csv'], opts);
MortalityTotal = readtable([FolderK67dendrometry 'K67.mortality.total.2021.csv'], opts);
RecruitmentTotal = readtable([FolderK67dendrometry 'K67.recruitment.total.2021.csv'], opts);
GrowthTotal = readtable([FolderK67dendrometry 'K67.grow.total.2021.csv'], opts);

% rowname expected q025 q075
opts = delimitedTextImportOptions("NumVariables", 13);
opts.DataLines = [3, Inf];              opts.Delimiter = ",";
opts.VariableNames = ["Var1", "yr1999a", "yr1999b", "yr1999c", "yr1999d", "yr2001a", "yr2001b", "yr2001c", "yr2001d", "yr2005a", "yr2005b", "yr2005c", "yr2005d", "yr2008a", "yr2008b", "yr2008c", "yr2008d", "yr2009a", "yr2009b", "yr2009c", "yr2009d", "yr2010a", "yr2010b", "yr2010c", "yr2010d", "yr2011a", "yr2011b", "yr2011c", "yr2011d", "yr2012a", "yr2012b", "yr2012c", "yr2012d", "yr2013a", "yr2013b", "yr2013c", "yr2013d", "yr2014a", "yr2014b", "yr2014c", "yr2014d", "yr2015a", "yr2015b", "yr2015c", "yr2015d", "yr2017a", "yr2017b", "yr2017c", "yr2017d"];
opts.SelectedVariableNames = ["yr1999a", "yr1999b", "yr1999c", "yr1999d", "yr2001a", "yr2001b", "yr2001c", "yr2001d", "yr2005a", "yr2005b", "yr2005c", "yr2005d", "yr2008a", "yr2008b", "yr2008c", "yr2008d", "yr2009a", "yr2009b", "yr2009c", "yr2009d", "yr2010a", "yr2010b", "yr2010c", "yr2010d", "yr2011a", "yr2011b", "yr2011c", "yr2011d", "yr2012a", "yr2012b", "yr2012c", "yr2012d", "yr2013a", "yr2013b", "yr2013c", "yr2013d", "yr2014a", "yr2014b", "yr2014c", "yr2014d", "yr2015a", "yr2015b", "yr2015c", "yr2015d", "yr2017a", "yr2017b", "yr2017c", "yr2017d"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
AGBclass    = readtable([FolderK67dendrometry 'K67.AGB.size.class.2021.csv'], opts);
MortalityClass     = readtable([FolderK67dendrometry 'K67.mortality.size.class.2021.csv'], opts);
RecruitmentClass   = readtable([FolderK67dendrometry 'K67.recruitment.size.class.2021.csv'], opts);
GrowthClass = readtable([FolderK67dendrometry 'K67.grow.size.class.2021.csv'], opts);
%              a        b         c        d   classes
ClassSize = {'<20cm';'20-35cm';'35-55cm';'>55cm'};
ColorFaceSize = [ .2 .2 .9;.2 .9 .2;.9 .9 .2;.9 .2 .2]; ColorEdgeSize = [ 0 0 .2;.0 .2 0; .2 .2 0; .2 0 0];
Transect = {'T01';'T02';'T03';'T04'};
TransectTick = {'o';'d';'s';'^'};
Probability = {'Expectation';'Q25';'Q75'};
clear opts

% Figures
AGB.year = [1999;2001;2005;2008;2009;2010;2011;2012;2013;2014;2015;2017];
AGB.total = [];
AGB.total(:,1) = (AGBtotal{1,:})';        AGB.total(:,2) = (AGBtotal{2,:})';
AGB.total(:,3) = (AGBtotal{3,:})';        AGB.total(:,4) = (AGBtotal{4,:})';
AGB.total(:,5) = (nanmean(AGB.total(:,1:4),2));
AGB.total(:,6) = (nanstd(AGB.total(:,1:4),[],2));
for ix = 1:length(AGB.year)
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'a'];
    AGB.total(ix,7:10) = (AGBclass{:,DemographyYearName})';
    AGB.total(:,11) = (nanmean(AGB.total(6:9),2));
    AGB.total(:,12) = (nanstd(AGB.total(6:9),[],2));
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'b'];
    AGB.total(ix,13:16) = (AGBclass{:,DemographyYearName})';
    AGB.total(:,17) = (nanmean(AGB.total(11:14),2));
    AGB.total(:,18) = (nanstd(AGB.total(11:14),[],2));
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'c'];
    AGB.total(ix,19:22) = (AGBclass{:,DemographyYearName})';
    AGB.total(:,23) = (nanmean(AGB.total(16:19),2));
    AGB.total(:,24) = (nanstd(AGB.total(16:19),[],2));
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'd'];
    AGB.total(ix,25:28) = (AGBclass{:,DemographyYearName})';
    AGB.total(:,29) = (nanmean(AGB.total(21:24),2));
    AGB.total(:,30) = (nanstd(AGB.total(21:24),[],2));
end

%
Mortality.total = [];
Mortality.total(:,1) = (MortalityTotal{1,:})';        Mortality.total(:,2) = (MortalityTotal{2,:})';
Mortality.total(:,3) = (MortalityTotal{3,:})';
for ix = 1:length(AGB.year)
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'a'];
    Mortality.total(ix,4:6) = (MortalityClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'b'];
    Mortality.total(ix,7:9) = (MortalityClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'c'];
    Mortality.total(ix,10:12) = (MortalityClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'd'];
    Mortality.total(ix,13:15) = (MortalityClass{:,DemographyYearName})';
end

Growth.total = [];
Growth.total(:,1) = (GrowthTotal{1,:})';        Growth.total(:,2) = (GrowthTotal{2,:})';
Growth.total(:,3) = (GrowthTotal{3,:})';
for ix = 1:length(AGB.year)
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'a'];
    Growth.total(ix,4:6) = (GrowthClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'b'];
    Growth.total(ix,7:9) = (GrowthClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'c'];
    Growth.total(ix,10:12) = (GrowthClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'd'];
    Growth.total(ix,13:15) = (GrowthClass{:,DemographyYearName})';
end

Recruitment.total = [];
Recruitment.total(:,1) = (RecruitmentTotal{1,:})';        Recruitment.total(:,2) = (RecruitmentTotal{2,:})';
Recruitment.total(:,3) = (RecruitmentTotal{3,:})';
for ix = 1:length(AGB.year)
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'a'];
    Recruitment.total(ix,4:6) = (RecruitmentClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'b'];
    Recruitment.total(ix,7:9) = (RecruitmentClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'c'];
    Recruitment.total(ix,10:12) = (RecruitmentClass{:,DemographyYearName})';
    DemographyYearName = ['yr' num2str(AGB.year(ix)) 'd'];
    Recruitment.total(ix,13:15) = (RecruitmentClass{:,DemographyYearName})';
end

%% .........................................................................
K67.year_avgAGB = interp1(AGB.year,AGB.total(:,5),K67.year);
K67.year_avgAGBT1 = interp1(AGB.year,AGB.total(:,1),K67.year);
K67.year_avgAGBT2 = interp1(AGB.year,AGB.total(:,2),K67.year);
K67.year_avgAGBT3 = interp1(AGB.year,AGB.total(:,3),K67.year);
K67.year_avgAGBT4 = interp1(AGB.year,AGB.total(:,4),K67.year);

v1 = 15;                         v2 = 20;       v3 = v1:((v2-v1)/5):v2;
mortalityYdate = datenum(AGB.year,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
for ix = 1:length(mortalityYdate)
    low = AGB.total(ix,5)-AGB.total(ix,6);    high = 2.*AGB.total(ix,6);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.0 .0 .0]);
    set(h(1),'FaceColor','none','EdgeColor','none');
    set(h,'BaseValue',0);   hold on;
    
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[AGB.total(ix,5);AGB.total(ix,5)],'-k','MarkerSize',30)
    hl2 = plot(mortalityYdate(ix)+160,AGB.total(ix,1),'.','MarkerSize',15,'color',[.9 0 0])
    hl3 = plot(mortalityYdate(ix)+160,AGB.total(ix,2),'.','MarkerSize',15,'color',[.7 0 0])
    hl4 = plot(mortalityYdate(ix)+160,AGB.total(ix,3),'.','MarkerSize',15,'color',[.5 0 0])
    hl5 = plot(mortalityYdate(ix)+160,AGB.total(ix,4),'.','MarkerSize',15,'color',[.5 .0 .5])
    
    count = count+1;
end
plot(K67.year_dateloc,(K67.year_dateloc.*0)+nanmean(AGB.total(:,5)),':','LineWidth',1,'color',[0.4 0.4 0.6]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);    XlabelYear(XlabelYear>datenum(K67.year(end)+1,1,1)) = [];
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
hl1 = plot(K67.year_dateloc+140,K67.year_avgAGB,'-','LineWidth',1,'color',[0.0 0.0 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT1,'--','LineWidth',1,'color',[0.9 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT2,'--','LineWidth',1,'color',[0.7 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT3,'--','LineWidth',1,'color',[0.5 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT4,'--','LineWidth',1,'color',[0.5 0. 0.5]);
grid on;                box on;
ylabel ({'AGB (kgC m^-^2)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;
legend([hl1,hl2,hl3,hl4,hl5],'total','T1','T2','T3','T4','orientation','horizontal');
legend box off;

%% ........................................................................
K67.year_avgMortality = interp1(AGB.year,Mortality.total(:,1),K67.year);
K67.year_avgMortalitySM = interp1(AGB.year,Mortality.total(:,4),K67.year);
K67.year_avgMortalityMD = interp1(AGB.year,Mortality.total(:,7),K67.year);
K67.year_avgMortalityLG = interp1(AGB.year,Mortality.total(:,10),K67.year);
K67.year_avgMortalityXL = interp1(AGB.year,Mortality.total(:,13),K67.year);

v1 = 0;                         v2 = 1.2;       v3 = v1:((v2-v1)/4):v2;
mortalityYdate = datenum(AGB.year,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
for ix = 1:length(mortalityYdate)
    %     low = Mortality.total(ix,1)-Mortality.total(ix,2);    high = 2.*Mortality.total(ix,2);
    low = Mortality.total(ix,2);        high = Mortality.total(ix,3)-Mortality.total(ix,2);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+100,mortalityYdate(ix)+260];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.2 .2 .2]);
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);
    hold on;
    
    im = 1;                 % small
    %     low = Mortality.total(ix,1)-.1;    high = 2.*.1;
    low = Mortality.total(ix,5);    high = Mortality.total(ix,6)-Mortality.total(ix,5);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,4);Mortality.total(ix,4)],'-k','MarkerSize',30)
    
    im = 2;                 % small
    low = Mortality.total(ix,8);    high = Mortality.total(ix,9)-Mortality.total(ix,8);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,7);Mortality.total(ix,7)],'-k','MarkerSize',30)
    
    im = 3;                 % small
    low = Mortality.total(ix,11);    high = Mortality.total(ix,12)-Mortality.total(ix,11);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,10);Mortality.total(ix,10)],'-k','MarkerSize',30)
    
    im = 4;                 % small
    low = Mortality.total(ix,14);    high = Mortality.total(ix,15)-Mortality.total(ix,14);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,13);Mortality.total(ix,13)],'-k','MarkerSize',30)
    
    plot([mortalityYdate(ix)+100;mortalityYdate(ix)+260],[Mortality.total(ix,1);Mortality.total(ix,1)],'-k','MarkerSize',30)
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,4);Mortality.total(ix,4)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(1,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,7);Mortality.total(ix,7)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(2,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,10);Mortality.total(ix,10)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(3,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Mortality.total(ix,13);Mortality.total(ix,13)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(4,:))
    
    count = count+1;
end
plot(K67.year_dateloc+150,(K67.year_dateloc.*0)+nanmean(Mortality.total(:,1)),':','LineWidth',1,'color',[0.4 0.4 0.4]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);    XlabelYear(XlabelYear>datenum(K67.year(end)+1,1,1)) = [];
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
hl1 = plot(K67.year_dateloc+140,K67.year_avgMortality,'-','LineWidth',2,'color',[0.4 0.4 0.4]);
hl2 = plot(K67.year_dateloc+140,K67.year_avgMortalitySM,'-','LineWidth',2,'color',ColorFaceSize(1,:));
hl3 = plot(K67.year_dateloc+140,K67.year_avgMortalityMD,'-','LineWidth',2,'color',ColorFaceSize(2,:));
hl4 = plot(K67.year_dateloc+140,K67.year_avgMortalityLG,'-','LineWidth',2,'color',ColorFaceSize(3,:));
hl5 = plot(K67.year_dateloc+140,K67.year_avgMortalityXL,'-','LineWidth',2,'color',ColorFaceSize(4,:));
grid on;                box on;     xtickangle(0);
ylabel ({'Mortality (kgC m^-^2 yr^-^1)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;
legend([hl1,hl2,hl3,hl4,hl5],'total','<20','20-35','35-55','>55 cm','orientation','horizontal');
legend box off;     xtickangle(0);

%% ........................................................................
K67.year_avgRecruitment   = interp1(AGB.year,Recruitment.total(:,1),K67.year);
K67.year_avgRecruitmentSM = interp1(AGB.year,Recruitment.total(:,4),K67.year);
K67.year_avgRecruitmentMD = interp1(AGB.year,Recruitment.total(:,7),K67.year);
K67.year_avgRecruitmentLG = interp1(AGB.year,Recruitment.total(:,10),K67.year);
K67.year_avgRecruitmentXL = interp1(AGB.year,Recruitment.total(:,13),K67.year);

v1 = 0;                         v2 = 0.6;       v3 = v1:((v2-v1)/3):v2;
mortalityYdate = datenum(AGB.year,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
for ix = 1:length(mortalityYdate)
    % %     low = Recruitment.total(ix,2);    high = Recruitment.total(ix,3)-Recruitment.total(ix,2);
    % %     y = [low,high;low,high];
    % %     x = [mortalityYdate(ix)+100,mortalityYdate(ix)+260];
    % %     h = area(x,y);      set(gca,'Layer','top');    hold on;
    % %     set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.2 .2 .2]);
    % %     set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);
    hold on;
    
    im = 1;                 % small
    %     low = Recruitment.total(ix,1)-.1;    high = 2.*.1;
    low = Recruitment.total(ix,5);    high = Recruitment.total(ix,6)-Recruitment.total(ix,5);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,4);Recruitment.total(ix,4)],'-k','MarkerSize',30)
    
    im = 2;                 % small
    low = Recruitment.total(ix,8);    high = Recruitment.total(ix,9)-Recruitment.total(ix,8);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,7);Recruitment.total(ix,7)],'-k','MarkerSize',30)
    
    im = 3;                 % small
    low = Recruitment.total(ix,11);    high = Recruitment.total(ix,12)-Recruitment.total(ix,11);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,10);Recruitment.total(ix,10)],'-k','MarkerSize',30)
    
    im = 4;                 % small
    low = Recruitment.total(ix,14);    high = Recruitment.total(ix,15)-Recruitment.total(ix,14);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,13);Recruitment.total(ix,13)],'-k','MarkerSize',30)
    
    %     plot([mortalityYdate(ix)+100;mortalityYdate(ix)+260],[Recruitment.total(ix,1);Recruitment.total(ix,1)],'-k','MarkerSize',30)
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,4);Recruitment.total(ix,4)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(1,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,7);Recruitment.total(ix,7)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(2,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,10);Recruitment.total(ix,10)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(3,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Recruitment.total(ix,13);Recruitment.total(ix,13)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(4,:))
    
    count = count+1;
end
% plot(K67.year_dateloc+150,(K67.year_dateloc.*0)+nanmean(Recruitment.total(:,1)),':','LineWidth',1,'color',[0.4 0.4 0.4]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);    XlabelYear(XlabelYear>datenum(K67.year(end)+1,1,1)) = [];
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
hl2 = plot(K67.year_dateloc+140,K67.year_avgRecruitmentSM,'-','LineWidth',2,'color',ColorFaceSize(1,:));
hl3 = plot(K67.year_dateloc+140,K67.year_avgRecruitmentMD,'-','LineWidth',2,'color',ColorFaceSize(2,:));
hl4 = plot(K67.year_dateloc+140,K67.year_avgRecruitmentLG,'-','LineWidth',2,'color',ColorFaceSize(3,:));
hl5 = plot(K67.year_dateloc+140,K67.year_avgRecruitmentXL,'-','LineWidth',2,'color',ColorFaceSize(4,:));
hl1 = plot(K67.year_dateloc+140,K67.year_avgRecruitment,'-','LineWidth',2,'color',[0.4 0.4 0.4]);
grid on;                box on;
ylabel ({'\Delta_{size class} (kgC m^-^2 yr^-^1)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;
% legend([hl1,hl2,hl3,hl4,hl5],'total','<20cm','20-35cm','35-55cm','>55cm','orientation','horizontal');
% legend box off;

%% ........................................................................
K67.year_avgGrowth = interp1(AGB.year,Growth.total(:,1),K67.year);
K67.year_avgGrowthSM = interp1(AGB.year,Growth.total(:,4),K67.year);
K67.year_avgGrowthMD = interp1(AGB.year,Growth.total(:,7),K67.year);
K67.year_avgGrowthLG = interp1(AGB.year,Growth.total(:,10),K67.year);
K67.year_avgGrowthXL = interp1(AGB.year,Growth.total(:,13),K67.year);

v1 = 0;                         v2 = 0.6;       v3 = v1:((v2-v1)/4):v2;
mortalityYdate = datenum(AGB.year,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
for ix = 1:length(mortalityYdate)
    low = Growth.total(ix,2);       high = Growth.total(ix,3)-Growth.total(ix,2);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+100,mortalityYdate(ix)+260];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.2 .2 .2]);
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);
    hold on;
    
    im = 1;                 % small
    low = Growth.total(ix,5);       high = Growth.total(ix,6)-Growth.total(ix,5);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
    h = area(x,y);      set(gca,'Layer','top');         hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,4);Growth.total(ix,4)],'-k','MarkerSize',30)
    
    im = 2;                 % small
    low = Growth.total(ix,8);       high = Growth.total(ix,9)-Growth.total(ix,8);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,7);Growth.total(ix,7)],'-k','MarkerSize',30)
    
    im = 3;                 % small
    low = Growth.total(ix,11);    high = Growth.total(ix,12)-Growth.total(ix,11);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,10);Growth.total(ix,10)],'-k','MarkerSize',30)
    
    im = 4;                 % small
    low = Growth.total(ix,14);    high = Growth.total(ix,15)-Growth.total(ix,14);
    y = [low,high;low,high];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',ColorFaceSize(im,:),'EdgeColor',ColorEdgeSize(im,:));
    set(h(1),'FaceColor','none','EdgeColor','none');    set(h,'BaseValue',0);         alpha(.5);
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,13);Growth.total(ix,13)],'-k','MarkerSize',30)
    
    
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,1);Growth.total(ix,1)],'-k','MarkerSize',30)
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,4);Growth.total(ix,4)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(1,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,7);Growth.total(ix,7)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(2,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,10);Growth.total(ix,10)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(3,:))
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[Growth.total(ix,13);Growth.total(ix,13)],'-',...
        'MarkerSize',30,'color',ColorEdgeSize(4,:))
    
    count = count+1;
end
plot(K67.year_dateloc+150,(K67.year_dateloc.*0)+nanmean(Growth.total(:,1)),':','LineWidth',1,'color',[0.4 0.4 0.4]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);    XlabelYear(XlabelYear>datenum(K67.year(end)+1,1,1)) = [];
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
hl1 = plot(K67.year_dateloc+140,K67.year_avgGrowth,'-','LineWidth',2,'color',[0.4 0.4 0.4]);
hl2 = plot(K67.year_dateloc+140,K67.year_avgGrowthSM,'-','LineWidth',2,'color',ColorFaceSize(1,:));
hl3 = plot(K67.year_dateloc+140,K67.year_avgGrowthMD,'-','LineWidth',2,'color',ColorFaceSize(2,:));
hl4 = plot(K67.year_dateloc+140,K67.year_avgGrowthLG,'-','LineWidth',2,'color',ColorFaceSize(3,:));
hl5 = plot(K67.year_dateloc+140,K67.year_avgGrowthXL,'-','LineWidth',2,'color',ColorFaceSize(4,:));
grid on;                box on;
ylabel ({'Growth (kgC m^-^2 yr^-^1)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;
legend([hl1,hl2,hl3,hl4,hl5],'total','<20cm','20-35cm','35-55cm','>55cm','orientation','horizontal');
legend box off;


%% .........................................................................
% % K67.year_avgAGB = interp1(AGB.year,AGB.total(:,5),K67.year);
% % K67.year_avgAGBT1 = interp1(AGB.year,AGB.total(:,1),K67.year);
% % K67.year_avgAGBT2 = interp1(AGB.year,AGB.total(:,2),K67.year);
% % K67.year_avgAGBT3 = interp1(AGB.year,AGB.total(:,3),K67.year);
% % K67.year_avgAGBT4 = interp1(AGB.year,AGB.total(:,4),K67.year);

v1 = 15;                         v2 = 20;       v3 = v1:((v2-v1)/5):v2;
mortalityYdate = datenum(AGB.year,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;        xtickangle(0);
% for ix = 1:length(mortalityYdate)
% %     low = AGB.total(ix,5)-AGB.total(ix,6);    high = 2.*AGB.total(ix,6);
% %     y = [low,high;low,high];
% %     x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
% %     h = area(x,y);      set(gca,'Layer','top');    hold on;
% %     set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.0 .0 .0]);
% %     set(h(1),'FaceColor','none','EdgeColor','none');
% %     set(h,'BaseValue',0);   hold on;

%     plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[AGB.total(ix,5);AGB.total(ix,5)],'-k','MarkerSize',30)
hl2 = plot(mortalityYdate+160,AGB.total(:,1),'Marker',TransectTick{1},'MarkerSize',5,'color',[.9 .0 .0],'MarkerFaceColor',[.9 .0 .0]);
hl3 = plot(mortalityYdate+160,AGB.total(:,2),'Marker',TransectTick{2},'MarkerSize',5,'color',[.7 .0 .0],'MarkerFaceColor',[.7 .0 .0]);
hl4 = plot(mortalityYdate+160,AGB.total(:,3),'Marker',TransectTick{3},'MarkerSize',5,'color',[.5 .0 .0],'MarkerFaceColor',[.5 .0 .0]);
hl5 = plot(mortalityYdate+160,AGB.total(:,4),'Marker',TransectTick{4},'MarkerSize',5,'color',[.5 .0 .5],'MarkerFaceColor',[.5 .0 .5]);

% % end
% for ix = 1:4
%     plot(mortalityYdate+200,AGB.total(:,ix*3-1),'.','MarkerSize',30,'color',ColorEdgeSize(ix,:))
% end

% plot(K67.year_dateloc,(K67.year_dateloc.*0)+nanmean(AGB.total(:,5)),':','LineWidth',1,'color',[0.4 0.4 0.6]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);    XlabelYear(XlabelYear>datenum(K67.year(end)+1,1,1)) = [];
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
hl1 = plot(K67.year_dateloc+140,K67.year_avgAGB,'-','LineWidth',1,'color',[0.0 0.0 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT1,'--','LineWidth',1,'color',[0.9 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT2,'--','LineWidth',1,'color',[0.7 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT3,'--','LineWidth',1,'color',[0.5 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT4,'--','LineWidth',1,'color',[0.5 0. 0.5]);
grid on;                box on;
ylabel ({'AGB (kgC m^-^2)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;        xtickangle(0);
legend([hl1,hl2,hl3,hl4,hl5],'total','T1','T2','T3','T4','orientation','horizontal');
legend box off;

%% .........................................................................
K67.year_avgAGB = interp1(AGB.year,AGB.total(:,5),K67.year);
K67.year_avgAGBT1 = interp1(AGB.year,AGB.total(:,1),K67.year);
K67.year_avgAGBT2 = interp1(AGB.year,AGB.total(:,2),K67.year);
K67.year_avgAGBT3 = interp1(AGB.year,AGB.total(:,3),K67.year);
K67.year_avgAGBT4 = interp1(AGB.year,AGB.total(:,4),K67.year);

v1 = 15;                         v2 = 20;       v3 = v1:((v2-v1)/5):v2;
mortalityYdate = datenum(AGB.year,1,1);
figure('color','white','Position',[50 50 600 400],'units','inches');
ix = find(K67.monthTRMM<100);    iy = ones(length(K67.month),1);	iy = iy.*v1;
iz = iy;                         iz(ix) = (v2+v1);
y  = [iy,iz];
h  = bar(K67.month',y,4.0);      set(gca,'Layer','top');
set(h(2),'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);       set(h(1),'FaceColor','none','EdgeColor','none');
set(h,'BaseValue',v1);          hold on;
for ix = 1:length(mortalityYdate)
    low = AGB.total(ix,5)-AGB.total(ix,6);    high = 2.*AGB.total(ix,6);
    y = [low,high;low,high];
    x = [mortalityYdate(ix)+140,mortalityYdate(ix)+220];
    h = area(x,y);      set(gca,'Layer','top');    hold on;
    set(h(2),'FaceColor',[.7 .7 .7],'EdgeColor',[.0 .0 .0]);
    set(h(1),'FaceColor','none','EdgeColor','none');
    set(h,'BaseValue',0);   hold on;
    
    plot([mortalityYdate(ix)+140;mortalityYdate(ix)+220],[AGB.total(ix,5);AGB.total(ix,5)],'-k','MarkerSize',30)
    % % %     hl2 = plot(mortalityYdate(ix)+160,AGB.total(ix,1),'.','MarkerSize',15,'color',[.9 0 0])
    % % %     hl3 = plot(mortalityYdate(ix)+160,AGB.total(ix,2),'.','MarkerSize',15,'color',[.7 0 0])
    % % %     hl4 = plot(mortalityYdate(ix)+160,AGB.total(ix,3),'.','MarkerSize',15,'color',[.5 0 0])
    % % %     hl5 = plot(mortalityYdate(ix)+160,AGB.total(ix,4),'.','MarkerSize',15,'color',[.5 .0 .5])
    
    count = count+1;
end
plot(K67.year_dateloc,(K67.year_dateloc.*0)+nanmean(AGB.total(:,5)),':','LineWidth',1,'color',[0.4 0.4 0.6]);
hold off;
hold on;
XlabelYear = datenum(K67.year(1:3:end),1,1);    XlabelYear(XlabelYear>datenum(K67.year(end)+1,1,1)) = [];
set(gca,'XTick',XlabelYear,'XTickLabel',datestr(XlabelYear,'yyyy'),'YTick',v3,...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)],'YLim',[v1 v2],'FontSize',14);
hl1 = plot(K67.year_dateloc+140,K67.year_avgAGB,'-','LineWidth',1,'color',[0.0 0.0 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT1,'--','LineWidth',1,'color',[0.9 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT2,'--','LineWidth',1,'color',[0.7 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT3,'--','LineWidth',1,'color',[0.5 0. 0.]);
plot(K67.year_dateloc+140,K67.year_avgAGBT4,'--','LineWidth',1,'color',[0.5 0. 0.5]);
grid on;                box on;
ylabel ({'AGB (kgC m^-^2)'},'FontSize',12);
ax1=gca;
ax2 = axes('YAxisLocation','right','Position',get(ax1,'Position'),'XAxisLocation','bottom',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick',datenum(K67.year,1,1),'XTickLabel',['','',''],...
    'XLim',[datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
set(ax2,'Layer','top','Ylim',[v1 v2],'YTick',['',''],'FontSize',14,'Position',get(gca,'Position'));
xlim([datenum(K67.year(1),1,1) datenum(K67.year(end)+1,1,1)]);
grid on;                box off;        xtickangle(0);
% legend([hl1,hl2,hl3,hl4,hl5],'total','T1','T2','T3','T4','orientation','horizontal');
% legend box off;
