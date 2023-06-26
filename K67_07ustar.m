% %Follows Richardson, A.D. et al. 2007
% %Calculates Re minimal use of environmental covariates (good thing)
% %*Annual* second order Fourier regression
% % .......................................................................

% .........................................................................
% .........................................................................
%Night-time values
% .........................................................................
% .........................................................................
flag_unix=0;
flag_ust_uncertain=1;       %run K67.ustar uncertainity =1 no=0
flag_calc=1;                flag_figure_ps=1;

if flag_unix==1
    load 'C:\Users\ncoupe\Documents\Amazon\AllSites\matlabBFstruc.mat'
    FileName =  'C:\Users\ncoupe\Documents\Amazon\AllSites\';
elseif flag_unix==0
    load '/home/ncoupe/Documents/Amazon/AllSites/matlabBFstruc.mat'
    FileName = '/home/ncoupe/Documents/Amazon/AllSites/';
end
addpath /home/ncoupe/Documents/OZ/;
addpath /home/ncoupe/Documents/Amazon/;

[site_name_BF, lat_BF, lon_BF, site_long_BF, year_ini, year_end, siteelevation_BF,...
    height_BF, canopy_BF, ust_tresh_BF, sunrise_BF, sunset_BF, utm_BF, ust_min_BF, ust_max_BF,...
    LAI_BF]  =  textread([FolderFlux 'sites_loc.txt'],...
    '%s	%f	%f	%s	%f	%f	%f	%f %f	%f	%f	%f	%f	%f	%f	%f',...
    'delimiter',',','headerlines', 1);
% load('matlabUSTuncertain.mat');

if flag_calc==1
    radiation(K67.So==0)=0;
    a = length(K67.dateloc);        g = length(K67.year);       e = length(K67.week);
    sunrise = 6;                    sunset = 18;
    % ..............................................................
    % K67.nee spike removal and night-time calculations
    % ..............................................................
    %     [NEE_flag4,NEE_spike4_free]=AM_spike_detection_filled(K67.nee,K67.dateloc,4,sunrise,sunset);
    [NEE_flag6,NEE_spike6_free]=AM_spike_detection_filled(K67.nee,K67.dateloc,6,sunrise,sunset);
    NEE_spike6_free(NEE_spike6_free<3) = NaN;
    %     [NEE_flag7,NEE_spike7_free]=AM_spike_detection_filled(K67.nee,K67.dateloc,7,sunrise,sunset);
    %     [NEE_flag9,NEE_spike9_free]=AM_spike_detection_filled(K67.nee,K67.dateloc,9,sunrise,sunset);
    %     length(K67.nee(~isnan(K67.nee)))
    %     length(NEE_spike4_free(~isnan(NEE_spike4_free)))
    %     length(NEE_spike6_free(~isnan(NEE_spike6_free)))
    %     length(NEE_spike7_free(~isnan(NEE_spike7_free)))
    %     length(NEE_spike9_free(~isnan(NEE_spike9_free)))
    
    NEE_night = NEE_spike6_free;	ta_night = K67.tasp;        ust_night = K67.ustar;
    ind = find((K67.So>0)|(K67.radiation>0));
    NEE_night(ind)=NaN;             ta_night(ind)=NaN;          ust_night(ind)=NaN;
    
    NEE_night=AM_rm_outlier(NEE_night,3);
    
    NEE_night_ust=NEE_night;         NEE_night_ust((isnan(K67.ustar)))=NaN;
    NEE_night_ust=AM_rm_outlier(NEE_night_ust,3);
    
    if flag_ust_uncertain==1
        % ..............................................................
        %% u* correction
        % ..............................................................
        %         figure ('color','white');
        %         subplot(2,1,1);
        %         plot(K67.nee);      hold on;        %[K67.nee,ix]=AM_rm_outlier(K67.nee,3);
        %         K67.nee((K67.nee<-40)|(K67.nee>40))=NaN;    plot(K67.nee,':c');     legend('NEE','NEE_{rm 3std}');
        %
        %         subplot(2,1,2); hold on;
        %         plot(K67.nee);
        %         plot(NEE_spike6_free,'c');           plot(NEE_spike7_free,':m');    plot(NEE_spike4_free,':g');
        %         legend('NEE','NEE_{spk6}','NEE_{spk7}','NEE_{spk4}');
        
        % ..............................................................
        %% night-time values
        % ..............................................................
        ind_night=find((K67.radiation>0)|(K67.So>0));
        NEE_night_all = NEE_spike6_free;              NEE_night_all(NEE_night_all<3) = NaN;
        NEE_night_all(ind_night)=NaN;
        ta_night_all  = AM_fill_plus(K67.tasp);       ta_night_all(ind_night)=NaN;
        ust_night_all = K67.ustar;                      ust_night_all(ind_night)=NaN;
        %             NEE_night_all=AM_rm_outlier(NEE_night,2);   ta_night_all=AM_rm_outlier(ta_night,2);
        %             ust_night_all=AM_rm_outlier(ust_night,2);
        % ..............................................................
        %% u* correction
        % ..............................................................
        for im = 1:4
            NEE_night = NEE_spike6_free;              NEE_night(ind_night) = NaN;
            ta_night = AM_fill_plus(K67.tasp);              ta_night(ind_night) = NaN;
            ust_night = K67.ustar;                          ust_night(ind_night) = NaN;
            id_month = (im*3)-2:(im*3);
            
            ta_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];
            NEE_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];
            ust_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];
            
            % split 7 Tair classes of ~equal sample size
            [ta_class,NEE_class_avg,NEE_class_std,NEE_class_count]=AM_class_eq_size(ta_night,NEE_night,7);
            if isnan(ta_class)
                NEE_night = NEE_spike6_free;                    NEE_night(ind_night) = NaN;
                ta_night  = AM_fill_plus(K67.tasp);             ta_night(ind_night)  = NaN;
                ust_night = K67.ustar;                          ust_night(ind_night) = NaN;
                if im>1;  id_month=((im-1)*3)-2:((im-1)*3);     end
                if im==1; id_month=((im+1)*3)-2:((im+1)*3);     end
                
                ta_night((K67.M<id_month(1))|(K67.M>id_month(end)))  = [];
                NEE_night((K67.M<id_month(1))|(K67.M>id_month(end))) = [];
                ust_night((K67.M<id_month(1))|(K67.M>id_month(end))) = [];
                [ta_class,NEE_class_avg,NEE_class_std,NEE_class_count] = AM_class_eq_size(ta_night,NEE_night,7);
            end
            
            ta_class = [ta_class(1,1);ta_class(:,2)];
            min_ust  = (floor((min(ust_night*100))))/100;  max_ust=(ceil(max(ust_night*100)))/100;
            ust_class = (min_ust:((max_ust-min_ust)/20):max_ust)';
            x3        = (min_ust:((max_ust-min_ust)/5):max_ust)';
            %             [ust_class]=AM_class_eq_size(ust_night,NEE_night,10);
            %             ust_class=[ust_class(1,1);ust_class(:,2)];
            min_sample_size = 3;
            [LUT_ta_ust,LUT_ta_ust_std,LUT_ta_ust_count] = AM_LUTgeneral(ta_night,ust_night,NEE_night,...
                ta_class,ust_class,1,'Ta(\circC)','u_{*}(m s^-^1)',...
                'NEE (\mumol im^-2 s^-^1)',min_sample_size);
            ind = find((LUT_ta_ust_count<40)&(LUT_ta_ust_std>nanmean(nanmean(LUT_ta_ust_std))*2));
            LUT_ta_ust(ind) = NaN;      LUT_ta_ust_std(ind) = NaN;      LUT_ta_ust_count(ind) = NaN;
            % ..............................................................
            % %             %        figure ('color','white');
            % %             for ik=1:(length(ust_class))-1
            % %                 %             plot(LUT_ta_ust(:,ik),'o','color',[0 1./ik 0.01*ik]); hold on;
            % %                 ix = AM_rm_outlier(LUT_ta_ust(:,ik),2);   LUT_ta_ust(:,ik) = ix;
            % %                 %             plot(LUT_ta_ust(:,ik),'.','color',[0 1./ik 0.01*ik]);
            % %             end
            % %             %         close all;
            %%
            color_mtx = jet(length(ta_class)-1);
            figure('color','white');
            subplot(2,2,1);
            for ik=1:length(ta_class)-1
                plot(ust_class(2:end),((LUT_ta_ust(ik,:))'),'color',color_mtx(ik,:)); hold on;
                plot(ust_class(2:end),((LUT_ta_ust(ik,:))'),'.','MarkerSize',12,'color',color_mtx(ik,:));
                LUT_ta_ust(ik,:) = AM_rm_outlier(LUT_ta_ust(ik,:),3);
                ix = (LUT_ta_ust(ik,:));                iy=1:length(ust_class)-1;
                ind = find((isnan(ix))|(isnan(iy)));    ix(ind)=[];         iz=iy;    iz(ind)=[];
                if length(ix)>2
                    iq=interp1(iz,ix,iy);                LUT_ta_ust(ik,:)=iq;
                end
%                 plot(ust_class(2:end),((LUT_ta_ust(ik,:))'),'o','MarkerSize',6,'color',color_mtx(ik,:));
            end
            title([month_name(im*3-2,:) '-' month_name(im*3,:)],'FontWeight','normal');
            ylabel('NEE_{night}(\mumol m^-^2 s^-^1)');   xlabel('u_{*night}(m s^-^1)');
            set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',x3,...
                'XAxisLocation','bottom','XTick',x3,'Xlim',([ust_class(1) ust_class(end)]),...
                'FontSize',12);          axis tight;
            subplot(2,2,2);
            for ik=1:length(ta_class)-1
                plot(ust_class(2:end),((LUT_ta_ust_std(ik,:))'),'color',color_mtx(ik,:)); hold on;
                plot(ust_class(2:end),((LUT_ta_ust_std(ik,:))'),'.','MarkerSize',12,'color',color_mtx(ik,:));
            end
            title(['Ta - u_{*} - NEE'],'FontWeight','normal');
            ylabel('NEE_{night_{std}}(\mumol m^-^2 s^-^1)');      xlabel('u_{* night}(m s^-^1)');
            set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',ust_class(1:2:end),...
                'XAxisLocation','bottom','XTick',x3,'Xlim',([ust_class(1) ust_class(end)]),...
                'FontSize',12);                                  axis tight;
            subplot(2,2,3);
            for ik=1:length(ta_class)-1
                plot(ust_class(2:end),((LUT_ta_ust_count(ik,:))'),'color',color_mtx(ik,:)); hold on;
                plot(ust_class(2:end),((LUT_ta_ust_count(ik,:))'),'.','MarkerSize',12,'color',color_mtx(ik,:));
            end
            title(['T_{air} - u_{*} - NEE'],'FontWeight','normal');
            ylabel('NEE_{night_{count}}(\mumol m^-^2 s^-^1)'); xlabel('u*_{night}(m s^-^1)');
            set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',ust_class(1:2:end),...
                'XAxisLocation','bottom','XTick',x3,'Xlim',([ust_class(1) ust_class(end)]),...
                'FontSize',12);          axis tight;
            
            % ..............................................................
            max_ust_ta_tresh_flag=[];
            for ik=1:(length(ta_class))-1
                max_ust_ta_tresh_flag=[max_ust_ta_tresh_flag;(ust_class(1:end-1))'];
            end
            max_ust_ta_tresh_flag=max_ust_ta_tresh_flag.*(LUT_ta_ust./LUT_ta_ust);
            
            for ik=1:(length(ta_class))-1
                ix=LUT_ta_ust(ik,:);
                ix(isnan(ix))=[];
                if length(ix)>1
                    iz=max(ix);
                    mp=find(ix==nanmax(ix));  mp=mp(1); if mp>1; mp=mp-1; end
                    jp=nanmean(ix(mp:end));
                    iz=1:length(ix);
                    max_ust_ta_tresh_max=jp*0.95;%0.95 at Papale's paper
                    ind=find((ix<max_ust_ta_tresh_max)&(iz<=mp));
                    max_ust_ta_tresh_flag(ik,ind)=NaN;
                end
            end
            
            subplot(2,2,4);
            for ik=1:length(ta_class)-1
                plot(ust_class(2:end),((max_ust_ta_tresh_flag(ik,:))'),'color',color_mtx(ik,:)); hold on;
                plot(ust_class(2:end),((max_ust_ta_tresh_flag(ik,:))'),'.','MarkerSize',12,'color',color_mtx(ik,:));
            end
            title('flag treshold T_{air} - u_{*} - NEE','FontWeight','normal');
            ylabel('NEE_{night}(\mumol m^-^2 s^-^1)'); xlabel('u_{*night}(m s^-^1)');
            set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',x3,...
                'XAxisLocation','bottom','XTick',x3,'Xlim',([ust_class(1) ust_class(end)]),...
                'FontSize',12);          axis tight;
            
            max_ust_ta=NaN(1,((length(ta_class))-1));
            for ik=1:(length(ta_class))-1
                ix=max_ust_ta_tresh_flag(ik,:);
                ix(isnan(ix))=[];     %if you have only one u* value missing, usually at low temp
                if ((~isempty(ix))&(length(ix))>1)
                    max_ust_ta(ik)=nanmin(ix);
                else
                    max_ust_ta(ik)=NaN;
                end
            end
            nanmean(max_ust_ta(max_ust_ta>0))
            ust_thr(im).K67.tasp=nanmean(max_ust_ta(max_ust_ta>0));
            
            ta_night=ta_night_all;      NEE_night=NEE_night_all;	ust_night=ust_night_all;
        end
        % ..............................................................
        % ..............................................................
        %% Nigth time respiration vs temperature
        % ..............................................................
        % ..............................................................
        %LUT month(as proxy for soil moisture) Ta night-time K67.nee
        [ta_class,NEE_class_avg,NEE_class_std,NEE_class_count]=...
            AM_class_eq_size(ta_night,NEE_night_ust,7);
        ta_class=[ta_class(1,1);ta_class(:,2)];
        M_class=(1:1:13)';
        [LUT_ta_month,LUT_ta_month_std,LUT_ta_month_count] =...
            AM_LUTgeneral(K67.tasp,K67.M,NEE_night,...
            ta_class,M_class,1,'T_{air}(C)','month','NEE (\mumol m^-2 s^-^1)',min_sample_size);
        
        %% ..............................................................
                x3 = (ta_class(1):(ta_class(end)-ta_class(1))/5:ta_class(end));
                color_mtx = jet(length(ta_class)-1);
        figure('color','white');
        subplot(2,2,1);
        for ik=1:length(ta_class)-1
            plot(ta_class(2:end),((LUT_ta_month(:,ik))'),'color',color_mtx(ik,:)); hold on;
        end
        title('T_{air}	NEE','FontWeight','normal');
        ylabel({'NEE_{night}';'(\mumol m^-^2 s^-^1)'});  xlabel('T_{air night}(\circC)');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0],...
            'XAxisLocation','bottom','Xlim',([ta_class(1) ta_class(end)]),'XTick',x3,...
            'FontSize',12);                          axis tight;
        
        % ..............................................................
        % ..............................................................
        %% Nigth time K67.ustar vs. month
        % ..............................................................
        % ..............................................................
        % month(as proxy for soil moisture) Ta night-time K67.nee
        
        [ust_class]=...
            AM_class_eq_size(ust_night,NEE_night,12);
        ust_class=[ust_class(1,1);ust_class(:,2)];
        month_class=(1:1:13)';
        [LUT_month_ust,LUT_month_ust_std,LUT_month_ust_count] =...
            AM_LUTgeneral(K67.M,ust_night,NEE_night,month_class,...
            ust_class,1,'month','u_{*}(m s^-^1)','NEE (\mumol m^-2 s^-^1)',min_sample_size);
        
        %% ..............................................................
        color_mtx = jet(length(month_class)-1);
        x3 = (ust_class(1):(ust_class(end)-ust_class(1))/5:ust_class(end));
        
        figure('color','white');
        subplot(2,2,1);
        for ik=1:length(month_class)-1
            plot(ust_class(2:end),((LUT_month_ust(ik,:))'),'color',color_mtx(ik,:)); hold on;
            plot(ust_class(2:end),((LUT_month_ust(ik,:))'),'.','MarkerSize',7,'color',color_mtx(ik,:));
            ix=(LUT_month_ust(ik,:));           iy=(1:length(ust_class)-1);
            ind=find((isnan(ix))|(isnan(iy)));  ix(ind)=[];         iz=iy;         iz(ind)=[];
            if length(iz)>2
                iq=interp1(iz,ix,iy);            LUT_month_ust(ik,:)=iq;
            end
%             plot(ust_class(2:end),((LUT_month_ust(ik,:))'),'o','MarkerSize',5,'color',color_mtx(ik,:));
        end
        title('month - u_{*} - NEE','FontWeight','normal');
        ylabel('NEE_{night}(\mumol m^-^2 s^-^1)');              xlabel('u_{*night}(m s^-^1)');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',x3,...
            'XAxisLocation','bottom','XTick',x3,'Xlim',([ust_class(1) ust_class(end)]),...
            'FontSize',12);          axis tight;
        subplot(2,2,2);
        for ik=1:length(month_class)-1
            plot(ust_class(2:end),((LUT_month_ust_std(ik,:))'),'color',color_mtx(ik,:)); hold on;
        end
        title('month - u_{*} - NEE','FontWeight','normal');
        ylabel('NEE_{night_{std}}(\mumol m^-^2 s^-^1)');           xlabel('u_{*night}(m s^-^1)');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',x3,...
            'XAxisLocation','bottom','Xlim',([ust_class(1) ust_class(end)]),...
            'FontSize',12);          axis tight;
        subplot(2,2,3);
        for ik=1:length(month_class)-1
            plot(ust_class(2:end),((LUT_month_ust_count(ik,:))'),'color',color_mtx(ik,:)); hold on;
        end
        title('month - u_{*} - NEE','FontWeight','normal');
        ylabel('NEE_{night_{count}}(\mumol m^-^2 s^-^1)');         xlabel('u_{*night}(m s^-^1)');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',x3,...
            'XAxisLocation','bottom','Xlim',([ust_class(1) ust_class(end)]),...
            'FontSize',12);          axis tight;
        
        % ..............................................................
        max_ust_month_tresh_flag=[];
        for ik=1:length(month_class)-1
            max_ust_month_tresh_flag=[max_ust_month_tresh_flag;(ust_class(2:end))'];
        end
        max_ust_month_tresh_flag=max_ust_month_tresh_flag.*(LUT_month_ust./LUT_month_ust);
        
        for ik=1:length(month_class)-1
            ix=LUT_month_ust(ik,:);   ix(isnan(ix))=[];
            if length(ix)>1
                max_ust_month_tresh_max=(mean(ix(end-1:end)))*0.95;%0.95 at Papale's paper
                max_ust_month_tresh_flag(ik,(LUT_month_ust(ik,:)<max_ust_month_tresh_max))=NaN;
            end
        end
        
        subplot(2,2,4);
        for ik=1:length(month_class)-1
            plot(ust_class(2:end),((max_ust_month_tresh_flag(ik,:))'),'color',color_mtx(ik,:)); hold on;
            plot(ust_class(2:end),((max_ust_month_tresh_flag(ik,:))'),'.','MarkerSize',12,'color',color_mtx(ik,:));
        end
        title('flag treshold month - u_{*} - NEE','FontWeight','normal');
        ylabel('NEE_{night}(\mumol m^-^2 s^-^1)'); xlabel('u_{*night}(m s^-^1)');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0],'XTick',x3,...
            'XAxisLocation','bottom','FontSize',12);          axis tight;
        
        max_ust_month=NaN(1,(length(month_class)-1));
        for ik=1:length(month_class)-1
            ix=max_ust_month_tresh_flag(ik,:); ix(isnan(ix))=[];     %if you have only one u* value missing, usually at low temp
            if length(ix)>1
                max_ust_month(ik)=min(ix);
            else
                max_ust_month(ik)=NaN;
            end
        end
        nanmean(max_ust_month(max_ust_month>0))
        ust_thr_month = nanmean(max_ust_month(max_ust_month>0));
        
        % %         if flag_figure_ps==1
        % %             for ik=1:13
        % %                 print('-dpsc', ['BFflux.ustr.ps' ], '-append');     close;
        % %             end
        % %         else
        % %             close all;
        % %         end
        
        %             close all;
        % ..............................................................
        %% Bootstrap this thing. 1) join both time series 2) calculate u*
        % ..............................................................
        ust_thr_btp10_per95   = [];   ust_thr_btp100_per95  = [];
        ust_thr_btp1000_per95 = [];   ust_thr_btp100_per99  = [];
        ust_thr_btp1000_per99 = [];   %ust_thr_btp100_per91 = [];
        %ust_thr_btp1000_per91 = [];
        
        ta_night = ta_night_all;      NEE_night = NEE_night_all;	
        ust_night = ust_night_all;    ust_night(ust_night==0) = NaN;
        for im=1:12
%             id_month=(im*3)-2:(im*3);
%             ta_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];	
%             NEE_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];
%             ust_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];

            ta_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];	
            NEE_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];
            ust_night((K67.M<id_month(1))|(K67.M>id_month(3)))=[];
            
            ta=ta_night;     ust=ust_night;        nee=NEE_night;
            ind=find((isnan(ta_night))|(isnan(ust_night))|(isnan(NEE_night)));
            ta(ind)=[];      ust(ind)=[];          nee(ind)=[];
            
            if (~isempty(ta))
%                 AM_ustar_papale(ta_night,ust_night,NEE_night,no_class_ta,...
%                     no_class_ust,flag_print,perc_flux_capture,flag_site)
                ust_thr_btp10(im).per95   = bootstrp(10,'AM_ustar_papale',ta,ust,nee,7,10,0,0.95,ip);
                ust_thr_btp100(im).per95  = bootstrp(100,'AM_ustar_papale',ta,ust,nee,7,10,0,0.95,ip);
                ust_thr_btp1000(im).per95 = bootstrp(10000,'AM_ustar_papale',ta,ust,nee,7,10,0,0.95,ip);
                ust_thr_btp100(im).per99  = bootstrp(100,'AM_ustar_papale',ta,ust,nee,7,10,0,0.99,ip);
                ust_thr_btp1000(im).per99 = bootstrp(10000,'AM_ustar_papale',ta,ust,nee,7,10,0,0.99,ip);
                %                     ust_thr_btp100(im).per91= bootstrp(100,'AM_ustar_papale',ta,ust,z,7,10,0,0.91,ip);
                %                     ust_thr_btp1000(im).per91= bootstrp(1000,'AM_ustar_papale',ta,ust,z,7,10,0,0.91,ip);
                %                     ust_thr_btp100.per85= bootstrp(100,'AM_ustar_papale',ta_night,ust_night,NEE_night,sunrise,sunset,1,0.85,ip);
                %                     ust_thr_btp1000.per85= bootstrp(1000,'AM_ustar_papale',ta_night,ust_night,NEE_night,sunrise,sunset,1,0.85,ip);
                close all;
                
                ust_thr_btp10_per95   = [ust_thr_btp10_per95;ust_thr_btp10(im).per95];
                ust_thr_btp100_per95  = [ust_thr_btp100_per95;ust_thr_btp100(im).per95];
                ust_thr_btp1000_per95 = [ust_thr_btp1000_per95;ust_thr_btp1000(im).per95];
                %                                         ust_thr_btp100_per91=[ust_thr_btp100_per91;ust_thr_btp100(im).per91];
                %                                         ust_thr_btp1000_per91=[ust_thr_btp1000_per91;ust_thr_btp1000(im).per91];
                ust_thr_btp100_per99  = [ust_thr_btp100_per99;ust_thr_btp100(im).per99];
                ust_thr_btp1000_per99 = [ust_thr_btp1000_per99;ust_thr_btp1000(im).per99];
                
            end
            ta_night = ta_night_all;      NEE_night = NEE_night_all;	ust_night = ust_night_all;
        end
        %%
        figure('color','white');
        subplot(2,2,1);
        histfit(ust_thr_btp10_per95);    h = get(gca,'Children');
        set(h(2),'FaceColor',[.8 .8 1]);     axis ([0 0.4 0 10]);
        xlabel('u_{* treshhold}');           ylabel({'incidence';'n=100, 95% max flux'});
        [muhat.btp10.per95,sigmahat.btp10.per95,muci.btp10.per95,...
            sigmaci.btp10.per95]= normfit(ust_thr_btp10_per95);
        %returns 95% confidence intervals for the rgsameter estimates on the mean
        title([site_name_BF]);
        
        subplot(2,2,2);
        histfit(ust_thr_btp100_per95,15);   h = get(gca,'Children');
        set(h(2),'FaceColor',[.8 .8 1]);        axis ([0.0 0.4 0 40]);
        xlabel('u_{* treshhold}');              ylabel({'incidence';'n=1000, 95% max flux'});
        [muhat.btp100.per95,sigmahat.btp100.per95,muci.btp100.per95,...
            sigmaci.btp100.per95]= normfit(ust_thr_btp100_per95);
        title([site_name_BF]);
        
        subplot(2,2,3);
        histfit(ust_thr_btp100_per99,15);	h = get(gca,'Children');
        set(h(2),'FaceColor',[.8 .8 1]);        axis ([0.0 0.4 0 40]);
        xlabel('u_{* treshhold}');              ylabel({'incidence';'n=100, 99% max flux'});
        [muhat.btp100.per99,sigmahat.btp100.per99,muci.btp100.per99,...
            sigmaci.btp100.per99]= normfit(ust_thr_btp100_per99);
        title([site_name_BF]);
        
        %             subplot(2,2,4);
        %             histfit(ust_thr_btp100_per91,15);	h = get(gca,'Children');
        %             set(h(2),'FaceColor',[.8 .8 1]);    axis ([0.0 04 0 40]);
        %             xlabel('u_{* treshhold}'); ylabel({'incidence';'n=100, 91% max flux'});
        %             [muhat.btp100.per91,sigmahat.btp100.per91,muci.btp100.per91,...
        %                 sigmaci.btp100.per91]= normfit(ust_thr_btp100_per91);
        %             %%
        figure('color','white');
        subplot(2,2,1);
        histfit(ust_thr_btp1000_per95,15);	h = get(gca,'Children');
        set(h(2),'FaceColor',[.8 .8 1]);        axis ([0 0.4 0 500]);
        xlabel('u_{* treshhold}');                 ylabel({'incidence';'n=1000, 95% max flux'});
        [muhat.btp1000.per95,sigmahat.btp1000.per95,muci.btp1000.per95,...
            sigmaci.btp1000.per95]= normfit(ust_thr_btp1000_per95);
        title([site_name_BF]);
        
        subplot(2,2,2);
        histfit(ust_thr_btp1000_per99,15);	h = get(gca,'Children');
        set(h(2),'FaceColor',[.8 .8 1]);        axis ([0 0.4 0 500]);
        xlabel('u_{* treshhold}');              ylabel({'incidence';'n=1000, 99% max flux'});
        [muhat.btp1000.per99,sigmahat.btp1000.per99,muci.btp1000.per99,...
            sigmaci.btp1000.per99]= normfit(ust_thr_btp1000_per99);
        title([site_name_BF]);
        
        %             subplot(2,2,3);
        %             histfit(ust_thr_btp1000_per91,15);	h = get(gca,'Children');
        %             set(h(2),'FaceColor',[.8 .8 1]);        axis ([0.1 0.8 0 500]);
        %             xlabel('u_{* treshhold}');                 ylabel({'incidence';'n=1000, 91% max flux'});
        %             [muhat.btp1000.per91,sigmahat.btp1000.per91,muci.btp1000.per91,...
        %                 sigmaci.btp1000.per91]= normfit(ust_thr_btp1000_per91);
        %
        %                 subplot(2,2,3);
        %                 histfit(ust_thr_btp100.per85,15);     h = get(gca,'Children');
        %                 set(h(2),'FaceColor',[.8 .8 1]);        axis ([0.1 0.4 0 50]);
        %                 xlabel('u_{* treshhold}'); ylabel({'incidence';'n=100, 85% max flux'});
        %                 [muhat.ust_thr_btp100.per85,sigmahat.ust_thr_btp100.per85,muci.ust_thr_btp100.per85,...
        %                     sigmaci.ust_thr_btp100.per85]= normfit(ust_thr_btp100.per85);
        %
        %         subplot(2,2,4);
        %         histfit(ust_thr_btp1000.per85,15);	h = get(gca,'Children');
        %         set(h(2),'FaceColor',[.8 .8 1]);        axis ([0.1 0.4 0 500]);
        %         xlabel('u_{* treshhold}');                 ylabel({'incidence';'n=1000, 85% max flux'});
        %         [muhat.ust_thr_btp1000.per85,sigmahat.ust_thr_btp1000.per85,muci.ust_thr_btp1000.per85,...
        %             sigmaci.ust_thr_btp1000.per85]= normfit(ust_thr_btp1000.per85);
    end
    %         save('matlabUSTuncertain.mat','ust_thr','muhat','sigmahat','muci','sigmaci')
end

%% ........................................................................
figure('color','white');
subplot(2,2,1);
ix = ust_thr_btp1000_per95;    ix(isnan(ix)) = [];
histfit(ust_thr_btp1000_per95,15);	h = get(gca,'Children');
[muhat.btp1000.per95,sigmahat.btp1000.per95,muci.btp1000.per95,...
    sigmaci.btp1000.per95] = normfit(ix);
set(h(1),'Color',[0 0 0]);
set(h(2),'FaceColor',[.8 .8 1],'EdgeColor',[.7 .7 1]);  axis ([0 0.2 0 800]);
text(.25, 250,char(site_BF));
ylabel({'incidence';'n=1000, 95% max flux'});           xlabel('u_{* treshhold}');

ix=[];           iz=[];             iw=[];
ix=([ix;muhat.btp1000.per95]);
iz=([iz;muci.btp1000.per95(:)]);
iw=([iw;sigmahat.btp1000.per95(:)]);
open ix;            open iz;            open iw;

subplot(2,2,2);
ix=ust_thr_btp1000_per99;    ix(isnan(ix))=[];
histfit(ust_thr_btp1000_per99,15);	h = get(gca,'Children');
[muhat.btp1000.per99,sigmahat.btp1000.per99,muci.btp1000.per99,...
    sigmaci.btp1000.per99] = normfit(ix)
set(h(1),'Color',[0 0 0]);
set(h(2),'FaceColor',[.8 .8 1],'EdgeColor',[.7 .7 1]);        axis ([0 0.2 0 800]);
text(.25, 250,char(site_BF));
ylabel({'incidence';'n=1000, 99% max flux'});   xlabel('u_{* treshhold}');

ix=[];   iz=[];   iw=[];
ix=([ix;muhat.btp1000.per99]);
iz=([iz;muci.btp1000.per99(:)]);
iw=([iw;sigmahat.btp1000.per99(:)]);
open ix;            open iz;            open iw;


%% ........................................................................
x1 = 0;         x2 = .6;        x3 = x1:((x2-x1)/4):x2;     x4 = x1:((x2-x1)/20):x2;
y1 = 0;         y2 = 20;        y3 = y1:((y2-y1)/4):y2;     y4 = y1:((y2-y1)/(length(K67.year)+2)):y2;

figure('color','white');
for im=1:12
    subplot(3,4,im);            hold on;
    title([month_name(im,:)],'FontName','arial','fontweight','normal');
    ta_night=K67.tasp;          NEE_night=K67.nee;         ust_night=K67.ustar;
    ta_night(K67.So>0)=NaN;     NEE_night(K67.So>0)=NaN;   ust_night(K67.So>0)=NaN;
                                NEE_night(K67.nee<-3)=NaN; 

                                
    ta_night(K67.M~=im)  = NaN;
    NEE_night(K67.M~=im) = NaN;
    ust_night(K67.M~=im) = NaN;
    
    plot(ust_night,NEE_night,'.','color',[.5 .5 .5],'MarkerSize',4);
    
    for in = 1:length(K67.year)
        ix = find(K67.Y==K67.year(in));
        x = ust_night(ix); y = NEE_night(ix);
        ind = find((isnan(x))|(isnan(y)));    x(ind) = [];y(ind) = [];
        %         [p,r2,~,~,sm] = AM_compute_regression(ust_night(ix),NEE_night(ix),color_mtx(im+2,:),0);
        %         y_model = polyval(p,ust_night);
        %         plot(ust_night,y_model,'-','color',color_mtx(in,:));
        [x_bin,y_bin] = AM_bin_noplot(x,y,20,5);
        f = @(b,x_bin) (b(1).* b(2).* x_bin)./(b(2) + (b(1).*x_bin));           % Rectangular hyperbola
        %         f = @(b,x_bin) b(1) + b(2)./(x_bin + b(3));
        %         f = @(b,x_bin) b(1).*exp(b(2).*x_bin);                                     % Objective Function
        B = fminsearch(@(b) norm(y_bin - f(b,x_bin)), [1; 1]);                  % Estimate Parameters
            plot(x_bin, f(B,x_bin), '-','color',color_mtx(in,:),'LineWidth',2)
        if im==8
            plot(x4(17),y4(in+1),'.','color',color_mtx(in,:),'MarkerSize',10);
            text (x4(18),y4(in+1),num2str(K67.year(in)),'FontSize',8);                            box on;
        end
    end
    box on;             grid on;
    set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],...
        'YTick',y3,'YTickLabel',y3,'Ylim',[y1 y2],'FontSize',12);
    if (im==1)|(im==5)|(im==9);  ylabel ({'NEE_{night} (\mumolCO_2 m^-^2 s^-^1)'},'FontSize',12); end
    if (im>8); xlabel ({'u_{* nighttime} (m s^-^1)'},'FontSize',12); end
end

figure('color','white');
for im=1:12
    subplot(3,4,im);            hold on;
    title([month_name(im,:)],'FontName','arial','fontweight','normal');
    ta_night=K67.tasp;          NEE_night=K67.nee;         ust_night=K67.ustar;
    ta_night(K67.So>0)=NaN;     NEE_night(K67.So>0)=NaN;   ust_night(K67.So>0)=NaN;
    NEE_night(K67.nee<-3)=NaN;
    
    
    ta_night(K67.M~=im)  = NaN;
    NEE_night(K67.M~=im) = NaN;
    ust_night(K67.M~=im) = NaN;
    
    plot(ust_night,NEE_night,'.','color',[.5 .5 .5],'MarkerSize',4);
    
    for in = 1:length(K67.year)
        ix = find(K67.Y==K67.year(in));
        x = ust_night(ix); y = NEE_night(ix);
        ind = find((isnan(x))|(isnan(y)));    x(ind) = [];y(ind) = [];
        %         [p,r2,~,~,sm] = AM_compute_regression(ust_night(ix),NEE_night(ix),color_mtx(im+2,:),0);
        %         y_model = polyval(p,ust_night);
        %         plot(ust_night,y_model,'-','color',color_mtx(in,:));
        [x_bin,y_bin] = AM_bin_noplot(x,y,20,5);
        plot(x_bin,y_bin,'.','color',color_mtx(in,:),'MarkerSize',12);
        f = @(b,x_bin) (b(1).* b(2).* x_bin)./(b(2) + (b(1).*x_bin));           % Rectangular hyperbola
        %         f = @(b,x_bin) b(1) + b(2)./(x_bin + b(3));
        %         f = @(b,x_bin) b(1).*exp(b(2).*x_bin);                                     % Objective Function
        B = fminsearch(@(b) norm(y_bin - f(b,x_bin)), [1; 1]);                  % Estimate Parameters
            plot(x_bin, f(B,x_bin), '-','color',color_mtx(in,:),'LineWidth',1)
        if im==8
            plot(x4(17),y4(in+1),'.','color',color_mtx(in,:),'MarkerSize',10);
            text (x4(18),y4(in+1),num2str(K67.year(in)),'FontSize',8);
        end
        box on;
    end
    box on;             grid on;
    set(gca,'XTick',x3,'XTickLabel',x3,'Xlim',[x1 x2],...
        'YTick',y3,'YTickLabel',y3,'Ylim',[y1 y2],'FontSize',12);
    if (im==1)|(im==5)|(im==9);  ylabel ({'NEE_{night} (\mumolCO_2 m^-^2 s^-^1)'},'FontSize',12); end
    if (im>8); xlabel ({'u_{* nighttime} (m s^-^1)'},'FontSize',12); end
end

