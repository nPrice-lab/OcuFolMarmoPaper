function expt1Plot(subj,subjS,pd,ext,filePath,figPath,p)
    
    fg = gcf;
    fg.Position = [0 0 1680 1050];
    set(fg,'Color', [1 1 1]);
    fontS = 6;
    p.pack('h',{1/3 1/3 1/3}); % 3 columns for 3 subjects +1 for legend
    for ss = 1:numel(subj)
        p(ss).pack('v',{6 47 47}) % 3 rows (names, eye trace, two axes plot)
    end

    for ss = 1:numel(subj)
        fname = strjoin({pd,subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        keep = d.keepMat(:,end);
        
        % convert negative speeds
        IdsLeft = d.th == 180;
        d.eyeData.dx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1;         

        dt = round(d.eyeData.dt);
        pSubj = p(ss); % select subject panel
        
        % average eye traces
        lineCr = parula;
        nline = numel(unique(d.ps));
        lineCr = lineCr(round(linspace(50,200,nline)),:);
        win = -50:dt:150;
        clear dta
        dta = eyeTrace(d,keep,win,'dx');
        
        pSubj(2).select();
        for ii = 1:numel(dta)
            src.shadedErrorBar(dta(ii).t,dta(ii).mu,dta(ii).SEM,'lineprops',{'color',lineCr(ii,:)});
            hold on
        end
        xlabel('Time (ms)');
        ylabel('Eye Speed (deg/s)');
        xline(0);
        yline(0);
        ylim([-6 12]);
        yticks(-4:2:10);
        xticks([0 50 100]);
        xlim([win(1) win(end)]);
%         title(['Monkey ',subjS{ss}]);
%         set(gca,'FontSize',12);
        
        
        % two axes plots
        
        % latency estimate
        [~,par,R2] = d.LatencyEst(d.eyeData,keep,30:dt:150,'dx',d.conds);
        lat = par(:,3);
        
        % latency bootstrap CI
        latSEM = d.LatencyEstCI(d.eyeData,keep,30:dt:150,'dx',d.conds);

        % Time window for extracting speed/direction in each trial
        % [latency+10 latency+60]
        % use the same time window (based on median latency across
        % conditions) for all conditions
        if mod(round(median(lat)),2)==1
            med = round(median(lat))+1;
        else
            med = round(median(lat));
        end
        AnWin = round(med)+10:dt:round(med)+60;
%         AnWin = round(med):dt:round(med)+50;
        
        % add dashed lines to figures
        hold on
        xline(AnWin(1),':')
        xline(AnWin(end),':')
%         fill([AnWin(1) AnWin(end) AnWin(end) AnWin(1)],reshape(repmat(ylim,2,1),4,1),[0.5 0.5 0.5],'FaceAlpha',0,'LineStyle','--')
        hold off
        
        
        clear dta
        for ii = 1:numel(lat)
            
%             % use different time window for different conditions
%             a = round(lat(ii))+10;
%             b = round(lat(ii))+60;
%             AnWin = a:dt:b;
            tmp =  eyeTrace(d,keep,AnWin,'dx');
            dta(ii).y = mean(tmp(ii).y,2);
            dta(ii).mu = mean(tmp(ii).mu,2);
            dta(ii).SEM = mean(tmp(ii).SEM);
            dta(ii).btSEM = std(bootstrp(1000,@mean,dta(ii).y));
        end
        
        mu = cat(1,dta.mu);
        SEM = cat(1,dta.SEM);
        btSEM = cat(1,dta.btSEM);
        
        pSubj(3).select();  
%         left_color = [.5 .5 .5];
%         right_color = [0 0 0];
        yyaxis left;
%         pp = plot(d.conds.value,lat,'o-','LineWidth',1,'MarkerSize',4);
%         pp.MarkerFaceColor = pp.Color;
        errorbar(d.conds.value,lat,latSEM,'o-','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','auto');
        ylabel('Latency (ms)');
        ylim([30 100]);
        yticks(40:10:90)
        yyaxis right;
%         errorbar(d.conds.value,mu,SEM*1.96,'o-','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','auto');
        errorbar(d.conds.value,mu,btSEM,'o-','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','auto');
        % plot(conds.value,y2);
        xlabel('Post-saccadic delay (ms)');
        ylabel('Eye Speed (deg/s)')
        ylim([0 6.5]);
        xlim([0 320]);
%         title(['Monkey ',subjS{ss}]);
%         ax = gca();
%         ax.YAxis(1).Color = left_color;
%         ax.YAxis(2).Color = right_color;
%         set(gca,'FontSize',12);
        
%         grid = grid+1;
    end
    
    % add panel numbers
    panel = {'A' 'B' 'C' 'D' 'E' 'F'};
    for ss = 1:numel(subj)
        for ii = 1:2
            p(ss,ii+1).select()
            idx = ss+numel(subj)*(ii-1);
            hold on
            yy = ylim;
            xx = xlim;
            text(xx(1),yy(end)+diff(yy)/12,panel{idx},'FontWeight','bold','HorizontalAlignment','right');
        end
    end

    % change background color to indicate different subjects
%     ax = subplot(nrow, ncol, 2);
%     left = ax.Position(1) -0.04;
%     width = ax.Position(3) + 0.07;
%     bottom = 0;
%     height = 1;
%     annotation(gcf,'rectangle',[left,bottom,width,height],'LineStyle','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha', 0.1);
%     
%     
    % add subj names
    for ss = 1:numel(subj)
        p(ss,1).select()
        yy = ylim;
        xx = xlim;
        text(diff(xx)*0.5,mean(yy),['Monkey ', subjS{ss}],'FontWeight','bold','HorizontalAlignment','center');
        axis off
    end
    
    % change margin to indicate different subjects
    p.de.margin = 4;
    for ss = 1:numel(subj)
        p(ss).marginleft = 15;
        p(ss,2).margintop = 7;
        p(ss,2).marginbottom = 10;
        p(ss,3).margintop = 10;
        p(ss,3).marginbottom = 7;
    end
    p.margin = [10 10 10 5];
    
    % global font size
    p.de.fontsize = fontS;
    
    pause(5);
    p.export(fullfile(figPath,'Expt1Results.tif'),'-rp');

end