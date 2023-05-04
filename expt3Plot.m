function expt3Plot(subj,subjS,pd,ext,filePath,figPath,p)
    % upper 6 panels: 2*3 trial-averaged eye velocity traces for hor/ver 
    % congruent vs incongruent
    % bottom 6 panels: 2*3 violin plots of eye speed
    
    fg = gcf;
    fg.Position = [0 0 1680 1050];
    set(fg,'Color', [1 1 1]);
    fontS = 6;
    p.pack('h',{1/3 1/3 1/3}); % 3 columns for 3 subjects +1 for legend
    for ss = 1:numel(subj)
        p(ss).pack('v',{5 40 40 15}) % 4 rows (names, eye trace, speed violin, latency)
        p(ss,2).pack('h',{50 50}) % 2 columns (hor/ver)
        p(ss,3).pack('h',{50 50}) % 2 columns (hor/ver)
        p(ss,4).pack('h',{50 50}) % 2 columns (hor/ver)
    end

    for ss = 1:numel(subj)
        fname = strjoin({pd,subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        keep = d.keepMat(:,end);   
        
        dt = round(d.eyeData.dt);
        
        % plot eye traces
        lineCr = jet; % hot;
        nline = 2;
        lineCr = lineCr(round(linspace(50,200,nline)),:);
        win = -50:dt:150;
        dtaX = eyeTrace(d,keep,win,'dx');
        dtaY = eyeTrace(d,keep,win,'dy');
        HVtitle = {'Horizontal' 'Vertical'};
        
        % select subject 
        pSubj = p(ss);
        
        for hv = 1:2
            % select panel
            pSubj(2,hv).select()
            
            for cc = 1:2
                cng = hv + (cc>1)*2;
                cngId = find(d.conds.value(:,3)==cng);
                flipId = find(d.conds.value(:,3)==cng & ismember(d.conds.value(:,1),[180 270]));

                for ii = 1:numel(flipId)
                    dtaY(flipId(ii)).y = dtaY(flipId(ii)).y*-1;
                    dtaX(flipId(ii)).y = dtaX(flipId(ii)).y*-1;
                end

                switch cng
                    case {1 3}
                        dv = cat(1,dtaX(cngId).y);
                    case {2 4}
                        dv = cat(1,dtaY(cngId).y);
                end

                mu = mean(dv,1);
                SEM = std(dv,[],1)/sqrt(size(dv,1));
                src.shadedErrorBar(dtaY(1).t,mu,SEM,'lineprops',{'color',lineCr(cc,:)});
                hold on
            end
            
            xline(0);
            yline(0);
            ylim([-6 12]);
            yticks(-4:4:10);
            xticks([0 50 100])
            xlim([win(1) win(end)]);
            title({HVtitle{hv},' '},'FontWeight','bold');
            
            if hv == 1                
                xlabel('Time (ms)');
                ylabel('Eye Speed (deg/s)');
            else
%                 xticks([])
                yticks([])
            end
            
        end
        
        % violin plot
        DV = cell(1,4);
%         projMu = cell(1,4);
%         projSEM = zeros(1,4);

        % Time window for extracting speed/direction in each trial
        % [latency+10 latency+60]
        % use the same time window (based on median latency across
        % conditions) for all conditions
        win = 30:dt:150;
        dtaX = eyeTrace(d,keep,win,'dx');
        dtaY = eyeTrace(d,keep,win,'dy');
        flipId = find(ismember(d.conds.value(:,1),[180 270]));
        for ii = 1:numel(flipId)
            dtaY(flipId(ii)).mu = dtaY(flipId(ii)).mu*-1;
            dtaX(flipId(ii)).mu = dtaX(flipId(ii)).mu*-1;
        end
        lat = zeros(4,1);
        for cng = 1:4
             idx = find(ismember(d.conds.value(:,3),cng));
             if mod(cng,2) == 1
                mu = cat(1,dtaX(idx).mu);
             else
                mu = cat(1,dtaY(idx).mu);
             end
             mu  = mean(mu,1);
             t = dtaX(1).t;
             [~,par,~] = src.FitPieceWiseLinearInv(mu,[-1 mean(mu) median(t)],t);
             lat(cng) = par(3);
        end
        
        if mod(round(median(lat)),2)==1
            med = round(median(lat))+1;
        else
            med = round(median(lat));
        end
        AnWin = round(med)+10:dt:round(med)+60;
        
        
        % go back and add two dashed lines to eye traces plots
        for hv = 1:2
            % select horizontal/vertical column
            pSubj(2,hv).select();
            hold on
            xline(AnWin(1),':');
            xline(AnWin(end),':');
%             fill([AnWin(1) AnWin(end) AnWin(end) AnWin(1)],reshape(repmat(ylim,2,1),4,1),[0.5 0.5 0.5],'FaceAlpha',0,'LineStyle','--')
            hold off
        end
        dtaX = eyeTrace(d,keep,AnWin,'dx');
        dtaY = eyeTrace(d,keep,AnWin,'dy');
        
        for cng = 1:4
            cngId = find(d.conds.value(:,3)==cng);
            flipId = find(d.conds.value(:,3)==cng & ismember(d.conds.value(:,1),[180 270]));

            for ii = 1:numel(flipId)
                dtaY(flipId(ii)).y = dtaY(flipId(ii)).y*-1;
                dtaX(flipId(ii)).y = dtaX(flipId(ii)).y*-1;
            end

            switch cng
                case {1 3}
                    dv = cat(1,dtaX(cngId).y);
                case {2 4}
                    dv = cat(1,dtaY(cngId).y);
            end

            dv = mean(dv,2);
            DV{cng} = dv;
%             projMu(cng) = mean(y);
%             projSEM(cng) = std(y,[],1)/sqrt(size(y,1)); 
        end


        vio = nan(max(cellfun(@numel,DV)),4);
        for ii = 1:4
            vio(1:numel(DV{ii}),ii) = DV{ii};
        end
        vio = [vio(:,1) vio(:,3) vio(:,2) vio(:,4)]; % swap columns
        type = {'Con' 'Incon'};
        
        pSubj = p(ss);
        for hv = 1:2
            % select panel
            pSubj(3,hv).select();
            idx = (hv-1)*2+1:(hv-1)*2+2;
            dv = vio(:,idx);
            src.violinplot(dv, type,'ViolinColor',lineCr,'MarkerSize',7);
%             xlabel('Condition');
            ylim([-2 12]);
            yticks(-2:2:10);
            if hv==1
                ylabel('Eye Speed (deg/s)');
            else
                yticks([]);
            end
            
%             title(HVtitle{hv});
%             set(gca,'FontSize',12);
            set(gca,'box','off')

            % ttest
            [~,pval,ci,stats] = ttest(vio(:,idx(1)),vio(:,idx(2)));
            if pval<.001; as = '***'; xloc = 1.5;
            elseif pval<.01; as = '**'; xloc = 1.5;
            elseif pval<.05; as = '*'; xloc = 1.5;
            end
            if pval<.05
                text(xloc,10,as,'FontSize',fontS-1)
            end
            save(fullfile(figPath,[subj{ss},'_',HVtitle{hv},'_ttest.mat']),'pval','ci','stats');
        end
        
        
        % latency bar plots
        XL = xlim; % use xlim of violin plots
        latBar = [lat(1) lat(3) lat(2) lat(4)]; % swap columns
        for hv = 1:2
            % select panel
            pSubj(4,hv).select();
            idx = (hv-1)*2+1:(hv-1)*2+2;
            dv = latBar(idx);
            b = bar(1:2,dv,0.5,'EdgeColor',[.5 .5 .5]);
            b.FaceColor = 'flat';
            b.CData = lineCr;
            xticklabels(type);
            xlim(XL);
            ylim([0 60]);
            yticks(0:20:40);
            if hv==1
                ylabel('Latency (ms)');
            else
                yticks([]);
            end
        end
    end
    
    %     add legend
%     p(2,3,1,2).select()
%     lgd = legend({'Cong.','Incong.'});
%     lgd.Location = 'southeast';
%     lgd.Position = [1.2 0.5924 0.02 0.019];
%     lgd.FontSize = fontS-1;    
%     p(3,2,2).select()
%     for cc = 1:2
%         plot(dtaY(1).t, nan,'color',lineCr(cc,:)) % plot nans (hack to generate correct legend but plot no data)
%     end
%     lgd = legend({'Cong.','Incong.'});
%     lgd.Title.String = {'Direction','Congruency'};
%     lgd.Location = 'southeast';
%     lgd.FontSize = fontS-1;
%     p(4,2).select(lgd);
%     lgd.Title.FontSize = fontS-1;
%     axis off    
%     p(2,4).marginleft = 3;
%     p(2,4).margintop = 3;
    
    % add panel numbers
    panel = {'A' 'D' 'G'; 'B' 'E' 'H'; 'C' 'F' 'I'};
    for ss = 1:numel(subj)
        for ii = 2:3
            p(ss,ii,1).select()
            hold on
            yy = ylim;
            xx = xlim;
            text(xx(1),yy(end)+diff(ylim)/12,panel{ss,ii-1},'FontSize',fontS+1,'FontWeight','bold','HorizontalAlignment','right');
        end
        
        ii = 4;
        p(ss,ii,1).select()
        hold on
        yy = ylim;
        xx = xlim;
        text(xx(1),yy(end)+diff(ylim)/5,panel{ss,ii-1},'FontSize',fontS+1,'FontWeight','bold','HorizontalAlignment','right');
    end
    
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
    p.export(fullfile(figPath,'Expt3Results.tif'),'-rp');
end