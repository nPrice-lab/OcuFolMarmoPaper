function otherLatEst(repoPath)
%%  try other latency estimation methods in experiment 1

%% General parameters
    if nargin<1
        repoPath = cd;
    end
    
    subj = {'butch' 'brick' 'm1899'};
    subjS = {'Bu' 'Br' 'Ni'};
    pd = 'postsac';
    ext = 'condsCB';
    filePath = fullfile(repoPath,'data');
    figPath = fullfile(repoPath,'figures/');
    if ~exist(figPath, 'dir')
        mkdir(figPath)
    end
    
    for ss = 1:numel(subj)
        %% Original method (piecewise regression)

        fname = strjoin({pd,subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        IdsLeft = d.th == 180;
        d.eyeData.dx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1; 
        d.eyeData.ddx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1; 
        keep = d.keepMat(:,end);
        dt = round(d.eyeData.dt);

        latOG = zeros(1,size(d.conds.value,1));
        [~,par,R2] = d.LatencyEst(d.eyeData,keep,30:dt:150,'dx',d.conds);
        lat = par(:,3);

        latOG(:) = lat;

        %% First wave of acceleration
        latAC = zeros(1,size(d.conds.value,1));

        win = 0:dt:150;
        dtaDDX = d.eyeTrace(keep,win,'ddx');
        thd = [40 20 10];
        for ii = 1:size(d.conds.value,1)

    %         mu = mean(dtaDDX(ii).mu(dtaDDX(ii).t<=40));
    %         SD = std(dtaDDX(ii).mu(dtaDDX(ii).t<=40));
    % 
    %         tmp = dtaDDX(ii).t(dtaDDX(ii).mu>(mu+2*SD));
    %         latAC(ss,ii) = tmp(1);    
    %         xline(tmp(1),'-.');

            tmp = dtaDDX(ii).t(dtaDDX(ii).mu>thd(ss));
            latAC(ii) = tmp(1); 
        end    

        %% 1.5 SD deviation from steady state

        latDE = zeros(1,size(d.conds.value,1));

        win = 20:dt:150;
        eos = 40; % end of steady state
        dtaDX = d.eyeTrace(keep,win,'dx');
        for ii = 1:size(d.conds.value,1)
            mu = mean(dtaDX(ii).mu(dtaDX(ii).t<=eos));
            SD = std(dtaDX(ii).mu(dtaDX(ii).t<=eos));

            tmp = dtaDX(ii).t(dtaDX(ii).mu>(mu+1.5*SD));
            tmp = tmp(tmp>eos);
            latDE(ii) = tmp(1);
        end
        
        %% plot figure
        p = panel();
        p.pack('v',{5,95});
        % subject Name
        p(1).select();
        yy = ylim;
        xx = xlim;
        text(diff(xx)*0.5,mean(yy),['Monkey ', subjS{ss}],'FontWeight','bold','HorizontalAlignment','center');
        axis off
        
        p(2).pack(2,3);
%         fg.Position = [0 0 1680 1050];
        
        win = -50:150;
        dtaDX = d.eyeTrace(keep,win,'dx');
        dtaDDX = d.eyeTrace(keep,win,'ddx');
        for ii = 1:size(d.conds.value,1)
            p(2,(ii>3)+1,ii-(ii>3)*3).select();

            yyaxis left
            p1 = plot(dtaDX(ii).t,dtaDX(ii).mu);
            ylim([-6 12])
            if ii == 1
                ylabel('Eye Speed (deg/s)')
            end
            
            yyaxis right
            p2 = plot(dtaDDX(ii).t,dtaDDX(ii).mu);
            ylim([-50 150])
            if ii == 1
                ylabel('Eye Acc (deg/s^2)')
            end
            
            p3 = xline(latOG(ii),'k');
            p4 = xline(latAC(ii),'-.');
            p5 = xline(latDE(ii),'--');
            
            xlim([win(1) win(end)])
            xlabel('Time (ms)')
%             set(findall(gcf,'-property','FontSize'),'FontSize',12)
            if ii == 3
                legend([p3 p4 p5],{'Piece-wise' 'Acc. Threshold' 'Deviation'},'Position',[0.8 0.9 0.1 0.05])
            end
                
        end
        p.de.margin = 12;
        p.de.fontsize = 6;
        p.export(fullfile(figPath,[subj{ss},'_otherLatEst.jpg']),'-rp');
%         exportgraphics(fg,fullfile(figPath,[subj{ss},'_otherLatEst.jpg']));
    end
    

end