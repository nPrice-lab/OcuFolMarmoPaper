function methodPlot(subj,subjS,pd,ext,filePath,figPath,p)

    fname = strjoin({pd,subj{1},ext},'_');
    load(fullfile(filePath,[fname,'.mat']),'d');
    keep = d.keepMat(:,end);
    % convert negative speeds
    IdsLeft = d.th == 180;
    d.eyeData.dx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1;   
    
    keep = keep & d.ps == 50; % use only postsac = 50 ms trials
    win = [-150 150];
    t = round(d.eyeData.t);
    tWin = t>=round(win(1)) & t<=round(win(end));
    t = t(tWin);
    dv = d.eyeData.dx;
%     fg = figure;
%     fg.Position = [0 0 1680/2 1050/2];
%     set(fg,'Color', [1 1 1]);
    fontS = 8;

    p.pack('v',{50,50})
    p(2).pack('h',{5 95})
    p(2,2).pack('v',{30 70})
    % stimulus
    p(2,2,1).select()
    plot([t(1) 0 0 t(end)],[0 0 30 30],'k','LineWidth',2);
    ylim([-5 35]);
    yticks([0 30]);
    xlim([min(win) max(win)]);
    ylabel('Stimulus');
%     title('Stimulus','FontSize',fontS,'FontWeight','bold');
%     set(gca,'FontSize',12);
    
    p(2,2,2).select()
    hold on
    numTrials = size(dv,1);
    for tr = 1:numTrials
        if keep(tr)
            plot(t,dv(tr,tWin), 'Color', [0.5 0.5 0.5 0.5]);
        end
    end
    
    mu = mean(dv(keep,:),1);
    plot(t,mu(tWin),'Color',[0 0 0],'LineWidth',2)
    hold off
    ylim([-5 15]);
    yticks([0:5:10]);
    yline(0);
    xline(0,':'); 
%     ylabel('Speed (deg/s)');
    ylabel('Eye')
    xlabel('Time From Motion Onset (ms)');
    xlim([min(win) max(win)]);
%     title('Eye','FontSize',fontS,'FontWeight','bold');
%     set(gca,'FontSize',12);
    
    p(2,1).select()
    yy = ylim;
    xx = xlim;
    tt = text(xx(1),mean(yy),'Speed (deg/s)','FontSize',fontS,'HorizontalAlignment','center');
    tt.Rotation = 90;
    axis off
    
    p.de.margin = 4;
    p(2,2,2).margintop = 10;
    p(2,2,1).marginbottom = 10;
    p.margin = [10 10 10 5];
    
    % global font size
    p.de.fontsize = fontS;
    
    pause(5);
    p.export(fullfile(figPath,'Method.tif'),'-rp');
end