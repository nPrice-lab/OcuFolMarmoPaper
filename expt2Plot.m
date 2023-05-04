function expt2Plot(subj,subjS,pd,ext,filePath,figPath,p)
    % row 1,2: 2*3 figures of eye speed vs time, grouped by tf
    % row 3: heatmap of eye speed at 3 time intervals
    % row 4: heatmap of gain at 3 time intervals
    
    fg = gcf;
    fg.Position = [0 0 1680 1050];
    set(fg,'Color', [1 1 1]);
    fontS = 6;
    p.pack('h',numel(subj)); % 2 columns for 2 subjects
    for ss = 1:numel(subj)
%         p(ss).pack('v',{6 47 47}) % 3 rows (names, eye trace, two axes plot)
        p(ss).pack('v',{4 32 32 32})
        p(ss,2).pack(2,{5 5 30 30 30})
        p(ss,3).pack(2,{5 5 30 30 30})
        p(ss,4).pack(2,{5 5 30 30 30})
    end

    for ss = 1:numel(subj)
        fname = strjoin({pd,subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        keep = d.keepMat(:,end);
        
        % convert negative speeds
        IdsLeft = d.th == 180;
        d.eyeData.dx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1;    
        
        dt = round(d.eyeData.dt);
        pSubj = p(ss);
        
        % average eye traces
        lineCr = colormap(colorcet('L05'));
        nline = numel(unique(d.sf));
        lineCr = lineCr(round(linspace(50,256,nline)),:);
        win = -50:dt:150;
        dta = eyeTrace(d,keep,win,'dx');
        sfCond = unique(d.conds.value(:,strcmp(d.conds.var,'sf')));
        tfCond = unique(d.conds.value(:,strcmp(d.conds.var,'tf')));
        
        for itf = 1:numel(tfCond)
            prow = 1 + 1*(itf>3);
            pcol = itf - 3*(prow-1) + 2;
            pSubj(2,prow,pcol).select();

            for isf = 1:numel(sfCond)
                idx = find(d.conds.value(:,1)==sfCond(isf) & d.conds.value(:,2)==tfCond(itf));
                src.shadedErrorBar(dta(1).t,dta(idx).mu,dta(idx).SEM,'lineprops',{'color',lineCr(isf,:)});
                hold on
            end
%             xlabel('Time (ms)');
%             ylabel('Eye Speed (deg/s)');
            title([num2str(tfCond(itf)),' Hz'])
            xline(0);
            yline(0);
            ylim([-6 12]);
            xlim([win(1) win(end)]);
            yticks(-4:4:10);
            xticks([0 150]);
            xtickangle(0);
            if ss == 1 && prow == 2 && pcol == 3
                xlabel('Time (ms)');
%                 ylabel('Eye Speed (deg/s)');           
            end
            if prow == 1
                xticks([])
            end
%             set(gca,'FontSize',12);
        end
        
        % heat maps & condition avg plot
        win = {60:dt:110,80:dt:130,100:dt:150};

        dvType = {'Speed','Gain'};
        for ii = 1:numel(win)
            for kk = 1: numel(dvType)
            AnWin = win{ii};
            dta = eyeTrace(d,keep,AnWin,'dx');
            map = zeros(numel(sfCond),numel(tfCond));
            for itf = 1:numel(tfCond)
                for isf = 1:numel(sfCond)
                    idx = find(d.conds.value(:,1)==sfCond(isf) & d.conds.value(:,2)==tfCond(itf));
                    if strcmp(dvType{kk},'Gain')
                        div = tfCond(itf)/sfCond(isf);
                    else
                        div = 1;
                    end
                    map(isf,itf) = mean(dta(idx).mu,'all')/div;
                end    
            end
            
            prow = kk;
            pcol = ii+2;
            pSubj(3,prow,pcol).select();
            colormap(hot)
            pcc = src.pcolorcent(sfCond',tfCond',map);
            
            %smooth
%             pcc.FaceColor = 'interp';
%             pcc.EdgeColor = 'none';
            
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
            xlim([pcc.XData(1) pcc.XData(end)])
            if prow == 2
                xticks(sfCond([2,4,6]));
                xticklabels({num2str(sfCond(2)) num2str(sfCond(4)) num2str(sfCond(6))});
                xtickangle(90)
            else
                xticks([]);
            end

            if pcol == 3
                yticks(tfCond([2,4]));
                yticklabels({num2str(tfCond(2)) num2str(tfCond(4))});
                ytickangle(0)
            else
                yticks([]);
            end
%             cr = colorbar;
            if strcmp(dvType{kk},'Gain')
                caxis([0 1]);
            else
                caxis([0 12]);
            end
            
%             if strcmp(dvType,'gain')
%                 cr.Label.String = "Mean Horizontal Gain";
%             else
%                 cr.Label.String = "Mean Horizontal Velocity (deg/s)";
%             end
%             set(gca,'FontSize',12);
            if kk ==1
                title([num2str(AnWin(1)),' - ',num2str(AnWin(end)),' ms']);
            end
            
%             yticks([1 10]);
%             xticks([0.1 1]);
            
            if ss == 1 && kk == 2 && ii == 1
%                 ylabel('Temporal Freq (Hz)');
                xlabel('Spatial Freq (cpd)');
            end
            
            %%%%%%%%%%%%%%%%%%%%%% condition average error bar plot
            pSubj(4,prow,pcol).select();
            lineCr = colormap(colorcet('L08'));
            nline = numel(unique(d.tf));
            lineCr = lineCr(round(linspace(50,256,nline)),:);
            
            for itf = 1:nline
                hold on
                plot(sfCond,map(:,itf),'Color',lineCr(itf,:));
            end
            set(gca, 'XScale', 'log')
            
            if strcmp(dvType{kk},'Gain')
                ylim([0 1.2]);
                yticks(0:0.5:1);
            else
                ylim([0 12]);
                yticks(0:5:10);
            end
            
            if prow == 2
                xticks(sfCond([2,4,6]));
                xticklabels({num2str(sfCond(2)) num2str(sfCond(4)) num2str(sfCond(6))});
                xtickangle(90)
            else
                xticks([]);
            end
            
            if ss == 1 && kk == 2 && ii == 1
%                 ylabel('Temporal Freq (Hz)');
                xlabel('Spatial Freq (cpd)');
            end
            
            end
        end
    end
    
    % add panel numbers
    panel = {'A' 'B'};
    for ss = 1:numel(subj)
        p(ss,2,1,3).select();
        hold on
        yy = ylim;
        xx = xlim;
        text(xx(1)-diff(xx)*0.15,yy(end)+diff(yy)*0.08,panel{ss},'FontSize',fontS+1,'FontWeight','bold');
    end
    
    panel = {'C' 'E'; 'D' 'F'};
    for ss = 1:numel(subj)
        for ii = 1:2
            p(ss,3,ii,3).select();
            hold on
            yy = ylim;
            xx = xlim;
            text(exp(log(xx(1))-diff(log(xx))*0.15),yy(end),panel{ss,ii},'FontSize',fontS+1,'FontWeight','bold');
        end
    end

    panel = {'G' 'I'; 'H' 'J'};
    for ss = 1:numel(subj)
        for ii = 1:2
            p(ss,4,ii,3).select();
            hold on
            yy = ylim;
            xx = xlim;
            text(exp(log(xx(1))-diff(log(xx))*0.15),yy(end),panel{ss,ii},'FontSize',fontS+1,'FontWeight','bold');
        end
    end    
    % change background color to indicate different subjects
%     ax = subplot(nrow, ncol, 4);
%     left = 0.51;
%     width = 0.5;
%     bottom = 0;
%     height = 1;
%     annotation(gcf,'rectangle',[left,bottom,width,height],'LineStyle','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha', 0.1);
%         
    
    % add subj names
    for ss = 1:numel(subj)
        p(ss,1).select()
        yy = ylim;
        xx = xlim;
        text(diff(xx)*0.5,mean(yy),['Monkey ', subjS{ss}],'FontSize',fontS+2,'FontWeight','bold','HorizontalAlignment','center');
        axis off
    end
    
    % add speed gain label
    lab  = {'Speed' 'Gain'};
    for ii = 1:2
        for jj = 1:2
            p(1,3+jj-1,ii,1).select()
            yy = ylim;
            xx = xlim;
            tt = text(xx(1),mean(yy),lab(ii),'FontSize',fontS+1,'FontWeight','bold','HorizontalAlignment','center');
            tt.Rotation = 90;
            axis off
        end
    end
 
    % add y-axis label
    lab  = {'Eye Speed (deg/s)' 'Temporal Freq (Hz)'};
    for ii = 1:2
        p(1,ii+1,2,2).select()
        yy = ylim;
        xx = xlim;
        if ii == 1
            xpos = xx(1);
        else
            xpos = xx(1)-mean(xx);
        end
        tt = text(xx(1)-mean(xx),mean(yy),lab(ii),'FontSize',fontS,'HorizontalAlignment','left');
        tt.Rotation = 90;
        axis off
    end    
    % change margin
    p.de.margin = 5;
    for ss = 1:numel(subj)
%         p(ss).marginleft = 15;
%         p(ss,2).margintop = 7;
        p(ss,2).marginbottom = 10;
        p(ss,3).margintop = 10;
        p(ss,3).marginbottom = 10;
        p(ss,4).margintop = 10;
%         p(ss,3).marginbottom = 7;
    end
%     p(1,3,2,2).marginright = 8;
    p.margin = [5 10 20 5];
    
    % global font size
    p.de.fontsize = fontS;
    
    colormap(hot)
    pause(5);
    p.export(fullfile(figPath,'Expt2Results.tif'),'-rp','-a1.2');
end