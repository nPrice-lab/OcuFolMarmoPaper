function randTest(repoPath)
%% Randomization test for latency difference
%   Input: repoPath -path for the OcuFolMarmoPaper repo. Assume current
%   directory is the repo if no input.

    %% General parameters
    if nargin<1
        repoPath = cd;
    end
    
    subj = {'butch' 'brick' 'm1899'};
    subjS = {'Bu' 'Br' 'Ni'};
    pd = {'postsac' 'sftf' 'direction'};
    ext = 'condsCB';
    filePath = fullfile(repoPath,'data');
    figPath = fullfile(repoPath,'figures/');
    if ~exist(figPath, 'dir')
        mkdir(figPath)
    end
      %% randomization test  

      expt1(subj,pd{1},ext,filePath,figPath);
      expt3(subj,pd{3},ext,filePath,figPath);
end

function expt3(subj,pd,ext,filePath,figPath)

    for ss = 1:numel(subj)
        fname = strjoin({pd,subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        keep = d.keepMat(:,end);

        dt = round(d.eyeData.dt);
        d = flip(d); % flip values in 180 and 270 conditions

        % empirical latencies
        latEm = expt3LatencyEst(d,keep,30:dt:150,d.conds);   
        latEm = [latEm(1) latEm(3) latEm(2) latEm(4)]; % swap columns
        diffEm{1} = latEm(2)-latEm(1); % horizontal difference
        diffEm{2} = latEm(4)-latEm(3); % vertical difference

        % shuffle cond labels, randomized estimates of latency
        nr = 1000;
        latRn = zeros(nr,4);
        clear dta
        for hv = 1:2
            for ii = 1:nr
                conds = d.conds;
                if hv == 1
                    idx = find(ismember(conds.value(:,3),[1 3])); 
                else
                    idx = find(ismember(conds.value(:,3),[2 4]));
                end
                tmpIds = conds.Ids(ismember(conds.Ids,idx));
                idx2 = find(ismember(conds.Ids,idx));
                ir = randperm(numel(tmpIds));
                conds.Ids(idx2) = tmpIds(ir);    
                tmp = expt3LatencyEst(d,keep,30:dt:150,conds);
                latRn(ii,:) = tmp;
            end
            latRn = [latRn(:,1) latRn(:,3) latRn(:,2) latRn(:,4)]; % swap columns
            diffRn{1} = latRn(:,2)-latRn(:,1); % horizontal difference
            diffRn{2} = latRn(:,4)-latRn(:,3); % vertical difference
            
            
            figure;
            h = histfit(diffRn{hv});
            hold on
            xline(diffEm{hv},':');
            text(diffEm{hv}+0.05,max(ylim(gca))-10,'data');
            xline(prctile(diffRn{hv},97.5),':r');
            text(prctile(diffRn{hv},97.5)+0.05,max(ylim(gca))-10,'97.5%');

            % get p value
            dist = fitdist(diffRn{hv},'Normal');
            pVal = 1-cdf(dist,diffEm{hv});

            dta.lat(hv).diffEm = diffEm;
            dta.lat(hv).diffRn = diffRn;
            dta.lat(hv).pVal = pVal;

        end
        
        % for speed
        if mod(round(median(latEm)),2)==1
            med = round(median(latEm))+1;
        else
            med = round(median(latEm));
        end
        AnWin = round(med)+10:dt:round(med)+60;

        spdEm = expt3SpeedEst(d,keep,AnWin,d.conds);
        spdEm = [spdEm(:,1) spdEm(:,3) spdEm(:,2) spdEm(:,4)]; % swap columns
        
        % t test
        for hv = 1:2
            idx = (hv-1)*2+1:(hv-1)*2+2;
            [~,pval,ci,stats] = ttest(spdEm(:,idx(1)),spdEm(:,idx(2)));
            dta.spd(hv).tpval = pval;
            dta.spd(hv).ci = ci;
            dta.spd(hv).stats = stats;
            dta.spd(hv).AnWin = AnWin;
        end
        
        idx = ~isnan(spdEm);
        spdEm = [mean(spdEm(idx(:,1),1)), mean(spdEm(idx(:,2),2)), mean(spdEm(idx(:,3),3)), mean(spdEm(idx(:,4),4))];
        diffEm{1} = spdEm(1)-spdEm(2); % horizontal difference
        diffEm{2} = spdEm(3)-spdEm(4); % vertical difference
        
        spdRn = zeros(nr,4);
        for hv = 1:2
            for ii = 1:nr
                conds = d.conds;
                if hv == 1
                    idx = find(ismember(conds.value(:,3),[1 3])); 
                else
                    idx = find(ismember(conds.value(:,3),[2 4]));
                end
                tmpIds = conds.Ids(ismember(conds.Ids,idx));
                idx2 = find(ismember(conds.Ids,idx));
                ir = randperm(numel(tmpIds));
                conds.Ids(idx2) = tmpIds(ir);    
                tmp = expt3SpeedEst(d,keep,AnWin,conds);
                tmp = [tmp(:,1) tmp(:,3) tmp(:,2) tmp(:,4)]; % swap columns
                idx = ~isnan(tmp);
                tmp = [mean(tmp(idx(:,1),1)), mean(tmp(idx(:,2),2)), mean(tmp(idx(:,3),3)), mean(tmp(idx(:,4),4))];
                spdRn(ii,:) = tmp;
            end
            diffRn{1} = spdRn(:,1)-spdRn(:,2); % horizontal difference
            diffRn{2} = spdRn(:,3)-spdRn(:,4); % vertical difference
        
        
            figure;
            h = histfit(diffRn{hv});
            hold on
            xline(diffEm{hv},':');
            text(diffEm{hv}+0.05,max(ylim(gca))-10,'data');
            xline(prctile(diffRn{hv},97.5),':r');
            text(prctile(diffRn{hv},97.5)+0.05,max(ylim(gca))-10,'97.5%');

            % get p value
            dist = fitdist(diffRn{hv},'Normal');
            pVal = 1-cdf(dist,diffEm{hv});

            dta.spd(hv).diffEm = diffEm;
            dta.spd(hv).diffRn = diffRn;
            dta.spd(hv).pVal = pVal;

        end
        
        save(fullfile(figPath,[subj{ss},'_Expt3randRes.mat']),'dta');

    end
end

function expt1(subj,pd,ext,filePath,figPath)

    for ss = 1:numel(subj)
        fname = strjoin({pd,subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        keep = d.keepMat(:,end);

        dt = round(d.eyeData.dt);
        
        % shuffle cond labels, randomized estimates of latency
        nr = 1000;
        latRn = zeros(nr,size(d.conds.value,1));
        for ii = 1:nr
            numTrials = size(d.eyeData.x2,1);
            ir = randperm(numTrials);
            conds = d.conds;
            conds.Ids = conds.Ids(ir);
            [~,par,R2] = d.LatencyEst(d.eyeData,keep,30:dt:150,'dx',conds);
            latRn(ii,:) = par(:,3);            
        end
        
        % fit empirical latency
        [~,par,R2] = d.LatencyEst(d.eyeData,keep,30:dt:150,'dx',d.conds);
        latEm = par(:,3)';
%         [fitEm,parEm,R2Em] = FitReg(latEm,[-1 mean(lat,'all')],d.conds.value');
        aaEm = polyfit(d.conds.value,latEm,1);
        yyEm = polyval(aaEm,d.conds.value);
        rrEm = corrcoef(yyEm,latEm);
        rrEm = rrEm(2);
        
        
        figure;
        plot(d.conds.value,latEm);
        hold on
%         plot(d.conds.value,fitEm);
        plot(d.conds.value,yyEm);
        
        % fit randomized latency
%         [~,par,R2] = FitReg(lat,[-1 mean(lat,'all')],d.conds.value');
%         a = par(:,1);
        aaRn = zeros(size(latRn,1),2);
        yyRn = zeros(size(latRn));
        rrRn = zeros(size(latRn,1),1);
        for ii = 1:size(latRn,1)
            aaRn(ii,:) = polyfit(d.conds.value,latRn(ii,:),1);
            yyRn(ii,:) = polyval(aaRn(ii,:),d.conds.value);
            tmp = corrcoef(yyRn(ii,:),latRn(ii,:));
            rrRn(ii) = tmp(2);
        end
        figure;
        h = histfit(aaRn(:,1));
        hold on
        xline(aaEm(:,1),':');
        text(aaEm(:,1)+0.05,max(ylim(gca))-10,'data');
        xline(prctile(aaRn(:,1),97.5),':r');
        text(prctile(aaRn(:,1),97.5)+0.05,max(ylim(gca))-10,'97.5%');
        
        % get p value
        dist = fitdist(aaRn(:,1),'Normal');
        pVal = 1-cdf(dist,aaEm(:,1));
        
        clear dta
        dta.lat.slopeEm = aaEm(:,1);
        dta.lat.rEm = rrEm;
        dta.lat.slopeRn = aaRn(:,1);
        dta.lat.rRn = rrRn;
        dta.lat.pVal = pVal;
%         save(fullfile(figPath,[subj{ss},'_randRes.mat']),'dta');
        
        % for speed
        % different window for different condition
        clear eyeT tmp
        for ii = 1:size(d.conds.value,1)
            tmp = eyeTraceWithCond(d,keep,[latEm(ii)+10:latEm(ii)+60],d.conds,'dx');
            eyeT(ii) = tmp(ii);
        end
        % normal linear regression
        for ii = 1:numel(eyeT)
            eyeT(ii).resX = ones(size(eyeT(ii).y,1),1)*d.conds.value(ii)/1000;
            eyeT(ii).resY = mean(eyeT(ii).y,2);
        end
        resX = cat(1,eyeT.resX);
        resX = [ones(size(resX)) resX];
        resY = cat(1,eyeT.resY);
        [b,bint,~,~,stats] = regress(resY,resX);
        dta.spd.b2 = b;
        dta.spd.bint2 = bint;
        dta.spd.stat2 = stats;

        % randomization
        spdEm = zeros(1,numel(eyeT));
        for ii = 1:numel(eyeT)
            spdEm(ii) = mean(eyeT(ii).mu);
        end
        
        spdRn = zeros(nr,size(d.conds.value,1));
        for ii = 1:nr
            numTrials = size(d.eyeData.x2,1);
            ir = randperm(numTrials);
            conds = d.conds;
            conds.Ids = conds.Ids(ir);
            [~,par,R2] = d.LatencyEst(d.eyeData,keep,30:dt:150,'dx',conds);
            latRn = par(:,3);  
            for jj = 1:size(d.conds.value,1)
                tmp = eyeTraceWithCond(d,keep,[latRn(jj)+10:latRn(jj)+60],conds,'dx');
                spdRn(ii,jj) = mean(tmp(jj).mu);
            end
        end
        
        aaEm = polyfit(d.conds.value/1000,spdEm,1);
        yyEm = polyval(aaEm,d.conds.value/1000);
        rrEm = corrcoef(yyEm,spdEm);
        rrEm = rrEm(2);
        
        aaRn = zeros(size(spdRn,1),2);
        yyRn = zeros(size(spdRn));
        rrRn = zeros(size(spdRn,1),1);
        for ii = 1:size(spdRn,1)
            aaRn(ii,:) = polyfit(d.conds.value/1000,spdRn(ii,:),1);
            yyRn(ii,:) = polyval(aaRn(ii,:),d.conds.value/1000);
            tmp = corrcoef(yyRn(ii,:),spdRn(ii,:));
            rrRn(ii) = tmp(2);
        end
        
        % get p value
        dist = fitdist(aaRn(:,1),'Normal');
        pVal = 1-cdf(dist,aaEm(:,1));

        dta.spd.slopeEm2 = aaEm(:,1);
        dta.spd.rEm2 = rrEm;
        dta.spd.slopeRn2 = aaRn(:,1);
        dta.spd.rRn2 = rrRn;
        dta.spd.pVal2 = pVal;        
        %% use the median window
        if mod(round(median(latEm)),2)==1
            med = round(median(latEm))+1;
        else
            med = round(median(latEm));
        end
        AnWin = round(med)+10:dt:round(med)+60;
        tmp = eyeTraceWithCond(d,keep,AnWin,d.conds,'dx');
        
        % normal linear regression
        for ii = 1:numel(tmp)
            tmp(ii).resX = ones(size(tmp(ii).y,1),1)*d.conds.value(ii)/1000;
            tmp(ii).resY = mean(tmp(ii).y,2);
        end
        resX = cat(1,tmp.resX);
        resX = [ones(size(resX)) resX];
        resY = cat(1,tmp.resY);
        [b,bint,~,~,stats] = regress(resY,resX);
        dta.spd.b = b;
        dta.spd.bint = bint;
        dta.spd.stat = stats;
        dta.spd.AnWin = AnWin;
        
        spdEm = mean(cat(1,tmp.mu),2);
        
        spdRn = zeros(nr,size(d.conds.value,1));
        for ii = 1:nr
            numTrials = size(d.eyeData.x2,1);
            ir = randperm(numTrials);
            conds = d.conds;
            conds.Ids = conds.Ids(ir);
            tmp = eyeTraceWithCond(d,keep,AnWin,conds,'dx');
            spdRn(ii,:) = mean(cat(1,tmp.mu),2);
        end
        
        aaEm = polyfit(d.conds.value/1000,spdEm,1);
        yyEm = polyval(aaEm,d.conds.value/1000);
        rrEm = corrcoef(yyEm,spdEm);
        rrEm = rrEm(2);
        
        aaRn = zeros(size(spdRn,1),2);
        yyRn = zeros(size(spdRn));
        rrRn = zeros(size(spdRn,1),1);
        for ii = 1:size(spdRn,1)
            aaRn(ii,:) = polyfit(d.conds.value/1000,spdRn(ii,:),1);
            yyRn(ii,:) = polyval(aaRn(ii,:),d.conds.value/1000);
            tmp = corrcoef(yyRn(ii,:),spdRn(ii,:));
            rrRn(ii) = tmp(2);
        end
        
        % get p value
        dist = fitdist(aaRn(:,1),'Normal');
        pVal = 1-cdf(dist,aaEm(:,1));
        
        dta.spd.slopeEm = aaEm(:,1);
        dta.spd.rEm = rrEm;
        dta.spd.slopeRn = aaRn(:,1);
        dta.spd.rRn = rrRn;
        dta.spd.pVal = pVal;
        save(fullfile(figPath,[subj{ss},'_Expt1randRes.mat']),'dta');
    end
end

function spd = expt3SpeedEst(d,keep,win,conds)
        dtaX = eyeTraceWithCond(d,keep,win,conds,'dx');
        dtaY = eyeTraceWithCond(d,keep,win,conds,'dy');
        
        for cng = 1:4
            cngId = find(conds.value(:,3)==cng);
%             flipId = find(d.conds.value(:,3)==cng & ismember(d.conds.value(:,1),[180 270]));
% 
%             for ii = 1:numel(flipId)
%                 dtaY(flipId(ii)).y = dtaY(flipId(ii)).y*-1;
%                 dtaX(flipId(ii)).y = dtaX(flipId(ii)).y*-1;
%             end

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


        spd = nan(max(cellfun(@numel,DV)),4);
        for ii = 1:4
            spd(1:numel(DV{ii}),ii) = DV{ii};
        end
        
end

function lat = expt3LatencyEst(d,keep,win,conds)
        dtaX = eyeTraceWithCond(d,keep,win,conds,'dx');
        dtaY = eyeTraceWithCond(d,keep,win,conds,'dy');
%         flipId = find(ismember(conds.value(:,1),[180 270]));
%         for ii = 1:numel(flipId)
%             dtaY(flipId(ii)).mu = dtaY(flipId(ii)).mu*-1;
%             dtaX(flipId(ii)).mu = dtaX(flipId(ii)).mu*-1;
%         end
        lat = zeros(4,1);
        for cng = 1:4
             idx = find(ismember(conds.value(:,3),cng));
             if mod(cng,2) == 1
                mu = cat(1,dtaX(idx).mu);
             else
                mu = cat(1,dtaY(idx).mu);
             end
             mu  = mean(mu,1);
             t = dtaX(1).t;
             [~,par,~] = analysis.src.FitPieceWiseLinearInv(mu,[-1 mean(mu) median(t)],t);
             lat(cng) = par(3);
        end
end

function d = flip(d)
    flipDir = [180 270];
    for ii = 1:numel(flipDir)
        flipId = find(ismember(d.conds.value(:,1),flipDir(ii)));
        flipTrial = ismember(d.conds.Ids,flipId);
        if ii == 1
            dv = {'x' 'dx' 'ddx'};
        else
            dv = {'y' 'dy' 'ddy'};
        end
        for kk = 1:numel(dv)
            tmp = d.eyeData.(dv{kk});
            tmp(flipTrial,:) = tmp(flipTrial,:).*-1;
            d.eyeData.(dv{kk}) = tmp;
        end
    end
end

function dta = eyeTraceWithCond(d,keep,win,conds,key)
    t = round(d.eyeData.t);
    tWin = t>=round(win(1)) & t<=round(win(end));
    t = t(tWin);
    
    dv = d.eyeData.(key);
    for kk = 1:size(conds.value,1)
        ix = keep & conds.Ids == kk;
        y = dv(ix,tWin);
        dta(kk).t = t;
        dta(kk).y = y;
        dta(kk).mu = mean(y,1); 
        dta(kk).SEM = std(y,[],1)/sqrt(size(y,1)); 
    end
end