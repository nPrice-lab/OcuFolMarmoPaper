function GratingDotDiff(repoPath)
%% examine the difference of eye speeds for grating and random dots with similar speed
%% General parameters
    if nargin<1
        repoPath = cd;
    end
    
    subj = {'butch' 'brick'};
    subjS = {'Bu' 'Br' 'Ni'};
    pd = {'postsac', 'sftf'};
    ext = 'condsCB';
    filePath = fullfile(repoPath,'data');
    figPath = fullfile(repoPath,'figures/');
    if ~exist(figPath, 'dir')
        mkdir(figPath)
    end
    
    
    for ss = 1:numel(subj)
        
        %% random dots
        fname = strjoin({pd{1},subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        IdsLeft = d.th == 180;
        d.eyeData.dx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1; 
        d.eyeData.ddx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1; 
        keep = d.keepMat(:,end);
        dt = round(d.eyeData.dt);
        
        % two axes plots
        
        % latency estimate
        [~,par,R2] = d.LatencyEst(d.eyeData,keep,30:dt:150,'dx',d.conds);
        lat = par(:,3);

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
         
        clear dta
            
        dta =  eyeTrace(d,keep,AnWin,'dx');
        dotRes = mean(dta(3).y,2); % extract ps = 50 condition
        
        %% random dots
        fname = strjoin({pd{2},subj{ss},ext},'_');
        load(fullfile(filePath,[fname,'.mat']),'d');
        IdsLeft = d.th == 180;
        d.eyeData.dx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1; 
        d.eyeData.ddx(IdsLeft,:) = d.eyeData.dx(IdsLeft,:)*-1; 
        keep = d.keepMat(:,end);
        dt = round(d.eyeData.dt);       
        
        win = 60:dt:110;
        dta = eyeTrace(d,keep,win,'dx');
        
        idx = d.conds.value(:,1)==0.62 & d.conds.value(:,2)==18.75; % speed = 30.24 deg
        graRes = mean(dta(idx).y,2); % extract 
        
        disp('mean speed for dots')
        disp(mean(dotRes));
        disp('SD speed for dots')
        disp(std(dotRes));
        disp('mean speed for gratings')
        disp(mean(graRes));
        disp('SD speed for gratings')
        disp(std(graRes));
        [~,p,~,stats] = ttest2(dotRes,graRes);
        disp(stats.tstat);
        disp(stats.df);
        disp(p);
    end
end