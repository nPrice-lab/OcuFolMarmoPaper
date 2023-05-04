classdef OFMarmo
    
  properties 
      filtOrder = 3;        % filter order 
%       filtFramelen = 51;    % filter threshold 
      filtFramelen = 21;    % filter threshold (for 500Hz)
  end 

  properties
    % stimulus properties 
    fixX;           % fix location X
    fixY;           % fix location Y
    th;             % direction 
    Bth;            % bandwidth of direction
    sf;             % spatial frequency
    tf;             % temporal frequency
    Bsf;            % Bandwidth of spatial frequency
    V;              % speed
    Bv;             % Bandwidth of speed    
    ps;             % post saccadic delay 
    contrast;       % contrast      
    
    % stimulus timing 
    fOffset;        % fix offset time (seconds) 
    tOnset;         % tar onset time (seconds) 
    mOnset;         % motion onset time (seconds) 

    % experimental conditions
    conds;
    
    % trials to keep after selection
    keepMat;
    
    % paradigm: postsac, direction, grating, or motion cloud
    pd;

    % eye movement data 
    eyeData;        % eyelink Data, aligned to motion onset
     
    % screening trials 
    fDropRaw;          % true|false - trials with framdrop
    fDropAfterOnset;    % true|false - trials with framdrop after stimulus motion onset
    started;        % true|false - trials that were started 
    complete;       % true|false - trials that were completed 
    intact;         % true|false - trials that pass sanity check of timing
    
  end     
  
  %% data quality / filtering
  methods (Access = public)

      function keepMat = TrialSelection(d,thdWin,posThd,sacThd,fDropType,maxRT)
        % Detect trials exceed position and saccade threshold
        if nargin < 6
            maxRT = inf;
        end
        AbovePosThd = d.PosDetect(posThd,thdWin);
        AboveSacThd = d.SacDetect(sacThd,thdWin);  
        
        keepMat = d.complete;
        keepMat(:,end+1) = keepMat(:,end) & d.intact;
        keepMat(:,end+1) = keepMat(:,end) & ~d.(fDropType);
        keepMat(:,end+1) = keepMat(:,end) & d.diff.SacDur<maxRT; % response time < 0.5s
        keepMat(:,end+1) = keepMat(:,end) & ~AbovePosThd;
        keepMat(:,end+1) = keepMat(:,end) & ~AboveSacThd;
      end

      function AbovePosThd = PosDetect(d,posThd,win)
        AbovePosThd = zeros(d.numTrials,1);
        t = round(d.eyeData.t);
        tWin = t>=round(win(1)) & t<=round(win(end));
        t = t(tWin);
        x = d.eyeData.x2;
        y = d.eyeData.y2;
        numTrials = size(d.eyeData.x2,1);
        for tr = 1:numTrials 
            if any(abs(x(tr,tWin))>posThd)|any(abs(y(tr,tWin))>posThd)
                AbovePosThd(tr) = 1;
            end
        end 
        AbovePosThd = logical(AbovePosThd);
      end
      
      function AboveSacThd = SacDetect(d,sacThd,win) 
          sp = d.eyeData.sp;
          t = round(d.eyeData.t);
          tWin = t>=round(win(1)) & t<=round(win(end));
          t = t(tWin);
          AboveSacThd = any(sp(:,tWin)>sacThd,2);  
      end          
  end
  
  %% additional properties
  methods (Access = public)
      function [est,par,R2] = LatencyEst(d,eyeData,keep,win,key,conds)
          t = round(d.eyeData.t);
          tWin = t>=round(win(1)) & t<=round(win(end));
          t = t(tWin);
          dv = eyeData.(key);
          numConds = size(conds.value,1);
            
          est = zeros(numConds,numel(t));
          par = zeros(numConds,3);
          R2 = zeros(1,numConds);
          % loop for each condition
          for i = 1:numConds
            ix = keep & conds.Ids == i;
            y = dv(ix,tWin);
            mu = mean(y,1); 
            [est(i,:),par(i,:),R2(i)] = analysis.src.FitPieceWiseLinearInv(mu,[-1 mean(mu) median(t)],t);
          end
      end      
      
        function dta = eyeTrace(d,keep,win,key)
            t = round(d.eyeData.t);
            tWin = t>=round(win(1)) & t<=round(win(end));
            t = t(tWin);

            dv = d.eyeData.(key);
            for kk = 1:size(d.conds.value,1)
                ix = keep & d.conds.Ids == kk;
                y = dv(ix,tWin);
                dta(kk).t = t;
                dta(kk).y = y;
                dta(kk).mu = mean(y,1); 
                dta(kk).SEM = std(y,[],1)/sqrt(size(y,1)); 
            end
        end      
        
        function SEM = LatencyEstCI(d,eyeData,keep,win,key,conds)
          t = round(eyeData.t);
          tWin = t>=round(win(1)) & t<=round(win(end));
          t = t(tWin);
          dv = eyeData.(key);
          numConds = size(conds.value,1);

          btsize = 1000;
          SEM = zeros(numConds,1);
        %   est = zeros(numConds,numel(t));
        %   par = zeros(numConds,3);
        %   R2 = zeros(1,numConds);
          % loop for each condition
          for ii = 1:numConds
            ix = keep & conds.Ids == ii;
            y = dv(ix,tWin);
            bt = bootstrp(btsize,@mean,y);
            dist = zeros(1,btsize);
            for jj = 1:size(bt,1)
                [~,par,~] = analysis.src.FitPieceWiseLinearInv(bt(jj,:),[-1 mean(bt(jj,:)) median(t)],t);
                dist(jj) = par(3);
            end
            SEM(ii) = std(dist);
          end
        end
  end
  
  
end