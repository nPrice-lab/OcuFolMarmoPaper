function [RespFit, Params, R2,p,F] = gauss2DSFTF(sf,tf,y,type,x0)
    %% doc string here
    %
    %   sf: ntrials X 1 vector of spatial freq conditions
    %   tf: ntrials X 1 vector of temporal freq conditions
    %   y: ntrials X 1 vector of measurements (e.g. eye speeds/spikes)
    %   x0: initial values of 7 parameters


%     xdata = [sf tf];
   
    % parameters (see equations 1 and 2 in Miura et al 2014
    % x(1) = A
    % x(2) = b
    % x(3) = sf0
    % x(4) = tf0
    % x(5) = sd_sf
    % x(6) = sd_tf
    % x(7) = Q
    
%     sf1 = logspace(-2,4,7);
%     tf1 = logspace(-1,5,6); 

%% stimulation
%     sf1 = [0.0400    0.0800    0.1600    0.3100    0.6200    1.2400    2.4800];
%     tf1 = [1.5600    3.1300    6.1300   12.2500   18.7500   25.0000];
% 
%     [sfm,tfm] = meshgrid(sf1,tf1);
% 
%     x(1) = 1;
%     x(2) = 0;
%     x(3) = 0.37;
%     x(4) = 19.38;
%     x(5) = 3; % SD
%     x(6) = 10;
%     x(7) = 1;    
% %     out = Gauss2dfitFun(x,[sfm(:) tfm(:)]);
% %     out = reshape(out,size(sfm));
%     xdata(:,:,1) = sfm;
%     xdata(:,:,2) = tfm;
%     out = Gauss2dfitFun(x,xdata);
%     
%     figure
%     surf(sfm,tfm,out)
% %     plot3(sf(:),tf(:),out,'k.'), 
%     set(gca,'xscale','log','yscale','log')
%     xlabel('sf')
%     ylabel('tf')
% 
% keyboard
% 
%     err = normrnd(0,0.05,size(out));
%     y = out+err;
%     figure
%     surf(sfm,tfm,y)
%     set(gca,'xscale','log','yscale','log')
% keyboard

%     sf = sfm(:);
%     tf = tfm(:);
%     xdata = [sf tf];
%     y = y(:);

%% calculate condition means

sf1 = unique(sf);
tf1 = unique(tf);
[tfm,sfm] = meshgrid(tf1,sf1);
xdata(:,:,1) = sfm;
xdata(:,:,2) = tfm;
ydata = prepData(sf,tf,y,type);
%% model fitting
%     if nargin < 4, x0 =  [max(y) mean(y) mean(sf) mean(tf) std(y) std(y) 1];  end
    
%     LB = [-inf -inf min(sf) min(tf) 0 0 0];
%     UB = [inf inf max(sf) max(tf) inf inf 1];

%       x0 = x; % using x works
%       x0 = [x(1) x(2) mean(sf1) mean(tf1) x(5) x(6) x(7)]; % works  
%       x0 = [max(y(:)) mean(y(:)) mean(sf1) mean(tf1) x(5) x(6) x(7)]; % works  
%       x0 = [max(y(:)) mean(y(:)) mean(sf1) mean(tf1) x(5) x(6) 0.5]; % works  
%       x0 = [max(y(:)) mean(y(:)) mean(sf1) mean(tf1) std(y(:)) std(y(:)) 0.5]; % not work 
%       x0 = [max(y(:)) mean(y(:)) mean(sf1) mean(tf1) std(sf1) std(tf1) 0.5]; % not work
      if nargin < 5, x0 = [max(y(:)) mean(y(:)) mean(sf) mean(tf) 50 50 1]; end% works 

      LB = [-inf -inf min(sf) min(tf) 0 0 0];
      UB = [inf inf max(sf) max(tf) inf inf 1];

    LSQoptions = optimset('Display', 'off'); %, 'TolFun',0.001); %'MaxIter', 400); %
    [Params,Resnorm, Resid, exitflag, output, lambda, Jacob] = ...
        lsqcurvefit(@Gauss2dfitFun, x0, xdata, ydata, LB, UB, LSQoptions);

    RespFit = Gauss2dfitFun(Params, xdata);
    R = corrcoef(RespFit, ydata);
    R2 = R(2)^2;
    df = numel(ydata)-7;
    
    % set Q = 0, don't estimate it
    x0 = x0(1:end-1);
    LB = LB(1:end-1);
    UB = UB(1:end-1);
    [Params0,Resnorm, Resid, exitflag, output, lambda, Jacob] = ...
        lsqcurvefit(@Gauss2dfitFun0, x0, xdata, ydata, LB, UB, LSQoptions);
    RespFit0 = Gauss2dfitFun0(Params0, xdata);
    
    R0 = corrcoef(RespFit0, ydata);
    R20 = R0(2)^2;
    df0 = numel(ydata)-6;
    
    [p, F] = ftest(ydata(:), RespFit0(:),RespFit(:), df0, df);
    
    
%     figure
% %     RespFit = reshape(RespFit,size(out));
%     surf(sfm,tfm,RespFit);
%     title(['free Q, R2=',num2str(R2)])
%     set(gca,'xscale','log','yscale','log')
% 
%         figure
% %     RespFit = reshape(RespFit,size(out));
%     surf(sfm,tfm,RespFit0);
%     title(['fix Q=0, R2=',num2str(R20)])
%     set(gca,'xscale','log','yscale','log')
% 
%         figure
% %     RespFit = reshape(RespFit,size(out));
%     surf(sfm,tfm,ydata);
%     set(gca,'xscale','log','yscale','log')
%     keyboard
%     
%     figure
% %     y = reshape(y,size(out));
%     surf(sfm,tfm,ydata);
%     set(gca,'xscale','log','yscale','log')
% keyboard
end

function y = Gauss2dfitFun(x,xdata)
% parameters (see equations 1 and 2 in Miura et al 2014
A = x(1);
b = x(2);
sf0 = x(3);
tf0 = x(4);
SDsf = x(5);
SDtf = x(6);
Q = x(7);

sf = xdata(:,:,1);
tf = xdata(:,:,2);

% sf = xdata(:,1);
% tf = xdata(:,2);
    
tfs = 2.^(Q*(log2(sf)-log2(sf0))+log2(tf0));
y = A*exp(-(log2(sf)-log2(sf0)).^2/SDsf^2).*exp(-(tf-tfs).^2/SDtf^2)+b;
end

function y = Gauss2dfitFun0(x,xdata)
% parameters (see equations 1 and 2 in Miura et al 2014, fix Q = 0
A = x(1);
b = x(2);
sf0 = x(3);
tf0 = x(4);
SDsf = x(5);
SDtf = x(6);
Q = 0;

sf = xdata(:,:,1);
tf = xdata(:,:,2);

% sf = xdata(:,1);
% tf = xdata(:,2);
    
tfs = 2.^(Q*(log2(sf)-log2(sf0))+log2(tf0));
y = A*exp(-(log2(sf)-log2(sf0)).^2/SDsf^2).*exp(-(tf-tfs).^2/SDtf^2)+b;
end

function map = prepData(sf,tf,y,type)
    
    sfCond = unique(sf);
    tfCond = unique(tf);

    map = zeros(numel(sfCond),numel(tfCond));
    for t = 1:numel(tfCond)
        for s = 1:numel(sfCond)
            ix = sf == sfCond(s) & tf == tfCond(t);
            if strcmp(type,'gain')
                div = tfCond(t)/sfCond(s);
            else
                div = 1;
            end
            map(s,t) = mean(y(ix))/div;
        end    
    end
end