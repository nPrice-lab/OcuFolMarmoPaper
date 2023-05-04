function [eyeData, params] = percPurs_Diff(eyeData, varargin)
%
% given input eyeData positions (.x/.y) this function calculates
% - smoothed eye position (.x2/.y2). Units = deg
% - velocity (.dx/.dy). Units = deg/s
% - speed (.sp)
% - acceleration (.acc). This is a scalar. Units deg/s/s
%
% Smoothing and differentiation are performed using Savitzky-Golay filter
% with parameters specified as named pairs
% e.g.
% ord = 5; fl=101;
% [eyeData, params] = percPurs_Diff(eyeData,'order',ord,'framelen',fl);
%
% The filter frequency response can be visualised with freqz(SGg(:,1))
% ord=5;fl=101; gives a 3db cutoff freq at ~34 Hz  
%
% This function replaces part of percPurs_smooth
% See also percPurs_FileList, percPurs_Merge
% NP 11-03-2021

%   file = 'c:\data\percPurs\NP.percPurs_Discrim250.mat';
%   load(file,'eyeData')

p = inputParser;
addRequired(p,'eyeData',@isstruct);
addParameter(p,'order',5,@isnumeric); % Differentiating filter - Order of polynomial fit
addParameter(p,'framelen',51,@(x) mod(x,2)==1); % Differentiating filter - Half window length
p.KeepUnmatched = true;
parse(p,eyeData,varargin{:})
params = p.Results;

% define Savitzky-Golay differentiating filter
[~,SGg] = sgolay(params.order, params.framelen);   % Calculate S-G coefficients

eyeData.dt = eyeData.t(2)-eyeData.t(1); % step size in ms. 

if iscell(eyeData.x)
    eyeData.x2 = cellfun(@(x) conv(x,SGg(:,1), 'same'),eyeData.x, 'UniformOutput',false); % position
    eyeData.y2 = cellfun(@(x) conv(x,SGg(:,1), 'same'),eyeData.y, 'UniformOutput',false);
    eyeData.dx = cellfun(@(x) -conv(x, SGg(:,2)', 'same')*(1000/eyeData.dt), eyeData.x, 'UniformOutput', false);
    eyeData.dy = cellfun(@(x) -conv(x, SGg(:,2)', 'same')*(1000/eyeData.dt), eyeData.y, 'UniformOutput', false);
    
    eyeData.ddx = cellfun(@(x) conv(x, 2*SGg(:,3), 'same')*(1000/eyeData.dt)^2, eyeData.x, 'UniformOutput', false);
    eyeData.ddy = cellfun(@(x) conv(x, 2*SGg(:,3), 'same')*(1000/eyeData.dt)^2, eyeData.y, 'UniformOutput', false);
    
    % absolute value of speed/accel
    eyeData.sp = cellfun(@(x1,x2) hypot(x1, x2), eyeData.dx, eyeData.dy, 'UniformOutput', false);
    eyeData.acc = cellfun(@(x1,x2) hypot(x1, x2), eyeData.ddx, eyeData.ddy , 'UniformOutput', false);
    
else
    eyeData.x2 = conv2(eyeData.x, SGg(:,1)', 'same'); % position
    eyeData.y2 = conv2(eyeData.y, SGg(:,1)', 'same');
    eyeData.dx = -conv2(eyeData.x, SGg(:,2)', 'same')*(1000/eyeData.dt); % velocity
    eyeData.dy = -conv2(eyeData.y, SGg(:,2)', 'same')*(1000/eyeData.dt);
    
    eyeData.ddx = conv2(eyeData.x, 2*SGg(:,3)', 'same')*(1000/eyeData.dt)^2; % acceleration
    eyeData.ddy = conv2(eyeData.y, 2*SGg(:,3)', 'same')*(1000/eyeData.dt)^2;
    
    % absolute value of speed/accel
    eyeData.sp = hypot(eyeData.dx, eyeData.dy);
    eyeData.acc = hypot(eyeData.ddx, eyeData.ddy);
end
% if 0
% figure
% plot(eyeData.t, eyeData.x(1,:),'b')
% hold on, 
% plot(eyeData.t, eyeData.x2(1,:),'r','linewidth',3)
% end