function [RespFit, ParamFinal, R2] = FitPieceWiseLinearInv(y,ParamInit,x);
%
%[RespFit, ParamFinal, R2] = FitPieceWiseLinear(y,ParamInit,x);
%finds least-square fits using a two-part piece-wise linear function with a
%plateau in the rightmost-section
%i.e.   y = ax+b if x>=c
%       y = ac+b if x<c
%if y is 2D, it does fits to each ROW in turn
%
% e.g.
% y = [10 9 8 7 6 5.1 5.05 5 5 5];
% [R,P,R2] = FitPieceWiseLinear(y)

if nargin < 3, x = 1:size(y,2); end
if nargin < 2, ParamInit = [-1 mean(y) median(x)]; end

PWLfun = inline('(x(1)*tt+x(2)).*(tt>=x(3)) + (x(1)*x(3)+x(2)).*(tt<x(3))', 'x','tt');
%x = [slope, y-offset, change-point];
LB = [-inf -inf -inf];
UB = [inf inf inf];
xin = x;

for a = 1:size(y,1), %loop through rows
    %need to ignore NaN values for fitting, but still output appropriate value
    %at specified x value
    yNaN = isnan(y(a,:));
    Thisy = y(a,~yNaN);
    Thisx = x(~yNaN);

    LSQoptions = optimset('Display', 'off'); %, 'TolFun',0.001); %'MaxIter', 400); %
    [Params,Resnorm, Resid, exitflag, output, lambda, Jacob] = ...
        lsqcurvefit(PWLfun, ParamInit, Thisx, Thisy,LB,UB,LSQoptions);
    ParamFinal(a,:) = Params;
    RespFit(a,:) = PWLfun(Params, xin);
    
    R = corrcoef(PWLfun(Params, Thisx), Thisy);
    R2(a) = R(2)^2;
    
end