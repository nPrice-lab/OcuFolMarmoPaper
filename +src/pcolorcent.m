function H = pcolorcent(varargin)
% PCOLORCENT(X,Y,C) - pseudocolor or "checkerboard" plot with colors in matrix
% C plotted centered on the coordinates in X,Y
% This is trying to replicate the behavior of image/imagesc, which plot
% each value in C into a single pixel, while allowing the variable X and Y
% positions as in pcolor.
% My main gripe with pcolor is that patch centers are not aligned on the
% X,Y coordinates 
%
% This also flips the x/y conventions, so that C is (x,y), rather than
% (row,column)

% Functionality is achieved just by recalculating X/Y values and passing
% them directly to pcolor
%
% See also pcolor, imagesc

if nargin==1,
    c = varargin{1};
    x = 1:size(c,1);
    y = 1:size(c,2);
elseif nargin>=3,
    [x,y,c] = deal(varargin{1:3});
end

% modify vectors
% Should implement a check that x,y are vectors, not matrices...

if length(x)==1,
    x2 = x*[0.75 1.25];
else
    dx = diff(x)/2;
    x2 = [x(1)-dx(1) x+[dx dx(end)]];
end
if length(y)==1,
    y2 = y*[0.75 1.25];
else
    dy = diff(y)/2;
    y2 = [y(1)-dy(1) y+[dy dy(end)]];
end

c2 = [c' nan(length(y),1); nan(1,length(x)+1)];

if nargin>=4
    H = pcolor(varargin{4},x2,y2,c2);
else
    H = pcolor(x2,y2,c2);    
end
