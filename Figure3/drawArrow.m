function [m] = drawArrow(x0,y0,xf,yf,scale,color,thickness)
%
% A function which draws an arrow starting at (x0,y0) and pointing to
% (xf,yf). 
%   'scale' = size of the arrowhead 
%   'color' = color of the lines
%   'thickness' = thickness of lines
%

if nargin < 5
    scale = 1;
    color = 'k';
    thickness = 1.5;
end
rotate = @(s) [cos(s) -sin(s); sin(s) cos(s)];
hold on;
dist = sqrt((xf-x0)^2+(yf-y0)^2);
tri = scale*(1/100)*[-1 2 -1 -1;
        2 0 -2 2];
t = 0:0.001:(1-scale*(1/50)/dist);
if isempty(t)
   t = 0:0.001:1; 
end
lx = (1-t)*x0 + t*xf;
ly = (1-t)*y0 + t*yf;

if xf >= x0
   ang = atan((yf-y0)/(xf-x0)); 
else
   ang = atan((yf-y0)/(xf-x0)) + pi;
end
tri = rotate(ang)*tri;
   
plot(lx,ly,'-','color',color,'linewidth',thickness);
F = fill(tri(1,:)+lx(end),tri(2,:)+ly(end),'k'); F.FaceColor = color; F.EdgeColor = 'none';

m=0;
end

