function [ h ] = drawCar( x,y,d,R,s,col,alp)
% inputs:
%   (x,y) - location of center of mass
%   d - dist from rear axle to center of mass
%   R - half width of rear axle
%   s - orientation of car (s = 0 is facing due right)
%   col - color of car (RGB triplet)
%   alp - face alpha of car (alp = 1 means not transparent, alp = 0 means fully transparent)


if nargin < 6
   col = 'r';
   alp = 1;
elseif nargin < 7
    alp = 1;
end

load car.mat;

rotate = @(theta) [cos(theta), -sin(theta); sin(theta),  cos(theta)];

% xH1 = d*[xCar{4}(40:150),xCar{4}(150),xCar{4}(40)];
% yH1 = R*[yCar{4}(40:150),0.8*yCar{4}(150),yCar{4}(40)];
xH1 = d*[xCar{4}(40:150),xCar{4}(150),xCar{4}(40)];
yH1 = R*[yCar{4}(40:150),yCar{4}(40),yCar{4}(40)];
H1 = [xH1;yH1];
H1 = rotate(s)*H1;

% xH2 = d*[xCar{3}(173:283),xCar{3}(173),xCar{3}(173)];
% yH2 = R*[yCar{3}(173:283),0.8*yCar{3}(173),yCar{3}(173)];
xH2 = d*[xCar{3}(173:283),xCar{3}(173),xCar{3}(173)];
yH2 = R*[yCar{3}(173:283),yCar{3}(283),yCar{3}(173)];
H2 = [xH2;yH2];
H2 = rotate(s)*H2;

xC = d*[xCar{1},xCar{2},xCar{3},xCar{4},xCar{1}(1)];
yC = R*[yCar{1},yCar{2},yCar{3},yCar{4},yCar{1}(1)];
C = [xC;yC];
C = rotate(s)*C;
k=1;
for i = windows
   xW{i}(end+1) = xW{i}(1); yW{i}(end+1) = yW{i}(1);
   xW{i} = d*xW{i}; yW{i} = R*yW{i};
   W{k} = [xW{i};yW{i}];
   W{k} = rotate(s)*W{k};
   k=k+1;
end

gold = [218,165,32]/255;
hold on; 
F = fill(C(1,:)+x,C(2,:)+y,col,'FaceAlpha',alp);
fill(H1(1,:)+x,H1(2,:)+y,gold,'facealpha',alp,'edgecolor','none');
fill(H2(1,:)+x,H2(2,:)+y,gold,'facealpha',alp,'edgecolor','none');
for k = 1:length(W)
    fill(W{k}(1,:)+x,W{k}(2,:)+y,'k','facealpha',alp);
end


h = 0;
end

