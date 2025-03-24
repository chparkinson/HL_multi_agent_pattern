clear; 

%
% In this example, we have 3 agents exhibiting isotropic motion 
% They will travel to their final points while trying to maintain
% formation as an equilateral triangle. 
%

%rng seed for exact results from paper when run on first author's desktop
% rng(12345);

% choose side length for the triangle
d = 1/2;

% choose initial points and final points
x_init = [2;-2;0;-2;-2;-2]; 
xf = [0; 1; 0; 1; 0; 1] + (d/sqrt(3))*[cos(pi/2);sin(pi/2);cos(pi/2 + 2*pi/3);sin(pi/2 + 2*pi/3);cos(pi/2+4*pi/3);sin(pi/2 + 4*pi/3)];

% obstacles (xC = centers, r = radii)
xC = [0 0;-1/2 1/2];
r = [1/4,1/4];
%xC = []; r = []; % for no obstacles leave these empty 


% choose weights which balance penalty for breaking formation or spending
% time away from final point
W1 = 0.5; % importance of minimizing travel time 
W2 = 3; % importance of keeping formation


% choose horizon time and dt
T = 6.1; 
dt = 0.1;
sig = 1; tau = 0.25/sig; tol = 5e-4; max_iter = 50000;
timer = tic;
[u,x,p,how_many_iter] =  SolveHJB(x_init,xf,T,dt,xC,r,d,W1,W2,sig,tau,tol,max_iter);
fprintf('Algorithm converged to tolerance %.4e after %i iterations and %.2f seconds\n',tol,how_many_iter,toc(timer));

%% plot results
num_agents = length(x_init)/2;
F = figure(11);
J = size(x,2);
t = 0:dt:T;
for i = 1:num_agents
   COLOR{i} = 0.3 + 0.7*rand(1,3); COLOR{i}(randi(3))= 0; 
end

v=@(x,y) 1+0.25*sin(2*x).*sin(2*y);
[X,Y] = ndgrid(-3:0.03:3,-3:0.03:3);
V = v(X,Y);

% Predefined colors for 3-agent case
COLOR(1:3) = {[0.9 0.1 0.1], [0.6 0.1 0.9], [0.4 0.4 0.4]};

s = 0:2*pi/100:2*pi;
count = 1;
% for j = [32,28,22,17,9]
for j = J:-1/(10*dt):2
    clf; hold on;
    contourf(X,Y,V,50,'edgecolor','none');
    xCt = [cos(2*t(j)) -sin(2*t(j)); sin(2*t(j)) cos(2*t(j))]*xC;
    for i = 1:length(r)
        fill(xCt(1,i)+r(i)*cos(s),xCt(2,i)+r(i)*sin(s),'k','EdgeColor','none');
        drawArrow(xCt(1,i)-0.9*r(i)*cos(2*t(j)+mod(i,2)*pi),xCt(2,i)-0.9*r(i)*sin(2*t(j)+mod(i,2)*pi),xCt(1,i)+r(i)*cos(2*t(j)+mod(i,2)*pi),xCt(2,i)+r(i)*sin(2*t(j)+mod(i,2)*pi),5,'g',4);
    end
    for i = 1:num_agents
        ind = (2*(i-1)+1):2*i;
        plot(x(ind(1),J),x(ind(2),J),'.','markersize',25,'color',COLOR{i});
        plot(x(ind(1),j:J),x(ind(2),j:J),'-','linewidth',2.5,'color',COLOR{i});
        plot(x(ind(1),j),x(ind(2),j),'.','markersize',25,'color',COLOR{i});
        plot(xf(ind(1)),xf(ind(2)),'.','markersize',25,'color',COLOR{i});
    end
    % TITLE=title(sprintf('$t=%.2f$',t(J-j+1)));TITLE.FontSize=20;TITLE.Interpreter = 'latex';
    TE = text(-1.6,2,sprintf('$t=%.2f$',t(J-j+1))); TE.FontSize = 20; TE.Interpreter = 'latex';
    TE.FontWeight = 'bold';
    axis([-2.25 2.25 -2.25 2.25]);
    axis square;
    axis off;
    pause(0.1);

    %%% uncomment this to print pictures
    %%%   WARNING: if you are doing the animation it will print every frame
    %
    % print(sprintf('ex3n%i',count),'-dpng');
    % count = count + 1;
end