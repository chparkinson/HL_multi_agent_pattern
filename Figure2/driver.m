clear;

%
% In this example, we have 4 agents. Agents 1 and 2 exhibit isotropic
% motion while agents 3 and 4 are Reeds-Shepp cars. The agents will attempt
% to navigate to their final locations while maintaining formation in a
% square, where the cars are opposite corners.
%
%

% This is the rng seed that led to the figures from the paper 
% on the first author's desktop computer
% rng(10101);

% choose side length of square
d = 1/2;


% choose initial points and final points
%       -  initial and final points for the cars consist
%          of (x,y,theta) where theta is orientation
x_init = [-1; -2; 1; -2; 0; -2; 0; 0; -1.5; pi];
xf = [-0.25; 2.25; 0.25; 1.75; -0.25; 1.75; pi/2; 0.25; 2.25; 3*pi/2];

% choose turning radius for cars (min turning radius = 1/w)
w = 2;

% obstacles (xC = centers, r = radii) - leave empty for no obs: xC = [], r = []
% the ones from the paper are saved in obs.mat
load obs.mat;

% UNCOMMENT THIS TO RANDOMLY GENERATE DIFFERENT OBSTACLES
% USING THE INCLUDED OBSTACLE GENERATION FUNCTION
%[xC,r] = generateDisjointCircles(-1,1,-0.375,0.375,0.05,0.15,0.05,12); 



% choose weights which balance penalty for breaking formation or spending
% time away from final point
W1 = 1; % importance of minimizing travel time
W2 = 1; % importance of keeping formation


% choose horizon time and dt and hyperparameters
T = 7.1;
dt = 0.1;
sig = 1; tau = 0.25/sig; tol = 5e-4; max_iter = 50000;
timer = tic;
[u,x,p,how_many_iter] =  SolveHJB(x_init,xf,T,dt,xC,r,w,d,W1,W2,sig,tau,tol,max_iter);
fprintf('Algorithm converged to tolerance %.4e after %i iterations and %.2f seconds\n',tol,how_many_iter,toc(timer));

%% plot results
F = figure(11);
J = size(x,2);
t = 0:dt:T;
s = 0:2*pi/100:2*pi;

% Predefined colors for 4-agent case
COLOR(1:4) = {[0.7597, 0, 0.8282], [0.5995, 0.4870, 0], [0,0.5848, 0.9568], [0.9027,0.1123,0.1603]};

count = 1;
% for snapshots from paper
% for j = [J, J-12, J-25, J-37,J-50]
for j = J:-1/(10*dt):22
    clf; hold on;
    for i = 1:length(r)
        fill(xC(1,i)+r(i)*cos(s),xC(2,i)+r(i)*sin(s),'k','EdgeColor','none');
    end
    %agent 1
    i=1;
    plot(x(1,J),x(2,J),'.','markersize',25,'color',COLOR{i});
    plot(x(1,j:J),x(2,j:J),'-','linewidth',2.5,'color',COLOR{i});
    plot(xf(1),xf(2),'.','markersize',25,'color',COLOR{i});
    plot(x(1,j),x(2,j),'.','markersize',25,'color',COLOR{i});

    %agent 2
    i=2;
    plot(x(3,J),x(4,J),'.','markersize',25,'color',COLOR{i});
    plot(x(3,j:J),x(4,j:J),'-','linewidth',2.5,'color',COLOR{i});
    plot(xf(3),xf(4),'.','markersize',25,'color',COLOR{i});
    plot(x(3,j),x(4,j),'.','markersize',25,'color',COLOR{i});

    %agent 3
    i=3;
    plot(x(5,J),x(6,J),'.','markersize',25,'color',COLOR{i});
    plot(x(5,j:J),x(6,j:J),'-','linewidth',2.5,'color',COLOR{i});
    plot(xf(5),xf(6),'.','markersize',25,'color',COLOR{i});
    drawCar(x(5,j),x(6,j),0.2,0.1,x(7,j),COLOR{i},1);

    %agent 4
    i=4;
    plot(x(8,J),x(9,J),'.','markersize',25,'color',COLOR{i});
    plot(x(8,j:J),x(9,j:J),'-','linewidth',2.5,'color',COLOR{i});
    plot(xf(8),xf(9),'.','markersize',25,'color',COLOR{i});
    drawCar(x(8,j),x(9,j),0.2,0.1,x(10,j),COLOR{i},1);

    %TITLE=title(sprintf('$t=%.2f$',t(J-j+1)));TITLE.FontSize=15;TITLE.Interpreter = 'latex';
    TE = text(-2.25,2.25,sprintf('$t=%.2f$',t(J-j+1))); TE.FontSize = 20; TE.Interpreter = 'latex';
    axis([-2.5 2.5 -2.5 2.5]);
    axis square;
    axis off;
    pause(0.1);

    % Uncomment to print images to .png
    % WARNING: if you do this for the animation, it will print every single
    %          frame
    % print(sprintf('ex2n%i',count),'-dpng');
    % count = count+1;
end