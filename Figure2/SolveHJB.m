function [u,x,p,how_many_iter] = SolveHJB(x_init,xf,T,dt,xC,r,w,d,W1,W2,sig,tau,tol,max_iter)
%INPUT:
%   x_init   - initial position of the agents
%   xf       - desired final location for agents
%   T        - travel time they are allowed
%   dt       - time discretization parameter
%   xC       - the locations of the centers of the obstacles
%   r        - the radii of the obstacles
%   w        - maximum angular acceleration for Reeds-Shepp vehicles
%   d        - the side length of the square the agents form
%   W1,W2    - weights on travel time vs pattern formation
%   sig,tau  - PDHG hyperparameters
%   tol      - the convergence tolerance
%   max_iter - maximum allowable iteration count
%
%OUTPUT:
%   u - value of the value function at (x_init, 0) given travel time T
%   x - approx optimal paths for agents
%   p - approx optimal costate paths
%   how_many_iter - the final iteration count 

% initial gradient descent rate
gd_rate = 0.1;

% time interval
dim = length(x_init);
t = 0:dt:T;
J = length(t);

%random initialization
x = xf+ 0.1*randn(dim,J);
x(:,round((J+1)/2):J) = x(:,round((J+1)/2):J) - xf  + x_init;
x(:,J) = x_init;
p = 0.1*randn(size(x)); p(:,1) = 0;

% hyperparameters
A = 1;
B = 100;

% split agents into individuals variables
x1 = x(1:2,:);
x2 = x(3:4,:);
x3 = x(5:7,:);
x4 = x(8:10,:);
z1 = x(1:2,:);
z2 = x(3:4,:);
z3 = x(5:7,:);
z4 = x(8:10,:);
p1 = p(1:2,:);
p2 = p(3:4,:);
p3 = p(5:7,:);
p4 = p(8:10,:);

% initialize obstacle functions
O1 = zeros(1,J);
O2 = zeros(1,J);
O3 = zeros(1,J);
O4 = zeros(1,J);
gO1 = zeros(2,J);
gO2 = zeros(2,J);
gO3 = zeros(2,J);
gO4 = zeros(2,J);



%iteration to solve HJB
for k = 1:max_iter
    x1_old = x1; p1_old = p1;
    x2_old = x2; p2_old = p2;
    x3_old = x3; p3_old = p3;
    x4_old = x4; p4_old = p4;
    


    for j = 2:J
        %resolve obstacles
        if isempty(r)
                O1(j) = 1;
                gO1(:,j) = [0;0];
                O2(j) = 1;
                gO2(:,j) = [0;0];
                O3(j) = 1;
                gO3(:,j) = [0;0];
                O4(j) = 1;
                gO4(:,j) = [0;0];
            else
                [D,m] = min((x1(1,j)-xC(1,:)).^2 + (x1(2,j)-xC(2,:)).^2-r.^2);
                O1(j) = 1/2 + (1/2)*tanh(B*D);
                gO1(:,j) = 2*B*sech(B*D)^2*(x1(:,j)-xC(:,m));
                [D,m] = min((x2(1,j)-xC(1,:)).^2 + (x2(2,j)-xC(2,:)).^2-r.^2);
                O2(j) = 1/2 + (1/2)*tanh(B*D);
                gO2(:,j) = 2*B*sech(B*D)^2*(x2(:,j)-xC(:,m));
                [D,m] = min((x3(1,j)-xC(1,:)).^2 + (x3(2,j)-xC(2,:)).^2-r.^2);
                O3(j) = 1/2 + (1/2)*tanh(B*D);
                gO3(:,j) = 2*B*sech(B*D)^2*(x3(1:2,j)-xC(:,m));
                [D,m] = min((x4(1,j)-xC(1,:)).^2 + (x4(2,j)-xC(2,:)).^2-r.^2);
                O4(j) = 1/2 + (1/2)*tanh(B*D);
                gO4(:,j) = 2*B*sech(B*D)^2*(x4(1:2,j)-xC(:,m));
        end
        
        % update costat variables
        b = p1(:,j) + sig*(z1(:,j)-z1(:,j-1));
        p1(:,j) = max(0,1-sig*dt*O1(j)/norm(b,2))*b;

        b = p2(:,j) + sig*(z2(:,j)-z2(:,j-1));
        p2(:,j) = max(0,1-sig*dt*O2(j)/norm(b,2))*b;

        g = [cos(x3(3,j));sin(x3(3,j))];
        b = p3(:,j) + sig*(z3(:,j)-z3(:,j-1));
        gtb = g'*b(1:2);
        p3(1:2,j) = (max(0,1-sig*dt*O3(j)/abs(gtb))-1)*gtb*g+b(1:2);
        p3(3,j) = max(0,1-w*dt*sig*O3(j)/abs(b(3)))*b(3);

        g = [cos(x4(3,j));sin(x4(3,j))];
        b = p4(:,j) + sig*(z4(:,j)-z4(:,j-1));
        gtb = g'*b(1:2);
        p4(1:2,j) = (max(0,1-sig*dt*O4(j)/abs(gtb))-1)*gtb*g+b(1:2);
        p4(3,j) = max(0,1-w*dt*sig*O4(j)/abs(b(3)))*b(3);
    end

    % update state variables
    x1(:,1) = x1(:,1)+tau*p1(:,2);
    x2(:,1) = x2(:,1)+tau*p2(:,2);
    x3(:,1) = x3(:,1)+tau*p3(:,2);
    x4(:,1) = x4(:,1)+tau*p4(:,2);
    for j = 2:J-1

        % with these variable definitions, the penalty function for pattern
        % formation is P = L1^2 + L2^2 + norm(Mp)^2 + Ip^2
        % This will zero iff the agents for a square with agent 1 and 2 in
        % opposite corners and agents 3 and 4 in opposite corners
        L1 = norm(x1(:,j)-x2(:,j),2)^2 - 2*d^2;
        L2 = norm(x3(1:2,j)-x4(1:2,j),2)^2 - 2*d^2;
        Mp = x1(:,j) + x2(:,j) - (x3(1:2,j)+x4(1:2,j));
        Ip = (x1(:,j)-x2(:,j))'*(x3(1:2,j)-x4(1:2,j));
        
        % gradients of the penalty functions P and C with respect to each
        % agent
        gP1 =  4*L1*(x1(:,j)-x2(:,j)) + 2*Ip*(x3(1:2,j)-x4(1:2,j)) + 2*Mp;
        gP2 =  4*L1*(x2(:,j)-x1(:,j)) - 2*Ip*(x3(1:2,j)-x4(1:2,j)) + 2*Mp;
        gP3 =  4*L2*(x3(1:2,j)-x4(1:2,j)) + 2*Ip*(x1(:,j)-x2(:,j)) - 2*Mp; 
        gP4 =  4*L2*(x4(1:2,j)-x3(1:2,j)) - 2*Ip*(x1(:,j)-x2(:,j)) - 2*Mp; 
        gC1 = 2*A*exp(-A*norm(x1(:,j)-xf(1:2),2)^2)*(x1(:,j)-xf(1:2));
        gC2 = 2*A*exp(-A*norm(x2(:,j)-xf(3:4),2)^2)*(x2(:,j)-xf(3:4));
        gC3 = 2*A*exp(-A*norm(x3(:,j)-xf(5:7),2)^2)*(x3(:,j)-xf(5:7));
        gC4 = 2*A*exp(-A*norm(x4(:,j)-xf(8:10),2)^2)*(x4(:,j)-xf(8:10));
        
        % update first two agents
        XI = x1(:,j) - tau*(p1(:,j)-p1(:,j+1));
        x1(:,j) = x1(:,j) - gd_rate*(-dt*tau*(norm(p1(:,j),2)*gO1(:,j) - W1*gC1 - W2*gP1) + (x1(:,j)-XI));
        XI = x2(:,j) - tau*(p2(:,j)-p2(:,j+1));
        x2(:,j) = x2(:,j) - gd_rate*(-dt*tau*(norm(p2(:,j),2)*gO2(:,j) - W1*gC2 - W2*gP2) + (x2(:,j)-XI));
        
        % update agents 3 and 4
        XI = x3(:,j) - tau*(p3(:,j)-p3(:,j+1));
        g = [cos(x3(3,j));sin(x3(3,j))];
        gprime = [-sin(x3(3,j));cos(x3(3,j))];
        gtp = g'*p3(1:2,j); gptp = gprime'*p3(1:2,j);
        x3(1:2,j) = x3(1:2,j) - gd_rate*(-dt*tau*((abs(gtp)+w*abs(p3(3,j)))*gO3(:,j) - W1*gC3(1:2) - W2*gP3) + (x3(1:2,j)-XI(1:2)));
        x3(3,j) = x3(3,j) - gd_rate*(-dt*tau*(O3(j)*sign(gtp)*gptp - W1*gC3(3)) + x3(3,j) - XI(3));
        XI = x4(:,j) - tau*(p4(:,j)-p4(:,j+1));
        g = [cos(x4(3,j));sin(x4(3,j))];
        gprime = [-sin(x4(3,j));cos(x4(3,j))];
        gtp = g'*p4(1:2,j); gptp = gprime'*p4(1:2,j);
        x4(1:2,j) = x4(1:2,j) - gd_rate*(-dt*tau*((abs(gtp)+w*abs(p4(3,j)))*gO4(:,j) - W1*gC4(1:2) - W2*gP4) + (x4(1:2,j)-XI(1:2)));
        x4(3,j) = x4(3,j) - gd_rate*(-dt*tau*(O4(j)*sign(gtp)*gptp-W1*gC4(3)) + x4(3,j) - XI(3));
    end
    
    %update relaxation variables
    z1 =  x1 + (x1-x1_old);
    z2 =  x2 + (x2-x2_old);
    z3 =  x3 + (x3-x3_old);
    z4 =  x4 + (x4-x4_old);

    % check convergence criterion
    change1 = max(max(max(abs(x1-x1_old))), max(max(p1-p1_old)));
    change2 = max(max(max(abs(x2-x2_old))), max(max(p2-p2_old)));
    change3 = max(max(max(abs(x3-x3_old))), max(max(p3-p3_old)));
    change4 = max(max(max(abs(x4-x4_old))), max(max(p4-p4_old)));
    change = max([change1,change2,change3,change4]);
    if change < tol
        break
    end
    
    % print progress
    if mod(k,1000) == 0
        fprintf('Completed iteration %i. change = %.4e\n',k,change);
        % decrement gradient descent rate every 1000 iterations (after 5000)
        if  k>=5000
            gd_rate = gd_rate*0.5;
        end
    end

end

% compute value function
u = 0;
for j = 2:J
    L1 = norm(x1(:,j)-x2(:,j),2)^2 - 2*d^2;
        L2 = norm(x3(1:2,j)-x4(1:2,j),2)^2 - 2*d^2;
        Mp = norm(x1(:,j) + x2(:,j) - (x3(1:2,j)+x4(1:2,j)));
        Ip = (x1(:,j)-x2(:,j))'*(x3(1:2,j)-x4(1:2,j));
    u = u + p1(:,j)'*(x1(:,j)-x1(:,j-1)) + p2(:,j)'*(x2(:,j)-x2(:,j-1))+ p3(:,j)'*(x3(:,j)-x3(:,j-1)) + p4(:,j)'*(x4(:,j)-x4(:,j-1))...
        + dt*W1*(4 - exp(-A*norm(x1(:,1)-xf(1:2))^2)-exp(-A*norm(x2(:,1)-xf(3:4))^2)...
                   - exp(-A*norm(x3(:,1)-xf(5:7))^2)-exp(-A*norm(x4(:,1)-xf(8:10))^2) )...
        + dt*W2*(L1^2 + L2^2 +Mp^2 + Ip^2)...
        + dt*(O1(j)*norm(p1(:,j))+O2(j)*norm(p2(:,j))+...
              O3(j)*(abs(p3(1,j)*cos(x3(3,j)) + p3(2,j)*sin(x3(3,j))) + w*abs(p3(3,j)))+...
              O4(j)*(abs(p4(1,j)*cos(x4(3,j)) + p4(2,j)*sin(x4(3,j))) + w*abs(p4(3,j))));
end


x(1:2,:)  = x1;
x(3:4,:)  = x2;
x(5:7,:)  = x3;
x(8:10,:) = x4;
p(1:2,:)  = p1;
p(3:4,:)  = p2;
p(5:7,:)  = p3;
p(8:10,:) = p4; 
how_many_iter = k;

end

