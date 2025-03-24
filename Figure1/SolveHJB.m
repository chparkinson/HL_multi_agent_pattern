function [u,x,p,how_many_iter] = SolveHJB(x_init,xf,T,dt,xC,r,d,W1,W2,sig,tau,tol,max_iter)
%INPUT:
%   x_init   - initial position of the agents
%   xf       - desired final location for agents
%   T        - travel time they are allowed
%   dt       - time discretization parameter
%   xC       - the locations of the centers of the obstacles
%   r        - the radii of the obstacles
%   d        - the side length of the equilateral triangle the agents form
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

% gradient descent rate
gd_rate = 0.1;

% time interval
dim = length(x_init);
t = 0:dt:T;
J = length(t);

% random initialization
x = xf + 0.1*randn(dim,J);
x(:,round((J+1)/2):J) = x(:,round((J+1)/2):J) - xf + x_init;
x(:,J) = x_init;
p = 0.1*randn(size(x)); p(:,1) = 0;

% velocity function and gradient
v = @(x) 1;
gv = @(x) zeros(2,1);

%hyperparameters and splitting up into agents
A = 1;
B = 100;
x1 = x(1:2,:);
x2 = x(3:4,:);
x3 = x(5:6,:);
z1 = x(1:2,:);
z2 = x(3:4,:);
z3 = x(5:6,:);
p1 = p(1:2,:);
p2 = p(3:4,:);
p3 = p(5:6,:);
O1 = zeros(1,J);
O2 = zeros(1,J);
O3 = zeros(1,J);
gO1 = zeros(2,J);
gO2 = zeros(2,J);
gO3 = zeros(2,J);



% iteration to solve the HJB equation
for k = 1:max_iter
    x1_old = x1; p1_old = p1;
    x2_old = x2; p2_old = p2;
    x3_old = x3; p3_old = p3;
    


    for j = 2:J
        %resolve obstacles
        if isempty(r)
                O1(j) = 1;
                gO1(:,j) = [0;0];
                O2(j) = 1;
                gO2(:,j) = [0;0];
                O3(j) = 1;
                gO3(:,j) = [0;0];
            else
                [D,m] = min((x1(1,j)-xC(1,:)).^2 + (x1(2,j)-xC(2,:)).^2-r.^2);
                O1(j) = 1/2 + (1/2)*tanh(B*D);
                gO1(:,j) = 2*B*sech(B*D)^2*(x1(:,j)-xC(:,m));
                [D,m] = min((x2(1,j)-xC(1,:)).^2 + (x2(2,j)-xC(2,:)).^2-r.^2);
                O2(j) = 1/2 + (1/2)*tanh(B*D);
                gO2(:,j) = 2*B*sech(B*D)^2*(x2(:,j)-xC(:,m));
                [D,m] = min((x3(1,j)-xC(1,:)).^2 + (x3(2,j)-xC(2,:)).^2-r.^2);
                O3(j) = 1/2 + (1/2)*tanh(B*D);
                gO3(:,j) = 2*B*sech(B*D)^2*(x3(:,j)-xC(:,m));
        end
        
        % update costate variables
        b = p1(:,j) + sig*(z1(:,j)-z1(:,j-1));
        p1(:,j) = max(0,1-sig*dt*O1(j)*v(x1(:,j))/norm(b,2))*b;
        b = p2(:,j) + sig*(z2(:,j)-z2(:,j-1));
        p2(:,j) = max(0,1-sig*dt*O2(j)*v(x2(:,j))/norm(b,2))*b;
        b = p3(:,j) + sig*(z3(:,j)-z3(:,j-1));
        p3(:,j) = max(0,1-sig*dt*O3(j)*v(x3(:,j))/norm(b,2))*b;
    end

    x1(:,1) = x1(:,1)+tau*p1(:,2);
    x2(:,1) = x2(:,1)+tau*p2(:,2);
    x3(:,1) = x3(:,1)+tau*p3(:,2);
    for j = 2:J-1
        
        % this is the marginal penalty function for pattern formation
        % (we only actually need its gradient)
        % P = (norm(x1(:,j)-x2(:,j),2)^2 - d^2)^2+...
        %    (norm(x1(:,j)-x3(:,j),2)^2 - d^2)^2+...
        %    (norm(x2(:,j)-x3(:,j),2)^2 - d^2)^2;

        % gradients of pattern function with respect to each agent
        gP1 = 2*((norm(x1(:,j)-x2(:,j),2)^2 - d^2)*(x1(:,j)-x2(:,j)) + (norm(x1(:,j)-x3(:,j),2)^2 - d^2)*(x1(:,j)-x3(:,j)));
        gP2 = 2*((norm(x2(:,j)-x1(:,j),2)^2 - d^2)*(x2(:,j)-x1(:,j)) + (norm(x2(:,j)-x3(:,j),2)^2 - d^2)*(x2(:,j)-x3(:,j)));
        gP3 = 2*((norm(x3(:,j)-x1(:,j),2)^2 - d^2)*(x3(:,j)-x1(:,j)) + (norm(x3(:,j)-x2(:,j),2)^2 - d^2)*(x3(:,j)-x2(:,j)));
        
        % gradients of travel time function with respect to each agent
        gC1 = 2*A*exp(-A*norm(x1(:,j)-xf(1:2),2)^2)*(x1(:,j)-xf(1:2));
        gC2 = 2*A*exp(-A*norm(x2(:,j)-xf(3:4),2)^2)*(x2(:,j)-xf(3:4));
        gC3 = 2*A*exp(-A*norm(x3(:,j)-xf(5:6),2)^2)*(x3(:,j)-xf(5:6));
    
        % gradient of the velocity function
        gV = [norm(p1(:,j),2)*gv(x1(:,j));norm(p2(:,j),2)*gv(x2(:,j));norm(p3(:,j),2)*gv(x3(:,j))];

        % update state variables via gradient descent
        XI = x1(:,j) - tau*(p1(:,j)-p1(:,j+1));
        x1(:,j) = x1(:,j) - gd_rate*(-dt*tau*(O1(j)*gV(1:2) + v(x1(:,j))*norm(p1(:,j),2)*gO1(:,j) - W1*gC1 - W2*gP1) + (x1(:,j)-XI));
        XI = x2(:,j) - tau*(p2(:,j)-p2(:,j+1));
        x2(:,j) = x2(:,j) - gd_rate*(-dt*tau*(O2(j)*gV(3:4) + v(x2(:,j))*norm(p2(:,j),2)*gO2(:,j) - W1*gC2 - W2*gP2) + (x2(:,j)-XI));
        XI = x3(:,j) - tau*(p3(:,j)-p3(:,j+1));
        x3(:,j) = x3(:,j) - gd_rate*(-dt*tau*(O3(j)*gV(5:6) + v(x3(:,j))*norm(p3(:,j),2)*gO3(:,j) - W1*gC3 - W2*gP3) + (x3(:,j)-XI));

    end
    
    % update relaxation variables
    z1 =  x1 + (x1-x1_old);
    z2 =  x2 + (x2-x2_old);
    z3 =  x3 + (x3-x3_old);

    % check convergence criterion
    change1 = max(max(max(abs(x1-x1_old))), max(max(p1-p1_old)));
    change2 = max(max(max(abs(x2-x2_old))), max(max(p2-p2_old)));
    change3 = max(max(max(abs(x3-x3_old))), max(max(p3-p3_old)));
    change = max([change1,change2,change3]);
    if change < tol
        break
    end
    
    % print progress
    if mod(k,1000) == 0
        fprintf('Completed iteration %i. change = %.4e\n',k,change);
        % after 5000 iterations, decrement gradient descent rate each 1000
        % iterations
        if  k>=5000
            gd_rate = gd_rate*0.5;
        end
    end
end

%compute the value function
u = 0;
for j = 2:J
    P = (norm(x1(:,j)-x2(:,j),2)^2 - d^2)^2+...
           (norm(x1(:,j)-x3(:,j),2)^2 - d^2)^2+...
           (norm(x2(:,j)-x3(:,j),2)^2 - d^2)^2;
   u = u + p1(:,j)'*(x1(:,j)-x1(:,j-1)) + p2(:,j)'*(x2(:,j)-x2(:,j-1))+ p3(:,j)'*(x3(:,j)-x3(:,j-1))...
        + dt*W1*(3 - exp(-A*norm(x1(:,1)-xf(1:2))^2)-exp(-A*norm(x2(:,1)-xf(3:4))^2)-exp(-A*norm(x3(:,1)-xf(5:6))^2))...
        + dt*W2*P...
        + dt*(O1(j)*v(x1(:,j))*norm(p1(:,j))+O2(j)*v(x2(:,j))*norm(p2(:,j))+O3(j)*v(x3(:,j))*norm(p3(:,j)));
end

% store all agents in single variables. 
x(1:2,:) = x1;
x(3:4,:) = x2;
x(5:6,:) = x3;
p(1:2,:) = p1;
p(3:4,:) = p2;
p(5:6,:) = p3;

% store the final iteration count
how_many_iter = k;
end

