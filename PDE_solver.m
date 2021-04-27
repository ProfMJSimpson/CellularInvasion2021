% Filename: PDE_solver.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, April 2021
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2021) 
% Travelling wave analysis of cellular invasion into surrounding tissues
% This file script contains:
%   - two calls to function InvasionModel that solves Equations(3)-(4)
%   with IC from Equations (5)-(6) using Newton-Raphson method. It displays
%   time-dependent solutions u(x,t) and v(x,t) as in Figures 2(b) and (h).
%	- function InvasionModel
%   - function tridia
% This function generates and displays the solutions u(x,t) and v(x,t)
% of invasion model with following parameters gamma = {1, 100}, tmax = 120,  
% dx = 0.01, dt = 0.01, L = 250, Vinit={0.5, 1}

% Generating Figure 2(b)
% Calling function InvasionModel(gamma, Vinit, L, tmax, dx, dt,
% timeToPrint) with following inputs: gamma = 1, Vinit=0.5, tmax = 120,
% dx = 0.01, dt = 0.01, and array of times at which we must print
% solutions u(x,t) and v(x,t), t=[0, 40, 80, 120]
InvasionModel(1, 0.5, 250, 120, 0.01, 0.01, [0 40 80 120]);

% Generating Figure 2(h)
% Calling function InvasionModel(gamma, Vinit, L, tmax, dx, dt,
% timeToPrint) with following inputs: gamma = 100, Vinit = 1, tmax = 120,
% dx = 0.01, dt = 0.01, and array of times at which we must print
% solutions u(x,t) and v(x,t), t=[0, 100, 110, 120]
InvasionModel(100, 1, 250, 120, 0.01, 0.01, [0 100 110 120]);

% Function InvasionModel
% This function solves Equations (3)-(4) with IC from Equations (5)-(6)
% using Newton-Raphson method. It displays time-dependent solutions 
% u(x,t) and v(x,t) 
% at requested times.
% INPUT ARGUMENTS:
% -- gamma: parameter gamma, a positive constant not equal to zero
% -- Vinit: constant to which we set initial density v(x,0), such as
% 0 < Vinit <=1
% -- L: total length of domain
% -- tmax: final time of solution, starting at t=0
% -- dx: spacing between two nodes in spatial grid
% -- dt: time step
% -- timeToPrint : an array of desired times when to display solutions
function InvasionModel(gamma, Vinit, L, tmax, dx, dt, timeToPrint)

    % colors used to display density profiles
    % u(x,t) is in brown and v(x,t) is in blue
    colors = [0 0.4470 0.7410;0.9290 0.6940 0.1250;];
    % tolerance factor
    tol = 1e-08;
    % total number of time steps
    ts = round(tmax/dt);
    % spatial grid going from 0 to L with spacing dx
    x = 0:dx:L;
    % total number of nodes in spatial grid
    nodes = length(x);
    
    % initialisation of variables used in Newton-Raphson algorithm
    % previous density profiles u(x,t) and v(x,t)
    u_p = zeros(1,nodes);
    v_p = zeros(1,nodes);
    % correction vector  for u and v
    delta_u = ones(1,nodes);
    delta_v = ones(1,nodes);
    % function F
    F = zeros(1,nodes);
    % coefficients a b c of the tridiagonal matrix
    % Jacobian
    coeffA = zeros(1,nodes);
    coeffB = zeros(1,nodes);
    coeffC = zeros(1,nodes);

    % array to store all current positions x* where u(x,t)=0.5
    PostionArray = zeros(1,ts);
    
    % initialisation of density profiles using Equations (5)-(6)
    % where beta = 10 and alpha = 1
    for i=1:nodes
        if (i*dx < 10)
            u_p(1,i) = 1;
        end
        if (i*dx >= 10)
            u_p(1,i) = 0;
        end
        v_p(1,i) = Vinit;
    end
    u = u_p;
    v = v_p;
    
    % opening a new figure
    figure
    hold on
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 18);
    %% Newton-Raphson algorithm - main loop
    % for all time steps
    for j = 1:ts
        condition = 1;
        % while tolerance is not reached
        while (condition)
            % Neumann boundary condition for u at x = 0
            coeffA(1,1) = 0;
            coeffB(1,1) = -1.0;
            coeffC(1,1) = 1.0;
            F(1) = -1.0*(u(1,2) - u(1,1));

            % Neumann boundary condition for u at x = L
            coeffA(1,nodes) = -1.0;
            coeffB(1,nodes) = 1.0;
            coeffC(1,nodes) = 0;
            F(1,nodes) = -1.0*(u(1,nodes) - u(1,nodes-1));

            % J(u) delta u = -F(u)
            % Equation (3)
            for i = 2:nodes-1
                coeffA(1,i) = (2-(v(1,i)+v(1,i-1)))/(2*dx^2);
                coeffB(1,i) = -((2-(v(1,i+1)+v(1,i)))...
                    +(2-(v(1,i)+v(1,i-1))))/(2*dx^2)+(1-2*u(1,i)-v(1,i))...
                    -1/dt;
                coeffC(1,i) = (2-(v(1,i+1)+v(1,i)))/(2*dx^2);
                F(1,i) = -(((2-(v(1,i+1)+v(1,i)))*(u(1,i+1)-u(1,i))...
                    - (2-(v(1,i)+v(1,i-1)))*(u(1,i)-u(1,i-1)))/(2*dx^2)...
                    + (u(1,i)*(1-u(1,i)-v(1,i)))-(u(1,i)-u_p(1,i))/dt);
            end         
            delta_u = tridia(coeffA, coeffB, coeffC, F, nodes);

            % correction of u(x,t)
            u(1,:) = u(1,:) + delta_u(1,:);

            % Equation (4)
            % J(v) delta v = -F(v)
            for i = 1:nodes
                coeffA(1,i) = 0;
                coeffB(1,i) = -gamma*u(1,i)  - 1.0/dt;
                coeffC(1,i) = 0;
                F(1,i) = -(-gamma*u(1,i)*v(1,i)- 1.0/dt*(v(1,i) - v_p(1,i)));
            end
            delta_v = tridia(coeffA, coeffB, coeffC, F, nodes);
            
            % correction of v(x,t)
            v(1,:) = v(1,:) + delta_v(1,:);

            % checking if tolerance is reached
             if (norm(delta_u,Inf) <= tol && norm(delta_v,Inf) <= tol)
                condition = 0;
            end           
        end

        % updating current u(x,t) and v(x,t)
        u_p(1,:) = u(1,:);
        v_p(1,:) = v(1,:);

        % finding and storing current position of u(x,t) = 0.5
		for inode = 1:nodes-1
			if (u(1,inode) >= 0.5 && u(1,inode+1) <= 0.5)
                x0 = inode*dx;
                x1 = (inode+1)*dx;
                y0 = u(1,inode);
                y1 = u(1,inode+1);
                PostionArray(1,j) = x0 + (0.5-y0)*(x1-x0)/(y1-y0); 
			end
        end

       % displaying current solutions u(x,t) and v(x,t) at requested 
       % times
       if (isempty(find(round(timeToPrint/dt) == j,1)) == 0)
            plot(x, u_p, 'm-', 'LineWidth', 2,'Color', colors(1,1:3));
            plot(x, v_p, 'y-', 'LineWidth', 2,'Color', colors(2,1:3));
       end
        
       % display current time at each 100 time steps
       if (mod(j,100)==0)
          disp(strcat('t=', num2str(j*dt,5)));
        end
    end
    
    % computing wave speed at each time step
    WaveSpeedArray = (PostionArray(1,2:end)-PostionArray(1,1:length(PostionArray)-1))/(dt);
    % computing final wave speed using an average over 101 last values
    WaveSpeed = sum(WaveSpeedArray(end-100:end))/101;
    
    % printing value of gamma
    textg = strcat(strcat('$\gamma = ',num2str(round(gamma,2))), '$');
    text(200,0.6,textg,'interpreter','latex','fontsize',18);
    % printing wave speed
    textc = strcat(strcat('$c = ',num2str(round(WaveSpeed,2))), '$');
    text(200,0.7,textc,'interpreter','latex','fontsize',18);
    hold off
    % setting limits of x-axis and y-axis
    xlim([0 L])
    ylim([0 1])
    box on
    % printing axis labels
    ylabel('$u(x,t) \quad v(x,t)$','interpreter','latex');
    xlabel('$x$','interpreter','latex');
end

%% Function tridia
% This function implements Thomas algorithm that solves a system
% of equations Ax = d, where A is a tridiagonal matrix. The parameters 
% a,b and c are three diagonals of matrix A. N is size of 
% vector solution x.
function x = tridia(a,b,c,d,N)
    x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end

    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end 
