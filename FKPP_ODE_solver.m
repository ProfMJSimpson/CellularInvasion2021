% Filename: FKPP_ODE_solver.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, April 2021
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2021) 
% Travelling wave analysis of cellular invasion into surrounding tissues
% This file script contains:
%   - two calls to function plotPhasePlane to generate travelling wave
%     solution of Fisher-KPP ODE in the phase plane. It displays green 
%     dashed curve as seen in Figures 4(b) and (e).
%	- function plotPhasePlane
%   - function HeunSolver

% Generating Figure 4(b)
% Calling function plotPhasePlane(c,dz,z_begin,z_end) with 
% following inputs: c = 1.24, dz = 0.001, z_begin = -50, z_end = 50
plotPhasePlane(1.24,0.001,-50,50);

% Generating Figure 4(e)
% Calling function plotPhasePlane(c,dz,z_begin,z_end) with 
% following inputs: c = 1.75, dz = 0.001, z_begin = -50, z_end = 50
plotPhasePlane(1.75,0.001,-50,50);

% Function plotPhasePlane
% This function solves Equations (20) and (21) in the phase plane
% by Heun's method and plots solution W(z) versus U(z) on the phase plane.
% The same plot shows also the equilibrium points (0,0) and (1,0).
% INPUT ARGUMENTS:
% -- c, positive wave speed
% -- dz, step size used to discretise domain of z.
% -- z_begin and z_end, lower and upper limits of domain of z, used
% to integrate numerically Equations (11) and (12) by Heun's method,
% such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
function plotPhasePlane(c,dz,z_begin,z_end)

    % green color used to display travelling wave solution
    colors = [0.4660 0.6740 0.1880;];

    % calculating eigenvalues and eigenvectors of solution 
    % around saddle (1,0)
    Us = 1;
    A = [0 1;1 -c];
    [eigenvec,eigenval] = eig(A);
    % determining initial conditions close to equilibrium point
    % (1,0) along eigenvector of solution
    % setting initial conditions to start on unstable manifold
    % to obtain invading travelling wave solution
    IC = [0;0];
    if (eigenval(1,1) > 0)
        IC = eigenvec(:,1);
    end
    if (eigenval(2,2) > 0)
        IC = eigenvec(:,2);
    end
    % initial conditions to obtain invading travelling wave solution
    IC = [Us+IC(2,1)*0.0001; IC(1,1)*0.0001;];
    
    % solving Equations (20) and (21) in the phase plane with
    % Heun's method
    [U, W] = heunSolver(c, dz, z_begin, z_end, IC);
    
    % opening a new figure
    figure
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 18);
    % plotting axis at U(z) and W(z)
    line([-0.5 1.05],[0 0],'Color','k','LineStyle','-','LineWidth',1);
    hold on
    line([0 0],[-0.5 0.05],'Color','k','LineStyle','-','LineWidth',1);
    % plotting equilibrium points (0,0) and (1,0)
    plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    plot(1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    % plotting solution 
    plot(U, W, '--', 'LineWidth', 2,'color',colors(1,:));
    % displaying value of wave speed c
    textg = strcat(strcat('$c = ',num2str(round(c,2))), '$');
    text(0.2,-0.1,textg,'interpreter','latex','fontsize',18);   
    box on
    % printing axis labels
    ylabel('$W(z)$','interpreter','latex');
    xlabel('$U(z)$','interpreter','latex');
    % setting limits of x-axis and y-axis
    % minimum value of W(z) used to limit y axis
    ymin = min(W);
    ylim([ymin-0.05 0.05])
    xlim([-0.2 1.05])
end

% Function heunSolver
% This function solves Equations (20) and (21) by Heun's method 
% INPUT ARGUMENTS:
% -- c, positive wave speed
% -- dz, step size used to discretise domain of z
% -- z_begin and z_end, lower and upper limits of numerical domain 
% of z such as z_begin <= z <= z_end. Initial conditions are applied at
% z = z_begin.
% -- IC: an array of values of initial conditions U(0) and W(0)
% OUTPUT ARGUMENTS:
% ** U : The solution U(z)
% ** W : The solution V(z)
function [U, W] = heunSolver(c,dz,z_begin,z_end,IC)
    % domain of z
    z = z_begin:dz:z_end;
    % number of nodes in z-domain
    sz = length(z);

    % initialisation 
    W = zeros(sz,1);
    U = zeros(sz,1);
    U(1) = IC(1);
    W(1) = IC(2);

    % integrating Equations (20) and (21) with Heun's method
    for i = 1:sz-1
        Ubar = dz * W(i) + U(i); 
        Wbar = dz * (-c*W(i)- U(i)*(1-U(i))) + W(i);
        U(i+1) = dz/2 * (W(i)+Wbar) + U(i); 
        W(i+1) = dz/2 * ((-c*W(i) - U(i)*(1-U(i))) + (-c*Wbar - Ubar *(1-Ubar))) + W(i); 
    end
 
end
