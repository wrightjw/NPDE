clear all, close all

%% Question 1:

% This script finds a numerical approximation for the solution to the
% initial value problem:
%   u_t = 0.5 u_xx  for  x \in (-1,1),  t \in (0,T),
%   u(x,t=0) = sin(pi*x) -3sin(3*pi*x),
%   u(x=-1,t) = u(x=1,t) = 0,
% using the Crank-Nicolson method.

% The grid spacing in the x-dimension
dx = 0.01;

% The x-dimension grid locations
x = (-1:dx:1)';

% The number of x-dimension grid divisions
M = length(x)-1; %M = 2/dx instead of 1/dx??? wat do??

% The time step size
dt = 0.25*dx;
% dt = C*dx^2; % Uncomment for Exercises 1.3 and 1.4

% The time integration length and coefficient k
T = 0.25;
k = 0.5;

% The number of timesteps
N = floor(T*(1+1e-10)/dt);

% The initial condition
U_n = [0; sin(pi*x(2:end-1))-3*sin(3*pi*x(2:end-1)); 0];
U_np1 = zeros(size(x)); % "np1" denotes "n plus 1"

% The global error associated with the initial condition is exactly zero.
% Use this property to initialise the norm calculations.
norm_2 = 0;
norm_2_dx = 0;

% Assemble matrix which applies the Crank-Nicolson method - BU(t+1) = AU(t)
A = sparse(M-1,M-1); %Sparse matrix for LHS (t-1) of C-N rearrangement
B = sparse(M-1, M-1); %Sparse matrix for RHS (t+1) of C-N rearrangement
C = sparse(M-1, M-1); %Sparse matrix to compute A/B to solve for t+1
mu = 0.5*k*dt/dx^2;

%Inserting matrix entries for A and B for AU=B
for i = 1:M-1
  A(i,i) = 1-2*mu;
  B(i,i) = 1+2*mu;
end
for i = 1:M-2
  A(i+1,i) = mu; A(i,i+1) = mu;
  B(i+1,i) = -mu; B(i,i+1) = -mu;
end

%Matrix that solves C-R expression explicitly for time t+1
C = mldivide(B,A);

for n = 0:N-1
  % A single Crank-Nicolsons timestep
  U_np1(2:end-1) = C*U_n(2:end-1);
  
  % The exact solution after (n+1) timesteps
  u_np1 = exp(-pi^2*k*(n+1)*dt).*sin(pi*x) ...
      -3*exp(-9*pi^2*k*(n+1)*dt).*sin(3*pi*x);
  
  % The global error after (n+1) timesteps
  e_np1 = U_np1-u_np1;
  
  % Update the error norms norm_2 and norm_2_dx
  norm_2 = max(norm_2,norm(e_np1));
  norm_2_dx = max(norm_2_dx,sqrt(dx)*norm_2);
  
  % Plot the solution at the (n+1)th timestep
  figure(1)
  set(gcf, 'Position', [100 300 700 700])
  plot(x,U_np1,'k-')
  title('Numerical Solution to Heat Equation Using Crank-Nicolson')
  xlabel('$x$','FontSize',16,'Interpreter','latex')
  ylabel('Numerical solution','FontSize',16)
  xlim([-1,1])
  ylim([-5,5])
  text(0.55,2,['$t=',num2str((n+1)*dt,'%.4f'),',n=',num2str(n+1),'$'], ...
      'FontSize',16,'Interpreter','latex')
  drawnow()
  
  U_n = U_np1;
end

% Display the error norms
disp(['\max_{n \in \{0,...,N\}} \| \underline{e}^n \|_2 = ', ...
    num2str(norm_2,'%.16e')])
disp(['\max_{n \in \{0,...,N\}} \| \underline{e}^n \|_{2,\Delta x} = ', ...
    num2str(norm_2_dx,'%.16e')])

% The maximum integers a and b asked for in the question are 
% a = 3
% b = 4

%% Question 2:

% This script finds a numerical approximation for the solution to the
% initial value problem:
% -u_xx - u_yy = -2x(x-1)-2y(y-1)+2pi^2*sin(pi*x)sin(pi*y)  for  x,y \in (0,1)^2,
% u(x=0,y)=u(x=1,y)=u(x,y=0)=u(x,y=1)=0,
% using 5-point Laplacian finite difference scheme.
% Function at bottom of script.
 
% Asking user for the grid spacing in the x-dimension
prompt = 'Give dx: ';
dx = input(prompt);

%Running finite difference scheme function for user defined dx.
[M, U, u] = FS(dx);

%Computing and displaying error norms
e=U-u; %the difference between numerical and exact solution for each m and p.
e_squared = e.^2;
norm_2 = sqrt(sum(e_squared));
norm_2_dx = dx*norm_2;
disp(['\max_{n \in \{0,...,N\}} \| \underline{e}^n \|_2 = ', ...
    num2str(norm_2,'%.16e')])
disp(['\max_{n \in \{0,...,N\}} \| \underline{e}^n \|_{2,\Delta x} = ', ...
    num2str(norm_2_dx,'%.16e')])

%Formatting numerical solution into matrix rows and columns for each 
% m and p for plotting, where dx = 0.005.
U_surf = zeros(M-1);
c=1; %initialising count for loop
if dx == 0.005
    for m = 1:M-1
        for p = 1:M-1
        U_surf(m,p)=U(c);
        c = c+1;
        end
    end
else
    %Setting dx=0.005 for plotting and rerunning algorithm
    dx=0.005;
    [M, U, u] = FS(0.005);
    for m = 1:M-1
        for p = 1:M-1
        U_surf(m,p)=U(c);
        c = c+1;
        end
    end
end

%Appending boundary conditions - noting these are all zero
U_plot = zeros(M+1); %initialising matrix with boundaries
U_plot(2:M,2:M) = U_surf; %embedding solution matrix


%Plotting numerical solution - we ignore the boundary as these
%are all 0 anyway.
figure(2)
set(gcf, 'Position', [800 300 700 700])
[X,Y]=meshgrid(0:dx:1,0:dx:1);
surf(X,Y,U_plot);
shading interp;
colormap cool;
title('Finite Difference Solution of Poisson Equation, $dx = 0.005$', 'Interpreter', 'Latex');
xlabel('$x$','FontSize',16,'Interpreter','latex');
ylabel('$y$','FontSize',16, 'Interpreter','latex');
zlabel('Numerical Solution');

%Function that performs finite difference scheme on poisson's equation for
%given dx
%inputs: dx
%outputs: [M, U, u]
function [M, U, u] = FS(dx)
    % The number of x-dimension grid divisions
    M = floor((1+1e-10)/dx);

    % Initialise matrices A and B, and the vector b
    A = sparse((M-1)^2,(M-1)^2); B = sparse(M-1,M-1); b = zeros((M-1)^2,1);
 
    % Assemble matrix B
    for i = 1:M-1
        if i == 1
            B(i,i) = 4;
            B(i, i+1) = -1;
        elseif i == M-1
            B(i,i-1) = -1;
            B(i,i) = 4;
        else
            B(i,i) = 4;
            B(i, i+1) = -1;
            B(i,i-1) = -1;
        end
    end
 
    % Assemble matrix A
    for i = 1:M-1
        A((i-1)*(M-1)+1:i*(M-1),(i-1)*(M-1)+1:i*(M-1)) = B; %Matrix B on diagonal
    end
    for i = 1:M-2
        A((i-1)*(M-1)+M:i*(M-1)+(M-1),(i-1)*(M-1)+1:i*(M-1)) = -1*eye(M-1); %-I above B
        A((i-1)*(M-1)+1:i*(M-1), (i-1)*(M-1)+M:i*(M-1)+(M-1)) = -1*eye(M-1); %-I below B
    end
    A = A/dx^2;

    c=1; %Initialising count to put values into b vector
    %Assemble vector b
    for m = 1:(M-1)
        for p = 1:(M-1)
            b(c)=-2*(m*dx)*((m*dx)-1)-2*(p*dx)*((p*dx)-1)+2*pi^2*sin((m*dx)*pi)*sin((p*dx)*pi);
            c=c+1;
        end
    end

    %Solving matrix system Au=b for u
    U = A\b;

    %Computing exact solution, ignoring boundaries as the global error is 0.
    u = zeros(size(U)); c=1; %initialising states for loop
    for m = 1:M-1
        for p = 1:M-1
            u(c) = (m*dx)*((m*dx)-1)*(p*dx)*((p*dx)-1)+sin(pi*(m*dx))*sin(pi*(p*dx)); 
            c=c+1;
        end
    end
end

% The maximum integers c and d asked for in the question are 
% c = 2
% d = 4