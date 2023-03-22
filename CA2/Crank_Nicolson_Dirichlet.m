% Slide 6 in Week 6. We solve a 1D convection equation using the Lax-Wendroff
% method. An implementation of the factorisation step is also shown.

clear
clc

% initialisation 
dt=0.02; Nt=100;
N=202;
x=linspace(0,1,N);
x=x(1:N); dx = x(2)-x(1); x=x(:);

c=58;
tau = ((58/10000) * dt)/(2*(dx^2));
t=0;
lambda = -58/10000 * 4 * pi^2 * (3/4)^2;

% initial condition
u = sin(2*pi*x*(3/4));
u(1) = 0;


% set up LHS matrix B with Dirichlet boundary condition
B = toeplitz([1 + 2 * tau; -tau; zeros(N-2,1)], [1 + 2 * tau  -tau zeros(1,N-2)]);
LHS = B(2:end-1, 2:end-1);
LHS(end, end) = 467 / 3; % LHS(N-1. N-1)
LHS(end, end-1) = -464 / 3; % LHS(N-2. N-1)

% set up RHS matrix A with Dirichlet boundary condition
A = toeplitz([1 - 2 * tau; tau;  zeros(N-2,1)], [1 - 2 * tau tau zeros(1,N-2)]);
RHS = A(2:end-1, 2:end-1);
RHS(end, end) = -461 / 3; % RHS(N-1, N-1)
RHS(end,end-1) = 464 / 3; % RHS(N-2, N-1)

% creation of solution matrix
for i=1:Nt
    % solve
    u_sol_1 = LHS\RHS * u(2:end-1);

    u_sol = zeros(N, 1);
    u_sol(2:end-1) = u_sol_1;
    u = u_sol;

    u_ana1 = exp(lambda * (t+dt)) * sin(2*pi*x*(3/4));
    u_ana = u_ana1(2:end-1);

    t=t+dt;
   
    plot(x(2:end-1),u_sol_1,'-*b', x(2:end-1), u_ana, '-or')
    ylim([-1,1])
    legend('Numerical','Analytical')
    title(['Time = ' num2str(t)])
    set(gca,'FontSize',30)

    pause(0.001)

end





