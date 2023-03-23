% Slide 6 in Week 6. We solve a 1D convection equation using the Lax-Wendroff
% method. An implementation of the factorisation step is also shown.

clear
clc

% initialisation 
dt=0.01; Nt=100;
N=100;
t=0;
x=linspace(0,1,N+1);
x=x(1:N); dx = x(2)-x(1); x=x(:);

nu=0.001;

for V=1:10
    % set up LHS matrix B
    LHS = toeplitz([((1/dt) + (2*nu)/(2*(dx^2))); -(nu/(2*(dx^2))); zeros(N-2,1)], [((1/dt) + (2*nu)/(2*(dx^2))) -(nu/(2*(dx^2)))  zeros(1,N-2)]);
    
    % set up RHS matrix A
    RHS = toeplitz([((1/dt) + (V/dx) - (2*nu)/(2*(dx^2))); (nu/(2*(dx^2)));  zeros(N-2,1)], [((1/dt) + (V/dx) - (2*nu)/(2*(dx^2))) (-(V/dx) + (nu/(dx^2))) zeros(1,N-2)]);

    G = LHS\RHS;

    [eigvec, eigval] = eig(G);

    figure
    plot(diag(eigval),'-*b')
    ylim([-1,1])
    title(['V = ' num2str(V)])
    set(gca,'FontSize',30)

    pause(0.0001)

end






