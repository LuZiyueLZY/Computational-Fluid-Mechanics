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
lambda = -58/10000 * 4 * pi^2;

% initial condition
u = sin(2*pi*x);
u(1) = 0;


% % set up LHS matrix B
% B = toeplitz([1 + 2 * tau; -tau; zeros(N-2,1)], [1 + 2 * tau  -tau zeros(1,N-2)]);
% LHS = B(2:end-1, 2:end-1);
% 
% % set up RHS matrix A
% A = toeplitz([1 - 2 * tau; tau;  zeros(N-2,1)], [1 - 2 * tau tau zeros(1,N-2)]);
% RHS = A(2:end-1, 2:end-1);

% % G = LHS\RHS
% G = LHS\RHS;
% [eigvec, eigval] = eig(G);

% checking sigma_i (eigen values)
for i=1:1000
    for j=1:1000
        tau = ((58/100) * i* dt)/(2*((j * dx)^2));

        % set up LHS matrix B
        B = toeplitz([1 + 2 * tau; -tau; zeros(N-2,1)], [1 + 2 * tau  -tau zeros(1,N-2)]);
        LHS = B(2:end-1, 2:end-1);
        
        % set up RHS matrix A
        A = toeplitz([1 - 2 * tau; tau;  zeros(N-2,1)], [1 - 2 * tau tau zeros(1,N-2)]);
        RHS = A(2:end-1, 2:end-1);

        % G = LHS\RHS
        G = LHS\RHS;
        [eigvec, eigval] = eig(G);

%         disp(eigval(100,100))
%         disp(000000)

        for k=1:200
            if abs(eigval(k,k)) >= 1
                disp(i)
                disp(j)
                disp(k)
            end
        end
    end
end

