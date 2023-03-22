% Implementation of a compact scheme for the first-order spatial
% derivative. The compact scheme is high-order FD method.

clear;

% initialisation
N=100;
x=linspace(0,1,N+1);
x=x(1:N); dx = x(2)-x(1); x=x(:);

y = sin(2*pi*x);

% set up the derivative matrix. This is a compact scheme, so you need to
% invert a matrix
D2R = toeplitz([30;-16;1;zeros(N-3,1)], [30 -16 1 zeros(1,N-3)]);
D2R(1, N-1)=1; D2R(N-1, 1)=1; D2R(2, N)=1; D2R(N, 2)=1;
D2R(1, N)= -16; D2R(N, 1) = -16; % periodic bc
%{
D2L = toeplitz([4;1;zeros(N-2,1)], [4 1 zeros(1,N-2)]);
D2L(1,N) = 1; D2L(N,1)=1; % periodic bc
%}
D2 = -1/(12 * dx^2) * D2R ;

% analytical solution of dy/dx
dy_ana = -2*pi*2*pi*sin(2*pi*x);

% numerical solution of dy/dx
dy_num = D2*y;

% plot the result together
figure
plot(x,dy_ana,'-*b')
hold on
plot(x,dy_num,'-or')
title('Compact center FD method')
legend('Analytical','Numerical')
set(gca,'FontSize',30)
xlabel('x');ylabel('y')