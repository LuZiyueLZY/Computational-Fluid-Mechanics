% Verify the order of accuracy (OoA) of the compact center finite-difference 
% method for the first-order spatial derivative dy/dx

clear;


% In a loop, we consider a series of values of N (or dx)
for N=50:50:300

% initialisation
x=linspace(0,1,N+1);
x=x(1:N); dx = x(2)-x(1); x=x(:);

y = -232 * pi^2 * sin(2*pi*x) + 116 * pi * cos(2*pi*x);

% set up the derivative matrix

% D1
D1R = toeplitz([0;-1;zeros(N-2,1)], [0 1 zeros(1,N-2)]);
D1R(1,N) = -1; D1R(N,1)=1; % periodic bc
D1L = toeplitz([4;1;zeros(N-2,1)], [4 1 zeros(1,N-2)]);
D1L(1,N) = 1; D1L(N,1)=1; % periodic bc
D1 = D1L\D1R*3/dx;

% D2
D2R = toeplitz([-2;1;zeros(N-2,1)], [-2 1 zeros(1,N-2)]);
D2R(1,N) = 1; D2R(N,1)=1; % periodic bc
D2L = toeplitz([120;12;zeros(N-2,1)], [120 12 zeros(1,N-2)]);
D2L(1,N) = 12; D2L(N,1)=12; % periodic bc
D2 = D2L\D2R*144/(dx^2);

% LHS matrices
A = D2 + D1;

% Solution without fixed point
%u = A\y;

% solution with fixed point
u = zeros(length(A), 1);
u(2:end) = A(2:end, 2:end)\y(2:end);

end

% plot the solution
figure
plot(x, u,'-*b')
set(gca,'FontSize',30)
xlabel('x');ylabel('u')
