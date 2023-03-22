% Verify the order of accuracy (OoA) of the compact center finite-difference 
% method for the first-order spatial derivative dy/dx

clear;

% initialize the arrays for the error and dx. In the end, we will plot the
% error in terms of dx
er_array=[];
dx_array=[];
diff = [];

% In a loop, we consider a series of values of N (or dx)
for N=50:50:300

% initialisation
x=linspace(0,1,N+1);
x=x(1:N); dx = x(2)-x(1); x=x(:);

y = cos(2*pi*x);

% set up the derivative matrix
D1R = toeplitz([0;-1;zeros(N-2,1)], [0 1 zeros(1,N-2)]);
D1R(1,N) = -1; D1R(N,1)=1; % periodic bc
D1L = toeplitz([4;1;zeros(N-2,1)], [4 1 zeros(1,N-2)]);
D1L(1,N) = 1; D1L(N,1)=1; % periodic bc
D1 = D1L\D1R*3/dx;

% analytical solution of dy/dx
dy_ana = - 2 * pi * sin(2 * pi * x);

% numerical solution of dy/dx
dy_num = D1*y;

% save the RMS error and dx
dx_array=[dx_array dx];
er_array=[er_array sqrt(sum((dy_ana-dy_num).^2)/N)];
end

% plot in a log-log figure
figure
loglog(dx_array,er_array,'-*b')
hold on
loglog(dx_array,dx_array.^4*20,'--r')
legend('RMS error of compact center FD','Er = C dx^4')
set(gca,'FontSize',30)
xlabel('dx');ylabel('RMS Error')
