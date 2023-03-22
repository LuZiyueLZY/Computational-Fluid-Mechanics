% Verify the order of accuracy (OoA) of the compact center finite-difference 
% method for the first-order spatial derivative dy/dx

clear;

% initialize the arrays for the error and dx. In the end, we will plot the
% error in terms of dx
er_array=[];
dx_array=[];

% In a loop, we consider a series of values of N (or dx)
for N=50:50:300

% initialisation
x=linspace(0,1,N+1);
x=x(1:N); dx = x(2)-x(1); x=x(:);

y = sin(2*pi*x);

% set up the derivative matrix
D2R = toeplitz([-2;1;zeros(N-2,1)], [-2 1 zeros(1,N-2)]);
D2R(1,N) = 1; D2R(N,1)=1; % periodic bc
D2L = toeplitz([120;12;zeros(N-2,1)], [120 12 zeros(1,N-2)]);
D2L(1,N) = 12; D2L(N,1)=12; % periodic bc
D2 = D2L\D2R*144/(dx^2);

% analytical solution of dy/dx
dy_ana = -2*pi*2*pi*sin(2*pi*x);

% numerical solution of dy/dx
dy_num = D2*y;

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
