% Verify the order of accuracy (OoA) of the compact center finite-difference 
% method for the first-order spatial derivative dy/dx

clear;

% initialize the arrays for the error and dx. In the end, we will plot the
% error in terms of dx
er_array=[];
dx_array=[];

% Time loop
for i=1:200

    er_array=[];
    dx_array=[];

    % initialisation 
    dt=0.02; Nt=100;
    N=202;
    x=linspace(0,1,N);
    x=x(1:N); dx = x(2)-x(1); x=x(:);
    
    c=58;
    tau = 232;
    t=0;
    lambda = -58/100 * 4 * pi^2;
    
    u = sin(2*pi*x);
    
    % set up LHS matrix B
    B = toeplitz([1 + 2 * tau; -tau; zeros(N-2,1)], [1 + 2 * tau  -tau zeros(1,N-2)]);
    LHS = B(2:end-1, 2:end-1);
    
    % set up RHS matrix A
    A = toeplitz([1 - 2 * tau; tau;  zeros(N-2,1)], [1 - 2 * tau tau zeros(1,N-2)]);
    RHS = A(2:end-1, 2:end-1);


    % analytical solution of dy/dx
    dy_ana = exp(lambda * (t+dt)) * sin(2*pi*x);

    % numerical solution of dy/dx
    u_sol_1 = LHS\RHS * u(2:end-1);

    u_sol = zeros(N, 1);
    u_sol(2:end-1) = u_sol_1;
    dy_num = u_sol;

    % Space Loop
    for j=1:N
       
        % save the RMS error and dx
        dx_array=[dx_array x(j)];
        er_array=[er_array sqrt(sum((dy_ana(j)-dy_num(j)).^2)/N)];
    end
    % plot in a log-log figure 
    if mod(i, 50) == 0
        figure
        loglog(dx_array,er_array,'-*b')
        hold on
        loglog(dx_array,dx_array.^2,'--r')
        legend('RMS error of compact center FD') %,'Er = C dx^2')
        set(gca,'FontSize',30)
        xlabel('dx');ylabel('RMS Error')
        
    end

    t=t+dt;
end