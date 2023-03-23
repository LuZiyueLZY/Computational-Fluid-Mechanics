%(a)
% Using a finite-difference method to solve the Poission equation
% d^2u/dx^2 + d^2u/dy^2 = 1 with x and y from 0 to 1
% boundary conditions are 0 on the three edges and a sine wave on one edge
% inital conditions can also be varied
% we use iterative methods to solve the numerical equations

clear
close all

%%%% initialisation
C = (019 + 020 + 015)/3
Ni_x=60;Ni_y=50;                        % Grid resolution and domain length
%Ni_x=20;Ni_y=10;                        % Grid resolution and domain length

Lx=1;Ly=1;

x=linspace(0,Lx,Ni_x+1);dx=x(2)-x(1);   % Set up the grids
y=linspace(0,Ly,Ni_y+1);dy=y(2)-y(1);

[A,b]=assembleAb(Ni_x,Ni_y,dx,dy,x,y,C);  % Assemble the A matrix and b vector
% the boundary condition here is smoother. Note the implementaiton of the
% boundary condition in this subfunction

%% question b LU decomposition
N=size(A,1);

% Lp matrix starts with an identity matrix. Initialize the U matrix (which
% is redundant).
Lp=eye(N,N);

% make a copy of A in U_Lu and we do the row operation on U
U_Lu=A;

% the loop to go from the 1st column to the second last column (N-1)
auxL_start=eye(N,N);    % auxiliary L matrix starts with an identity matrix
for i=1:N-1
    auxL=auxL_start;    % auxiliary L matrix starts with an identity matrix
    
    auxL(i+1:N,i) = - U_Lu(i+1:N,i)/U_Lu(i,i);  % the equation to calculate the coefficients in the L matrix
    U_Lu(i+1:N,:) = U_Lu(i+1:N,:) + auxL(i+1:N,i)*U_Lu(i,:);  % updating the U matrix using the row operations. 
                                                     % You are using U(j,:) row to do these row operations.
    Lp=sparse(auxL)*sparse(Lp);   % save the auxiliary L matrix. 
                  % Remember that the last operation will left-multiply the previous operations. 
                  % So auxL appears on the left in auxL*Lp.
end

% We got A=LU. Now we can use it to solve Au=b.

% b vector. When boundary conditions change or the source terms change, you
% can still use the LU decomposition of A above to solve Au=b (for
% different b).

% get the intermediate vector y
y_Lu=Lp*b;

% solve for u. Note that here I just use inv(U_Lu) directly. 
u_Lu=U_Lu\y_Lu;
plot_results_solution_only(u_Lu,Ni_x,Ni_y,x,y,"(a) u_Lu")
%QR decomposition
N=size(A,1);

Q = zeros(N,N);

% v1=u1
Q(:,1)=A(:,1);

% get v2,v3...,vj,vN in the orthogonal matrix
for j=2:N
    Proj=zeros(N,1);
    % to get vj, you first do projection of uj on v1--v_k---v_{j-1}
    for k=1:j-1
        Proj = Proj + ( Q(:,k)'*A(:,j)  )/( Q(:,k)'*Q(:,k)   )*Q(:,k);
        %"Hello-1"+num2str(k)
    end
    % and then subtract the projection from uj
    Q(:,j)  =  sparse(A(:,j))  - Proj;
    "Hellorere"+num2str(j)
end

% normalisation
for j=1:N
    Q(:,j) = Q(:,j) /   sqrt( Q(:,j)'*Q(:,j) );
    "Hellokekeke"+num2str(j)
end

R=zeros(N,N);

% get R
for j=1:N
    % you do the projection of uj on e1,e2,...,e_j
    for k=1:j
        R(k,j) = Q(:,k)'*sparse(A(:,j));
    end
    "Hellojejeje"+num2str(j)
end

y_QR=Q'*b;
u_QR=R\y_QR;
plot_results_solution_only(u_QR,Ni_x,Ni_y,x,y,"(a) u_QR")

%%%% (c) Start the iterative procedure

% initial guess and specify the objective of residual
u0 = rand((Ni_x-1)*(Ni_y-1),1);
%u0 = -0.5*ones((Ni_x-1)*(Ni_y-1),1);
resobj = 10^-7;

% Jacobi method (method=1), Gauss-Seidel (method=2), SOR (method=3)
method = 1;
[solution_1,k_1,resarray_1]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u0,resobj,u_QR,nan);

method = 2;
[solution_2,k_2,resarray_2]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u0,resobj,u_QR,nan);

method = 3;
[solution_3,k_3,resarray_3]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u0,resobj,u_QR,1.5);

%%%% plot the results for each method
plot_results(k_1,resarray_1,solution_1,Ni_x,Ni_y,x,y,"*r","(c) Jacobi")
plot_results(k_2,resarray_2,solution_2,Ni_x,Ni_y,x,y,"*b","(c) Gauss-Seidel")
plot_results(k_3,resarray_3,solution_3,Ni_x,Ni_y,x,y,"*g","(c) SOR, w = 1.5") % omega = 1.5

%as omega increase, the number of iterations fall
method = 3;
[solution_4,k_4,resarray_4]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u0,resobj,u_QR,1.1);
[solution_5,k_5,resarray_5]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u0,resobj,u_QR,1.3);
[solution_6,k_6,resarray_6]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u0,resobj,u_QR,1.7);

plot_results(k_4,resarray_4,solution_4,Ni_x,Ni_y,x,y,"^c","(c) SOR, w = 1.1") % omega = 1.1
plot_results(k_5,resarray_5,solution_5,Ni_x,Ni_y,x,y,"^m","(c) SOR, w = 1.3") % omega = 1.3
plot_results(k_6,resarray_6,solution_6,Ni_x,Ni_y,x,y,"^y","(c) SOR, w = 1.7") % omega = 1.7


% (d) 
u_provided=zeros((Ni_x-1)*(Ni_y-1),1);
for j=1:Ni_y-1
    for i=1:Ni_x-1
        po=i+(j-1)*(Ni_x-1);
        u_provided(po) = C * cos(2*pi*(2*x(i)+3*y(j)));
    end
end

u_guess_7 = rand((Ni_x-1)*(Ni_y-1),1);
u_guess_8 = rand((Ni_x-1)*(Ni_y-1),1);
u_guess_9 = rand((Ni_x-1)*(Ni_y-1),1);

method = 3;
[solution_provided,k_provided,resarray_provided]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u_provided,resobj,u_QR,1.5);
[solution_7,k_7,resarray_7]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u_guess_7,resobj,u_QR,1.5);
[solution_8,k_8,resarray_8]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u_guess_8,resobj,u_QR,1.5);
[solution_9,k_9,resarray_9]=iterativemethods(method,A,b,x,y,Ni_x,Ni_y,u_guess_9,resobj,u_QR,1.5);

%residual is lower for the provided initial guess as the initial guess is
%close to the u_true
%but number of iterations for provided initial guess is the higher
plot_results(k_provided,resarray_provided,solution_provided,Ni_x ,Ni_y,x,y,"r","(d) provided initial guess")
plot_results(k_7,resarray_7,solution_7,Ni_x,Ni_y,x,y,"b","(d) random initial guess 1")
plot_results(k_8,resarray_8,solution_8,Ni_x,Ni_y,x,y,"g","(d) random initial guess 2")
plot_results(k_9,resarray_9,solution_9,Ni_x,Ni_y,x,y,"k","(d) random initial guess 3")

% (e)
Ny=[11:10:61];Nx = Ny+10;                        % Grid resolution and domain length
number_of_iterations = size(Ny,2)
Lx=1;Ly=1;
[result_u_Lu, result_ana] = order_accuracy(Nx,Ny,Lx,Ly,C)
