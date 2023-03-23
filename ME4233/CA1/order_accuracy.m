function [result_u_Lu, result_ana]=order_accuracy(xarray,yarray,Lx,Ly,C)

% initialisation

error = zeros(length(xarray),1);
dxarray= zeros(length(xarray),1);
dyarray= zeros(length(yarray),1);
harray= zeros(length(xarray),1);

result_u_Lu = []
result_ana = []

for iM=1:length(xarray)
    M = xarray(iM);
    P = yarray(iM);

    x=linspace(0,Lx,M+1);dx=x(2)-x(1);   % Set up the grids
    y=linspace(0,Ly,P+1);dy=y(2)-y(1);
    
    [A,b]=assembleAb_e(M,P,dx,dy,x,y,C);  % Assemble the A matrix and b vector
                                            % the boundary condition here is smoother. Note the implementaiton of the
                                            % boundary condition in this subfunction

    %%LU decomposition
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


    % get the analytical solution and calculate the error
    ana=zeros((M-1)*(P-1),1);
    for j=1:P-2
        for i=1:M-1
            po=i+(j-1)*(M-1);
            ana(po) = C*cos(2*pi*(2*x(i)+3*y(j)));
        end
    end
    
    error(iM) = sqrt(sum(( u_Lu - ana ).^2)/(M*N)) % RMS
    size(u_Lu)
    size(ana)
    dxarray(iM) = dx;
    dyarray(iM) = dy;
    harray(iM) = sqrt(dx * dy)
    result_u_Lu{end+1} = [u_Lu]
    result_ana{end+1} = [ana]

end

% plot for verifying the order of accuracy
set(0,'DefaultFigureWindowStyle','normal') 
figure
loglog(harray,error,'-*b')
hold on
loglog(harray, 2*harray.^2,'-r')
legend('Your result','f = c dx^2','Location','northwest')
xlabel('h');ylabel('error')
set(gca,'FontSize',10)