function [u0,k,resarray]=iterativemethods(method,A,b,x,y,Nx,Ny,u0,resobj,u_true, omega)
A=sparse(A);

if method ==1 & isnan(omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobi

D=diag(diag(A));
R=A-D;

k=0;

resarray=[];
uks=[];

while 1
    uks=[uks u0];
    
    u1=D\(b-R*u0);
    

    res=norm(u1 - u_true);       
    resarray=[resarray res];
        
    if res < resobj
         break
    end   
 
    disp(['Finish iteration ' num2str(k) '. res = ' num2str(res)])

    
    u0=u1;
    k=k+1;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif method ==2 & isnan(omega)
% GS

L=tril(A);

U=A-L;

k=0;
resarray=[];

uks=[];
while 1
    uks=[uks u0];
    
    u1=L\(b-U*u0);
    
    
    res=norm(u1 - u_true);
    resarray=[resarray res];
    if res<resobj
        break
    end   
    disp(['Finish iteration ' num2str(k) '. res = ' num2str(res)])

      
    u0=u1;
    k=k+1;
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif method == 3

% SOR
L=tril(A,-1);
D=diag(diag(A));
U=A-L-D;

k=0;
resarray=[];

uks=[];
while 1
    uks=[uks u0]; 
    
    u1=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0);  % SOR method
    
    res=norm(u1 - u_true);
    resarray=[resarray res];
    if res<resobj
        break
    end   
    disp(['Finish iteration ' num2str(k) '. res = ' num2str(res)])

    
    u0=u1;
    k=k+1;
end



end

