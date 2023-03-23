function [u0]=SOR(A,b,u0,resobj)

omega=1.5;

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
    
    res=norm(b-A*u1);
    resarray=[resarray res];
    if res<resobj
        break
    end   
    disp(['Finish iteration ' num2str(k) '. res = ' num2str(res)])

    
    u0=u1;
    k=k+1;
end
end
