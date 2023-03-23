function [LHS,RHS_add_to_vort_at_timestep_n] = assembleLHS(Nx,Ny,stmfunc,vort,Re,dx,dy,f,t,C)
%ASSEMBLERHS Summary of this function goes here
%   Detailed explanation goes here

% the boundary conditions at the boundaries remains the same 
U_south = zeros(Nx-1,1);
U_north = ones(Nx-1,1)*sin(2*pi*f*t); 
U_west  = zeros(Ny-1,1);
U_east  = zeros(Ny-1,1);

nu=4/Re;
dx2=dx*dx;
dy2=dy*dy;

LHS= zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1));
RHS_add_to_vort_at_timestep_n = zeros((Nx-1)*(Ny-1),1);
stmfunc = reshape(stmfunc,Nx-1,Ny-1); % make streamfunction into 2D

fac1 =  nu*-2*((dx2)+(dy2))-4*C;
% interior
for j=2:Ny-2
    for i=2:Nx-2

        fac2 =  (stmfunc(i,j-1) - stmfunc(i,j+1))*dx*dy + nu*dy2;
        fac3 =  (stmfunc(i,j+1) - stmfunc(i,j-1))*dx*dy + nu*dy2;
        fac4 =  (stmfunc(i-1,j) - stmfunc(i+1,j))*dx*dy + nu*dx2;
        fac5 =  (stmfunc(i+1,j) - stmfunc(i-1,j))*dx*dy + nu*dx2;    

        po=i+(j-1)*(Nx-1);
        LHS(po,po)= fac1;                           % (i,j)
        LHS(po,po+1)=fac2;                          % (i+1,j)
        LHS(po,po-1)=fac3;                          % (i+-1,j)
        LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
        LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)

                       
    end
end

% south
j=1;

for i=2:Nx-2
    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2); %remains the same as page 16 week 11 as taylor series expansion of vorticity is same
        
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (0              - stmfunc(i,j+1))*dx*dy + nu*dy2; % stream function (i, j-1) = 0
    fac3 =  (stmfunc(i,j+1) - 0             )*dx*dy + nu*dy2; % stream function (i, j-1) = 0
    fac4 =  (stmfunc(i-1,j) - stmfunc(i+1,j))*dx*dy + nu*dx2;
    fac5 =  (stmfunc(i+1,j) - stmfunc(i-1,j))*dx*dy + nu*dx2;    

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    LHS(po,po+1)=fac2;                          % (i+1,j)
    LHS(po,po-1)=fac3;                          % (i+-1,j)
    %LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)

    RHS_add_to_vort_at_timestep_n(po)= vortsouthbc*-fac4;     
end

% North
j=Ny-1; % j+1
for i=2:Nx-2
    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2); %remains the same as page 17 week 11 as taylor series expansion of vorticity is same
        
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (stmfunc(i,j-1) - 0             )*dx*dy + nu*dy2; % stream function (i, j+1) = 0
    fac3 =  (0              - stmfunc(i,j-1))*dx*dy + nu*dy2; % stream function (i, j+1) = 0
    fac4 =  (stmfunc(i-1,j) - stmfunc(i+1,j))*dx*dy + nu*dx2;
    fac5 =  (stmfunc(i+1,j) - stmfunc(i-1,j))*dx*dy + nu*dx2;    

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    LHS(po,po+1)=fac2;                          % (i+1,j)
    LHS(po,po-1)=fac3;                          % (i+-1,j)
    LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    %LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)
    
    RHS_add_to_vort_at_timestep_n(po)= vortnorthbc*-fac5; 
end

% West
i=1; % i-1
for j=2:Ny-2
    vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx)/(0.5*dx2); %remains the same as page 18 week 11 as taylor series expansion of vorticity is same
    
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (stmfunc(i,j-1) - stmfunc(i,j+1))*dx*dy + nu*dy2;
    fac3 =  (stmfunc(i,j+1) - stmfunc(i,j-1))*dx*dy + nu*dy2;
    fac4 =  (0              - stmfunc(i+1,j))*dx*dy + nu*dx2; % stream function (i-1, j) = 0
    fac5 =  (stmfunc(i+1,j) - 0             )*dx*dy + nu*dx2; % stream function (i-1, j) = 0

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    LHS(po,po+1)=fac2;                          % (i+1,j)
    %LHS(po,po-1)=fac3;                          % (i+-1,j)
    LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)
    
    RHS_add_to_vort_at_timestep_n(po)= vortwestbc*-fac3; 
end

% East
i=Nx-1; % i+1
for j=2:Ny-2
    vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx)/(0.5*dx2); %remains the same as page 19 week 11 as taylor series expansion of vorticity is same
    
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (stmfunc(i,j-1) - stmfunc(i,j+1))*dx*dy + nu*dy2;
    fac3 =  (stmfunc(i,j+1) - stmfunc(i,j-1))*dx*dy + nu*dy2;
    fac4 =  (stmfunc(i-1,j) - 0             )*dx*dy + nu*dx2; % stream function (i+1, j) = 0
    fac5 =  (0              - stmfunc(i-1,j))*dx*dy + nu*dx2; % stream function (i+1, j) = 0

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    %LHS(po,po+1)=fac2;                          % (i+1,j)
    LHS(po,po-1)=fac3;                          % (i+-1,j)
    LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)
    
    RHS_add_to_vort_at_timestep_n(po)= vorteastbc*-fac2; 
end

% South-west
i=1;j=1;

    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2); %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
    vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx)/(0.5*dx2);   %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
    
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (0              - stmfunc(i,j+1))*dx*dy + nu*dy2; % stream function (i, j-1) = 0
    fac3 =  (stmfunc(i,j+1) - 0             )*dx*dy + nu*dy2; % stream function (i, j-1) = 0
    fac4 =  (0              - stmfunc(i+1,j))*dx*dy + nu*dx2; % stream function (i-1, j) = 0
    fac5 =  (stmfunc(i+1,j) - 0             )*dx*dy + nu*dx2; % stream function (i-1, j) = 0

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    LHS(po,po+1)=fac2;                          % (i+1,j)
    %LHS(po,po-1)=fac3;                          % (i+-1,j)
    %LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)
    
    RHS_add_to_vort_at_timestep_n(po)= vortsouthbc*-fac4 + vortwestbc*-fac3 ;  
                    
% South-east          
i=Nx-1;j=1;

    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2); %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
    vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx)/(0.5*dx2);   %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
    
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (0              - stmfunc(i,j+1))*dx*dy + nu*dy2; % stream function (i, j-1) = 0
    fac3 =  (stmfunc(i,j+1) - 0             )*dx*dy + nu*dy2; % stream function (i, j-1) = 0
    fac4 =  (stmfunc(i-1,j) - 0             )*dx*dy + nu*dx2; % stream function (i+1, j) = 0
    fac5 =  (0              - stmfunc(i-1,j))*dx*dy + nu*dx2; % stream function (i+1, j) = 0 

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    %LHS(po,po+1)=fac2;                          % (i+1,j)
    LHS(po,po-1)=fac3;                          % (i+-1,j)
    %LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)

    RHS_add_to_vort_at_timestep_n(po)= vortsouthbc*-fac4 + vorteastbc*-fac2 ;  
                    
% North-east         
i=Nx-1;j=Ny-1;

    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2); %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
    vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx)/(0.5*dx2);   %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
        
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (stmfunc(i,j-1) - 0             )*dx*dy + nu*dy2; % stream function (i, j+1) = 0
    fac3 =  (0              - stmfunc(i,j-1))*dx*dy + nu*dy2; % stream function (i, j+1) = 0
    fac4 =  (stmfunc(i-1,j) - 0             )*dx*dy + nu*dx2; % stream function (i+1, j) = 0
    fac5 =  (0              - stmfunc(i-1,j))*dx*dy + nu*dx2; % stream function (i+1, j) = 0 

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    %LHS(po,po+1)=fac2;                          % (i+1,j)
    LHS(po,po-1)=fac3;                          % (i+-1,j)
    LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    %LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)

    RHS_add_to_vort_at_timestep_n(po)= vortnorthbc*-fac5 + vorteastbc*-fac2 ;  

% North-west          
i=1;j=Ny-1;

    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2); %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
    vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx)/(0.5*dx2);   %remains the same as page (nil) week 11 as taylor series expansion of vorticity is same
        
    %fac1 =  nu*-2*((dx2)+(dy2))+4*C;
    fac2 =  (stmfunc(i,j-1) - 0             )*dx*dy + nu*dy2; % stream function (i, j+1) = 0
    fac3 =  (0              - stmfunc(i,j-1))*dx*dy + nu*dy2; % stream function (i, j+1) = 0
    fac4 =  (0              - stmfunc(i+1,j))*dx*dy + nu*dx2; % stream function (i-1, j) = 0
    fac5 =  (stmfunc(i+1,j) - 0             )*dx*dy + nu*dx2; % stream function (i-1, j) = 0   

    po=i+(j-1)*(Nx-1);
    LHS(po,po)= fac1;                           % (i,j)
    LHS(po,po+1)=fac2;                          % (i+1,j)
    %LHS(po,po-1)=fac3;                          % (i+-1,j)
    LHS(po,po-(Nx-1))=fac4;                     % (i,j-1)
    %LHS(po,po+(Nx-1))=fac5;                     % (i,j+1)
          
    RHS_add_to_vort_at_timestep_n(po)= vortnorthbc*-fac5 + vortwestbc*-fac3 ;  
end

