function [u,v] = get_uv(stmfunc,Nx,Ny,dx,dy,f,t)
%GET_UV Summary of this function goes here
%   Detailed explanation goes here

stmfunc = reshape(stmfunc,Nx-1,Ny-1); %make streamfunction into 2D
stmfunc=[0             zeros(1,Ny-1)   0
         zeros(Nx-1,1) stmfunc         zeros(Nx-1,1)
         0             zeros(1,Ny-1)   0 ]; 
     
u =  (stmfunc(2:end-1,3:end) - stmfunc(2:end-1,1:end-2))/2/dy;  %interior grid points streamfunctionat 2nd index to streanfunction at end-1 index and using central difference scheme+formula relating streamfunction and velocity, 
v = -(stmfunc(3:end,2:end-1) - stmfunc(1:end-2,2:end-1))/2/dx;

u=[0             zeros(1,Ny-1)  0
   zeros(Nx-1,1) u              ones(Nx-1,1)*sin(2*pi*f*t)
   0             zeros(1,Ny-1)  0 ]; 

v=[0             zeros(1,Ny-1)  0
   zeros(Nx-1,1) v              zeros(Nx-1,1)
   0             zeros(1,Ny-1)  0 ]; 

end

