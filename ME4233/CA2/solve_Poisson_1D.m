function [stmfunc] = solve_Poisson_1D(vort,A,Nx,Ny)
%SOLVE_POISSON Summary of this function goes here
%   Detailed explanation goes here

% RHS of the Poisson equation. Note that the bc for streamfunction is zero.
b = -vort; % vort is already a 1D vector, this is the b matrix for that timestep (the b matrix in the for loop)

Ngs=(Nx-1)*(Ny-1);
u0=zeros(Ngs,1); %initial guess of stream function
resobj=10^-4;

A=sparse(A);
stmfunc=SOR(A,b,u0,resobj);

% When we solve the Poisson equation, we stack the 2D grid points into a 1D
% vector. But when we solve the time-dependent equation, we can work with a
% 2D mesh. So here, we reshape the 1D vector back to the 2D mesh

end

