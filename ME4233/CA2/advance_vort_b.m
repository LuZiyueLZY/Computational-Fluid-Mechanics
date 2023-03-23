function [vortnew,vortnew_exact] = advance_vort_b(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,f,t) 
%SOLVE_ Summary of this function goes here

%Detailed explanation goes here

%getting the RHS of the vorticity equation i.e. getting the b matrix for
%the vorticity equation i.e. getting the the sum of the (-4 * C * vorticity at timestep n) and (the boundary conditions at the edges and corners) 

C = ((dx*dx)*(dy*dy))/dt;
RHS_at_timestep_n = -4 * C * vort; %vort here is vorticity at timestep n
%RHS = assembleRHS_b(Nx,Ny,stmfunc,vort,Re,dx,dy,f,t);

[LHS,RHS_add_to_vort_at_timestep_n] = assembleLHS(Nx,Ny,stmfunc,vort,Re,dx,dy,f,t,C); % LHS is the A matrix in the page 32 week 11 notes 
RHS = RHS_at_timestep_n + RHS_add_to_vort_at_timestep_n; % RHS is the sum of the (-4 * C * vorticity at timestep n) and (the boundary conditions at the edges and corners) 

Ngs=(Nx-1)*(Ny-1);
u0=zeros(Ngs,1); %initial guess of vorticity at timestep n+1
resobj=10^-4;

A=sparse(LHS);
vortnew = SOR(A,RHS,u0,resobj);  
vortnew_exact = A\RHS; % this is the u matrix, i.e. vortecity at timestep n+1%use SOR here instead of just inversing
%vortnew = vort + dt*RHS; %formula to get vortnew using explicit euler method
                                
end









