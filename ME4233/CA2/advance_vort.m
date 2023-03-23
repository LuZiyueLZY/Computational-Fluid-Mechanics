function [vortnew] = advance_vort(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,f,t)
%SOLVE_ Summary of this function goes here
%   Detailed explanation goes here

RHS = assembleRHS(Nx,Ny,stmfunc,vort,Re,dx,dy,f,t);

vortnew = vort + dt*RHS; %formula to get vortnew using explicit euler method
                                
end


