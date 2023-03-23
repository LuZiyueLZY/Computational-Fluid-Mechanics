function plot_results(k,resarray,u0,Nx,Ny,x,y, color,name)

% plot 1 for the residue
figure('NumberTitle', 'off', 'Name', 'Residual plot for '+name)
semilogy(0:k,resarray,color)
set(gca,'FontSize',40)
xlabel('k')
ylabel('res')

solu=u0;

soluf=reshape(solu,Nx-1,Ny-1);
soluf=[0             zeros(1,Ny-1) 0
       zeros(Nx-1,1) soluf         zeros(Nx-1,1)
       0             zeros(1,Ny-1) 0 ];   
   
% % plot 2 for a 2D illustration
% figure
% contourf(x,y,soluf')
% set(gca,'FontSize',40)
% 
% % plot 3 for a 3D illustration
% figure
% surf(x,y,soluf')
% set(gca,'FontSize',40)

end

