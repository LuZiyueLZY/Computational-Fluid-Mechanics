function plot_results_solution_only(u0,Nx,Ny,x,y,name)
solu=u0;

soluf=reshape(solu,Nx-1,Ny-1);
soluf=[0             zeros(1,Ny-1) 0
       zeros(Nx-1,1) soluf         zeros(Nx-1,1)
       0             zeros(1,Ny-1) 0 ];   
   
% plot 2 for a 2D illustration
figure('NumberTitle', 'off', 'Name', '2D plot for '+name);
contourf(x,y,soluf')
set(gca,'FontSize',40)

% plot 3 for a 3D illustration
figure('NumberTitle', 'off', 'Name', '3D plot for '+name)
title("3D plot for "+name)
surf(x,y,soluf')
set(gca,'FontSize',40)

end

