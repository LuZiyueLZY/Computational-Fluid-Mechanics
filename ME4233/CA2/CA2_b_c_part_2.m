%clear;

% specify Reynolds number
Re=18;

% specify frequency value
f = 5;

% set up the grids
Nx=50;Ny=20;   % no. of intervals e.g. the index in x-direction is 1-51 and the index in y-direction is 1 to 21
Lx=3;Ly=1;

x=linspace(0,Lx,Nx+1);dx=x(2)-x(1);
y=linspace(0,Ly,Ny+1);dy=y(2)-y(1);

[x2,y2]=meshgrid(x,y);    % generate a 2D mesh

% specify the value of dt
%dt=0.01433;   % dt_max = 0.01327 for implicit euler forward
dt = 0.002;
tf=2; % Total time
NMax=tf/dt;   % maximum of time steps. Total time = NMax*dt

% assemble A matrix to solve the Poisson equation Au=b
[A]=assembleA(Nx,Ny,dx,dy);

%%% now, this whole part is to get the vorticity and this is the b matrix

% initial condition for the vorticity
vort = zeros((Nx-1)*(Ny-1),1); % b matrix initial guess which has to be a 1D vector

t=0;

energyarray=[];

figure('units','normalized','outerposition',[0 0 1 1])


%Part c initiation


%index of u matrix for any coordinate = ((x or y) coordinate) / (dx or dy) + 1

x_index_for_point_C = 1.5/dx + 1; %West to East on diagram (1.5 units distance) i.e. top to bottom on u matrix (row number i.e. x index)
y_index_for_point_C = 0.5/dy + 1; %South to North on diagram (0.5 units distance) i.e. left to right on u matrix (column number i.e. y index)
u_value_for_point_C_array = [];

x_index_for_point_boundary = 1.5/dx + 1;
y_index_for_point_boundary = Ly/dy + 1;
u_value_for_point_boundary_array = [];

%End of Part c initiation



% time evoluion using SOR method to get stream function

for i=1:NMax
stmfunc = solve_Poisson_1D(vort,A,Nx,Ny); % solving the Poisson equation for streamfunction, implicit euler does not affect this function because this function just solves the poisson equation, we are doing the implicit method for the vorticity equation
                                     
[vort, vortnew_exact] = advance_vort_b(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,f,t); % advance the vorticity transport equation, this will get affected by the implicit euler method
%round(vort,2) == round(vortnew_exact,2)

disp(['Finish timestep' num2str(i)])

t=t+dt;

% for plotting
figure(1)
subplot(1,2,1)
[u,v]=get_uv(stmfunc,Nx,Ny,dx,dy,f,t);
quiver(x2,y2,u',v')                     % the aprostaphy means transpose
%axis equal
xlim([0 Lx]);ylim([0 Ly+0.1]);
title(['At time step = ' num2str(i)])
set(gca,'FontSize',20)
xlabel('x');ylabel('y')
pause(0.0000001)
% drawnow

energyarray = [energyarray sqrt(u(end-1,end-1)^2+v(end-1,end-1)^2)];

subplot(1,2,2)
plot(dt:dt:dt*i,energyarray,'-*b')
title(['The kinetic energy at the grid point (' num2str(x(end-1)) ',' num2str(y(end-1)) ')'])
set(gca,'FontSize',20)
xlabel('time');ylabel('energy')

%part c


u_value_for_point_C = u(x_index_for_point_C,y_index_for_point_C);
u_value_for_point_C_array = [u_value_for_point_C_array u_value_for_point_C];

figure(2)
subplot(1,1,1)
plot(dt:dt:dt*i,u_value_for_point_C_array,'-m')
title(['The u value at the grid point (' num2str(x(x_index_for_point_C)) ',' num2str(y(y_index_for_point_C)) ')'])
set(gca,'FontSize',10)
xlabel('time');ylabel('u value (velocity along x direction)')     

u_value_for_point_boundary = u(x_index_for_point_boundary,y_index_for_point_boundary);
u_value_for_point_boundary_array = [u_value_for_point_boundary_array u_value_for_point_boundary];

figure(3)
subplot(1,1,1)
plot(dt:dt:dt*i,u_value_for_point_boundary_array,'-k')
title(['The u value at the open boundary (' num2str(x(x_index_for_point_boundary)) ',' num2str(y(y_index_for_point_boundary)) ')'])
set(gca,'FontSize',10)
xlabel('time');ylabel('u value (velocity along x direction)')    

end

% post-processing
% stream function without boundary condition
stmfunc = reshape(stmfunc,Nx-1,Ny-1); %make streamfunction into 2D
plot_results(stmfunc,vort,Nx,Ny,dx,dy,x,y,x2,y2,f,t) %plot stream function with boundary condition which is 0

%% Part C (doing fourier series of the u value at c)

length_of_signal = length(u_value_for_point_C_array); %number of sample points on the TOTAL fourier series
Fs = length_of_signal / t                             %sampling frequency

%Compute the Fourier transform of the u at grid point C
Y_gridpoint_C = fft(u_value_for_point_C_array);

P2 = abs(Y_gridpoint_C/length_of_signal);
P1 = P2(1:length_of_signal/2+1);
P1(2:end-1) = 2*P1(2:end-1);
frequency_gridpoint_C_array = Fs*(0:(length_of_signal/2))/length_of_signal;

[amplitude_of_gridpoint_C,index_of_gridpoint_C] = max(P1)
frequency_gridpoint_C = frequency_gridpoint_C_array(1,index_of_gridpoint_C)

figure
plot(frequency_gridpoint_C_array,P1) 
title("Single-Sided Amplitude Spectrum of gridpoint C")
xlabel("frequency of gridpoint C (Hz)")
ylabel("|P1(f)|")

%Compute the Fourier transform of the u at grid point boundary
Y_gridpoint_boundary = fft(u_value_for_point_boundary_array);

P2_b = abs(Y_gridpoint_boundary/length_of_signal);
P1_b = P2_b(1:length_of_signal/2+1);
P1_b(2:end-1) = 2*P1_b(2:end-1);
frequency_gridpoint_boundary_array = Fs*(0:(length_of_signal/2))/length_of_signal;

[amplitude_of_gridpoint_boundary,index_of_gridpoint_boundary] = max(P1_b)
frequency_gridpoint_boundary = frequency_gridpoint_boundary_array(1,index_of_gridpoint_boundary)

figure
plot(frequency_gridpoint_boundary_array,P1_b) 
title("Single-Sided Amplitude Spectrum of boundary")
xlabel("frequency of gridpoint boundary (Hz)")
ylabel("|P1_b (f)|")

% Determine the max value and max point.

[mag_C idx_C] = max(abs(Y_gridpoint_C));
[mag_boundary idx_boundary] = max(abs(Y_gridpoint_boundary));


% determine the phase difference at the maximum point.
pY_gridpoint_C = angle(Y_gridpoint_C(idx_C));
pY_gridpoint_boundary = angle(Y_gridpoint_boundary(idx_boundary));

phase_lag = pY_gridpoint_C - pY_gridpoint_boundary %output wave - input wave % units in radians % the wave at C is lagging behind the wave at boundary

phase_lag_in_degrees = phase_lag *180/pi 