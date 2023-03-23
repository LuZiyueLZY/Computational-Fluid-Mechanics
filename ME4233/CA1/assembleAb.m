function [A,b]=assembleAb(Nx,Ny,dx,dy,x,y,C)

dx2=dx*dx;
dy2=dy*dy;
dx2dy2=dx*dx*dy*dy;

A=zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1));
b=zeros((Nx-1)*(Ny-1),1);

% interior grid points without boundary conditions
for j=2:Ny-2
    for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        % the following lines are slightly different than those in the
        % slides because we don't assume dx=dy in this code
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        A(po,po+1)=1/4;             % (i+1,j)
        A(po,po-1)=1/4;             % (i+-1,j)
        A(po,po-(Nx-1))=dx2/(9*dy2);        % (i,j-1)
        A(po,po+(Nx-1))=dx2/(9*dy2);        % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;
    end
end

% South
j=1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        A(po,po+1)=1/4;                        % (i+1,j)
        A(po,po-1)=1/4;                        % (i+-1,j)
        % A(po,po-(Nx-1))=dx2/(9*dy2);         % (i,j-1)
        A(po,po+(Nx-1))=dx2/(9*dy2);           % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;
end

% North
j=Ny-1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        A(po,po+1)=1/4;                        % (i+1,j)
        A(po,po-1)=1/4;                        % (i+-1,j)
        A(po,po-(Nx-1))=dx2/(9*dy2);           % (i,j-1)
        % A(po,po+(Nx-1))=dx2/(9*dy2);         % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;
end

% West
i=1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        A(po,po+1)=1/4;                        % (i+1,j)
        % A(po,po-1)=1/4;                      % (i+-1,j)
        A(po,po-(Nx-1))=dx2/(9*dy2);           % (i,j-1)
        A(po,po+(Nx-1))=dx2/(9*dy2);           % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;
end

% East
i=Nx-1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        % A(po,po+1)=1/4;                      % (i+1,j)
        A(po,po-1)=1/4;                        % (i+-1,j)
        A(po,po-(Nx-1))=dx2/(9*dy2);           % (i,j-1)
        A(po,po+(Nx-1))=dx2/(9*dy2);           % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;
end

% South-west
i=1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        A(po,po+1)=1/4;                        % (i+1,j)
        % A(po,po-1)=1/4;                      % (i+-1,j)
        % A(po,po-(Nx-1))=dx2/(9*dy2);           % (i,j-1)
        A(po,po+(Nx-1))=dx2/(9*dy2);         % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;

% South-east        
i=Nx-1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        % A(po,po+1)=1/4;                      % (i+1,j)
        A(po,po-1)=1/4;                        % (i+-1,j)
        % A(po,po-(Nx-1))=dx2/(9*dy2);         % (i,j-1)
        A(po,po+(Nx-1))=dx2/(9*dy2);           % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;

% North-east         
i=Nx-1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        % A(po,po+1)=1/4;                      % (i+1,j)
        A(po,po-1)=1/4;                        % (i+-1,j)
        A(po,po-(Nx-1))=dx2/(9*dy2);           % (i,j-1)
        % A(po,po+(Nx-1))=dx2/(9*dy2);         % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;

% North-west           
i=1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -((2*dx2)/(9*dy2) + 0.5);    % (i,j)
        A(po,po+1)=1/4;                        % (i+1,j)
        % A(po,po-1)=1/4;                      % (i+-1,j)
        A(po,po-(Nx-1))=dx2/(9*dy2);           % (i,j-1)
        % A(po,po+(Nx-1))=dx2/(9*dy2);         % (i,j+1)
        
        b(po)=-8 * C * pi^2 * cos(2*pi*(2*x(i)+3*y(j))) * dx2;
        
end