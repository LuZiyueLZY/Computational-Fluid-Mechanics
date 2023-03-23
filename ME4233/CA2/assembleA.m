function [A]=assembleA(Nx,Ny,dx,dy) % assemble A matrix in the poisson equation

dx2=dx*dx;
dy2=dy*dy;

A=zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1)); % if there are 51 points on the grid in x direction, there are 50 intervals, so, there are 49 interior points. Now, this Nx is number of intervals

% interior
for j=2:Ny-2
    for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        A(po,po+1)=1/(dx2);                          % (i+1,j)
        A(po,po-1)=1/(dx2);                          % (i+-1,j)
        A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)
    end
end

% South
j=1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        A(po,po+1)=1/(dx2);                          % (i+1,j)
        A(po,po-1)=1/(dx2);                          % (i+-1,j)
        %A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)
end

% North
j=Ny-1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        A(po,po+1)=1/(dx2);                          % (i+1,j)
        A(po,po-1)=1/(dx2);                          % (i+-1,j)
        A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        %A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)
end

% West
i=1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        A(po,po+1)=1/(dx2);                          % (i+1,j)
        %A(po,po-1)=1/(dx2);                          % (i+-1,j)
        A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)
end

% East
i=Nx-1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        %A(po,po+1)=1/(dx2);                          % (i+1,j)
        A(po,po-1)=1/(dx2);                          % (i+-1,j)
        A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)
end

% South-west
i=1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        A(po,po+1)=1/(dx2);                          % (i+1,j)
        %A(po,po-1)=1/(dx2);                          % (i+-1,j)
        %A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)

% South-east          
i=Nx-1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        %A(po,po+1)=1/(dx2);                          % (i+1,j)
        A(po,po-1)=1/(dx2);                          % (i+-1,j)
        %A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)

% North-east         
i=Nx-1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        %A(po,po+1)=1/(dx2);                          % (i+1,j)
        A(po,po-1)=1/(dx2);                          % (i+-1,j)
        A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        %A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)

% North-west          
i=1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*((1/dx2)+(1/dy2));            % (i,j)
        A(po,po+1)=1/(dx2);                          % (i+1,j)
        %A(po,po-1)=1/(dx2);                          % (i+-1,j)
        A(po,po-(Nx-1))=1/(dy2);                     % (i,j-1)
        %A(po,po+(Nx-1))=1/(dy2);                     % (i,j+1)

end