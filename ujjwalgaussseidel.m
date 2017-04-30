%UJJWALGAUSSSEIDEL
% This is the code for solving the Poisson Equation APc1-1 using
% Gauss-seidel method

% Setting of number of interior nodes in x-direction
% equal number of grids necessary

clear all; 
clc;

M=input('M=');

% Setting up number of interior points in y-direction
N= input('N=');


%% setting up increments in each direction with the paranmeters provided in
% the problem statement
%first_xincrement=M+1;
%first_yincrement=N+1;
%second_xincrement=M+2;
%second_yincrement=N+2;


%Creation of x and y values or descretization
x = linspace(-pi,pi,M+2);
y = linspace(-pi,pi,N+2);

% U is the solution of the given problem
U = ones(M+2,N+2); % The solution grid set up  
F = rightside(x,y);

%% Associated boundary conditions

% bottom bounddary condition
 ubottom = (x.*(pi-x).^2); %this is the boundary condition for y=-pi 
                       % and all x's
                       
 
% top boundary condition at y=pi and all the x's
  utop = (cos(x).*(pi-x).^2);
  
 %left hand side boundary condition
 uleft = -(4*pi^3+((y+pi)*2*pi*(pi-1)));%boundary condition evaluated at x=-pi
                                  %for all y's

  % Setting up increments along with x and y increments
  dx = (2*pi)/M;
  dy = (2*pi)/N;

 % placing the BC's on the top and bottom of the solution grid
  U(1,:)   = ubottom/dy^2;
  U(end,:) = utop/dy^2; 
 
 % Placing BC's on the left side of the solution grid
  U(:,1) = uleft/dx^2;
  
  %% Applying Neuman Condition in the right side of the Solution-Grid
  %% Applying Neuman Condition in the right side of the Solution-Grid
  % right side of the grid U will be computed using preset values
  
  for j = 2:N+1
      U(j,end) = F(j,end) - ((2*U(j,end-1))/dx^2 - (U(j-1,end))/dy^2 - (U(j+1,end))/dy^2);
  end
 
  %% Solving the grids or internal nodes of the solution vector
  
  for k = 2:M+1
    for j = 2:N+1
        U(j,k) =   F(j,k) - (U(j,k-1))/dx^2 - (U(j,k+1))/dx^2- (U(j-1,k))/dy^2 - (U(j+1,k))/dy^2 ;
    end
 end
  
