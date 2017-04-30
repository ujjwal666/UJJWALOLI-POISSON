%UJJWALGAUSSSEIDEL
% This is the code for solving the Poisson Equation APc1-1 using
% Gauss-seidel method

% Setting of number of interior nodes in x-direction
% equal number of grids necessary

M=input('M=');

% Setting up number of interior points in y-direction
N= input('N=');


%% setting up increments in each direction with the paranmeters provided in



%Creation of x and y values or descretization
x = linspace(-pi,pi,M+2);
y = linspace(-pi,pi,N+2);

% U is the solution of the given problem
U = ones(M+2,N+2); % The solution grid set up 
                   % added 2 to account for initital and final point 
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
 
  
 % Multipliers to be used while solving the equation
  
 E = 1/dx^2;
 R = 1/dy^2;
 T = -((2*E)+(2*R));
 
 % placing the BC's on the top and bottom of the solution grid
  U(1,:)   = R*ubottom;
  U(end,:) = R*utop; 
 
 % Placing BC's on the left side of the solution grid
  U(:,1) = E*uleft;
  
%% Applying Neuman Condition in the right side of the Solution-Grid
  % right side of the grid U will be computed using preset values
  err = 10; % setting up error constraint
  
  
  
  while err > 1E-5  % Setting up loop for error calculation
  B=U; % Setting up matrix for error calculation
  for j = 2:N+1
      U(j,end) = 1/T*(F(j,end) - (2*E*U(j,end-1) -R*U(j-1,end) - R*U(j+1,end)));
  end
 
  %% Solving the grids or internal nodes of the solution vector
  
   for k = 2:M+1
    for j = 2:N+1
        U(j,k) =   1/T*(F(j,k) - E*U(j,k-1) - E*U(j,k+1)- R*U(j-1,k) - R*U(j+1,k)) ;
   
     end
   end
   err = abs(max(max(((B-U)./B)))); 
  end
  
  surf(U)
  
  
