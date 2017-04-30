
% this function file creates 2D F matrix 
% x = values in x direction 
% y = values in y direction

function F = rightside(x,y)
  % values of ax, ay, bx, by as defined on the problem statement
  
  ax = -pi;
  ay = -pi;
  bx =  pi;
  by =  pi;
  
  [Y,X]=meshgrid(y,x);
  f=cos(pi/2*(((X-ax)/(bx-ax))+1)).*sin(pi*(Y-ay)/(by-ay));
  F=-f;
  
end

  
