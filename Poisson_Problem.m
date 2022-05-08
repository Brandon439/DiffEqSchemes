% MT2 -- solve the Poisson problem u_{xx} + u_{yy} = 0
% on [0,5] x [0,2].  
%
% with BC's u(x,0) = 0, u(x,2) = 1
% u(0,y) = y/2, u(5,y) = 0
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.

a = 0;  
b = 5; 
c = 0;
d = 2;
m = 99;
n = 39; % chosen so that step-size is equal in both directions
h = (b-a)/(m+1);
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(c,d,n+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:n+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

rhs = zeros(m,n);          % rhs f = 0, is modified below for BC's

% set boundary conditions around edges of usoln array:
% warning: entry (1,1) corresponds to bottom left corner of domain
% turn head 90 deg clockwise to view entries of usoln as points on plot
usoln = zeros(m+2,n+2);

for j = 1:n+2
usoln(1,j) = (d-c)/(n+1)*(j-1)/2; %u(0,y) = y/2
end

usoln(:,1) = 0; %u(x,0) = 0
usoln(:,n+2) = 1; %u(x,2) = 1
usoln(m+2,:) = 0; %u(5,y) = 0


% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
rhs(:,n) = rhs(:,n) - usoln(Iint,n+2)/h^2;
rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;


% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,m*n,1);

% form matrix A:
I = speye(m);
J = speye(n);
e = ones(m,1);
f = ones(n,1);
T = spdiags([e -4*e e],[-1 0 1],m,m);
S = spdiags([f f],[-1 1],n,n);
%  kron(where to copy to make big matrix, the little matrix to copy)
A = (kron(J,T) + kron(S,I)) / h^2;

% Solve the linear system:
uvec = A\F;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usoln(Iint,Jint) = reshape(uvec,m,n);

clf
hold on
% plot solution:
contourf(X,Y,usoln,30,'k')
axis([a b c d])
daspect([1 1 1])
title('Contour plot of computed solution')
hold off
figure(2);
plot(Y, usoln(95,:),'-x')
title('Line plot of computed solution at x = 4.7')
fprintf('u(4.85,1.9) = %7.5f \n', usoln(98,39) )