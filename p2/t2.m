clear all
close all

geometry = @circleg;
hmax = 1/32;
[p,e,t] = initmesh(geometry ,'hmax',hmax);

[A,unused,b] = assema(p,t,1,1,1); % assemble
np = size(p,2); % total number of nodes

fixed = unique([e(1,:) e(2,:)]); % boundary nodes
free = setdiff([1:np],fixed); % interior nodes
b = b(free)-A(free,fixed)*myExacFun(fixed); % modify stiffness matrix
A = A(free,free); % modify load vector
xi = zeros(np,1); % allocate solution vector
xi(fixed) = myExacFun(fixed); % insert fixed node values
xi(free) = A\b; % solve for free node values