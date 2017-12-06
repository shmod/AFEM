close all
clear all

L = 300;                % number of time levels
T = 30;                  % final time
t = linspace(0,T,L+1);   % time grid
h = 5;
alph = 0.01;

geometry = @circleg;
hmax = 1/2^h;
[p,e,t2] = initmesh(geometry ,'hmax',hmax);
U = labfun(p(1,:), p(2,:))';          % inital condition

A = StiffMat2D(p,t2,1);
M = MassMat2D(p, t2);
b0 = LoadVec2D(p, t2, @labfun);
b1 = LoadVec2D(p, t2, @labfun2);
np = size(p,2);

fixed = unique([e(1,:) e(2,:)]);            % boundary nodes
free = setdiff([1:np],fixed);               % interior nodes
g = zeros(size(fixed,2),1);

b0 = b0(free);             % modify load vector
b1 = b1(free);
A = A(free,free);
M = M(free,free);
U(fixed) = 0;



for l = 1:L
    1
    k = t(l+1) - t(l);
    U(free) = (M+k/2*A*alph)\((M- k*alph/2*A)*U(free)+k/2*(b1+b0));    %nota that b should be zero...
    b0 = 0;
        pdesurf(p,t2,U)
        drawnow
        pause(T/100);
end



