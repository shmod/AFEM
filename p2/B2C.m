close all
clear all
tic
L = 300;                % number of time steps
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

Mass0 = 0;
UR = U;
for K = 1:size(t2, 2);
    nodes = t2(1:3,K);
    area = polyarea(p(1,nodes), p(2,nodes));
    Mass0 = Mass0 + 1/3*sum(U(nodes))*area;
end

pdesurf(p,t2,U)
MassT = zeros(1,L);
zlabel('Concentration [mmol/mm^3]', 'fontsize', 16);

for l = 1:L
    k = t(l+1) - t(l);
    U(free) = (M+k/2*A*alph)\((M- k*alph/2*A)*U(free)+k/2*(b1+b0));    %nota that b should be zero...
    b0 = 0;
    %           pdesurf(p,t2,U)
    %           drawnow
    %           pause(T/100);
    for K = 1:size(t2, 2);
        nodes = t2(1:3,K);
        area = polyarea(p(1,nodes), p(2,nodes));
        MassT(l) = MassT(l) + 1/3*sum(U(nodes))*abs(area);
    end
end
figure
pdesurf(p,t2,U);
zlabel('Concentration [mmol/mm^3]', 'fontsize', 16)

figure
MassLoss = Mass0 - MassT;
plot(0:300, [0 MassLoss]);
ylabel('Massloss', 'fontsize', 16);
xlabel('Timestep', 'fontsize', 16);
toc

