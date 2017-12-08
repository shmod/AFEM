close all
clear all
%comment
geometry = @circleg;
error2 = zeros(5,1);
step = zeros(5,1);
for h = 1:5
    step(h) = 1/2^h;
    hmax = 1/2^h;
    [p,e,t] = initmesh(geometry ,'hmax',hmax);
    np = size(p,2);
    A = sparse(np, np);
    bk = zeros(np,1);
    r = zeros(np,1);
    I = eye(length(p));
    
    for K = 1:size(t,2);                    % loop over the triangles
        nodes = t(1:3,K);                     % find triangle K?s nodes
        % compute the (3 x 3) stiffness matrix AK
        nodeCoordinates = p(:,nodes);
        [area, b ,c] = Gradients(nodeCoordinates(1,:), nodeCoordinates(2,:));
        CoG = sum(nodeCoordinates, 2)/3;
        AK = area*[b*b' + c*c'];
        A(nodes,nodes) = A(nodes,nodes)+AK;   % add AK(i,j), i,j=1,2,3,
        % to A(nodes(i),nodes(j))
        
        %Compute bik
        bk(nodes,1) = bk(nodes,1)+ [1; 1; 1].*myFun(CoG(1), CoG(2))*area/3;
        
    end
    
    %Compute rik by iterating through the edges
%     for E = 1:size(e,2)
%         nodes = e(1:2, E);
%         x = p(1, nodes);
%         y = p(2, nodes);
%         length_E = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
%         r(nodes,1) = r(nodes,1) + length_E*myExacFun(x,y)'/2;
%     end

    A(e(1,:),:) = I(e(1,:),:);
    bk(e(1,:)) = 0;
    
    fixed = unique([e(1,:) e(2,:)]);            % boundary nodes

    g = myExacFun(p(1, fixed),p(2, fixed));
    free = setdiff([1:np],fixed);               % interior nodes
    bk = bk(free)-A(free,fixed)*g';             % modify load vector
    A = A(free,free);                           % modify stiffness matrix
    Z = zeros(np,1);                            % allocate solution vector
    Z(fixed) = g';                              % insert fixed node values
    Z(free) = A\bk;                             % solve for free node values
    
    %Error
    exactSol = myExacFun(p(1, free),p(2, free));
    error1 = exactSol'-Z(free);
    error2(h) = sqrt(error1'*A*error1);
    if (h == 1 | h == 5)
         pdesurf(p,t,Z)
         xlabel('x','fontsize',16)
         ylabel('y','fontsize',16)
         zlabel('u_h','fontsize',16)
         figure
    end
    %pdemesh(p,e,t) 
    %     figure
end

loglog(step,error2,'*')
gamma = abs((error2(end)-error2(1))/(1/2-1/32));
hold on
loglog(step,step.^gamma);
xlabel('h_{max}','fontsize',16)
ylabel('Energy norm of the error','fontsize',16)
legend('Energynorm', 'h_{max}^{gamma}')