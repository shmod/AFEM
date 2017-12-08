function z = labfun(x,y)
rho = 10;
R = 0.5;
r = 0.3;
l = sqrt(x.^2+y.^2);
for i = 1:size(x,2)
    if (abs(R-l(1,i)) <= r)
        z(i) = rho;
    else
        z(i) = 0;
    end
end