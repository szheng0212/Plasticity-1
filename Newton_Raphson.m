x0=0.5;
x_m0=0.5;
x1=1;
phi=0;
count=0;
while abs(x1-x_m0)>1.0e-3
    x1=x0-(3*x0*x0-3)/(6*x0);
    x_m0=x0;
    x0=x1;
    count=count+1;
end
