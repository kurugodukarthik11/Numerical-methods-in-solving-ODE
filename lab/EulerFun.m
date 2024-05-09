%%%Calling Function for Euler Method

function funcVal = EulerFun(m)

h=1/3;
b=1;
a=0;
N=(b-a)/h;


x0=0;
y0=1;
z0=m;
f=@(x,y,z) z;
g=@(x,y,z) 6*y^2-x;

fprintf('The value of slope (alpha) = %12.12f\n\n',m);

for n=1:1:N
y1=y0+h*f(x0,y0,z0);
z1=z0+h*g(x0,y0,z0);
y0=y1;
z0=z1;
x0=x0+h;
fprintf('x = %2.3f, y = %12.12f\n\n',x0,y1);
%fprintf('x = %2.3f, z = %12.12f\n\n',x0,z1);
end
funcVal = y1;
end