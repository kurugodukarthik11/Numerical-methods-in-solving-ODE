%%%%%Second order Implicit RK Method%%%%%%
h=0.2;
b=0.4;
a=0;
N=(b-a)/h;
L = linspace(a,b,N+1);
S = zeros(1,N+1);
S(1)=1;
for n=1:N                                          
     x=L(n);
     y=S(n);       
 eps = 1; tol = 10^(-3); total = 100; j = 0; format long; 
 z=y;
 while ((eps > tol)&&(j < total))
   f = z+2*(x+h/2)*(y+(h*z)/2)^2;
   f1 = 1+2*h*(x+h/2)*(y+(h*z)/2);
   zz = z-f/f1;
   eps = abs(zz-z); z = zz;
   fprintf('j= %2.0f, The root = %12.12f\n\n',j,z);
   j = j+1;
 end
     y=y+h*z; 
     k=n+1;
     S(k)=y;
     fprintf('n= %2.0f, Solution= %12.12f\n\n',n,S(k));
end
syms u(t)
usol(t)=dsolve(diff(u,t)==-2*(t*u^2),u(0)==1);
fprintf('Exact solution= %12.12f\n\n',usol(b));
Error=abs(usol(b)-S(N+1));
fprintf('Error= %12.12f\n\n',Error);