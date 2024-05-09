%%%%%Euler Method%%%%%% %(First Problem)%
h=0.1;
b=0.4;
a=0;
N=(b-a)/h;
L = linspace(a,b,N+1);
S = zeros(1,N+1);
S(1)=1;
for n=1:N                                          
     x=L(n);
     y=S(n);
     y=y+(h*(x/y));
     k=n+1;
     S(k)=y;
     fprintf('n= %2.0f, Solution= %12.12f\n\n',n,S(k))
end
syms u(t)
usol(t)=dsolve(diff(u,t)==(t/u),u(0)==1);
Error=abs(usol(b)-S(N+1));
fprintf('Error= %12.12f\n\n',Error);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%Euler Method%%%%%% %(2nd Problem)%
% h=0.1;
% b=1;
% a=0;
% N=(b-a)/h;
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% S(1)=1;
% for n=1:N                                          
%      x=L(n);
%      y=S(n);
%      y=y+(h*(-2*x*y^2));
%      k=n+1;
%      S(k)=y;
%      fprintf('n= %2.0f, Solution= %12.12f\n\n',n,S(k))
% end
% 
% syms u(t)
% 
% usol(t)=dsolve(diff(u,t)==-2*t*u^2,u(0)==1);
% 
% Error=abs(usol(b)-S(N+1));
% fprintf('Error= %12.12f\n\n',Error);
