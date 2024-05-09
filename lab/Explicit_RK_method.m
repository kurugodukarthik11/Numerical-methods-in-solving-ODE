% %%%%%Modified Euler Method%%%%%% 
% h=0.2;
% b=0.4;
% a=0;
% N=(b-a)/h;
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% S(1)=1;
% for n=1:N                                          
%      x=L(n);
%      y=S(n);
%      K1=-2*x*y^2;
%      K2=-2*(x+h/2)*(y+(h/2)*K1)^2;
%      y=y+(h*K2);
%      k=n+1;
%      S(k)=y;
%      fprintf('n= %2.0f, Solution= %12.12f\n\n',n,S(k))
% end
% 
% syms u(t)
% 
% usol(t)=dsolve(diff(u,t)==-2*(t*u^2),u(0)==1);
% 
% Error=abs(usol(b)-S(N+1));
% fprintf('Error= %12.12f\n\n',Error);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%Heun's Method%%%%%% 
% h=0.2;
% b=0.4;
% a=0;
% N=(b-a)/h;
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% S(1)=1;
% for n=1:N                                          
%      x=L(n);
%      y=S(n);
%      K1=-2*x*y^2;
%      K2=-2*(x+h)*(y+h*K1)^2;
%      y=y+(h/2)*(K1+K2);
%      k=n+1;
%      S(k)=y;
%      fprintf('n= %2.0f, Solution= %12.12f\n\n',n,S(k))
% end
% 
% syms u(t)
% 
% usol(t)=dsolve(diff(u,t)==-2*(t*u^2),u(0)==1);
% 
% Error=abs(usol(b)-S(N+1));
% fprintf('Error= %12.12f\n\n',Error);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%RK Method of Fourth order%%%%%% 
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
     K1=-2*x*y^2;
     K2=-2*(x+h/2)*(y+(h/2)*K1)^2;
     K3=-2*(x+h/2)*(y+(h/2)*K2)^2; 
     K4=-2*(x+h)*(y+h*K3)^2;
     y=y+(h/6)*(K1+2*K2+2*K3+K4);
     k=n+1;
     S(k)=y;
     fprintf('n= %2.0f, Solution= %12.12f\n\n',n,S(k))
end
syms u(t)
usol(t)=dsolve(diff(u,t)==-2*(t*u^2),u(0)==1);
fprintf('Exact solution= %12.12f\n\n',usol(b));
Error=abs(usol(b)-S(N+1));
fprintf('Error= %12.12f\n\n',Error);




