% a = [1,2;3,4;5,6];
% b = [1,1;2,2];
% disp(a*b);
% % disp(a.*b);
% disp(length(a));
% disp(size(a));
% x = 20;
% if x==20
%     disp(x);
% end
% sum = 0;
% for i=0:1:10
%     if rem(i,2) == 0
%         disp(i)
%         sum = sum + i;
%     end
% end
% disp(sum)

% x = 1;
% f = @(x) x.^2 + 5*x;
% I = integral(f,0,1);
% disp(I)

% f = @(x,y) 1./(sqrt(x+y) .* (1+x+y).^2);
% % ymax = @(x) 1-x;
% I2 = integral2(f,0,1,0,@(x) 1-x);
% disp(I2)

% f = @(r,t,p) r.^2.*sin(t).*sin(p);
% I3 = integral3(f,0,pi,0,pi,0,@(r) r);
% I4 = integral()

% ode = diff(y,t)==t*y^2;
% sol(t) = dsolve(ode,y(0)==1);
% disp(sol(t));

% ode = diff(y,x,2)==exp(2*x)-2*y-x;
% Dy = diff(y);
% sol = dsolve(ode,y(0)==1,Dy(0)==0);
% disp(sol);

%%%Implicit
% x0 = 3;
% f(x)= x^3 - 5*x + 1;
% Df = diff(f,x);
% % disp(Df(x0))
% x1 = 1;
% tol = 2;
% while tol > 10^(-14)
%     x1 = x0 - f(x0)/Df(x0);
%     tol = abs(x0 - x1);
%     x0 = x1;
%     fprintf("%12.12f\n",x0);
% end
% fprintf("Root of the eqn: %12.12f\n",x0);

% syms y(x);
% ode = diff(y,x) ==x/y;
% cond = y(0)==1;
% sol(x) = dsolve(ode,cond);
% z = sol(0.4);
% fprintf('\nThe value of y at x=0.4 is y=%2.4f',z);

% syms y(x);
% ode = diff(y,x)==x/y;
% sol(x) = dsolve(ode,y(0)==1);
% z = sol(0.4);
% fprintf("%2.4f\n", z);

% h1 = 0.2;
% x0=0;
% y0=1;
% f = @(x,y) x/y;
% while x0<0.4
%     y1 = y0+ h1*f(x0,y0);
%     x0 = x0 + h1;
%     y0=y1; 
% end
% fprintf("%2.6f\n", y0);

% h2 = 0.1;
% x0=0;
% y0=1;
% f = @(x,y) x/y;
% while x0<0.4
%     % fprintf("%.3f\n",f(x0,y0))
%     y1 = y0+ h2*f(x0,y0);
%     x0 = x0 + h2;
%     y0=y1; 
% end
% fprintf("%2.6f", y0);

% syms y(x);
% 
% ode = diff(y,x) == -2*x*y^2;
% sol(x) = dsolve(ode,y(0)==1);
% fprintf("%2.4f %2.4f\n",sol(0.2), sol(0.4));

% %Backward Euler method
% a=0;
% b=0.4;
% h=0.2;
% N= (b-a)/h;
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% S(1)=1;
% F = @(xn1,yn1,yn,h) yn1-yn+2*h*xn1*(yn1).^2;
% F1 = @(xn1,yn1,h) 1+4*h*xn1*yn1;
% tol=10^-5;
% for n=1:N
%     x=L(n+1);
%     y=S(n);
%     eps=1;
%     z=y;
%     while(eps>tol)
%         zz = z - (F(x,z,y,h))/(F1(x,z,h));
%         eps = abs(z-zz);
%         z = zz;
%         fprintf("%d %12.12f\n",n,z);
%     end
%     k=n+1;
%     S(k)=z;
%     fprintf("%12.12f\n",S(k));
% end


% % h = 0.2;
% % x0=0;
% % y0=1;
% % f = @(x,y) -2*x*(y).^2;
% % % S = zeros(1,4);
% % % S(1) = 1;
% % i=1;
% % while x0<0.4
% %     % S(i) = y0;
% %     k1 = f(x0,y0);
% %     k2 = f(x0 + h/2, y0+((h/2)*k1));
% %     k3 = f(x0 + h/2, y0+((h/2)*k2));
% %     k4 = f(x0+h,y0+(h*k3));
% %     y1 = y0 + (h/6 *(k1 + 2*k2 + 2*k3 + *k4));
% %     x0 = x0 + h;
% %     y0=y1;
% %     % i = i+1;
% % end
% % fprintf("%2.6f\n", y0);
% % disp(S);

% RK 4th order
% a=0;
% b=0.4;
% h=0.2;
% N= (b-a)/h;
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% S(1)=1;
% f = @(x,y) -2*x.*(y).^2;
% for n=1:N
%     x0 = L(n);
%     y0= S(n);
%     k1 = f(x0,y0);
%     k2 = f(x0 + h/2, y0+((h/2)*k1));
%     k3 = f(x0 + h/2, y0+((h/2)*k2));
%     k4 = f(x0+h,y0+(h*k3));
%     y1 = y0 + (h/6 *(k1 + 2*k2 + 2*k3 + k4));
%     k=n+1;
%     S(k) = y1;
%     fprintf("n=%2.0f, Solution = %12.12f\n\n",n,S(k));
% end
% fprintf("Error : %12.12f", abs(sol(0.4) - S(N+1)));

% Implicit RK 2nd order
% a=0;
% b=0.4;
% h=0.2;
% N = (b-a)/h;
% L = linspace(a,b,N+1);
% S1 = zeros(1,N+1);
% S1(1)=1;
% f = @(x,y) -2*x.*(y).^2;
% for n=1:N
%     x0 = L(n);
%     y0= S1(n);
%     tol = 2;
%     % j=0;
%     % total = 100;
%     z=y0;
%     while (tol>(10^(-3)))
%         F = z - f(x0+(h/2),(y0+(h*z)/2));
%         % F = z+2*(x0+h/2)*(y0+(h*z)/2)^2;
%         F1 = 1+2*h*(x0+h/2)*(y0 + (h*z)/2);
%         zz = z - F/F1;
%         tol = abs(zz-z);
%         z = zz;
%         % j=j+1;
%     end
%     % fprintf("z: %f\n",z);
%     y1 = y0 + h*z;
%     k=n+1;
%     S1(k) = y1;
%     fprintf("n=%2.0f, Solution = %12.12f\n\n",n,S1(k));
% end
% fprintf("Error : %12.12f", abs(sol(0.4) - S1(N+1)));

% syms x y(x);
% ode = diff(y,x) == x.^2 + y.^2;
% sol(x) = dsolve(ode,y(0)==1);
% ex(1) = sol(0.3);
% ex(2) = sol(0.4);
% fprintf("%2.4f\n", ex(1));
% fprintf("%2.4f\n", ex(2));
% 
% f1 = @(x,y) x.^2 + y.^2;
% % f1 = x^2+y^2;
% % disp(f1(1,2));
% % syms(f1,x)
% % f2 = diff(f1,x);
% f2 = @(x,y) 2*x + 2*y*f1(x,y);
% % disp(f2);
% % f3 = diff(f2,x);
% f3 = @(x,y) 2+2*(y*f2(x,y)+(f1(x,y).^2));
% 
% Adam Bashforth of order 3 using power series
% a=0;
% b=0.2;
% h=0.1;
% N = (b-a)/h;
% L = linspace(a,b,N+1);
% S1 = zeros(1,N+1);
% S1(1)=1;
% for n=1:N
%     x0 = L(n);
%     y0 = S1(n);
%     y1 = y0 + h*f1(x0,y0) + (h.^2)*f2(x0,y0)/2 + (h.^3)*f3(x0,y0)/6;
%     k=n+1;
%     S1(k) = y1;
% end
% S1(N+2) = S1(N+1) + h*(23*f1(L(N+1),S1(N+1)) - 16*f1(L(N),S1(N)) + 5*f1(L(N-1),S1(N-1)))/12;
% fprintf("%12.12f\n",S1(N+2));
% 
% Adam Bashforth of order 4 using euler method
% a=0;
% b=0.3;
% h=0.1;
% N = int16((b-a)/h);
% L = linspace(a,b,N+1);
% S1 = zeros(1,N+1);
% S1(1)=1;
% for n=1:N
%     x0 = L(n);
%     y0 = S1(n);
%     y1 = y0 + h*f1(x0,y0);
%     k=n+1;
%     S1(k) = y1;
% end
% % disp(S1);
% % disp(L);
% S1(N+2) = S1(N+1) + h*(55*f1(L(N+1),S1(N+1)) - 59*f1(L(N),S1(N)) + 37*f1(L(N-1),S1(N-1)) - 9*f1(L(N-2),S1(N-2)))/24;
% fprintf("%12.12f",S1(N+2));


% a=0;
% b=0.5;
% h=0.1;
% N = (b-a)/h;
% L = linspace(a,b,N+1);
% S1 = zeros(1,N+1);
% S1(1)=1;
% f1 = @(x,y) x+y;
% f2 = @(x,y) 1+f1(x,y);
% f3 = @(x,y) f2(x,y);
% m=3;
% for j=m:N
% for n=1:j-1
%     x0 = L(n);
%     y0 = S1(n);
%     y1 = y0 + h*f1(x0,y0) + (h.^2)*f2(x0,y0)/2 + (h.^3)*f3(x0,y0)/6;
%     k=n+1;
%     S1(k) = y1;
% end
% S1(j+1) = S1(j) + h*(23*f1(L(j),S1(j)) - 16*f1(L(j-1),S1(j-1)) + 5*f1(L(j-2),S1(j-2)))/12;
% end
% fprintf("%12.12f\n", S1(6))

% syms y(x);
% ode = diff(y,x) == x+y;
% sol(x) = dsolve(ode,y(0)==1);
% ex(1) = sol(0.5);
% fprintf("%2.4f\n", ex(1));
% 
% fprintf("%12.12f", abs(ex(1)-S1(6)));

% a=0;
% b=0.4;
% h=0.2;
% N=(b-a)/h;
% L = linspace(a,b,N+1);
% S1 = zeros(1,N+1);
% S2 = zeros(1,N+1);
% S1(1) = 1;
% S2(1) = -1;
% f = @(x,y,z) x+(y*z);
% g = @(x,y,z) y+(x*z);
% for n=1:N
%     x0 = L(n);
%     y0 = S1(n);
%     z0 = S2(n);
%     y1 = y0 + h*f(x0,y0,z0);
%     z1 = z0 + h*g(x0,y0,z0);
%     k = n+1;
%     S1(k) = y1;
%     S2(k) = z1;
% end
% disp(S1);
% disp(S2);

% Higher order IVP using RK 4th order
% a=0;
% b=0.4;
% h=0.2;
% N=(b-a)/h;
% L = linspace(a,b,N+1);
% S1 = zeros(1,N+1);
% S2 = zeros(1,N+1);
% S1(1) = 1;
% S2(1) = 0;
% K1 = @(x,y,z) z;
% L1 = @(x,y,z) cos(x) - 4*y;
% K2 = @(x,y,z) K1(x+(h/2),y+(h/2)*K1(x,y,z),z+(h/2)*L1(x,y,z));
% L2 = @(x,y,z) L1(x+(h/2),y+(h/2)*K1(x,y,z),z+(h/2)*L1(x,y,z));
% K3 = @(x,y,z) K1(x+(h/2),y+(h/2)*K2(x,y,z),z+(h/2)*L2(x,y,z));
% L3 = @(x,y,z) L1(x+(h/2),y+(h/2)*K2(x,y,z),z+(h/2)*L2(x,y,z));
% K4 = @(x,y,z) K1(x+h,y+h*K3(x,y,z),z+h*L3(x,y,z));
% L4 = @(x,y,z) L1(x+h,y+h*K3(x,y,z),z+h*L3(x,y,z));
% for n=1:N
%     x0 = L(n);
%     y0 = S1(n);
%     z0 = S2(n);
%     y1 = y0 + (h/6)*(K1(x0,y0,z0) + 2*K2(x0,y0,z0) + 2*K3(x0,y0,z0) + K4(x0,y0,z0));
%     z1 = z0 + (h/6)*(L1(x0,y0,z0) + 2*L2(x0,y0,z0) + 2*L3(x0,y0,z0) + L4(x0,y0,z0));
%     k = n+1;
%     S1(k) = y1;
%     S2(k) = z1;
% end
% fprintf("Numerical: %12.12f\n",S1(3));
% ex = (2*cos(2*0.4) + cos(0.4))/3;
% fprintf("Exact: %12.12f\n",ex);

% a=0;
% b=1;
% h=1/3;
% g = @(x,y) 6*y^2 - x;
% f = @(z) z;
% N = (b-a)/h;
% L = linspace(a,b,N+1);
% S11 = zeros(1,N+1);
% S12 = zeros(1,N+1);
% S21 = zeros(1,N+1);
% S22 = zeros(1,N+1);
% S11(1) = 1;
% S12(1) = 1.2;
% S21(1) = 1;
% S22(1) = 1.5;
% for n=1:N
%     x10 = L(n);
%     y10 = S11(n);
%     z10 = S12(n);
%     x20 = L(n);
%     y20 = S21(n);
%     z20 = S22(n);
%     y11 = y10 + h*f(z10);
%     z11 = z10 + h*g(x10,y10);
%     y21 = y20 + h*f(z20);
%     z21 = z20 + h*g(x20,y20);
%     k = n+1;
%     S11(k) = y11;
%     S12(k) = z11;
%     S21(k) = y21;
%     S22(k) = z21;
% end

% a=0;
% b=1;
% h=1/3;
% g = @(x,y) 6*y^2 - x;
% f = @(z) z;
% tar = 5;
% ph = @(a) a - tar;
% N = (b-a)/h;
% L = linspace(a,b,N+1);
% S11 = zeros(1,N+1);
% S12 = zeros(1,N+1);
% S11(1) = 1;
% 
% 
% S12(1) = 1.2;
% S11 = euler(S11,S12,L,f,g,N,h);
% disp(S11);
% ph1 = ph(S11(N+1));
% 
% 
% S12(1)=1.5;
% S12 = euler(S11,S12,L,f,g,N,h);
% disp(S12);
% ph2 = ph(S12(N+1));
% 
% 
% a0 = 1.2;
% a1 = 1.5;
% eps = 0.05;
% ph3 = ph2;
% 
% while(abs(ph3)>eps)
%     a2 = a1 - ((a1-a0)/(ph2-ph1)) * ph2;
%     S12(1) = a2;
%     S12 = euler(S11,S12,L,f,g,N,h);
%     ph3 = ph(S12(N+1));
%     ph1 = ph2;
%     ph2 = ph3;
% end
% disp(S12)
% disp(ph3)
% function S11 = euler(S11,S12,L,f,g,N,h)
% for n=1:N
%     x10 = L(n);
%     y10 = S11(n);
%     z10 = S12(n);
%     y11 = y10 + h*f(z10);
%     z11 = z10 + h*g(x10,y10);
%     k = n+1;
%     S11(k) = y11;
%     S12(k) = z11;
% end
% end

u0 = 2;
t0 = 1;
tn = 1.2;
h = 0.1;

u = method(@g,u0,t0,tn,h);
for i = 1:3
    fprintf("y(%.1f) = %.6f\n",(9+i)/10,u(i));
end

syms u(t)
usol(t) = dsolve(diff(u,t) == t*t + u*u, u(1) == 2);

function u = method(g,u0,t0,tn,h)
    n = (tn-t0)/h;
    u = 1:n+2;
    t = 1:n+2;
    t(1) = t0;
    for i = 1:n+1
        t(i+1) = t(i) + h;
        u(i+1) = 0;
    end
    u(1) = u0;
    for i = 1:n+1
        e = 1;
        u(i+1) = u(i) + h*((t(i)+h/2)^2 + (u(i) + (1/2)*(h*g(u(i),t(i))))^2);
        while e>0.00005
            l = u(i+1) - u(i) - h*((t(i)+h/2)^2 + ((u(i+1)+u(i))/2)^2);
            l1 = 1 - h*(u(i+1) + u(i));
            c = u(i+1) - l/l1;
            e = c - u(i+1);
            u(i+1) = c;
        end
    end
end

function y = g(u,t)
    y = t^2 + u^2;
end


% ⁠u0 = 0;
% t0 = 1;
% tn = 1.4;
% h = 0.1;
% 
% u = predCorr(@g,u0,t0,tn,h);
% for i = 1:5
%     fprintf("y(%.1f) = %.6f\n",(9+i)/10,u(i))
% end
% 
% 
% function u = predCorr(g,u0,t0,tn,h)
%     n = (tn-t0)/h;
%     u = 1:n+2;
%     t = 1:n+2;
%     u(1) = u0;
%     t(1) = t0;
%     for i=1:n+1
%         t(i+1) = t(i) + h;
%         u(i+1) = 0;
%     end
%     j = 4;
%     for i = 1:j-1
%         u(i+1) = u(i) + h*g(u(i),t(i));
%     end
% 
%     for i = j:n+1
%         u(i+1) = u(i-3) + ((4*h)/3)*( 2*g(u(i),t(i)) - g(u(i-1),t(i-1))+ 2*g(u(i-2),t(i-2)));
%         for k = 1:2
%             u(i+1) = u(i-1) + (h/3)*(g(u(i+1),t(i+1))+4*g(u(i),t(i))+g(u(i-1),t(i-1)));
%         end
%     end
% end
% function y = g(u,t)
%     y = (1+t^2)*(u-1);
% end ⁠