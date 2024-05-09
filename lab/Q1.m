u0=2;
t0=1;
tn=1.2;
h=0.1;

u=q1(@g,u0,t0,tn,h);
for i = 1:3
    fprintf("Implicit y(%.1f) = %.6f\n",(9+i)/10,u(i));
end

syms u1(t)
usol(t) = dsolve(diff(u1,t) == t*t + u1*u1, u1(1) == 2);
exact=zeros(1,4);
for i = 1:3
    exact(i) = usol((9+i)/10);
    fprintf("Exact y(%.1f) = %.6f\n",(9+i)/10,usol((9+i)/10));
end

Err = zeros(1,4);
for i = 1:3
    Err(i) = exact(i)-u(i);
    fprintf("Error at x = %.1f is %.6f\n",(9+i)/10,abs(usol((9+i)/10)-u(i)));
end 

function u = q1(g,u0,t0,tn,h)
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
        tol = 1;
        u(i+1) = u(i) + h*((t(i)+h/2)^2 + (u(i) + (1/2)*(h*g(u(i),t(i))))^2);
        while tol>0.0005
            l = u(i+1) - u(i) - h*((t(i)+h/2)^2 + ((u(i+1)+u(i))/2)^2);
            l1 = 1 - h*(u(i+1) + u(i));
            u1 = u(i+1) - l/l1;
            tol = abs(u1-u(i+1));
            u(i+1) = u1;
        end
    
    end
end

function y = g(u,t)
    y = t^2 + u^2;
end