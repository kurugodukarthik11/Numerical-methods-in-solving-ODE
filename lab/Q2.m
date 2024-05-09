u0 = 0;
t0 = 1;
tn = 1.4;
h = 0.1;

u = q2(@g,u0,t0,tn,h);
for i = 1:5
   fprintf("u(%.1f) = %.6f\n",(9+i)/10,u(i))
end

function u = q2(g,u0,t0,tn,h)
    n = (tn-t0)/h;
    u = 1:n+2;
    t = 1:n+2;
    u(1) = u0;
    t(1) = t0;
    for i=1:n+1
        t(i+1) = t(i) + h;
        u(i+1) = 0;
    end
    for i = 1:3
        u(i+1) = u(i) + h*g(u(i),t(i));
    end
    for i=1:4
        fprintf("Euler method u(%.1f) = %.6f\n",(9+i)/10,u(i))
    end
    for i = 4:n+1
        u(i+1) = u(i-3) + ((4*h)/3)*( 2*g(u(i),t(i)) - g(u(i-1),t(i-1))+ 2*g(u(i-2),t(i-2)));
        fprintf("Predicted method");
        disp(u(i+1));
        for k = 1:2
            u(i+1) = u(i-1) + (h/3)*(g(u(i+1),t(i+1))+4*g(u(i),t(i))+g(u(i-1),t(i-1)));
            fprintf("Corrected method, Iteration number %d:%12.12f\n",k,u(i+1));
        end
    end
end

function y = g(u,t)
    y = (1+t^2)*(u-1);
end