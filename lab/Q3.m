q3()

function q3()
    h = 0.25;
    x0 = 0;
    u0 = 1/3;
    s1 = 0.005;
    s2 = 0.2;
    tar = 0.5;
    u11= euler(u0, s1, h);
    u12 = euler(u0, s2, h);
    fprintf("Initial slope: %.2f u(1): %.6f\n",s1,u11)
    fprintf("Initial slope: %.2f u(1): %.6f\n",s2,u12)
    for i= 1:1:2  
        s3 = s2 - (u12 - tar) * (s2 - s1) / (u12 - u11);
        s1 = s2;
        s2 = s3;
        u11 = u12;
        u12 = euler(u0, s3, h);
        fprintf("Iteration: %d Initial slope: %12.12f u(1): %12.12f\n",i,s3,u12);
    end
    fprintf("Error %12.12f:\n",abs(u12-0.5));
end
function u111 = euler(u0, v0, h)
    x = 0;
    u = u0;
    v = v0;
    while x < 1
        u1 = u + h * v;
        v1 = v + h * 2 * u * v;
        u = u1; v = v1;
        x = x + h;
    end
    u111 = u;
end