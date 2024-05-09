% h=1/2;
% f=@(x,y,z) z;
% g=@(x,y,z) 6*y^2-x;
% a=0;
% b=1;

yb=5;
m0=1.2;%%%%%Initial guess%%%%%%
m1=1.5;%%%%%Initial guess%%%%%%
x=m1;
y=m0;

err=0.0005;%%%%tolerance%%%%%%

iter=100;

FinalAlpha=0;
for i=1:1:iter
fprintf('iteration\n=%2.0f\n\n',i);
phi0=EulerFun(m0)-yb;
phi1=EulerFun(m1)-yb;
% phi1-phi0
m2=m1-((m1-m0)/(phi1-phi0))*phi1;%%%%Secant method%%%%%%%
yf=EulerFun(m2);
if(abs(yf-yb)<=err)
    FinalAlpha=m2;
    fprintf('Final value of Alpha\n=%12.12f\n\n',FinalAlpha);
    %fprintf('Final values of y\n=%12.12f\n\n',yf);
    break
else
    m0=m1;
    m1=m2;
end
end