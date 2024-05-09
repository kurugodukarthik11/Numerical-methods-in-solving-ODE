
% %Find only one root of the equation x^3-5x+1=0; See the following code;

c=zeros(1,60);
x=0.5;

  eps = 1; tol = 10^(-14); total = 100; k = 0; format long; 
 while ((eps > tol)&&(k < total))
   f = x^3-5*x+1;
   f1=3*x^2-5;
   xx = x-f/f1;
   eps = abs(xx-x); x = xx;
   fprintf('k= %2.0f, The root = %12.12f\n\n',k,x);
   k = k+1;
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%Find all roots of the equation x^3-5x+1=0; See the following code;


% c=zeros(1,2);
% c(1)=0.5;
% 
% 
% x=c(1);
% for n=1:1:2
%   eps = 1; tol = 10^(-14); total = 100; k = 0; format long; 
%  while ((eps > tol)&&(k < total))
%    f = x^3-5*x+1;
%    f1=3*x^2-5;
%    xx = x-f/f1;
%    eps = abs(xx-x); x = xx;
%    %fprintf('k= %2.0f, The root = %12.12f\n\n',k,x);
%    k = k+1;
%  end
% c(n)=x;
% x=c(n)+2;
% fprintf('n= %2.0f, The root = %12.12f\n\n',n,c(n));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Another Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % ROOT FINDING %%%%%%%%%%%%%
% 
% F=@(x) x^3-5*x+1;
% interval = [0, 3];
% N = 50;
% start_pts = linspace(interval(1), interval(2), N);
% found_roots = [];
% for i = 1:numel(start_pts)-1
%     try
%         found_roots(end+1) = fzero(F, [start_pts(i), start_pts(i+1)]);
%     end
% end
% mu_n = found_roots
% mu_0 = found_roots(1);
% Nmax = length(found_roots);



