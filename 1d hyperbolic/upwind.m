function sol = upwind(u0,g)
%
%UPWIND scheme solving the avection equation u_t + au_x = 0 
%   with a = 1 therefore only using FTBS 
%   with lambda = CLF number / dt/dx
%
%   only the last value of the solution is saved, hence only this one is 
%   required by the problem
global Nx CFL dt a Tmax

lam = a*CFL;
aa = lam;
bb = 1 -lam;
A = gallery('tridiag',Nx-1,aa,bb,0);
b = zeros(Nx-1,1);
u=u0(2:end);
 
for k = dt:dt:Tmax
    b(1) = aa*g(k);
    u_time_next = A*u +b;
    u = u_time_next;
end
sol = u;
sol = [g(k);sol];
end

