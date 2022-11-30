function sol = lax_wendroff(u0,g)
%LAX_WENDROFF Method to solve the avection equation u_t + au_x = 0 
%
%   only the last value of the solution is saved, hence only this one is 
%   required by the problem

global Nx CFL a dt Tmax

lam = a*CFL;
aa = lam/2 + lam^2/2;       % u^n_{j-1}
bb = 1 - lam^2;             % u^n_j
cc = -lam/2 + lam^2/2;      % u^n_{j+1}

A = gallery('tridiag',Nx-1,aa,bb,cc);
b = zeros(Nx-1,1);
u=u0(2:end);

for k = dt:dt:Tmax
    %extra polation
    uN = 2*u(end-1) - u(end-2);
    b(1) = lam*g(k);
    b(end) = uN*cc;
    u_time_next = A*u +b;
    u = u_time_next;
end

sol = u;
sol = [g(k);sol];


end

