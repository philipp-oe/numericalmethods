% Program for 2D hyperbolic PDE solving waterwaves in a small channel
% governed by the linearized Saint–Venant equations
%
% Philipp Oelze, as of 17th Nov. 2022
%% Parametes
global A Tmax dx Nx CFL dt
g = 9.81; V0 = 2; H0 = 3;
L = 3; L0 = -1; L1 = 2;
w = 0.15;

%% 
A = [V0 H0; g V0];
disp(eig(A)) % check if system is parabolic -> no zero eigenvalues
% how far the waves propagate in t=0.2sec

Tmax = 0.2;
Nx = 1000;
dx = L/(Nx+1);
dt = Tmax/(Nx);
CFL = dt/dx;

vx0 = zeros(2,Nx+1);

f= @(x,t) -g*((sin(30*pi*t +pi/6) >0.5).*(abs(x)<w).*sin(pi*x/w));

sol_uw = upwind_sys(vx0,f);
h_uw = sol_uw{1};
u_uw = sol_uw{2};

sol_lw = Lax_wendroff2(vx0,f);
h_lw = sol_lw{1};
u_lw = sol_lw{2};

x = linspace(L0,L1,Nx+1);
t = [0:dt:Tmax];
[X T] = meshgrid(x,t);

figure
mesh(X,T,h_uw)

figure
mesh(X,T,h_lw)

figure
subplot(2,1,1)
plot(x,h_uw(end,:),x,h_lw(end,:))
legend('h_{uw}','h_{lw}')
subplot(2,1,2)
plot(x,u_uw(end,:),x,u_lw(end,:))
legend('u_{uw}','u_{lw}')


