% Program for hyperbolic PDEs solving the 1D avection equation 
% du/dt + a du/dx = 0 
% using Upwind, Lax-Friedrichs and Lax-Wendroff methods
%
% as of 14th Nov. 2022

%% Parameters 
global Nx CFL dt a Tmax

tau=2.5;    a=1;    D=6;

CFL = 0.9;
Nx=100;
Tmax = 5;

dx=D/(Nx-1);
dt=(CFL*dx/a);

gsin = @(t) (sin(2*pi*t/tau));
gsqr = @(t) (square(2*pi*t/tau));
ux0 = zeros(Nx,1);

%% computation
% upwind sin
sol_uw_sin = upwind(ux0,gsin);
% upwind square
sol_uw_sqr = upwind(ux0,gsqr);

if CFL == 1
    exact_sol =[sol_uw_sin sol_uw_sqr];
    save('exact_sol');
end

% lax-friedrichs sin
sol_LF_sin = lax_friedrichs(ux0,gsin);
% lax-friedrichs sqr
sol_LF_sqr = lax_friedrichs(ux0,gsqr);

% lax-wendroff sin
sol_LW_sin = lax_wendroff(ux0,gsin);
% lax-wendroff sqr
sol_LW_sqr = lax_wendroff(ux0,gsqr);
%% plot
%x=0:dx:D;
x=linspace(0,D,Nx+1);
load('exact_sol.mat')
%
figure('Name','sin-wave')
plot(x,sol_uw_sin)
hold on
plot(x,sol_LF_sin)
plot(x,sol_LW_sin)
legend('Upwind','L-Friedrichs','L-Wendroff')

figure("Name",'square-wave')
plot(x,sol_uw_sqr)
hold on
plot(x,sol_LF_sqr)
plot(x,sol_LW_sqr)
legend('Upwind','L-Friedrichs','L-Wendroff')
%}

figure
hold on
subplot(2,1,1)
plot(x,sol_uw_sin,x,sol_LF_sin,x,sol_LW_sin,x,exact_sol(:,1),'k')
title('sine-wave')
xlabel('x \rightarrow')
ylabel('u \rightarrow')
legend('Upwind','L-Friedrichs','L-Wendroff','exact Sol.')

subplot(2,1,2)
plot(x,sol_uw_sqr,x,sol_LF_sqr,x,sol_LW_sqr,x,exact_sol(:,2),'k')
xlabel('x \rightarrow')
ylabel('u \rightarrow')
legend('Upwind','L-Friedrichs','L-Wendroff','exact Sol.')
title('square-wave')