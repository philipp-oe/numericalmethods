% Program for Numerical Methods
% Exercise 3 - Problem 1
% Parabolic PDE's in 1D with explicit Euler, Methode of Lines 
% and comparison to Matlab's build in functions
%
% Philipp Oelze, as of 03st Nov. 2022
clearvars
%% Parameters
L = 1;      Tmax = 2;
Nx = 50;    hx = L/Nx 
% ht = hx^2; Nt = floor(Tmax/ht);

Nt = 5000;  % stable solution
%Nt = 2900;  % unstable solution
ht =Tmax/Nt

t = [0:ht:Tmax];    x = [0:hx:L];
u0 = zeros(Nx+1,1);

%% part a 
sol = mol_sol(Nx,hx,Nt,ht);

u0t = zeros(Nt+1,1);
for i=1:Nt/2
    u0t(i) = sin(pi*t(i)); end

sol = [u0t sol];

% line below plots the figures for part a
plot1(sol,x,t,Nt)

%% part b compare to the MATLAB inbuild function ode23(s/sJ)
% parameters
N=400       
hx = L/N;
a=0.3;      
c1 = a/hx^2;
tspan = [0 Tmax];
u0 = zeros(1,N);        u(:,1) = u0;

A=gallery('tridiag',N,1,-2,1);
A(N,N-1) = 2;
A = c1*A;

% comutation ode23
opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
tic;
[t u] = ode23(@(t,u) rhs_heat(u,A,t,N,c1),tspan, u0,opts);
u0t = zeros(length(t),1);
i = 1;
while t(i) <1
    u0t(i) = sin(pi*t(i)); 
    i = i+1;
end
u = [u0t u];    % adding the boundary condtions for u(0,t)

disp('ode23')
timesteps = length(t)
compute_time = toc 

%computation ode23s
tic;
[t u] = ode23s(@(t,u) rhs_heat(u,A,t,N,c1),tspan, u0,opts);
u0t = zeros(length(t),1);
i = 1;
while t(i) <1
    u0t(i) = sin(pi*t(i)); 
    i = i+1;
end
u = [u0t u];    % adding the boundary condtions for u(0,t)

disp('ode23s')
timesteps = length(t)
compute_time = toc 

%comtutation ode23s Jacobian 
opts = odeset(opts,'Jacobian',A);
tic;
[t u] = ode23s(@(t,u) rhs_heat(u,A,t,N,c1),tspan, u0,opts);
u0t = zeros(length(t),1);
i = 1;
while t(i) <1
    u0t(i) = sin(pi*t(i)); 
    i = i+1;
end
u = [u0t u];    % adding the boundary condtions for u(0,t)

disp('ode23sJ')
timesteps = length(t)
compute_time = toc

% remove comment for plot
%{
u1 = u(:,1)';
x=[0:hx:L];
[X T] =meshgrid(x,t);
figure 
mesh(X,T,u)
xlabel('x \rightarrow')
ylabel('\leftarrow t')
zlabel('Temp. (arb. unit) \rightarrow')
%}

%% functions
function sol = mol_sol(Nx,hx,Nt,ht)
% solves du/dt = rhs_heat with explicit euler, as part of methods of lines
N=Nx;
a= 0.3;
c1 = a/hx^2;

u0 = zeros(N,1);
u(:,1) = u0;

A=gallery('tridiag',N,1,-2,1);
A(N,N-1) = 2;

A = c1*A;

% Euler Method u(x,t+1) = u(x,t) + ht( A*u(x,t) + b(t))
for i = 1:Nt
    t = (i-1)*ht;
    u(:,i+1) = u(:,i) + ht*( rhs_heat(u(:,i),A,t,N,c1) );
end

sol =[u]';

end

function fun=rhs_heat(u,A,t,N,c1)
% right hand side of the PDE du/dt = a*d^2u/dt^2
    b=zeros(N,1);
    if t <= 1;
        b(1) = c1*sin(pi*t); end
    fun = A*u+b;
end

function plot1(sol,x,t,Nt)
    % 3D Plot 
    [X T] = meshgrid(x,t);

    figure
    mesh(X,T,sol);
    xlabel('x \rightarrow')
    ylabel('\leftarrow t')
    zlabel('Temp. (arb. unit) \rightarrow')

    % 2D Plot at t=1
    t1 = Nt/2;
    u_t1 = sol(t1,:);

    figure
    plot(x,u_t1,'r*')
    xlabel('x \rightarrow')
    ylabel('Temp. (arb. unit) \rightarrow')

    figure 
    plot(t,sol(:,1),'r',t,sol(:,end),'k')
    xlabel('t \rightarrow')
    ylabel('Temp. (arb. unit) \rightarrow')
    legend('left side', 'right side')

end