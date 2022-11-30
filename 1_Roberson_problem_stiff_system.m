% Program to solve a stiff problem, here done with the Robertson Problem
% modeling the reactions of three chemicals A, B and C,
% computation time for Runge-Kutta 3 and implicit Euler method are compared
% 
% 
% as of 15th Sep. 2022

clearvars
close all
%% Solving with RK-3

x0 = [1;0;0];
N_array = 1E3*[2,4,8.5,10,12,16];
T = 10;
k = length(N_array); 
sols = cell(k,1);

for i=1:k
    N = N_array(i);
    h = T/N;
    sols{i} = RK(x0,N,h);
end

% plot of solution for N =12000
i=5;
t = [0:T/N_array(i):T];

figure('Name','Robertson Problem using RK-3 and N=12000')
plot(t,sols{i}(1,:),t,sols{i}(2,:),t,sols{i}(3,:))
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')

figure('Name','semilogy-plot RP using RK-3 and N=12000')
semilogy(t,sols{i}(1,:),t,sols{i}(2,:),t,sols{i}(3,:))
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')

% plot for N=10000, bc there is an unstable region at the end
i=4;
t = [0:T/N_array(i):T];

figure('Name','semilogy-plot RP using RK-3 and N=10000')
semilogy(t,sols{i}(1,:),t,sols{i}(2,:),t,sols{i}(3,:))
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')

clear h i k N t 

%% Stability and eigenvalues of the Jacobi matrix
X = sols{4};

maxEV = max(abs(JacobiEV(X)));
t = [0:T/N_array(4):T];

figure('Name','max ev of the jacobi matrix for a nearly stable solution N=10k')
plot(t,maxEV)
xlabel('Time')
ylabel('max. EV of the Jacobi Matrix')
legend('N=10k')

X = sols{5};

maxEV = max(abs(JacobiEV(X)));
t = [0:T/N_array(5):T];

figure('Name','max ev of the jacobi matrix for a stable solution N=12k')
plot(t,maxEV)
xlabel('Time')
ylabel('max. EV of the Jacobi Matrix')
legend('N=12k')

X = sols{3};

maxEV = max(abs(JacobiEV(X)));
t = [0:T/N_array(3):T];

figure('Name','max ev of the jacobi matrix for a unstable solution N=8.5k')
plot(t,maxEV)
xlabel('Time')
ylabel('max. EV of the Jacobi Matrix')
legend('N=8.5k')

%% Solving for T:[0,1000] with RK-3
T = 1000;
N = 5E6;
h=T/N;
t = [0:h:T];


sol_RK = RK(x0,N,h);    
tic;
RKend = RK_end(x0,N,h);     % adapted RK3 that only saves the last value
T_rk1000 = toc;

disp('time and endvalues for t=0...1000 with RK3')
format short
h
T_rk1000
format long
RKend 
format short 

% Plot to check if the solution is stable
figure('Name','RK3 plot for t=0...1000 with N=5E6')
plot(t,sol_RK)
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')

figure('Name','RK3 semilogy plot for t=0...1000 with N=5E6')
semilogy(t,sol_RK)
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')
exact = [0.33687453061;
         0.0000020137023;
         0.6631234557];

disp('error for RK3')
%error = abs(RKend-exact)./exact
error = norm(RKend -exact)/norm(exact)
%% solving using implicit euler IE

% for T = 10
T = 10;
h = 0.5;
q = 1;      % for q != to 1 only the last value of impeuler is returned

[t, sol_IE] = impeuler(T,x0,h,q);

figure('Name','solution IE with h=1')
plot(t,sol_IE)
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')

% for T = 1000
T = 1000;
h = 5;
q = 1;      % for q != to 1 only the last value of impeuler is returned

[t, sol_IE] = impeuler(T,x0,h,q);

figure('Name','solution IE with h=5')
plot(t,sol_IE)
xlabel('Time')
ylabel('scaled concentrations')
legend('x_1','x_2','x_3')



%% errors using implicit euler IE
T = 1000;
h = 1E-3;
q = 0;      % for q != to 1 only the last value of impeuler is returned

tic;
[t, sol] = impeuler(T,x0,h,q);
toc
sol(:,end);
disp('error for IE')
error = norm(sol-exact)/norm(exact)

% this calcualtion uses a changend version of impeuler, where only the last
% value of the solution is returned
% h = [1, 0.1, 0.01, 0.001, 1E-4];
% k =length(h);
% Time = zeros(1,k);
% sol_ie = zeros(3,k);
% 
% for i=1:k
%     start = tic;
%     [t, sol] = impeuler(T,x0,h(i));
%     Time(i) = toc;
%     sol_ie(:,1) = sol;
% end
% 
% error = abs(sol_ie-exact)./exact

%% Function to compute the Runge Kutta Methode
function result = RK(x0,N,h)
    % Inditional conditions
    x_array = x0;
    x = x0;

    % computation till i = N therefore all time steps are done
    for i=1:N
        x_next = step_RK(h,x);
        x_array(:,end+1) = x_next;
        x = x_next;
    end
    
    result = x_array;
end

%% adapted RK Methode with only solving the last point of the solution
function result = RK_end(x0,N,h)
    % Inditional conditions
    x_array = x0;
    x = x0;

    % computation till i = N therefore all time steps are done
    for i=1:N
        x_next = step_RK(h,x);
        x = x_next;
    end
    
    result = x_next;
end

%% Function to calulate the next step with RK

function u_next=step_RK(h,u)
    k1 = Rfun(u);
    k2 = Rfun(u +h*k1);
    k3 = Rfun(u +h*k1/4 +h*k2/4);
    
    u_next= u +h/6*(k1 +k2 +4*k3);
end

%% Function to describe the Robertson Problem

function fun = Rfun(x)
    r1=0.04;
    r2=1E4;
    r3=3E7;
    f = [-r1*x(1)+r2*x(2)*x(3);
         r1*x(1)-r2*x(2)*x(3)-r3*x(2)^2;
         r3*x(2)^2;
         ];
    fun=f;
end

%% Jacobi Matrix Eigenvalues

function A = JacobiEV(X)
    r1=0.04;
    r2=1E4;
    r3=3E7;

    J = zeros(3,3);

    for i=1:length(X)
        x = X(:,i);
        J = [r1, r2*x(3), r2*x(2);
             -r1, -r2*x(3)-2*r3*x(2), -r2*x(2);
             0, 2*r3*x(2), 0];
        A(:,i) = eig(J);
    end
end

function [t y]=impeuler(Tend, y0, h, q)
%IMPEULER  Solve Robertson problem with the Implicit Euler method.
%   [TOUT,YOUT] = IMPEULER(T,Y0,h) computes the concentrations xA, xB, xC from 
%   time t=0 to t=T using timestep h, and initial data Y0=[xA(0),xB(0),xC(0)]. 
%   The solution is returned in the array YOUT, where each row contains the 
%   concentrations xA, xB, xC at at the time given in the accompanying column 
%   vector TOUT.

r1 = 0.04;  % Rate constants
r2 = 1e4;
r3 = 3e7;

f = @(u) [-r1*u(1)+r2*u(2)*u(3); r1*u(1)-r2*u(2)*u(3)-r3*u(2)^2; r3*u(2)^2];  % RHS function
Jf = @(u) [-r1 r2*u(3) r2*u(2); r1 -r2*u(3)-2*r3*u(2) -r2*u(2); 0 2*r3*u(2) 0]; % Jacobian of f
I = eye(3,3);

tol = 1e-10;  % Tolerance in Newton

t = 0:h:Tend;
u = y0(:);    % Make sure u is a column vector
y = u';


for t1=t(2:end)
    d = 1;  % Dummy value
    v = u;  % Starting guess for u_{n+1} taken as u_n
    
    while (norm(d)>tol)    % Solve F(v) = v - h*f(v) - u_n = 0, with Newton's method
        F = v - h*f(v) - u;
        J = I - h*Jf(v);   % Jacobian of F
        d = -J\F;
        v = v + d;
    end
    
    u = v;        % u_{n+1} = v, solution to F(v)=0.
    if q == 1
        y = [y; u'];  % Save values for plotting
    end
end
    if q == 1
        y = y;
    else
        y = v;
    end

t = t';
end
