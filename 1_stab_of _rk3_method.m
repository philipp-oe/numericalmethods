% Program regarding the Accuracy and Stabiliy of Runge-Kutta method 
% for the Landau-Lifschitz-equation for an initial magnetization 
%
% m(0) = (1,0,0)', a = 0.5*[1;1;sqrt(2)], alpha = 0.1
% dm/dt = a cross m + alpha a cross (a cross m)
% 
% as of 14th Sep. 2022

clearvars
close all

%% Starting Parameters for one solution
m0 = [1; 0; 0];
T = 40;             % End time
N = 640;            % Number of steps
h = T/N;            % Step size
Time = [0:h:T];     % Time array to plot solution over time

%% Solution and Plot for one solution

sol = RK(m0,N,h);

% x,y,z plot 

figure('Name','xyz-plot of Landau Lifshitz equation with 640 steps')
plot(Time,sol(1,:),Time,sol(2,:),Time,sol(3,:))
ylabel({'m_{x/y/z}'})
xlabel('time t')
legend({'x','y','z'})

% 3d plot

figure('Name','3d plot of Landau Lifshitz')
axis equal
plot3(sol(1,:),sol(2,:),sol(3,:))
xlabel('x')
ylabel('y')
zlabel('z')
hold on
grid on
plot3([0,0.5],[0,0.5],[0,sqrt(2)/2])
legend('Trajectory','$\vec{a}$','Interpreter','latex','Location','best')

clear a sol h N

%% Comutation for diffent number of steps to check the order of accuracy

T = 40;     % chosing the same total time as before for simplicity reasons
N_array = [20, 40, 80, 160, 320, 640];
solutions = {};
for i=1:length(N_array)
    h = T/N_array(i);
    solutions{i} = RK(m0,N_array(i),h);
end

% calculating the difference beween N and 2N

for i=1:length(N_array)-1
    diff(i) = norm(solutions{i}(:,end)-solutions{i+1}(:,end));
end

h_array = T./N_array;
h_plot = h_array(1:end-1);

figure('Name','abs. error vs h, h^2, h^3')
loglog(h_plot,diff,h_plot,h_plot,h_plot,h_plot.^2,h_plot,h_plot.^3)
xlabel('h')
ylabel('abs. error')
legend('abs error','h','h^2','h^3')
grid on

clear diff h h_array h_plot N_array solutions Time

%% Stability check for h=2 and h=2,5
T = 40;
h_stab = 2;
N = ceil(T/h_stab);
sol_stab = RK(m0,N,h_stab);
t_stab = [0:h_stab:T];

h_unstab = 2.5;
N = ceil(T/h_unstab);
sol_unstab = RK(m0,N,h_unstab);
t_unstab = [0:h_unstab:T];

figure('name','Stablility of the solution for h=2 and h=2.5')
plot(t_stab,sol_stab(1,:),t_unstab,sol_unstab(1,:))
hold on
grid on

legend('stable','unstable')
axis([0 40 -0.5 1.5])


%% Function to compute the Runge Kutta Methode
function result = RK(m0,N,h)
    % Inditional conditions
    m_array = m0;
    m = m0;

    % computation till i = N therefore all time steps are done
    for i=0:N-1
        m_next = step_RK(h,m);
        m_array(:,end+1) = m_next;
        m = m_next;
    end
    
    result = m_array;
end

%% Function to calulate the next step with RK

function u_next=step_RK(h,u)
    k1 = LLfun(u);
    k2 = LLfun(u +h*k1);
    k3 = LLfun(u +h*k1/4 +h*k2/4);
    
    u_next= u +h/6*(k1 +k2 +4*k3);
end

%% Landau-Lifshitz function impelmentation (with parameters)

function f = LLfun(m)
    % parameters:
    alpha = 0.1;
    a = 0.5*[1;1;sqrt(2)];

    axm = cross(a,m);

    f = axm + alpha*cross(a,axm);
end