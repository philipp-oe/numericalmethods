% Program for a three-body problem (Earth-Moon-Satellite)
% using Adam Bashford 4 and ode23 in comparsion to solve the problem
% 
% as of 16th Sep. 2022

clearvars
close all

%% Solving the ODE using Adam Bashford 4

% starting contitions 
R0 = [1.2; 0; 0; -1];       % p1 = x, p2=y, p3=v_x, p4=v_y
T = 7;      % max time  
h = 1E-4;   % step size
N = T/h;    % step count

t=[0:h:T];      % time array for using in the plot

R = AB4(R0,N,h);        % computation using the Adam Bashford 4 method

figure('Name','position-time plot')
plot(t,R(1,:),t,R(2,:))
grid on
xlabel('time')
ylabel('postion')
legend('x','y')


figure('Name','velocity-time-plot')
plot(t,R(3,:),t,R(4,:))
grid on
xlabel('time')
ylabel('velocity')
legend('x','y')

figure('name','Trajectory')
plot(R(1,:),R(2,:),'r')
grid on
hold on
mu=1/82.45;
plot(mu,0,'bX','LineWidth',2)     % Marker Earth position
plot(1-mu,0,'kX','LineWidth',2)   % Marker Moon position
xlabel('r_x')
ylabel('r_y')
legend('Trajectory','Earth','Moon')
axis equal

figure('name','v - Trajectory')
plot(R(3,:),R(4,:),'r')
grid on
hold on
xlabel('v_x')
ylabel('v_y')
axis equal

%% 2nd part calculating T_acc

tol = 0.015;
h = 0.0005;

method = 3;
% first step to calculate an initial value set
res = Result(R0,N,h,method);
% halving of the step size
res2help =Result(R0,2*N,h/2,method);
% deleting every second value, so res and res2 are the same lenght
res2=res2help(:,1:2:end);
%calculating the differenc between both
diffres = res(1:2,:)-res2(1:2,:);   % only the positions are from interest
for kk=1:N
    if norm(diffres(1:2,kk)) > tol
        %Tacc is time when difference is too big 
        nmin=kk;
        Tacc = kk*h-h
        break
    else 
        Tacc = 20;
    end

end

tt = [0:h:Tacc];
N = Tacc/h;
RR = Result(R0,N,h,method);

figure('Name','AB4 xy t = 0...T_{acc}')
plot(RR(1,:),RR(2,:),'r')
grid on
hold on
mu=1/82.45;
plot(mu,0,'bX','LineWidth',2)     % Marker Earth position
plot(1-mu,0,'kX','LineWidth',2)   % Marker Moon position
xlabel('r_x')
ylabel('r_y')
legend('Trajectory','Earth','Moon')
axis equal



%% Step size plot of ode23

res = Result(R0,N,h,4);
t = res(5,:);

figure
plot(t,res(1,:),t,res(2,:));
xlabel('time')
ylabel('position')
legend('x','y')

figure
plot(res(1,:),res(2,:),'r');
grid on
hold on
mu=1/82.45;
plot(mu,0,'bX','LineWidth',2)     % Marker Earth position
plot(1-mu,0,'kX','LineWidth',2)   % Marker Moon position
xlabel('r_x')
ylabel('r_y')
legend('Trajectory','Earth','Moon')
axis equal

figure('name','v - Trajectory')
plot(res(3,:),res(4,:),'r')
grid on
hold on
xlabel('v_x')
ylabel('v_y')
axis equal

%calulate the stepsize
h = diff(t);

figure 
plot(t(2:end),h,'*')
xlabel('time')
ylabel('step size')

%% Result function for easier handeling of the second question
function res = Result(R0,N,h,method)
    % Explicit Euler
    if method == 1 
        res =exEuler(R0,N,h);

    % Runge Kutta 3
    elseif method == 2
        res =RK(R0,N,h);

    % Adam Bashford 4
    elseif method == 3
        res = AB4(R0,N,h);

    % ODE23
    elseif method == 4
        R(:,1) = R0;
        options = odeset(RelTol=1E-5);
        [ti, R] = ode23(@dfun,[0 19.7725],R(:,1),options);
        res = [R'; ti'];
    else
        disp('wrong methode counter')
    end

end

%% Function for Adam Bashford 4
function res = AB4(R0,N,h)

    % calulate the first 4 steps using RK3
    U = RK(R0,3,h);
    t=0;
    % Adam Bashford 
    for i=4:N
        u_next = U(:,i) + h/24*( 55*dfun(t,U(:,i)) -59*dfun(t,U(:,i-1)) ...
            +37*dfun(t,U(:,i-2)) -9*dfun(t,U(:,i-3)) );
        U(:,end+1) = u_next;
    end
    res = U;
end

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
    t=0;
    k1 = dfun(t,u);
    k2 = dfun(t+h, u +h*k1);
    k3 = dfun(t+h/2,u +h*k1/4 +h*k2/4);
    
    u_next= u +h/6*(k1 +k2 +4*k3);
end

%% Explicit Euler 
function res = exEuler(R0,N,h)
    res = R0;
    u = R0;
    t = 0;
    for i=1:N
        u_next =  u + dfun(t,u)*h;
        res(:,end+1) = u_next;
        u = u_next;
    end
end

%% function with x,y,v_x,v_y
function d = dfun(t,R)
    mu = 1/82.45;
    r0 = [-mu;0];
    r1 = [1-mu;0];
    rr0 = R(1:2);

    dearth = norm(rr0-r0);
    dmoon = norm(rr0-r1);

    % (1,2,3,4)^T = (x,y,v_x,v_y)^T

    % postion derevatives
    d(1) = R(3);
    d(2) = R(4);
    
    % velocity derivatives
    d(3) = -(1-mu)*(R(1)-r0(1))/dearth^3 ...
           -mu*(R(1)-r1(1))/dmoon^3 ...
           +2*R(4)+R(1);
    d(4) = -(1-mu)*(R(2)-r0(2))/dearth^3 ...
           -mu*(R(2)-r1(2))/dmoon^3 ...
           -2*R(3)+R(2);
    d = d';
end
