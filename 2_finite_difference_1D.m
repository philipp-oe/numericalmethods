% Program for Finite Difference aproximation in 1D 
% simulating a fluid flow through a pipe with a heating element placed
% inside of the pipe, comparing different stepsizes and velocities 
%
% as of 26th Sep. 2022

close all
clearvars
%% assuming v=1 = const.
N = [10, 20, 40, 80, 160];
h = 1./N;
v=1;
temps = {};
z = {};

for i=1:length(N)
    temps{i} = Temp(N(i),v);
    z{i} = [0:h(i):1];
end

figure
plot(z{1},temps{1},'--',z{2},temps{2},':',z{3},temps{3},'-.'...
    ,z{4},temps{4},'-')
grid on
legend('N = 10','N = 20', 'N = 40', 'N = 80')
xlabel('z \rightarrow')
ylabel('T \rightarrow')

% differnce between N and 2N

for i=1:length(N)-1
    dd(i) = abs(temps{i}(N(i)/2 +1) - temps{i+1}(N(i+1)/2 +1) );
    diff(i) = abs(temps{i}(end)-temps{i+1}(end));
end
hplot = h(1:end-1);

figure
loglog(hplot,diff,hplot,dd,hplot,hplot,hplot,hplot.^2)
legend('error @ z=1','error @ z=0.5','h','h^2')
grid on
xlabel('h \rightarrow')
ylabel('error \rightarrow')

RES = [ N(3), temps{3}(N(3)/2+ 1), temps{3}(end);
        N(4), temps{4}(N(4)/2+ 1), temps{4}(end);
        N(5), temps{5}(N(5)/2+ 1), temps{5}(end);]

%% Calculations for different velocities 
v = [1,3,10,30,100];
N = v*1E2;
h = 1./N;

for i=1:length(v)
    temps{i} = Temp(N(i),v(i));
    z{i} = [0:h(i):1];
end

figure
plot(z{1},temps{1},'-',z{2},temps{2},'-',z{3},temps{3},'-'...
    ,z{4},temps{4},'-',z{5},temps{5},'-')
grid on
legend('1','3','10','30','100')
xlabel('z \rightarrow')
ylabel('T \rightarrow')


%% Function to solve the Temperature equation

function sol = Temp(N,v)
%TEMP solving the convection diffusion equation
% sol = TEMP(steps, velocity)
% A*T = Q(z)

% parameters
L = 1;
h = L/N;

% constants
a = 0.1;
b = 0.5;
Q0 = 4000;
alpha0 = 100;
Tout = 20;
T0 = 50;

alpha = sqrt(v^2/4 + alpha0^2) - v/2;

% elements of the tridiagonal matrix A
aa = -v/(2*h) - 1/h^2;
bb = (2/h^2);
cc = v/(2*h) -1/h^2;

% vectors to build the tridiagonal martix using matlabs diag command
aaa = aa*ones(1,N-1);   % lower diagonal
bbb = bb*ones(1,N);     % main diagonal
ccc = cc*ones(1,N-1);   % upper diagonal

% matrix A initialisation (w/o boundary conditions)
A=sparse(diag(aaa,-1))+sparse(diag(bbb))+sparse(diag(ccc,+1));

% boundary condtions for z = L 
bc = 2*h*alpha;

A(N,N) = bb - cc*bc;
A(N,N-1) = aa+cc; 

% Initialising the driving force Q(z)
Q = zeros(N,1);
f = zeros(N,1);
z = [0:h:1];

% Nb = b/h +1;

for i=1:N
    if a <= z(i) && z(i) <= b
        Q(i) = Q0 * sin(pi*(z(i)-a)/(b-a));
        f(i-1) = Q(i);
    end
end

f(1) = f(1) - aa*T0;
f(N) = f(N) - cc*bc*Tout;

sol = [T0, (A\f)'];

end
