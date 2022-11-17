% Program for Numerical Methods
% Exercise 3 - part 2
% Crank Nicelson method for 2D PDE
%
% Philipp Oelze as of 9th Nov. 2022 
% running version as of 09.11.2020 16:13
clearvars 
close all
%% Parameters and main

% for animation set q = 0 and chose h = 0.1
q = 1;

Lx = 10;    Ly = 6;     L = [Lx Ly];
h = 0.05; ht = h;
M = Lx/h; N = Ly/h;      % #total inner points (w/o ghost points)
M = M+1;
Tmax = 50;
tsteps = Tmax/ht;
T=[0:ht:Tmax];

% Matrix Implementation
A = getA(M,N,h);
k = length(A);
B = sparse(eye(k) - ht/2*A);
C = sparse(eye(k) + ht/2*A);
clear k

fun = get_f(L(1),L(2),h);

T0 = 20*ones(M,N+1)';
temps = cell(4,1);
temps{1} = T0;
k = 20*ones(1,M);

u = reshape(T0(2:end,:),[],1);
for i=2:tsteps+1
    u_next = cn_step(B,C,h,u,fun);
    u = u_next;
    % saveing t=1,10,50
    if i==(tsteps/50)
        temps{2} = [k; reshape(u,M,N)'];
    elseif i==(tsteps/5)
        temps{3} = [k; reshape(u,M,N)'];
    elseif i==tsteps
        temps{4} = [k; reshape(u,M,N)'];
    end

    %saving all (5,3) points for all t 
    uu = [k; reshape(u,M,N)'];
    t53(i-1) = uu((N+2)/2,(M+1)/2);
     
    % saving all time stept for the animation
    if q==1
        temps{i,1} = [k; reshape(u,M,N)']; end
end  

% plots for t=0,1,10,50
plots(L,h,temps)

t53 = [20 t53];
%plots temperature at(5,3) over time
figure('Name','plot1')
plt1 = plot(T,t53);
hold on
datatip(plt1,'DataIndex',length(t53), ...
    'location','northwest');
xlabel('t \rightarrow')
ylabel('Temp. \rightarrow')
hold off

%% Animation for t=0...50
if q == 1
    figure 
    [X,Y] = meshgrid([0:h:Lx],[0:h:Ly]);
    for i=1:length(temps)
        mesh(X,Y,temps{i})
        axis([0 10 0 6 20 40])
        xlabel('x \rightarrow')
        ylabel('y \rightarrow')
        zlabel('Temp. \rightarrow')
        pause(0.05)
    end
end

%% function for A Matrix 
function A = getA(M,N,h)

    SM = gallery('tridiag',M,1,-2,1);
    SM(1,2) = 2;
    SM(end,end-1) = 2;
    SM = SM/h^2;

    SN = gallery('tridiag',N,1,-2,1);
    SN(end,end-1) = 2;
    SN = SN/h^2;

    A = sparse( kron(eye(N),SM) + kron(SN,eye(M)));
end
%% Funtion for crank nicelson
function step = cn_step(B,C,ht,u,fun)
    step = B\(C*u + ht*fun);
end
%% Function f
function fun = get_f(Lx,Ly,h)
    T0 = 20;
    x = [0:h:Lx];
    y = [h:h:Ly];

    for i=1:length(x)
        for j=1:length(y)
            f(i,j) = 20*exp(-0.1*(x(i)-3)^2 -2*(y(j)-1)^2);
        end
    end
    f(:,1) = f(:,1) + T0/h^2;
    fun =reshape(f,[],1);
end

%% Plot function
function plots(L,h,temps)

[X,Y] = meshgrid([0:h:L(1)],[0:h:L(2)]);

figure 
mesh(X,Y,temps{1}) 
axis([0 10 0 6 20 40])
xlabel('x \rightarrow')
ylabel('y \rightarrow')
zlabel('Temp. \rightarrow')


figure 
mesh(X,Y,temps{2}) 
axis([0 10 0 6 20 40])
xlabel('x \rightarrow')
ylabel('y \rightarrow')
zlabel('Temp. \rightarrow')


figure
mesh(X,Y,temps{3})  
axis([0 10 0 6 20 40])
xlabel('x \rightarrow')
ylabel('y \rightarrow')
zlabel('Temp. \rightarrow')


figure 
mesh(X,Y,temps{4}) 
axis([0 10 0 6 20 40])
xlabel('x \rightarrow')
ylabel('y \rightarrow')
zlabel('Temp. \rightarrow')


end

