% Program for Finite Difference aproximation in 2D
% heat flow through a solid metall block, described by -laplace T = f
% 
% as of 4th Oct. 2022

%% assuming a constant function f

Lx = 10;    Ly = 6;
h = 0.2;
M = Lx/h-1;   N = Ly/h-1; % total inner points( w/o ghost points)

fun = get_f_a(M,N,h);
temps = temp2(M,N,h,fun);
[X,Y] = meshgrid([0:h:Lx],[0:h:Ly]);
temps = [20*ones(1,M); temps; temps(end,:)];
temps = [temps(:,1), temps, temps(:,end)];

mesh(X,Y,temps)
xlabel('x')
ylabel('y')

disp('Temperature at (5,3) for f=const.') 
disp(temps((N+3)/2,(M+3)/2))

%% calculation with a none constant function
temps = cell(3,1);
h = [0.2, 0.1, 0.05, 0.025, 0.0125];
M = Lx./h -1;
N = Ly./h -1;

for i=1:5
    tic;
    fun = get_f_c(Lx,Ly,h(i));
    temps{i} = temp2(M(i),N(i),h(i),fun);
    temps{i} =[20*ones(1,M(i)); temps{i}; temps{i}(end,:)];
    temps{i} = [temps{i}(:,1), temps{i}, temps{i}(:,end)];
    toc
    disp('step done')
end

figure
mesh(X,Y,temps{1})
xlabel('x')
ylabel('y')

[X,Y] = meshgrid([0:h(2):Lx],[0:h(2):Ly]);
figure
contour(X,Y,temps{2},'ShowText','on')
xlabel('x')
ylabel('y')
 
figure
imagesc([0:h(2):Lx],[0:h(2):Ly],temps{2})
xlabel('x')
ylabel('y')
colorbar

for i=1:3
    k = ['Temperature at (5,3) with h=',num2str(h(i))];
    disp(k)
    disp(temps{i}((N(i)+3)/2,(M(i)+3)/2)) 
end

%% check for second order accuracy 
for i=1:length(h)-1
    diff(i) = abs(temps{i}((N(i)+3)/2,(M(i)+3)/2)...
        -temps{i+1}((N(i+1)+3)/2,(M(i+1)+3)/2) );
end
hplot = h(2:end);

figure
loglog(hplot,diff,hplot,hplot,hplot,hplot.^2)
legend('error','h','h^2')
grid on
%% Function for Finite Differences in 2D
function res = temp2(M,N,h,fun)
% with M - x step, N - y steps, h - stepsize

    a = 2*eye(M);
    b = -ones(M-1,1);
    SM = 1/h^2*(a+diag(b,1)+diag(b,-1));
    SM(1,2) = -2/h^2;
    SM(end,end-1) = -2/h^2;
    SM =sparse(SM);

    c = 2*eye(N);
    d = -ones(N-1,1);
    SN = 1/h^2*(c+diag(d,1)+diag(d,-1));
    SN(end,end-1) = -2/h^2;
    SN = sparse(SN);

clear a b c d

    A = kron(eye(N),SM) + kron(SN,eye(M));
    A = sparse(A);

clear SN SM

    reshelp =A\fun; 
    res = reshape(reshelp,M,N)';
end

%% constant function f_a
function fun = get_f_a(M,N,h)
    T0 = 20;
    f = 5*ones(M,N);
    f(:,1) = f(:,1)+ T0/h^2;
    fun = reshape(f,[],1);
end

%% none constant function f_c
function fun = get_f_c(Lx,Ly,h)
    T0 = 20;
    x = [h:h:Lx-h];
    y = [h:h:Ly-h];

    for i=1:length(x)
        for j=1:length(y)
            f(i,j) = 20*exp(-0.1*(x(i)-3)^2 -2*(y(j)-1)^2);
        end
    end
    f(:,1) = f(:,1)+ T0/h^2;
    fun =reshape(f,[],1);
end
