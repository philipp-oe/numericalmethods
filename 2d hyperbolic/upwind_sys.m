function sols = upwind_sys(v0,f)
%UPWIND fuction solving PDE's with upwind method
%   Detailed explanation goes here

% A is the matrix infront of dv/dx
% dv/dt - Adv/dx = F
% v = (h,u)^T
global A dt Tmax dx Nx CFL

lam = CFL;

[S L] =eig(A); 
Lp = L.*(L>0);  Lm = L.*(L<0);
Ap = S*Lp*inv(S);   Am = S*Lm*inv(S);

sols = cell(2,1);
sols{1} = v0(1,:);
sols{2} = v0(2,:);

v = v0;

x = linspace(-1,2,Nx+1);

j = 0;


while j < Tmax
    for i=2:Nx             
        k = x(i);
        v_next(:,i) = v(:,i) -lam*Ap*(v(:,i)-v(:,i-1))...
            -lam*Am*(v(:,i+1)-v(:,i)) + dt*[0;f(k,j)];  
    end
    v = v_next;
    expol_left = (2*v_next(:,1) -v_next(:,2));
    expol_right = (2*v_next(:,end) -v_next(:,end-1));
    v = [expol_left v(:,2:end) expol_right];
    sols{1} = [sols{1}; v(1,:)];
    sols{2} = [sols{2}; v(2,:)];
    j = j + dt;
end



end

