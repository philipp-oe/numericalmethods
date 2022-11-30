function sols = Lax_wendroff2(v0,f)

global Nx dx dt Tmax A CFL
lam =CFL;

sols = cell(2,1);
sols{1} = v0(1,:);
sols{2} = v0(2,:);

v = v0;

x = linspace(-1,2,Nx+1);

j = 0;
while j < Tmax
    for i =2:Nx;
        % adjusted source term f
        fun = [0; f(x(i),j)] -lam/4*A*[0; f(x(i+1),j)-f(x(i-1),j)]...
            +1/2*[0; f(x(i),j+dt)- f(x(i),j)];

        v_next(:,i) = v(:,i) -lam/2*A*( v(:,i+1)-v(:,i-1) )...
            + lam^2/2 *A^2 *( v(:,i+1) -2*v(:,i) +v(:,i-1) )...
            + dt*fun;
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

