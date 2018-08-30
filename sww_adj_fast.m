function [u_new, v_new, h_new] = sww_adj_fast(u, v, h)

the_globals
the_parameters
the_const_and_init

% create 9x9 matrix for each (i,j) of the grid
%
% computes model update as
% u_new = u - dudt
% v_new = v - dvdt
% h_new = h - dhdt


u_new = zeros(size(u,1), size(u,2));
v_new = zeros(size(v,1), size(v,2));
h_new = zeros(size(h,1), size(h,2));

mmtl = zeros(9);
mmadj = zeros(9);

for i=2:(ni-1)
    for j=2:(nj-1)
        %
        U = 0.25*(u(i+1,j) + u(i-1,j) + u(i,j-1) + u(i,j+1));
        V = 0.25*(v(i+1,j) + v(i-1,j) + v(i,j-1) + v(i,j+1));
        H = 0.25*(h(i+1,j) + h(i-1,j) + h(i,j-1) + h(i,j+1));
        %
        mmtl(4,1) = -U;
        mmtl(4,3) = -H;
        mmtl(5,2) = -U;
        
        mmtl(6,1) = -g;
        mmtl(6,3) = -U;
        mmtl(7,2) = -V;
        mmtl(8,1) = -V;
        mmtl(8,3) = -H;
        
        mmtl(9,2) = -g;
        mmtl(9,3) = -V;
        
        mmtl(4,4) = 1;
        mmtl(5,5) = 1;
        mmtl(6,6) = 1;
        mmtl(7,7) = 1;
        mmtl(8,8) = 1;
        mmtl(9,9) = 1;
        %
        %mmadj = mmtl';    %ADj is TL transposed
        %
        dudx = 0.5*(u(i+1,j) - u(i-1,j))/dx;
        dudy = 0.5*(u(i,j+1) - u(i,j-1))/dy;
        dvdx = 0.5*(v(i+1,j) - v(i-1,j))/dx;
        dvdy = 0.5*(v(i,j+1) - v(i,j-1))/dy;
        dhdx = 0.5*(h(i+1,j) - h(i-1,j))/dx;
        dhdy = 0.5*(h(i,j+1) - h(i,j-1))/dy;
        %
        dudt = U*dt;
        dvdt = V*dt;
        dhdt = H*dt;
        %
        %updating the state increment
        dstate = [dudt dvdt dhdt dudx dvdx dhdx dudy dvdy dhdy]';
        dstate = mmadj * dstate;
        %
        u_new(i,j) = u(i,j) - dstate(1);
        v_new(i,j) = v(i,j) - dstate(2);
        h_new(i,j) = h(i,j) - dstate(3);
    end
end

%
%new boundary: reflective
h_new(1 , : )  = h(2   ,    :);  % along y-axis
h_new(ni, : )  = h(ni-1,    :);
h_new(: ,  1)  = h(:   ,    2);  % along x-axis
h_new(: , nj)  = h(:   , nj-1);

u_new(1 , : )  = -u(2   ,    :);
u_new(ni, : )  = -u(ni-1,    :);
u_new(: ,  1)  =  u(:   ,    2);
u_new(: , nj)  =  u(:   , nj-1);

v_new(1 , : )  =  v(2   ,    :);
v_new(ni, : )  =  v(ni-1,    :);
v_new(: ,  1)  = -v(:   ,    2);
v_new(: , nj)  = -v(:   , nj-1);

%