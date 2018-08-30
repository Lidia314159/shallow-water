function [u_new, v_new, h_new] = bc(u, v, h, bc_choice)
% computes boundary conditions: 
% (choice=0) original code, 
% (choice=1) reflective, 
% (choice=2) outbound
% (choice=3) periodic in x


the_globals

u_new = u;
v_new = v;
h_new = h; 

switch bc_choice
    
    case 0   
    % (original code) Define Boundary Conditions
    u_new(1,:) = 2.5*u(2,:) - 2*u(3,:) + 0.5*u(4,:);
    u_new(length(x),:) = 2.5*u(ni-1,:) - 2*u(ni-2,:) + 0.5*u(ni-3,:);
    u_new(:,1) = 2.5*u(:,2) - 2*u(:,3) + 0.5*u(:,4);
    u_new(:,length(y)) = 2.5*u(:,nj-1) - 2*u(:,nj-2) + 0.5*u(:,nj-3);
    
    v_new(1,:) = 2.5*v(2,:) - 2*v(3,:) + 0.5*v(4,:);
    v_new(length(x),:) = 2.5*v(ni-1,:) - 2*v(ni-2,:) + 0.5*v(ni-3,:);
    v_new(:,1) = 2.5*v(:,2) - 2*v(:,3) + 0.5*v(:,4);
    v_new(:,length(y)) = 2.5*v(:,nj-1) - 2*v(:,nj-2) + 0.5*v(:,nj-3);
    
    h_new(1,:) = 2.5*h(2,:) - 2*h(3,:) + 0.5*h(4,:);
    h_new(length(x),:) = 2.5*h(ni-1,:) - 2*h(ni-2,:) + 0.5*h(ni-3,:);
    h_new(:,1) = 2.5*h(:,2) - 2*h(:,3) + 0.5*h(:,4);
    h_new(:,length(y)) = 2.5*h(:,nj-1) - 2*h(:,nj-2) + 0.5*h(:,nj-3);

    case 1
    %new boundary: reflective
    h_new(1 , : )  = h(2   ,    :);
    h_new(ni, : )  = h(ni-1,    :);
    h_new(: ,  1)  = h(:   ,    2);
    h_new(: , nj)  = h(:   , nj-1);
    
    u_new(1 , : )  = -u(2   ,    :);
    u_new(ni, : )  = -u(ni-1,    :);
    u_new(: ,  1)  =  u(:   ,    2);
    u_new(: , nj)  =  u(:   , nj-1);

    v_new(1 , : )  =  v(2   ,    :);
    v_new(ni, : )  =  v(ni-1,    :);
    v_new(: ,  1)  = -v(:   ,    2);
    v_new(: , nj)  = -v(:   , nj-1);

    case 2
    %new boundary: outbound
    h_new(1 , : )  = h(2   ,    :);
    h_new(ni, : )  = h(ni-1,    :);
    h_new(: ,  1)  = h(:   ,    2);
    h_new(: , nj)  = h(:   , nj-1);

    u_new(1 , : )  =  u(2   ,    :);
    u_new(ni, : )  =  u(ni-1,    :);
    u_new(: ,  1)  =  u(:   ,    2);
    u_new(: , nj)  =  u(:   , nj-1);

    v_new(1 , : )  =  v(2   ,    :);
    v_new(ni, : )  =  v(ni-1,    :); 
    v_new(: ,  1)  =  v(:   ,    2);
    v_new(: , nj)  =  v(:   , nj-1);
    
    case 3
    %new boundary: periodic
    h_new(1 , : )  = h(ni-1,    :);
    h_new(ni, : )  = h(2   ,    :);
    h_new(: ,  1)  = h(:   , nj-1);
    h_new(: , nj)  = h(:   ,    2);
    
    u_new(1 , : )  = -u(ni-1,    :);
    u_new(ni, : )  = -u(2   ,    :);
    u_new(: ,  1)  =  u(:   , nj-1);
    u_new(: , nj)  =  u(:   ,    2);

    v_new(1 , : )  =  v(ni-1,    :);
    v_new(ni, : )  =  v(2,       :);
    v_new(: ,  1)  = -v(:   , nj-1);
    v_new(: , nj)  = -v(:   ,    2);
    %
end