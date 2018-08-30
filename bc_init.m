function [u_new, v_new, h_new] = bc_init(u, v, h)
% initialize boundary conditions: 
% (choice=0) original code does none, 
% (choice=1) reflective, 
% (choice=2) outbound
% (choice=3) periodic in x


the_globals
the_parameters
the_const_and_init

u_new = u;
v_new = v;
h_new = h; 

switch bc_choice
    
    case 0   
    % in the original code there was no initializing the boundary conditions
    return

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
    
    u_new(1 , : )  =  u(ni-1,    :);
    u_new(ni, : )  =  u(2   ,    :);
    u_new(: ,  1)  = -u(:   , nj-1);
    u_new(: , nj)  = -u(:   ,    2);

    v_new(1 , : )  = -v(ni-1,    :);
    v_new(ni, : )  = -v(2,       :);
    v_new(: ,  1)  =  v(:   , nj-1);
    v_new(: , nj)  =  v(:   ,    2);
    %

end