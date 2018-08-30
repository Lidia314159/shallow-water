%% TANGENT LINEAR SOLVE
clear all;

g = 9.8; % units of m/s^2

ni = 100; % number of gridpoints in x
nj = 100; % number of gridpoints in y
nt = 10000; % number of timesteps

dt = 0.01; % timestep
dx = 0.5; % xgrid step
dy = 0.5; % ygrid step

u_initial = zeros(ni,nj); % initialize grids
v_initial = zeros(ni,nj); % initialize grids
h_initial = zeros(ni,nj); % initialize grids
u = zeros(ni,nj); u_old = zeros(ni,nj);
v = zeros(ni,nj); v_old = zeros(ni,nj);
h = zeros(ni,nj); h_old = zeros(ni,nj);

% initial condition - some random lines
v_initial(3,:) = 1;
h_initial(:,24) = 1;



u_old = u_initial;
v_old = v_initial;
h_old = h_initial;

for t = 1:nt
for i=2:(ni-1)
    for j=2:(nj-1)        
        u(i,j) = 0.25*(u_old(i+1,j) + u_old(i-1,j) + u_old(i,j+1) + u_old(i,j-1)) ...
            - (0.5/2)*(dt/dx) * u_initial(i,j) * (u_old(i+1,j) - u_old(i-1,j)) ...
            - 0.5*(dt/dy) * v_initial(i,j) * (u_old(i,j+1) - u_old(i,j-1)) ...
            - 0.5*g*(dt/dx) * (h_old(i+1,j) - h_old(i-1,j));

        v(i,j) = 0.25*(v_old(i+1,j) + v_old(i-1,j) + v_old(i,j+1) + v_old(i,j-1)) ...
            - (0.5/2)*(dt/dy) * v_initial(i,j) * (v_old(i,j+1) - v_old(i,j-1)) ...
            - 0.5*(dt/dx) * u_initial(i,j) * (v_old(i+1,j) - v_old(i-1,j)) ...
            - 0.5*g*(dt/dy) * (h_old(i,j+1) - h_old(i,j-1));

        h(i,j) = 0.25*(h_old(i+1,j) + h_old(i-1,j) + h_old(i,j+1) + h_old(i,j-1)) ...
            - 0.5*(dt/dx)* u_initial(i,j) * (h_old(i+1,j) - h_old(i-1,j)) ...
            - 0.5*(dt/dy)* v_initial(i,j) * (h_old(i,j+1) - h_old(i,j-1)) ...
            - 0.5*(dt/dx) * h_initial(i,j) * (u_old(i+1,j) - u_old(i-1,j)) ...
            - 0.5*(dt/dy) * h_initial(i,j) * (v_old(i,j+1) - v_old(i,j-1));

    end
end

    u_old = u;
    v_old = v;
    h_old = h;
end

first_h = h;

%% FULL NONLINEAR SOLVE
clear all;

g = 9.8; % units of m/s^2

ni = 100; % number of gridpoints in x
nj = 100; % number of gridpoints in y
nt = 10000; % number of timesteps

dt = 0.01; % timestep
dx = 0.5; % xgrid step
dy = 0.5; % ygrid step

u_initial = zeros(ni,nj); % initialize grids
v_initial = zeros(ni,nj); % initialize grids
h_initial = zeros(ni,nj); % initialize grids
u = zeros(ni,nj); u_old = zeros(ni,nj);
v = zeros(ni,nj); v_old = zeros(ni,nj);
h = zeros(ni,nj); h_old = zeros(ni,nj);

% initial condition - some random lines
v_initial(3,:) = 1;
h_initial(:,24) = 1;



u_old = u_initial;
v_old = v_initial;
h_old = h_initial;

for t = 1:nt
for i=2:(ni-1)
    for j=2:(nj-1)        
        u(i,j) = 0.25*(u_old(i+1,j) + u_old(i-1,j) + u_old(i,j+1) + u_old(i,j-1)) ...
            - (0.5/2)*(dt/dx) * u_old(i,j) * (u_old(i+1,j) - u_old(i-1,j)) ...
            - 0.5*(dt/dy) * v_old(i,j) * (u_old(i,j+1) - u_old(i,j-1)) ...
            - 0.5*g*(dt/dx) * (h_old(i+1,j) - h_old(i-1,j));

        v(i,j) = 0.25*(v_old(i+1,j) + v_old(i-1,j) + v_old(i,j+1) + v_old(i,j-1)) ...
            - (0.5/2)*(dt/dy) * v_old(i,j) * (v_old(i,j+1) - v_old(i,j-1)) ...
            - 0.5*(dt/dx) * u_old(i,j) * (v_old(i+1,j) - v_old(i-1,j)) ...
            - 0.5*g*(dt/dy) * (h_old(i,j+1) - h_old(i,j-1));

        h(i,j) = 0.25*(h_old(i+1,j) + h_old(i-1,j) + h_old(i,j+1) + h_old(i,j-1)) ...
            - 0.5*(dt/dx)* u_old(i,j) * (h_old(i+1,j) - h_old(i-1,j)) ...
            - 0.5*(dt/dy)* v_old(i,j) * (h_old(i,j+1) - h_old(i,j-1)) ...
            - 0.5*(dt/dx) * h_old(i,j) * (u_old(i+1,j) - u_old(i-1,j)) ...
            - 0.5*(dt/dy) * h_old(i,j) * (v_old(i,j+1) - v_old(i,j-1));

    end
end

    u_old = u;
    v_old = v;
    h_old = h;
end

second_h = h;
%%
%s = surf(h,'edgecolor','none');

figure
s = surf(first_h,'edgecolor','none');
figure
s = surf(second_h,'edgecolor','none');
figure
s = surf(first_h-second_h,'edgecolor','none');


