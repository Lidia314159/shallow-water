%Steven McHale
%Tsunami Model
%Shallow-Water Wave Equation
%Crank-Nicholson Discretization

diary('tsunami.out')

clear;
clf;
clc;

% Constants
g  = 9.81;
u0 = 0.0;
v0 = 0.0;
b  = 0.0;
h_level = 5000.0;
h0 = 5030.0;

% Define the x domain
ni = 151;
xmax = 100000.0;
dx = xmax/(ni-1.0);
x  = [0:dx:xmax];

% Define the y domain
nj = 151;
ymax = 100000.0;
dy = ymax/(nj-1.0);
y  = [0:dy:ymax];

% Define the wavespeed
wavespeed = u0 + sqrt(g*(h0 - b));

% Define time-domain
dt = 0.68*dx/wavespeed;
tmax = 1500.0;
t=[0:dt:tmax];
%courant = wavespeed*dt/dx;

% Build empty u, v, b matrices
u=zeros(length(x), length(y));
v=zeros(length(x), length(y));
b=zeros(length(x), length(y));

% Define h
h=zeros(length(x), length(y));
h(:,:) = h_level;
xc = floor(max(x) / 2);
yc = floor(max(y) / 2);
sigma = floor(max(x) / 20);
for i = 1:ni
  for j = 1:nj
      dsqr = ((x(i) - xc)^2 + (y(j) - yc)^2);
      h(i,j) = h_level + exp(-dsqr / sigma^2) * (h0 - h_level);
  end
end

dtdx = dt / dx;
dtdy = dt / dy;

figure(1)

% Employ Lax
for n=1:(length(t)-1)

%    mesh(x,y,h(:,:))
%    axis ([0 100000 0 100000 4990 5030])
%    title ('Gaussian Pulse with Reflective Boundaries')
%    xlabel('X Domain [m]')
%    ylabel('Y Domain [m]')
%    zlabel('Height [m]')
%    pause(0.02)

    fprintf('%7d t=%15.6f u(75,75) = %9.6f v(75,75) = %9.6f h(75,75) = %11.6f\n', n-1, t(n), u(75,75), v(75,75), h(75,75))

    % Update boundary conditions: reflective
    h_new(:,1)  = h(:, 2);
    h_new(:,nj) = h(:,nj-1);
    h_new(1,:)  = h(2,:);
    h_new(ni,:) = h(ni-1,:);

    u_new(:,1)  =  u(:,2);
    u_new(:,nj) =  u(:,nj-1);
    u_new(1,:)  = -u(2,:);
    u_new(ni,:) = -u(ni-1,:);

    v_new(:,1)  = -v(:,2);
    v_new(:,nj) = -v(:,nj-1);
    v_new(1,:)  =  v(2,:);
    v_new(ni,:) =  v(ni-1,:);

    for i=2:(ni - 1)
        for j=2:(nj - 1)
            u_new(i,j) = ((u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)) / 4)  ...
                       - 0.5 * dtdx * ((u(i+1,j)^2) / 2 - (u(i-1,j)^2) / 2) ...
                       - 0.5 * dtdy * (v(i,j)) * (u(i,j+1) - u(i,j-1))      ...
                       - 0.5 * g * dtdx * (h(i+1,j) - h(i-1,j));

            v_new(i,j) = ((v(i+1,j) + v(i-1,j) + v(i,j+1) + v(i,j-1)) / 4)  ...
                       - 0.5 * dtdy * ((v(i,j+1)^2) / 2 - (v(i,j+1)^2) / 2) ...
                       - 0.5 * dtdx * (u(i,j)) * (v(i+1,j) - v(i-1,j))      ...
                       - 0.5 * g * dtdy * (h(i,j+1) - h(i,j-1));

            h_new(i,j) = ((h(i+1,j) + h(i-1,j) + h(i,j+1) + h(i,j-1)) / 4)                       ...
                       - 0.5 * dtdx * (u(i,j)) * ((h(i+1,j) - b(i+1,j)) - (h(i-1,j) - b(i-1,j))) ...
                       - 0.5 * dtdy * (v(i,j)) * ((h(i,j+1) - b(i,j+1)) - (h(i,j-1) - b(i,j-1))) ...
                       - 0.5 * dtdx * (h(i,j) - b(i,j)) * (u(i+1,j) - u(i-1,j))                  ...
                       - 0.5 * dtdy * (h(i,j) - b(i,j)) * (v(i,j+1) - v(i,j-1));
        end
    end

    u = u_new;
    v = v_new;
    h = h_new;

end

quit;