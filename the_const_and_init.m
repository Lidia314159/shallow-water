% constants, initial condition, etc

% Constants
g  = 9.81;
u0 = 0;                               
v0 = 0;
b0  = 0;
h0 = 5030;   

%initialize h0 shape, pick one
h_init_shape = 1;  %=0 prism, =1 Gaussian, in the middle
%boundary condition choice, pick one
bc_choice = 1;     %=0 old code (none), =1 reflective, =2 oubound, =3 periodic 

% Define the x domain
ni = 151;                               
xmax = 100000;                      
dx = xmax/(ni-1);
x  = (0:dx:xmax);

% Define the y domain
nj = 151;
ymax = 100000;                     
dy = ymax/(nj-1);
y  = (0:dy:ymax);

% Define the wavespeed
wavespeed = u0 + sqrt(g*(h0 - b0));

% Define time-domain
dt = 0.68*dx/wavespeed;  % smaller is better, =1 makes bad animation
%tmax = 1500;
%t = [0:dt:tdomain];
tmax = dt*2500;  %old value: 1500
t=(0:dt:tmax);
%courant = wavespeed*dt/dx;

%b 
b = zeros(length(x), length(y));

% we ignore the terrain profile
% Define b: 
% for i = 1:length(x)
%     if x(i) > 20001
%         b(:,i) = 0;
%     elseif x(i) < 20000
%         b(:,i) = 5000/20000*(20000-x(i));
%     end
% end