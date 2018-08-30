% make sw_video

%% Making the video

clear
close

the_globals
the_parameters

% ----------
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
ni = 301; %151;                               
xmax = 200000;                      
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

% ----------

% Build empty u, v, h, b matrices
u_init = zeros(length(x), length(y));
v_init = zeros(length(x), length(y));
h_init = zeros(length(x), length(y)); 
%b = zeros(length(x), length(y));

% Initialize u, v
u = u_init;
v = v_init;
% Initialize h: h_init_shape=0 (column), h_init_shape=1 (Gaussian) 
h = h_start(h_init);
%
%check the look of the init h0
%figure
mesh(y, x, h) % mesh requires y first
axis ([0 200000 0 100000 4970 5030])
title ('Shallow Water Wave - SWW')
ylabel('X Domain [m]')
xlabel('Y Domain [m]')
zlabel('Height [m]')


% initialize bondary conditions
[u_new, v_new, h_new] = bc_init(u, v, h);
u = u_new;
v = v_new; 
h = h_new;

%plot for check
figure
%
for index=1:length(t)-1 %plotting loop
%
    [u_new, v_new, h_new] = sww(u, v, h);
    %either 
    %mesh(x,y,h);
    %or
    ss = surf(y,x,h, 'FaceAlpha',0.5);
    set(ss,'LineStyle','none')
    %
    axis ([0 200000 0 100000 4970 5030])
    title ('Shallow Water Wave - SWW')
    xlabel('X Domain [m]')
    ylabel('Y Domain [m]')
    zlabel('Height [m]')
    pause(0.02)
    %Check statistics of hmin, hmax, hmean
    if save_stats
        if ismember(index, [50, 100, 500, 1000, 1500, 2000, 2500, 5000, 10000])
            disp('index')
            index
            disp('h_max')
            max(max(h))
            disp('h_min')
            min(min(h))
            disp('h_mean')
            mean(mean(h))
        end
    end
    % old=new and continue loop
    u = u_new;
    v = v_new;
    h = h_new;

end % end plotting loop

% -----------------

writerObj = VideoWriter('ShallowWater_2.avi');
writerObj.FrameRate = 30;
open(writerObj); 
%
for index=1:length(t)-1 %plotting loop
    [u_new, v_new, h_new] = sww(u, v, h);
    %either 
    %mesh(x,y,h);
    %or
    ss = surf(y,x,h, 'FaceAlpha',0.5);  %surf is like mesh, (y,x,h)
    set(ss,'LineStyle','none')
    axis ([0 200000 0 100000 4970 5030])
    title ('Shallow Water Wave - SWW')
    xlabel('X Domain [m]')
    ylabel('Y Domain [m]')
    zlabel('Height [m]')
    %pause(0.01)
    %
    frame = getframe(gcf); 
    writeVideo(writerObj, frame);
    %
    % old=new and continue loop
    u = u_new;
    v = v_new;
    h = h_new;
end
close(writerObj);
    %

