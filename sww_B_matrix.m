%From Steven McHale: Tsunami Model
%Shallow-Water Wave Equation
%Crank-Nicholson Discretization

%%
clear
close
%clf;
%clc;

the_globals
%the_parameters
the_const_and_init

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

% initialize bondary conditions
[u_new, v_new, h_new] = bc_init(u, v, h);
u = u_new;
v = v_new; 
h = h_new;

% -----------------

% Employ Lax-Friedrichs method or Crank-Nicolson

% init 
%u;
%v;
%h;

%% Build B-matrix: from t=500 to 600 and from t=500 to 1500
T1 = 1000;
Tmax = 2000;

saveBu = zeros(ni, nj, Tmax-T1);
saveBv = zeros(ni, nj, Tmax-T1);
saveBh = zeros(ni, nj, Tmax-T1);

for index = 1:T1
    [u_new, v_new, h_new] = sww(u, v, h);
    % old=new and continue loop
    u = u_new;
    v = v_new;
    h = h_new;
end

for index = 1:(Tmax-T1)
    [u_new, v_new, h_new] = sww(u, v, h);
    %
    saveBu(:,:,index) = u_new;
    saveBv(:,:,index) = v_new;
    saveBh(:,:,index) = h_new;

    % old=new and continue loop
    u = u_new;
    v = v_new;
    h = h_new;
end

Bu = zeros(ni, nj);
Bv = Bu; 
Bh = Bu;

temp = zeros(1,Tmax-T1);

for i=1:ni
    for j=1:nj
        temp = squeeze(saveBu(i,j,:));
        Bu(i,j) = cov(temp);
        temp = squeeze(saveBv(i,j,:));
        Bv(i,j) = cov(temp);
        temp = squeeze(saveBh(i,j,:));
        Bh(i,j) = cov(temp);
    end
end

figure
ss = surf(Bu, 'FaceAlpha',0.5);
set(ss,'LineStyle','none')
title('Covariance Bu in time from t=1000 to 2000')

figure
ss = surf(Bv, 'FaceAlpha',0.5);
set(ss,'LineStyle','none')
title('Covariance Bu in time from t=1000 to 2000')

figure
ss = surf(Bh, 'FaceAlpha',0.5);
set(ss,'LineStyle','none')
title('Covariance Bu in time from t=1000 to 2000')

%load train
%sound(y, Fs)
