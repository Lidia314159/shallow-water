% Printing RMSE and timing of model, TL, and ADj

%%
clear
close

the_globals
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

% if check_plot_h0  %check the look of the init h0
%     %figure
%     mesh(x, y, h)
%     axis ([0 100000 0 100000 4970 5030])
%     title ('Shallow Water Wave - SWW')
%     xlabel('X Domain [m]')
%     ylabel('Y Domain [m]')
%     zlabel('Height [m]')
% end

% initialize bondary conditions
[u_new, v_new, h_new] = bc_init(u, v, h);
u = u_new;
v = v_new; 
h = h_new;

%% Testing RMSE of the TL and ADj

%rmse_tl_u = zeros(1, length(t)-2);
%rmse_tl_v = rmse_tl_u;

loww = 1;  %1500
highh = 1500; %5000

rmse_tl_h = zeros(1, highh);
%rmse_adj_u = rmse_tl_u;
%rmse_adj_v = rmse_tl_u;
rmse_adj_h = zeros(1, highh);

% %init
% for index=1:loww
%     [u_new, v_new, h_new] = sww(u, v, h);
%     %update, cont loop
%     u = u_new;
%     v = v_new;
%     h = h_new;
% end

u_tm1 = u;
v_tm1 = v; 
h_tm1 = h;  %advance model from t=1 to t=loww

delta_hTL_min = zeros(1, highh);
delta_hTL_max = zeros(1, highh);
delta_hADj_min = zeros(1, highh);
delta_hADj_max = zeros(1, highh);

for index=1:highh     %loop t=loww to t=highh
    
    [u_t, v_t, h_t] = sww(u_tm1, v_tm1, h_tm1); %tm1 == time minus one
    
    [u_tp1, v_tp1, h_tp1] = sww(u_t, v_t, h_t); %tp1 == time plus one
    
    [u_TL, v_TL, h_TL] = sww_tl(u_t, v_t, h_t);
    
    %
    [u_ADj, v_ADj, h_ADj] = sww_adj_fast(u_t, v_t, h_t);
 
    %rmse_tl_u(index-1) = sqrt(sum(dot(u_tp1 - u_TL, u_tp1 - u_TL)));
    %rmse_tl_v(index-1) = sqrt(sum(dot(v_tp1 - v_TL, v_tp1 - v_TL)));
    rmse_tl_h(index) = sqrt(sum(dot(h_tp1 - h_TL, h_tp1 - h_TL))/(ni*nj)); 
    delta_hTL_max(index) = max(max(h_tp1 - h_TL));
    delta_hTL_min(index) = min(min(h_tp1 - h_TL));
    %
    %rmse_adj_u(index-1) = sqrt(sum(dot(u_tm1 - u_ADj, u_tm1 - u_ADj)));
    %rmse_adj_v(index-1) = sqrt(sum(dot(v_tm1 - v_ADj, v_tm1 - v_ADj)));
    rmse_adj_h(index) = sqrt(sum(dot(h_tm1 - h_ADj, h_tm1 - h_ADj))/(ni*nj));
    delta_hADj_max(index) = max(max(h_tm1 - h_ADj));
    delta_hADj_min(index) = min(min(h_tm1 - h_ADj));
    
    %move time by one time step
    u_tm1 = u_t;
    v_tm1 = v_t;
    h_tm1 = h_t;
    %
end

figure
plot(loww:1:highh, rmse_tl_h)
hold on; grid on
plot(loww:1:highh, rmse_adj_h)
title('RMSE of TL vs State(k+1) and ADj vs State(k-1)')
xlabel('time [k]')
ylabel('RMSE')
legend('RMSE TL', 'RMSE ADj')
%


figure
plot(loww:1:highh, delta_hTL_max)
hold on; grid on
plot(loww:1:highh, delta_hTL_min)
plot(loww:1:highh, delta_hADj_max)
plot(loww:1:highh, delta_hADj_min)
title('Max and Min diff: h(TL) vs h(k+1) and h(ADj) vs h(k-1)')
xlabel('time [k]')
ylabel('Max and Min diff of h')
legend('max h TL', 'min h TL', 'max h ADj', 'min h ADj')

%% timing

timer_model = 0; 
%timer_tl = 0;
%timer_adj_fast = 0;

% for cc = 1:200
% 
%     f1 = @() sww(u_t, v_t, h_t);
%     timer_model = timer_model + timeit(f1);  %in sec
% 
%     f2 = @() sww_tl(u_t, v_t, h_t);
%     timer_tl = timer_tl + timeit(f2); 
% 
%     f3 = @() sww_adj_fast(u_t, v_t, h_t);
%     timer_adj_fast = timer_adj_fast + timeit(f3);
% 
% end

tstart = tic;
for cc = 1:100000
    
    %[u_n, v_n, h_n] = sww(u, v, h);
    %[u_n, v_n, h_n] = sww_tl(u, v, h);
    %[u_n, v_n, h_n] = sww_adj_fast(u, v, h);
    [u_n, v_n, h_n] = sww_adj(u, v, h);
end
tend = toc(tstart);


fprintf('\n %15.8f', tend/cc)
%fprintf('%15.8f', timer_tl/cc)
%fprintf('%15.8f', timer_adj_fast/cc)

% end of the RMSE and timing tests <<

