%From Steven McHale: Tsunami Model
%Shallow-Water Wave Equation
%Crank-Nicholson Discretization

%%
clear
close
%clf;
%clc;

the_globals
the_parameters
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
%
if check_plot_h0  %check the look of the init h0
    %figure
    mesh(x, y, h)
    axis ([0 100000 0 100000 4970 5030])
    title ('Shallow Water Wave - SWW')
    xlabel('X Domain [m]')
    ylabel('Y Domain [m]')
    zlabel('Height [m]')
end

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

%% Building TL

if build_TL
    %advance the model: get the sww wave state at T
    T = floor(length(t)*0.25);
    %
    % run sww from time=1 to time=T+1, save off at T and at T+1
    for index=1:T+1  
        [u_new, v_new, h_new] = sww(u, v, h);
        % old=new and continue loop
        u = u_new;
        v = v_new;
        h = h_new;
        %
        if index==T
            u_T = u;
            v_T = v;
            h_T = h;
        end
    end
    u_T1 = u;
    v_T1 = v;
    h_T1 = h;
    % 
    mesh(x, y, (h_T1 - h_T))
    %
    %get TL of the state (u, v, h) at time t=T
    [u_Tnew, v_Tnew, h_Tnew] = sww_tl(u_T, v_T, h_T);
    %
    mesh(x, y, (h_Tnew - h_T1));
    %
    [u_T1new, v_T1new, h_T1new] = sww_adj(u_T1, v_T1, h_T1);
    %
    mesh(x, y, (h_T1new - h_T));
    %
end

%% Testing RMSE of the TL and ADj

if rmse_TL_ADj
    %rmse_tl_u = zeros(1, length(t)-2);
    %rmse_tl_v = rmse_tl_u;
    rmse_tl_h = zeros(1, length(t)-1);
    %rmse_adj_u = rmse_tl_u;
    %rmse_adj_v = rmse_tl_u;
    rmse_adj_h = zeros(1, length(t)-1);
    
    %init
    u_tm1 = u;
    v_tm1 = v; 
    h_tm1 = h;  %move from t=1 to t=2
    %
    for index=2:length(t)-1     %loop over all time
        %
        [u_t, v_t, h_t] = sww(u_tm1, v_tm1, h_tm1); %tm1 == time minus one
        [u_tp1, v_tp1, h_tp1] = sww(u_t, v_t, h_t); %tp1 == time plus one
        %
        [u_TL, v_TL, h_TL] = sww_tl(u_t, v_t, h_t);
        [u_ADj, v_ADj, h_ADj] = sww_adj(u_t, v_t, h_t);
        %    
        %rmse_tl_u(index-1) = sqrt(sum(dot(u_tp1 - u_TL, u_tp1 - u_TL)));
        %rmse_tl_v(index-1) = sqrt(sum(dot(v_tp1 - v_TL, v_tp1 - v_TL)));
        rmse_tl_h(index-1) = sqrt(sum(dot(h_tp1 - h_TL, h_tp1 - h_TL))/(ni*nj));
        %
        %rmse_adj_u(index-1) = sqrt(sum(dot(u_tm1 - u_ADj, u_tm1 - u_ADj)));
        %rmse_adj_v(index-1) = sqrt(sum(dot(v_tm1 - v_ADj, v_tm1 - v_ADj)));
        rmse_adj_h(index-1) = sqrt(sum(dot(h_tm1 - h_ADj, h_tm1 - h_ADj))/(nj*nj));
            
        %move time by one time step
        u_tm1 = u_t;
        v_tm1 = v_t;
        h_tm1 = h_t;
        %
    end
    
    figure
    plot(rmse_tl_h)
    hold on; grid on
    plot(rmse_adj_h)
    %
end
% end RMSE code <<


%% Write sww state to files

if write_to_files
    cd '/Users/Lidija/Documents/MATLAB/data'
    for index=1:length(t)-1 %plotting loop
        % 
        [u_new, v_new, h_new] = sww(u, v, h);
        %
        save_sww = zeros(size(u,1), size(u,2), 3);
        save_sww(:,:,1) = u_new;
        save_sww(:,:,2) = v_new;
        save_sww(:,:,3) = h_new;
        %
        %write to file
        save_file_name = ['save_sww_' num2str(index) '.mat'];
        save(save_file_name,'save_sww');
        %
        txt_file_name = ['sww_' num2str(index) '.txt'];
        fid = fopen(txt_file_name, 'w+');  % or wt
        %
        fprintf(fid, '%s\n', 'u');
        for ii = 1:size(u,1)
            fprintf(fid, '%5.12f ', u(ii,:));   % '%g ' 
            fprintf(fid, '\n');
        end
        fprintf(fid, '%s\n', 'v');
        for ii = 1:size(v,1)
            fprintf(fid, '%5.12f ', v(ii,:));   % '%g ' 
            fprintf(fid, '\n');
        end
        fprintf(fid, '%s\n', 'h');
        for ii = 1:size(h,1)
            fprintf(fid, '%5.12f ', h(ii,:));   % '%g ' 
            fprintf(fid, '\n');
        end
        fclose(fid);
        %
        %update: old=new and continue loop
        u = u_new;
        v = v_new;
        h = h_new;
        %
    end
    cd '/Users/Lidija/Documents/MATLAB'
end


%% Save off one sww state at time T to a file and restart sww

if save_off_state
    T = floor(length(t)*0.6);
    %write ctate to a file, note: add time stamp to the file name
    saved_file = 'save_sww.mat';
    [saved_file, T] = sww_save(saved_file, T);
    %read state from a file 
    [u, v, h, T] = sww_read(saved_file, T);
    %
    if check_save_off
        %check from time=T+1, new init: [u, v, h, T]
        for index=(T+1):(length(t)-1)
            [u, v, h] = sww(u, v, h); 
            mesh(x,y,h);
            axis ([0 100000 0 100000 4970 5030])
            title ('Shallow Water Wave - SWW')
            xlabel('X Domain [m]')
            ylabel('Y Domain [m]')
            zlabel('Height [m]')
            pause(0.02)
            % old=new and continue loop
            u_old = u;
            v_old = v;
            h_old = h;
        end   
    end %check_save_off
end %save_off_state


%% Making the video

if make_video
    writerObj = VideoWriter('ShallowWater_1.avi');
    writerObj.FrameRate = 30;
    open(writerObj); 
    %
    for index=1:length(t)-1 %plotting loop
        [u_new, v_new, h_new] = sww(u, v, h);
        %either 
        %mesh(x,y,h);
        %or
        ss = surf(x,y,h, 'FaceAlpha',0.5);
        set(ss,'LineStyle','none')
        axis ([0 100000 0 100000 4970 5030])
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
end


%% General check with plotting -- if it works

if check_plot
    figure
    %
    for index=1:length(t)-1 %plotting loop
    %
        [u_new, v_new, h_new] = sww(u, v, h);
        %either 
        %mesh(x,y,h);
        %or
        ss = surf(x,y,h, 'FaceAlpha',0.5);
        set(ss,'LineStyle','none')
        %
        axis ([0 100000 0 100000 4970 5030])
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
end

%load train
%sound(y, Fs)
