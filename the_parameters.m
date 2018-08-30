% the test flags&parameters: flags used for testing, plotting, etc

%turn on for plotting h0: =0 no plot, =1 do plot
check_plot_h0 = 0;   %for plotting h0

build_TL = 0;         %TL and ADj code
rmse_TL_ADj = 0;      %plot RMSE of TL and ADj 

% save off one file 
save_off_state = 0;   %save off a sww state to a file for restart, =0 skip, =1 do it
check_save_off = 0;   %plot-check by plotting from t=T onward

%write sww states to files, for all t
write_to_files = 0; 

make_video = 0;   %=1 only if making the video

check_plot = 0;   %plot wave in time
save_stats = 0;   %saves off stats: hmin, hmax, hmean




