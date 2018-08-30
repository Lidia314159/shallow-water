function [saved_file, T] = sww_save(saved_file, T)

%Save off a sww state at time T (time step) to a file and restart sww
%T = floor(length(t)*0.25);
%write to file: 'save_sww'

the_globals

for index=1:T  % run sww from time=1 to time=T
    [u_new, v_new, h_new] = sww();
    % old=new and continue loop
    u = u_new;
    v = v_new;
    h = h_new;
end
%
%save off to file save_sww.mat
% size(u)
% size(v)
% size(h)

save_sww = zeros(size(u,1), size(u,2), 3);

save_sww(:,:,1) = u;
save_sww(:,:,2) = v;
save_sww(:,:,3) = h;

save('save_sww.mat','save_sww');


