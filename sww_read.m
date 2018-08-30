function [u_new, v_new, h_new, tt] = sww_read(saved_file, T)

%saved_file

%read file save_sww.mat, reconstruct state [u,v,h,tt]
reconstruct = matfile(saved_file);
ssaved = reconstruct.save_sww;
%
u_new = ssaved(:,:,1);
v_new = ssaved(:,:,2);
h_new = ssaved(:,:,3);
tt = T;

