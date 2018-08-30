function h = h_start(h_init)
%initialize h 
% choice=0, original code, column 
% choice=1, Gaussian pulse

the_globals
the_parameters
the_const_and_init

h = h_init;

h_mat = h_init;
h_level = 5000;
h0 = 5030;
H = h0 - h_level;

switch h_init_shape
    case 0
    % old code: a column in the middle
    h(:,:) = h_level;                              
    mid_pt = floor(ni/2);
    mid_start = mid_pt - floor(15/2);
    mid_end = mid_pt + ceil(15/2);
    h(mid_start:mid_end, mid_start:mid_end) = h0;

    case 1
    % 2D Gaussian in the middle
    xc = floor(max(x)/2);
    yc = floor(max(y)/2);
    %sig = floor(max(x)/20);
    sig = floor(max(x)/40);

    % Gaussian at mid-point 
    for i=2:(ni-1)
        for j=2:(nj-1)
            expon = ((x(i)-xc)^2 + (y(j)-yc)^2) / sig^2;
            h_mat(i,j) = exp(-expon);      
        end
    end
    h = h_mat*(h0 - h_level) + h_level;

    %matlab tricks that work
%     [X, Y] = meshgrid(x, y);
%     hh = exp(-1/sig^2*((X-xc).^2 + (Y-yc).^2));
%     h = h0*hh;
end

