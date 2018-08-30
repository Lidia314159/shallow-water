% collecting data

% block
index = [500, 1000, 1500, 2000, 5000, 10000];
h_max = [5.0044e+03, 5.0045e+03, 5.0028e+03, 5.0028e+03, 5.0031e+03, 5.0006e+03];
h_min = [4.9976e+03, 4.9977e+03, 4.9965e+03, 4.9970e+03, 4.9997e+03, 4.9997e+03];
h_mean = [5.0003e+03, 5.0004e+03, 5.0003e+03, 5.0003e+03, 5.0004e+03, 5.0003e+03];

% Gauss
indexG = [50, 100, 500, 1000, 1500, 2000, 2500, 5000, 10000];
h_maxG = [5.0041e+03, 5.0038e+03, 5.0029e+03, 5.0029e+03, 5.0019e+03, 5.0019e+03, 5.0015e+03, 5.0020e+03, 5.0004e+03];
h_minG = [4.9975e+03, 4.9985e+03, 4.9985e+03, 4.9986e+03, 4.9977e+03, 4.9980e+03, 4.9995e+03, 4.9998e+03, 4.9998e+03];
h_meanG = [5.0002e+03, 5.0003e+03, 5.0002e+03, 5.0003e+03, 5.0002e+03, 5.0002e+03, 5.0002e+03, 5.0002e+03, 5.0002e+03];

figure

plot(index, h_max-5000, '-gs', indexG, h_maxG-5000, '-ks', ...
    index, h_min-5000, '-yo', indexG, h_minG-5000, '-bo', ...
    index, h_mean-5000, '-m^', indexG, h_meanG-5000, '-r^')
grid on
title('hmax, hmin, hmean with hInitBlock vs. hInitGauss')
xlabel('Number of time steps')
ylabel('(hmax, hmin, hmean) - h0, where h0=5000')
legend('hmaxBlock', 'hmaxGauss', 'hminBlock', 'hminGauss', 'hmeanBlock', 'hmeanGauss')


