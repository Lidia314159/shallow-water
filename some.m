%% Define Function and Variables
clear, close all
clc
c_0 = 5;
D_ab = 1E-4;
num_x = 100;
num_y = 100;
num_t = 100;
t_max = 50000;
x_range = linspace(-5,5,num_x);
y_range = linspace(-5,5,num_y);
t_range = linspace(1,t_max,num_t);
c_val = zeros(num_x,num_y,num_t);
for i=1:num_x;
    for j=1:num_y
        for t=1:num_t
            c_val(i,j,t)=(c_0/sqrt(pi*D_ab*t_range(t)))*exp(-((x_range(i).^2)+(y_range(j).^2))/(4*D_ab*t_range(t)));
        end
    end
end
%% Record Movie
writerObj = VideoWriter('myfirst3Dmovie.avi');
writerObj.FrameRate = 30;
open(writerObj); 
for t=1:num_t    
    mesh(x_range,y_range,c_val(:,:,t));
      xlim([-5,5])
      ylim([-5,5])
      zlim([0,15])
      xlabel('Position-X [cm]')
      ylabel('Position-Y [cm]')
      zlabel('Concentration [mol / m^2 ]')
      title('Diffusion of a point source over time')
      frame = getframe(gcf); 
      writeVideo(writerObj, frame);
   end
close(writerObj);
