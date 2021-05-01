%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MATLAB script that calculates the time-varying load in the shock 
%   chord using Newton's second law. The rocket is modelled as a point mass and
%   the shock cord is inextensible. The variation of the parachute force as it 
%   opens is modelled by a polynomial opening period followed by a constant 
%   force, with the parachute opening time as a input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

file_name = 'nose_click_test1';        

    %   inputs
Cd = 1.4;           %   parachute drag coefficient
L = 3;              %   shock cord length                               [m]
S = pi*0.305^2;     %   parachute area                                  [m2]
m = 0.9;            %   rocket mass                                     [kg]
t_delay = 2.6;      %   parachute delployment delay (after apogee)      [s]

    %   simulation time parameters
dt = 0.01;
total_time = 5;
sim_time = 0;

V0 = sqrt(2*9.81*L+t_delay^2*9.81^2);   
v = [];
v(1) = -V0;
v(2) = V0+9.81*dt;
k = 0.5*1.225*Cd*S/m;
grad = [];
grad(1:2) = 0;

error = 10;

    %   parachute force scale factor
F = ones(1,total_time/dt+3);
F_open = 0.3+dt;                  %   time taken to open parachute fully [s]

F_ticks = F_open/dt;

disp('Assume parachute opening can be modelled by an nth order polynomial')
n = input('n = ');

if n ~= 0
    
    for i = 1:F_ticks
   
        F(i) = (i*dt)^n/(F_open)^n;
    
    end
end

i = 2;

    %   solver loop
while sim_time<total_time+dt
    v(i) = (k*v(i-1)^2*F(i)-9.81)*dt+v(i-1);
    %error = (v(i)-v(i-1))^2;
    sim_time = sim_time+dt;
    i = i+1;
    
end

time = 0:dt:(length(v)-1)*dt;

    %   plot to show parachute force scaling
figure
hold on
title('Chute Opening Profile')
plot(time,F)
xlabel('time / s')
ylabel('Chute diameter / percentage of max')
ylim([0 1.2])
xlim([0 F_open*1.2])
grid on
grid minor
box on
hold off
exportgraphics(gcf, sprintf('%s_chute_opening.png',  file_name), 'Resolution', 600);

for i = 1:length(v)-1
    grad(i) = (v(i+1)-v(i))/dt;
    
end

figure
plot(time,v)
title('Descent Velocity After Parachute Deployment')
xlabel('time / s')
ylabel('Velocity / ms^{-1}')
grid minor
grid on
exportgraphics(gcf, sprintf('%s_velocity.png', file_name), 'Resolution', 600);

figure
plot(time(1:length(time)-1),m*grad)
title(['Shock Load due to Parachute Deployment with t_{free} = ', num2str(t_delay), 's'])
grid minor
grid on
text(total_time*0.3,max(m*grad)/2,['Max shock load = ', num2str(double(max(m*grad))), ' N'])
text(total_time*0.3,max(m*grad)/2.4,['Average load = ', num2str(sum(m*grad)*dt/total_time), ' N'])
xlabel('time / s')
ylabel('Force / N')
exportgraphics(gcf, sprintf('%s_shock_loads.png',  file_name), 'Resolution', 600);

disp(['Total impulse delivered = ', num2str(sum(grad)*dt*m), ' Ns'])
disp(['Max shock load = ', num2str(max(m*grad)), ' N'])
disp(['Max acceleration = ', num2str(max(grad)), ' m/s2'])
disp(['Average load = ', num2str(sum(m*grad)*dt/total_time), ' N'])
