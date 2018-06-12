clearvars

%%%%%%%%%%%%%%%%%
wheel_base = 118.11; %3 meters in inches
%EQUATION FOR STEERIN COLUMN
hole_r = [0.4, 2.5440+.12, 1.59+.03]'; %VARY THIS

%rotate around vector V
V = [0 -5.92274, 9.47638]';

below = -1.15620; %OR VARY THIS, how far attachments are below holes
attach_r = (below.*V)./norm(V)+hole_r; 

%%%%%%%%%%%%%%%%%
%EQUATIONS FOR CONTROL ARMS
p_right = [7.44904, 1.73588, 0.75816]';
c_right = [6.979, -2.60568, 0.85416]';

z_height = p_right(3);

point_list_rx =[];
point_list_ry =[];
theta_list = [];

for iter = linspace(-pi/2,pi/2,5000)
    Rot_z_r = [cos(iter) -sin(iter) 0;...
     sin(iter) cos(iter) 0;...
     0 0 1];
 
    pr = Rot_z_r*(p_right-c_right)+c_right;
    
    point_list_rx = [point_list_rx pr(1)];
    point_list_ry = [point_list_ry pr(2)];
    theta_list = [theta_list iter];
end

theta_list(theta_list>pi) = theta_list(theta_list>pi) - 2*pi;
point_list_mat = [point_list_rx' point_list_ry' theta_list'];

%%%%%%%%%%%%%%%%%
length_r = norm(p_right-attach_r);

%find locus of turning radii
steering_list = linspace(-pi/4,pi/4,5000);
turning_mat = [];

for iter = steering_list
    [attach_2] = rotate_column(V,attach_r, iter);
    wait = iter;
    
    [new_r, theta] = minimize_distance(attach_2, length_r, point_list_mat,z_height);
    turning_mat = [turning_mat; theta];
end

turning_radii = wheel_base./sin(turning_mat);

plot(steering_list, turning_mat,'.','MarkerSize',5)
grid on

function [new_r] = rotate_column(V, hole_r, steer_theta)
     %rotate so vector V is along z axis
     theta = pi/365;
     theta_sum = 0;
     Rot_x = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
     while(abs(V(2)) > 0.1)
        V = Rot_x*V;
        hole_r = Rot_x*hole_r;
        theta_sum = theta_sum + theta;
     end

     %rotate around transformed steering column
     Rot_z = [cos(steer_theta) -sin(steer_theta) 0;...
         sin(steer_theta) cos(steer_theta) 0;...
         0 0 1];
     hole_r = Rot_z*hole_r;

     %inverse previous x rotation
     Rot_x = [1 0 0; 0 cos(-theta_sum) -sin(-theta_sum); 0 sin(-theta_sum) cos(-theta_sum)];
     new_r = Rot_x*hole_r;
end
 
function [new_r, theta] = minimize_distance(hole_r, length_r, point_list_mat,z_height)
    %1,2 are left brake mount, 3,4 are right brake mount
    tol_r = zeros(5000,1);
    
    z = zeros(5000,1);
    point = [point_list_mat(:,1:2) z];
    test = hole_r' - point;
    for i = 1:5000
        test3(i) = norm(test(i,:));
    end
    test3 = test3';
    tol_r = abs(test3 - length_r);
    
    
    [val,min_r] = min(tol_r);
    
    new_r = [point_list_mat(min_r,1) point_list_mat(min_r,2) z_height];
    theta = [point_list_mat(min_r,3)];
end