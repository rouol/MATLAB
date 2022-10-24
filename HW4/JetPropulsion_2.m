%Task of hitting a target with a rocket. But now the missile has a jet engine.
%Find the angle at which the missile will hit the target with the given delta distance.
%Realize the problem using the firing method.

%Set parameters
g = 9.81; %acceleration of gravity
k = 0.5; %air drag coefficient
m = 100; %the mass of the missile
M0 = 10; %initial mass of fuel
FCR = 1; %fuel consumption factor
F_jet = 100; %engine thrust
delta = 5; %delta distance
L = 150; %distance to target
v0 = 50; % initial speed
start_point = [0,0]; %start_point
end_point = [L,0]; %end_point

% r_jet = GetJetTrajectory(start_point, v0, pi/3, M0, F_jet, FCR, g, k, m);
% r = GetTrajectory(start_point, v0, pi/3, g, k, m);
% plot(r(:,1), r(:,2), 'b', r_jet(:,1), r_jet(:,2), 'r');

%Solve the problem
% The system of equations of motion of the rocket with allowance for air resistance, gravity, and jet thrust has the form:
%(m+M(t))*x'(t) + k*x'(t) = F_jet_x
%(m+M(t))*y''(t) + k*y'(t) + g*(m+M(t)) = F_jet_y
%M'(t) = FCR (Fuel Consumption Rate)
% Let's reduce the system to a first-order system:
%x'(t) = v_x(t)
%y'(t) = v_y(t)
%v_x'(t) = (F_jet_x - k*v_x(t))/(m+M(t))
%v_y'(t) = (F_jet_y - k*v_y(t))/(m+M(t)) - g
%M'(t) = FCR
%Initial conditions:
%M(0) = M0
%x(0) = start_point(1)
%y(0) = start_point(2)
%v_x(0) = v0*cos(alpha)
%v_y(0) = v0*sin(alpha)
% Using the firing method, find the angle alpha at which the missile hits the target at a given delta distance.
%We will use Newton's method to find the angle alpha.

%Implementation of the firing method

%Set initial values
alpha_min = 0;
alpha_max = pi/2;
alpha = NaN;
distance = NaN;
r = [];
is_found = false;
N = 0;

%Till we find the exact angle alpha, we change it to alpha_step
while true
    %Increase the number of iterations
    N=N+1;

    % Find two rocket trajectories.
    r1 = GetJetTrajectory(start_point, v0, (alpha_max+2*alpha_min)/3, M0, F_jet, FCR, g, k, m);
    r2 = GetJetTrajectory(start_point, v0, (2*alpha_max+alpha_min)/3, M0, F_jet, FCR, g, k, m);

    %Determination of the distances between the end point and the hit point of the missile
    distance1 = abs(r1(end,1) - end_point(1));
    distance2 = abs(r2(end,1) - end_point(1));

    %If the distance is less than the delta distance, then the angle alpha is found
    if distance1 < delta
        alpha = (alpha_max+2*alpha_min)/3;
        distance = distance1;
        r = r1;
        is_found = true;
        break
    end
    if distance2 < delta
        alpha = (2*alpha_max+alpha_min)/3;
        distance = distance2;
        r = r2;
        is_found = true;
        break
    end

    %If both distances are greater than the delta distance, then change the boundaries of the angle alpha
    if distance1 < distance2
        alpha_max = (2*alpha_max+alpha_min)/3;
    else
        alpha_min = (alpha_max+2*alpha_min)/3;
    end

    %If the difference between the bounds of the angle alpha is less than the specified accuracy, then the angle alpha is not found
    if abs(alpha_min - alpha_max) < 0.001
        alpha = (alpha_max + alpha_min)/2;
        distance = abs(r1(end,1) - end_point(1));
        r = r1;
        is_found = false;
        break
    end
end

%Result display
if ~is_found
    disp('Alpha not found, but the closest values:')
end
disp(['Alpha = ', num2str(alpha)])
disp([['Distance from target = ', num2str(distance)]])
disp([['Number of iterations = ', num2str(N)]])

%Draw the trajectory of the rocket's flight, but interpolate it beforehand with a Bezier curve of third order

%To remove repetitive points
r = unique(r, 'rows');

%Interpolating the missile's trajectory with a 3rd-order Bézier curve
x_interp = linspace(min(r(:,1)), max(r(:,1)), 1000);
y_interp = interp1(r(:,1), r(:,2), x_interp, 'spline');

%Drawing the trajectory of the missile's 2 thickness
plot(x_interp, y_interp, 'LineWidth', 2)
hold on
%Draw target point as a green flag
plot(end_point(1), end_point(2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
%Draw the hit point of the missile as a cross
plot(r(end,1), r(end,2), 'kx', 'MarkerSize', 10, 'LineWidth', 2)
grid on
%Sign the axes with LaTeX and set the font size to 14
xlabel('x, m', 'FontSize', 14)
ylabel('y, m', 'FontSize', 14)
title('Trajectory of the rocket')
legend('Trajectory', 'Target', 'Hit point', 'Location', 'best')