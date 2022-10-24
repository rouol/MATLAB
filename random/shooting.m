% task 0
% shoot object with a given velocity and angle
% and determine the time of flight and place of impact ground considering air resistance
% and the maximum height reached
% and plot the trajectory

% clear the workspace
close all
clc

% define the initial velocity and angle
v0 = 100;
alpha = 45;
% define the acceleration due to gravity
g = 9.81;
% define the drag coefficient
c = 0.5;
% define the mass of the object
m = 1;
% define the initial position
x0 = 0;
y0 = 0;
% define the time step
dt = 0.01;
% define the initial time
t0 = 0;
% make system of equations
f = @(t, x) [x(2); -c * x(2) / m; x(4); -g - c * x(4) / m];
% define impact position
x_impact = 10;
y_impact = 0;
% t_flight function without air resistance
t_flight = @(v0, theta) 2 * v0 * sin(theta) / g;

% find alpha for which the object hits the ground at x_impact
% using fzero
% alpha = fzero(@(alpha) x_impact - x0 - v0 * cosd(alpha) * t_flight(v0, alpha), alpha);


% print the angle
% fprintf('The angle is %f degrees', alpha);

% plot the trajectory
% using ode45
% [t, x] = ode45(f, [t0, t_flight(v0, alpha)], [x0, v0 * cosd(alpha), y0, v0 * sind(alpha)]);
% plot(x(:, 1), x(:, 3), 'r');
% hold on;

% solve considering air resistance
% using ode45
% for different values of alpha
T = 10;
for alpha = 15:5:75
    % 10 attempts
    for it = 1:1:10
        [t, x] = ode45(f, [t0, T], [x0, v0 * cosd(alpha), y0, v0 * sind(alpha)]);
        max_height = max(x(:, 3));
        max_height_index = find(x(:, 3) == max_height);
        t_impact = interp1(x(max_height_index:end, 3), t(max_height_index:end), 0);
        % t_impact = T
        % check if the t_impact is NaN
        if isnan(t_impact)
            T = T * 10;
        else
            % print the time of flight
            % fprintf('The time of flight is %f seconds', t_impact);
            [t, x] = ode45(f, [t0, t_impact], [x0, v0 * cosd(alpha), y0, v0 * sind(alpha)]);
            % plot the trajectory
            plot(x(:, 1), x(:, 3), 'r');
            hold on;
            break;
        end
    end
end
% [t, x] = ode45(f, [t0, 10], [x0, v0 * cosd(alpha), y0, v0 * sind(alpha)]);
% plot(x(:, 1), x(:, 3), 'b');
% hold on;
xlabel('x')
ylabel('y')
title('Trajectory of the object')
grid on


% find alpha for which the object hits the ground at x_impact
% using newton's method
alpha = 45;
% define the tolerance
tol = 1e-6;
% define the maximum number of iterations
max_iter = 2;

% newton's method
for it = 1:1:max_iter
    % calculate the function value
    x = x_impact_calc(alpha, f, t0, T, x0, y0, v0);
    % calculate the derivative
    dx = (x_impact_calc(alpha + tol, f, t0, T, x0, y0, v0) - x_impact_calc(alpha, f, t0, T, x0, y0, v0)) / tol;
    x
    dx
    alpha
    % update the value of alpha
    alpha = alpha - x / dx * alpha;
    % check if the value of alpha is close enough to the root
    if abs(x-x_impact) < tol
        break;
    end
end