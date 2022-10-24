%Task of hitting the target with a missile.
%Find the angle at which the missile will hit the target with a given delta distance.
%Solve the problem using the firing method.

%Solve the problem
% The system of equations of motion of the missile, taking into account air resistance and gravity, has the form:
%x''(t) + k/m*x'(t) = 0
%y''(t) + k/m*y'(t) + g = 0
% Let's reduce the system to a first-order system:
%x'(t) = v_x(t)
%y'(t) = v_y(t)
%v_x'(t) = -k/m*v_x(t)
%v_y'(t) = -k/m*v_y(t) - g
%Initial conditions:
%x(0) = start_point(1)
%y(0) = start_point(2)
%v_x(0) = v0*cos(alpha)
%v_y(0) = v0*sin(alpha)
% Using the firing method, find the angle alpha at which the missile hits the target at a given delta distance.
%We will use Newton's method to find the angle alpha.

L = 100; %Distance to target
g = 9.81; %acceleration of gravity
k = 0.5; %air drag coefficient
m = 1; %the mass of the missile
delta = 5; %delta distance
v0 = 102; % initial velocity
start_point = [0,0]; %start_point
end_point = [L,0]; %end point

%Fire method implementation
alpha = Alpha(L);

%Draw target point as a green flag
plot(end_point(1), end_point(2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold on

for i = 1:length(alpha)
    r = GetTrajectory(start_point, v0, alpha(i), g, k, m);
    %Interpolation of the missile's trajectory with a 3rd order Bézier curve
    x_interp = linspace(min(r(:,1)), max(r(:,1)), 1000);
    y_interp = interp1(r(:,1), r(:,2), x_interp, 'spline');
    %Draw the trajectory of the missile's 2 thickness
    plot(x_interp, y_interp, 'LineWidth', 2)
    hold on
    %Draw the impact point of the missile as a cross
    plot(r(end,1), r(end,2), 'kx', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
    grid on
    %Sign the axes with LaTeX and set the font size to 14
    xlabel('x, m', 'FontSize', 14)
    ylabel('y, m', 'FontSize', 14)
    title('Trajectory of the rocket')
    legend('Target', 'Trajectory 1', 'Hit point', 'Trajectory 2', 'Location', 'best')
end