%Task of hitting the target with a missile.
%Find the angle at which the missile will hit the target with a given delta distance.
%Realize the problem using the firing method.

%Set parameters
g = 9.81; %acceleration of gravity
k = 0.5; %air drag coefficient
m = 1; %mass of the missile
delta = 1; %delta distance
L = 150; %distance to target
v0 = 102; % initial velocity
start_point = [0,0]; %start_point
end_point = [L,0]; %end_point

%Solution
% The system of equations of motion of the rocket with allowance for air resistance and gravity has the form:
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

angles = [];
is_found = false;
rs = [];

%Firing method implementation
for i = 1:3
    %Set initial values
    alpha_min = pi/6*(i-1);
    alpha_max = pi/6*i;
    distance = NaN;
    r = [];
    N = 0;

    %Till we find the exact angle alpha, we change it to alpha_step
    while true
        %Increase number of iterations
        N=N+1;

        % Find two trajectories of the rocket.
        r1 = GetTrajectory(start_point, v0, (alpha_max+2*alpha_min)/3, g, k, m);
        r2 = GetTrajectory(start_point, v0, (2*alpha_max+alpha_min)/3, g, k, m);

        %Draw the trajectory of the rocket, but before that we interpolate it with a 3rd-order Bézier curve

        %Interpolation of the missile's trajectory with a 3rd-order Bézier curve
        x_interp = linspace(min(r1(:,1)), max(r1(:,1)), 1000)
        y_interp = interp1(r1(:,1), r1(:,2), x_interp, 'spline');

        %Drawing the trajectory of the missile's 2 thickness
        % plot(x_interp, y_interp, 'LineWidth', 2)
        hold on

        %Interpolating the missile's trajectory with a 3rd order Bézier curve
        x_interp = linspace(min(r2(:,1)), max(r2(:,1)), 1000)
        y_interp = interp1(r2(:,1), r2(:,2), x_interp, 'spline');

        %Drawing the trajectory of a missile with a thickness of 2
        % plot(x_interp, y_interp, 'LineWidth', 2)

        %Determine the distance between the end point and the hit point of the missile
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

        %If the difference between the limits of the angle alpha is less than the specified accuracy, then the angle alpha is not found
        if abs(alpha_min - alpha_max) < 0.1
            alpha = (alpha_max + alpha_min)/2;
            distance = abs(r1(end,1) - end_point(1));
            r = r1;
            is_found = false;
            break
        end
    end

    if is_found
        angles = [angles, alpha];
        %Draw the missile's trajectory, but interpolate it with a 3rd-order Bézier curve beforehand

        %Interpolating the missile's trajectory with a 3rd order Bézier curve
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
    end
end

%Result display
disp(angles)
% disp(['Distance from target = ', num2str(distance)])
% disp(['Number of iterations = ', num2str(N)])