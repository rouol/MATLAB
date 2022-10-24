%Task of hitting a target with a missile.
% Plotting the distance vs. angle of approach.

%Setting parameters
g=9.81; %acceleration of gravity
v0=10; %initial velocity
start_point=[0,0]; %start point
k=linspace(0, 10, 50); %resistance coefficient
m=1; %mass of rocket
alpha=linspace(0, pi/2, 20); %angle of launch
L = zeros(length(alpha)); %Define distance vector

% array for marker coordinates
% marker_coordinates = zeros(0, 2);
for n = 1:length(k)
    for i=1:length(alpha)
        %Getting the flight path
        r = GetTrajectory(start_point, v0, alpha(i), g, k(n), m);
        %Getting the flight distance
        L(i) = r(end, 1);
    end

    %Interpolation with cubic spline L(alpha)
    alpha_int = linspace(0, pi/2, 1000);
    L_spline = interp1(alpha, L, alpha_int, 'spline');

    %Build the graph and the maximal distance point on it
    plot(alpha_int, L_spline, 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5]);
    hold on;
    [maxL, maxL_index] = max(L_spline);
    marker_x = alpha_int(maxL_index);
    marker_y = maxL;
    plot(marker_x, marker_y, 'black.', 'MarkerSize', 25-n);
    % save the coordinates of the marker
    % marker_coordinates(n, :) = [marker_x(1), marker_y(1)];
    xlabel('Angle, rad', 'FontSize', 14);
    ylabel('Distance, m', 'FontSize', 14);
    hold on;
end
%Plot the markers
% plot(marker_coordinates(:, 1), marker_coordinates(:, 2), 'red');