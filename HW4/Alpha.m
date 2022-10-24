% Alpha(L) function, which will return the angles at which the missile hits the point with coordinates L.
% To do this, search between which angles the point with coordinates L is located and interpolate between them.
function alpha = Alpha(L)
    %Setting parameters
    g = 9.81; %acceleration of gravity
    k = 0.5; %air drag coefficient
    m = 1; %mass of the missile
    delta = 5; %delta distance
    v0 = 102; % initial velocity
    start_point = [0,0]; %start_point
    end_point = [L,0]; %end_point
    % Set maximum and minimum angles
    alpha_min = 0;
    alpha_max = pi/2;
    %Insert an array of angles
    alpha = [];
    % Initialization of the curve
    curve = [];
    %Shoot N times and interpolate the curve of angle versus distance
    N = 50;
    angle = linspace(alpha_min, alpha_max, N);
    for i=1:N
        r = GetTrajectory(start_point, v0, angle(i), g, k, m);
        curve = [curve; [r(end,1), angle(i)]];
    end
    %Graph of the angle versus distance
    plot(curve(:,1), curve(:,2));
    %Linear interpolation
    for i=1:size(curve,1)-1
        if curve(i,1) <= L && L <= curve(i+1,1) || curve(i,1) >= L && L >= curve(i+1,1)
            alpha = [alpha, curve(i,2) + (curve(i+1,2) - curve(i,2)) * (L - curve(i,1)) / (curve(i+1,1) - curve(i,1))];
        end
    end
end