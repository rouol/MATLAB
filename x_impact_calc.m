function x_impact = x_impact_calc(alpha, f, t0, T, x0, y0, v0)
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
            fprintf('The time of flight is %f seconds', t_impact);
            [t, x] = ode45(f, [t0, t_impact], [x0, v0 * cosd(alpha), y0, v0 * sind(alpha)]);
            x_impact = x(end, 1);
            break;
        end
    end
end