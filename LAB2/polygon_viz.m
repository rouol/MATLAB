% define function that get x and y vertical vectors with coordinates of the points of the polygon
% plot the polygon
% return the perimeter of the polygon
% return the area of the polygon
function [perimeter, area] = polygon_viz(x, y)
    % plot the polygon
    plot(x, y, 'r');
    hold on;
    % get the perimeter of the polygon
    perimeter = 0;
    for i = 1:length(x)
        if i == length(x)
            perimeter = perimeter + sqrt((x(i) - x(1))^2 + (y(i) - y(1))^2);
        else
            perimeter = perimeter + sqrt((x(i) - x(i+1))^2 + (y(i) - y(i+1))^2);
        end
    end
    % get the area of the polygon
    area = 0;
    for i = 1:length(x)
        if i == length(x)
            area = area + (x(i) + x(1)) * (y(i) - y(1));
        else
            area = area + (x(i) + x(i+1)) * (y(i) - y(i+1));
        end
    end
    area = abs(area) / 2;
    % print the perimeter and the area
    fprintf('Perimeter: %f\n', perimeter);
    fprintf('Area: %f\n', area);
end