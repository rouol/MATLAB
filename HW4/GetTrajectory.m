%The function that returns the trajectory of the missile
function r = GetTrajectory(start_point, v0, angle, g, k, m)
    %If angle is 0, the trajectory will be a point
    if angle == 0
        r = start_point;
        return;
    end
    
    T = 2*v0/g; %time of flight
    
    % Solve a system of equations
    [~, rv] = ode45(@(t,rv) [rv(3); rv(4); -k/m*rv(3); -k/m*rv(4)-g], [0, T], [start_point(1), start_point(2), v0*cos(angle), v0*sin(angle)]);
    
    %Determination of the missile's trajectory
    r = rv(:,1:2);

    %Determining the point of intersection of a trajectory with the OX axis
    for i = 1:length(r)-1
        if r(i,2)*r(i+1,2) < 0
            %Determination of the intersection point coordinate
            x_cross = r(i,1) + (r(i+1,1)-r(i,1))*abs(r(i,2))/(abs(r(i,2))+abs(r(i+1,2)));
            r = [r(1:i,:); [x_cross, 0]];
            break
        end
    end
end
