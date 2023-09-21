function [value, isterminal, direction] = eventFunction(t, s, LU, MU)
    % Event function to detect both Earth and Sun impacts

    % Earth's radius in meters
    radius_earth = 6.4e3/LU;

    % Sun's radius in meters
    radius_sun = 6.96e5/LU; % 7e5

    % Extract the position components from the state vector y
    x = s(1);
    y = s(2);
    z = s(3);

    % Compute the distance to the center of Earth and the Sun
    distance_to_earth = sqrt((x-(1-MU))^2 + y^2 + z^2) - radius_earth;
    distance_to_sun   = sqrt((x-(-MU))^2 + y^2 + z^2) - radius_sun;

    % The event function should return a vector of values, isterminal, and direction
    value = [distance_to_earth; distance_to_sun];
    isterminal = [1; 1];  % Stop the integration when any of the events are triggered
    direction = [0; 0];   % Detect events regardless of the direction (both crossing upwards and downwards).
end