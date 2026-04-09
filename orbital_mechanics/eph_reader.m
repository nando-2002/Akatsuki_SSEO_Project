%% Parse ephemeris data
[utc_time, mission_info] = read_ephemeris_file('horizons_results (4).txt');

% Extract Keplerian elements in correct order: [a, e, i, Ω, ω, θ]
% Note: mission_info(:,10) is semi-major axis, mission_info(:,2) is periapsis distance
% We need to compute eccentricity from periapsis distance and semi-major axis
a_eph = mission_info(:,10);      % Semi-major axis [km]
e_eph = mission_info(:,1);       % Eccentricity (already in column 1)
i_eph = deg2rad(mission_info(:,3));  % Inclination [rad]
OM_eph = deg2rad(mission_info(:,4)); % RAAN [rad]
w_eph = deg2rad(mission_info(:,5));  % Argument of periapsis [rad]
TA_eph = deg2rad(mission_info(:,9)); % True anomaly [rad]

keplerian_elements_eph = [a_eph, e_eph, i_eph, OM_eph, w_eph, TA_eph];

% Plot the six Keplerian elements from ephemeris vs UTC time (in datetime)
elem_names = {'Semi-major axis a [km]','Eccentricity e','Inclination i [deg]',...
    'RAAN \Omega [deg]','Arg. of perigee \omega [deg]','True anomaly \theta [deg]'};
elem_vals = [keplerian_elements_eph(:,1), keplerian_elements_eph(:,2), ...
    rad2deg(keplerian_elements_eph(:,3)), rad2deg(keplerian_elements_eph(:,4)), ...
    rad2deg(keplerian_elements_eph(:,5)), rad2deg(keplerian_elements_eph(:,6))];

% Ensure utc_time is datetime; if it's numeric MJD or similar, try to convert gracefully
if ~isdatetime(utc_time)
    try
        utc_time = datetime(utc_time,'ConvertFrom','datenum');
    catch
        % If conversion fails, create simple time vector in seconds relative to first entry
        utc_time = seconds(utc_time - utc_time(1));
    end
end

for k = 1:6
    figure;
    plot(utc_time, elem_vals(:,k), '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
    grid on;
    xlabel('UTC Time');
    ylabel(elem_names{k});
    title(sprintf('Ephemeris: %s vs UTC', elem_names{k}));
    if k==1
        % give minor formatting for semi-major axis
        ylim([min(elem_vals(:,1))*0.995, max(elem_vals(:,1))*1.005]);
    end
end