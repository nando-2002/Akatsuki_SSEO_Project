function [time_utc, keplerian_elements] = read_ephemeris_file(filename)
% READ_EPHEMERIS_FILE Reads Chandra ephemeris file and extracts orbital elements
%   [time_utc, keplerian_elements] = read_ephemeris_file(filename) reads the
%   JPL Horizons ephemeris file and returns:
%   - time_utc: datetime array of observation times
%   - keplerian_elements: matrix containing Keplerian orbital elements
%
%   Keplerian elements columns:
%   1: Eccentricity (EC)
%   2: Periapsis distance (QR) [km]
%   3: Inclination (IN) [deg]
%   4: Longitude of Ascending Node (OM) [deg]
%   5: Argument of Periapsis (W) [deg]
%   6: Time of periapsis (Tp) [Julian Date]
%   7: Mean motion (N) [deg/day]
%   8: Mean anomaly (MA) [deg]
%   9: True anomaly (TA) [deg]
%   10: Semi-major axis (A) [km]
%   11: Apoapsis distance (AD) [km]
%   12: Orbital period (PR) [sec]

    % Open and read the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Read all lines
    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = lines{1};
    fclose(fid);
    
    % Find the start of data section
    data_start = find(contains(lines, '$$SOE'), 1) + 1;
    if isempty(data_start)
        error('Could not find data start marker ($$SOE) in file');
    end
    
    % Find the end of data section
    data_end = find(contains(lines, '$$EOE'), 1) - 1;
    if isempty(data_end)
        % If no end marker, use end of file
        data_end = length(lines);
    end
    
    % Initialize arrays
    num_points = data_end - data_start + 1;
    time_utc = NaT(num_points, 1);
    keplerian_elements = zeros(num_points, 12);
    
    % Parse each data line
    data_count = 0;
    for i = data_start:data_end
        current_line = strtrim(lines{i});
        
        % Skip empty lines or lines that don't contain data
        if isempty(current_line) || startsWith(current_line, '$$')
            continue;
        end
        
        data_count = data_count + 1;
        
        % Split the line by commas
        tokens = strsplit(current_line, ',');
        
        if length(tokens) < 14
            warning('Line %d has insufficient data, skipping', i);
            continue;
        end
        
        % Parse Julian Date (first column)
        jd = str2double(tokens{1});
        
        % Parse calendar date (second column)
        date_str = strtrim(tokens{2});
        
        % Extract just the date part (remove 'A.D.' and 'TDB')
        % The format is: 'A.D. 2025-Nov-14 00:00:00.0000'
        date_parts = strsplit(date_str, ' ');
        
        % Find the parts that contain the actual date and time
        date_time_str = '';
        for k = 1:length(date_parts)
            part = strtrim(date_parts{k});
            if ~isempty(part) && ~strcmp(part, 'A.D.') && ~strcmp(part, 'TDB')
                if isempty(date_time_str)
                    date_time_str = part;
                else
                    date_time_str = [date_time_str ' ' part]; %#ok<AGROW>
                end
            end
        end
        
        % Parse the date and time
        try
            % Try parsing with the expected format
            time_utc(data_count) = datetime(date_time_str, ...
                'InputFormat', 'yyyy-MMM-dd HH:mm:ss.SSSS', ...
                'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
        catch
            try
                % Alternative parsing approach
                % Extract year, month, day, time components
                date_part = extractBetween(date_time_str, 1, 11);
                time_part = extractBetween(date_time_str, 13, length(date_time_str));
                
                % Convert month name to number
                month_names = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
                month_str = extractBetween(date_part, 6, 8);
                month_num = find(strcmp(month_str, month_names), 1);
                
                year_str = extractBetween(date_part, 1, 4);
                day_str = extractBetween(date_part, 10, 11);
                
                % Create datetime from components
                time_utc(data_count) = datetime(str2double(year_str{1}), ...
                    month_num, str2double(day_str{1}), ...
                    'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
                
                % Add time component if available
                if ~isempty(time_part)
                    time_parts = strsplit(time_part{1}, ':');
                    if length(time_parts) >= 3
                        hours = str2double(time_parts{1});
                        minutes = str2double(time_parts{2});
                        seconds_parts = strsplit(time_parts{3}, '.');
                        seconds = str2double(seconds_parts{1});
                        if length(seconds_parts) > 1
                            milliseconds = str2double(seconds_parts{2});
                        else
                            milliseconds = 0;
                        end
                        
                        time_utc(data_count).Hour = hours;
                        time_utc(data_count).Minute = minutes;
                        time_utc(data_count).Second = seconds;
                        time_utc(data_count).Millisecond = milliseconds;
                    end
                end
                
            catch ME
                warning('Failed to parse date: %s. Using Julian date conversion. Error: %s', ...
                    date_time_str, ME.message);
                % Fallback: use Julian date conversion
                time_utc(data_count) = datetime(jd - 1721058.5, 'ConvertFrom', 'modifiedjuliandate');
            end
        end
        
        % Parse Keplerian elements (columns 3-14)
        for j = 1:12
            keplerian_elements(data_count, j) = str2double(tokens{j+2});
        end
    end
    
    % Trim arrays to actual data count
    time_utc = time_utc(1:data_count);
    keplerian_elements = keplerian_elements(1:data_count, :);
    
    fprintf('Successfully read %d ephemeris points\n', data_count);
    fprintf('Time range: %s to %s\n', ...
        char(time_utc(1)), char(time_utc(end)));
end