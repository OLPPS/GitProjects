function GateDataFile = ExtractMCDataTimeSlice(MCDataMatFileName, t1, t2)
% Efficiently extracts the rows from a .mat file where the time
% equals or exceeds t1 but is less than t2. Assumes variable is stored in
% -v7.3 format. The .mat file should have the data needed stored in an
% array named "EventsData", with the 4th column storing the time data.

TimeColNumber = 4;

    % Create matfile object
    mfile = matfile(MCDataMatFileName);

    % Access target column
    times = mfile.('EventsData')(:, TimeColNumber);  % Loads only column number TimeColNumber into RAM

    % Find start index
    startIdx = find(times >= t1, 1, 'first');

    % Find end index
    endRel = find(times(startIdx:end) >= t2, 1, 'first');
    if isempty(endRel)
        endIdx = size(times,1);  % Go to end of array
    else
        endIdx = startIdx + (endRel-1) - 1;
    end

    % Extract only required rows
    GateDataFile = mfile.('EventsData')(startIdx:endIdx, :);
end

