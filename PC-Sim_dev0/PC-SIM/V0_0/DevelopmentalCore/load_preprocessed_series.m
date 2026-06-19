function TOTAL = load_preprocessed_series(baseName, nFiles)
% baseName: string with a placeholder for the index, e.g. 'FILENAME_%dns_preprocessed'
% nFiles: number of files, e.g. 35

    % --- Load the first file fully ---
    firstFile = sprintf(baseName, 1);
    data = dlmread(firstFile);
    TOTAL(:,1:4) = data(:,1:4);

    % --- Loop over remaining files ---
    for ii = 2:nFiles
        ii
        thisFile = sprintf(baseName, ii);
        data = dlmread(thisFile);
        TOTAL(:, 3+ii) = data(:,4);
    end
end
