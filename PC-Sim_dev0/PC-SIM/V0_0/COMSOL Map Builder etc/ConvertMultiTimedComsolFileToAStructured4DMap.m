function [Structured4DMap_AllTimes, xVals, yVals, zVals] = ConvertMultiTimedComsolFileToAStructured4DMap(rawDataArray)
% Constructs a 4D map from raw COMSOL file of all times being mapped.
%
% Inputs:
%   COMSOL file of 2D format with columns (x,y,z,CIE@t_1...CIE@t_end)
%
% Outputs:
%   -Structured4DMap - 4D array of format [xDimension, yDimension, zDimension, timeDimension]
%   -[xVals, yVals, zVals] - ordered list of unique values for each spatial dimension

    % Extract unique grid values (non-uniform spacing allowed)
    xVals = unique(rawDataArray(:,1));
    yVals = unique(rawDataArray(:,2));
    zVals = unique(rawDataArray(:,3));

    % Map each coordinate to its index in the grid
    [~, xIdx] = ismember(rawDataArray(:,1), xVals);
    [~, yIdx] = ismember(rawDataArray(:,2), yVals);
    [~, zIdx] = ismember(rawDataArray(:,3), zVals);

    % Preallocate map
    numberOfTimePoints = size(rawDataArray,2) - 3;
    Structured4DMap_AllTimes = zeros(length(xVals), length(yVals), length(zVals), numberOfTimePoints);

    % Populate map
    for t = 1:numberOfTimePoints
        for n = 1:length(xIdx)
            Structured4DMap_AllTimes(xIdx(n), yIdx(n), zIdx(n), t) = rawDataArray(n, t+3);
        end
    end
end

