function [SelectedSignalStructured3DMap, SelectedSignalTimeIndices3DMap,xVals, yVals, zVals] = GenerateCIEMap(ComsolFileAddress,SelectionCriteria)
%Load and construct 4D array of CIE maps for all times, arranged in ordered and normalised 3D space
[Structured4DMap_AllTimes,xVals, yVals, zVals] = ConvertMultiTimedComsolFileToAStructured4DMap(ComsolFileAddress);
%Make array equal to the same 3D structure but with maximum signal across
%simulated time.
if SelectionCriteria == 'Max'
    [SelectedSignalStructured3DMap, SelectedSignalTimeIndices3DMap] = max(Structured4DMap_AllTimes, [], 4);
else
    SelectedSignalStructured3DMap = 0;
    SelectedSignalTimeIndices3DMap = 0;
end
end