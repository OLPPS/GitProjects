function [sourceIndices,targetIndices] = FindMaskBasedOnCScaseCodes(GateDataFile_plusCS,TargetCScode)

mask = ismember(GateDataFile_plusCS,TargetCScode);
sourceIndices = find(mask);
targetIndices = find(mask) + 1;

end