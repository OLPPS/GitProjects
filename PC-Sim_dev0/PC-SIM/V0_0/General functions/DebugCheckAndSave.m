function CurrentWorkspaceVariables = DebugCheckAndSave(DEBUGGING,Label,Suffix,varsBefore)

 if DEBUGGING == 1
            CurrentWorkspaceVariables = who;
            newVars = setdiff(CurrentWorkspaceVariables, varsBefore);
            outputFileName = strcat(Label,'_DDD_PixBrdAllProccessed');
            save(outputFileName, newVars{:})

        elseif DEBUGGING == 2
            varsAfter = who;
            newVars = setdiff(CurrentWorkspaceVariables, varsBefore);
            outputFileName = strcat(Label,'_DDD_PixBrdAllProccessed_v7pt3');
            save(outputFileName, newVars{:}, '-v7.3')

        end

end