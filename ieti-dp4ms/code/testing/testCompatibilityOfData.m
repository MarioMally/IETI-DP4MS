function bool = testCompatibilityOfData(folderPath,fileName,paperName,tol)
    fileRef = strcat(folderPath,'referenceData/',paperName,fileName);
    fileTest = strcat(folderPath,'computedData/',fileName);

    tableRef = readtable(fileRef,"Delimiter",',',VariableNamingRule='preserve');
    tableTest = readtable(fileTest,"Delimiter",',',VariableNamingRule='preserve');

    bool = true;

    for iVar=tableRef.Properties.VariableNames
        refVals = tableRef(:,iVar).(1);
        testVals = tableTest(:,iVar).(1);
        
        if iscell(refVals)
            if ~all(cellfun(@(x,y) isequal(x,y),refVals,testVals))
                bool = false;
                fprintf('There are not negligiable differences in column  %s!\n',iVar{1});
                return;
            end
        else
            filter = ~isnan(refVals) & ~(refVals > 1e10);
            refVals = refVals(filter);
            testVals = testVals(filter);

            if sum(abs(refVals-testVals)./abs(refVals))/numel(refVals)>tol
                bool = false;
                fprintf('There are not negligiable differences in column  %s!\n',iVar{1});
                return;
            end
        end
    end
end