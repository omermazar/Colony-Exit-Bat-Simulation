function [AllParams] = ReadDefaultParamsTable(filename)

% reading default parameters table
% filename = 'DATA\DefaultParamsTable.xlsx';

% general read
opts = detectImportOptions(filename);
shts = sheetnames(filename);
opts.VariableTypes{strcmp(opts.VariableNames, 'Value')} = 'char';

AllParams= struct();

for iSheet = shts'
%     display(iSheet);
    t2 = readtable(filename, opts, 'Sheet', iSheet, 'ReadRowNames',true);
    % convert rowNames to fieldnames
    fieldNames = t2.Properties.RowNames;
    for k = 1:numel(fieldNames)
        iField = fieldNames{k};
        numValue = str2double(t2.Value{k});
            
        if ~isnan(numValue)
            % numeric values
            AllParams.(iSheet).(iField) = numValue;
            %
        elseif ~contains(t2.Value{k}, '[') 
            % Char values
            AllParams.(iSheet).(iField) = t2.Value{k};
        else 
            % array
            str1 = erase(t2.Value{k}, '[');
            str1 = erase(str1, ']');
            Str = sprintf('%s,', str1);
            AllParams.(iSheet).(iField) = sscanf(Str, '%g,', [3, inf]).'; 

        end % if

    end % for k
end % for iSheet = shts'