
function [IndVec] = getManuverInd(ManueverTypeCell,ManToSearch, MaxIndex)

%     IndexCell = strfind(ManueverTypeCell, ManToSearch);
%     FullIndex = find(~(cellfun('isempty', IndexCell)));
    IndexCell = strcmp(ManueverTypeCell, ManToSearch);
    FullIndex = find(IndexCell);
    
    MaxVecIndex = find(FullIndex>=MaxIndex,1)-1;
    
    if ~isempty(MaxVecIndex)
        IndVec = FullIndex(1:MaxVecIndex);
    else % if ~isempty(MaxVecIndex)
        IndVec = FullIndex;
    end  % if ~isempty(MaxVecIndex)


end % function [IndVec] = getManuverInd
