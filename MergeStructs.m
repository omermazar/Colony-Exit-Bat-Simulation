 function [merged_struct] = MergeStructs(struct_a,struct_b)
% struct_b is the added structure to struct_a

%%% if one of the structres is empty do not merge
if isempty(struct_a)
    merged_struct=struct_b;
    return
end
if isempty(struct_b)
    merged_struct=struct_a;
    return
end

%%% concatenate
f = fieldnames(struct_b);
 for i = 1:length(f)
    struct_a.(f{i}) = struct_b.(f{i});
 end
 
 merged_struct = struct_a;