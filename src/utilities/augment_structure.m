% Function to generate additional TF fields on the fly and add to
% regression data structure
%augment_structure(table,pd_fields)
% VARIABLES:
%   table: table structure containing base TF fields
%   pd_fields: list of variables used for regression analysis
% RETURNS: 
%   table fo same length as original with additional predictor fileds added
function out_table = augment_structure(table,base_fields, new_pd_fields)
out_table = table;
out_table = removevars(out_table,base_fields);
% for each base predictor, find fields that are functions of that predictor
% and add derivative field to structure
for b = 1:numel(base_fields)
    bp = base_fields{b};
    match_indices = find(contains(new_pd_fields, bp));
    match_indices = match_indices(1:end);
    tf_vec = table.(bp);
    for mi = match_indices
        dp = new_pd_fields{mi};
        space_ind = strfind(dp,'_');
        if ~isempty(space_ind)
            order = str2num(dp(space_ind+1:end));        
        else
            order = 1;
        end
        tf_new = tf_vec.^order;
        if order > 1
            norm_term = nanmean(tf_new(:));
            tf_new = (tf_new - min(tf_new(:))) / norm_term;                
        end
        out_table.(dp) = tf_new; 
    end
end