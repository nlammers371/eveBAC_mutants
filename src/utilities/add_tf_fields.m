% Function to join in vivo regression table with fixed TF fields

% add_tf_fields(regression_set,tf_data,ap_shift)

function reg_set_out = add_tf_fields(table,tf_data,ap_field_eve2)

reg_set_out = table;
% extract useful fields
ap_vec = tf_data(1).ap_vec;
InterpGrid = tf_data(1).InterpGrid;
ap_shift = table.AP(1)-table.APOrig(1);

for f = 1:numel(tf_data)
    TF = tf_data(f).TF;
    tf_shift = 0;
    % extract tf level vectors
    pt_vec = reshape(tf_data(f).pt_v1_normed,[],1);    
    % adjust Bcd separately, since it also used in vivo AP coordinates
    if strcmp(TF,'Bcd') && strcmp(ap_field_eve2,'AP')
        tf_shift = ap_shift;
    end
    % get mapping indices
    pt_index_array = [repelem(ap_vec'+tf_shift,numel(InterpGrid)),repmat(InterpGrid',numel(ap_vec),1)];    
    [~, aa] = ismember([reg_set_out.(ap_field_eve2) reg_set_out.Time],pt_index_array,'row');
    % add to reg_set_out       
    reg_set_out.(TF) = pt_vec(aa);    
end