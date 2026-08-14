%Use auditory oddball tutorial dataset because it has both atlases

s600_17_ind = 14;
s600_7_ind = 15;

s600_17 = surface_file.Atlas(s600_17_ind);
s600_7 = surface_file.Atlas(s600_7_ind);

num_scouts=600;
label_v_dict_7 = cell(num_scouts, 2);
label_v_dict_17 = cell(num_scouts, 2);
%Key idea: all sets of vertices for one match up exactly to the other

%Build region - vertex dictionaries for each atlas (FIRST VERTEX only)
for i=1:num_scouts
    rgn_name_7 = s600_7.Scouts(i).Label;
    first_vert_7 = s600_7.Scouts(i).Vertices(1);
    label_v_dict_7(i,1) = {rgn_name_7};
    label_v_dict_7(i,2) = {first_vert_7};

    rgn_name_17 = s600_17.Scouts(i).Label;
    first_vert_17 = s600_17.Scouts(i).Vertices(1);
    label_v_dict_17(i,1) = {rgn_name_17};
    label_v_dict_17(i,2) = {first_vert_17};

end

%Sort both vertex lists and line them up
col2_7 = cell2mat(label_v_dict_7(:,2));
[~,sort_ind_7] = sort(col2_7);
label_v_dict_7_sort = label_v_dict_7(sort_ind_7,:);

col2_17 = cell2mat(label_v_dict_17(:,2));
[~,sort_ind_17] = sort(col2_17);
label_v_dict_17_sort = label_v_dict_17(sort_ind_17,:);

%Output: a containers.Map() dictionary of 17 net --> 7 net label
s600_17_to_7 = containers.Map();
for i=1:num_scouts
    name_17 = string(label_v_dict_17_sort(i,1));
    name_7 = string(label_v_dict_7_sort(i,1));
    s600_17_to_7(name_17) = name_7;
end

save('s600_17_to_7.mat', 's600_17_to_7');






