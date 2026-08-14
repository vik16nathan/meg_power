function v_area_dict=create_tri_area_dict(triangles, vertices_xyz)
    v_area_dict = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    %cellWrapper = struct('data', {}); %allow insert into v_t_dict
    for i = 1:size(triangles,1)
        %Iterate through each vertex of each triangle
        %Find the area using Heron's formula
        tri = triangles(i,:);
        v1_xyz = vertices_xyz(tri(1),:);
        v2_xyz = vertices_xyz(tri(2),:);
        v3_xyz = vertices_xyz(tri(3),:);
        area = heronsFormula(v1_xyz, v2_xyz, v3_xyz);
        for j = 1:3
            ver = triangles(i,j);
            if ~isKey(v_area_dict, ver)
               v_area_dict(ver) = {area};
            else
               v_area_dict(ver) = [v_area_dict(ver), area];
            end
        end
    end
end