%Vectorized version of previous implementation
function avg_areas = oneThirdSumAreas(num_v, v_area_dict)
    avg_areas = zeros(num_v, 1);
    for v=1:num_v
        attached_areas = v_area_dict(v);
        avg_areas(v) = 1/3*sum(cell2mat(attached_areas));
    end
end