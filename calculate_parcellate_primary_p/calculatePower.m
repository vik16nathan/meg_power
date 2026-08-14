
function outputCellArray=calculatePower(curr_xyz, thickness_in_m, R)
%Power dissipated by an equivalent current dipole moment flowing through a
%cortical column of resistance R (Ohms): P = (I_dip/thickness)^2 * R
%(Fig. 1A). R is supplied by the caller (e.g. the uniform grey-matter
%resistivity column resistance from
%process_brainstorm_data/parcellate_cortical_column_res_vol.m) rather than
%computed here, since different callers derive it differently.
    I_x = curr_xyz(:, :, 1) ./ thickness_in_m;
    I_y = curr_xyz(:, :, 2) ./ thickness_in_m;
    I_z = curr_xyz(:, :, 3) ./ thickness_in_m;

    I_equiv = sqrt(I_x.^2 + I_y.^2 + I_z.^2);
    power_dip_mat_full = I_equiv.^2 .* R;
    power_vs_t = sum(power_dip_mat_full, 1);
    outputCellArray={};
    outputCellArray{1} = power_vs_t;
    outputCellArray{2} = power_dip_mat_full;

end

