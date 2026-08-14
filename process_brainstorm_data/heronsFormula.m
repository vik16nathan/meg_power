function area = heronsFormula(v1_xyz, v2_xyz, v3_xyz)
    d_12 = pdist2(v1_xyz, v2_xyz);
    d_23 = pdist2(v2_xyz, v3_xyz);
    d_13 = pdist2(v1_xyz, v3_xyz);
    s = (d_12 + d_23 + d_13)/2;
    area = sqrt(s*(s-d_12)*(s-d_23)*(s-d_13));
end
