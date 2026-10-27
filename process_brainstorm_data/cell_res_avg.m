bst_headless_init();

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

PBD = '/export02/data/vikramn/hbm_manuscript_code/process_brainstorm_data/'; % code dir: only for input reference data (cell_props_filtered_redo.csv)
OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
output_dir = '/export02/data/vikramn/hbm_manuscript_code/outputs/resistances/';
thickness_areas = load(fullfile(OUTPUTS_DIR, 'thickness_areas_omega_fem.mat'));

dend_r=1.816*10^7;
gm_res=3.5; %ohm.m, Abboud et al. 2021
%Load the conversion from DKT regions to areal neuron densities
%Load surface file from brainstorm (oddball?)
rgn_neur_dens = readtable(fullfile(PBD, 'cell_props_filtered_redo.csv'));
for s=1:numel(SubjectNames)
    subject = SubjectNames{s};
    surface_file_path = sprintf('/export02/data/vikramn/OmegaSubset/anat/%s/tess_cortex_central_low.mat', subject);
    surface_file = in_bst_data(surface_file_path);
    thickness = thickness_areas.subject_thicknesses(subject);
    support_areas = thickness_areas.subject_support_areas(subject);
    % Uniform-resistivity per-vertex map, precomputed once in
    % parcellate_cortical_column_res_vol.m -- load it instead of
    % recalculating the same gm_res*thickness/support_areas formula here.
    resGm = load(fullfile(output_dir, sprintf('%s_vertex_res_gm.mat', subject)));
    vertex_res_map_old = resGm.vertex_res_map_gm;
    % Each subject's own low-res cortex surface has a slightly different
    % actual vertex count (15002-15011 across the cohort) -- size this to
    % match thickness/support_areas, not a hardcoded constant.
    num_v = numel(thickness);
    vertex_resistance_map=zeros(num_v,1);
    for dkt_region=surface_file.Atlas(2).Scouts
        region=dkt_region.Label;
        vertices_in_rgn = dkt_region.Vertices;

        % Find the row index where "brainstorm_name" is equal to "entorhinal L"
        rowIndex = find(strcmp(rgn_neur_dens.brainstorm_name, region));
        if(isempty(rowIndex))
            % No neuron-density data for this region -- fall back to the
            % uniform-resistivity formula per vertex (previously used a
            % stray, never-looped 'v' left over from whichever region was
            % last matched, silently overwriting the wrong vertex instead
            % of this region's own).
            for i=1:numel(vertices_in_rgn)
                v=vertices_in_rgn(i);
                vertex_resistance_map(v) = gm_res*thickness(v)/support_areas(v);
            end

            disp(region);
            continue
        end

        % Extract the corresponding value in the "neuron_area_density" column
        neuronDensityValue = rgn_neur_dens.neuron_area_density(rowIndex);

        for i=1:numel(vertices_in_rgn)
            v=vertices_in_rgn(i);
            N=neuronDensityValue*(1000)^2*support_areas(v); %convert from neurons/mm^2 --> neurons/m^2
            R=dend_r/N;
            vertex_resistance_map(v)=R;


        end
    end
    save(sprintf('%s%s_region_res_cell_column.mat', output_dir, subject), 'vertex_res_map_old', ...
                'vertex_resistance_map');
end