%Exports each subject's per-vertex support_areas (same vertex ordering as
%cell_res_avg.m's vertex_res_map_old/vertex_resistance_map, since both are
%computed from the identical per-vertex indexing in that script) to CSV, so
%res_correlations.R can remove the shared support-area (A) term from the
%denominator of both resistance formulations without needing to read a
%MATLAB containers.Map from R (R.matlab cannot parse containers.Map objects).

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
output_dir = '/export02/data/vikramn/hbm_manuscript_code/outputs/resistances/';
thickness_areas = load(fullfile(OUTPUTS_DIR, 'thickness_areas_omega_fem.mat'));

listing = dir(fullfile(output_dir, '*_region_res_cell_column.mat'));
SubjectNames = regexprep({listing.name}, '_region_res_cell_column\.mat$', '');

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};
    outFile = fullfile(output_dir, sprintf('%s_support_areas.csv', subject));
    if exist(outFile, 'file')
        fprintf('Skipping %s: already exported.\n', subject);
        continue;
    end
    support_areas = thickness_areas.subject_support_areas(subject);
    csvwrite(outFile, support_areas);
    fprintf('%s: exported support_areas (%d/%d).\n', subject, s, numel(SubjectNames));
end
