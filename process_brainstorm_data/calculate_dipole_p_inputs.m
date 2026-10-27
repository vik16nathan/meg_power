%Export each subject's cortical surface vertices and thickness/support-area
%maps -- inputs consumed downstream by parcellate_cortical_column_res_vol.m
%(resistance) and calculate_parcellate_primary_p/rms_idip_rest.m (which
%applies process_meg_compute_mn_source_fem.m's kernel to per-second data
%blocks; this script does not compute or duplicate that step).

bst_headless_init();

% sub-0013 excluded: DUNEuro forward-model computation failed even after
% increasing the BEM brain-inner margin to 3/6/9mm (see SUPPLEMENTARY_METHODS.md
% Section 1) -- final analyzed sample is n=31.
SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';

%Save surface files for all individuals as subject_vertices.csv
for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};
    outFile = fullfile(OUTPUTS_DIR, sprintf('%s_vertices.csv', subject));
    if exist(outFile, 'file')
        fprintf('Skipping %s vertices: already exported.\n', subject);
        continue;
    end
    surface_file = in_bst_data(sprintf(['/export02/data/vikramn/OmegaSubset/anat/' ...
        '%s/tess_cortex_central_low.mat'], subject));
    csvwrite(outFile, surface_file.Vertices);
    fprintf('%s: exported vertices (%d/%d).\n', subject, s, numel(SubjectNames));
end

%Build per-subject support-area / cortical-thickness maps used downstream.
%The CAT12 processing timestamp is discovered automatically from the results
%file already on disk (matched on its 15002V_02 suffix), rather than
%hardcoded per subject, since it's just a processing-run timestamp, not
%something meaningful to record by hand for 32 subjects.
thicknessAreasFile = fullfile(OUTPUTS_DIR, 'thickness_areas_omega_fem.mat');
if exist(thicknessAreasFile, 'file')
    % Resume from a prior (possibly stall-killed) attempt instead of
    % redoing create_tri_area_dict for every subject -- this loop produces
    % no incremental console output, so a restart with no checkpoint would
    % always redo the same slow work and get killed at the same point again.
    prev = load(thicknessAreasFile);
    subject_v_area_dicts = prev.subject_v_area_dicts;
    subject_support_areas = prev.subject_support_areas;
    subject_thicknesses = prev.subject_thicknesses;
else
    subject_v_area_dicts = containers.Map();
    subject_support_areas = containers.Map();
    subject_thicknesses = containers.Map();
end
for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};
    if isKey(subject_thicknesses, subject)
        fprintf('Skipping %s thickness/area maps: already computed.\n', subject);
        continue;
    end
    %surface_file = in_bst_data(sprintf([ '/export02/data/vikramn/OmegaSubset/' ...
    %    'anat/%s/tess_cortex_pial_low.mat'], subject)); %pial
    surface_file = in_bst_data(sprintf([ '/export02/data/vikramn/OmegaSubset/' ...
        'anat/%s/tess_cortex_central_low.mat'], subject)); %NOT pial
    % Each subject's own low-res cortex surface has a slightly different
    % actual vertex count (15002-15011 observed across the cohort -- a
    % normal side effect of Brainstorm's surface downsampling, not an
    % error), so support_areas must be sized to THIS subject's real vertex
    % count rather than a fixed constant. Using a hardcoded 15002 here
    % silently truncated support_areas for any subject whose real surface
    % had more vertices, producing a length mismatch against the
    % correctly-per-subject-sized projected thickness vector below.
    subject_num_v = size(surface_file.Vertices, 1);
    v_area_dict = create_tri_area_dict(surface_file.Faces, surface_file.Vertices);
    subject_v_area_dicts(subject) = v_area_dict;
    subject_support_areas(subject) = oneThirdSumAreas(subject_num_v, v_area_dict);

    cat12_dir = sprintf('/export02/data/vikramn/OmegaSubset/data/%s/CAT12/', subject);
    thickness_listing = dir(sprintf('%sresults_surface_thickness_*.mat', cat12_dir));
    % bst_project_sources (below) saves its projected output back into this
    % SAME folder, named "<original>_<N>V.mat" -- exclude it here so a
    % re-run doesn't see its own prior output as a second candidate "original".
    isProjectedOutput = ~cellfun(@isempty, regexp({thickness_listing.name}, '_\d+V\.mat$', 'once'));
    thickness_listing = thickness_listing(~isProjectedOutput);
    if numel(thickness_listing) ~= 1
        error('Expected exactly one results_surface_thickness_*.mat file for %s, found %d', ...
            subject, numel(thickness_listing));
    end

    % CAT12 imports thickness onto the high-res native cortex
    % (tess_cortex_central_high.mat: 215667 vertices for sub-0002), not the
    % low-res surface everything else here uses (support_areas, the MN
    % source space) -- confirmed by loading both directly: thickness came
    % back (215667,1) vs. support_areas (15002,1) for the same subject.
    % Project it onto the subject's own low-res cortex first so it shares
    % the same per-vertex indexing as support_areas.
    sSubject = bst_get('Subject', subject);
    iCortexLow = find(strcmpi({sSubject.Surface.SurfaceType}, 'Cortex') & ...
        contains({sSubject.Surface.FileName}, 'tess_cortex_central_low'));
    if isempty(iCortexLow)
        error('%s: could not find tess_cortex_central_low.mat among registered Cortex surfaces.', subject);
    end
    lowResCortexFile = sSubject.Surface(iCortexLow(1)).FileName;

    thicknessFileRel = fullfile(subject, 'CAT12', thickness_listing(1).name);
    projectedFiles = bst_project_sources({thicknessFileRel}, lowResCortexFile, 0, 0);
    if isempty(projectedFiles)
        error('%s: bst_project_sources failed to project thickness onto the low-res cortex (check the exported report).', subject);
    end
    thickness = in_bst_results(projectedFiles{1}, 1);
    thickness_in_m = thickness.ImageGridAmp(:,1) ./ 1000;
    if numel(thickness_in_m) ~= subject_num_v
        error('%s: projected thickness has %d vertices but support_areas has %d -- projection landed on the wrong surface.', ...
            subject, numel(thickness_in_m), subject_num_v);
    end
    subject_thicknesses(subject) = thickness_in_m;

    save(thicknessAreasFile, ...
        'subject_support_areas', 'subject_v_area_dicts', 'subject_thicknesses');
    fprintf('%s: computed thickness/area maps (%d/%d).\n', subject, s, numel(SubjectNames));
end
