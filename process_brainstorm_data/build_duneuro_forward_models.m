%Builds a 3-layer (scalp/outer-skull/inner-skull) FEM head mesh via Iso2mesh,
%then computes a DUNEuro FEM forward model (head model) for MEG on the
%cortex-surface source space, for every subject in SubjectNames. Mirrors:
%https://neuroimage.usc.edu/brainstorm/Tutorials/Duneuro
%
%Prerequisites:
% 1) preprocess_meg_notch_highpass.m must have already run for a subject
%    (this script looks up that subject's filtered resting condition from
%    preproc_condition_manifest.mat and skips subjects not yet preprocessed).
% 2) Each subject must already have scalp/outer-skull/inner-skull BEM
%    surfaces registered in its Brainstorm anatomy (the tutorial generates
%    these from a segmented MRI before running "Generate FEM mesh" -- check
%    that batch_cat12.m's extramaps=1 output already produces them; if not,
%    a BEM-surface-generation step needs to be added before this script).
%
%Meshing method: 8 subjects (sub-0002/8/19/22/29/32/36/37) already have a
%valid head model computed with 'iso2mesh-2021' + 'mergemesh' -- this is
%the method documented in the manuscript's methods section, so those 8 are
%left alone (skipped via the existingHeadmodel check below). All remaining
%subjects use 'iso2mesh-2026' instead: iso2mesh-2021 failed on them with a
%tissue-labeling error from its mergemesh/mergesurf merge step, while
%iso2mesh-2026 uses an unrelated code path (surfboolean + tetgen) that was
%confirmed by hand in the GUI to succeed on several of them.

bst_headless_init();

% A subject preprocessed for the first time in one MATLAB session isn't
% reliably visible to a later session without a full reload -- confirmed
% for sub-0013 (the one subject preprocessed fresh today, vs. the rest,
% which were already done in earlier sessions): db_reload_subjects(iSubject)
% alone did NOT pick up its newly-created filtered conditions, but this
% protocol-wide reload did. Same class of staleness as the "No head model
% available" bug fixed the same way in process_meg_compute_mn_source_fem.m,
% just one step earlier in the pipeline (preprocessing -> forward modeling
% here, vs. forward modeling -> source estimation there).
db_reload_database('current');

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

sFiles = [];
bst_report('Start', sFiles);

ProtocolInfo = bst_get('ProtocolInfo');

%RestConditionManifest maps subject -> preprocessed (notch+high-pass)
%resting-state condition name, produced by preprocess_meg_notch_highpass.m.
manifest = load(fullfile(OUTPUTS_DIR, 'preproc_condition_manifest.mat'));
RestConditionManifest = manifest.RestConditionManifest;

%Use Brainstorm's current DUNEuro defaults (all layers, default
%conductivities) -- matches the tutorial, which keeps the default
%conductivity values for scalp/skull/brain rather than customizing them.
DuneuroOptions = bst_get('DuneuroOptions');

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};

    if ~isKey(RestConditionManifest, subject)
        fprintf('Skipping %s: not preprocessed yet (run preprocess_meg_notch_highpass.m first).\n', subject);
        continue;
    end
    condition = RestConditionManifest(subject);

    existingHeadmodel = dir(fullfile(ProtocolInfo.STUDIES, subject, condition, '*duneuro*.mat'));
    if ~isempty(existingHeadmodel)
        fprintf('Skipping %s: DUNEuro head model already exists.\n', subject);
        continue;
    end

    % process_fem_mesh (Iso2mesh) needs scalp/outer-skull/inner-skull BEM
    % surfaces already registered on the subject's anatomy -- CAT12's
    % extramaps output does not include these, so generate them here if
    % missing ("Method requires three surfaces... Create them with process
    % 'Generate BEM surfaces' first").
    existingBem = dir(fullfile(ProtocolInfo.SUBJECTS, subject, '*innerskull*.mat'));
    if isempty(existingBem)
        fprintf('=== %s: generating BEM surfaces (scalp/outer-skull/inner-skull) ===\n', subject);
        try
            bst_process('CallProcess', 'process_generate_bem', sFiles, [], ...
                'subjectname', subject, ...
                'nscalp',      1922, ...
                'nouter',      1922, ...
                'ninner',      1922, ...
                'thickness',   4, ...
                'method',      'brainstorm');
        catch err
            fprintf('  FAILED to generate BEM surfaces for %s: %s\n', subject, err.message);
            continue;
        end
        existingBem = dir(fullfile(ProtocolInfo.SUBJECTS, subject, '*innerskull*.mat'));
        if isempty(existingBem)
            fprintf('  FAILED to generate BEM surfaces for %s: no inner-skull surface on disk afterward (check the exported report).\n', subject);
            continue;
        end
    end

    % 'iso2mesh-2021' (with either 'mergemesh' or 'mergesurf') failed for
    % most subjects with "Problem with tissue labels: Brainstorm cannot
    % understand the output labels (1 3)" -- a real iso2mesh/tetgen
    % tissue-labeling failure in that method's merge step. 'iso2mesh-2026'
    % uses a completely different code path (surfboolean + s2m/tetgen1.5,
    % process_fem_mesh.m:628-689) that doesn't go through mergemesh/
    % mergesurf at all, and was confirmed by hand in the GUI to succeed on
    % several subjects that failed under iso2mesh-2021 (sub-0009/10/11/14).
    %
    % IMPORTANT: unlike iso2mesh-2021/iso2mesh, the iso2mesh-2026 code path
    % in process_fem_mesh.m has NO fallback for an empty BemFiles list, and
    % process_fem_mesh's Run() function (the bst_process('CallProcess',...)
    % entry point) never populates OPTIONS.BemFiles at all -- there is no
    % bst_process option for it. Calling it via bst_process('CallProcess',
    % 'process_fem_mesh', ..., 'method','iso2mesh-2026', ...) therefore
    % crashes 100% of the time with "Unrecognized function or variable
    % 'TissueLabels'" (process_fem_mesh.m:647) -- confirmed by inspecting
    % the exported report after a full-cohort run. Brainstorm's own error
    % reporting swallows this (report system, not a MATLAB exception), so
    % the crash doesn't get caught by our try/catch, and the subsequent
    % "find any FEM surface" verification below can silently pick up an
    % unrelated LEFTOVER FEM surface from an earlier iso2mesh-2021 attempt
    % instead of catching the failure -- confirmed: sub-0039/40/41 "passed"
    % this way using their original single iso2mesh-2021 mesh from Jul 29,
    % not a real iso2mesh-2026 mesh at all. Call process_fem_mesh's Compute()
    % directly (bypassing bst_process, matching what the GUI's
    % ComputeInteractive does when the user manually picks BEM surfaces) so
    % BemFiles is always explicitly supplied.
    existingFemMesh = dir(fullfile(ProtocolInfo.SUBJECTS, subject, 'tess_fem_iso2mesh-2026_*V.mat'));
    if ~isempty(existingFemMesh)
        fprintf('=== %s: reusing existing iso2mesh-2026 FEM mesh (%s) ===\n', subject, existingFemMesh(end).name);
    else
        fprintf('=== %s: generating FEM mesh (Iso2mesh-2026, 3 layers) ===\n', subject);
        sSubjectBem = bst_get('Subject', subject);
        if isempty(sSubjectBem.iInnerSkull) || isempty(sSubjectBem.iOuterSkull) || isempty(sSubjectBem.iScalp)
            fprintf('  FAILED %s: missing scalp/outer-skull/inner-skull registration needed to build BemFiles.\n', subject);
            continue;
        end
        femOptions = process_fem_mesh('GetDefaultOptions');
        femOptions.Method        = 'iso2mesh-2026';
        femOptions.BemFiles      = { ...
            sSubjectBem.Surface(sSubjectBem.iInnerSkull).FileName, ...
            sSubjectBem.Surface(sSubjectBem.iOuterSkull).FileName, ...
            sSubjectBem.Surface(sSubjectBem.iScalp).FileName};
        femOptions.MaxVol            = 0.1;
        femOptions.KeepRatio         = 1.0; % fraction (Compute() expects 0-1, not 0-100)
        femOptions.VertexDensity     = 0.5;
        femOptions.NbVertices        = 15000;
        femOptions.isEegCaps         = 0;
        femOptions.NodeShift         = 0.3;
        femOptions.Downsample        = 3;
        femOptions.Zneck             = -115;
        femOptions.ZefMeshResolution = 3;
        femOptions.ZefUseGPU         = 0;
        femOptions.ZefAdvancedInterface = 0;
        [~, iSubjectBem] = bst_get('Subject', subject);
        try
            [isFemOk, femErrMsg] = process_fem_mesh('Compute', iSubjectBem, [], 0, femOptions);
        catch err
            fprintf('  FAILED to generate FEM mesh for %s: %s\n', subject, err.message);
            continue;
        end
        if ~isFemOk
            fprintf('  FAILED to generate FEM mesh for %s: %s\n', subject, femErrMsg);
            continue;
        end
    end

    % Find the FEM surface matching the tess_fem_iso2mesh-2026_*V.mat name
    % pattern specifically (NOT "any FEM surface", which -- as above --
    % could silently be an unrelated leftover from a prior iso2mesh-2021
    % attempt) and set it as this subject's default FEM mesh, matching the
    % GUI's explicit "set as default" step that process_fem_mesh itself
    % does not perform (without it, process_headmodel doesn't know which
    % FEM mesh to use, and DUNEuro exports an empty grid). The in-memory
    % subject record can lag briefly after Compute() writes to disk, so
    % poll/retry a few times rather than give up on the first miss.
    [~, iSubject] = bst_get('Subject', subject);
    iFemSurfaces = [];
    for iAttempt = 1:10
        db_reload_subjects(iSubject);
        sSubject = bst_get('Subject', subject);
        isNewFemMesh = ~cellfun(@isempty, regexp({sSubject.Surface.FileName}, 'tess_fem_iso2mesh-2026_\d+V\.mat$', 'once'));
        iFemSurfaces = find(strcmpi({sSubject.Surface.SurfaceType}, 'FEM') & isNewFemMesh);
        if ~isempty(iFemSurfaces)
            break;
        end
        pause(2);
    end
    if isempty(iFemSurfaces)
        fprintf('  FAILED %s: no iso2mesh-2026 FEM surface found in sSubject.Surface after process_fem_mesh (waited 10x2s).\n', subject);
        continue;
    end
    db_surface_default(iSubject, 'FEM', iFemSurfaces(end));
    [sSubject, iSubject] = bst_get('Subject', subject); % re-fetch: iFEM just changed
    if isempty(sSubject.iFEM)
        fprintf('  FAILED %s: no default FEM surface after db_surface_default (check the exported report).\n', subject);
        continue;
    end

    % panel_duneuro.m's GUI computes FemNames/FemCond/FemSelect fresh for
    % whichever mesh is open (FemSelect for MEG = only the innermost layer,
    % since MEG lead fields only depend on the innermost compartment) --
    % reusing one global bst_get('DuneuroOptions') for every subject leaves
    % these sized/ordered for whatever mesh they last applied to, which
    % silently filters out every element of a mismatched subject's mesh
    % (bst_duneuro.m: iRemove = ~ismember(FemMat.Tissue, find(FemSelect))).
    % Replicated here from panel_duneuro.m's CreatePanel/GetDefaultCondutivity
    % (MEG-only branch, 'simbio' reference conductivities).
    %
    % iso2mesh-2026 meshes have purely numeric TissueLabels ('1','2','3',
    % raw tetgen region IDs) instead of names -- confirmed by direct
    % inspection (and by hand in the GUI) that this numeric ordering has NO
    % consistent meaning across subjects (e.g. label "1" was the
    % largest-volume compartment for one subject and the smallest for
    % another), so neither GetDefaultConductivitySimbio's name matching nor
    % a fixed FemSelect(1)=1 can be trusted blindly here. Reclassify by
    % testing which of the subject's own BEM surfaces each tissue region's
    % boundary vertices coincide with (verified reliable on 4 subjects):
    % scalp touches {scalp, outer-skull} only, skull touches {outer-skull,
    % inner-skull} only, and the innermost compartment touches only
    % inner-skull -- labeled 'brain' (there is no separate CSF surface in
    % this 3-layer scalp/outer-skull/inner-skull setup; the innermost
    % compartment includes CSF elements but is treated as one brain
    % compartment), matching process_fem_mesh.m's own default TissueLabels
    % = {'brain','skull','scalp'} for iso2mesh-2021 when BemFiles isn't
    % explicitly passed (our case), so conductivity assignment stays
    % consistent with the 8 iso2mesh-2021 subjects already processed.
    femFile = sSubject.Surface(sSubject.iFEM(1)).FileName;
    FemNames = ClassifyFemTissueLabels(file_fullpath(femFile), subject, ProtocolInfo);
    DuneuroOptions.FemNames = FemNames;
    DuneuroOptions.FemCond  = GetDefaultConductivitySimbio(FemNames);
    iInnermost = find(contains(lower(FemNames), 'brain'), 1);
    if isempty(iInnermost)
        fprintf('  FAILED %s: could not identify the innermost FEM tissue layer for MEG (labels: %s).\n', ...
            subject, strjoin(FemNames, ', '));
        continue;
    end
    DuneuroOptions.FemSelect = zeros(1, numel(FemNames));
    DuneuroOptions.FemSelect(iInnermost) = 1; % MEG: innermost layer only

    sFilesRest = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname', subject, ...
        'condition',   condition);
    if isempty(sFilesRest)
        fprintf('  SKIPPED head model for %s: no resting-state recordings found in %s.\n', ...
            subject, condition);
        continue;
    end

    fprintf('=== %s: computing DUNEuro FEM head model (MEG, cortex surface) ===\n', subject);
    try
        bst_process('CallProcess', 'process_headmodel', sFilesRest, [], ...
            'sourcespace', 1, ...   %Cortex surface (matches compute_mn_source_fem.m)
            'meg',         5, ...   %DUNEuro FEM
            'duneuro',     DuneuroOptions);
    catch err
        fprintf('  FAILED to compute DUNEuro head model for %s: %s\n', subject, err.message);
        continue;
    end

    % process_fem_mesh/process_headmodel are "import"-type processes that
    % report failures through Brainstorm's report system rather than a
    % MATLAB exception or an empty return value, so verify the real output
    % file exists rather than trusting the calls above completed silently.
    if isempty(dir(fullfile(ProtocolInfo.STUDIES, subject, condition, '*duneuro*.mat')))
        fprintf('  FAILED %s: no DUNEuro head model file on disk afterward (check the exported report for the real cause).\n', subject);
        continue;
    end
end

db_save(); % flush before the next stage's fresh MATLAB session starts

%Save and export report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile, '/export02/data/vikramn/brainstorm3/reports/');


function TissueLabels = ClassifyFemTissueLabels(femFile, subject, ProtocolInfo)
% For iso2mesh-2026 output only: the saved TissueLabels are raw numeric
% tetgen region IDs ('1','2','3'), not names, and their ordering carries no
% consistent anatomical meaning across subjects. Reclassify each label as
% 'brain' (innermost)/'skull'/'scalp' by testing which of the subject's own
% scalp/outer-skull/inner-skull BEM surfaces (generated earlier in this
% script) its boundary vertices spatially coincide with: the scalp
% compartment's vertices lie on {scalp, outer-skull} only, skull on
% {outer-skull, inner-skull} only, and the innermost compartment on
% {inner-skull} only (this compartment includes CSF elements, but there is
% no separate CSF surface in this 3-layer setup, so it's one lumped brain
% compartment, matching process_fem_mesh.m's own default TissueLabels for
% iso2mesh-2021 -- see build_duneuro_forward_models.m for detail). Verified
% against 4 subjects by direct inspection before relying on it here.
% Non-numeric TissueLabels (iso2mesh-2021 etc., already named) are passed
% through unchanged.
femInfo = load(femFile, 'Vertices', 'Elements', 'Tissue', 'TissueLabels');
if any(isnan(str2double(femInfo.TissueLabels)))
    TissueLabels = femInfo.TissueLabels;
    return;
end

scalpFile = dir(fullfile(ProtocolInfo.SUBJECTS, subject, 'tess_head_bem_*V.mat'));
outerFile = dir(fullfile(ProtocolInfo.SUBJECTS, subject, 'tess_outerskull_bem_*V.mat'));
innerFile = dir(fullfile(ProtocolInfo.SUBJECTS, subject, 'tess_innerskull_bem_*V.mat'));
if isempty(scalpFile) || isempty(outerFile) || isempty(innerFile)
    error('%s: missing scalp/outer-skull/inner-skull BEM surface(s) needed to classify iso2mesh-2026 tissue labels.', subject);
end
scalp = load(fullfile(scalpFile(1).folder, scalpFile(1).name), 'Vertices');
outer = load(fullfile(outerFile(1).folder, outerFile(1).name), 'Vertices');
inner = load(fullfile(innerFile(1).folder, innerFile(1).name), 'Vertices');

tol = 1e-5; % 10 micrometers: FEM mesh boundary nodes coincide near-exactly with input BEM surface nodes
uniqueLabels = unique(femInfo.Tissue);
if ~isequal(uniqueLabels(:)', 1:numel(uniqueLabels))
    error('%s: unexpected FEM Tissue values %s (expected 1:N with no gaps).', subject, mat2str(uniqueLabels));
end
TissueLabels = cell(1, numel(uniqueLabels));
for i = 1:numel(uniqueLabels)
    vertIdx = unique(reshape(femInfo.Elements(femInfo.Tissue == uniqueLabels(i), 1:4), [], 1));
    pts = femInfo.Vertices(vertIdx, :);
    [~, dScalp] = knnsearch(scalp.Vertices, pts);
    [~, dOuter] = knnsearch(outer.Vertices, pts);
    [~, dInner] = knnsearch(inner.Vertices, pts);
    touchesScalp = mean(dScalp < tol) > 0.01;
    touchesOuter = mean(dOuter < tol) > 0.01;
    touchesInner = mean(dInner < tol) > 0.01;
    if touchesInner && ~touchesOuter && ~touchesScalp
        TissueLabels{i} = 'brain'; % innermost compartment
    elseif touchesInner && touchesOuter && ~touchesScalp
        TissueLabels{i} = 'skull';
    elseif touchesOuter && touchesScalp && ~touchesInner
        TissueLabels{i} = 'scalp';
    else
        error('%s: could not classify FEM tissue label %d (touches: scalp=%d outer-skull=%d inner-skull=%d).', ...
            subject, uniqueLabels(i), touchesScalp, touchesOuter, touchesInner);
    end
end
end


function FemCond = GetDefaultConductivitySimbio(FemNames)
% Replicates panel_duneuro.m's GetDefaultCondutivity(FemNames, 'simbio'):
% name-matches each tissue label and assigns the corresponding simbio
% reference conductivity (S/m). Default (unmatched) = gray matter.
conductivity = [0.14, 0.33, 1.79, 0.025, 0.008, 0.43]; % white,gray,csf,?,skull,scalp
FemCond = conductivity(2) * ones(1, numel(FemNames));
for i = 1:numel(FemNames)
    name = lower(FemNames{i});
    if contains(name, 'white') || contains(name, 'wm')
        FemCond(i) = conductivity(1);
    elseif contains(name, 'brain') || contains(name, 'grey') || contains(name, 'gray') || ...
            contains(name, 'gm') || contains(name, 'cortex')
        FemCond(i) = conductivity(2);
    elseif contains(name, 'csf') || contains(name, 'inner')
        FemCond(i) = conductivity(3);
    elseif contains(name, 'spong') || contains(name, 'bone') || contains(name, 'skull') || contains(name, 'outer')
        FemCond(i) = conductivity(5);
    elseif contains(name, 'skin') || contains(name, 'scalp') || contains(name, 'head')
        FemCond(i) = conductivity(6);
    end
end
end
