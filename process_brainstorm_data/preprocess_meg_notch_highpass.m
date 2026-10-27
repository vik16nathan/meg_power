%Applies the OMEGA tutorial's MEG preprocessing (notch filter, high-pass
%filter, ECG artifact cleaning) to each subject's resting-state recording
%and its matched empty-room noise recording, before FEM forward modeling /
%source estimation. Mirrors the "PRE-PROCESSING" and "ARTIFACT CLEANING"
%sections of brainstorm3/toolbox/script/tutorial_omega.m (the same script
%Brainstorm's own OMEGA tutorial is generated from), applied per-subject
%instead of to the whole database at once, and starting from data that is
%already imported (does not repeat process_import_bids/headpoints/ctf_convert,
%since the OmegaSubset recordings are already imported).
%
%Each subject's raw recordings start at an inconsistent session number
%(e.g. sub-0002 starts at ses-02, sub-0018 starts at ses-01) and use
%different task-rest RUN numbers (e.g. sub-0011's usable rest recording is
%run-02, sub-0013's is run-02 -- neither has a run-01 recording at all).
%So this script picks, per subject, the EARLIEST session that has BOTH a
%task-rest recording (any run number) AND a task-noise recording --
%self-contained, so noise covariance for source estimation can always be
%computed from data in the exact same session, with no cross-session or
%cross-subject matching needed. Confirmed by inspecting every subject's
%BIDS _scans.tsv that all 32 subjects have at least one such session, once
%the rest run number isn't hardcoded to 01.
%
%NOT YET RUN/VERIFIED against a live Brainstorm session.

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

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';

sFiles = [];
bst_report('Start', sFiles);

ProtocolInfo = bst_get('ProtocolInfo');

%RestConditionManifest maps subject -> the resulting filtered (notch +
%high-pass) resting-state condition name, so downstream scripts
%(build_duneuro_forward_models.m, process_meg_compute_mn_source_fem.m) don't
%need to guess Brainstorm's auto-generated folder name.
%NoiseConditionManifest maps subject -> the filtered task-noise condition
%name from that SAME session, so process_meg_compute_mn_source_fem.m can
%compute noise covariance directly from same-session data without any
%cross-session/cross-subject date matching.
RestConditionManifest = containers.Map();
NoiseConditionManifest = containers.Map();

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};

    % Skip subjects whose rest recording has already been preprocessed --
    % find the matching same-session preprocessed noise recording too.
    existingRest = dir(fullfile(ProtocolInfo.STUDIES, subject, '*task-rest_run-*_meg_notch_high'));
    if ~isempty(existingRest)
        restName = existingRest(1).name;
        fprintf('Skipping %s: already preprocessed (%s).\n', subject, restName);
        RestConditionManifest(subject) = restName; %#ok<NASGU>
        tok = regexp(restName, 'ses-(\d+)', 'tokens', 'once');
        ses_num = str2double(tok{1});
        existingNoise = dir(fullfile(ProtocolInfo.STUDIES, subject, sprintf('*ses-%02d_task-noise_*_meg_notch_high', ses_num)));
        if ~isempty(existingNoise)
            NoiseConditionManifest(subject) = existingNoise(1).name; %#ok<NASGU>
        else
            fprintf('  WARNING: %s has no preprocessed same-session (ses-%02d) noise recording -- source estimation will skip this subject.\n', subject, ses_num);
        end
        continue;
    end

    % Find the earliest session with BOTH a task-rest recording (any run
    % number) AND a task-noise recording.
    restListing = dir(fullfile(ProtocolInfo.STUDIES, subject, ...
        sprintf('*%s_ses-*_task-rest_run-*_meg', subject)));
    if isempty(restListing)
        fprintf('Skipping %s: no task-rest recording found.\n', subject);
        continue;
    end
    ses_nums = nan(1, numel(restListing));
    for i = 1:numel(restListing)
        tok = regexp(restListing(i).name, 'ses-(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            ses_nums(i) = str2double(tok{1});
        end
    end
    [sortedSes, sortIdx] = sort(ses_nums);
    restListingSorted = restListing(sortIdx);

    ses_num = NaN;
    restTag = '';
    for i = 1:numel(restListingSorted)
        thisSes = sortedSes(i);
        hasNoise = ~isempty(dir(fullfile(ProtocolInfo.STUDIES, subject, sprintf('*ses-%02d_task-noise_*_meg', thisSes))));
        if hasNoise
            ses_num = thisSes;
            tok = regexp(restListingSorted(i).name, '(task-rest_run-\d+)', 'tokens', 'once');
            restTag = tok{1};
            break;
        end
    end
    if isnan(ses_num)
        fprintf('Skipping %s: no session has both a task-rest recording and a task-noise recording.\n', subject);
        continue;
    end

    fprintf('=== %s: preprocessing ses-%02d (%s + task-noise) ===\n', subject, ses_num, restTag);

    % Select this subject's raw files, then narrow to this session's rest/noise recordings
    sFilesSubj = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname', subject);
    sFilesRestRaw = bst_process('CallProcess', 'process_select_tag', sFilesSubj, [], ...
        'tag',    sprintf('ses-%02d_%s', ses_num, restTag), ...
        'search', 1, ...
        'select', 1);
    sFilesNoiseRaw = bst_process('CallProcess', 'process_select_tag', sFilesSubj, [], ...
        'tag',    sprintf('ses-%02d_task-noise', ses_num), ...
        'search', 1, ...
        'select', 1);
    if isempty(sFilesRestRaw)
        fprintf('  SKIPPED %s: could not re-select the raw rest recording.\n', subject);
        continue;
    end
    if isempty(sFilesNoiseRaw)
        fprintf('  SKIPPED %s: could not re-select the raw noise recording.\n', subject);
        continue;
    end
    sFilesRaw = [sFilesRestRaw, sFilesNoiseRaw];

    try
        % Convert to continuous (matches tutorial_omega.m): these @raw
        % recordings are imported but registered as native epoched/averaged,
        % not continuous, and process_notch/process_bandpass refuse to run
        % on anything but continuous data ("Impossible to process native
        % epoched/averaged files... convert them to continuous").
        sFilesRaw = bst_process('CallProcess', 'process_ctf_convert', sFilesRaw, [], ...
            'rectype', 2);  % Continuous

        % Notch filter: 60/120/180/240/300 Hz (matches tutorial_omega.m)
        % 'MEG' only (not 'MEG, EEG' like tutorial_omega.m): this CTF system
        % has zero EEG channels (confirmed from the channel file), and
        % requesting a sensor type with no matching channels makes these
        % processes return empty rather than just processing MEG.
        sFilesNotch = bst_process('CallProcess', 'process_notch', sFilesRaw, [], ...
            'freqlist',    [60, 120, 180, 240, 300], ...
            'sensortypes', 'MEG', ...
            'read_all',    1);

        % High-pass: 0.3 Hz, strict attenuation (matches tutorial_omega.m)
        sFilesBand = bst_process('CallProcess', 'process_bandpass', sFilesNotch, [], ...
            'sensortypes', 'MEG', ...
            'highpass',    0.3, ...
            'lowpass',     0, ...
            'attenuation', 'strict', ...
            'mirror',      0, ...
            'useold',      0, ...
            'read_all',    1);
    catch err
        fprintf('  FAILED preprocessing for %s: %s\n', subject, err.message);
        continue;
    end
    if isempty(sFilesBand)
        % process_notch/process_bandpass report failures through Brainstorm's
        % own report system rather than a MATLAB exception, so a failure
        % here doesn't get caught above -- it just comes back empty.
        fprintf('  FAILED preprocessing for %s: notch/bandpass returned no files (check the exported report for the real cause).\n', subject);
        continue;
    end

    % sFilesRaw was [rest, noise], and bst_process preserves input order
    sFilesRestBand  = sFilesBand(1);
    sFilesNoiseBand = sFilesBand(2);

    % Artifact cleaning on the rest recording only: detect + regress out cardiac artifact
    try
        bst_process('CallProcess', 'process_evt_detect_ecg', sFilesRestBand, [], ...
            'channelname', 'ECG', ...
            'timewindow',  [], ...
            'eventname',   'cardiac');
        bst_process('CallProcess', 'process_ssp_ecg', sFilesRestBand, [], ...
            'eventname',   'cardiac', ...
            'sensortypes', 'MEG', ...
            'usessp',      1, ...
            'select',      1);
    catch err
        fprintf('  WARNING: ECG artifact cleaning failed for %s (continuing without it): %s\n', ...
            subject, err.message);
    end

    [~, restCondition]  = fileparts(fileparts(sFilesRestBand.FileName));
    [~, noiseCondition] = fileparts(fileparts(sFilesNoiseBand.FileName));
    RestConditionManifest(subject)  = restCondition;
    NoiseConditionManifest(subject) = noiseCondition;
    fprintf('  %s preprocessed -> rest condition "%s", noise condition "%s"\n', ...
        subject, restCondition, noiseCondition);
end

save(fullfile(OUTPUTS_DIR, 'preproc_condition_manifest.mat'), ...
    'RestConditionManifest', 'NoiseConditionManifest');

% Force the protocol database to disk before this session exits -- without
% this, the next stage's fresh MATLAB session can start before the last few
% subjects' additions are flushed, and fail with "No study found" for them
% even though they were just successfully preprocessed.
db_save();

ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile, '/export02/data/vikramn/brainstorm3/reports/');
