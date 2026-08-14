%STANDALONE analysis script -- NOT part of run_primary_p_pipeline.sh.
%
%Compares "primary" dipole P (calculate_fem_column_power.m: I_dip through
%the cortical column resistance R=rho*L/A, Fig. 1) against "secondary"
%dipole P (a separate DUNEuro/Python volume-conductor simulation of the
%power dissipated in the surrounding tissue by each dipole's own field --
%see meg_power/secondary_p/eeg_dipole_power_computation.py, Fig. 3).
%Secondary P is currently only computed for the original n=4 subject set
%(sub-0002/0004/0006/0007), so this only runs for whichever of those
%subjects have both files on disk -- NOT the full 32-subject list.

SnrFixed = 3; %primary P's SnrFixed value to compare against secondary P
              %(secondary P has no SNR/regularization parameter of its own)

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
POWER_DIR = fullfile(OUTPUTS_DIR, 'primary_p', 'fem_column_power'); % calculate_fem_column_power.m output
SECONDARY_DIR = '/export02/data/vikramn/duneuro_malte/';            % eeg_dipole_power_computation.py output

candidateSubjects = {'sub-0002', 'sub-0004', 'sub-0006', 'sub-0007'};
SubjectNames = {};
for s = 1:numel(candidateSubjects)
    subject = candidateSubjects{s};
    primFile = fullfile(POWER_DIR, sprintf('%s_snr%d_dip_p_i_rms.mat', subject, SnrFixed));
    secFile = fullfile(SECONDARY_DIR, sprintf('%s_dip_fem_p_rms.csv', subject));
    if exist(primFile, 'file') && exist(secFile, 'file')
        SubjectNames{end+1} = subject; %#ok<SAGROW>
    else
        fprintf('Skipping %s: missing %s and/or %s.\n', subject, primFile, secFile);
    end
end

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};
    fprintf('Subject %s\n', subject);

    primData = load(fullfile(POWER_DIR, sprintf('%s_snr%d_dip_p_i_rms.mat', subject, SnrFixed)));
    sub_pri_p = primData.dip_p;
    sub_sec_p = load(fullfile(SECONDARY_DIR, sprintf('%s_dip_fem_p_rms.csv', subject)));

    fprintf('  Correlation between primary and secondary P (120 s RMS): %s\n', num2str(corr(sub_pri_p, sub_sec_p)));
    fprintf('  Total secondary P: %s W\n', num2str(sum(sub_sec_p)));
    fprintf('  Total primary P (SnrFixed=%d): %s W\n', SnrFixed, num2str(sum(sub_pri_p)));
    fprintf('**************************************************\n');
end
