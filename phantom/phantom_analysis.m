%%% DIPOLE SCANNING %%%

%Elekta phantom analysis
amplitudes={elekta_dipole_scan.Dipole.Amplitude};
goodness={elekta_dipole_scan.Dipole.Goodness};
amp_vec=[];
num_dipoles=32;
for i=1:num_dipoles
  amp_vec = [amp_vec, norm(amplitudes{i})];
end

disp(mean(amp_vec)/10^-9);
disp(std(amp_vec/10^-9)/sqrt(length(amp_vec)));
disp(mean(cell2mat(goodness)));

%%% MIN NORM %%%

%CTF: note, in microA.m not nanoA.m

%% 200 µA
ctf_mn_200 = in_bst_results(['link|PhantomCTF-ds/phantom_200uA_20150709_01/' ...
   'results_MN_MEG_KERNEL_231018_1231.mat|PhantomCTF-ds/' ...
   'phantom_200uA_20150709_01/data_stim_average_231014_1930.mat'], 1);

t0_cols_200 = [1, 43, 85]; % NOTE: 0 will be ignored (MATLAB is 1-based)
idips_200 = totalIDip_at_cols(ctf_mn_200, t0_cols_200);

fprintf('200 µA: mean = %.6f, SEM = %.6f (n = %d)\n', ...
    mean(idips_200), std(idips_200)/sqrt(numel(idips_200)), numel(idips_200));

%% 20 µA   generalized to match the 200 µA logic
ctf_mn_20 = in_bst_results( ...
  ['link|PhantomCTF-ds/phantom_20uA_20150603_03/' ...
   'results_MN_MEG_KERNEL_231018_1231.mat|PhantomCTF-ds/' ...
   'phantom_20uA_20150603_03/data_stim_average_231014_1930.mat'], 1);

t0_cols_20 = [13, 65, 117, 169, 221, 273, 325]; % amplitude peak columns
idips_20 = totalIDip_at_cols(ctf_mn_20, t0_cols_20);

fprintf('20 µA:  mean = %.6f, SEM = %.6f (n = %d)\n', ...
    mean(idips_20), std(idips_20)/sqrt(numel(idips_20)), numel(idips_20));

%% Elekta
%ELEKTA (load diff protocol... don't run at same time as lines above)

elekta_idips = [];
t60_col = 8; %column corresponding to 60 ms "peak" in Elekta dipole activity
for i=1:num_dipoles
    dip_num = sprintf('%02d', i);
    elekta_dip = in_bst_results(sprintf(['link|Kojak/kojak_all_200nAm_pp_no_chpi_no_ms_raw/' ...
        'results_MN_MEG_GRAD_KERNEL_231020_1214.mat|' ...
        'Kojak/kojak_all_200nAm_pp_no_chpi_no_ms_raw/' ...
        'data_%s_average_231016_1441.mat'], dip_num), 1);
    elekta_amp = elekta_dip.ImageGridAmp(:, t60_col);
    elekta_total_idip = calculateTotalIDip(elekta_amp);
    elekta_idips = [elekta_idips, elekta_total_idip];

end

disp(mean(elekta_idips));
disp(std(elekta_idips)/sqrt(32));

%%%%%%HELPER FUNCTION FOR CTF PHANTOM%%%%%%%%%%%%
%% Helper: compute Total IDip at specific time columns
function total_idips = totalIDip_at_cols(bst_results, t0_cols)
    % bst_results: struct from in_bst_results(...)
    % t0_cols: vector of column indices (1-based)
    nT = size(bst_results.ImageGridAmp, 2);
    t0_cols = t0_cols(t0_cols >= 1 & t0_cols <= nT); % drop invalid cols (e.g., 0)
    total_idips = zeros(1, numel(t0_cols));
    for k = 1:numel(t0_cols)
        amp = bst_results.ImageGridAmp(:, t0_cols(k));
        total_idips(k) = calculateTotalIDip(amp);
    end
end  
