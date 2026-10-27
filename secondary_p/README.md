The instructions to derive the three-layer FEM meshes were provided by Eleanor Hill (see attached PDF).
The code to calculate secondary currents at each point in a FEM mesh are largely due to Malte Hoelterschinken's efforts.

Pipeline (run via `run_pipeline.sh`, or each script individually in this order):
1. `convert_fem_mat_to_vc.py` -- any python3. Reads each subject's FEM mesh from
   `../outputs/mergemesh_fem_backup/<subject>/` (build_duneuro_forward_models.m's
   output) and writes `../outputs/secondary_p/vc_3layer/<subject>_vc_3layer.npz`.
2. `eeg_dipole_power_computation.py` -- needs `/usr/bin/python3.8` specifically
   (must match the Python ABI duneuropy.so was built against). For each subject,
   solves the EEG forward problem for all ~15,000 cortical dipoles and writes each
   dipole's total dissipated power to
   `../outputs/secondary_p/fem_p_rms/<subject>_snr3_fem_p_rms.csv`.
3. `eeg_dipole_vsum_sec_curr.py` -- same python3.8 requirement. For a handful of
   exemplar subjects, vector-sums the per-element secondary current density across
   all dipoles and derives the volume/power/current breakdown by tissue
   (brain/skull/scalp), writing to `../outputs/secondary_p/vsum_sec_curr/`.
4. `generate_figures.py` -- any python3. Builds every figure under `figures/` from
   whatever subjects have data so far (degrades gracefully if a subject is only
   partway through the pipeline).

Both duneuropy-dependent stages (2 and 3) reuse one `MEEGDriver3d` per subject
across all of that subject's dipoles instead of rebuilding it per dipole, and
parallelize dipoles across worker processes (`N_WORKERS` in
`eeg_dipole_power_computation.py`) -- each worker builds its own driver once and
reuses it for every dipole assigned to it.

`SECONDARY_P_SUBJECTS=sub-0002,sub-0008` (comma-separated, no spaces) restricts
stages 2-3 to specific subjects, e.g. for a quick single-subject test run.
