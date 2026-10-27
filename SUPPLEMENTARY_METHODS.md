# Supplementary Methods: MEG Forward/Inverse Modeling Pipeline Details

This document records methodological details of the scripted (headless Brainstorm)
MEG processing pipeline in `process_brainstorm_data/`, at the level of detail
needed for a supplementary methods section. It covers decisions that affect the
computed results; purely technical/software fixes that do not change what was
computed (e.g., stale-database-cache workarounds, file-naming corrections) are
omitted.

## 0. Sample

Of 32 subjects with usable MEG resting-state and empty-room recordings, 31
were successfully processed through the full forward/inverse modeling
pipeline described below. One subject (sub-0013) was excluded due to an
unresolved FEM forward-modeling failure (see Section 1). Final analyzed
sample: n=31 (17 male / 14 female; age 21.0–47.0 years, mean 26.9, SD 6.1;
30 right-handed / 1 left-handed; see `outputs/demographics_table.csv`).

## 1. FEM head model construction

Three-layer boundary element (BEM) surfaces (scalp, outer skull, inner skull;
1922 vertices each) were generated from each subject's segmented anatomy using
Brainstorm's template-based method. A tetrahedral finite-element mesh was then
built from these three surfaces using Iso2mesh, and a DUNEuro FEM forward model
was computed for MEG on the cortical surface source space (`tess_cortex_central_low`,
mid-cortical surface — not pial — to keep the source space within gray matter for
FEM analysis).

**Mesh-generation algorithm.** Two Iso2mesh implementations were used across the
cohort:

- **8 of 32 subjects** (sub-0002, 0008, 0019, 0022, 0029, 0032, 0036, 0037): the
  originally-used Iso2mesh algorithm, merging the three BEM surfaces via
  `mergemesh`. This is the method documented in the main Methods.
- **Remaining subjects**: `mergemesh` failed for these subjects' anatomy with a
  tissue-labeling error during volumetric meshing (self-intersecting or
  near-touching BEM surfaces). A different Iso2mesh code path was used instead,
  merging the BEM surfaces via boolean surface union followed by Tetgen
  tetrahedralization. Both approaches produce a 3-compartment (scalp / skull /
  brain) tetrahedral mesh from the same BEM surfaces and are treated as
  methodologically equivalent constructions of the same forward model; the
  choice between them was a per-subject meshing-robustness fallback, not a
  substantive change in head-model definition.

**Cortex-to-inner-skull margin.** For subjects whose cortical surface vertices
fell partly outside the reconstructed inner-skull boundary (causing the FEM
solver to reject affected source locations as outside the mesh volume — subjects
sub-0012 and sub-0034), the BEM generation margin between the cortical surface
and the inner-skull boundary was increased from the default 3 mm to 6 mm before
regenerating the BEM surfaces and FEM mesh for those subjects. This is a
geometric safety margin in surface reconstruction, not a change to tissue
conductivity or forward-model physics.

For one subject (sub-0013), the same issue could not be resolved: increasing
the margin to 6 mm and then 9 mm did not eliminate the out-of-grid source
vertex, indicating a genuine anatomical constraint (cortical surface geometry
too close to the reconstructed inner-skull boundary in this subject even with
a substantially enlarged margin) rather than a numerical/parameter issue.
Sub-0013 was therefore excluded from all downstream analyses. Final analyzed
sample: **31 of 32 subjects**.

## 2. Tissue conductivity assignment

Isotropic conductivities were assigned per SimBio reference values: brain 0.33
S/m, skull 0.008 S/m, scalp 0.43 S/m (no separate white/gray/CSF compartments —
the innermost FEM compartment lumps brain tissue and CSF together, consistent
with the 3-surface BEM construction above).

For subjects processed with the boolean-merge/Tetgen mesh path, tissue
compartments in the output mesh are labeled with arbitrary numeric region IDs
(not anatomical names), and this numeric ordering was confirmed to carry no
consistent meaning across subjects. Tissue identity was therefore determined
per subject by testing which of that subject's own scalp / outer-skull /
inner-skull BEM surfaces each FEM tissue compartment's boundary vertices
spatially coincide with (the innermost compartment touches only the
inner-skull surface; the skull compartment touches inner- and outer-skull; the
scalp compartment touches outer-skull and scalp), then assigning the
corresponding conductivity. This reclassification does not apply to subjects
processed with the `mergemesh` path, whose tissue labels are already
anatomically named by construction.

## 3. DUNEuro forward model computation

The MEG lead field was computed using only the innermost (brain) FEM
compartment (MEG lead fields depend only on conductivity within the innermost
conducting boundary; skull and scalp compartments do not contribute).
Source model: Venant (fitted, CG). Source space and conductivity values were
identical across all subjects regardless of which mesh-generation path was
used (verified post hoc by comparing the computed forward-model parameters —
conductivities, source model, source space file — across all subjects).

## 4. Noise covariance

Noise covariance was computed from each subject's own empty-room recording,
acquired in the same recording session as their resting-state data (never
borrowed across sessions or subjects). The noise recording received the same
preprocessing as the resting-state recording (notch filter at 60/120/180/240/300
Hz; 0.3 Hz high-pass, strict attenuation) before covariance estimation.

## 5. Minimum-norm source estimation

Source estimation used `process_inverse_2018` (minimum-norm, amplitude
measure) with unconstrained (free) source orientation (3 orthogonal dipole
components per cortical vertex, matching the original hand-computed analysis).
Depth weighting was enabled (weight exponent 0.5, weight limit 10) with a loose
orientation constraint of 0.2. Noise covariance was regularized (regularization
factor 0.1). A sensitivity analysis on the MNE regularization parameter was run
by sweeping the fixed-SNR parameter (`SnrFixed`) over 1, 3, and 5.

## 6. One-second epoching for RMS current-dipole calculation

For the primary-current-dipole RMS analysis (see `calculate_parcellate_primary_p/`),
each subject's preprocessed (notch + high-pass filtered) continuous
resting-state recording was segmented into contiguous, non-overlapping 1-second
blocks (cardiac-artifact SSP projectors applied during segmentation; a trailing
partial block shorter than 1 s was discarded). The minimum-norm inverse kernel
was then applied to each of the first 120 one-second blocks to obtain
per-vertex, per-second primary current dipole moments, from which the RMS
across the 120 s window was computed.

## 7. Cortical thickness resampling

CAT12 cortical thickness is computed on each subject's native, high-resolution
central cortical surface. Because the source space and per-vertex support-area
calculations use the lower-resolution cortical surface
(`tess_cortex_central_low`), thickness values were projected from the native
high-resolution surface onto the same low-resolution surface used for source
modeling (same-subject surface-to-surface projection/interpolation) before
being combined with support-area and resistance calculations, so that all
per-vertex quantities share the same vertex indexing.
