# xd-qalas

QALAS T1/T2 mapping using DMv4 for free-breathing liver data. The inputs are scanner DICOMs converted to NIfTI with dcm2niix (not XD-reconstructed yet). B1 correction uses the Siemens `tfl_b1map`.

## script_qalas_mapping_xd_v1.m

Multi-echo QALAS (4 echoes, bipolar readout): T1, T2, inversion efficiency (IE) and R2* mapping.

1. **B0 map** [Hz]: ROMEO on the 4 echoes of the last QALAS contrast, with bipolar phase offset correction.
2. **Bipolar readout distortion**: the readout (voxel axis i, L-R) shifts each echo by B0/bw_px voxels, with alternating sign. Echoes 2 and 4 are resampled by 2·B0/bw_px onto the echo 1/3 geometry (Hanning-windowed sinc, 7 voxels). Echoes 1 and 3 are not resampled.
3. **T1, T2, IE**: root-sum-of-squares over the corrected echoes, then one dictionary fit with B1 binning. The TFL B1 map is resliced onto QALAS using the NIfTI headers (same session, no registration).
4. **R2\*** [1/s]: root-sum-of-squares over the 5 QALAS contrasts of the corrected echoes, then an ARLO fit over the echoes (`r2star_arlo` from Chi-separation Toolbox v1.1.2, in `third_party/chisep/`).

Outputs, in `exp102326_MEQALAS_NIfTIITIs/me_qalas_fb_nii_output/`:

| File | Content |
|---|---|
| `B0_c5_Hz.nii` | B0 map from contrast 5 [Hz] |
| `T1_map_use_rsos_over_echoes.nii` | T1 [ms] |
| `T2_map_use_rsos_over_echoes.nii` | T2 [ms] |
| `IE_map_use_rsos_over_echoes.nii` | inversion efficiency |
| `R2star_map_use_rsos_over_contrasts.nii` | R2* [1/s] |
| `romeo_c5/` | ROMEO inputs and outputs |
| `qalas_rsos.nii.gz` | rsos over the contrasts of echo 1; header and geometry reference for all maps |
| `tfl_b1map2qalas_flirt.nii.gz`, `tfl_img2qalas_flirt.nii.gz` | TFL B1 map and TFL image resliced onto QALAS |

Notes:

- Fat is not corrected: its B0 aliases with the 3.23 ms echo spacing.
- All echoes keep the odd-echo readout displacement of B0/bw_px (about ±0.5 voxel in liver), as in a single-echo acquisition.

## Other scripts

- `script_qalas_mapping_xd_v0.m`: T1, T2, IE mapping from echo 1 only, with TFL B1 correction.
- `script_load_invivo_FB.m`: loads raw twix data (needs mapVBVD, not included).

## Requirements

- Run the scripts from the repository root.
- FSL in the path (`fslmaths`, `fslroi`, `fslmerge`, `flirt`).
- v1 also needs the ROMEO binary (set `romeo_path`) and the Parallel Computing Toolbox.
