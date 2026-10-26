%%
%
% script_qalas_mapping.m : uses an R=2-fold accelerated 3d-qalas
% acquisition at 1.15 mm iso resolution as well as an lower resolution
% 3d-AFI b1 mapping sequence for T1, T2, PD (and inversion efficiency)
% mapping. data are from a healthy volunteer.
%
% the script expects names of the nifti (and json to read parameters) files
% obtained by converting the dicoms from the scanner using e.g. dcm2niix
%
% the script also assumes that FSL is in the path since it uses BET, FLIRT
% and ROBUSTFOV functions
%
% the example data are in vivo, thus BET is used to accelerate the mapping
% to eliminate non-brain tissue. for phantom experiments, the BET
% functionality will need to be bypassed by (setting use_bet=0).
%
% the script loads a b1 map acquired using the AFI sequence, but it is compatible
% with the Siemens tfl_b1map -> this can be enabled using the "b1_is_afi" flag


% ## [v4] - 2024-08-01

% ### Added
%  - Inversion-efficiency-simulation option
%  - Simulates inversion pulse for a range of T1 and T2 values
%  - Results in smaller dictionary and reduced computation

% ### Updated
%  - Echo spacing import from JSON file
%  - Adapted to accommodate DICOM header changes from VE to XA baseline
%
%
% update: 2025-01-07
% b1_is_afi argument was not being passed to the function in "script_qalas_mapping.m"
% and the function was using the default value. this is now fixed -> reported by
% Maksim Slivka

