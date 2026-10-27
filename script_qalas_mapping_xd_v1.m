%--------------------------------------------------------------------------
%% multi-echo qalas T1, T2, IE and R2* mapping for free-breathing liver data
%% NOTE THAT THE SCRIPT REQUIRES FSL TO BE IN THE PATH
%--------------------------------------------------------------------------

% 1. B0 map with romeo from the last qalas contrast
% 2. bipolar readout: use B0 to move the even echoes onto the odd echo geometry
% 3. T1, T2, IE: rsos over echoes, then dictionary fit with tfl b1 correction
% 4. R2*: rsos over qalas contrasts, then fit over echoes

addpath Qalas_Inversion_Efficiency_Toolbox_v4/utils/
addpath third_party/chisep/     % r2star_arlo (chi-separation toolbox v1.1.2 p-code)

nii_folder = 'exp102326_MEQALAS_NIfTIITIs/me_qalas_fb/';
output_path = 'exp102326_MEQALAS_NIfTIITIs/me_qalas_fb_nii_output/';    % where to save the output maps

% dicoms have been converted to nifti and json using dcm2niix
% qalas magnitude (series 41) and phase (series 42): one nifti per echo, 5 qalas contrasts along dim 4
qalas_mag_prefix = [nii_folder, '_tfl_3d_qalas_liver_pt_4echoes_20261023183146_41_e'];
qalas_phs_prefix = [nii_folder, '_tfl_3d_qalas_liver_pt_4echoes_20261023183146_42_e'];

% siemens tfl_b1map: first series is the raw image, second series is the b1 map
tfl_img_nifti_fname = [nii_folder, '_BH_TB1TFL_2x2x5_middlepos_20261023183146_28.nii'];
tfl_b1_nifti_fname = [nii_folder, '_BH_TB1TFL_2x2x5_middlepos_20261023183146_29.nii'];

romeo_path = '/autofs/cluster/berkin/yuting/MATLAB/demo/mritools_ubuntu-20.04_4.5.3/bin/romeo';


num_echoes = 4;
b0_contrast = 5;                % qalas contrast used for the romeo B0 map: the last one, its B0 map looked cleanest

num_b1_bins = 25;               % number of discrete bins to use to quantize the b1 map (in the range [0.65, 1.35] at 3T)
inv_eff = 0.75:0.025:1;         % inversion efficiency range to use for the Bloch sim dictionary

use_inv_eff_lkp = 1;            % set to 0 to use inv_eff = 0.75:0.025:1
estimate_pd_map = 0;            % set to 0 to avoid estimating proton density -> estimation becomes much faster
use_parfor = 1;                 % use parallel processing with parfor
num_workers = 4;                % number of parfor workers: recommended total_ram / 64 gb -> e.g. use 4 workers for 256 gb ram

DEBUG_MODE = 1;                 % set to 1 to plot images


%--------------------------------------------------------------------------
%% load multi-echo magnitude and json header
%--------------------------------------------------------------------------

dt_img = load_nifti([qalas_mag_prefix, '1.nii']);
N = s(dt_img.vol(:,:,:,1));

% magnitude echoes: [N, 5 qalas contrasts, num_echoes]
img_me = zeross([N, 5, num_echoes]);

TE_s = zeros(1, num_echoes);    % echo times in sec

for e = 1:num_echoes
    dt_echo = load_nifti([qalas_mag_prefix, num2str(e), '.nii']);
    img_me(:,:,:,:,e) = single(reshape(dt_echo.vol, [N, 5]));

    json_echo = jsondecode(fileread([qalas_mag_prefix, num2str(e), '.json']));
    TE_s(e) = json_echo.EchoTime;
end

json_qalas = jsondecode(fileread([qalas_mag_prefix, '1.json']));

esp = json_qalas.RepetitionTime             % echo spacing in sec: XA60 stores the correct echo spacing here
turbo_factor = json_qalas.EchoTrainLength   % ETL
bw_px = json_qalas.PixelBandwidth           % readout bandwidth [Hz/px]


%--------------------------------------------------------------------------
%% rsos qalas image using fsl: geometry reference and header for the output maps
%--------------------------------------------------------------------------

system(['fslmaths ', qalas_mag_prefix, '1.nii -sqr ', output_path, 'qalas_rsos'])
system(['fslmaths ', output_path, 'qalas_rsos', ' -Tmean ', output_path, 'qalas_rsos'])
system(['fslmaths ', output_path, 'qalas_rsos', ' -sqrt ', output_path, 'qalas_rsos'])

hdr = load_nifti([output_path, 'qalas_rsos.nii.gz']);
qalas_rsos = hdr.vol;

hdr.datatype = 16;
hdr.bitpix = 32;


%--------------------------------------------------------------------------
%% B0 map with romeo from the last qalas contrast
%--------------------------------------------------------------------------

% romeo uses the 4 echoes of qalas contrast b0_contrast, with mcpc-3d-s bipolar phase offset correction
romeo_dir = [output_path, 'romeo_c', num2str(b0_contrast)];
[~, ~] = mkdir(romeo_dir);

for e = 1:num_echoes
    system(['fslroi ', qalas_phs_prefix, num2str(e), '_ph.nii ', romeo_dir, '/phase_e', num2str(e), ' ', num2str(b0_contrast-1), ' 1'])
    system(['fslroi ', qalas_mag_prefix, num2str(e), '.nii ', romeo_dir, '/mag_e', num2str(e), ' ', num2str(b0_contrast-1), ' 1'])
end
system(['fslmerge -t ', romeo_dir, '/phase ', sprintf([romeo_dir, '/phase_e%d '], 1:num_echoes)])
system(['fslmerge -t ', romeo_dir, '/mag ', sprintf([romeo_dir, '/mag_e%d '], 1:num_echoes)])
delete([romeo_dir, '/phase_e*.nii.gz'], [romeo_dir, '/mag_e*.nii.gz'])

te_str = ['[', strjoin(arrayfun(@num2str, TE_s*1e3, 'UniformOutput', false), ','), ']'];   % [ms]

% matlab's LD_LIBRARY_PATH makes romeo segfault -> unset for this call; -B filename is relative to the -o folder
system(['env -u LD_LIBRARY_PATH ', romeo_path, ' -p ', romeo_dir, '/phase.nii.gz -m ', romeo_dir, '/mag.nii.gz -t "', te_str, '" -B B0.nii --phase-offset-correction bipolar -o ', romeo_dir])

dt_b0 = load_nifti([romeo_dir, '/B0.nii']);
B0 = single(dt_b0.vol);         % [Hz]

hdr.vol = B0;
save_nifti(hdr, [output_path, 'B0_c', num2str(b0_contrast), '_Hz.nii'])

if DEBUG_MODE
    imagesc3d2( B0, N/2, 30, [180,-0,90], [-150,150]), setGcf(.5), colormap jet, colorbar
end


%--------------------------------------------------------------------------
%% bipolar readout: move e2, e4 onto odd echo geometry
%--------------------------------------------------------------------------

% bipolar readouts shift each echo by B0/bw_px voxels along the readout (voxel axis i, L-R), with alternating sign
% -> sample e2, e4 at i + 2*B0/bw_px so they match the e1, e3 geometry; e1, e3 are not resampled
%    (direction checked on this data: shifted e2 matches sqrt(e1.*e3) better, the opposite shift makes it worse)
% fat is left as is: with 3.23 ms echo spacing its frequency aliases, so its B0 and displacement are wrong
% even echoes: hanning-windowed sinc along i, full width 7 voxels (flirt -interp sinc defaults), weights normalized to sum 1
% edge samples are replicated
shift_px = 2 * B0 / bw_px;
sinc_hw = 3.5;

[I, J, K] = ndgrid(1:N(1), 1:N(2), 1:N(3));
x = I + shift_px;
x0 = floor(x);

img_me_interp = img_me;

for e = 2:2:num_echoes
    for c = 1:5
        vol = img_me(:,:,:,c,e);
        num = zeross(N);
        den = zeross(N);

        for t = -3:4
            n = x0 + t;
            d = x - n;
            w = sin(pi*d) ./ (pi*d) .* (0.5 + 0.5*cos(pi*d/sinc_hw)) .* (abs(d) < sinc_hw);
            w(d == 0) = 1;
            num = num + w .* vol(sub2ind(N, min(max(n, 1), N(1)), J, K));
            den = den + w;
        end

        img_me_interp(:,:,:,c,e) = num ./ den;
    end
end

if DEBUG_MODE
    for e = 1:num_echoes
        imagesc3d2( img_me(:,:,:,1,e), N/2, e, [180,-0,90], [0,400]), setGcf(.5)
        imagesc3d2( img_me_interp(:,:,:,1,e), N/2, 10+e, [180,-0,90], [0,400]), setGcf(.5)
    end
end


%--------------------------------------------------------------------------
%% reslice tfl b1 map onto qalas
%--------------------------------------------------------------------------

% tfl and qalas are from the same session and share scanner coordinates -> reslice using nifti headers, no registration
% -paddingsize 1 stops flirt from filling the outermost qalas planes with its background value
system(['flirt -interp sinc -applyxfm -usesqform -paddingsize 1 -in ', tfl_img_nifti_fname, ' -ref ', output_path, 'qalas_rsos -out ', output_path, 'tfl_img2qalas_flirt'])
system(['flirt -interp sinc -applyxfm -usesqform -paddingsize 1 -in ', tfl_b1_nifti_fname, ' -ref ', output_path, 'qalas_rsos -out ', output_path, 'tfl_b1map2qalas_flirt'])

dt_b1_flirt = load_nifti([output_path, 'tfl_b1map2qalas_flirt.nii.gz']);

% no brain extraction for liver data: fit all voxels
msk = ones(N);


%--------------------------------------------------------------------------
%% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

thre_high = 1.35;
thre_low = 0.65;

% tfl_b1map stores the flip angle x 10 for a nominal 80 deg -> divide by 800 for relative b1
img_b1 = dt_b1_flirt.vol / 800 .* msk;

img_b1(img_b1 > thre_high) = thre_high;
img_b1(img_b1 < thre_low) = thre_low;

img_b1 = img_b1 .* msk;

if DEBUG_MODE
    dt_tfl2qalas = load_nifti([output_path, 'tfl_img2qalas_flirt.nii.gz']);

    imagesc3d2( qalas_rsos, N/2, 31, [180,-0,90], [0,400]), setGcf(.5)
    imagesc3d2( dt_tfl2qalas.vol, N/2, 32, [180,-0,90], [0,400]), setGcf(.5)
    imagesc3d2( img_b1, N/2, 33, [180,-0,90], [0.,1.35]), setGcf(.5), colormap jet
end


%--------------------------------------------------------------------------
%% create masks for each b1 value
%--------------------------------------------------------------------------

b1_val = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins );

if length(b1_val) == 1
    % do not use b1 correction
    msk_b1 = msk;
else
    msk_b1 = zeross([N,length(b1_val)]);

    for t = 1:length(b1_val)
        if t > 1
            msk_b1(:,:,:,t) = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
        else
            % first bin: contains all pixels less than b1_val(1)
            msk_b1(:,:,:,t) = msk.*(img_b1 <= b1_val(t));
        end

        if t == length(b1_val)
            % last bin: contains all pixels greater than b1_val(end-1)
            msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
        end

        msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
    end
end


%--------------------------------------------------------------------------
%% create look up table
%--------------------------------------------------------------------------

t1_entries = [5:5:3000, 3100:100:5000];
t2_entries = [1:1:700, 710:20:1000, 1100:100:2500];

T1_entries = repmat(t1_entries.', [1,length(t2_entries)]).';
T1_entries = T1_entries(:);

T2_entries = repmat(t2_entries.', [1,length(t1_entries)]);
T2_entries = T2_entries(:);

t1t2_lut = cat(2, T1_entries, T2_entries);

% remove cases where T2>T1
keep_entry = t1t2_lut(:,1) >= t1t2_lut(:,2);
t1t2_lut_prune = single(t1t2_lut(keep_entry,:));

if use_inv_eff_lkp
    load('Qalas_Inversion_Efficiency_Toolbox_v4/utils/inversion-efficiency-simulation/data/ielookup_toolbox_v3.mat');

    ie_lkp = reshape((ielookup.ies_mtx),[length(t1t2_lut),1]);
    ie_lkp_prune = single(ie_lkp(keep_entry,:));

    length_inv_eff = 1;
else
    length_inv_eff = length(inv_eff);
end

length_dict = length(t1t2_lut_prune);

disp(['dictionary entries: ', num2str(length_dict)])


TR = 4500e-3;       % excluding dead time at the end of the sequence
alpha_deg = 4;

num_reps = 5;       % number of repetitions to simulate to reach steady state
echo2use = 1;

gap_between_readouts = 900e-3;
time2relax_at_the_end = 0;  % actual TR on the console = TR + time2relax_at_the_end


%--------------------------------------------------------------------------
%% simulate Mz and Mxy for T1 T2 values in the look up table
%--------------------------------------------------------------------------

length_b1_val = length(b1_val);

signal = zeross([length_dict, 5, length_b1_val, length_inv_eff]);

delete(gcp('nocreate'))
pc = parcluster('local');
parpool( min(ceil(pc.NumWorkers*.5), length_b1_val) )

parfor b1 = 1:length_b1_val
    for ie = 1:length_inv_eff
        if use_inv_eff_lkp == 1
            [Mz, Mxy] = sim_qalas_pd_b1_eff_T2_ie_lkp(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1),  ie_lkp_prune(:,1));
        else
            [Mz, Mxy] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1), inv_eff(ie));
        end

        temp = abs( Mxy(:,:,end).' );     % transverse magnitude signal

        % normalize every entry of dict.signal
        for n = 1:size(temp,1)
            temp(n,:) = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
        end

        signal(:,:,b1,ie) = temp;
    end
end

delete(gcp('nocreate'))

% concat dictionary across inv_eff dimension
dict = zeross([length_dict * length_inv_eff, 5, length_b1_val]);

for t = 1:length_inv_eff
    dict(1 + (t-1)*length_dict : t*length_dict, :, :) = signal(:,:,:,t);
end


%--------------------------------------------------------------------------
%% dictionary fit to rsos over echoes -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------

img2use = sum(abs(img_me_interp).^2,5).^.5;

T1_map = zeross(N);
T2_map = zeross(N);

PD_map = zeross(N);     % proton density-> won't be estimated if estimate_pd_map is set to 0
IE_map = zeross(N);     % inversion efficiency

delete(gcp('nocreate'))

if use_parfor
    parpool( num_workers )
    max_workers = num_workers;
else
    max_workers = 0;    % parfor runs serially on the client
end

tic
parfor (slc_select = 1:N(3), max_workers)
    disp(num2str(slc_select))

    for b = 1:length(b1_val)
        msk_slc = msk_b1(:,:,slc_select,b);
        num_vox = sum(msk_slc(:)~=0);

        if num_vox > 0
            img_slc = zeross([5, num_vox]);

            for t = 1:5
                temp = sq(img2use(:,:,slc_select,t));
                img_slc(t,:) = temp(msk_slc~=0);
            end

            res = dict(:,:,b) * img_slc; % dot product

            % find T1, T2 values
            [~, max_idx] = max(abs(res), [], 1);

            max_idx_t1t2 = mod(max_idx, length_dict);
            max_idx_t1t2(max_idx_t1t2==0) = length_dict;

            res_map = t1t2_lut_prune(max_idx_t1t2,:);

            % find inv_eff
            if use_inv_eff_lkp
                ie_to_use = ie_lkp_prune(max_idx_t1t2,1);
            else
                max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;

                ie_to_use = inv_eff(max_idx_ie);
            end

            if estimate_pd_map
                if use_inv_eff_lkp
                    [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2_ie_lkp(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use);
                else
                    [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                end
            end

            t1_map = zeross(N(1:2));
            t1_map(msk_slc==1) = res_map(:,1);

            t2_map = zeross(N(1:2));
            t2_map(msk_slc==1) = res_map(:,2);

            ie_map = zeross(N(1:2));
            ie_map(msk_slc==1) = ie_to_use;

            if estimate_pd_map
                Mxy_sim_use = abs(Mxy_sim(:,:,end));

                scl = zeross([num_vox,1]);

                for idx = 1:size(Mxy_sim_use,2)
                    scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                end

                pd_map = zeross(N(1:2));
                pd_map(msk_slc~=0) = scl;
                PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
            end

            T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
            T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
            IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;
        end
    end
end
toc

delete(gcp('nocreate'))


%--------------------------------------------------------------------------
%% save T1, T2, IE (and PD) maps
%--------------------------------------------------------------------------

hdr.vol = T1_map;
save_nifti(hdr, [output_path, 'T1_map_use_rsos_over_echoes.nii'])

hdr.vol = T2_map;
save_nifti(hdr, [output_path, 'T2_map_use_rsos_over_echoes.nii'])

hdr.vol = IE_map;
save_nifti(hdr, [output_path, 'IE_map_use_rsos_over_echoes.nii'])

if estimate_pd_map
    hdr.vol = PD_map;
    save_nifti(hdr, [output_path, 'PD_map_use_rsos_over_echoes.nii'])
end


%--------------------------------------------------------------------------
%% R2* mapping: rsos over qalas contrasts, then fit over echoes
%--------------------------------------------------------------------------

mag_rsos_contrasts = squeeze(sum(abs(img_me_interp).^2,4).^.5);   % [N, num_echoes]

tic
    r2star = r2star_arlo(mag_rsos_contrasts, TE_s*1e3, msk);       % TE in [ms] -> R2* in [1/s]
toc

hdr.vol = r2star;
save_nifti(hdr, [output_path, 'R2star_map_use_rsos_over_contrasts.nii'])


%--------------------------------------------------------------------------
%% display maps
%--------------------------------------------------------------------------

imagesc3d2( T1_map, N/2, 41, [180,-0,90], [0,2500]), colormap jet, colorbar
imagesc3d2( T2_map, N/2, 42, [180,-0,90], [0,200]), colormap jet, colorbar
imagesc3d2( IE_map, N/2, 43, [180,-0,90], [0,1]), colormap jet, colorbar
imagesc3d2( r2star, N/2, 44, [180,-0,90], [0,150]), colormap jet, colorbar
imagesc3d2( 1e3./(eps+r2star), N/2, 45, [180,-0,90], [0,50]), colormap jet, colorbar      % T2* [ms]
