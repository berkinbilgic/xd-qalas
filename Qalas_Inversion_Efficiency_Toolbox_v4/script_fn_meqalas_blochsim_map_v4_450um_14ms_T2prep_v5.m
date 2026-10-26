

load_data_type = 1; % 1: grappa, 2: cs, 3: zs-ssl


%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------

addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/qalas_7T/Qalas_Inversion_Efficiency_Toolbox_v4'))

data_path = '/autofs/space/baldur_001/users/yohan/data/qalas/2026_04_20_bay2_invivo_450um_vnav_me_qalas/';
file_name = 'meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V';

% load online recon without mppca denoising
qalas_nifti_filename = [data_path,'nifti/_meqalas_vnav_3e_450um_Vref280V_Inv394V_20260420121712_22_e1.nii'];

if load_data_type ~= 3
    % load grappa or cs recon
    for select_echo = 1:1 % load first echo
        for select_contrast = 1:5
            tic
            fprintf('load data... contrast: %d/%d ',select_contrast,5);

            if load_data_type == 1
                load([data_path,'output/',file_name,'_img_grappa_vol_espirit_rec_contrast_', num2str(select_contrast), '_echo_', num2str(select_echo), '.mat']);
                img(:,:,:,select_contrast) = abs(img_grappa_vol); % magnitude only
            elseif load_data_type == 2
                load([data_path,'output/',file_name,'_Img_CS_vol_espirit_rec_contrast_', num2str(select_contrast), '_echo_', num2str(select_echo), '_lambda_', num2str(1e-07), '.mat']);
                img(:,:,:,select_contrast) = abs(img_cs); % magnitude only
            end
            toc
        end
    end
    img = padarray(img,[0,2,0,0]);
    img = flip(flip(permute(img,[2,1,3,4]),1),2);
end

% load zs-ssl recon
if load_data_type == 3
    tic
    fprintf('load data... ');

    dt_nifti_e1 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e1.nii']);
    img = dt_nifti_e1.vol;
    toc
end

dt_img = load_nifti(qalas_nifti_filename);

N = s(img(:,:,:,1));

if load_data_type == 3
    imagesc3d2( img(:,:,:,1), N/2, 1, [180,0,90], [0, 50])
    imagesc3d2( rsos(img,4), N/2, 2, [180,0,90], [0,100])
else
    imagesc3d2( img(:,:,:,1), N/2, 1, [180,0,90], [0.,3e-5])
    imagesc3d2( rsos(img,4), N/2, 2, [180,0,90], [0.,1e-4])
end

% save rsos qalas image for registration
dt_rsos = dt_img;
dt_rsos.vol = single((flip(flip(permute(rsos(padarray(img,[0,2,0,0]),4),[2,1,3]),1),2)))./max(vec(single(rsos(img,4)))).*max(dt_img.vol(:));
% dt_rsos.vol = single((((rsos(img,4)))))./max(vec(single(rsos(img,4)))).*max(dt_img.vol(:));

% save_nifti(dt_rsos, [data_path,'nifti/rsos_img_zsssl_',file_name,'.nii']);


%--------------------------------------------------------------------------
%% set flags
%--------------------------------------------------------------------------

num_b1_bins = 100;
estimate_pd_map = 1;     % set to 1 to estimate PD map (makes it slower)

use_parfor = 0;
use_inv_eff_lkp = 1;


%--------------------------------------------------------------------------
%% load brain mask
%--------------------------------------------------------------------------

% msk = single(BET(rsos(img,4), N, [1,1,1].*0.45, 0.3));
% genNii(msk, [1,1,1]*.45, [data_path, 'nifti/bet_mask_', file_name, '.nii'])

dt_bet = load_nifti([data_path, 'nifti/bet_mask_zsssl_', file_name, '.nii']);

msk = dt_bet.vol;

if load_data_type ~= 3    
    imagesc3d2( img(:,:,:,1) .* msk, N/2, 11, [180,0,90], [0.,3e-5])
    imagesc3d2( rsos(img,4) .* msk, N/2, 12, [180,0,90], [0.,1e-4])
else
    imagesc3d2( img(:,:,:,1) .* msk, N/2, 11, [180,0,90], [0,50]) % zs-ssl
    imagesc3d2( rsos(img,4) .* msk, N/2, 12, [180,0,90], [0,100]) % zs-ssl
end


%--------------------------------------------------------------------------
%% register afi onto qalas: load itk results
%% in itk-snap: Use mutual information, coarsest level 1x, finest level 1x
%--------------------------------------------------------------------------

% afi acquired in sag direction
dt_afi_std_sag_map_mf0 = load_nifti([data_path, 'nifti/_gre_afi_b1_mgh_wtr_Ref220V_20260420121712_17_reg2zsssl_e1.nii']);      % map

% tic
% dt_afi_std_sag_map_mf9 = dt_afi_std_sag_map_mf0;
% dt_afi_std_sag_map_mf9.vol = medfilt3(dt_afi_std_sag_map_mf0.vol, [9,9,9]);
% dt_afi_std_sag_map_mf9.data_type = 16;
% dt_afi_std_sag_map_mf9.bitpix = 32;
% toc
% save_nifti(dt_afi_std_sag_map_mf9, [data_path, 'nifti/_gre_afi_b1_mgh_wtr_Ref220V_20260420121712_17_reg2zsssl_e1_mf9.nii'])
dt_afi_std_sag_map_mf9 = load_nifti([data_path, 'nifti/_gre_afi_b1_mgh_wtr_Ref220V_20260420121712_17_reg2zsssl_e1_mf9.nii']);

% tic
% dt_afi_std_sag_map_mf9_pfit = dt_afi_std_sag_map_mf0;
% dt_afi_std_sag_map_mf9_pfit.vol = polyfit3D_NthOrder(dt_afi_std_sag_map_mf9.vol./60,msk,5).*60;
% dt_afi_std_sag_map_mf9_pfit.data_type = 16;
% dt_afi_std_sag_map_mf9_pfit.bitpix = 32;
% toc
% save_nifti(dt_afi_std_sag_map_mf9_pfit, [data_path, 'nifti/_gre_afi_b1_mgh_wtr_Ref220V_20251021101600_8a_reg_mf9_pfit3.nii'])

% b1_map_use = flip(flip(permute(dt_afi_std_sag_map_mf9.vol(3:end-2,:,:),[2,1,3]),1),2);
b1_map_use = dt_afi_std_sag_map_mf9.vol;

% imagesc3d2( flip(flip(permute(dt_afi_std_sag_map_mf0.vol(3:end-2,:,:),[2,1,3]),1),2) .* msk / 60, s(img)/2, 7, [180,0,0], [0.,1.5], [], 'afi mf0'), colormap jet
imagesc3d2( dt_afi_std_sag_map_mf0.vol .* msk / 60, s(img)/2, 7, [180,0,90], [0.,1.5], [], 'afi mf0'), colormap jet
imagesc3d2( b1_map_use .* msk / 60, s(img)/2, 8, [180,0,90], [0.,1.5], [], 'afi mf9'), colormap jet
% imagesc3d2( dt_afi_std_sag_map_mf9_pfit.vol .* msk /60, s(img)/2, 9, [180,0,-90], [0.,1.5], [], 'afi mf9 pfit order 3'), colormap jet


%--------------------------------------------------------------------------
%% load json header
%--------------------------------------------------------------------------

qalas_json_filename = '_meqalas_vnav_3e_450um_Vref280V_Inv394V_20260420121712_22_e1.json';

json = fileread([data_path,'nifti/',qalas_json_filename]);

jsn = jsondecode(json);

% esp = res.RepetitionTimeInversion                % echo spacing in sec, VE12U
esp = jsn.RepetitionTime                % echo spacing in sec, XA60
turbo_factor = jsn.EchoTrainLength      % ETL


%--------------------------------------------------------------------------
%% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

b1_val = b1_map_use(msk==1) / 60;

b1_val_nonzero = b1_val(b1_val>0);

figure(100), histogram(b1_val_nonzero, 100), title('afi medfilt 9')

mean(b1_val_nonzero)
[prctile(b1_val_nonzero,.1), prctile(b1_val_nonzero,99.9)]


%--------------------------------------------------------------------------
%% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

b1_val = b1_map_use(msk==1) / 60;

thre_high = prctile(b1_val, 99.9);
thre_low = prctile(b1_val, .1);

% img_b1_load = dt_afi_std_sag_map_mf3.vol / afi_scale_factor;
img_b1_load = b1_map_use / 60;

temp = img_b1_load .* msk;

temp(temp > thre_high) = thre_high;
temp(temp < thre_low) = thre_low;

temp = temp .* msk;
img_b1 = temp .* msk;    

imagesc3d2( img_b1, s(img_b1)/2, 5, [180,0,90], [0.,1.5]), colormap jet

rmse(img_b1, img_b1_load.*msk)


%--------------------------------------------------------------------------
%% create masks for each b1 value
%--------------------------------------------------------------------------

b1_val = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins );

sum_msk = sum(msk(:));

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


idx = 1;

for t = 1:20:length(b1_val)
    imagesc3d2( sum(msk_b1(:,:,:,1:t),4) .* img_b1, N/2, idx, [180,0,0], [0,thre_high]), colormap jet
    
    idx = idx+1;
end
        
imagesc3d2( sum(msk_b1(:,:,:,:),4) .* img_b1, N/2, idx, [180,0,0], [0,thre_high]), colormap jet

msk_b1 = single(msk_b1);


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
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) < t1t2_lut(t,2)
        idx = idx+1;
    end
end


t1t2_lut_prune = zeross( [length(t1t2_lut) - idx, 2] );

if use_inv_eff_lkp ==1
    % FA=500 inversion pulse
    % load /autofs/cluster/berkin/berkin/Matlab_Code_New/qalas_7T/Qalas_Inversion_Efficiency_Toolbox_v4/utils/inversion-efficiency-simulation/data/ielookup_toolbox_v3_FA500.mat
    
    % FA=750 inversion pulse
    load /autofs/cluster/berkin/berkin/Matlab_Code_New/QALAS/qalas_7T/Qalas_Inversion_Efficiency_Toolbox_v4/utils/inversion-efficiency-simulation/data/ielookup_toolbox_v3_FA750.mat

    % FA=1500 inversion pulse
    % load('utils/inversion-efficiency-simulation/data/ielookup_toolbox_v3.mat');

    ie_lkp = reshape((ielookup.ies_mtx),[length(t1t2_lut),1]);
    ie_lkp_prune = zeross( [length(t1t2_lut) - idx, 1] );

    idx = 0;
    for t = 1:length(t1t2_lut)
        if t1t2_lut(t,1) >= t1t2_lut(t,2)
            idx = idx+1;
            t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
            ie_lkp_prune(idx,:) = ie_lkp(t,:);
        end
    end

    length_inv_eff = 1;

else
    idx = 0;
    for t = 1:length(t1t2_lut)
        if t1t2_lut(t,1) >= t1t2_lut(t,2)
            idx = idx+1;
            t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
        end
    end

    length_inv_eff = length(inv_eff);
end


length_dict = length(t1t2_lut_prune);

disp(['dictionary entries: ', num2str(length_dict)])


%--------------------------------------------------------------------------
%%  simulate Mz and Mxy for T1 T2 values in the look up table
%--------------------------------------------------------------------------

addpath '/autofs/cluster/berkin/yohan/matlab_code/qalas/Qalas_Inversion_Efficiency_Toolbox_v4/utils'

TR = 5000e-3;       % vNav TR excluding dead time at the end of the sequence 
alpha_deg = 4;

num_reps = 5;       % number of repetitions to simulate to reach steady state
echo2use = 1;

gap_between_readouts = 900e-3;
time2relax_at_the_end = 0;  % actual TR on the console = TR + time2relax_at_the_end


length_b1_val = length(b1_val);

signal = zeross([length_dict, 5, length_b1_val, length_inv_eff]);


delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; 
parpool( min(ceil(total_cores*.5), length_b1_val) ) 


parfor b1 = 1:length_b1_val                    
    for ie = 1:length_inv_eff
        [Mz, Mxy] = sim_qalas_pd_b1_eff_T2_ie_lkp_YJ(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1),  ie_lkp_prune(:,1));

        temp = abs( Mxy(:,:,end).' );     % transverse magnitude signal

        % normalize every entry of dict. signal
        for n = 1:size(temp,1)
            temp(n,:) = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
        end

        signal(:,:,b1,ie) = temp;
    end    
end

delete(gcp('nocreate'))


%--------------------------------------------------------------------------
% concat dictionary across inv_eff dimension
%--------------------------------------------------------------------------


dict = zeross([length_dict * length_inv_eff, 5, length_b1_val]);

for t = 1:length_inv_eff
    dict(1 + (t-1)*length_dict : t*length_dict, :, :) = signal(:,:,:,t);
end


% output_path = ['/autofs/space/marduk_001/users/yohan/data/qalas/2025_09_04_bay2_qalas_knee_invivo/output/'];
% 
% if ~isfolder(output_path)
%     mkdir(output_path)
% end

tic
    save([data_path, 'output/dict_b1_100bins_medfilt9_fa750.mat'], 'dict', '-v7.3')
toc


%--------------------------------------------------------------------------
%% load dictionary  
%--------------------------------------------------------------------------

% output_path = '/autofs/space/marduk_001/users/yohan/data/qalas/2025_09_04_bay2_qalas_knee_invivo/output/';

load([data_path, 'output/dict_b1_100bins_medfilt9_fa750.mat'])


%--------------------------------------------------------------------------
%% dictionary fit -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------

T1_map = zeross(N);
T2_map = zeross(N);

PD_map = zeross(N);     % proton density-> won't be estimated if estimate_pd_map is set to 0
IE_map = zeross(N);     % inversion efficiency

num_workers = 2;

if use_parfor
    delete(gcp('nocreate'))

    parpool( num_workers )

    parfor slc_select = 1:N(3)
        disp(num2str(slc_select))

        for b = 1:length(b1_val)
            msk_slc = msk_b1(:,:,slc_select,b);
            num_vox = sum(msk_slc(:)~=0);

            if num_vox > 0
                img_slc = zeross([5, num_vox]);

                for t = 1:5
                    temp = sq(img(:,:,slc_select,t));
                    img_slc(t,:) = temp(msk_slc~=0);
                end

                tic
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
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2_ie_lkp_YJ(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use);
                    else
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
                end
                toc

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

    delete(gcp('nocreate'))

else
    for slc_select = 1:N(3)
        tic
        % don't use parfor
        disp(num2str(slc_select))

        for b = 1:length(b1_val)
            msk_slc = msk_b1(:,:,slc_select,b);
            num_vox = sum(msk_slc(:)~=0);

            if num_vox > 0
                img_slc = zeross([5, num_vox]);

                for t = 1:5
                    temp = sq(img(:,:,slc_select,t));
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
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2_ie_lkp_YJ(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use);
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

                toc
            end
        end
    end
end

%%

% save as .mat
if load_data_type == 1
    save([data_path,'output/',file_name,'_T1_map_grappa_b1_100bins_medfilt9_fa750.mat'], 'T1_map')
    save([data_path,'output/',file_name,'_T2_map_grappa_b1_100bins_medfilt9_fa750.mat'], 'T2_map')
    save([data_path,'output/',file_name,'_IE_map_grappa_b1_100bins_medfilt9_fa750.mat'], 'IE_map')
elseif load_data_type == 2
    save([data_path,'output/',file_name,'_T1_map_cs_b1_100bins_medfilt9_fa750.mat'], 'T1_map')
    save([data_path,'output/',file_name,'_T2_map_cs_b1_100bins_medfilt9_fa750.mat'], 'T2_map')
    save([data_path,'output/',file_name,'_IE_map_cs_b1_100bins_medfilt9_fa750.mat'], 'IE_map')
elseif load_data_type == 3
    save([data_path,'output/',file_name,'_T1_map_zsssl_b1_100bins_medfilt9_fa750.mat'], 'T1_map')
    save([data_path,'output/',file_name,'_T2_map_zsssl_b1_100bins_medfilt9_fa750.mat'], 'T2_map')
    save([data_path,'output/',file_name,'_IE_map_zsssl_b1_100bins_medfilt9_fa750.mat'], 'IE_map')
end

if estimate_pd_map
    if load_data_type == 1
        save([data_path,'output/',file_name,'_PD_map_grappa_b1_100bins_medfilt9_fa750.mat'], 'PD_map')
    elseif load_data_type == 2
        save([data_path,'output/',file_name,'_PD_map_cs_b1_100bins_medfilt9_fa750.mat'], 'PD_map')
    elseif load_data_type == 3
        save([data_path,'output/',file_name,'_PD_map_zsssl_b1_100bins_medfilt9_fa750.mat'], 'PD_map')
    end
end


%% T2*

addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/chi-separation_v2/chi-separation-main/Chisep_Toolbox_v1.1.2'))

if load_data_type ~= 3
    % load grappa or cs recon
    for select_echo = 1:3 % load first echo
        for select_contrast = 1:5
            tic
            fprintf('load data... echo: %d/%d, contrast: %d/%d ',select_echo,3,select_contrast,5);

            if load_data_type == 1
                load([data_path,'output/',file_name,'_img_grappa_vol_espirit_rec_contrast_', num2str(select_contrast), '_echo_', num2str(select_echo), '.mat']);
                img_me(:,:,:,select_contrast,select_echo) = abs(img_grappa_vol); % magnitude only
            elseif load_data_type == 2
                load([data_path,'output/',file_name,'_Img_CS_vol_espirit_rec_contrast_', num2str(select_contrast), '_echo_', num2str(select_echo), '_lambda_', num2str(1e-07), '.mat']);
                img_me(:,:,:,select_contrast,select_echo) = abs(img_cs); % magnitude only
            end
            toc
        end
    end
    img_me = padarray(img_me,[0,2,0,0,0]);
    img_me = flip(flip(permute(img_me,[2,1,3,4,5]),1),2);
end

% load zs-ssl recon
if load_data_type == 3
    % load zs-ssl recon
    tic
    fprintf('load data... ');

    % dt_nifti_e1 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e1.nii']);
    % dt_nifti_e2 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e2.nii']);
    % dt_nifti_e3 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e3.nii']);

    dt_nifti_e1 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e1.nii']);
    dt_nifti_e2 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e2.nii']);
    dt_nifti_e3 = load_nifti([data_path, 'nifti/img_zsssl_meas_MID00172_FID23506_meqalas_vnav_3e_450um_Vref280V_Inv394V_e3.nii']);

    img_me(:,:,:,:,1) = dt_nifti_e1.vol;
    img_me(:,:,:,:,2) = dt_nifti_e2.vol;
    img_me(:,:,:,:,3) = dt_nifti_e3.vol;
    toc
end


TE = [0.00284    0.00646    0.01098];

R2_s = r2star_arlo(sq(abs(rsos(img_me,4))), TE*1e3, msk);
T2s_map = 1./R2_s.*1e3.*msk;

if load_data_type == 1
    save([data_path,'output/',file_name,'T2s_map_grappa_rsos_contrast.mat'], 'T2s_map')
elseif load_data_type == 2
    save([data_path,'output/',file_name,'T2s_map_cs_rsos_contrast.mat'], 'T2s_map')
elseif load_data_type ==3
    save([data_path,'output/',file_name,'T2s_map_zsssl_rsos_contrast.mat'], 'T2s_map')
end


R2_c(:,:,:,1) = r2star_arlo(sq(abs(img_me(:,:,:,1,:))), TE*1e3, msk);
R2_c(:,:,:,2) = r2star_arlo(sq(abs(img_me(:,:,:,2,:))), TE*1e3, msk);
R2_c(:,:,:,3) = r2star_arlo(sq(abs(img_me(:,:,:,3,:))), TE*1e3, msk);
R2_c(:,:,:,4) = r2star_arlo(sq(abs(img_me(:,:,:,4,:))), TE*1e3, msk);
R2_c(:,:,:,5) = r2star_arlo(sq(abs(img_me(:,:,:,5,:))), TE*1e3, msk);
T2s_map = 1./mean(R2_c,4).*1e3.*msk;

if load_data_type == 1
    save([data_path,'output/',file_name,'T2s_map_grappa_each_contrast_mean.mat'], 'T2s_map')
elseif load_data_type == 2
    save([data_path,'output/',file_name,'T2s_map_cs_each_contrast_mean.mat'], 'T2s_map')
elseif load_data_type ==3
    save([data_path,'output/',file_name,'T2s_map_zsssl_each_contrast_mean.mat'], 'T2s_map')
end


%--------------------------------------------------------------------------
%% save results as nifti
%--------------------------------------------------------------------------

output_path = '/autofs/space/marduk_001/users/yohan/data/qalas/2025_09_04_bay2_qalas_knee_invivo/output/';

% load([output_path, 'T1_map_zsssl_b1_100bins_medfilt9_fa750.mat'])
% load([output_path, 'T2_map_zsssl_b1_100bins_medfilt9_fa750.mat'])
% load([output_path, 'IE_map_zsssl_b1_100bins_medfilt9_fa750.mat'])

% save([output_path, 'T1_map_zsssl_b1_100bins_medfilt9_fa750_ssl_qalas+.mat'], 'T1_map')
% save([output_path, 'T2_map_zsssl_b1_100bins_medfilt9_fa750_ssl_qalas+.mat'], 'T2_map')
% save([output_path, 'PD_map_zsssl_b1_100bins_medfilt9_fa750_ssl_qalas+.mat'], 'PD_map')
% save([output_path, 'IE_map_zsssl_b1_100bins_medfilt9_fa750_ssl_qalas+.mat'], 'IE_map')

% load the qalas nifti to use its header
qalas_nifti_filename = '/autofs/space/marduk_001/users/yohan/data/qalas/2025_09_04_bay2_qalas_knee_invivo/nifti/img_zsssl_correctHeader_meas_MID00115_FID41751_qls_500um_Vinv537_Vref295_goiaT2RF.nii';

hdr = load_nifti(qalas_nifti_filename);

hdr.datatype = 16;
hdr.bitpix = 32;


hdr.vol = T1_map;
% save_nifti(hdr, [output_path, 'T1_map_zsssl_b1_100bins_medfilt9_polyfit5_fa750_dict.nii'])
save_nifti(hdr, [output_path, 'T1_map_zsssl_b1_100bins_medfilt9_polyfit5_fa750_ssl_qalas+_test.nii'])

hdr.vol = T2_map;
save_nifti(hdr, [output_path, 'T2_map_zsssl_b1_100bins_medfilt9_polyfit5_fa750_ssl_qalas+_test.nii'])

hdr.vol = IE_map;
save_nifti(hdr, [output_path, 'IE_map_zsssl_b1_100bins_medfilt9_polyfit5_fa750_ssl_qalas+_test.nii'])

if estimate_pd_map
    hdr.vol = PD_map;
    save_nifti(hdr, [output_path, 'PD_map_zsssl_b1_100bins_medfilt9_polyfit5_fa750_ssl_qalas+_test.nii'])
end

 