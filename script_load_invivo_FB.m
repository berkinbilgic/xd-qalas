%--------------------------------------------------------------------------
%%  look at Cemre's email, use specific mapvbvd
%--------------------------------------------------------------------------

restoredefaultpath

addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/QALAS/abdominal/mapVBVD-PMU/

data_path = '/autofs/space/marduk_004/users/berkin/2026_03_11_bay5_qalas_freebreathing_invivo/';

file_name = 'meas_MID00577_FID09192_qalas_350FOV_R1_2mm';

tic
    dt = mapVBVD([data_path, file_name, '.dat'], 'removeOS');

    kspace = dt{end}.image();
toc

%--------------------------------------------------------------------------
%%   
%--------------------------------------------------------------------------

addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'))

kspace = squeeze(kspace);

kspace = permute(kspace,[1,3,4,2,5]);

img = ifft3call(kspace);

img_rsos = rsos(img,4);

for t = 1:5
    imagesc3d2(img_rsos(:,:,:,:,t),size(img_rsos)/2,t,[90,90,90],[0,2e-4])
end

%--------------------------------------------------------------------------
%%    bin PT signal
%--------------------------------------------------------------------------

addpath('/autofs/cluster/berkin/yuting/data/20250802_abd_invivo/utils')

tabPT = transpose(abs(dt{end}.PMUdata.RESP));


figure(11), subplot(211), plot(tabPT(1:2000,1))
figure(11), subplot(212), plot(tabPT(1:2000,2))


%%

num_state = 4;

PT = tabPT(:,2);

% According to the length of PT, it's acquired per readout
pt_bin = binningRespSignal_cyt(PT,num_state,1); 


state_ind = zeros(size(pt_bin));
for ii = 1:num_state
    for jj = 1:size(pt_bin,1)
        state_ind(jj,ii) = find(PT==sq(pt_bin(jj,ii)));
    end
end

ind_sort = sort(state_ind,1);

tt = unique(ind_sort(:));
assert(length(tt)==length(PT))

% show example of state #1
TR_index   = floor((sq(ind_sort(:,1))-1)/(nETL*5)) + 1;   
echo_index = mod(sq(ind_sort(:,1))-1, (nETL*5)) + 1;     

[unique_vals, ~, idx] = unique(TR_index);

counts = accumarray(idx, 1);

figure(12); scatter(TR_index(:),echo_index(:),5,'filled')
xlabel('index of TR')
ylabel('t')

%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

restoredefaultpath
addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'))


tic
    dt2 = mapVBVD2([data_path, file_name, '.dat'], 'removeOS');

    kspace2 = dt2{end}.image();
toc

%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

addpath /cluster/berkin/berkin/Matlab_Code_New/QALAS/abdominal/read_log/


tic
    Resp = Read_respdata_VIDA([data_path, 'qalas_R1_2mm.resp']);
toc

%%

figure(20), subplot(211), plot(Resp.Rsr_preprcd_HeadFirst)
figure(20), subplot(212), plot(Resp.Rsr_raw_HeadFirst)
