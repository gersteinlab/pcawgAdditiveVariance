%% converts null and obs keys into a single key using only common SNVs (with coding/driver/funseq information)
%% inputs and outputs to keys folder 

input_tag = ['../keys/' cohortName];

load([input_tag '.null.key.mat']);
snv_ids_key_rand = snv_ids_key;
snv_cd_key_rand = snv_cd_key;
snv_drv_key_rand = snv_drv_key;
snv_fsq_key_rand = snv_fsq_key;
snv_prm_key_rand = snv_prm_key;
load([input_tag '.obs.key.mat']);
load(strrep([input_tag '.obs.null.merged.mat'],'keys','SNVstats'));

filt = (snv_shared>1);
ordKey_missing = zeros(1,sum(filt));
ordKey_cd = zeros(1,sum(filt));
ordKey_drv = zeros(1,sum(filt));
ordKey_prm = zeros(1,sum(filt));
ordKey_fsq = zeros(1,sum(filt));
ordKey_snv_ids = cell(1,sum(filt));
ordKey_input_snvs = zeros(1,sum(filt));

count = 1;
for i = 1:N_snv
    if mod(i,10000)==0
        display([num2str(i) '/' num2str(N_snv)]);
    end
    if filt(i)==0
        continue;
    end
    
    for j = 1:(N_samp)
        if isempty(strfind(samp_ids{j},'null'))
            if ismember(i,sampXsnv_cell{j})
                ordKey_input_snvs(count) = 1;
                break;
            end
        end
    end
    
    snv_id = snv_ids{i};
    
    if ordKey_input_snvs(count)==1
        idx = find(strcmp(snv_id,snv_ids_key),1);
        if isempty(idx)
            ordKey_missing(count) = 1;
        else
            ordKey_snv_ids{count} = snv_id;
            ordKey_cd(count) = snv_cd_key(idx);
            ordKey_drv(count) = snv_drv_key(idx);
            ordKey_fsq(count) = snv_fsq_key(idx);
		ordKey_prm(count) = snv_prm_key(idx);
        end
    else
        idx = find(strcmp(snv_id,snv_ids_key_rand),1);
        if isempty(idx)
            ordKey_missing(count) = 1;
        else
            ordKey_snv_ids{count} = snv_id;
            ordKey_cd(count) = snv_cd_key_rand(idx);
            ordKey_drv(count) = snv_drv_key_rand(idx);
            ordKey_fsq(count) = snv_fsq_key_rand(idx);
		ordKey_prm(count) = snv_prm_key_rand(idx);
        end        
    end
    
    count = count + 1;
            
end

save([input_tag '.orderedKey.mat'],'ordKey_missing','ordKey_snv_ids',...
    'ordKey_cd','ordKey_drv','ordKey_fsq','ordKey_prm','ordKey_input_snvs');

