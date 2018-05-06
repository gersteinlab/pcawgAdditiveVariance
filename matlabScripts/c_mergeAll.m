%% merges all chromosome level SNVstats into one file

fname_all = ['../SNVstats/' cohortName '.obs.null.merged.mat'];

snv_ids_all = cell(0,0);
snv_refs_all = [];
snv_alts_all = [];
samp_ids_all = cell(0,0);
sampXsnv_cell_all = cell(0,0);
N_samp_all = 0;
N_snv_all = 0;

for cChr = 1:22
    
    cChr
    
    fname = strrep(fname_all,'.mat',['.chr' num2str(cChr) '.mat']);
    
    load(fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
        'N_samp','N_snv');
    
    N_snv0 = length(snv_ids_all);
    N_samp0 = length(samp_ids_all);
    snv_ids_all = {snv_ids_all{1:N_snv0} snv_ids{1:N_snv}};
    snv_refs_all = [snv_refs_all snv_refs];
    snv_alts_all = [snv_alts_all snv_alts];
    if isempty(samp_ids_all)
        samp_ids_all = samp_ids;
        sampXsnv_cell_all = sampXsnv_cell;
        continue;
    end
    
    for i = 1:N_samp
        sampXsnv_cell{i} = sampXsnv_cell{i} + N_snv0;
    end
    
    for i = 1:N_samp
        idx = find(strcmp(samp_ids{i},samp_ids_all));
        sampXsnv_cell_all{idx} = [sampXsnv_cell_all{idx} sampXsnv_cell{i}];
    end
    
end

snv_ids = snv_ids_all;
snv_refs = snv_refs_all;
snv_alts = snv_alts_all;
samp_ids = samp_ids_all;
sampXsnv_cell = sampXsnv_cell_all;

% update stats
N_snv = length(snv_ids);
N_samp = length(samp_ids);
display(['# snvs: ' num2str(N_snv)]);
display(['# samps: ' num2str(N_samp)]);
snv_shared = zeros(1,N_snv);
for i = 1:N_samp
    snv_shared(sampXsnv_cell{i}) = snv_shared(sampXsnv_cell{i}) + 1;
end
h_shared = hist(snv_shared,1:max(snv_shared));
% display(h_shared);
% figure(1); plot(h_shared);

save(fname_all,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
    'N_samp','N_snv','snv_shared','h_shared');


