%% merges null and obs SNVstats per chromosome

%% NB can be parallelized across chromosomes

for cChr = 1:22
    
    input_fname = ['../SNVstats/' cohortName '.obs.chr' num2str(cChr) '.mat'];
    null_fname = ['../SNVstats/' cohortName '.null.chr' num2str(cChr) '.mat'];
    
    load(null_fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
        'N_samp','N_snv');
    snv_ids2 = snv_ids;
    snv_refs2 = snv_refs;
    snv_alts2 = snv_alts;
    samp_ids2 = samp_ids;
    sampXsnv_cell2 = sampXsnv_cell;
    N_samp2 = N_samp;
    N_snv2 = N_snv;
    
    load(input_fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
        'N_samp','N_snv');
    
    snv_ids = {snv_ids{1:N_snv} snv_ids2{1:N_snv2}};
    snv_refs = [snv_refs snv_refs2];
    snv_alts = [snv_alts snv_alts2];
    samp_ids = {samp_ids{1:N_samp} samp_ids2{1:N_samp2}};
    sampXsnv_cell = {sampXsnv_cell{1:N_samp} sampXsnv_cell2{1:N_samp2}};
    
    for i = 1:N_samp2
        samp_ids{N_samp+i} = [samp_ids{N_samp+i} '-null'];
        sampXsnv_cell{N_samp+i} = sampXsnv_cell{N_samp+i} + N_snv;
    end
    
    cCmp = {snv_ids{N_snv+1:length(snv_ids)}};
    for i = 1:N_snv
        vec = strcmp(snv_ids{i},cCmp);
        if sum(vec)>0
            display([num2str(i) '/' num2str(N_snv)]);
            idx = find(vec==1);
            idx1 = i; idx2 = idx + N_snv;
            snv_ids = {snv_ids{1:idx2-1} snv_ids{idx2+1:length(snv_ids)}};
            snv_refs = [snv_refs(1:idx2-1) snv_refs(idx2+1:end)];
            snv_alts = [snv_alts(1:idx2-1) snv_alts(idx2+1:end)];
            for j = 1:N_samp2
                vec2 = sampXsnv_cell{N_samp+j};
                vec2(vec2==idx2) = idx1;
                vec2(vec2>idx2) = vec2(vec2>idx2) - 1;
                sampXsnv_cell{N_samp+j} = vec2;
            end
            cCmp = {snv_ids{N_snv+1:length(snv_ids)}};
        end
    end
    
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
    
    fname2 = strrep(input_fname,'obs.chr','obs.null.merged.chr');
    save(fname2,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
        'N_samp','N_snv','snv_shared','h_shared');
    
end
