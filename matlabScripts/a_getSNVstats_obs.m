%% calculates chromosome level SNV statistics for obs samples
%% input from bedfile, output to SNVstats folder

%% NB can be parallelized across chromosomes

inFold = '../bedFiles/';
outFold = '../SNVstats/';
fname = [cohortName '.obs.bed'];

for cChr = 1:22
    cChr = num2str(cChr);
    
    snv_ids = cell(0,0);
    snv_refs = [];
    snv_alts = [];
    samp_ids = cell(0,0);
    sampXsnv_cell = cell(0,0);
    
    [dum nLine] = system(['wc -l ' inFold fname]);
    idx = strfind(nLine,' ');
    nLine = nLine(1:idx(1)-1);
    
    fid = fopen([inFold fname]);
    cLine = 0;
    cCount = 0;
    tic
    while(1)
        
        cLine = cLine+1;
        cCount = cCount+1;
        if cCount==10000
            display([num2str(cLine) '/' nLine]);
            cCount = 0;
        end
        
        ln = fgetl(fid);
        if ln==-1
            break;
        end
        if ln(1)=='#'
            continue;
        end
        
        tabidx = strfind(ln,sprintf('\t'));
        a = ln(1:tabidx(1)-1);
        b = ln(tabidx(2)+1:tabidx(3)-1);
        if length(a)==1||length(a)==2
            a = ['chr' a];
        end
        if ~(strcmp(a,['chr' cChr]))
            continue;
        end
        snv_id = [a ':' b];
        
        vec = strcmp(snv_id,snv_ids);
        if sum(vec)==0
            idx1 = length(snv_ids)+1;
            snv_ids{1,idx1} = snv_id;
            snv_refs = [snv_refs ln(tabidx(3)+1)];
            snv_alts = [snv_alts ln(tabidx(4)+1)];
        else
            idx1 = find(vec==1,1);
        end

        samp_id = ln(tabidx(5)+1:tabidx(6)-1);
        
        vec = strcmp(samp_id,samp_ids);
        if sum(vec)==0
            idx2 = length(samp_ids)+1;
            samp_ids{1,idx2} = samp_id;
        else
            idx2 = find(vec==1,1);
        end
        
        if idx2<=length(sampXsnv_cell)
            sampXsnv_cell{idx2} = [sampXsnv_cell{idx2} idx1];
        else
            sampXsnv_cell{idx2} = [idx1];
        end
        
    end
    toc
    fclose(fid);
    
    % get stats
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
    
    fname2 = strrep(fname,['.bed'],['.chr' cChr '.mat']);
    save([outFold fname2],'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
        'N_samp','N_snv','snv_shared','h_shared');
end
