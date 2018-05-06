%% calculates keys (coding/driver/funseq infomation) for null and obs SNVs
%% input from summary files, output to keys folder

for cIt = 1:2
    
    if cIt==1
        input_fname = [cohortName '.null.summary.txt'];
    else
        input_fname = [cohortName '.obs.summary.txt'];
    end
    output_fname = strrep(input_fname,'summary.txt','key.mat');
    inFold = '../summaryFiles/';
    outFold = '../keys/';
    
    fid = fopen([inFold input_fname],'r');
    nLine = 0;
    while(1)
        ln = fgetl(fid);
        if ln==-1
            break;
        end
        nLine = nLine + 1;
    end
    fclose(fid);
    
    snv_ids_key = cell(1,nLine);
    snv_cd_key = zeros(1,nLine);
    snv_drv_key = zeros(1,nLine);
    snv_fsq_key = zeros(1,nLine);
    snv_prm_key = zeros(1,nLine);
    snv_int_key = zeros(1,nLine);
    snv_chr_key = zeros(1,nLine);
    snv_gid_key = zeros(1,nLine);
    geneList = cell(1,nLine);
    
    fid = fopen([inFold input_fname]);
    ln = fgetl(fid);
    cLine = 0;
    cCount = 0;
    cGene = 0;
    while(1)
        
        if cCount==10000
            display([num2str(cLine) '/' num2str(nLine)]);
            cCount = 0;
        end
        
        ln = fgetl(fid);
        if ln==-1
            break;
        end
        
        tabidx = strfind(ln,sprintf('\t'));
        
        snv_id = ln(1:tabidx(1)-5);
        
        if ~(strcmp(snv_id(1:3),'chr'))
            continue;
        end
        
        cLine = cLine+1;
        cCount = cCount+1;
        
        snv_cd = ln(tabidx(1)+1);
        snv_drv = ln(tabidx(2)+1);
        snv_fsq = ln(tabidx(3)+1:tabidx(4)-1);
        snv_prm = ln(tabidx(4)+1);
        snv_int = ln(tabidx(5)+1);
        snv_gen = ln(tabidx(6)+1:end);
        spidx = strfind(snv_gen,'[');
        if ~isempty(spidx)
            spidx = spidx(1);
            if snv_gen(spidx-1)==' '
                spidx = spidx-1;
            end
            snv_gen = snv_gen(1:spidx-1);
        end
        clidx = strfind(snv_id,':'); clidx = clidx(1);
        if strcmp(snv_id(4:clidx-1),'X')
            snv_chr = 23;
        elseif strcmp(snv_id(4:clidx-1),'Y')
            snv_chr = 24;
        elseif strcmp(snv_id(4:clidx-1),'M')
            snv_chr = 25;
        else
            snv_chr = str2num(snv_id(4:clidx-1));
        end
        
        snv_ids_key{cLine} = snv_id;
        if strcmp(snv_cd,'1')
            snv_cd_key(cLine) = 1;
        end
        if strcmp(snv_drv,'1')
            snv_drv_key(cLine) = 1;
        end
        snv_fsq_key(cLine) = str2num(snv_fsq);
        if strcmp(snv_prm,'1')
            snv_prm_key(cLine) = 1;
        end
        if strcmp(snv_int,'1')
            snv_int_key(cLine) = 1;
        end
        snv_chr_key(cLine) = snv_chr;
        
        vec = strcmp(snv_gen,{geneList{1:cGene}});
        idx = find(vec==1,1);
        if ~isempty(idx)
            snv_gid_key(cLine) = idx;
        else
            cGene = cGene + 1;
            geneList{cGene} = snv_gen;
            snv_gid_key(cLine) = cGene;
        end
        
    end
    fclose(fid);
    
    snv_ids_key = {snv_ids_key{1:cLine}};
    snv_cd_key = snv_cd_key(1:cLine);
    snv_drv_key = snv_drv_key(1:cLine);
    snv_fsq_key = snv_fsq_key(1:cLine);
    snv_prm_key = snv_prm_key(1:cLine);
    snv_int_key = snv_int_key(1:cLine);
    snv_chr_key = snv_chr_key(1:cLine);
    snv_gid_key = snv_gid_key(1:cLine);
    geneList = {geneList{1:cGene}};
    
    chr_idx = cell(1,25);
    chr_ids = cell(1,25);
    for i = 1:25
        chr_idx{i} = find(snv_chr_key==i);
        chr_ids{i} = {snv_ids_key{chr_idx{i}}};
    end
    
    save([outFold output_fname],'snv_ids_key','snv_cd_key','snv_drv_key',...
        'snv_fsq_key','snv_prm_key','snv_int_key','snv_chr_key',...
        'snv_gid_key','geneList','chr_idx','chr_ids');
end
