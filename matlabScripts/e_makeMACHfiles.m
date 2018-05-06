%% writes files necessary for GCTA to gctaFiles folder (separate inputs for each funseq threshold)
%% inputs are from keys and bedfiles folders
%% additionally, an intermediate file is saved in the machMats folder, containing the genotype matrices saved as a matlab structure

input_tag = cohortName;
outFold = '../gctaFiles/';

bedFold = '../bedFiles/';
keyFold = '../keys/';
machFoldSNV = '../machMats/';
mergeFold = '../SNVstats/';

load([keyFold input_tag '.orderedKey.mat']);
load([mergeFold input_tag '.obs.null.merged.mat']);

nSNV = length(ordKey_fsq);
nSamp = N_samp;
filt = (snv_shared>1);
filt_idx = find(filt);
nFilt = sum(filt);
if ~(nFilt==length(ordKey_missing))
    display('filter test failed');
end

machMatsSNV = zeros(nSamp,nSNV,6);
machMat0SNV = zeros(nSamp,nSNV);
machMatsSNV_noncd = zeros(nSamp,nSNV,6);
machMat0SNV_noncd = zeros(nSamp,nSNV);
machMatsSNV_nonprm = zeros(nSamp,nSNV,6);
machMat0SNV_nonprm = zeros(nSamp,nSNV);
machMatsSNV_drv = zeros(nSamp,nSNV,6);
machMat0SNV_drv = zeros(nSamp,nSNV);
snv_chr_all = zeros(1,nSNV);
rej = 0;

for cSNV = 1:nSNV
    
    if mod(cSNV,1000)==0
        display([num2str(cSNV) '/' num2str(nSNV)]);
    end
    
    cIdx = filt_idx(cSNV);
    
    sampIdxs = [];
    for cSamp = 1:nSamp
        if ismember(cIdx,sampXsnv_cell{cSamp})
            sampIdxs = [sampIdxs cSamp];
        end
    end
    
    snv_chr = snv_ids{cIdx};
    clnIdx = find(snv_chr==':',1);
    snv_chr = snv_chr(4:clnIdx-1);
    snv_chr = str2num(snv_chr);
    if isempty(snv_chr)
        snv_chr = 23;
    end
    if snv_chr < 1
        snv_chr = 23;
    end
    snv_chr_all(cSNV) = snv_chr;
    
    snv_cd = ordKey_cd(cSNV);
    snv_drv = ordKey_drv(cSNV);
    snv_prm = ordKey_prm(cSNV);
    snv_fsq = ordKey_fsq(cSNV);
    lim = floor(snv_fsq);
    
    if lim>0
        if snv_drv==1
            machMatsSNV_drv(sampIdxs,cSNV,1:lim) = 1;
        else
            machMatsSNV(sampIdxs,cSNV,1:lim) = 1;
            if snv_cd==1
                machMatsSNV_noncd(sampIdxs,cSNV,1:lim) = 1;
                if snv_prm==0
                    machMatsSNV_nonprm(sampIdxs,cSNV,1:lim) = 1;
                end
            end
        end
    else
        if snv_drv==1
            machMat0SNV_drv(sampIdxs,cSNV) = 1;
        else
            machMat0SNV(sampIdxs,cSNV) = 1;
            if snv_cd==1
                machMat0SNV_noncd(sampIdxs,cSNV) = 1;
                if snv_prm==0
                    machMat0SNV_nonprm(sampIdxs,cSNV) = 1;
                end
            end
        end
    end
end

for i = 1:length(samp_ids)
    if ~isempty(strfind(samp_ids{i},'null'))
        samp_ids{i} = strrep(samp_ids{i},'-null','_null');
    end
end

singSamps = zeros(length(samp_ids),1);
phenVec = zeros(length(samp_ids),1);
for i = 1:length(samp_ids)
    if ~isempty(strfind(samp_ids{i},'null'))
        vec = strcmp(samp_ids{i}(1:end-5),samp_ids);
        if sum(vec)==0
            singSamps(i) = 1;
        end
        continue;
    else
        phenVec(i) = 1;
        vec = strcmp([samp_ids{i} '_null'],samp_ids);
        if sum(vec)==0
            singSamps(i) = 1;
        end
    end
end
phenVec = phenVec(singSamps==0);

fname = [machFoldSNV input_tag '.mat'];
save(fname,'machMatsSNV','machMat0SNV','machMatsSNV_noncd','machMat0SNV_noncd',...
    'machMatsSNV_nonprm','machMat0SNV_nonprm','machMatsSNV_drv','machMat0SNV_drv',...
    'samp_ids','rej','singSamps','phenVec','snv_chr_all');

% write out gcta files
output_tag = [outFold input_tag];
cMachMats = machMatsSNV;
cMachMat0 = machMat0SNV + machMatsSNV(:,:,1);

for cFsq = 0:6
    if cFsq==0
        cMat = cMachMat0;
    else
        cMat = cMachMats(:,:,cFsq);
    end
    cMat = cMat(singSamps==0,:);
    cBinMat = (cMat>0);
    colTots = sum(cBinMat,1);
    colInds = find((colTots>1)&(snv_chr_all<=22));
    cMat = cMat(:,colInds);
    nGene = size(cMat,2);
    N_samp = size(cMat,1);
    
    fid = fopen([output_tag '.fsq' num2str(cFsq) '.dose'],'w');
    for i = 1:N_samp
        display([num2str(i) '/' num2str(N_samp)]);
        
        nn = num2str(i);
        
        ln = ['samp' nn ' ALLELE_DOSE'];
        
        for j = 1:nGene
            str = [num2str(cMat(i,j),'%.4f')];
            ln = [ln ' ' str];
        end
        ln = [ln ' ' sprintf('\n')];
        
        fprintf(fid,ln);
    end
    fclose(fid);
    
    fid = fopen([output_tag '.fsq' num2str(cFsq) '.info'],'w');
    ln = ['SNP' sprintf('\t') 'Al1' sprintf('\t') 'Al2' sprintf('\t') ...
        'Freq1' sprintf('\t') 'MAF' sprintf('\t') 'Quality' sprintf('\t') ...
        'Rsq' sprintf('\n')];
    fprintf(fid,ln);
    for i = 1:nGene
        display([num2str(i) '/' num2str(nGene)]);
        
        gen_id = ['snp' num2str(i)];
        maf = 0.5;
        freq = 0.5;
        ln = [gen_id sprintf('\t') 'A' sprintf('\t') ...
            'T' sprintf('\t') num2str(freq,'%.4f') sprintf('\t') ...
            num2str(maf,'%.4f') sprintf('\t') '1.0000' sprintf('\t') '1.0000'...
            sprintf('\n')];
        fprintf(fid,ln);
    end
    fclose(fid);
    
    fid = fopen([output_tag '.fsq' num2str(cFsq) '.phen'],'w');
    nonSings = find(singSamps==0);
    for i = 1:N_samp
        display([num2str(i) '/' num2str(N_samp)]);
        
        ln = ['samp' num2str(i) sprintf('\t')];
        ln = [ln 'samp' num2str(i) sprintf('\t')];
        ii = nonSings(i);
        if isempty(strfind(samp_ids{ii},'null'))
            ln = [ln '1'];
        else
            ln = [ln '0'];
        end
        ln = [ln sprintf('\n')];
        
        fprintf(fid,ln);
    end
    fclose(fid);
end

