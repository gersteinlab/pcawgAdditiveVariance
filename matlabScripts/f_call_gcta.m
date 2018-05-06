%% calls GCTA for all funseq thresholds

cd('../gctaFiles');

tags = {cohortName};

for cTag = 1:length(tags)
    for cFsq = 0:6
        
        tag = [tags{cTag} '.fsq' num2str(cFsq)];
        
        cmdStr = ['./gcta64 --dosage-mach ' tag '.dose ' tag '.info ' ...
            '--imput-rsq 0.3 --maf 0 --make-grm --out ' tag];
        system(cmdStr);
        
        cmdStr = ['./gcta64 --reml --grm ' tag ' --pheno ' tag '.phen --out ' ...
            tag ' --prevalence 0.5 --reml-pred-rand'];
        system(cmdStr);
        
        cmdStr = ['./gcta64 --dosage-mach ' tag '.dose ' tag '.info --blup-snp ' tag '.indi.blp --out ' ...
            tag];
        system(cmdStr);
        
    end
end

cd('../matlabScripts');
