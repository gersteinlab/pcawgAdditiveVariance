%% reads the GCTA outputs, and writes a summary results file in the results folder

tags = {cohortName};
nTag = length(tags);
outFold = '../gctaFiles/';

pvalsAll = [];
advarsAll = [];
tagList = 1:nTag;

for cTag = tagList
    
    input_tag = tags{cTag};
    
    pvals = ones(1,7)*-1;
    addvars = ones(1,7)*-1;
    
    display([input_tag ':']);
    display(' ');
    
    for fsqThres = [0:6]
        
        tag = [input_tag '.fsq' num2str(fsqThres)];
        
        fname = [outFold tag '.hsq'];
        
        if ~exist(fname)
            continue;
        end
        
        fid = fopen(fname);
        while(1)
            ln = fgetl(fid);
            if length(ln)<4
                continue;
            end
            if strcmp(ln(1:4),'Pval')
                break;
            end
        end
        pval = ln(6:end);
        pvals(fsqThres+1) = str2num(pval);
%         if str2num(pval)<1e-4
%             display([num2str(cCondit) ',' num2str(fsqThres) ': ' pval]);
%         end
        fclose(fid);
        
        fid = fopen(fname);
        while(1)
            ln = fgetl(fid);
            if length(ln)<7
                continue;
            end
            if strcmp(ln(1:7),'V(G)/Vp')
                break;
            end
        end
        tabidx = strfind(ln,sprintf('\t'));
        addvars(fsqThres+1) = str2num(ln(tabidx(1)+1:tabidx(2)-1));
        fclose(fid);
        
    end
    
    % nSNPs
    pvals
    addvars
    
    resFile = ['../results/' input_tag '.txt'];
    fid = fopen(resFile,'w');
    fprintf(fid,'fsq_thres\tadditive_var\tp-val\n');
    for fsqThres = [0:6]
        fprintf(fid,[num2str(fsqThres) '\t']);
        fprintf(fid,[num2str(addvars(fsqThres+1)) '\t']);
        fprintf(fid,[num2str(pvals(fsqThres+1)) '\n']);
    end
    fclose(fid);
    
end

