%%% additive variance demo

cohortName = 'Breast-AdenoCa';

cd matlabScripts

%%% call matlab pipeline in order
% a_getSNVstats_null
% a_getSNVstats_obs
% a_makeKeys
% b_mergeSNVstats
% c_mergeAll
% d_makeOrderedKey
% e_makeMACHfiles
% f_call_gcta
f_summarize_results

cd ..
