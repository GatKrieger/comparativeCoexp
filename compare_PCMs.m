% compare pairwise correlation matrices (PCMs)
%% load data

exp_per_sp_table = ['expression_per_species_test.xlsx'];

load_expression_data;
F = fields(dataS);

%% parameters
global PCMs expCorr F Flabels intint ifFlt corrThrs dataLabels NCorrThrs
ifFlt = 0; % whether to filter the pairwise correlation matrix
corrThrs = 0.2; % if filtering, take only |R| > 'corrThrs'
geneThrs = 10; % filter out genes that are experssed in only 'geneThrs'% of the samples (or less)
NCorrThrs = 20; % min number of genes to compute correlation
expThrs = min(min(dataS.(F{1}))); % minimal expression level
expThrsList = repmat(expThrs, 1, length(F));
%% compare PCM of two datasets
disp('generating Pairwise correlation matrices...');
run PCMs_generate.m;

%% expected correlations (dataset-control)
disp('generating expected correlation values...');
run expCorr_generate.m;

disp(['# genes: ', num2str(length(intint))]);
%% divergece table
[obsP, expP, pvalP] = corr_of_corr(1,2, ifFlt);
ParentsDiv = (expP - obsP)./expP;
[obsH, expH, pvalH] = corr_of_corr(3,4, ifFlt);
HybridDiv = (expH - obsH)./expH;
[obsCerHyc] = corr_of_corr(1,3, ifFlt);
[obsParHyp] = corr_of_corr(2,4, ifFlt);
Ngenes = nan(length(intint), 4);
for i = 1:length(F)
    Ngenes(:, i) = sum(PCMs.(F{i}).PCMflt ~= 0, 2);
end
%% overlap of top correlated genes
NgenesToOL = [3, 10, 50, 100, 200, 500, 1000, 2000, 4000];
if ifFlt; measure = 'PCMflt'; else; measure = 'PCM'; end;
comps = {[1,2], [3,4], [1,3], [2,4]};
compsTitle = {'cerpar', 'hychyp', 'cerhyc', 'parhyp'};
OLs = nan(length(intint), length(NgenesToOL), length(comps));
for i = 1:4
    for j = 1:length(NgenesToOL)
        OLs(:, j, i) = OLngenes(PCMs.(F{comps{i}(1)}).(measure), PCMs.(F{comps{i}(2)}).(measure), NgenesToOL(j));
    end
end
%%
global divT;
divT = table(intint, orfNames(intint), ParentsDiv, HybridDiv, ...
    obsP, obsH, obsCerHyc, obsParHyp, ...
    nanmedian(expCorr.(F{1}).corrDiag, 2), nanmedian(expCorr.(F{2}).corrDiag, 2),...
    nanmedian(expCorr.(F{3}).corrDiag, 2), nanmedian(expCorr.(F{4}).corrDiag, 2), ...
    Ngenes(:, 1), Ngenes(:, 2), Ngenes(:, 3), Ngenes(:, 4), ...
    'variableNames', {'geneID', 'gene','div_cer_par', 'div_hyc_hyp',...
    'obs_cer_par', 'obs_hyc_hyp', 'obs_cer_hyc', 'obs_par_hyp',...
    'exp_cer', 'exp_par', 'exp_hyc', 'exp_hyp', 'Ngenes_cer', 'Ngenes_par', 'Ngenes_hyc','Ngenes_hyp'});

disp('divT is ready');
%% save PCM
name = datestr(datetime);
name = strrep(name(1:11), '-', '_');
ggc = struct;
ggc.PCMs = PCMs;
ggc.expCorr = expCorr;
ggc.intint = intint;
ggc.corrThrs = corrThrs;
ggc.Flabels = Flabels;
ggc.ifFlt = ifFlt;
ggc.F = F;
save(['divT', name, '.mat'], 'divT'); 
save(['ggc_', name, '.mat'], '-struct', 'ggc', '-v7.3');