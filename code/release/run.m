clear
load datasets.mat
addpath(genpath('pogs'))

%% Select expression profile of brain tissue. This is a vector with the same number of elements as the number of vertices. Each entry represents the expression value of a gene
tissue_id = find(strcmp(tissue_names, 'Brain'));
tissue_expression = GTEx_expression(:, tissue_id);

%% Run ActPro. POGS implementation ('ActPro_L1_POGS') is experimental. For better performance, install CVX/Mosek (http://cvxr.com/cvx/doc/mosek.html) and use ''ActPro_L1' method
[ TS_interactome ] = construct_TS_interactome( interactome, tissue_expression, 'method', 'ActPro_L1_POGS' );
