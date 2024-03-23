default = 1000; % fallback bounds
model = readSBML('../external_data/mitocore/mitocore_v1.01.xml', default);
changeCobraSolver('gurobi');
model = changeObjective(model, {'OF_ATP_MitoCore'}, 1);

%% read omics data
all_proteomics = readtable('../external_data/293parp_abundance_ratios_mapped.csv');
control = all_proteomics(:, [1,2]);
mito = all_proteomics(:, [1,3]);

%%
control.Properties.VariableNames{1} = 'gene';
control.Properties.VariableNames{2} = 'value';
mito.Properties.VariableNames{1} = 'gene';
mito.Properties.VariableNames{2} = 'value';

%% map omics data onto reactions
mapped_mito= mapExpressionToReactions(model, mito);
mapped_control = mapExpressionToReactions(model, control);

%% export mapped omics data
writetable(table(model.rxns, mapped_mito), '../generated_data/mapped_mito.csv', 'WriteRowNames', true);
writetable(table(model.rxns, mapped_control), '../generated_data/mapped_control.csv', 'WriteRowNames', true);

%% fill NaNs
fill_val = -1.; 
mapped_mito = fillmissing(mapped_mito, 'constant', fill_val);
mapped_control = fillmissing(mapped_control, 'constant', fill_val);

%% GIMME %%
frac = .8; % objective fraction

threshold = median([
    median(mapped_mito), 
    median(mapped_control)
    ]);

gimme_mito = make_gimme(model, mapped_mito, threshold, frac, '../generated_models/gimme_mito.xml');
gimme_wt= make_gimme(model, mapped_control, threshold, frac, '../generated_models/gimme_control.xml');


%% functions
function gimme_model = make_gimme(model, data, threshold, obj_frac, outfile)
    options = 'options_GIMME';
    options = load(['options_methods' filesep options]);
    options.solver = 'GIMME';
    options.threshold = threshold;
    options.obj_frac = obj_frac;
    options.expressionRxns = data;
    gimme_model = createTissueSpecificModel(model, options);

    writeSBML(gimme_model, outfile);
end