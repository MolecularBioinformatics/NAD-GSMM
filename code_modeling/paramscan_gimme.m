% scan using giime for different thresholds

%% read model
default = 1000; % fallback bounds
model = readSBML('../external_data/mitocore/mitocore_v1.01.xml', default);
changeCobraSolver('gurobi');
model = changeObjective(model, {'OF_ATP_MitoCore'}, 1);

%% read omics data
all_proteomics = readtable('../external_data/293parp_abundance_ratios_mapped.csv');
data_control = all_proteomics(:, [1,2]);
data_mito = all_proteomics(:, [1,3]);

data_mito.Properties.VariableNames{1} = 'gene';
data_mito.Properties.VariableNames{2} = 'value';
data_control.Properties.VariableNames{1} = 'gene';
data_control.Properties.VariableNames{2} = 'value';

%% map omics data onto reactions
mapped_mito = mapExpressionToReactions(model, data_mito);
mapped_control = mapExpressionToReactions(model, data_control);

fill_val = -1.; 
mapped_mito = fillmissing(mapped_mito, 'constant', fill_val);
mapped_control = fillmissing(mapped_control, 'constant', fill_val);
%% GIMME: scan objective fraction

threshold = median([
    median(mapped_mito) 
    median(mapped_control) 
    ]);
omics = mapped_control;
out = '../generated_models/models_paramscan_gimme/control_objfrac_%0.2f.xml';

for i=1:100
    of = 1 * (i/100); % current objective fraction
    
    make_gimme(model, omics, threshold, of, sprintf(out, of));
end

omics = mapped_mito;
out = '../generated_models/models_paramscan_gimme/mito_objfrac_%0.2f.xml';

for i=1:100
    of = 1 * (i/100); % current objective fraction
    
    make_gimme(model, omics, threshold, of, sprintf(out, of));
end

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
