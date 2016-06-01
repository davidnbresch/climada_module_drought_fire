% bf_TEST: bushfire TEST batch file 
% MODULE:
%   drought_fire
% NAME:
%   bf_TEST
% PURPOSE:
%   TEST bushfire module
%   
%   Uses TEST data within the module's data folder (and subfolders). For
%   operational use, use the standard climada data folder (see
%   climada_global.data_dir), i.e. copy centrpids, entities etc. to the
%   respective sub-folders there.
%
%   please refer to the main code climada_bf_hazard_set
% CALLING SEQUENCE:
%   bf_TEST
% EXAMPLE:
%   bf_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
% MODIFICATION HISTORY:
% beuschl@student.ethz.ch, 20160601, initial, key author
% horatc@student.ethz.ch, 20160601, initial, key author
% david.bresch@gmail.com, 20160601, climada-compatibility
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables


% Parameters
%
% find the local path to the module
module_data_dir=[fileparts(fileparts(which('bf_TEST'))) filesep 'data'];
%
% define the path and filename of the intermediate step data file
% (generated by bf_generator)
bf_file=[module_data_dir filesep 'hazards' filesep 'external_model_output' filesep 'AUS_BF_proto_data.mat'];
%
% define the path and filename of the centroids file
centroids_file = [module_data_dir filesep 'centroids' filesep 'AUS_BF_centroids.mat'];

entity_today_file=[module_data_dir filesep 'entities' filesep 'Victoria_today_adaptation.xlsx'];
entity_2030_file=[module_data_dir filesep 'entities' filesep 'Victoria_2030_adaptation.xlsx'];

if ~exist(bf_file,'file')
    [fP,fN]=fileparts(bf_file);
    fprintf('ERROR: %s not found\n',bf_file)
    fprintf('-> please locate intermediate step data file %s, edit %m\n',fN)
    return
else
    load(bf_file)
end

if ~exist(centroids_file,'file')
    [fP,fN]=fileparts(centroids_file);
        fprintf('ERROR: %s not found\n',centroids_file)
    fprintf('-> please locate the centroids file %s, edit %m\n',fN)
    return
else
    load(centroids_file)
end

if ~exist(entity_today_file,'file')
    [fP,fN]=fileparts(entity_today_file);
        fprintf('ERROR: %s not found\n',entity_today_file)
    fprintf('-> please locate the entity file %s, edit %m\n',fN)
    return
end

if ~exist(entity_2030_file,'file')
    [fP,fN]=fileparts(entity_2030_file);
        fprintf('ERROR: %s not found\n',entity_2030_file)
    fprintf('-> please locate the future entity file %s, edit %m\n',fN)
    return
end



% generate the hazard sets (_ in filename to easily find them)
hazard_today        = climada_bf_hazard_set(bf,centroids,1.00,'','_TEST_AUS_Australia_BF_today');
hazard_CC_f3        = climada_bf_hazard_set(bf,centroids,1.03,'','_TEST_AUS_Australia_BF_CC_f3');
hazard_CC_f65       = climada_bf_hazard_set(bf,centroids,1.65,'','_TEST_AUS_Australia_BF_CC_f65');

% read and encode the entities
entity_today = climada_entity_read(entity_today_file,hazard_today);
entity_2030  = climada_entity_read(entity_2030_file, hazard_today);

% calculate risk today, + ecnomic development, + moderate/high climate change (CC)
EDS_today = climada_EDS_calc(entity_today,hazard_today);
EDS_2030 = climada_EDS_calc(entity_2030,hazard_today);
EDS_2030_CC_f3 = climada_EDS_calc(entity_2030,hazard_CC_f3);
EDS_2030_CC_f65 = climada_EDS_calc(entity_2030,hazard_CC_f65);

figure;climada_waterfall_graph(EDS_today,EDS_2030,EDS_2030_CC_f3); title('moderate climate change')
figure;climada_waterfall_graph(EDS_today,EDS_2030,EDS_2030_CC_f65);title('high climate change')

% Adaptation
% ----------
% first to establish the baseline, today
res_today_adapt = climada_measures_impact(entity_today,hazard_today,'no');

% impact of adaptation on AED, just for understanding, as shown in
% adaptation cost curve
%
% impact of fire fighters
diff_AED_ff = res_today_adapt.ED(end)-res_today_adapt.ED(1);
% impact of education
diff_AED_edu = res_today_adapt.ED(end)-res_today_adapt.ED(2);
% impact of sprinkler system
diff_AED_sprinkler = res_today_adapt.ED(end)-res_today_adapt.ED(3);

%figure;climada_adaptation_cost_curve(res_today_adapt)

res_2030_adapt_CC_f3  = climada_measures_impact(entity_2030, hazard_CC_f3, res_today_adapt);
res_2030_adapt_CC_f65 = climada_measures_impact(entity_2030, hazard_CC_f65,res_today_adapt);

figure;climada_adaptation_cost_curve(res_2030_adapt_CC_f3)
figure;climada_adaptation_cost_curve(res_2030_adapt_CC_f65)

% for it to work needed to comment out row 324 in
% climada_adaptation_cost_curve.m: clear clear measures_impact_comparison
figure;climada_adaptation_cost_curve(res_2030_adapt_CC_f3,res_2030_adapt_CC_f65)