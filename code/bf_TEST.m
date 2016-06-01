% climada template
% MODULE:
%   drought_fire
% NAME:
%   bf_TEST
% PURPOSE:
%   TEST bushfire module
% CALLING SEQUENCE:
%   bf_TEST
% EXAMPLE:
%   bf_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160701
%-

% HALLO

% next version of our incredible bushfire ^^
%startup climada
%startup; %need to do in folder climada master

% % Parameters
% hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'BFXX_hazard.mat'];
% 
% if ~exist(hazard_set_file,'file')
%     % generate hazard file
%     hazard = climada_bf_hazard_set;
% end

% the folder we store intermediate steps' data
data_folder=[climada_global.data_dir filesep 'hazards' filesep '_data'];
if ~exist(data_folder,'dir'),[fP,fN]=fileparts(data_folder);mkdir(fP,fN);end

bf_file = [data_folder filesep 'AUS_bf.mat'];
load(bf_file)

centroids_file = [climada_global.data_dir filesep 'centroids' filesep 'AUS_bf_centroids.mat'];
load(centroids_file)

hazard        = climada_bf_hazard_set(bf,centroids,1.00,'','_hazard');
hazard_CC_f3  = climada_bf_hazard_set(bf,centroids,1.03,'','_hazard_CC_f3');
hazard_CC_f65 = climada_bf_hazard_set(bf,centroids,1.65,'','_hazard_CC_f65');


% generate the entities
entity_today       = climada_entity_read('Victoria_today_adaptation.xlsx',hazard_CC_f3); %choose Victoria_today..
entity_2030        = climada_entity_read('Victoria_2030_adaptation.xlsx',hazard_CC_f3); %choose Victoria_2030..

%clear entity_today.assets.VALNaN % doesn't work -> if want need to do by 
% hand!

EDS_today = climada_EDS_calc(entity_today,hazard);
EDS_2030 = climada_EDS_calc(entity_2030,hazard);
EDS_2030_CC_f3 = climada_EDS_calc(entity_2030,hazard_CC_f3);
EDS_2030_CC_f65 = climada_EDS_calc(entity_2030,hazard_CC_f65);
% waterfall plot only works if in climada_EDS_stats.m everything commented 
% out with historical data!!!!


figure
res_f3=climada_waterfall_graph(EDS_today,EDS_2030,EDS_2030_CC_f3);

figure
res_f65=climada_waterfall_graph(EDS_today,EDS_2030,EDS_2030_CC_f65);

% Adaptation
% adaptation measures -> implement all of them
res_today_adapt = climada_measures_impact(entity_today,hazard,'no');

%impact of adaptation on AED
% impact of fire fighters
diff_AED_ff = res_today_adapt.ED(end)-res_today_adapt.ED(1);
% impact of education
diff_AED_edu = res_today_adapt.ED(end)-res_today_adapt.ED(2);
% impact of sprinkler system
diff_AED_sprinkler = res_today_adapt.ED(end)-res_today_adapt.ED(3);

figure
climada_adaptation_cost_curve(res_today_adapt)

%show arrows that show which ones cost effective
%climada_adaptation_cost_curve(res_today_adapt,'','','','','','',1)

res_2030_adapt_CC_f3  = climada_measures_impact(entity_2030, hazard_CC_f3, res_today_adapt);
res_2030_adapt_CC_f65 = climada_measures_impact(entity_2030, hazard_CC_f65,res_today_adapt);

climada_adaptation_cost_curve(res_2030_adapt_CC_f3)
climada_adaptation_cost_curve(res_2030_adapt_CC_f65)

% for it to work needed to comment out row 324 in
% climada_adaptation_cost_curve.m: clear clear measures_impact_comparison
climada_adaptation_cost_curve(res_2030_adapt_CC_f3,res_2030_adapt_CC_f65)