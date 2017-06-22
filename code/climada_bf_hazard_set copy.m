function hazard=climada_bf_hazard_set(bf,centroids,fscrew,areascrew,hazard_set_file)
% climada template
% MODULE:
%   drought_fire
% NAME:
%   climada_bf_hazard_set
% PURPOSE:
%   Generate a BF (bushfire) hazard event set. 
%
%   The TEST case for average fires in Victoria, Australia, see code bf_TEST
%
%   Each fire is attributed to 1 centroid (valid for centroids down to 10km
%   resolution, revise if using 1km or finer).
%
%   previous: likely climada_bf_generator.m
%   next: e.g climada_EDS_calc
% CALLING SEQUENCE:
%   hazard = climada_bf_hazard_set(bf,centroids,fscrew,areascrew,hazard_set_file);
% EXAMPLE:
%   hazard = climada_bf_hazard_set;
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   bf: file containing lon, lat, intentsity of each wildfire and number
%       of years simulated (the 'proto-hazard set')
%   centroids:  centroids in my domain (to contribute intensity to)
%   fscrew:     to account for climate change: constant factor I multiply
%               frequency with -> if empty set to 1 by default -> input in
%               format 1+%increase -> eg for 30% increase fscrew = 1.3
%   areascrew:  to account for climate change: constant factor I multiply
%               affected area with -> if empty set to 1 by default -> input
%               in format 1+%increase -> eg for 30% increase areascrwe = 1.3
% OUTPUTS:
% MODIFICATION HISTORY:
% beuschl@student.ethz.ch, 20160601, initial, key author
% horatc@student.ethz.ch, 20160601, initial, key author
% david.bresch@gmail.com, 20160601, climada-compatibility
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('bf','var'),bf=[];end % in case we want to pass all parameters as structure
if ~exist('centroids','var'),centroids=[];end % in case we want to pass all parameters as structure
if ~exist('fscrew','var'),fscrew=[];end
if ~exist('areascrew','var'),areascrew=[];end
if ~exist('hazard_set_file','var'),hazard_set_file=[];end

% locate the module's (or this code's) data folder (usually  afolder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% set default value for fscrew and areascrew if not given
if isempty(fscrew),fscrew=1;end % default=1
if isempty(areascrew),areascrew=1;end % default=1
%
% the folder we store intermediate steps' data
data_folder=[climada_global.data_dir filesep 'hazards' filesep '_data'];
if ~exist(data_folder,'dir'),[fP,fN]=fileparts(data_folder);mkdir(fP,fN);end


% prompt for bf intermediate step data file if not given
if isempty(bf) % local GUI
    bf                   = [data_folder filesep '*.mat'];
    [filename, pathname] = uigetfile(bf, 'Select proto hazard set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        bf=fullfile(pathname,filename);
    end
    load(bf); % ADDED because couldn't load file?!?!
end

if isempty(centroids) % local GUI
    centroids             = [climada_global.data_dir filesep 'centroids' filesep '*.mat'];
    [filename, pathname] = uigetfile(centroids, 'Select centroids:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids=fullfile(pathname,filename);
    end
    load(centroids); % ADDED because couldn't load file?!?!
end

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep 'BFXX_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save bushfire (BF) hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
end

% complete path, if missing
[fP,fN,fE]=fileparts(hazard_set_file);
if isempty(fP),hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep fN fE];end

if isfield(centroids,'assets')
    % centroids contains in fact an entity
    entity=centroids; centroids=[]; % silly switch, but fastest
    centroids.lat =entity.assets.lat;
    centroids.lon=entity.assets.lon;
    centroids.centroid_ID=1:length(entity.assets.lon);
    % treat optional fields
    if isfield(entity.assets,'distance2coast_km'),centroids.distance2coast_km=entity.assets.distance2coast_km;end
    if isfield(entity.assets,'elevation_m'),centroids.elevation_m=entity.assets.elevation_m;end
    if isfield(entity.assets,'country_name'),centroids.country_name=entity.assets.country_name;end
    if isfield(entity.assets,'admin0_name'),centroids.admin0_name=entity.assets.admin0_name;end
    if isfield(entity.assets,'admin0_ISO3'),centroids.admin0_ISO3=entity.assets.admin0_ISO3;end
    if isfield(entity.assets,'admin1_name'),centroids.admin1_name=entity.assets.admin1_name;end
    if isfield(entity.assets,'admin1_code'),centroids.admin1_code=entity.assets.admin1_code;end
    clear entity
end

% template for-loop with waitbar or progress to stdout
t0       = clock;
n_events = size(bf,2);
msgstr   = sprintf('processing %i events @ %i centroids',n_events,length(centroids.lon));
mod_step = 10; % first time estimate after 10 events, then every 100

fprintf('%s\n',msgstr);
format_str='%s';

hazard.intensity = zeros(size(bf,2),size(centroids.centroid_ID,2));
hazard.intensity = sparse(hazard.intensity);

for i=1:n_events
    
    % your calculations here
    %for i=1:5000,sqrt(i)*exp(event_i);end % DUMMY
    mean_bf.lat(i) = mean(bf(i).lat(:));
    mean_bf.lon(i) = mean(bf(i).lon(:));
    mean_bf.intensity(i) =mean(bf(i).intensity(:))*length(bf(i).intensity)*0.1651*areascrew;
    % scaled with the area so that our bigger fire can actually generate
    % more damage ;)
    
    % attribute each fire to a specific centroid
    diff_c_h_lat = abs(centroids.lat-mean_bf.lat(i));
    [min_lat, loc_min_lat] = min(diff_c_h_lat);
    diff_c_h_lon = abs(centroids.lon-mean_bf.lon(i));
    [min_diff_lon, loc_min_lon] = min(diff_c_h_lon);
    
    all_min_lat = find(centroids.lat(:) == centroids.lat(loc_min_lat));
    all_min_lon = find(centroids.lon(:) == centroids.lon(loc_min_lon));
    
    loc_log = ismember(all_min_lat,all_min_lon);
    min_c_ID = all_min_lat(loc_log);
    
    hazard.intensity(i,min_c_ID) = mean_bf.intensity(i);
    hazard.frequency(i) = fscrew*1/bf(i).years_simulated;
    
    % the progress management
    if mod(i,mod_step)==0
        mod_step          = 100;
        t_elapsed_event   = etime(clock,t0)/i;
        events_remaining  = n_events-i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('est. %3.0f sec left (%i/%i events)',t_projected_sec,   i,n_events);
        else
            msgstr = sprintf('est. %3.1f min left (%i/%i events)',t_projected_sec/60,i,n_events);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
    
end % event_i
fprintf(format_str,''); % move carriage to begin of line
fprintf('after the loop\n')

% fill in the other hazard fields
hazard.lon = centroids.lon;
hazard.lat = centroids.lat;
hazard.centroid_ID = centroids.centroid_ID;

hazard.event_ID = 1:size(bf,2);
hazard.orig_event_flag = zeros(1,size(bf,2));
hazard.peril_ID = 'BF';
% hazard.filename = hazard_set_file;

if fscrew~=1
    hazard.comment = 'with climate change';
else
    hazard.comment =  'without climate change';
end

fprintf('saving BF hazard set as %s\n',hazard_set_file)

hazard.filename = hazard_set_file;
save(hazard_set_file,'hazard') % ADDED because didn't save before -> not sure wether right approach ^^

end % climada_bf_hazard_set