function hazard=bf_generator_large(csv_file,centroids,hazard_set_file)
% climada template
% MODULE:
%   drought_fire
% NAME:
%   bf_generator_large
% PURPOSE:
%   Process the MODIS5 firms data in order to generate large fires.
%   See also bf_generator_small to generate small fires (using a cellular
%   automat)
%
%   Please proceed as follows:
%   0. run bf=bf_generator_large('TEST') to test whether all works
%      This generates a _TEST_BF_hazard_large hazard set
%   1. Download the historic fire data (since Nov 2000-today) as .csv file
%      from https://firms.modaps.eosdis.nasa.gov/download (the full global
%      dataset since 2000 is too large for a single download) and store the
%      .csv file as ....
%   2. run bf_generate_large to process the data (calls firms_read) and to
%      generate the bushfire hazard event set
%
%   previous call: generate an entity, see e.g. climada_nightlight_entity,
%   climada_srtm_entity or climada_GDP_entity (all in
%   https://github.com/davidnbresch/climada_module_country_risk)
%   next call: e.g. climada_EDS_calc
%
%   Citation for both MCD14DL and MCD14ML: This data set was
%   provided by the University of Maryland and NASA FIRMS operated by
%   NASA/GSFC/ESDIS with funding provided by NASA/HQ
%
%   Developers hint: search for bf_generator in the code to integrate
%   cellular automat to (re)fine fires.
%
% CALLING SEQUENCE:
%   hazard=bf_generator_large(csv_file,param2)
% EXAMPLE:
%   hazard=bf_generator_large('TEST') % use the TEST dataset
% INPUTS:
%   csv_file: the name and path of the .csv file as downloaded from
%       https://firms.modaps.eosdis.nasa.gov/download
%       > promted for filename for if not given
%       for TESTS, there is a .csv file for Australia in the drought_fire
%       module's data folder, just set csv_file='TEST'. In this case, test
%       centroids are used, too.
%   centroids: a centroids structure (see climada_centroids_read) or an
%       entity structure (in which case entity.assets.lon and .lat are used)
%       > promted for if empty
%   hazard_set_file: the filename (and path, optional) of the hazard
%       event set. If no path provided, default path ../data/hazards is used
%       (and name can be without extension .mat)
%       > promted for if empty
% OPTIONAL INPUT PARAMETERS:
%   param2: as an example
% OUTPUTS:
%   hazard: a struct, the hazard event set, more for tests, since the
%       hazard set is stored as hazard_set_file, see code
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       intensity(event_i,centroid_i), sparse: the hazard intensity of
%           event_i at centroid_i
%       frequency(event_i): the frequency of each event
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril, e.g. 'TC' for
%       tropical cyclone or 'ET' for extratropical cyclone
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one
%       event_ID: a unique ID for each event
%       date: the creation date of the set
%       matrix_density: the density of the sparse array hazard.intensity
%       windfield_comment: a free comment, not in all hazard event sets
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
% MODIFICATION HISTORY:
% david.bresch@gmail.com, 20160703
% david.bresch@gmail.com, 20170508, hint to bf_generator added
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('csv_file','var'),csv_file=[];end % OR:
if ~exist('centroids','var'),centroids=[];end % OR:
if ~exist('hazard_set_file','var'),hazard_set_file=[];end

% locate the module's (or this code's) data folder (usually  afolder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% since we store the hazard as sparse array, we need an a-priory estimation
% of it's density
hazard_arr_density=0.03; % 3% sparse hazard array density (estimated)
% define the reference year for this hazard set
hazard_reference_year = climada_global.present_reference_year; % default for present hazard is normally 2015
%
% define the scenario name for this hazard set
% we assume no climate change when creating a hazard set from tc tracks
hazard_scenario = 'no climate change';
%
% whether we create the yearset (=1, grouping events into years) or not (=0)
% the yearset is only produced for original tracks, since probabilistic
% ones can be identified as following original indices +1:ens_size, with
% ens_size=(hazard.event_count/hazard.orig_event_count)-1; see climada_EDS2YDS
% Note: the yearset creation assumes tracks to be ordered by ascending year
% (that's the case for UNISYS tracks as read by climada_tc_read_unisys_database)
create_yearset=0; % default=1 (not implemented yet)
%
% TEST files
TEST_csv_file=[module_data_dir filesep 'hazards' filesep 'external_model_output' filesep 'firms.csv'];
TEST_centroids_file=[module_data_dir filesep 'centroids' filesep 'AUS_BF_centroids.mat'];
TEST_hazard_set_file='_TEST_AUS_BF_hazard_large';
TEST_event_i=3839;
%TEST_centroids_file=[module_data_dir filesep 'AUS_Australia_Victoria_entity'];


% prompt for csv_file if not given
if isempty(csv_file) % local GUI
    csv_file=[climada_global.data_dir filesep '*.csv'];
    [filename, pathname] = uigetfile(csv_file, 'Select firms .csv file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        csv_file=fullfile(pathname,filename);
    end
end

% special treatment for TEST case
if strcmp(csv_file,'TEST')
    TEST_mode=1;
    csv_file       =TEST_csv_file;
    centroids      =TEST_centroids_file;
    hazard_set_file=TEST_hazard_set_file;
    if exist(csv_file,'file')
        fprintf('TEST mode, using %s\n',TEST_csv_file);
    else
        fprintf('ERROR: %s not found, consider fetching from GitHub\n',TEST_csv_file);
        fprintf('ERROR: Test data %s not found\n',TEST_csv_file);
        fprintf(['Please download from <a href="https://github.com/davidnbresch/climada_module_drought_fire">'...
            'climada_module_drought_fire</a> from Github.\n'])
        return
    end
else
    TEST_mode=0;
end % strcmp(csv_file,'TEST')

% check (eventually prompt for) centroids
% (climada_centroids_load does also convert an entity into  centroids)
centroids=climada_centroids_load(centroids);

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'BF_hazard_large.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save hazard event set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file=fullfile(pathname,filename);
    end
end

% complete path, if missing
[fP,fN,fE]=fileparts(hazard_set_file);
if isempty(fP),fP=[climada_global.data_dir filesep 'hazards'];end
if isempty(fE),fE='.mat';end
hazard_set_file=[fP filesep fN fE];


% read the firms data
% -------------------
firms=firms_read(csv_file);
% lat: [92736x1 double]
% lon: [92736x1 double]
% brightness: [92736x1 double]
% datenum: [92736x1 double]
% datenum_unique: [6919x1 double]
% unique_IA: [6919x1 double]
% unique_IC: [92736x1 double]

hazard.source_filename=csv_file;

n_events=length(firms.datenum_unique);
n_centroids=length(centroids.lon);

% figure which fires are within the box around centroids
bbox=[min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)];
dx=min(diff(unique(centroids.lon)));dy=min(diff(unique(centroids.lat)));
x=[bbox(1)-dx bbox(1)-dx bbox(2)+dx bbox(2)+dx bbox(1)-dx];
y=[bbox(3)-dy bbox(4)+dy bbox(4)+dy bbox(3)-dy bbox(3)-dy];

% restrict firms to points within
firms_in=inpolygon(firms.lon,firms.lat,x,y);

if TEST_mode % TEST plot
    figure('Name','TEST','Position',[576 198 791 475]);subplot(1,2,1);
    event_i=TEST_event_i; % max number of points for event_i=357;
    event_pos= firms.datenum==firms.datenum_unique(event_i);
    plot(firms.lon,firms.lat,'.k','MarkerSize',1); hold on % all black
    xlim([bbox(1)-5*dx bbox(2)+5*dx]);ylim([bbox(3)-5*dy bbox(4)+5*dy]);axis equal
    plot(firms.lon(firms_in),firms.lat(firms_in),'.r','MarkerSize',1); % in red
    plot(centroids.lon,centroids.lat,'ob','MarkerSize',3) % centroids blue circels
    plot(firms.lon(event_pos),firms.lat(event_pos),'xg'); hold on
    plot(x,y,'-m');
    legend('fires','within','centroids','largest event','bounding box')
    climada_plot_world_borders(1,'','',1);
    title(sprintf('largest event (%i)',event_i));set(gcf,'Color',[1 1 1])
end % TEST plot

% fill the hazard structure
min_year = str2double(datestr(min(firms.datenum_unique),'yyyy'));
max_year = str2double(datestr(min(firms.datenum_unique),'yyyy'));
orig_years = max_year - min_year+1;
hazard.peril_ID         = 'BF';
hazard.reference_year   = hazard_reference_year;
hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
hazard.orig_years       = orig_years;
hazard.orig_event_count = n_events;
hazard.event_count      = n_events;
hazard.event_ID         = 1:hazard.event_count;
hazard.orig_event_flag  = ones(1,n_events);
hazard.datenum          = firms.datenum;
hazard.scenario         = hazard_scenario;

% allocate the hazard array (sparse, to manage MEMORY)
intensity = spalloc(n_events,n_centroids,...
    ceil(n_events*n_centroids*hazard_arr_density));
%intensity = zeros(n_tracks,n_centroids); % FASTER
category  = zeros(1,n_events); % storing sum of brightness
zero_event= zeros(1,n_centroids);

cos_centroids_lat = cos(centroids.lat/180*pi); % calculate once for speedup

% now restrict to points within 'reach' of centroids (speeds up interpolation)
firms.lon=firms.lon(firms_in);
firms.lat=firms.lat(firms_in);
firms.brightness=firms.brightness(firms_in);
firms.datenum=firms.datenum(firms_in);

% template for-loop with waitbar or progress to stdout
t0       = clock;
fprintf('processing %i events @ %i centroids\n',n_events,n_centroids);
mod_step = 10; % first time estimate after 10 events, then every 100
format_str='%s';

centroids_lon=centroids.lon; % for speedup
centroids_lat=centroids.lat; % for speedup
firms_datenum=firms.datenum; % for speedup
firms_datenum_unique=firms.datenum_unique; % for speedup
firms_brightness=firms.brightness; % for speedup
firms_lon=firms.lon;
firms_lat=firms.lat;

for event_i=1:n_events
    
    % find all points for a given event
    event_pos= firms_datenum==firms_datenum_unique(event_i);
    
    if sum(event_pos)>0 % event within centroids
        
        % make explicit for speedup (less indexing)
        lon=firms_lon(event_pos);
        lat=firms_lat(event_pos);
        brightness=firms_brightness(event_pos);
        category(event_i)=sum(brightness);
        
        event_intensity=zero_event; % init
        
        % assign burning points to centroids (currently one point to nearest centroid)
        for point_i=1:length(lon)
            
            % find closest centroid (these two lines MOST TIME CONSUMING)
            dd=((centroids_lon-lon(point_i)).*cos_centroids_lat).^2+(centroids_lat-lat(point_i)).^2; % in km^2
            [~,pos] = min(dd);
            pos=pos(1);
            event_intensity(pos)=max(event_intensity(pos),brightness(point_i));
            
        end % point_i
        
        intensity(event_i,:)= event_intensity;
        
        % here, one might consider to invoke bf_generator
        % in order to generate small (sub)scale fires and
        % higher resolution events using a cellular automat.
         
    end % sum(event_pos)>0
    
    % the progress management
    if mod(event_i,mod_step)==0
        mod_step          = 100;
        t_elapsed_event   = etime(clock,t0)/event_i;
        events_remaining  = n_events-event_i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('est. %3.0f sec left (%i/%i events)',t_projected_sec,   event_i,n_events);
        else
            msgstr = sprintf('est. %3.1f min left (%i/%i events)',t_projected_sec/60,event_i,n_events);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
    
end % event_i
fprintf(format_str,''); % move carriage to begin of line

hazard.intensity=intensity; clear intensity % for speedup above
hazard.category=category; clear category % for speedup above

if TEST_mode % TEST plot
    subplot(1,2,2);climada_hazard_plot(hazard,-1);
end

if create_yearset
    % not implemented yet
end % create_yearset

fprintf('saving hazard as %s\n',hazard_set_file)
save(hazard_set_file,'hazard');

end % bf_generator_large