function [entity,hazard,params]=climadacrop_init(params)
% climada crop setup
% MODULE:
%   drought_fire
% NAME:
%   climadacrop_init
% PURPOSE:
%   setup the climada crop module, generate the hazard set (hazard) the
%   asset and damage function database (entity) from netCDF and .csv files
%
%   DESCRIPTION PENDING
%
%   next call: climada_EDS_calc(entity,hazard)
% CALLING SEQUENCE:
%   [entity,hazard,params]=climadacrop_init(params)
% EXAMPLE:
%   [entity,hazard]=climadacrop_init; % set it all up
%   climada_EDS_DFC(climada_EDS_calc(entity,hazard)) % plot risk curve
%   for all damage function sampled:
%   for i=1:max(entity.damagefunctions.DamageFunID)
%     entity.assets.DamageFunID=entity.assets.DamageFunID*0+i;
%     EDS(i)=climada_EDS_calc(entity,hazard,'',0,2);
%   end
%   params=climadacrop_init('params') % get default parameters
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields (see also ISO3='params' above):
%    hazard_filename: filename of the .nc file with the crop data
%       if no path is provided, the default climadacrop path is assumed
%       (=[climada_global.data_dir,'climadacrop']).
%    hazard_varname: the varialbe name of in the netCFD file,
%       default='days above 32 degrees'
%    hazard_set_filename: the name the climada hazard set gets stored to
%       default=[climada_global.hazards_dir 'climadacrop_' params.peril_ID '.mat']
%    exposure_filename: filename of the .nc file with the exposure data
%       if no path is provided, the default climadacrop path is assumed
%       (=[climada_global.data_dir,'climadacrop']).
%    exposure_varname: the varialbe name of in the netCFD file,
%       default='days above 32 degrees'
%    damfun_filename: the filename of the .csv file containing the damage
%       function(s)
%    check_plot: whether show a check plot (=1), or not (=0, default)
%       if =2, run a test calculation, too, and show resulting damage
%       frequency curve (DFC). If=-2, only test calc and DFC, no other plot
%    add_country_ISO3: whether show add country ISO3 code to each gridcell
%       (default=1)
%    distance_to_coast: whether we calculate distant to coast (in km) of
%       call centroids (default=1), which speeds up later climada calculations
%       substantially (as coastal hazards need not to be evaluated at
%       inner-continental points). Set =0 in special cases, as initial
%       calculation might easily take some time, but since
%       climada_distance2coast_km listens to climada_global.parfor, set
%       climada_global.parfor=1 for speedup).
%    peril_ID: the 2-digit peril ID, default='AT'
%    hazard_units: the units of the hazard intensity (used in hazard.units
%       and entity.damagefunctions.Intensity_unit
% OUTPUTS:
%   entity: a climada entity structure, see climada_entity_read for a full
%       description of all fields
%    PLUS the fields
%       entity.assets.ISO3_pos: the country number (ISO3_list(:,2))
%       entity.assets.ISO3_list: the list linking isimip country numbers
%           in ISO3_list(:,2) with ISO names in ISO3_list(:,1)
%       Note that entity.assets.centroid_index is already pointing to the
%           correct centroids in the corresponding hazard set as generated
%           in this code, too.
%   hazard: a climada hazard set, see documentation
%    PLUS the field params, a copy of params for this code
%   params: the params structure, as on input, just completed.
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170210, initial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('params','var'), params = struct;end

if strcmpi(params,'params')
    return_params=1;clear params
    params = struct;
else
    return_params=0; % init
end

% check for some parameter fields we need
if ~isfield(params,'hazard_filename'),    params.hazard_filename='';end
if ~isfield(params,'hazard_varname'),     params.hazard_varname='';end
if ~isfield(params,'hazard_set_filename'),params.hazard_set_filename='';end
if ~isfield(params,'exposure_filename'),  params.exposure_filename='';end
if ~isfield(params,'exposure_varname'),   params.exposure_varname='';end
if ~isfield(params,'damfun_filename'),    params.damfun_filename='';end
if ~isfield(params,'check_plot'),         params.check_plot=[];end
if ~isfield(params,'add_country_ISO3'),   params.add_country_ISO3=[];end
if ~isfield(params,'distance_to_coast'),  params.distance_to_coast=[];end
if ~isfield(params,'centroids_file'),     params.centroids_file='';end % output only
if ~isfield(params,'peril_ID'),           params.peril_ID='';end
if ~isfield(params,'hazard_units'),       params.hazard_units='';end

% PARAMETERS
%
verbose=1; % default=1, to suppress output to stdout later
%
climadacrop_data_dir_name='climadacrop'; % without path
% define the defaut folder for climadacrop data
climadacrop_data_dir=[climada_global.data_dir filesep climadacrop_data_dir_name];
if ~isdir(climadacrop_data_dir)
    mkdir(climada_global.data_dir,climadacrop_data_dir_name); % create it
    fprintf('NOTE: store your clima(da)crop input data in %s\n',climadacrop_data_dir);
end
%
% the entity template to construct the entity from
entity_template=[climada_global.entities_dir filesep 'entity_template'];
%
% admin0 shape file (for add_country_ISO3 option):
admin0_shape_file=climada_global.map_border_file;
%
% populate default parameters in params
if isempty(params.check_plot),        params.check_plot=0;end
if isempty(params.add_country_ISO3),  params.add_country_ISO3=0;end  %=1
if isempty(params.distance_to_coast), params.distance_to_coast=1;end %=1
if isempty(params.peril_ID),          params.peril_ID='AT';end
if isempty(params.hazard_units),      params.hazard_units='Tcrit32';end

if isempty(params.hazard_set_filename),params.hazard_set_filename=...
        [climada_global.hazards_dir filesep 'climadacrop_' params.peril_ID '.mat'];end
%
if isempty(params.hazard_filename),  params.hazard_filename  =[climadacrop_data_dir filesep 'Tcrit32.nc'];end
if isempty(params.hazard_varname),   params.hazard_varname  ='days above 32 degrees';end
if isempty(params.exposure_filename),params.exposure_filename=[climadacrop_data_dir filesep 'Maize_exposure.nc'];end
if isempty(params.exposure_varname), params.exposure_varname='yield area';end
%if isempty(params.damfun_filename), params.damfun_filename='damage_function_best.csv';end
%if isempty(params.damfun_filename), params.damfun_filename='damage_functions.csv';end
if isempty(params.damfun_filename), params.damfun_filename='damage_function_sampled.csv';end

% complete path, if missing:
[fP,fN,fE]=fileparts(params.hazard_set_filename);
if isempty(fP),fP=climada_global.hazards_dir;end
if isempty(fE),fE='.mat';end
params.hazard_set_filename=[fP filesep fN fE];

% prepend climadacrop_data_dir in case only filenames are passed
if ~exist(params.hazard_filename,'file'),params.hazard_filename=[climadacrop_data_dir filesep params.hazard_filename];end
if ~exist(params.exposure_filename,'file'),params.exposure_filename=[climadacrop_data_dir filesep params.exposure_filename];end
if ~exist(params.damfun_filename,'file'),params.damfun_filename=[climadacrop_data_dir filesep params.damfun_filename];end

if return_params==1,entity=params;return;end % special case, return the full params structure

% read the hazard data and convert into a climada hazard event set
% ----------------------------------------------------------------

nc_hazard.info = ncinfo(params.hazard_filename);
nc_hazard.lon  = ncread(params.hazard_filename,'longitude')';
nc_hazard.lat  = ncread(params.hazard_filename,'latitude')';
nc_hazard.time = ncread(params.hazard_filename,'z')'; % for time in years, in fact
nc_hazard.data = ncread(params.hazard_filename,params.hazard_varname);

% create the grid
if verbose,fprintf('creating regular grid (meshgrid) ...');end
[gridlon,gridlat] = meshgrid(nc_hazard.lon,nc_hazard.lat);
%gridlon=gridlon';gridlat=gridlat'; % still as grid
vectlon=reshape(gridlon,[1 numel(gridlon)]); % as 1-D vect
vectlat=reshape(gridlat,[1 numel(gridlat)]);
if verbose,fprintf(' done\n');end

if params.distance_to_coast
    fprintf('> calculating distance to coast of %i land points, takes time...\n',length(vectlat))
    distance2coast_km=climada_distance2coast_km(vectlon,vectlat);
end

% no data in hazard
n_events=length(nc_hazard.time);
n_centroids=length(vectlon);
% make a local copy of all the small fields
hazard.peril_ID       = params.peril_ID;
hazard.units          = params.hazard_units;
hazard.reference_year = nc_hazard.time(1);
hazard.lon            = vectlon;
hazard.lat            = vectlat;
if exist('distance2coast_km','var'),hazard.distance2coast_km=distance2coast_km;end
hazard.centroid_ID    = 1:n_centroids;
hazard.orig_years     = n_events;
hazard.orig_event_count = n_events;
hazard.event_count    = length(nc_hazard.time);
hazard.event_ID       = 1:n_events;
hazard.orig_event_flag= ones(1,n_events);
hazard.yyyy           = nc_hazard.time; % comes as year
hazard.mm             = repmat(12,1,n_events);
hazard.dd             = repmat(31,1,n_events);
hazard.date           = datestr(now);
hazard.comment        = sprintf('from %s',params.hazard_filename);
hazard.filename       = params.hazard_set_filename;

% fill in the hazard intensity, one timestep at a time
for event_i=1:n_events
    hazard.intensity(event_i,:)=sparse(reshape(nc_hazard.data(:,:,event_i)',[1 numel(nc_hazard.data(:,:,event_i))]));
    hazard.name{event_i}=sprintf('%i',nc_hazard.time(event_i));
end
hazard.fraction=spones(hazard.intensity); % fraction 100%
hazard.frequency=ones(1,n_events)/n_events;
params.mfilename=mfilename;
hazard.params=params;
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);

if verbose,fprintf('saving hazard as %s\n',hazard.filename);end
save(hazard.filename,'hazard')

% read the exposure data and convert into a climada entity
% --------------------------------------------------------

nc_exposu.info = ncinfo(params.exposure_filename);
nc_exposu.lon  = ncread(params.exposure_filename,'longitude')';
nc_exposu.lat  = ncread(params.exposure_filename,'latitude')';
nc_exposu.data = ncread(params.exposure_filename,params.exposure_varname)';

if abs(sum(nc_exposu.lon-nc_hazard.lon))+abs(sum(nc_exposu.lat-nc_hazard.lat))>0
    fprintf('WARNING: hazard and exposure grids do not match\n');
end

entity=climada_entity_load(entity_template);

entity.assets.filename=[climada_global.entities_dir filesep 'climadacrop_entity'];
entity.assets.lon=vectlon; % constrain to land and convert to 1D
entity.assets.lat=vectlat; % constrain to land and convert to 1D
entity.assets.centroid_index=1:length(entity.assets.lat);
entity.assets.Value=reshape(nc_exposu.data,[1 numel(nc_exposu.data)]); % as 1-D vect

if params.add_country_ISO3
    % read admin0 (country) shape file
    admin0_shapes=climada_shaperead(admin0_shape_file);
    
    entity.assets.ISO3_pos=entity.assets.lon*0; % init
    
    next_shape=1;
    fprintf('assigning ISO3 country IDs to gridcells ...');
    for shape_i=1:length(admin0_shapes)
        ISO3_pos=climada_inpolygon(entity.assets.lon,entity.assets.lat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y,0);
        if ~isempty(ISO3_pos) % country found
            entity.assets.ISO3_list{next_shape,1}=admin0_shapes(shape_i).ISO_A3;
            entity.assets.ISO3_list{next_shape,2}=shape_i;
            entity.assets.ISO3_pos(ISO3_pos)=shape_i;
            next_shape=next_shape+1;
        end
    end % shape_i
    fprintf(' done\n');
end % params.add_country_ISO3

% complete entity
entity.assets.Cover      =entity.assets.Value;
entity.assets.Deductible =entity.assets.Value*0;
entity.assets.DamageFunID=ones(1,length(entity.assets.Value));
entity.assets.Category_ID=entity.assets.DamageFunID;
entity.assets.Region_ID  =entity.assets.DamageFunID;

if isfield(entity.assets,'Value_unit'),entity.assets=rmfield(entity.assets,'Value_unit');end
if isfield(entity.assets,'hazard'),entity.assets=rmfield(entity.assets,'hazard');end

% add damagefunction (imported from .csv file)
damfun_data=climada_csvread(params.damfun_filename);
n_points=length(damfun_data.Tcritdays);
if isfield(entity,'damagefunctions'),entity=rmfield(entity,'damagefunctions');end
if isfield(damfun_data,'damage'),damfun_data.sample_001=damfun_data.damage;end % old naming
entity.damagefunctions.filename=params.damfun_filename;
entity.damagefunctions.DamageFunID=ones(n_points,1);
entity.damagefunctions.Intensity=damfun_data.Tcritdays';
if max(entity.damagefunctions.Intensity)<max(max(hazard.intensity))
    fprintf('Warning: damage function intensity range (%2.2f..%2.2f) smaller than hazard intensity (0..%2.2f) range\n',...
        min(entity.damagefunctions.Intensity),max(entity.damagefunctions.Intensity),full(max(max(hazard.intensity))));
end
if max(damfun_data.sample_001)>2
    damfun_data_fact=1/100;
    fprintf('Note: damage function MDD converted to decimal from percent, new min/max: %2.3f/%2.3f\n',...
        damfun_data_fact*min(damfun_data.sample_001),damfun_data_fact*max(damfun_data.sample_001));
end
entity.damagefunctions.MDD=damfun_data_fact*damfun_data.sample_001';

n_samples=0; % init
DamageFunID=zeros(1,101); % to start with
% more than one sample, stack it up
% parse all names, figure the number of samples
% (to pre-allocate for speedup)
sample_names=fieldnames(damfun_data);
for sample_i=1:length(sample_names)
    if strfind(sample_names{sample_i},'sample')
        DamageFunID(sample_i)=str2double(strrep(sample_names{sample_i},'sample_',''));
        if ~isnan(DamageFunID(sample_i)),n_samples=n_samples+1;end
    end
end % sample_i
% allocate
entity.damagefunctions.DamageFunID=zeros(n_samples*n_points,1);
entity.damagefunctions.Intensity  =zeros(n_samples*n_points,1);
entity.damagefunctions.MDD        =zeros(n_samples*n_points,1);
entity.damagefunctions.PAA        = ones(n_samples*n_points,1);
pos1=1;pos2=n_points;
if params.check_plot,figure('Name','Damage function (samples)','Color',[1 1 1]);end
% 2nd round, now populating the data
if n_samples>1,fprintf('%i damage function samples, each %i points\n',n_samples,n_points);end
for ID_i=1:length(DamageFunID)
    if DamageFunID(ID_i)>0
        entity.damagefunctions.DamageFunID(pos1:pos2)=DamageFunID(ID_i);
        entity.damagefunctions.Intensity(pos1:pos2)=damfun_data.Tcritdays';
        entity.damagefunctions.MDD(pos1:pos2)=damfun_data_fact*damfun_data.(sample_names{ID_i})';
        entity.damagefunctions.MDD(pos1:pos2)=sort(entity.damagefunctions.MDD(pos1:pos2)); % sort ascending
        if params.check_plot
            plot(entity.damagefunctions.Intensity(pos1:pos2),entity.damagefunctions.MDD(pos1:pos2),'-b');hold on
        end

        pos1=pos2+1;pos2=pos2+n_points; % point to next free range
    end 
end % ID_i
% add further fields
entity.damagefunctions.peril_ID=repmat({hazard.peril_ID},size(entity.damagefunctions.MDD));
entity.damagefunctions.Intensity_unit=repmat({params.hazard_units},size(entity.damagefunctions.MDD));
entity.damagefunctions.name=repmat({'Tcritdays'},size(entity.damagefunctions.MDD));
entity.damagefunctions.datenum=repmat(now,size(entity.damagefunctions.MDD));

entity.assets = climada_assets_complete(entity.assets);

if exist('distance2coast_km','var'),entity.assets.distance2coast_km=distance2coast_km;end

% add the source files
entity.assets.exposure_filename =params.exposure_filename;

if verbose,fprintf('saving entity as %s\n',entity.assets.filename);end
save(entity.assets.filename,'entity','-v7.3'); % -v7.3 for size...

if params.check_plot>0,figure;climada_entity_plot(entity);end

if abs(params.check_plot)>1
    fprintf('sampling %i EDSs ...\n',n_samples);
    %n_samples=length(unique((entity.damagefunctions.DamageFunID));
    climada_progress2stdout    % init
    for sample_i=1:n_samples
        entity.assets.DamageFunID=entity.assets.DamageFunID*0+sample_i;
        EDS(sample_i)=climada_EDS_calc(entity,hazard,'',0,2);
        climada_progress2stdout(sample_i,n_samples,10,'samples');
    end % sample_i
    climada_progress2stdout(0) % terminate
    fprintf('Plotting...\n');
    climada_EDS_DFC(EDS,EDS(1)); % plot risk curve
end % params.check_plot==2

end % climadacrop_init