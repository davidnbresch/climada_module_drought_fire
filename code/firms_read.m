function firms=firms_read(csv_file,check_plot)
% climada template
% MODULE:
%   module name
% NAME:
%   firms_read
% PURPOSE:
%   Read the MODIS5 firms csv data
%
%   Please proceed as follows:
%   1. Download the historic fire data (since Nov 2000-today) as .csv file
%      from https://firms.modaps.eosdis.nasa.gov/download (the full global
%      dataset since 2000 is too large for a single download) and store the
%      .csv file as ....
%      download the file as .csv and select MODIS C6
%   2. run firms_read to read the data
%
%   next call: bf_generator_large
%
%   Citation for both MCD14DL and MCD14ML: This data set was
%   provided by the University of Maryland and NASA FIRMS operated by
%   NASA/GSFC/ESDIS with funding provided by NASA/HQ
% CALLING SEQUENCE:
%   firms=firms_read(csv_file,check_plot)
% EXAMPLE:
%   firms=firms_read('TEST') % use the TEST dataset
% INPUTS:
%   csv_file: the name and path of the .csv file as downloaded from
%       https://firms.modaps.eosdis.nasa.gov/download
%       > promted for filename for if not given
%       for TESTS, there is a .csv file for Australia in the drought_fire
%       module's data folder, just set csv_file='TEST'
% OPTIONAL INPUT PARAMETERS:
%   param2: as an example
% OUTPUTS:
%   firms: the data, a structure with fields
%       filename: the original .csv filename with path
%           reads files such as firms.csv, fire_nrt_M6_*.cvs, fire_archive_M6_*.cvs
%       lat(i): latitudes of burning point i
%       lon(i): longitudes of burning point i
%       brightness(i): brightness of burning point i
%       datenum(i): the matlab datenum for each point i, also
%           datenum_unique: the unique datenums
%           unique_IA: the positions of the first unique elements
%           unique_IC: the output of [~,IA,IC]=unique
% MODIFICATION HISTORY:
% david.bresch@gmail.com, 20160703, initial
% david.bresch@gmail.com, 20160715, formts archive_M6 and nrt_M6 added
%-

firms=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('csv_file','var'),csv_file=[];end % OR:
if ~exist('check_plot','var'),check_plot=[];end

% locate the module's (or this code's) data folder (usually  afolder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% set default value for check_plot if not given
if isempty(check_plot),check_plot=0;end
%
file_version='firms'; % default
%
% TEST csv-file
TEST_csv_file=[module_data_dir filesep 'hazards' filesep 'external_model_output' filesep 'firms.csv'];

% template to prompt for filename if not given
if isempty(csv_file) % local GUI
    csv_file=[climada_global.data_dir filesep '*.csv'];
    [filename, pathname] = uigetfile(csv_file, 'Select firms .csv file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        csv_file=fullfile(pathname,filename);
    end
end

if strcmp(csv_file,'TEST')
    csv_file=TEST_csv_file;
    if exist(csv_file,'file')
        fprintf('TEST mode, using %s\n',TEST_csv_file);
    else
        fprintf('ERROR: %s not found, consider fetching from GitHub\n',TEST_csv_file);
        return
    end
end % strcmp(csv_file,'TEST')


[fP,fN]=fileparts(csv_file);
mat_file=[fP filesep fN '.mat'];

% figure file verison (as there seem to be diverse ones
if findstr(fN,'archive_M6'),file_version='archive_M6';end
if findstr(fN,'nrt_M6'),file_version='nrt_M6';end
if findstr(fN,'nrt_V1'),file_version='nrt_V1';end

if climada_check_matfile(csv_file,mat_file)
    load(mat_file);
else
    
    firms.filename=csv_file;
    
    % the csv file header lone reads:
    % geom,latitude,longitude,brightness,scan,track,acq_date,acq_time,satellite,confidence,version,bright_t31,frp
    % a first data line reads e.g.
    %010100000023DBF97E6A4C62405A643BDF4F0D43C0,-38.104,146.388,316.5,1,1,2006-01-24,1309,Terra,91,5.1       ,286.9,22.5
    
    fprintf('reading ...');
    % note: reads version as %s, since sometimes just '-', not a number
    %[geom,firms.latitude,firms.longitude,firms.brightness,scan,track,acq_date,acq_time,...
    
    switch file_version
        case 'firms'
            [~,firms.lat,firms.lon,firms.brightness,~,~,acq_date,acq_time,~,~,~,~,~] = ...
                textread(csv_file,'%s%f%f%f%f%f%s%f%s%f%s%f%f','delimiter',',','headerlines',1);
        case 'archive_M6' % no first field, additional column instrument
            [firms.lat,firms.lon,firms.brightness,~,~,acq_date,acq_time,~,~,~,~,~,~] = ...
                textread(csv_file,'%f%f%f%f%f%s%f%s%s%f%s%f%f','delimiter',',','headerlines',1);
        case 'nrt_M6' % no first field, additional columns instrument and daynight
            [firms.lat,firms.lon,firms.brightness,~,~,acq_date,acq_time,~,~,~,~,~,~,~] = ...
                textread(csv_file,'%f%f%f%f%f%s%f%s%s%f%s%f%f%s','delimiter',',','headerlines',1);
        case 'nrt_V1' % no first field, additional columns instrument and daynight, confidence a string
            [firms.lat,firms.lon,firms.brightness,~,~,acq_date,acq_time,~,~,~,~,~,~,~] = ...
                textread(csv_file,'%f%f%f%f%f%s%f%s%s%s%s%f%f%s','delimiter',',','headerlines',1);
    end
    
    % process acq_date and acq_time
    fprintf(' converting date/time ...');
    hours=fix(acq_time/100);
    minutes=fix(acq_time-hours*100);
    firms.datenum=datenum(acq_date)+hours/24+minutes/(24*60);
    % find unique date/time, i.e. the (daily) events
    %[firms.datenum_unique,firms.unique_IA,firms.unique_IC]=unique(firms.datenum);
    firms.datenum_unique=unique(firms.datenum);
    fprintf('done\n');
    
    fprintf('data period %s ..%s, %i events (%i records)\n',...
        datestr(min(firms.datenum)),datestr(max(firms.datenum)),...
        length(firms.datenum_unique),length(firms.brightness))
    
    if check_plot
        plot(firms.lon,firms.lat,'.r','MarkerSize',1)
        hold on;climada_plot_world_borders(1,'','',1);
    end % check_plot
    
    fprintf('saving as %s\n',mat_file);
    save(mat_file,'firms');
end % climada_check_matfile

end % firms_read