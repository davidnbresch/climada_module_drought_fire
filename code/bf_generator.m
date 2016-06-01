function bf=bf_generator(number_years,store_intensity_field)
% BF event probabilistic
% MODULE:
%   drought_fire
% NAME:
%   bf_generator
% PURPOSE:
%   generate a hazard event set of bushfires for a specified number of
%   years and domain. It is currently set up for the region of Melbourne
%   with corresponding data. Next step: climada_bf_hazard_set.m
% CALLING SEQUENCE:
%   bf_generator(number_years);
% EXAMPLE:
%   bf_generator(5);
% INPUTS:
%   number_years: how many years of bushfire you want to simulate
%   store_intensity_field: if non-zero: matrix of size dx times dy, 
%   bf(i).intensity_field, is stored for all fires. for each fire, fire then 
%   can easily be plotted with for instance contourf(bf(1).intensity_field). 
%   if 0: matrix bf(i).intensity_field is not stored.
%   note: for many years, intensity_field must be set to 0, because not
%   enough memory capacity!
% OPTIONAL INPUT PARAMETERS:
%   none 
% OUTPUTS: 
%   a struct bf where bf(i) is one individual fire:
%   with bf(i).lon: array with longitudes of intensities
%   with bf(i).lat: array with latitudes of intensities
%   with bf(i).intensity: array with intensities at locations (bf(i).lon, bf(i).lat)
%   with bf(i).no_year: number of years simulated
%   with bf(i).intensity_field: matrix of size dx times dy which contains
%   intensities - note: only needed for plots. comment otherwise  needs a 
%   lot of memory (see below)
% MODIFICATION HISTORY:
% beuschl@student.ethz.ch, 20160601, initial, key author
% horatc@student.ethz.ch, 20160601, initial, key author
% david.bresch@gmail.com, 20160601, climada-compatibility
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables
if ~exist('bf_file','var'),bf_file=[];end

% DEFAULT
if ~exist('number_years','var'),number_years=0.01;end
if ~exist('store_intensity_field','var'),store_intensity_field=0;end
store_intensity_field = logical(store_intensity_field);

% DOMAIN
% lower left corner: (lat,lon) = (-38,143)
% upper right corner: (lat,lon) = (-34,147)

% domain area: 165'107 km^2
% current resolution: 1 cell == 0.1651 km^2 (with dx = dy = 1000)

domain_lat_ll = -38; % latitude value of lower left corner of domain
domain_lon_ll = 143; % longitude value of lower left corner of domain

domain_lat_ur = -34; % latitude value of upper right corner of domain
domain_lon_ur = 147; % longitude value of upper right corner of domain

dx = 1000; % number of cells in west-east direction
dy = 1000; % number of cells in south-north direction

scale_lat = dx/abs(domain_lat_ll - domain_lat_ur);
scale_lon = dy/abs(domain_lon_ll - domain_lon_ur);

% DATA
% Scaling data (1976-77 to 1995-96), Australian Bureau of Statistics
% (http://www.abs.gov.au/ausstats/abs@.nsf/0/ccb3f2e90ba779d3ca256dea00053977)

average_no = 424; % Average no. of fires: 424 fires/year in our domain 
average_si = 12; % Average fire destroys 12 cells

% CLIMATE CHANGE
% increase area burnt by factor area_screw due to climate change
area_screw = 1; % set to 1 for no climate change

% CALIBRATION
% total area burnt must be conserved - if cells_burnt is larger than 
% set_cells_burnt, computation is stopped

set_cells_burnt = average_no*average_si*area_screw*number_years;
cells_burnt = 0;

% tune with probability threshold to get right number of fires. you can
% also tune with number of time steps parameter time_steps inside for loop

p_threshold = 0.475; % if random number > p_threshold, fire is lighted

% init
n_o=1000000;
bf_lon=zeros(1,n_o);
bf_lat=zeros(1,n_o);
bf_intensity=zeros(1,n_o);
bf_no_year=zeros(1,n_o);

% START OF BUSHFIRE COMPUTATION
for o = 1:n_o
    if cells_burnt >= set_cells_burnt
        break;
    else
        
    % INITIALISATION
    % elements of tensor A can take on three different values: 
    % 0: unburnt cell - ready to burn
    % 1: fire         - certain probability that neighbouring cell is on fire 
    %                   as well in next time step, diagonal neighbours excluded
    % 2: ember (Glut) - still causing damage but does not ignite any
    %                   neighbouring cells

    % number of time steps is a random number, it can be played with it
    time_steps = max(1,round(rand*100)); 
    A = zeros(dx,dy,time_steps);

    % IGNITE FIRE
    % set spark somewhere randomly in the domain
    A(max(1,round(rand*dx)),max(1,round(rand*dy)),1) = 1; 

    % AREA CONSERVATION
    % set a maximum area in order to be able to better calibrate
    % area_max is taken for each fire from an exponential distribution with 
    % mean average_si, if it is reached, computation is stopped. ensures that
    % fires do not get too large

    area_max = max(1,round(exprnd(average_si)));
    count_cell = 0;

    for t = 1:(time_steps-1) % generation of one individual fire starts here
        if count_cell >= area_max
            break;
        else
            for i = 2:(dx-1)
                for j = 2:(dy-1)

                    % light new fire in the neighbourhood (if it is not yet
                    % burning or ember)
                    if A(i,j,t) == 1
                        if A(i-1,j,t) == 0 && rand > p_threshold;
                            A(i-1,j,t+1) = 1;
                        end
                        if A(i+1,j,t) == 0 && rand > p_threshold;
                            A(i+1,j,t+1) = 1;
                        end
                        if A(i,j-1,t) == 0 && rand > p_threshold;
                            A(i,j-1,t+1) = 1;
                        end
                        if A(i,j+1,t) == 0 && rand > p_threshold;
                            A(i,j+1,t+1) = 1;
                        end

                        A(i,j,t+1) = 2;
                    end

                    % ember remains ember - does not light new fire but is
                    % still increasing the damage
                    if A(i,j,t) == 2;
                            A(i,j,t+1) = 2;
                    end

                end
            end
            count_cell = nnz(A(:,:,t+1));
        end
    end

    % NUMERICAL INTENSITY
    % calculate intensity - sum up fields over all time steps
    B = zeros(dx,dy);
    for i = 1:time_steps
        B = B + A(:,:,i);
    end

    % PREPARE DATA FOR FURTHER USAGE BY climada_bf_hazard_set.m
    % create array of matrix B 
    [i,j,s] = find(B);

    % transform indizes into longitude and latitude, write it in struct array
    bf_lon(o) = domain_lon_ll + i./scale_lon;
    bf_lat(o) = domain_lat_ll + j./scale_lat;
    bf(o).intensity = s;
    bf(o).no_year = number_years;

    if store_intensity_field
        bf(o).intensity_field = B; 
    end
    
    cells_burnt = cells_burnt + nnz(B);
    %o % to print to stdout takes a lot of time
    end
end

if isempty(bf_file) % local GUI
    bf_file      = [climada_global.data_dir filesep 'hazards' filesep 'bf.mat'];
    hazard.filename = bf_file; 
    [filename, pathname] = uiputfile(bf_file, 'Save BF file as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        bf_file = fullfile(pathname,filename);
    end
    save(bf_file,'bf'); 
end

% store to output
bf.lon       = bf_lon;
bf.lat       = bf_lat;
bf.intensity = bf_intensity;
bf.no_year   = bf_no_year;

fprintf('Observed number of fires: ');
number_years*average_no
fprintf('Simulated number of fires: ');
length(bf)
fprintf('Observed number of cells burnt: ');
set_cells_burnt
fprintf('Simulated number of cells burnt: ');
cells_burnt

return
