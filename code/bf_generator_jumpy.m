function bf=bf_generator_jumpy(number_years,event_i,store_intensity_field)
% creates probabilistic bushfire events based on data from FIRMS
% MODULE:
%   drought_fire
% NAME:
%   bf_generator_jumpy
% PURPOSE:
%   generate a hazard event set of bushfires for a specified number of
%   years and domain. Next step: climada_bf_hazard_set.m
% TEST:
%   version: use bf_TEST_simple to get started on how to evaluate
%   the output of this function
% CALLING SEQUENCE:
%   bf_generator_jumpy(years_number, event number, store_output_field);
% EXAMPLE:
%   bf=bf_generator_jumpy(1, 100, 1);
%   bf_generator_jumpy(1,100,1): % fast test
%   when asked for firms file to select: choose a firms csv file
%   that can be downloaded from the FIRMS website
%   (please see PURPOSE of firms_read.m)
% INPUTS:
%   number_years: how many years of bushfire you want to simulate
%   store_intensity_field: if non-zero: matrix of size dx times dy,
%   bf(i).intensity_field, is stored for all fires. for each fire, fire then
%   can easily be plotted with for instance contourf(bf(1).intensity_field).
%   if 0: matrix bf(i).intensity_field is not stored.
%   note: for many years, intensity_field must be set to 0, because not
%   enough memory capacity
% OPTIONAL INPUT PARAMETERS:
%   none
% OUTPUTS:
%   a struct bf where bf(i) is one individual fire:
%   with bf(i).lon: array with longitudes of intensities
%   with bf(i).lat: array with latitudes of intensities
%   with bf(i).intensity: array with intensities at locations (bf(i).lon, bf(i).lat)
%   with bf(i).years_simulated: number of years simulated
%   with bf(i).no_cells: number of burned cells
%   with bf(i).intensity_field: matrix of size dx times dy which contains
%   intensities - note: only needed for plots. comment otherwise  needs a
%   lot of memory (see below)
%   with bf(i).A
% MODIFICATION HISTORY:
% beuschl@student.ethz.ch, 20160601, initial, key author
% horatc@student.ethz.ch, 20160601, initial, key author
% david.bresch@gmail.com, 20160601, climada-compatibility
% Roxanne Doerge, Dina Haenseler, Elisabeth Tschumi, 20170609, initial, key author
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables
if ~exist('number_years','var'),number_years=[];end
if ~exist('store_intensity_field','var'),store_intensity_field=[];end
if ~exist('csv_file','var'),csv_file=[];end % OR:
% Parameters
%
if isempty(number_years),         number_years         = 1;end
if isempty(store_intensity_field),store_intensity_field= 0;end
store_intensity_field = logical(store_intensity_field);
%
% define the filename for the intermediate step data file
bf_file = [climada_global.data_dir filesep 'hazards' filesep 'external_model_output' filesep 'AUS_BF_proto_data_TEST.mat'];
%
average_no = 12800; % Average no. of fires: 3200 fires/year in our domain
average_si = 8; % Average fire destroys 32 cells: 8km^2

% CLIMATE CHANGE
% increase area burnt by factor area_screw due to climate change
area_screw = 1; % set to 1 for no climate change

% CALIBRATION
% total area burnt must be conserved - if cells_burnt is larger than
% set_cells_burnt, computation is stopped

set_cells_burnt = average_no*average_si*area_screw*number_years;
cells_burnt = 0;
wiggle=0.05; % max. of cells by which onset of fire can be wiggled, 5km

% tune with probability threshold to get right number of fires. you can
% also tune with number of time steps parameter time_steps inside for loop

p_threshold = 0.475; % if random number > p_threshold, fire is lighted

n_o=event_i;   % simulated no of fires

% prompt for csv_file if not given
if ~exist('csv_file','var'),csv_file=[];end % OR:

if isempty(csv_file) % local GUI
    csv_file=[climada_global.modules_dir filesep 'climada_module_drought_fire-master/data/hazards/external_model_output' filesep '*.csv'];
    [filename, pathname] = uigetfile(csv_file, 'Select firms .csv file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        csv_file=fullfile(pathname,filename);
    end
end

firms = firms_read(csv_file);


% START OF BUSHFIRE COMPUTATION
dx = 100; % number of cells in west-east direction
dy = 100; % number of cells in south-north direction

for o = 1:n_o
    if cells_burnt >= set_cells_burnt
        disp('Break')
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
%         time_steps = 100;
        time_steps = 50;%max(1,round(rand*100));
        A = zeros(dx,dy,time_steps);
        
        % IGNITE FIRE
        % set spark somewhere randomly in the domain
        % Start new fire at position

        which=ceil(rand*length(firms.lat));
        fire_lat = firms.lat(which)+2*(rand-0.5)*wiggle;
        fire_lon = firms.lon(which)+2*(rand-0.5)*wiggle;
        clear which
        
        % square area of generator: 50x50km
        % 10000 cells, 500x500m
        lat.min = fire_lat-0.2577; %+25km
        lat.max = fire_lat+0.2577;
        lat.diff = lat.max-lat.min;
        lon.min = fire_lon-0.225;
        lon.max = fire_lon+0.225;
        lon.diff = lon.max-lon.min;

        scale_lat = dx/(lat.max - lat.min);
        scale_lon = dy/(lon.max - lon.min);

        ny = ceil((fire_lat-lat.min)*scale_lat);%+ceil(2*(rand-0.5)*wiggle);
        nx = ceil((fire_lon-lon.min)*scale_lon);%+ceil(2*(rand-0.5)*wiggle);

        A(nx,ny,1) = 1;
        
        % AREA CONSERVATION
        % set a maximum area in order to be able to better calibrate
        % area_max is taken for each fire from an exponential distribution with
        % mean average_si, if it is reached, computation is stopped.
        % ensures that fires do not get too large
        
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
            %C(:,:,t)=sum(A,3);
        end
        
        
        % NUMERICAL INTENSITY
        % calculate intensity - sum up fields over all time steps
        
         B=sum(A,3) ; % sum over 3rd dim of A
        
        % PREPARE DATA FOR FURTHER USAGE BY climada_bf_hazard_set.m
        % create array of matrix B
        [i,j,s] = find(B);
        
        % transform indizes into longitude and latitude, write it in struct array
        bf(o).lon = lon.min + i./scale_lon;
        bf(o).lat = lat.min + j./scale_lat;
        bf(o).intensity = s;
        bf(o).years_simulated = number_years; %no_year
        bf(o).no_cells = nnz(B);   
        
             if store_intensity_field
                bf(o).intensity_field = B;
             end
            
        bf(o).A=A;
        
        cells_burnt = cells_burnt + nnz(B);
    end
end

fprintf('saving bushfire proto-hazard set (bf) as %s\n',bf_file);
save(bf_file,'bf',climada_global.save_file_version) % for HDF5 format (portability)

fprintf('Observed number of fires : %i\n',number_years*average_no);
fprintf('Simulated number of fires: %i\n',length(bf));
fprintf('Observed number of cells burnt : %i\n',set_cells_burnt);
fprintf('Simulated number of cells burnt: %i\n',cells_burnt);

end % bf_generator