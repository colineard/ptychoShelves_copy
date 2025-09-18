%%  TEST TEMPLATE FOR FUNTIONALITY OF RECONSTRUCTION ON REAL DATA 
% 1) call standard template to get fresh settings defaults 
% 2) load example of measured data 
% 3) all c_solver engine to test quality 


%% set shared parameters for all test scripts 
run(fullfile( fileparts(mfilename('fullpath')), 'init_test.m'))

%% General

p.   z = 7.2045;                                            % Distance from object to detector 
p.   src_metadata = 'nexusCXS'; 

% Scan queue
p.   scan_number = [35 36];                                    % Multiple scan numbers for shared scans


% Data preparation
p.   asize = [1000 1000];                                     % Diffr. patt. array size   
p.   ctr = [509 845];                                       % Diffr. patt. center coordinates (y,x) (empty means middle of the array); e.g. [100 207;100+20 207+10];
p.   detector.check_2_detpos = [];                                   % = []; (ignores)   = 270; compares to dettrx to see if p.ctr should be reversed (for OMNY shared scans 1221122), make equal to the middle point of dettrx between the 2 detector positions
% p.   prepare.data_preparator = 'python';                                  % data preparator; 'python' or 'matlab' 
p.   src_metadata = 'nexusCXS';                    % load meta data from file; currently only 'spec' is supported;

% Scan positions
p.   src_positions = 'nexusCXS';                           % 'spec', 'orchestra', 'load_from_file', 'matlab_pos' (scan params are defined below)
p.   angular_correction_setup = 'omny';                         % if src_positions=='orchestra', choose angular correction for specific cSAXS experiment: 'flomni', 'omny', 'lamni', 'none', 
p.   orchestra.positions_file = [];%['/das/work/p16/p16812/data/nexus/specES1/scan_positions/scan_%05d.dat'];    % Filename pattern for position files, Example: ['../../specES1/scan_positions/scan_%05d.dat']; (the scan number will be automatically filled in)
p.   detector.name = 'Eiger1p5m'; 
% p.   detector.oldEiger1p5m = true;

p.   affine_matrix = [1 0;tan(0.46*pi/180) 1];                                    % Applies affine transformation (e.g. rotation, stretching) to the positions (ignore by = []). Convention [yn;xn] = M*[y;x]. For flOMNI we found in September 2018: = [1 0;tan(0.36*pi/180) 1]; for OMNY we found in October 2018: = [1 0;tan(0.4*pi/180) 1]; laMNI in June 2018  [1,0.0154;-0.0017,1.01]; laMNI in August [1.01 0.0031; -0.0018 1.00] 

p.   prepare.force_preparation_data = true; 

% I/O
p.   base_path = fullfile(ptycho_path, 'tests');
p.   specfile = [];%'/das/work/p16/p16812/data/eiger1p5_h5/';                                         % Name of spec file to get motor positions and check end of scan, defaut is p.spec_file == p.base_path;
p.   detector.name = 'eiger1p5M';                                  % 'pilatus' or 'eiger'
p.   raw_data_path{1} = '/das/work/p16/p16812/data/nexus/';                                 % Default using compile_x12sa_filename, used only if data should be prepared automatically

%% Reconstruction

% Initial iterate object
p.   model_object = true;                                   % Use model object
p.   model_object_type = 'rand';                            % specify how the object shall be created; use 'rand' for a random initial guess; use 'amplitude' for an initial guess based on the prepared data

p.   initial_iterate_object_file{1} = '';                   %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'

% Initial iterate probe
p.   model_probe = false;                                    % Use model probe
p.   initial_probe_file = '/das/work/units/csaxs/p16812/data/nexus/analysis/S00000-00999/S00012/S00012_512x512_b0_run_1_recons.h5';% Use probe from this mat-file (not used if model_probe is true)

% Shared scans - Currently working only for sharing probe and object
p.   share_probe = 0;                                       % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
p.   share_object = 0;                                      % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 

% Modes
p.   probe_modes  = 1;                                      % Number of coherent modes for probe
p.   object_modes = 1;                                      % Number of coherent modes for object
% Mode starting guess
p.   mode_start_pow = [0.02];                               % Normalized intensity on probe modes > 1. Can be a number (all higher modes equal) or a vector
p.   mode_start = 'herm';                                   % (for probe) = 'rand', = 'herm' (Hermitian-like base), = 'hermver' (vertical modes only), = 'hermhor' (horizontal modes only)
p.   ortho_probes = true;                                   % orthogonalize probes after each engine

%% Plot and save
p.   plot.prepared_data = false;                         % plot prepared data
p.   plot.calc_FSC = true;                                      % Calculate the Fourier Shell correlation for 2 scans
p.   plot.show_FSC = false;                                % Show the FSC plots, including the cropped FOV
p.   save.store_images = 0;                                       % Write nice jpegs in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
p.   plot.plot_interval = inf;                                       % Write nice jpegs in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho



%% ENGINES
% External C++ code
if 0
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    eng.  name = 'presolver';                   % A wrapper around c_solver engine which uses DM to presolve reconstruction first in low resolution. It can be used to refine the lowest spatial frequencies in a computationally cheaper way
    eng.  asize_presolve = [350 350];           % Diffraction size (p.asize value) for the downsampled dataset, IMPORTANT: this is asize before binning
    eng.  number_iterations = 500;              % Total number of iterations of DM method 
    eng.  probe_regularization = .1;            % DM: Weight factor for the probe update; 
    eng.  probe_change_start = 1;               % DM: Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % DM:   Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .05;                % DM: Relaxation in the Fourier domain projection, = 0  for full projection    
    
    % advanced setting unrelated to reconstruction quality 
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = [];                  % cSAXS only option: beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  use_gpu = true;                      % use GPU if available
    eng.  gpu_id = 3;%gpuid;%utils.find_unused_gpu_card(4); %1;                           % GPU device number; can use   utils.find_unused_gpu_card(4);   to automatically check; or check with nvidia-smi (+1 for Matlab indexing)
    
    eng.  ra_nodes = 2;                         % PSI only option:   number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file(19.83, 19.83) = '';    		% Name of the initial conditions file, default if empty. 
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. 
    eng.  solution_file = '';                   % Name of the solution file, default if empty. 
    eng.  force_prepare_h5_files = false;       % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
end
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    eng.  name = 'c_solver';
    eng.  number_iterations = 150;              % Total number of iterations
    eng.  opt_iter = 300;                       % Iterations for optimization     
    eng.  probe_regularization = .1;            % Weigth factor for the probe update; 
    eng.  probe_change_start = 1;               % Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .08;                % Relaxation in the Fourier domain projection, = 0  for full projection    
    eng.  background = 0;                       % [PARTIALLY IMPLEMENTED (not fully optimized)] Add background to the ML model in form:  |Psi|^2+B, B is in average counts per frame and pixel
    eng.  probe_support_fft = false;            % [PARTIALLY IMPLEMENTED (not fully optimized)] Apply probe support in Fourier space, ! uses model zoneplate settings to estimate support size 

    eng.  N_layer = 1;                          % Number of virtual object layers (slices)
    eng.  delta_z = 0e-6 * ones(1, eng.N_layer-1); % Separation between object slices 
    %eng.  ms_init_ob_fraction = [1 0];
    if eng.  N_layer>1
        p.sufix = [p.sufix '_N' num2str(eng. N_layer)];
        eng.  number_iterations = 0; % highly recommended
    end
    
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = {};                  % beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  use_gpu = true;                       % use GPU if available
    eng.  gpu_id = 3;
    eng.  num_gpus = 1;                         % number of GPUs
    eng.  ra_nodes = 2;                         % number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. Do not use ~ in the path
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. Do not use ~ in the path
    eng.  solution_file = '';                   % Name of the solution file, default if empty. Do not use ~ in the path
    eng.  force_prepare_h5_files = 0;           % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
   

if gpuDeviceCount
    
    eng.  name = 'c_solver';
    eng.  number_iterations = 100;              % Total number of iterations
    eng.  opt_iter = 300;                       % Iterations for optimization     
    eng.  probe_regularization = .1;            % Weigth factor for the probe update; 
    eng.  probe_change_start = 1;               % Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .08;                % Relaxation in the Fourier domain projection, = 0  for full projection    
    eng.  background = 0;                       % [PARTIALLY IMPLEMENTED (not fully optimized)] Add background to the ML model in form:  |Psi|^2+B, B is in average counts per frame and pixel
    eng.  probe_support_fft = false;            % [PARTIALLY IMPLEMENTED (not fully optimized)] Apply probe support in Fourier space, ! uses model zoneplate settings to estimate support size 

    eng.  N_layer = 1;                          % Number of virtual object layers (slices)
    eng.  delta_z = 0e-6 * ones(1, eng.N_layer-1); % Separation between object slices 
    %eng.  ms_init_ob_fraction = [1 0];
    if eng.  N_layer>1
        p.sufix = [p.sufix '_N' num2str(eng. N_layer)];
        eng.  number_iterations = 0; % highly recommended
    end
    
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = [];                  % beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  use_gpu = true;                       % use GPU if available
    eng.  gpu_id = 3;                           % GPU device number; check with nvidia-smi (+1 for Matlab indexing)
    eng.  num_gpus = 1;                         % number of GPUs
    eng.  ra_nodes = 2;                         % number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. Do not use ~ in the path
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. Do not use ~ in the path
    eng.  solution_file = '';                   % Name of the solution file, default if empty. Do not use ~ in the path
    eng.  force_prepare_h5_files = 0;           % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
    
% % --------- GPU engines  -------------   See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
%     eng = struct();                 % reset settings for this engine 
%     eng. name = 'GPU';    
%     eng. gpu_id = [];               % default GPU id, [] means choosen by matlab
%     eng. probe_modes = 1; 
%     eng. probe_support_radius = 0.9;      % Normalized radius of circular support, = 1 for radius touching the window    
%     eng. probe_support_fft = false;        % assume that there is not illumination intensity out of the central FZP cone 
%     eng. probe_support_fft_type = 'circ';   % choose the type of Fourier support; select 'circ' for a circular support based on FZP parameters, 'gaps' for detector gaps and 'mask' for the detector mask
%     eng. probe_support_fft_float = false;   % apply only a weak probe support instead of setting the values to 0
%     eng. probe_support_fft_relax = 0.9;     % probe support relaxation, 1 for full projection
% 
%     % basic recontruction parameters 
%     % PIE / ML methods              % See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
%     eng. beta_object = 1;           % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
%     eng. beta_probe = 1;            % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
%     eng. delta_p = 0.1;             % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, 
%     eng. momentum = 0.5;              % add momentum term to the MLc method, eng.momentum = multiplication gain for velocity
%     eng. accelerated_gradients_start = 2;   % iteration number from which the Nesterov gradient acceleration should be applied, this option is supposted only for MLc method 
% 
%     % DM
%     eng. pfft_relaxation = 0.05;     % Relaxation in the Fourier domain projection, = 0  for full projection 
%     eng. probe_regularization = 0.1; % Weight factor for the probe update (inertia)
% 
%     % other extensions 
%     eng. background = 0.001;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
% 
%     eng. method = 'DM';            % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
%     eng. number_iterations = 100;   % number of iterations for selected method 
% 
%     [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
% 
%     eng. method = 'MLc';            % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
%     eng. number_iterations = 500;   % number of iterations for selected method 
% 
%     [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end
    

%% Run the reconstruction
% matlab data prep
p.prepare.data_preparator = 'matlab';
run_recons_test(p, 'matlab_ps', 1);

% % matlab data prep WITH BINNING
% p.prepare.data_preparator = 'matlab';
% p.detector.binning = true; 
% run_recons_test(p, 'matlab_ps', 1);

% if gpuDeviceCount
%     % matlab data prep
%     p.detector.binning = false; 
%     p.prepare.data_preparator = 'matlab';
%     run_recons_test(p, 'matlab_ps', 2:3);
% end




function run_recons_test(p, arg, engine_ids )

% run only preselected engines 
p.engines = p.engines(engine_ids); 
% reconstruct 
out = core.ptycho_recons(p);

if p.detector.binning
    fprintf('Testing BINNED dataset "EIGER1p5M - OMNY - %s  engine %s"  .... resolution %3.2f\n', arg, p.engines{1}.name,  out.FSC.resolution(end))
    return
end
    

cmp = load('/das/work/p16/p16812/data/nexus/tests/reference.mat');
if all(cmp.ref.resolution-out.FSC.resolution>=1)
    fprintf('Testing real dataset "NeXus EIGER1p5M - OMNY - %s  engine %s"  .... OK\n', arg, p.engines{1}.name)
    fprintf('The resolution improved from [%3.2f %3.2f] to [%3.2f %3.2f]. Please consider updating the reference!\n', cmp.ref.resolution, out.FSC.resolution)
elseif all(cmp.ref.resolution-out.FSC.resolution<=-1)
    fprintf('Testing real dataset "NeXus EIGER1p5M - OMNY - %s  engine %s"  .... failed\n', arg, p.engines{1}.name)
    warning('The resolution dropped from %3.2f to %3.2f!\n', cmp.ref.resolution(1), out.FSC.resolution(1))
else
    fprintf('Testing real dataset "NeXus EIGER1p5M - OMNY - %s  engine %s"  .... OK\n', arg, p.engines{1}.name)
    fprintf('The resolution is the same as saved as reference [%3.2f %3.2f] nm.\n', cmp.ref.resolution)
end

% delete temporal data
for path = out.save_path
    rmdir(path{1}, 's')
end
end


%end
% 2011-11-24
% Parameter to autoposition windows on first display - p.windowautopos
% Replaced powerbound with countbound. countbound represents the mean
  % number of photons in a change below which no projection is taken. It
  % scales automatically with exposure time (number of photons in
  % measurement)
% Real axes option to show plots in microns
% Read parameters from spec
% Implement user suplied object_flat_region
% Implemented option for reconstructing when having 2 repeated scans in the prepared data file

% 2011-11-29
% Template seemed extracted from an AFS run, I modified directories for
% direct use on ../../
% Implemented test mode
% Added cutoff value at beginning
% Added auto settings for prepare data, scan numbers
% Implemented reading from spec. Note it will use the values from the first
% scan
% Added option for repeated scan, should be enabled for 2 detector positions

% 2012-08-23
% Replaced default prepare data function to prepare_data_2d
% In I/O section: added option for a sufix 
% Added default option for raw data path based on compile_x12sa_filename
% Added options to autoprepare data, with cutoff and burstmode detected if
  % the file does not exist. Also added the possiblity to override and
  % force a repreparation of data
% Added a data prefix option (for eaccount_1_) and defaults using
  % identify_eaccount

% 2012-10-29
% Added option for binning and some checks for OMNY detector position scans

% 2012-10-31
% Added options to use the external C-code for testing

% 2015-05-13
% Added option to queue file tasks from OMNI.
% For this I moved the default checks and generation of default names and
% paths to ptycho_recons. Que Dios se apiade de nosotros.

% 2016-02-11
% Removed old option for dump files
% Added p.store_images, if this flag is on and p.use_display it will open 
% figures in the background and write nice jpegs of the reconstruction and error metric anyway

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the PtychoShelves 
%   computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the PtychoShelves package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite 
% K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
% high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
% and for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for LSQ-ML:
% M. Odstrčil, A. Menzel, and M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express 26(3), 3108 (2018). 
%   (doi: 10.1364/OE.26.003108),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089),
% and/or for OPRP:
% M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation. 
% Opt. Express 24.8 (8360-8369) 2016. (doi: 10.1364/OE.24.008360).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland.

