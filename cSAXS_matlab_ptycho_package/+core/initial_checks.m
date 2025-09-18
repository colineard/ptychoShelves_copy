%INITIAL_CHECKS set default values for the most common variables
%
% ** p          p structure
%
% returns:
% ++ p          p structure
%
% see also: core.initialize_ptycho
%

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

function [p] = initial_checks(p)


import utils.*
import io.*

%%%%%%%%%%%%%
%% General %%
%%%%%%%%%%%%%

% check matlab version
check_matlab_version(9.3);

if ~isfield(p, 'use_display') || isempty(p.use_display)
    if verbose > 1
        p.use_display = true;
    else
        p.use_display = false;
    end
end

if ~usejava('desktop')
   % test if matlab was called with -nodisplay option, if yes then
   % use_display should be false
   p.use_display = false;
end

% check if fourier ptycho recon is needed
if ~isfield(p, 'fourier_ptycho')
    p.fourier_ptycho = (isfield(p, 'FP_focal_distance') && ~isempty(p.FP_focal_distance));
end

if ~isfield(p, 'sample_rotation_angles')
    p.sample_rotation_angles = [0,0,0];   % 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)  
end


%%% Derived quantities %%%
assert(~isempty(p.energy), 'Provide p.energy or source of metadata p.src_metadata')
p.lambda = 1.2398e-9/p.energy;                                % wavelength
if isscalar(p.asize); p.asize = [p.asize p.asize]; end

%%%%%%%%%%%%%%%%%%%%
%% Scan meta data %%
%%%%%%%%%%%%%%%%%%%%

% calculate fourier ptycho geometry
if p.fourier_ptycho
    if ~isfield(p, 'FP_focal_distance')
        error('For running Fourier Ptychography, please specify the focal length of your objective lens (p.FP_focal_distance)');
    end
    if ~get_option(p, 'z_lens')         
        p.z_lens = 1/(1/(p.FP_focal_distance)-1/(p.z));     
    end
end

%%%%%%%%%%%%%%%%
%% Scan queue %%
%%%%%%%%%%%%%%%%

% number of attempts to reconstruct the given dataset
if ~isfield(p.queue, 'max_attempts')
    p.queue.max_attempts = 5;
end

% lock files
if ~isfield(p.queue, 'lockfile') || isempty(p.queue.lockfile)
    p.queue.lockfile = false;
end

if ~isfield(p.queue, 'file_queue_timeout')
    p.queue.file_queue_timeout = 10; % time to wait for a new dataset in queue
end

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%


% suffix for prepared data file
if ~isfield(p.prepare, 'prep_data_suffix')
    p.prepare.prep_data_suffix = '';
end

if p.asize(1) ~= p.asize(2) && p.prepare.force_preparation_data == false
    verbose(1, 'Loading prepared data is not supported for asymmetric p.asize, enforce load from raw data ')
    p.prepare.force_preparation_data = true; 
end


% data preparator
if ~isfield(p.prepare, 'data_preparator')
    p.prepare.data_preparator = 'matlab_ps';
end
if any(strcmpi(p.prepare.data_preparator, {'python', 'libDetXR','json'}))
    error('libDetXR / python data preparator is deprecated.')
end

if strcmpi(p.prepare.data_preparator, 'matlab')
    p.prepare.data_preparator = 'matlab_ps';
end

% binning 
if ~isfield(p.detector,'binning')
    p.detector.binning = false; 
end

% upsampling 
if ~isfield(p.detector,'upsampling')
    p.detector.upsampling = false; 
end

% interpolate 
if ~isfield(p.detector,'interpolate')
    p.detector.interpolate = false; 
end




% prealignment for Fourier Ptychography
if check_option(p, 'FP_focal_distance')
    p.   fourier_ptycho = true;                                % set to true for Fourier Ptychography
else
    p.   fourier_ptycho = false;           
end

if ~isfield(p, 'prealign_FP')
    p.prealign_FP = false;
end

% Fourier Ptycho is only supported by Matlab data preparation
if p.fourier_ptycho
    % set prealign_data to true if not distortion correction is available
    if p.prealign_FP && ~p.prealign.prealign_data && isempty(p.prealign.distortion_corr)
        p.prealign.prealign_data = true;
    end
end

% set defaults for matlab_ps
if strcmpi(p.prepare.data_preparator, 'matlab_ps')
    if ~isfield(p.io, 'data_precision')
        p.io.data_precision = 'single';
    end
    if ~isfield(p.io, 'data_nthreads')
        p.io.data_nthreads = 2;
    end
    if ~isfield(p.prepare, 'pad_data')
        p.prepare.pad_data = false;
    end
end


if isfield(p, 'prop_regime') && ~ismember(p.prop_regime, {'nearfield', 'farfield'})
    error(['Nonexistent propagation regime ', p.prop_regime ])
end    

% store prepared data
if ~isfield(p.prepare, 'store_prepared_data')
    p.prepare.store_prepared_data = true; 
end

%%%%%%%%%%%%%%%%%%%%
%% Scan positions %%
%%%%%%%%%%%%%%%%%%%%

% load positions from prepared file
if ~isfield(p.io, 'load_prep_pos')
    p.io.load_prep_pos = false;
end

%%%%%%%%%
%% I/O %%
%%%%%%%%%

% file compression
if ~isfield(p.io, 'file_compression')
    p.io.file_compression = 0;
end
if ~isfield(p.io, 'data_compression')
    p.io.data_compression = 3;
end


% run name
if ~check_option(p, 'run_name')
    % check if prefix is defined 
    if isempty(p.prefix) 
        if iscell(p.scan_str)
            p.prefix =  p.scan_str{1};
        else
            p.prefix = p.scan_str;
        end
    end
    p.run_name = sprintf('%s_%s', p.prefix, datestr(now, 'yyyy_mm_dd'));
end
verbose(3, 'run_name = %s', p.run_name);

%%%%%%%%%%%%%%%%%%%%
%% Reconstruction %%
%%%%%%%%%%%%%%%%%%%%

% backward compatibilty for initial_iterate
if isfield(p, 'initial_iterate') && ~isfield(p, 'initial_iterate_object')
    p.initial_iterate_object = p.initial_iterate;
    p = rmfield(p, 'initial_iterate');
end
if isfield(p, 'initial_iterate_file') && ~isfield(p, 'initial_iterate_object_file')
    p.initial_iterate_object_file = p.initial_iterate_file;
    p = rmfield(p, 'initial_iterate_file');
end

% model probe
if ~isfield(p.model, 'probe_central_stop')
    p.model.probe_central_stop = false;
end

if ~isfield(p.model, 'probe_central_stop_diameter') && p.model.probe_central_stop
    p.model.probe_central_stop_diameter = 50e-6;
end

%%%%%%%%%%%%%%%%%%%
%% Plot and save %%
%%%%%%%%%%%%%%%%%%%

% plot prepared data
if ~isfield(p.plot, 'prepared_data') || (isfield(p.plot, 'prepared_data')&& isempty(p.plot.prepared_data))
    if p.verbose_level > 2
        p.plot.prepared_data = true;
    else
        p.plot.prepared_data = false;
    end
end

% plotting
if ~isfield(p.plot, 'interval') || isempty(p.plot.interval)
    if verbose > 2
        p.plot.interval = 10;
    else
        p.plot.interval = 200;
    end
end

% external call to save figures
if ~isfield(p.save, 'external')
    p.save.external = false;
end

% propagation and apodization
if ~isfield(p.plot, 'obj_apod')
    p.plot.obj_apod = false;
end
if ~isfield(p.plot, 'prop_obj')
    p.plot.prop_obj = 0;
end

% calculate FSC
if ~isfield(p.plot, 'calc_FSC')
    p.plot.calc_FSC = false;
end
if ~isfield(p.plot, 'show_FSC')
    p.plot.show_FSC = utils.verbose>2;
end
if ~isfield(p.plot, 'probe_spectrum')|| isempty(p.plot.probe_spectrum)
    p.plot.probe_spectrum = utils.verbose>2;
end
if ~isfield(p.plot, 'object_spectrum')|| isempty(p.plot.object_spectrum)
    p.plot.object_spectrum = utils.verbose>2;
end
if ~isfield(p.save, 'store_images_ids' )|| isempty(p.save.store_images_ids)
    p.save.store_images_ids = 1:4; 
end


%%%%%%%%%%%%%
%% Engines %%
%%%%%%%%%%%%%

% at least one engine has to be specified
if ~isfield(p, 'engines')
    error('At least one reconstruction engine has to be selected. Please check your template.')
end



% first engine is external
p.external_engine0 = strcmpi(p.engines{1}.name, 'c_solver') || ...
    (isfield(p.engines{1}, 'external') && p.engines{1}.external);

if ~isfield(p,'remove_scaling_ambiguity')
    % if true. try to keep norm(probe) constant during the reconstruction
    p.remove_scaling_ambiguity = true; 
end


end

