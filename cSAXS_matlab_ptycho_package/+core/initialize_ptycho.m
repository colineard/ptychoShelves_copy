%INITIALIZE_PTYCHO 
% everything that needs to be done before triggering the reconstruction. This includes 
% inter alia initial checks, intial guess preparations and loading the data. 
%
% ** p          p structure
%
% returns:
% ++ p          p structure
% ++ status     status flag
%
%
% see also: core.ptycho_recons
%
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



function [ p, status ] = initialize_ptycho( p )

import utils.*
import io.*

%%% read meta data %%%
if ~isfield(p, 'src_metadata')
    verbose(0,' p.src_metadata is not set, using default p.src_metadata = ''spec''')
    p.   src_metadata = 'spec';    % load meta data from file; currently only 'spec' is supported;
end

% check store_images flag
if ~isfield(p.save, 'store_images')
    p.save.store_images = true;
end
if ~p.save.store_images
    close all
end

% prepare container for meta data
assert( isnumeric(p.scan_number), 'p.scan_number has to contain an integer number')
p.numscans = length(p.scan_number); % Number of scans
p.meta = cell(1,length(p.scan_number));

p = scans.read_metadata(p);

for ii = 1:p.numscans
    p.   scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
end

% write procID
write_procID(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checks and defaults %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = core.initial_checks(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% prepare paths, note that it was already initialized in ptycho_recons %%%

p = core.ptycho_prepare_paths(p);


%%% load detector settings %%%
p = detector.load_detector(p);


if isfield(p, 'ds') && ~isempty(p.ds)
    warning(['Defining ds in the template is not supported anymore and ' ...
    'will not change the pixel size. Please make sure '...
    'that it is set correctly in +detector/+%s/%s.m and remove ds from your template.'], p.detector, p.detector)
end
for ii=1:length(p.detectors)
    assert(p.detectors(1).params.pixel_size==p.detectors(ii).params.pixel_size, 'Different detector pixel sizes are not supported at the moment.')
end
p.ds = p.detectors(1).params.pixel_size;


if check_option(p, 'prop_regime', 'nearfield')
    % nearfield ptychography 
    assert(check_option(p,'focus_to_sample_distance'), 'Undefined p.focus_to_sample_distance that is required for nearfield ptychography')
    p.nearfield_magnification = (p.z-p.focus_to_sample_distance)/p.focus_to_sample_distance; 
    verbose(1, 'Propagation in nearfield regime, magnification = %g', p.nearfield_magnification)
    p.dx_spec = [p.ds,p.ds] / p.nearfield_magnification;
    p.z = p.z / p.nearfield_magnification;
else
    % standard farfield ptychography
    p.dx_spec = p.lambda*p.z ./ (p.asize*p.ds);                   % resolution in the specimen plane
    p.dx_spec = p.dx_spec ./ cosd(p.sample_rotation_angles(1:2));      % account for a tilted sample ptychography      
end


%%% prepare positions %%%
p = scans.read_positions(p);

%%% find which positions belongs to each object 
p = core.find_shared_IDs(p); 


%%% prepare positions
% Prepare positions, note the output is already in probe positions which
% are different from object (scan) positions by a minus sign
p = core.ptycho_adjust_positions(p);
% p.positions_orig = p.positions;
% p.numpts_orig = p.numpts;
p.numpos = sum(p.numpts);


p.asize_nobin = p.asize;



%%%%%%%%%%%%%%%%%%%%%%
%%% prepare scans %%%%
%%%%%%%%%%%%%%%%%%%%%%


% make sure that all scans have a ctr
numctr = size(p.ctr,1);
if p.numscans > numctr
    for ii=numctr+1:p.numscans
        p.ctr(end+1,:) = p.ctr(numctr,:);
    end
elseif p.numscans < numctr
    p.ctr(p.numscans+1:end,:) = [];
end

%%%%  load data, mask and generate initial estimate of the probe
if p.prepare.auto_prepare_data
    [p, status]=core.ptycho_prepare_scans(p);
else
    if ~isa(p.prepare.prepare_data_function, 'function_handle')
        error(['Expected function handle as p.prepare.prepare_data_function. '... 
            'Please update p.prepare.auto_prepare_data or set p.prepare.auto_prepare_data=true.'])
    else
        [p, status] = p.prepare.prepare_data_function(p);
    end
end
if ~status
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot prepared data %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.plot.prepared_data && p.use_display
    core.analysis.plot_raw_data(p)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot initial guess %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define combined strings for figure title
p.plot.obtitlestring = '';
p.plot.prtitlestring = '';
p.plot.errtitlestring = '';
if p.share_object
    p.plot.obtitlestring =  [core.generate_scan_name(p) ' '];
end
if p.share_probe
    p.plot.prtitlestring = [core.generate_scan_name(p) ' '];
end
p.plot.errtitlestring = [core.generate_scan_name(p) ' '];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot initial guess %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.use_display
    p.plot.extratitlestring = sprintf(' (%dx%d) - Initial guess', p.asize(2), p.asize(1));
    core.analysis.plot_results(p, 'use_display', p.use_display);
end
p.plot.extratitlestring = sprintf(' (%dx%d)', p.asize(2), p.asize(1));

if ~isfield(p.plot, 'windowautopos')
    p.plot.windowautopos = false; % So resizing after first time display is respected
end

verbose(1, 'Finished data preparation and initialization.')

end




