%EIGER1P5M settings for Eiger 1.5M detector
%   This file is called if you specify p.detector = 'eiger1p5M'
%
% ** p          p structure
% returns:
% ++ det        detector structure; later accessible via p.detectors(ii).params
%
% see also: detector.load_detector


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

function [ det ] = eiger1p5M( p )


%% parameter
det.pixel_size = 75e-6;                 % detector pixel size
det.file_extension = 'h5';             % raw data file extension 
det.mask = fullfile(p.base_package_path, 'detector_masks', 'eiger_4_valid_mask.h5');    % detector mask
det.maskPath = '/entry/collection/mask';% H5 location for HDF5 mask files
if exist(fullfile(p.raw_data_path{1}, 'data'), 'dir')
    det.basepath_dir = 'data/';           % raw data directory
elseif exist(fullfile(p.raw_data_path{1}, 'eiger_4'), 'dir')
    det.basepath_dir = 'eiger_4/';           % raw data directory
elseif exist(fullfile(p.raw_data_path{1}, 'eigeromny'), 'dir')
    det.basepath_dir = 'eigeromny/';           % raw data directory, legacy, remove in future 
end
det.mask_saturated_value = [];          % mask pixels above given value
det.mask_below_value = [];              % mask pixels below given value
det.data_stored = true;                 % false == data are generated "onfly", no need to load / store (simulations)
det.detposmotor = 'dettr';              % name of spec motor used for checking detector position
det.readout = [];                       % detector readout, e.g. {[50:500], [50:500]}; leave empty for full readout


% det.orientation = [0 0 1];              % [<Transpose> <FlipLR> <FlipUD>]
% changed from 30/08/2018
det.orientation = [1 0 0];              % [<Transpose> <FlipLR> <FlipUD>]

% old Eiger1.5M data structure?
det.oldEiger1p5m = false;               
if isfield(p.detector, 'oldEiger1p5m') && p.detector.oldEiger1p5m
    det.oldEiger1p5m = true;
end

% H5 data path
if strcmpi(det.file_extension, 'h5')
    if det.oldEiger1p5m
        det.H5_detector_path = '/entry/data/eiger_4/';
    else
        det.H5_detector_path = '/entry/instrument/eiger_4/data/';
    end
end

%% define directory tree structure
% define a function handle for compiling the scan directory tree. Input
% should be the scan number.
det.read_path = @utils.compile_x12sa_dirname; 

%% filename
% specify funcion for compiling the filename. Input and output should be p and a 
% temporary structure, containing information about the currently prepared scan. 
% Cf. '+scans/get_filenames_cSAXS.m' 

det.get_filename = @scans.get_filenames_cSAXS;

% additional parameters for the filename function (det.get_filename)
% define filename pattern
% det.read_path_format = '%05d/';
% det.filename_pattern{1}.str = '%05d';
% det.filename_pattern{1}.del = '_';
% det.filename_pattern{1}.content = 'scan';
% det.filename_pattern{1}.start = 1;
% det.filename_pattern{2}.str = '%05d';
% det.filename_pattern{2}.del = '.';
% det.filename_pattern{2}.content = 'pos';
% det.filename_pattern{2}.start = 1;
% 
% det.data_prefix = 'run_';


% for scan = p.scan_number
% %     if ~isempty(dir(fullfile(p.raw_data_path{1}, 'eiger_4', utils.compile_x12sa_dirname(scan), '*.cbf')))
% %         det.file_extension = 'cbf';
% %     end
%     %%%% convert data from raw to cbf or hdf5
%     switch det.file_extension
%         case 'cbf'
%             detector.eiger1p5M.convert(scan, p.raw_data_path{1});
%         case 'h5'
%             detector.eiger1p5M.convert2hdf5(scan, p.raw_data_path{1});
%         otherwise
%             error('Unknown file extension for Eiger 1.5M.')
%     end
% 
% end


end

