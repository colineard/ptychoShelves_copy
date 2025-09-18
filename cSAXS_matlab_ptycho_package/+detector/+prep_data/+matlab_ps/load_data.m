%LOAD_DATA prepare filenames and load data
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: detector.prep_data.matlab_ps.prepare_data
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

function [ p ] = load_data( p )
import utils.find_files
import io.image_read

det.params = p.detectors(p.scanID).params;
detStorage = p.detectors(p.scanID).detStorage;

%% prepare filenames
if isempty(detStorage.files)
    [p] = det.params.get_filename(p);
    if numel(detStorage.files) == 0
        error('Did not find any files.')
    end
end
files = detStorage.files;

%% select data loading routine
current_version = version;
ver_str = strsplit(current_version, '.');
ver_num = str2double([ver_str{1} '.' ver_str{2} ver_str{3}]);
if ver_num >= 9.4 
    c_reader = false;
else
    utils.verbose(3, 'Fast data reader is not available for Matlab version %d. Switching to image_read.', ver_num)
    c_reader = false;
end
if c_reader && ~ismember(det.params.file_extension, {'h5', 'cbf', 'tiff', 'tif'})
    utils.verbose(3, 'Fast data reader is not available for selected file format %s. Switching to image_read.', det.params.file_extension)
    c_reader = false;
end

%% get data orientation
if isempty(det.params.orientation) 
    if strcmpi(det.params.file_extension, 'h5')
        det.params.orientation = utils.get_detector_orientation(det.params.H5_detector_path, 'filename', files{1}, 'orientType', 'transpose_fliplr_flipud');
    else
        error('Could not retrieve data orientation from raw data file. Please specify the orientation in your detector template (+detector/+<detector_name>/<detector_name>.m)')
    end
end

%% get data size
if strcmpi(det.params.file_extension, 'h5')
    dataInfo = h5info(files{1}, det.params.H5_detector_path);
    if isfield(dataInfo, 'Dataspace')
        det.params.data_size = dataInfo.Dataspace.Size(1:2);
        if isempty(detStorage.h5_group)
            detStorage.h5_group{1} = det.params.H5_detector_path;
        end
    elseif isfield(dataInfo, 'Datasets')
        jj = 1;
        while true
            if strcmpi(dataInfo.Datasets(jj).Name, 'data')
                det.params.data_size = dataInfo.Datasets(jj).Dataspace.Size(1:2);
                break
            end
            jj = jj+1;
        end
        if jj==numel(dataInfo.Datasets(jj))
            error('Could not find entry "data".')
        end
        if isempty(detStorage.h5_group)
            detStorage.h5_group{1} = fullfile(det.params.H5_detector_path, dataInfo.Datasets(jj).Name);
        end
    else
        error('Failed to retrieve data dimensions from H5 file.')
    end
    if det.params.orientation(1)
        det.params.data_size = fliplr(det.params.data_size);
    end
else
    % read first file to get orientation
    dataInfo = io.image_read(files{1}, 'orientation', [0 0 0]);
    det.params.data_size = size(dataInfo.data(:,:,1));
    if det.params.orientation(1)
        det.params.data_size = fliplr(det.params.data_size);
    end   
end

%% convert center from user orientation to raw data orientation
ctr = math.convert_location('pos', detStorage.ctr, 'framesz', ...
    det.params.data_size, 'vecIn', det.params.orientation, 'vecTypeIn', 'transpose_fliplr_flipud');

%% load data
if c_reader
    % use fast data reader
    if strcmpi(det.params.file_extension, 'tif')
        det.params.file_extension = 'tiff'; 
    end
    if strcmpi(det.params.file_extension, 'tiff')
        % keep results consistent with matlab's imread and io.image_read for tiff files 
        det.params.orientation(1) = ~det.params.orientation(1); 
    end
    arg.ctr = ctr - 1; % C convention
    
    % flip XY
    arg.ctr = fliplr(arg.ctr);
    
    % create structure for c_reader
    arg.data_path = files;
    arg.nthreads = p.io.data_nthreads;
    arg.precision = p.io.data_precision;
    arg.extension = det.params.file_extension;
    arg.asize = detStorage.read_size;
    arg.data_location = detStorage.h5_group;
    assert(all(arg.ctr>0),'Raw data center position has to be positive.')
    %% load data 
    utils.verbose(2, 'Loading raw data of scan S%05d.', p.scan_number(p.scanID))
    try
        data = io.read_measurement(arg);
    catch
        error('Data loading failed. Please check that the data exists and make sure that the detector parameters are set correctly in (%s).', ['+detector/+' p.detector.name '/+' p.detector.name '.m']);
    end
    
    %% apply orientation
    data = squeeze(data);
    if det.params.orientation(1)
        data = permute(data, [2 1 3]);
    end
    if det.params.orientation(2) && det.params.orientation(3)
        data = rot90(data,2);  % merge the fliplr and flipup operations 
    else
        if det.params.orientation(2)
            data = fliplr(data);
        end
        if det.params.orientation(3)
            data = flipud(data);
        end
    end
    
    
    
else
    
    %% prepare image_read_extraags
    det.params.image_read_extraargs{end+1} = 'RowFrom';
    det.params.image_read_extraargs{end+1} = detStorage.lim_inf(1);
    det.params.image_read_extraargs{end+1} = 'RowTo';
    det.params.image_read_extraargs{end+1} = detStorage.lim_sup(1);
    det.params.image_read_extraargs{end+1} = 'ColumnFrom';
    det.params.image_read_extraargs{end+1} = detStorage.lim_inf(2);
    det.params.image_read_extraargs{end+1} = 'ColumnTo';
    det.params.image_read_extraargs{end+1} = detStorage.lim_sup(2);
    det.params.image_read_extraargs{end+1} = 'Orientation';
    det.params.image_read_extraargs{end+1} = det.params.orientation;
    det.params.image_read_extraargs{end+1} = 'OrientByExtension';
    det.params.image_read_extraargs{end+1} = 0;   

    %% load data
    if isempty(detStorage.h5_group)
        if numel(files)==1
            % one directory contains single file
            dataaux = image_read(files{1}, det.params.image_read_extraargs);
            data = dataaux.data;
        else
            % multiple files per directory
            dataaux = image_read(files, det.params.image_read_extraargs);
            data = dataaux.data;
        end
    else
        data = zeros([p.asize numel(detStorage.h5_group)]);
        
        if numel(files)==1
            if numel(detStorage.h5_group)==1
                % one hdf5 file; one group
                dataaux = image_read(files{1}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{1});
                data = dataaux.data;
            else
                % one hdf5 file; multiple groups
                for ii=1:numel(detStorage.h5_group)
                    dataaux = image_read(files{1}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{ii});
                    data(:,:,ii) = dataaux.data;
                end
            end
        else
            if numel(detStorage.h5_group)==1
                % multiple files, single H5 group
                dataaux = image_read(files, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{1});
                data = dataaux.data;
            else
                % multiple files, but different H5 group
                for ii=1:numel(detStorage.h5_group)
                    dataaux = image_read(files{ii}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{ii});
                    data(:,:,ii) = dataaux.data;
                end
            end
        end
    end
end
    
detStorage.data = data;

end
