%PTYCHO_PREPARE_PATHS Prepare paths and check defaults

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

function [ p ] = ptycho_prepare_paths( p, varargin)
import utils.*


if nargin > 1
    init = varargin{1};
else
    init = false;
end

if init
    verbose(p.verbose_level);

    verbose(1, 'Preparing paths.')

    % base paths
    if isempty(p.base_path)
        p.base_path = './';
        verbose(3, 'Using default value for base_path: %s', p.base_path)
    end
    
    p.base_path = abspath(p.base_path);
    
    
    % specfile
    if isfield(p, 'specfile')
        if isempty(p.specfile) %&& isfield(p, 'src_metadata') && strcmp(p.src_metadata, 'spec')
            p.specfile = p.base_path;
            verbose(3, 'Using default value for specfile: %s', p.specfile )
        end
        p.specfile = abspath(p.specfile);
    end
    
    
    % ptycho path
    if ~isfield(p, 'ptycho_package_path') && isfield(p, 'ptycho_matlab_path')
        %legacy
        p.ptycho_package_path = p.ptycho_matlab_path;
    end
    if isempty(p.ptycho_package_path)
        % find template_ptycho and assume that the base path is there 
        p.ptycho_package_path = fileparts(which('template_ptycho.m'));
        verbose(3, 'Using default value for ptycho_package_path: %s', p.ptycho_package_path )
    end

    % base package path
    % legacy:
    if ~isfield(p, 'base_package_path') && isfield(p, 'cSAXS_matlab_path')
        p.base_package_path = p.cSAXS_matlab_path;
    end
    if isempty(p.base_package_path)
        if exist(fullfile(p.base_path,'cxs_software/base/+math'), 'dir')
            p.base_package_path = fullfile(p.base_path,'cxs_software/base');
            verbose(3, 'Using default value for base_package_path: %s', p.base_package_path )
        end
    end
    
    % raw data path (part I)
    if isempty(p.raw_data_path) || iscell(p.raw_data_path)&&isempty(p.raw_data_path{1})
        if isempty(p.scan_number)
            p.raw_data_path{1} = p.base_path;
        else
            for ii = 1:length(p.scan_number)
                p.raw_data_path{ii} = p.base_path;
            end
        end
    end
    % data prefix
    if isempty(p.detector.data_prefix)
        import beamline.identify_eaccount  %% not included in the ptychoshelves package
        eaccount = identify_eaccount;
        if ~isempty(eaccount) && eaccount(1) == 'e'
            % default setting for cSAXS beamline
            p.detector.data_prefix = [eaccount '_1_'];
        else
            verbose(0,'p.detector.data_prefix is not defined')
        end
    end

    
    % do some basic corrections of the paths 
    for path = {'base_path', 'ptycho_package_path', 'base_package_path', 'prepare_data_path', 'positions_file', 'raw_data_path'}
        if isfield(p, path{1}) && ~isempty(p.(path{1}))
            if iscell(p.(path{1}))
                for ii=1:length(p.(path{1}))
                    p.(path{1}){ii} = abspath(p.(path{1}){ii});
                end
            else
                p.(path{1}) = abspath(p.(path{1}));
            end
        end
    end

    
    %% add paths, but only if not included already 

    if exist(p.ptycho_package_path, 'dir') && ~exist(fullfile('+core', 'get_projections.m'),'file') 
        addpath(p.ptycho_package_path)
        % check if ptycho_package_path is already included 
    elseif ~exist( fullfile('+core', 'get_projections.m'), 'file')
        verbose(1,'Nonexistent ptycho_package_path: "%s"', p.ptycho_package_path)
    end
    
    if exist(fullfile(p.ptycho_package_path, 'utils'), 'dir') && ~exist(fullfile('aligned_FSC.m'),'file') 
        addpath(fullfile(p.ptycho_package_path, 'utils'))
    % check if ptycho_package_path/utils is already included 
    elseif ~exist(fullfile('aligned_FSC.m'),'file') 
        verbose(1,'Nonexistent ptycho_package_path: "%s/utils"', p.ptycho_package_path)
    end
    
    if exist(p.base_package_path, 'dir') && ~exist(fullfile('+math', 'argmax.m'), 'file')  
        addpath(p.base_package_path);
    % check if base_package_path is already included 
    elseif ~exist(fullfile('+math', 'argmax.m'), 'file')  
        verbose(1,'Nonexistent base_package_path: "%s"', p.base_package_path)
    end
    verbose(p.verbose_level);

    
else
    
    % base path
    if ~exist(p.base_path, 'dir')
        error('base_path = %s : Base directory does not exist.', p.base_path);
    end
    verbose(2, 'base_path = %s', p.base_path);
    
    % Save data path
    for ii = 1:length(p.scan_number)
        p.   scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
    end
    
    if isempty(p.save_path) || iscell(p.save_path)&&isempty(p.save_path{1})
        verbose(3, 'Using default save_path');
        for ii = 1:length(p.scan_number)
            p.save_path{ii} = fullfile(p.base_path, 'analysis',utils.compile_x12sa_dirname(p.scan_number(ii)),'');
            if ~exist(p.save_path{ii}, 'dir')
                mkdir(p.save_path{ii})
            end
            verbose(2, 'save_path = %s', p.save_path{ii});
        end
    else
        % not enough save paths; replicate
        if length(p.scan_number) > length(p.save_path)
            verbose(1, 'Number of save paths does not match number of scans. I will use only the first save path.')
            
            if contains(p.save_path{1}, '%')
                % fill in the scan number if needed
                container = p.save_path{1};
                for ii=1:length(p.scan_number)
                    p.save_path{ii} = sprintf(container, p.scan_number(ii));
                end
            else
                % if there is no variable to fill, append the scan number to
                % the given path
                container = p.save_path{1};
                for ii=1:length(p.scan_number)
                    p.save_path{ii} = fullfile(container, p.scan_str{ii});
                end
            end
            clear container

        else
            if contains(p.save_path{1}, '%')
                % fill in the scan number if needed
                for ii=1:length(p.scan_number)
                    p.save_path{ii} = sprintf(p.save_path{ii}, p.scan_number(ii));
                end
            else
                % if there is no variable to fill, append the scan number to
                % the given path
                for ii=1:length(p.scan_number)
                    p.save_path{ii} = fullfile(p.save_path{ii}, p.scan_str{ii});
                end
            end
        end
        
        for ii = 1:length(p.scan_number)
            p.save_path{ii} = rm_delimiter(p.save_path{ii});
            if ~exist(p.save_path{ii}, 'dir')
                mkdir(p.save_path{ii})
            end
            verbose(2, 'save_path = %s', p.save_path{ii});
        end

    end
    % Prepare data path
    if isempty(p.prepare_data_path)|| iscell(p.prepare_data_path)&&isempty(p.prepare_data_path{1})
        verbose(3, 'Using default prepared data path');
        p.prepare_data_path = p.save_path{1};
    else
        if iscell(p.prepare_data_path)
            p.prepare_data_path = cell2str(p.prepare_data_path);
        end
        if contains(p.prepare_data_path, '%')
            % fill in the scan number if needed
            p.prepare_data_path = sprintf(p.prepare_data_path, p.scan_number(1));
        else
            % if there is no variable to fill, append the scan number to
            % the given path
            p.prepare_data_path = fullfile(p.prepare_data_path, p.scan_str{1});
        end
    end
        p.prepare_data_path = replace(p.prepare_data_path,'\','\\');
        p.prepare_data_path = add_delimiter(p.prepare_data_path);
    if ~exist(p.prepare_data_path, 'dir')
        mkdir(p.prepare_data_path)
    end
    verbose(2, 'prepare_data_path = %s', p.prepare_data_path);
    
    
    
    % prepare data filename
    if isempty(p.prepare_data_filename)
        verbose(3, 'Using default prepared data filename');
        if ~p.detector.binning
            p.prepare_data_filename = [sprintf('S%05d_data_%03dx%03d',p.scan_number(1), p.asize(1), p.asize(2)) p.prepare.prep_data_suffix '.h5'];
        else
            p.prepare_data_filename = [sprintf('S%05d_data_%03dx%03d_b%i',p.scan_number(1), p.asize(1), p.asize(2),p.detector.binning) p.prepare.prep_data_suffix '.h5']; 
        end
        verbose(2, 'prepare_data_filename = %s', p.prepare_data_filename);
    end
    
    
    % raw data path
    if length(p.raw_data_path) ~= length(p.scan_number)
        if ~iscell(p.raw_data_path)
            p = str2cell(p, 'raw_data_path');
        end
        for ii = 2:length(p.scan_number)
            p.raw_data_path{ii} = p.raw_data_path{1};
        end
    end
    if p.prepare.auto_prepare_data
        % do some basic replacement to get the real path
        for ii = 1:length( p.raw_data_path)
            p.raw_data_path{ii} = abspath(p.raw_data_path{ii});
            p.raw_data_path{ii} = sprintf(replace(p.raw_data_path{ii},'\','\\') , p.scan_number(1));
        end
    else
        prepare_data_full_filename = fullfile(p.prepare_data_path, p.prepare_data_filename);
        if ~exist(prepare_data_full_filename, 'file'); error('prepared data file does not exist (%s).', prepare_data_full_filename); end
    end
        
    % do some basic corrections of the paths
    for path = {'base_path', 'specfile', 'ptycho_package_path', 'base_package_path', 'prepare_data_path'}
        p.(path{1}) = abspath(p.(path{1}));
    end
    
end
end

% convert single cell entry to string
function p = cell2str(p, fn)
    tmp = p.(fn){1};
    p = rmfield(p, fn);
    p.(fn) = tmp;

end

% convert string to cell
function p = str2cell(p, fn)
    tmp = p.(fn);
    p = rmfield(p, fn);
    p.(fn){1} = tmp;

end

% make sure that the path does not end with /
function path = rm_delimiter(path)
if ~isempty(path) && any(strcmp(path(end), {'\', '/'})) 
    path = path(1:end-1);
end
end

% make sure that the path ends with /
function path = add_delimiter(path)
if ispc
    delimiter = '\';
else
    delimiter = '/';
end
if ~isempty(path) && ~strcmp(path(end), delimiter)
    path = [path, delimiter];
end
end

