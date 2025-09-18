%SETUP_REMOTE_RECONS
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

function setup_remote_recons()
import beamline.*

eaccount = beamline.identify_eaccount;

pgroup = strrep(eaccount, 'e', 'p');

fprintf('#########################################\n')
fprintf('Setup for remote reconstruction.\n')
fprintf('#########################################\n\n')


chk = false;
use_exchange = false;

while ~chk
    str_chk = input(sprintf('Do you want to use exchange instead of sshfs? [y/N]:'), 's');
    if strcmpi(str_chk, 'n')
        use_exchange = false;
        break;
    elseif strcmpi(str_chk, 'y')
        use_exchange = true;
        break;
    elseif isempty(str_chk)
        use_exchange = false;
        break;
    else
        fprintf('Please type "y" or "n".\n')
    end
end

if use_exchange
    
    ptycho_package_path = fullfile('/sls/X12SA/Data10/', eaccount, 'cxs_software/ptycho');
    base_package_path = fullfile('/sls/X12SA/Data10/', eaccount, 'cxs_software/base');
    exchange_dir = fullfile('/sls/X12SA/exchange/', pgroup, 'remote_recons');
    
    system(sprintf('cd /sls/X12SA/exchange/p17960; mkdir %s; chmod 770 %s', exchange_dir, exchange_dir));
    fid = fopen(fullfile(utils.abspath(exchange_dir), [eaccount '_receiver_template.m']), 'w');
    fprintf(fid, '%% ptycho receiver template\n');
    fprintf(fid, '%% execute this script on the remote machine\n\n');
    
    fprintf(fid, 'p.  verbose_level = 2;\n');
    
    str = sprintf('p.  base_path = ''%s'';',  exchange_dir);
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Base path for the remote reconstruction.\n']);
    
    str = sprintf('p.  ptycho_package_path = ''%s'';',  ptycho_package_path);
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% PtychoShelves path for the remote reconstruction.\n']);
    
    str = sprintf('p.  base_package_path = ''%s'';',  base_package_path);
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Base package path for the remote reconstruction.\n']);
    
    str = sprintf('p.  queue.remote_path = ''%s'';',  exchange_dir);
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Path to the remote queue. Needs to be accessible from (local) primary and (remote) replica.\n']);
    
    str = sprintf('p.  queue.recon_latest_first = true;');
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Reconstruct latest dataset first.\n']);
    
    str = sprintf('p.  raw_data_path{1} = ''%s'';', '/sls/X12SA/Data10/');
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Path to the raw data on Data10.\n']);
    
    str = sprintf('p.  orchestra.positions_file = [''%s''];', fullfile('/sls/X12SA/Data10/',eaccount,'specES1/scan_positions/scan_%%05d.dat'));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Orchestra path on Data10.\n\n\n']);

    str = sprintf('p.  detector.data_prefix = %s_1_;', eaccount);
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Data prefix.\n\n\n']);
    
    fprintf(fid, 'addpath(genpath(p.ptycho_package_path));\n');
    fprintf(fid, 'addpath(genpath(p.base_package_path));\n\n\n');
    fprintf(fid, 'while true\n\t core.run_receiver(p);\nend\n');
    
    fclose(fid);
    
    fid = fopen(fullfile(utils.abspath(exchange_dir), 'start-remote-session.sh'), 'w');
    fprintf(fid, ['source ' fullfile(base_package_path, 'setup-environment.sh\n')]);
    fprintf(fid, 'matlab -nodisplay -nodesktop -nosplash\n');
    fclose(fid);
    system(sprintf('chmod 770 %s', fullfile(exchange_dir, '*')));
    fprintf(['You can now open NX/NoMachine or ssh to the daas cluster and execute\n\t<strong>cd %s; source start-remote-session.sh</strong>\nThis will start a matlab session in which you can launch the receiver. Just type\n\t<strong>'...
        eaccount '_receiver_template' '</strong>\n\n'], exchange_dir);
    fprintf('Your replica should now be ready to receive files. Let''s go back to your primary template.\n')
    fprintf('Please change/set the following parameters:\n')
    fprintf('\tp.queue.remote_recons = true;\n')
    fprintf('\tp.queue.tmp_dir_remote = ''%s'';\n', exchange_dir)
    fprintf('\tp.queue.remote_path = ''%s'';\n\n', fullfile(exchange_dir, 'remote_queue'));
    fprintf('That''s it! Good luck!\n')
    
    
else
    input('\nCurrently sshfs is only supported on the gpu node.\nPlease make sure you are running this script on x12sa-gpu-1,\nor abort now and open a new matlab there (Continue).','s');
    
    pgroup_path1 = fullfile('/das/work/', pgroup(1:3), pgroup);
    pgroup_path2 = fullfile('/das/work/units/csaxs/', pgroup);
    
    ptycho_package_path = fullfile('/sls/X12SA/Data10/', eaccount, 'cxs_software/ptycho');
    base_package_path = fullfile('/sls/X12SA/Data10/', eaccount, 'cxs_software/base');
    
    
%    fprintf('Please make sure that you have created a directory as mounting point for the p-group.\n')
%    fprintf(sprintf('If this has not been done yet, please do it now.\nFor an eaccount %s and p-group %s, a directory like ~/Data10/%s would be good.\n', eaccount, pgroup, pgroup))
    fprintf(sprintf('\nCreating mounting point for the p-group, or enter the directory if already created.\nFor an eaccount %s and p-group %s, a directory like ~/Data10/%s would be good.\n', eaccount, pgroup, pgroup))
    
    chk = false;
    while ~chk
        str = input(sprintf('\nHow shall the directory be called? Please specify the full path (~/Data10/%s).:', pgroup), 's');
        if isempty(str)
            str = sprintf('~/Data10/%s', pgroup);
            break;
        else
            while ~chk
                str_chk = input(sprintf('The directory will be called "%s", okay? [y/n]:', str), 's');
                if strcmpi(str_chk, 'n')
                    fprintf('Okay, lets try again.\n')
                    break;
                elseif strcmpi(str_chk, 'y')
                    chk=true;
                else
                    fprintf('Please type "y" or "n".\n')
                end
            end
        end
    end
    fprintf('Creating directory %s.\n', str);
    gpu_dir = str;
    mkdir(gpu_dir);
    
    while true
        str = input('\nTo find the correct path on DaaS, we have to know whether the pgroup has been moved to cSAXS internal storage.\n If you are not sure, login to DaaS, go to your pgroup and type "pwd -P". If the path starts with /das/work/units/, the pgroup is already been moved to internal storage.\n Did you move the pgroup to cSAXS internal storage? (n) [y/n]:', 's');
        if strcmpi(str, 'n')
            pgroup_default = pgroup_path1;
            break;
        elseif strcmpi(str, 'y')
            pgroup_default = pgroup_path2;
            break;
        elseif isempty(str)
            pgroup_default = pgroup_path1;
            break;
        else
            fprintf('Please type "y" or "n".\n')
        end
    end
    
    chk = false;
    while ~chk
        str = input(sprintf('\nAssume you are on DaaS, please specify the full path of your working directory inside the p group (%s):', pgroup_default), 's');
        if isempty(str)
            pgroup_path = pgroup_default;
            break;
        else
            while ~chk
                str_chk = input(sprintf('The path on DaaS is "%s", okay? [y/n]:', str), 's');
                if strcmpi(str_chk, 'n')
                    fprintf('Okay, lets try again.\n')
                    break;
                elseif strcmpi(str_chk, 'y')
                    chk=true;
                    pgroup_path = str;
                else
                    fprintf('Please type "y" or "n".\n')
                end
            end
        end
    end
    fprintf('\nOkay, now it is time to mount the p group.\n')
    fprintf(sprintf('This can be achieved by \n\t sshfs <daas-username>@ra-l-004.psi.ch:%s %s\nPlease ssh to x12sa-gpu-1 and mount the p group.\n', pgroup_path, gpu_dir))
    while true
        str = input('Done? [y/n]:', 's');
        if strcmpi(str, 'n')
            fprintf('Please do it now.\n')
        elseif strcmpi(str, 'y')
            break;
        else
            fprintf('Please type "y" or "n".\n')
        end
    end
    fprintf('Great!\n\n')
    
    
    fid = fopen(fullfile(utils.abspath(gpu_dir), [eaccount '_receiver_template.m']), 'w');
    fprintf(fid, '%% ptycho receiver template\n');
    fprintf(fid, '%% execute this script on the remote machine\n\n');
    
    fprintf(fid, 'p.  verbose_level = 2;\n');
    
    str = sprintf('p.  base_path = ''%s'';',  add_delimiter(pgroup_path));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Base path for the remote reconstruction.\n']);
    
    
    chk = false;
    while ~chk
        str = input(sprintf('Please specify the full PtychoShelves path / ptycho_package_path (%s):', ptycho_package_path), 's');
        if isempty(str)
            break;
        elseif strcmpi(str(1),'~') || strcmpi(str(1),'.')
            fprintf('Please do not use relative paths.\n')
        else
            while ~chk
                str_chk = input(sprintf('The new path is "%s", okay? [y/n]:', str), 's');
                if strcmpi(str_chk, 'n')
                    fprintf('Okay, lets try again.\n')
                    break;
                elseif strcmpi(str_chk, 'y')
                    chk=true;
                    ptycho_package_path = str;
                else
                    fprintf('Please type "y" or "n".\n')
                end
            end
        end
    end
    str = sprintf('p.  ptycho_package_path = ''%s'';',  add_delimiter(ptycho_package_path));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% PtychoShelves path for the remote reconstruction.\n']);
    
    chk = false;
    while ~chk
        str = input(sprintf('Please specify the full base package path / base_package_path (%s):', base_package_path), 's');
        if isempty(str)
            break;
        elseif strcmpi(str(1),'~') || strcmpi(str(1),'.')
            fprintf('Please do not use relative paths.\n')
        else
            while ~chk
                str_chk = input(sprintf('The new path is "%s", okay? [y/n]:', str), 's');
                if strcmpi(str_chk, 'n')
                    fprintf('Okay, lets try again.\n')
                    break;
                elseif strcmpi(str_chk, 'y')
                    chk=true;
                    base_package_path = str;
                else
                    fprintf('Please type "y" or "n".\n')
                end
            end
        end
    end
    str = sprintf('p.  base_package_path = ''%s'';',  add_delimiter(base_package_path));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Base package path for the remote reconstruction.\n']);
    
    remote_queue_path = fullfile(pgroup_path, 'remote_queue');
    str = sprintf('p.  queue.remote_path = ''%s'';',  add_delimiter(remote_queue_path));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Path to the remote queue. Needs to be accessible from (local) primary and (remote) replica.\n']);
    
    str = sprintf('p.  queue.recon_latest_first = true;');
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Reconstruct latest dataset first.\n']);
    
    str = sprintf('p.  raw_data_path{1} = ''%s'';', fullfile('/sls/X12SA/Data10/', eaccount));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Path to the raw data on Data10.\n']);
    
    str = sprintf('p.  orchestra.positions_file = [''%s''];', fullfile(utils.abspath('~/'), 'specES1/scan_positions/scan_%%05d.dat'));
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Orchestra path on Data10.\n\n\n']);
    
    str = sprintf('p.  detector.data_prefix = ''%s_1_'';', eaccount);
    fprintf(fid, [str repmat(' ', 1, max(75-length(str), 0)) '%% Data prefix.\n\n\n']);
    
    fprintf(fid, 'addpath(genpath(p.ptycho_package_path));\n');
    fprintf(fid, 'addpath(genpath(p.base_package_path));\n\n\n');
    fprintf(fid, 'while true\n\t core.run_receiver(p);\n\tpause(2);\nend\n');
    
    fclose(fid);
    
    fid = fopen(fullfile(utils.abspath(gpu_dir), 'start-remote-session.sh'), 'w');
    fprintf(fid, ['source ' fullfile(base_package_path, 'setup-environment.sh\n')]);
    fprintf(fid, 'matlab -nodisplay -nodesktop -nosplash\n');
    fclose(fid);
    
    fprintf(['You can now open NX/NoMachine or ssh to the daas cluster and execute\n\tsource start-remote-session.sh\nThis will start a matlab session in which you can launch the receiver. Just type\n\t'...
        eaccount '_receiver_template' '\n\n']);
    fprintf('Your replica should now be ready to receive files. Let''s go back to your primary template.\n')
    fprintf('Please change/set the following parameters:\n')
    fprintf('\tp.queue.remote_recons = true;\n')
    fprintf('\tp.queue.tmp_dir_remote = ''%s'';\n', add_delimiter(gpu_dir))
    fprintf('\tp.queue.remote_path = ''%s'';\n\n', fullfile(gpu_dir, 'remote_queue'));
    
    fprintf('Settings to run receiver in login node with c_solver and let the code allocate a compute node: \n')
    fprintf('\teng.  use_gpu = false;\n')
    fprintf('\teng.  gpu_id = [];\n')
    fprintf('\teng.  num_gpus = 0;  \n')  
    fprintf('\teng.  ra_nodes = 1; \n')
    fprintf('\tfor the path of p.initial_probe_fileprobe, do not use ~, use /sls/X12SA/Data10/e17966/\n\n')
    
    fprintf('If the code crashes and you want to clean up and try again, don''t forget to clean the remote queue folder ~/Data10/p17966/remote_queue/\n\n')
    
    fprintf('That''s it! Good luck!\n')
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
