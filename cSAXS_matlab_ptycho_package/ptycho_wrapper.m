%PTYCHO_WRAPPER launch multiple matlab sessions using screen and monitor
%their workload
% ptycho_wrapper()
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


function ptycho_wrapper()
import utils.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% modify the following section %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_nodes = 10;
file_queue = '/das/work/units/tomcat/p17371/offline/reconstruction/';
recon_file = 'moench_offline';
setup_environment = '/das/work/units/tomcat/p17371/offline/cSAXS_matlab_base/setup-environment.sh';
ptycho_path = '/das/work/units/tomcat/p17371/offline/cSAXS_matlab_ptycho';


fext = 'mat';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s = [];
s.name = '';
s.idle = true;
s.initialized = false;


% initialize
for ii=1:max_nodes
    s(ii).name = sprintf('ps_wrapper_%u', ii);
    s(ii).idle = true;
    s(ii).initialized = false;
    s(ii).connected = true;
    system(sprintf('screen -S %s -dm', s(ii).name));
    
end

finishup = utils.onCleanup(@(x) wrapper_exit(x), s);


for ii=1:max_nodes
    cmd = sprintf('stuff "source %s\n cd %s; matlab -nodisplay -nosplash -nodesktop\n"', setup_environment, ptycho_path);
    s(ii).initialized = true;
    s(ii).idle = true;
    finishup.update(s);
    send_cmd(s(ii).name, cmd);
end



% 'add worker' cmd
cmd = sprintf('stuff "%s\n"', recon_file);

% start first worker
s(1).idle = false;
send_cmd(s(1).name, cmd);
pause(2)

while true
    % check queue
    fqueue = dir(fullfile(file_queue,['*.' fext]));
    progQueue = dir(fullfile(file_queue, 'in_progress', ['*.' fext]));
    if isempty(fqueue) && isempty(progQueue)
        break
    end
    % more workers?
    if length(fqueue)>2
        % add new workers, if possible
        idle_workers = find([s.idle]~=0);
        if ~isempty(idle_workers)
            fprintf('Adding a new worker...\n');
            send_cmd(s(idle_workers(1)).name, cmd);
            s(idle_workers(1)).idle = false;
        end
    elseif numel(find([s.idle]==0))>length(progQueue)+2
        % stop a worker
        try
            fprintf('Removing a worker...\n');
            active_workers = find([s.idle]==0);
            stop_matlab(s(active_workers(end)).name);
            s(active_workers(end)).idle = true;
        catch
            fprintf('Failed to remove a worker.')
        end
    end
    pause(2)
    
end

end



function wrapper_exit(s)
% cleanup function for ptycho_wrapper - make sure that all screen sessions
% are terminated

fprintf('Closing sessions, please wait...\n')
for ii=1:numel(s)
    if s(ii).initialized || s(ii).connected
        fprintf('Closing session %u...\n', ii);
        stop_matlab(s(ii).name)
%         pause(0.5)
        send_cmd(s(ii).name, 'quit')
    else
        continue
    end
end
fprintf('Done.\n')

end

function send_cmd(session, cmd)
% send a command to a screen session
    system(sprintf('screen -x %s -X %s', session, cmd));
end

function stop_matlab(session)
% stop matlab in a screen session
    for jj=1:10 
        send_cmd(session, 'stuff $"\003"')
        pause(0.2)
    end
end
