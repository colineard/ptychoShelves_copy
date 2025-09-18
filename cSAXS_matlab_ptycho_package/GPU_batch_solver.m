%GPU_BATCH_SOLVER 
% multi GPU batch solver
%
% use parallel toolbox to distribute tasks into several asynchronously
% running workers


%% HOW TO USE BATCH SOLVER: 
% 1) prepare a template to load dat from filelist queue
% 2) dont use while loop in the template !!!
% 3) test that the template works ok, !! reconstruct at least one scan with it manually !!
% 4) modify parameters below for your scan 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%% SET PARAMETERS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nthread_per_GPU = 3;                        % how many workers will use one GPU, the workers will wait till there is enough memory to run them
GPU_ids = 1:4;                              % IDs of the used GPUs 
executed_template = 'ptycho_wrapper_GPU';   % name of the script to be run in parallel, TEST FIRST MANUALLY ONE RECONSTRUCTION 
base_path = '../../';                       % place to store temporal data, log files etc 
logging_folder = 'log_batch_solver';        % folder where the generated log files are stored 
final_cleanup = true;                      % clean the jobs when the solver is manually interrupted, if false, the worker will be running in background after CTRC+C interruption  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist(executed_template, 'file')~=0, sprintf('Missing script %s', executed_template))
cSAXS_base_path = core.find_base_package; 
addpath(cSAXS_base_path)
addpath('./utils')
utils.verbose(struct('prefix', 'batch_solver'))

% get p-structure but without data yet 
run(executed_template)

%%  set some default options %%%%%
p.use_display = false; % avoid plotting on screen
p.verbose_level = 1; % set the verbose level to be low but nonzero
p.save.external = false; % avoid opening another matlab process 
p.cSAXS_base_path = cSAXS_base_path; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


assert(~isempty(p.queue.path), 'Set p.queue.path value in template')
assert(~isempty(p.queue.name), 'Set p.queue.name value ''filelist'' in template')

% remove unimportant warning 
warning('off','parallel:task:DiaryIncomplete')
warning('off','parallel:batch:DiaryIncomplete')
warning('off','parallel:cluster:DepfunError')
warning('off','MATLAB:MKDIR:DirectoryExists')


%% PREPARE WORKERS 

cluster = parcluster('local');
if exist([base_path,'/local_cluster_jobs'], 'dir')
    rmdir([base_path,'/local_cluster_jobs'], 's')  % delete folder with jobs (prevent accumulation ) 
end
mkdir([base_path,'/local_cluster_jobs'])   % recreate the folder 
cluster.JobStorageLocation = [base_path,'/local_cluster_jobs'];


utils.verbose(0,'Starting jobs ... ')
[jobs, tasks, log_files, task_StartDateTime] =  init_jobs(cluster, p, executed_template, GPU_ids, Nthread_per_GPU,logging_folder);

%% RUN THE TASKS AND REPORTS  IN LOOP 

run_batch_solver(jobs, tasks, log_files, task_StartDateTime, p, final_cleanup)


%% manually clean all jobs if not done automatically (final_cleanup == false)
for job = jobs
    cancel(job)
    delete(job)
end


%%  USE FOLLOWING COMMANDS TO TERMINATE JOBS THAT ARE RUNNING ON DAAS CPU NODES
% !squeue
% % use if you need to cancel jobs on the Daas CPU nodes , fill in your user name 
% ! for id in  `squeue  | grep odstrcil | awk '{print $1}'`; do scancel $id ; done 
    






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FUNCTION DEFINITIONS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  p = worker(p)
    % worker is an infinite loop running ptychography for given template 
    
    warning off backtrace
    addpath('./utils')
    while true
        p.run_name = '';  % some bug in Klaus code
        fprintf('=== Started ===\n') % markers for the log parser 

        [~, status] = core.ptycho_recons(p);
        if ~status
            fprintf('=== Failed ===\n') % markers for the log parser 
            pause(10); % give a chance to the reporting script to detect this event 
        else
            fprintf('=== Finished ===\n')  % markers for the log parser 
            pause(5); % give a chance to the reporting script to detect this event 
        end
        
        try
            reset(gpuDevice( p.gpu_id ))
        end

    end
end


function [jobs, tasks, log_files, task_StartDateTime] =  init_jobs(cluster, p, executed_template, GPU_ids, Nthread_per_GPU,logging_folder)
    % initialize jobs , create log files

    N_GPU = length(GPU_ids);
    log_files = {};
    % reset all used GPUs, it may help if there was an error before 
    for id = GPU_ids
        reset(gpuDevice(id)); 
    end
    if ~exist(logging_folder, 'dir')
       mkdir(logging_folder) 
    end
    for j = 1:Nthread_per_GPU
        for i = 1:N_GPU
            % choose GPU to use 
            GPU_id = GPU_ids(i);
            job_id = i+(j-1)*N_GPU; 
            p.gpu_id = GPU_id;   % set the default GPU
            p = set_all_engines(p, 'gpu_id', GPU_id);
            utils.verbose(0,'Opening job %i on GPU %i thread %i',job_id, GPU_id, j )
            % run the worker in background 
            jobs(job_id) = createJob(cluster);
            tasks(job_id) = createTask(jobs(job_id),@worker,1,{p}, 'CaptureDiary', true);
            log_files{job_id} = sprintf('%s/%s_job_%02i.log', logging_folder, executed_template, job_id );
            task_StartDateTime(job_id) = now; 
        end
    end
    utils.verbose(0,'All jobs initialized ... ')
    Njobs = job_id;

    utils.verbose(0,'Log files: ')
    for job_id = 1:Njobs
        disp(log_files{job_id})
        if exist(log_files{job_id}, 'file'); delete(log_files{job_id}); end 
    end
end



function run_batch_solver(jobs, tasks, log_files, task_StartDateTime, p, final_cleanup)
    % iterate through all the threads, check their status and report if
    % something interesting have happened. If the thread failed, try to
    % resubmit new job 
    
    
    utils.verbose(0,'Submitting,  monitoring jobs and saving into logs ...    ')
    Njobs = length(tasks); 
    task_stopped = zeros(Njobs,1);
    task_waiting = false(Njobs,1);
    task_finished = false(Njobs,1);

    pause_time = 1; % time between the log updates 
    dead_time = 60; % time to report task dead 
    msg_status = 0; 
    t_start = tic; 
    use_gui = usejava('desktop'); % find if GUI or command line is used 

    % prepare automatic cleaner that will kill all jobs when function is
    % interruped 
    if final_cleanup
        finishup = onCleanup(@()myCleanupFun(jobs,p.queue.path)); 
    end
    
    while true 
        pause(pause_time)
        % Check RAM 
        %utils.check_available_memory();


        for job_id = 1:Njobs
            file = dir(log_files{job_id});
            if isempty(file)
                fileID = fopen(log_files{job_id},'w');
                fprintf(fileID,'Submitting ...  ');
                fclose(fileID);  
                continue 
            end

            if strcmp(jobs(job_id).State, 'pending')
                utils.verbose(0,'Submitting job %i',job_id); msg_status = 0; 
                try
                    submit(jobs(job_id));
                catch err
                    warning(err.message )
                end
                utils.verbose(0,'%s -- Job %i was submitted  ',datetime, job_id);msg_status = 0; 
                continue 
            end

            %%%%%%%%% update logs %%%%%%%%%%%%%%
            [log_size, log_size_0] = update_logs(file, tasks(job_id), log_files{job_id}); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if log_size == 0  && toc(t_start) > 20*Njobs
                warning('Log in task %i is empty', job_id); 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %% RESUBMIT TASKS IF FINISHED %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if strcmp(tasks(job_id).State, 'finished')  ||  task_finished(job_id)
                [Y,M,D,H,MN,S] = datevec(now -  task_StartDateTime(job_id)); 
                job_msg = sprintf('%s -- Job %i report: -  run time: %i:%i:%i, worker status: "%s"   ',datetime, job_id, H, MN,round(S), tasks(job_id).State);
                report_outputs(log_files{job_id}, job_msg); msg_status = 0; 
                
                
                if strcmp(tasks(job_id).State, 'finished') 
                    % if the task is really finished, e.g. crashed, delete is and open solver again 
                    if ~isempty(tasks(job_id).Error)
                        utils.verbose(0,'------ Error ------- '); msg_status = 0; 
                        disp(getReport( tasks(job_id).Error, 'extended', 'hyperlinks', 'on' ))
                    end
                    keyboard
                    p = tasks(job_id).InputArguments{1}; 
                    try; jobs(job_id).Tasks(1).delete; end
                    utils.verbose(0,'%s -- Resubmitting job %i on GPU %i  ',datetime,job_id, p.gpu_id);msg_status = 0; 
                    jobs(job_id) = jobs(job_id).recreate;
                    tasks(job_id) = createTask(jobs(job_id),@worker,1,{p}, 'CaptureDiary', true);
                end
                task_stopped(job_id) = false;
                task_finished(job_id) = false;
                task_StartDateTime(job_id) = now; 
                continue 
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% PARSE LOGS AND AUTOREPORT  %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,cmdout] =system(['tail -n2 ', log_files{job_id},  '| grep "Started solver" | wc -l']);
            if str2num(cmdout) == 0 && log_size_0 == log_size && task_stopped(job_id) <= dead_time
                task_stopped(job_id) =  task_stopped(job_id) + pause_time;
                if  task_stopped(job_id) > dead_time
                    utils.verbose(0,'%s -- Job %i is probably failed, see details in %s   ', datetime, job_id, log_files{job_id}); msg_status = 0; 
                end
            elseif strcmp(jobs(job_id).State, 'running') && log_size_0 ~= log_size && task_stopped(job_id) > dead_time
                utils.verbose(0,'%s -- Job %i continues, see details in %s', datetime, job_id, log_files{job_id});msg_status = 0; 
                task_stopped(job_id) = 0;
            elseif log_size_0 ~= log_size
                task_stopped(job_id) = 0; 
            end

            % check for serious error
            [~,cmdout]  =system(['tail -n10 ', log_files{job_id},  '| grep "Error in distcomp_evaluate_filetask_core" | wc -l']);

            if (str2num(cmdout))&& ~task_stopped(job_id)
                utils.verbose(0,'%s -- Job %i propably failed, see details in %s', datetime, job_id, log_files{job_id});msg_status = 0; 
                task_finished(job_id) = true;
            end

             % check for well handle error
            [~,cmdout1] =system(['tail -n10 ', log_files{job_id},  '| grep "Reconstruction stopped!" | wc -l']);
            [~,cmdout2] =system(['tail -n50 ', log_files{job_id},  '| grep "## PONG! ##" | wc -l']);
            if (str2num(cmdout1) || str2num(cmdout2))&& ~task_finished(job_id)
                utils.verbose(0,'%s -- Job %i reconstruction failed, see details in %s', datetime, job_id, log_files{job_id});msg_status = 0; 
                task_finished(job_id) = true;
            end
                        
            % check for queue
            [~,cmdout] =system(['tail -n5 ', log_files{job_id},  '| grep "Did not find enough files in queue" | wc -l']);
            if str2num(cmdout) && ~task_waiting(job_id)
                utils.verbose(0,'%s -- Job %i exitting because no data are available   ', datetime, job_id);msg_status = 0; 
                task_waiting(job_id) = true; 
            elseif str2num(cmdout) == 0
                task_waiting(job_id) = false;
            end

            % check for finished 
            [~,cmdout] =system(['tail -n1 ', log_files{job_id},  '| grep "=== Finished ===" | wc -l']);
            t = datevec(now -  task_StartDateTime(job_id)); 
            if str2num(cmdout)>0 && ~task_finished(job_id) && ~task_waiting(job_id) && (t(5) > 1)
                utils.verbose(0,'%s -- Job %i is finished  ',datetime, job_id);msg_status = 0; 
                task_finished(job_id) = true; 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% show current status %%  
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            jobs_stats = [sum(task_finished), sum(task_waiting), sum(task_stopped>dead_time)]; 
            if all(task_waiting)
                string = sprintf('[batch_solver] :  All jobs finished, waiting for new data ... ');
            else
                string = sprintf('[batch_solver] : Running jobs: %2i/%2i  --  Finished: %2i\tWaiting: %2i\tStopped: %2i -- Scans in queue: %4i ... ',...
                    Njobs-sum(task_finished|task_waiting|task_stopped>dead_time) , Njobs, jobs_stats(1), jobs_stats(2), jobs_stats(3), length(dir(fullfile(p.queue.path, '*.dat')))); 
            end
            
            if use_gui
                chain = '|/-\';
                if msg_status == 0
                    fprintf([string, '%s\n'], chain(1+mod(msg_status,4))); 
                else
                    % delete previous message and update the last line of
                    % the text , note \b = backspace, \r (carriage return)
                    % does not work in Matlab command window 
                    fprintf([repmat('\b',1,length(string)+2),string, '%s\n'], chain(1+mod(msg_status,4))); 
                end
                msg_status = msg_status + 1; 
            end
        end
    end

end

function [log_size, log_size_0] = update_logs(file, tasks, log_files)
    %  write new information to the log files 
    log_size_0 = file.bytes; 
    diary = tasks.Diary; 
    fileID = fopen(log_files,'w');
    diary = strrep(diary,'%','%%');
    fprintf(fileID,diary);
    fclose(fileID);
    pause(0.1)
    file = dir(log_files);
    log_size = file.bytes; 
end

function report_outputs(log_file, job_msg)
    % print a nice report engine output 
    utils.verbose(0,'=========================================================================================')
    utils.verbose(0,job_msg)
    utils.verbose(0,'=========================================================================================')
    
    [~,cmdout] =system(['tail -n5 ', log_file,  '| grep "=== Failed ===" | wc -l']);
    if str2num(cmdout)
       Nlines = 80; 
    else
       Nlines = 20;  
    end
    [~,out] = system(sprintf('tail -n%i %s',Nlines, log_file));
    out = splitlines(out); 
    for line = out'
       fprintf('[report]    %s \n', line{1}); 
    end
    utils.verbose(0,'=========================================================================================')

end

function myCleanupFun(jobs, queue_path)
    %% end and delete all jobs -> cleanup prevents tasks running on GPU without control 
    for job = jobs
        cancel(job)
        delete(job)
    end
    
    % move all files from in_progress back to queue_path so that they can
    % be processed when the code runs again
    files = dir(fullfile(queue_path,'in_progress','*.dat')); 
    for ii = 1:length(files)
        io.movefile_fast(fullfile(queue_path,'in_progress', files(ii).name),fullfile(queue_path, files(ii).name))
    end
    
    
    warning on
    warning off backtrace
    warning('========================================')
    warning('All batch processed jobs were terminated')
    warning('========================================')
    warning on backtrace
end


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

