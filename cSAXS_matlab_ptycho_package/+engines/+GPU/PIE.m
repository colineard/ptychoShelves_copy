% PIE - generalized version of the ptychographic iterative engine 
%
%[self, cache, fourier_error] =  PIE(self,par,cache,fourier_error,iter)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ** iter      current iteration number 
% returns:
% ++ self        self-like structure with final reconstruction
% ++ cache     structure with precalculated values to avoid unnecessary overhead
% ++ fourier_error  array [Npos,1] containing evolution of reconstruction error 
%
%
%   Publications most relevant to the Difference-Map implementation
%       Odstrcil, M., Baksh, P., Boden, S. A., Card, R., Chad, J. E., Frey, J. G., & Brocklesby, W. S
%       "Ptychographic coherent diffractive imaging with orthogonal probe relaxation."
%       Optics express 24.8 (2016): 8360-8369.

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


function [self, cache, fourier_error] =  PIE(self,par,cache,fourier_error,iter)
        import engines.GPU.GPU_wrapper.*
        import math.*
        import utils.*
        import engines.GPU.shared.*
        import engines.GPU.PIE.*
        
        %% MULTILAYER EXTENSION fix or remove in future
        par.probe_modes = max(par.Nlayers, par.probe_modes); 
        par.object_modes = max(par.Nlayers, par.object_modes); 
        
        
        if par.variable_probe
            % expand the probe to the full size
            probe_0 = self.probe{1}; 
            self.probe{1} = reshape(self.probe{1},prod(self.Np_p),[]);
            self.probe{1} = reshape(self.probe{1} * self.probe_evolution', self.Np_p(1), self.Np_p(2), []);
        end
        par.multilayer_object = par.Nlayers > 1; 
        par.multilayer_probe = false;  % not supported anymore 
        probe_amp_corr = [0,0]; 

        psi = cell(par.Nmodes, 1);
        Psi = cell(par.Nmodes, 1);
        probe_max = cell(par.probe_modes,1); 
        object_max = cell(par.object_modes,1);
        aprobe2 = abs(mean(self.probe{1},3)).^2; % update only once per iteration 
                
        if (is_method(par, 'ePIE') && ...
                ... % empirical estimation when the hybrid PIE method should be used 
                par.grouping > self.Npos/sqrt( pi^2 * mean(Ggather(cache.MAX_ILLUM)) / max(aprobe2(:)))) %  ||  ...% use hybrid ePIE in case of large grouping
%                 (~isempty(self.modes{end}.ASM_factor) && par.grouping > 1) || ...
%                 (is_method(par, 'ePIE')&&  par.multilayer_object)
            par.method = 'hPIE';  % hybrid ePIE
            if iter == 1;verbose(1,'Switching to hybrid PIE method '); end 
        end
        if is_method(par, {'hPIE'})
            for ll = 1:(par.object_modes*par.Nscans)
                for layer = 1:par.Nlayers
                    obj_illum_sum{ll,layer} = Gzeros(self.Np_o);
                    object_upd_sum{ll,layer} = Gzeros(self.Np_o, true);
                end
            end
            for ll = 1:par.probe_modes
                probe_upd_sum{ll}= (1+1i)*1e-8; 
                probe_illum_sum{ll} = 1e-8; 
            end
        end

        
        for ll = 1:par.Nlayers
            obj_proj{ll} = Gzeros([self.Np_p, par.grouping], true);
        end    
        if par.delta_p 
            grad = @get_grad_lsq;  % dumped least squares gradient (preconditioner)
        else
            grad = @get_grad_flat; % common PIE gradient 
        end

        %%%%%%%%%%%%%%%%%% ePIE algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        if par.grouping > 1
            for ll = 1:par.probe_modes
                if ~par.multilayer_object || ll == 1  % update only first layer probe 
                    probe_max{ll} = Ggather(max2(abs(self.probe{ll}))).^2;
                end
            end
            for ll = 1:par.object_modes
                if ~par.multilayer_probe || ll == 1  % update only first layer object 
                    object_max{ll} = Ggather(max2(abs(self.object{ll}(cache.object_ROI{:})))).^2;
                end
            end
        end
        
        if any(isnan(object_max{1})) ||  any(isnan(probe_max{1}))            
            error('Object or probe is nan, try smaller grouping')
        end
    

    
        % use already precalculated indices 
        rand_ind = randi(length( cache.preloaded_indices_sparse));
        indices = cache.preloaded_indices_sparse{rand_ind}.indices;
        scan_ids = cache.preloaded_indices_sparse{rand_ind}.scan_ids;

        for  ind_ii = randperm(length(indices))
            g_ind = indices{ind_ii}; 
            for ll = 1:max([par.probe_modes, par.object_modes,par.Nlayers]) 
                % generate indices of the used probes 
                if par.variable_probe && ll == 1
                    p_ind{ll} = g_ind;
%                 elseif par.multilayer_object && ll > 1
%                     p_ind{ll} = 1:length(ii);
                else  % single probe only
                     if par.share_probe %|| ll > 1  % share incoherent modes 
                        p_ind{ll} = 1;
                     else
                         if all(scan_ids{ind_ii} == scan_ids{ind_ii}(1))
                            p_ind{ll} =  scan_ids{ind_ii}(1);
                         else
                            p_ind{ll} =  scan_ids{ind_ii};
                         end
                     end
                end
            end

            %% load data to GPU (if not loaded yet)
            modF = get_modulus(self, cache, g_ind);
            mask = get_mask(self, cache.mask_indices, g_ind);
            noise = get_noise(self, par, g_ind); 
            

            % get objects projections 
            
            
            
            for layer = 1:par.Nlayers
                ll = 1; 
                obj_proj{layer} = get_views(self.object, obj_proj{layer},layer,ll,g_ind, cache, scan_ids{ind_ii},[]);
            end
            % get illumination probe 
            for ll = 1:par.probe_modes
                if ~par.multilayer_object || ismember(ll, [1,par.Nlayers:par.probe_modes])
                    probe{ll} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}));
                end
            end
                       
            for ll = 1:max([par.probe_modes, par.object_modes, par.Nlayers])  % Nlayers 
                %% fourier propagation  
               
                if (ll == 1 && (par.multilayer_object || par.multilayer_probe) )
                    probe{ll} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}));
                end
                if ll > 1 && par.multilayer_object ||  par.grouping == 1 % raw ePIE 
                    probe_max{ll} = max(Ggather(max2(abs(probe{ll})))).^2; 
                end
                if ll > 1 && par.multilayer_probe ||  par.grouping == 1 % raw ePIE 
                    object_max{ll} = max(Ggather(max2(abs(obj_proj{ll})))).^2;
                end
                if (ll == 1 && par.apply_subpix_shift) 
                    probe{ll}  = apply_subpx_shift(probe{ll}, self.modes{min(end,ll)}.sub_px_shift(g_ind,:) );
                end
                
                probe{ll} = apply_subpx_shift_fft(probe{ll}, self.modes{1}.probe_fourier_shift(g_ind,:)); 


                % get projection of the object and probe 
                psi{ll} = bsxfun(@times, obj_proj{min(ll,end)}, probe{min(ll,end)});
                Psi{ll} = fwd_fourier_proj(psi{ll} , self.modes{min(end, ll)}); 
                
                if ll < par.Nlayers && par.multilayer_object
                    probe{ll+1} =  Psi{ll};
                end
                
                if ll < par.Nlayers && par.multilayer_probe
                    obj_proj{ll+1} =  Psi{ll};
                end
            end
                        
            % get intensity (modulus) on detector including different corrections
            aPsi = get_reciprocal_model(self, Psi, modF, mask,iter, g_ind, par, cache);
                
                              
            %%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL CORRECTIONS END %%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%            
            if iter > 0 && par.get_error && (mod(iter,min(20, 2^(floor(2+iter/50)))) == 0 || (iter < 20)) || iter == par.number_iterations  % calculate only sometimes to make it faster
                [fourier_error(iter,g_ind)] = get_fourier_error(modF, aPsi, noise,mask, par.likelihood);
            end
            

        
            if (par.multilayer_object || par.multilayer_probe)
                [Psi(end),R] = modulus_constraint(modF,aPsi,Psi(end), mask, noise, par, 1); %  apply only on the last layer !!! 
            else
                [Psi,R] = modulus_constraint(modF,aPsi,Psi, mask, noise, par,1);
            end
            if iter == 0
                % in the first iteration only find optimal scale for the probe
                probe_amp_corr(1) = probe_amp_corr(1) + Ggather(sum(modF(:).^2));
                probe_amp_corr(2) = probe_amp_corr(2) + Ggather(sum(aPsi(:).^2));
                continue
            end
            if strcmp(par.likelihood, 'poisson')    % calculate only for the first mode
                %% automatically find  optimal step-size, note that for Gauss it is 1 !! 
                cache.beta_xi(g_ind)  =  gradient_descent_chi_solver(self,modF, aPsi2, R,mask, g_ind, mean(cache.beta_xi), cache);
            end
            
            if(par.multilayer_object || par.multilayer_probe)
                ind_modes = par.Nlayers:-1:1;
            else
                ind_modes = 1:max(par.probe_modes, par.object_modes);
            end
            

            for ll = ind_modes
                layer = ll; 
                chi = back_fourier_proj(Psi{min(end,ll)}, self.modes{min(end,ll)})-psi{min(end,ll)};
                               
                 
                %% get optimal gradient lenghts 
                 object_update=0; probe_update=0;m_probe_update= 0;
                
                 if iter >= par.object_change_start && (ll <= max(par.Nlayers, par.object_modes) || par.apply_multimodal_update)
                     object_update = Gfun(grad,chi, probe{min(ll,end)},...
                        probe_max{min(end,ll)}(1,1,min(end,p_ind{ll})),par.delta_p);
                 end
                
                
                if iter >= par.probe_change_start && ll <= max(par.probe_modes, par.Nlayers)
                    %% find optimal probe update !!! 
                        probe_update = Gfun(grad,chi,obj_proj{min(end,ll)},...
                            object_max{min(ll,end)}, par.delta_p);  
                       m_probe_update = mean(probe_update,3);
                end

                if ((ll == 1 && ~(par.multilayer_object || par.multilayer_probe)) || ...
                    (ll == par.Nlayers && (par.multilayer_object || par.multilayer_probe)))  && ...
                     (par.beta_LSQ ||  iter >= par.probe_position_search) 
                    %% variable step extension, apply only the first mode except the 3PIE case 
                    %% it will use the same alpha for the higher modes !!!!    
                                        
                    [cache.beta_probe(g_ind),cache.beta_object(g_ind)] =  gradient_projection_solver(self,chi,obj_proj{min(end,ll)},probe{ll},...
                        object_update, m_probe_update,p_ind{ll}, par, cache);
                     if any(g_ind ==1)
                         verbose(1,'Average alpha p:%3.3g  o:%3.3g ', mean(cache.beta_probe),mean(cache.beta_object));
                     end
                end
                     

                
                
                %%%%%%%%%%%%%%%%%%%%% PROBE UPDATE   %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%
                %% update probe first 
                if ((iter >= par.probe_change_start) && ll <= max(par.probe_modes,par.Nlayers) ...
                        &&  is_method(par, 'PIE')) || ...
                        (par.multilayer_object && ll > 1 && ll <= par.probe_modes)
                        % only in case of first layer probe, otherwise update interprobes
                    
                    beta_probe = get_vals(cache.beta_probe,g_ind) .* get_vals(cache.beta_xi,g_ind);                       
                            
                    if is_method(par, {'ePIE', 'hPIE'})
                        %%%%%%%%%%% update probe %%%%%%%%%%%%%%%%%%%%%%%%%%
                        probe{ll} = Gfun(@upd_probe_Gfun,probe{ll},probe_update, beta_probe); 

                        if (ll == 1 && par.apply_subpix_shift) 
                            probe{ll}  = apply_subpx_shift(probe{ll} , -self.modes{min(end,ll)}.sub_px_shift(g_ind,:));
                        end     
        
                     
                        if iter >= par.probe_change_start  
                            if (par.variable_probe && ll == 1)
                                self.probe{ll}(:,:,p_ind{ll}) = probe{ll};  % slowest line for large datasets !!!!!
                            elseif (~par.multilayer_object || ll == 1) && ll <= par.probe_modes
                                self.probe{ll} = mean(probe{ll},3); % merge information from all the shifted  probes if needed
                            end
                        end
                    end
                end
                             
                if iter > par.probe_fourier_shift_search && ll == 1
                    % search position corrections in the Fourier space, use
                    % only informatiom from the first mode, has to be after
                    % the probes updated 
                    self.modes{1} = gradient_fourier_position_solver(chi, obj_proj{1},probe{1},self.modes{1}, g_ind);
                end
                if  iter >= par.probe_position_search
                    % find optimal position shift that minimize chi{1} in current iteration 
                    [pos_update, cache] = gradient_position_solver(self, chi, obj_proj{1},probe{1,layer}, g_ind, iter, cache);
                    self.modes{1}.probe_positions(g_ind,:)=self.modes{1}.probe_positions(g_ind,:)+pos_update;
                end
                %%%%%%%%%%%%%%%%%%%%% OBJECT UPDATE   %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%
                if iter >= min([par.object_change_start]) && ....
                        ( ll <= max(par.Nlayers, par.object_modes) || par.apply_multimodal_update ) 
                    if ll ~= 1 && ~(par.multilayer_object || par.multilayer_probe) ; continue; end

                    if iter >= par.object_change_start  % && ~(par.multilayer_probe && ll > 1)  % the objects are just empty 
                    
                        beta_object =  get_vals(cache.beta_object,g_ind) .* get_vals(cache.beta_xi,g_ind);
                        if par.share_object
                            obj_ids = 1;  % update only the first object 
                        else
                            obj_ids = unique(scan_ids{ind_ii});  % update only the objects processed in this block 
                        end   
                        
                        if any(beta_object ~= 1)
                            object_update = bsxfun(@times, object_update, beta_object);
                        end
                                     

                        if  is_method(par, 'ePIE')  % use always in nearfield 
                            % classical ePIE, faster constraint application, but it will fail with too high grouping         
                            self.object = set_views(self.object, object_update,layer, obj_ids, g_ind, cache,  scan_ids{ind_ii},[]);
                        elseif is_method(par, 'hPIE')  %% hybrid PIE
                            if par.Nscans == 1 || par.share_object
                                ind_tmp = 1;
                            else
                                ind_tmp = 1+par.object_modes* ((1:par.Nscans)-1);
                            end
                            for kk = ind_tmp
                                obj_illum_sum{kk,layer}(:) = 0;
                                object_upd_sum{kk,layer}(:) = 0;
                            end
                            object_update = bsxfun(@times,  object_update,  aprobe2); % make is more like dumped LSQ solution 
                            [object_upd_sum,obj_illum_sum] = set_views_rc(object_upd_sum,obj_illum_sum, object_update,aprobe2,layer,obj_ids, g_ind, cache, scan_ids{ind_ii},[]);
                            for kk = ind_tmp
                                 self.object{kk,layer} = Gfun(@object_update_Gfun, self.object{kk,layer},object_upd_sum{kk,layer}, obj_illum_sum{kk,layer}, cache.MAX_ILLUM(min(kk,end)));
                            end
                        else
                            error('Unimplemented method %s ', par.method)
                        end
                    end                   
                end
                if ll > 1 && par.multilayer_object
%                     % apply rescaling to make scaling correction 
                    Psi{ll-1} =  probe{ll};
                end
                if ll > 1 && par.multilayer_probe
                    Psi{ll-1} =  obj_proj{ll} +  object_update;
                end  
            end
            
             if check_avail_memory < 0.2 || ~par.keep_on_gpu
                % slow step that is not needed if there is enough memory ,
                % it can slow down almost twice !!!
                 clear Psi R aPsi2 psi probe_update object_update chi 
                 if par.keep_on_gpu
                    warning('Low GPU memory')
                 end
             end
            
        end
       
        

       
       if par.multilayer_object
           self.probe = self.probe(1); 
       end
 
        if par.variable_probe && iter >= par.probe_change_start
            [self.probe{1}, self.probe_evolution] = apply_SVD_filter(self.probe{1}, par.variable_probe_modes+1,  self.modes{1});
        elseif par.variable_probe && iter < par.probe_change_start
            self.probe{1} = probe_0; 
        elseif ( ~isempty(self.probe_support)) &&  iter >= par.probe_change_start
            self.probe{1} = apply_probe_contraints(self.probe{1}, self.modes{1});
        end

       if iter == 0
           % apply initial correction for the probe intensity and return
           % it seems to be safer to underestimate the probe amplitude for variable probe method 
           probe_amp_corr = 0.5*sqrt(probe_amp_corr(1) / probe_amp_corr(2)); %% calculate ratio between modF^2 and aPsi^2
           
           for ii = 1:par.probe_modes
               self.probe{ii} = self.probe{ii}*probe_amp_corr;
           end
           verbose(2,'Probe amplitude corrected by %.3g',probe_amp_corr)
           return
       end
       
end

function probe = upd_probe_Gfun(probe,probe_update, alpha_p)
    probe =  probe + alpha_p.*probe_update;
end

function grad = get_grad_flat(chi,proj,max2, delta )
%   ePIE method
    grad =  (1/max2) * chi .* conj(proj) ;
end

function grad = get_grad_lsq(chi,proj,max2, delta )
    % dumped LSQ method  
    aproj = abs(proj) ;
    grad =  chi .*  conj(proj) .* aproj ./ (aproj.^2 + delta*max2)./sqrt(max2);
end

function object = object_update_Gfun(object,object_upd_sum, obj_illum_sum, max)
    object = object +  object_upd_sum ./ (obj_illum_sum+1e-9* max);
end

function array = get_vals(array, ind)
    if isscalar(array)
        return
    else
        array = reshape(array(ind),1,1,[]);
    end
end
