% LOAD_FROM_P load parameters from the p-structure to param and self structures for GPU
% engine 
% 
% [self, param] = load_from_p(self, param, p)
% 
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** param       structure containing parameters for the engines 
% ** p          ptychoshelves p structure 
%
% returns: 
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ param       structure containing parameters for the engines 


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



function [self, param, p] = load_from_p(param, p)

    import math.*
    import utils.*
    import engines.GPU.shared.*
    
    
    [Np_p(1),Np_p(2),Npos] = size( p.fmag);
    self.reconstruct_ind = p.scanidxs;

    
    %% load default variables with different name from the main ptycho code 
    param.Nmodes = p.probe_modes; 
    param.plot_results_every = p.plot.interval; 
    
    % if defined in p, use from p otherwise use defaults in param 
    try;param.object_regular = p.reg_mu; end
    try;param.probe_inertia = p.probe_regularization; end

    if check_option(p,'opt_errmetric','poisson')
        param.likelihood = 'poisson'; 
    else
        param.likelihood = 'L1'; 
    end
    if get_option(p,'background_width') &&  get_option(p,'binning')
        param.background_width = param.background_width / 2^p.binning;
    end

        
    %% load variables from the main ptycho code and merge it with the defaults 
    for field = fieldnames(param)'
        field = field{1}; 
        if isfield(p, field)
            param.(field)  = p.(field);
        end
    end
          
    % set verbosity for GPU engine 
    param.verbose_level = max(-2,p.verbose_level-2); % adjust verbosity for GPU code , verbose_level 0 is enough for commmon use
    verbose(param.verbose_level)


    % load additional reconstructed parameters , otherwise use default 
    for item = {{'background',[]}, {'intensity_corr',[]}, {'probe_fourier_shift',[]}, {'rotation',0},{'shear',0},{'relative_pixel_scale',1}}
        item = item{1}; 
        if isfield(p, item{1}) && ~isempty(p.(item{1}))
            self.(item{1}) = p.(item{1});
        else
            self.(item{1}) = item{2}; 
        end
    end
    if any(ismember(fieldnames(p), {'shear', 'rotation', 'relative_pixel_scale'}))  && isfield(p, 'positions_0')
        warning('Reseting probe positions to original values')
        p.positions = p.positions_0; 
        p = rmfield(p, 'positions_0'); 
    end
    
    if isempty(p.affine_matrix)
       p.affine_matrix = diag([1,1]); 
    end
    
    self.diffraction_deform_matrix = [];

    if ~check_option(p,'asize_presolve') 
        param.Np_p_presolve = []; 
    else
        param.Np_p_presolve = min(p.asize_presolve, p.asize); 
    end
  
    % Other parameters 
    param.fourier_ptycho = check_option(p,'fourier_ptycho'); 
    if check_option(p, 'z_lens')
        self.z_lens = p.z_lens;
    end
    param.upsampling_data_factor = p.detector.upsampling; 
    
    % Offaxis ptychography correction 
    if check_option(p, 'sample_rotation_angles')
        param. sample_rotation_angles = p.sample_rotation_angles;  % 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)
    else
        param. sample_rotation_angles = [0,0,0]; % conventional ptychography
    end
    
    self.pixel_size = p.dx_spec;
    self.Np_p = Np_p; 
    Nscans = length(self.reconstruct_ind); 

    %%%%%%%%%%%%%%%%%%%%%%
    %% load probes%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    if param.share_probe
        p.share_probe_ID(:) = 1; % enforce single object if sharing is requested
    end
    
    % normalization for consistency with the other CPU engines 
    probes = single(p.probes ./ (prod(sqrt(Np_p))*2*p.renorm)); 

  
    for i = 1:min(p.probe_modes, size(probes,4))
        % variable probe 
        if isfield(p, 'probe_PCA') && ~isempty(p.probe_PCA) && i == 1   && size(p.probe_PCA.eigen_vec,1) == p.asize(1) && param.variable_probe && is_method(param, 'PIE')
            verbose(1,'Loading PCA probe (%i)', i)
            self.probe{i} = reshape(p.probe_PCA.eigen_vec,prod(p.asize),[]) * p.probe_PCA.evolution' /(prod(sqrt(Np_p))*2*p.renorm);% normalization for consistency with the CPU code;
        elseif  isfield(p, 'probe_variable') && ~isempty(p.probe_variable) && i == 1   && size(p.probe_variable.eigen_vec,1) == p.asize(1) && param.variable_probe && is_method(param, 'ML')
             % ML methods, OPRP approx
            verbose(0,'Loading variable probe (%i)', i)
            % constant part 
            self.probe{i}(:,:,:,1) = probes(:,:,:,1); 
            % variable part  
            self.probe{i}(:,:,:,2) = p.probe_variable.eigen_vec /(prod(sqrt(Np_p))*2*p.renorm);
            % evolution of the variable part 
            self.probe_evolution =  p.probe_variable.evolution ;
            assert(length(self.probe_evolution)==Npos, 'Wrong size of variable probe coefficients')
        else
            % constant probe 
            verbose(1,'Loading constant probe (%i)', i)
            if param.share_probe
                self.probe{i} = mean(probes(:,:,:,i),3); 
            else
                if size(p.probes,3) == Nscans 
                    % one probe for each scan 
                    self.probe{i} = probes(:,:,:,i);
                else
                    % rather take only first to avoid issues 
                    self.probe{i} = probes(:,:,1,i);
                end
            end
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% prepare support contraints%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % provide estimate of the probe support
    if check_option(p, 'probe_mask') && check_option(p,'use_probe_support') && any(p.probe_mask)
        self.probe_support = ~p.probe_mask; 
    elseif check_option(p,'probe_support_radius') && p.probe_support_radius < sqrt(2)
         % very useful for DM code 
         [X,Y] = meshgrid((-p.asize(2)/2+1:p.asize(2)/2)/p.asize(2), (-p.asize(1)/2+1:p.asize(1)/2)/p.asize(1));
         self.probe_support = sqrt(X.^2+Y.^2) < p.probe_support_radius/2; 
    else
        self.probe_support = [];
    end
   
    % estimate of the probe support in detector plane 
    if check_option(p,'probe_support_fft') &&  ~check_option(p, 'prop_regime', 'nearfield')
        if ~check_option(p.model, 'probe_focal_length') && ~check_option(p.model, 'probe_outer_zone_width')
            error('Missing  model.probe_focal_length and model.probe_outer_zone_width of Fresnel zone plate' )
        end      
        if check_option(p.model, 'probe_support_fft_type')
            param.probe_support_fft_type = p.model.probe_support_fft_type;
        end

        if strcmpi(param.probe_support_fft_type, 'circ')
            if ~check_option(p.model, 'probe_outer_zone_width')
                p.model.probe_outer_zone_width = p.lambda * p.model.probe_focal_length / p.model.probe_diameter;
            end
            FZP_cone_diameter = p.lambda* p.z/(p.model.probe_outer_zone_width * p.ds);
            % add some extra space
            FZP_cone_diameter = FZP_cone_diameter * 1.2;
            [X,Y] = meshgrid(-p.asize(2)/2+1:p.asize(2)/2, -p.asize(1)/2+1:p.asize(1)/2);
            self.probe_support_fft = utils.imgaussfilt2_fft(sqrt(X.^2+Y.^2) < FZP_cone_diameter/2, FZP_cone_diameter/50);
            af_probe = sqrt(abs(fftshift(fft2(self.probe{1}(:,:,1)))));
            [cx, cy] = center(max(0,af_probe-0.1*max(af_probe(:))));
            self.probe_support_fft = imshift_fast(self.probe_support_fft, -cx, -cy,[], 'nearest');
            self.probe_support_fft = max(0, min(1, self.probe_support_fft));
            verbose(1, 'Using farfield probe support constraint')
        elseif strcmpi(param.probe_support_fft_type, 'gaps')
            [m, my, mx, self.probe_support_fft] = core.get_detector_gaps(p.fmask);
            if size(m,1)==1 % no module
                self.probe_support_fft = [];
            else
                if mx > my
                    shift = [0, m(2,4)-m(1,2)];
                else
                    shift = [m(2,3)-m(1,1), 0];
                end
                self.probe_support_fft = fftshift(self.probe_support_fft);
                if check_option(p, 'probe_support_fft_float')
                    self.probe_support_fft_shifted = ~((abs(pshift(self.probe_support_fft, shift))-self.probe_support_fft)==-1);
                    self.probe_support_fft_relax = p.probe_support_fft_relax;
                end
            end
        elseif strcmpi(param.probe_support_fft_type, 'mask')
            [m, my, mx, self.probe_support_fft] = core.get_detector_gaps(p.fmask);
            if size(m,1)==1 % no module
                self.probe_support_fft = [];
            else
                if mx > my
                    shift = [0, m(2,4)-m(1,2)];
                else
                    shift = [m(2,3)-m(1,1), 0];
                end
                self.probe_support_fft = mean(p.fmask,3);
                self.probe_support_fft = fftshift(self.probe_support_fft);
                if check_option(p, 'probe_support_fft_float')
                    self.probe_support_fft_shifted = ~((abs(pshift(self.probe_support_fft, shift))-self.probe_support_fft)==-1);
                    self.probe_support_fft_relax = p.probe_support_fft_relax;
                end
            end
        else
            error('Unknown probe support type.');
        end
    else
        self.probe_support_fft = []; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    %% load object%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%
    
    if param.share_object
        p.share_object_ID(:) = 1; % enforce single object if sharing is requested
    end
    
    % correct positions for sample tilt 
    positions = p.positions ; 
    if isfield(p, 'positions_0')
        positions_0 = p.positions_0 ;
    end
    
    verbose(1, 'Load and pad object and probe')
    % update current size of the object 
    for i = 1:length(p.object)
        for j = 1:size(p.object{i},4) % load multiple layers of the object 
            self.object{i,j} = p.object{i}(:,:,1,j);
        end
        p.object_size(i,:) = size(self.object{i,1});
    end
    
    % calculate optimal size !! find minimal object to fit all scans 
    Np_o = max(p.object_size,[],1); 
    % expand object size if the probe p
    
    if p.number_iterations > check_option(p, 'probe_position_search') && is_method(param, {'ML', 'PIE'})
        extra = 0.2; % add plenty of extra space for geometry refinement 
    else
        extra = 0.05; % do just a little of extra space 
    end

    

    
    % shift the positions to account for the expanded object size AND
    % center them !!! (GPU code assumes positions to be centered, better for scale / probe positions are unknown)

    %%%%%%%%%%%%%%%%%%%%%%
    %% load positions %%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    for i = unique(p.share_object_ID)
        ind = [self.reconstruct_ind{p.share_object_ID == i}];
        position_offset = 1+floor((max(positions(ind,:))-min(positions(ind,:)))/2 + min(positions(ind,:)) );
        if isfield(p, 'positions_0')
            self.probe_positions_0(ind,:) = positions_0(ind,:) - position_offset;
            self.probe_positions(ind,:)   = positions(ind,:)   - position_offset;
        else
            self.probe_positions_0(ind,:) = positions(ind,:)   - position_offset;
            self.probe_positions = [];
        end
    end
    
    % get object extent 
    self.Np_o = max(Np_o,  ceil((1+extra) * ( self.Np_p +  (max(self.probe_positions_0) - min(self.probe_positions_0)) ))); 
    % store object size without padding, useful for plotting
    p.object_size = max(p.object_size, ceil(( self.Np_p +  (max(self.probe_positions_0) - min(self.probe_positions_0)) ))); 

    
    self.probe_positions_0 = self.probe_positions_0(:,[2,1]);
    if ~isempty(self.probe_positions)
        self.probe_positions = self.probe_positions(:,[2,1]);
    end

    self.Npos = Npos;
    
    % only a relative correction with respect to the affine matrix already
    % applied in p-struct 
    for ii = 1:Nscans
        self.affine_matrix{ii} = diag([1,1]); 
    end

    
    %%%%%%%%%%%%%%%%%%%%%%
    %% adjust object %%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    for i = 1:size(self.object,1)
        for layer = 1:size(self.object,2)
            % if needed expand the object to allow position refinement 
            % and shift for consistency with the CPU code
            self.object{i,layer} = imshift_fast(self.object{i,layer},1,1,self.Np_o, 'nearest', mean(self.object{i,layer}(:)));
        end
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load data, mask noise %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    verbose(1, 'Preparing data and masks')
    assert(all(isfinite(p.fmag(:))), 'Provide p.fmag contains NaN/Inf')
    self.noise  = [];
    self.diffraction = (single(p.fmag .* p.fmask) / single(p.renorm) ).^2;
    self.mask = (~p.fmask);
    
    
    if check_option(p, 'damped_mask')
        % if relaxed mask is used, try to fill the smallest gaps (hot pixels) by neighbors 
        mask_ind = find(self.mask); 
        self.diffraction(mask_ind) = self.diffraction(min(mask_ind+1, numel(self.diffraction))); 
    end

    if any(sum(sum(self.diffraction)) / prod(p.asize) < 0.01)
        error('%i patterns has average photon count < 0.01', mean(sum(sum(self.diffraction)) / prod(p.asize) < 0.01))
    end
    
    %% automatic data centering / flipping / tilted plane correction 
    if check_option(p, 'auto_center_data') || check_option(p, 'custom_data_flip') || check_option(p, 'sample_rotation_angles')
        self.diffraction = fftshift_2D(self.diffraction);
        self.mask = fftshift_2D(self.mask);
        warning on 
        warning off backtrace
        
        if check_option(p, 'auto_center_data')
            warning('Enforcing automatic data centering')
            for ii = 1:Nscans
                [x0,y0]=math.center(abs(fftshift(fft2(fftshift(self.probe{1}(:,:,min(end,ii)))))));
                W = mean2(self.diffraction(:,:,self.reconstruct_ind{ii}));
                W = ((W - min(W)) / (max(W)-min(W))).^4;  % give more weight to the more transpared regions (air)
                avg_pattern = mean(W.*sqrt(max(0,single(self.diffraction(:,:,self.reconstruct_ind{ii})))),3);
                [x,y]=math.center(avg_pattern);
                x = round(x-x0); y = round(y-y0); 
                self.diffraction(:,:,self.reconstruct_ind{ii}) = imshift_fast(self.diffraction(:,:,self.reconstruct_ind{ii}),x,y);
                self.mask(:,:,self.reconstruct_ind{ii}) = imshift_fast(self.mask(:,:,self.reconstruct_ind{ii}), x,y);
                fprintf('Data in scan %i shifted by %i %i pixels\n', ii, x,y); 
            end
            if ~isempty(self.probe_support_fft )
                self.probe_support_fft  = imshift_fast(self.probe_support_fft , x,y);
            end
        end
        
        % apply custom flip of the diffraction data 
        if check_option(p, 'custom_data_flip')  && any(p.custom_data_flip)
            warning('Applying custom data flip: %i %i %i ', p.custom_data_flip(1), p.custom_data_flip(2), p.custom_data_flip(3))
            if p.custom_data_flip(1)
                self.diffraction =  flipud(self.diffraction);
                self.mask        =  flipud(self.mask); 
            end
            if p.custom_data_flip(2)
                self.diffraction =  fliplr(self.diffraction); 
                self.mask        =  fliplr(self.mask); 
            end
            if p.custom_data_flip(3)
                self.diffraction =  permute(self.diffraction, [2,1,3]); 
                self.mask        =  permute(self.mask, [2,1,3]); 
            end
        end

        %     
         if isfield(p, 'sample_rotation_angles') &&  any(p.sample_rotation_angles) && check_option(p, 'apply_tilted_plane_correction', 'diffraction') 
            %% OFFAXIS PTYCHOGRAPHY CORRECTION 
            if utils.verbose > -1
                warning('Applying tilted plane correction: %3.5g %3.3g %3.3g\n Note that current implementation assumes low NA illumination, if this is not true, the central diffraction cone can be malformed', p.sample_rotation_angles(1), p.sample_rotation_angles(2), p.sample_rotation_angles(3))
            end
            % create matrix of deformation to apply effects similar to 
            deform_mat =  get_tilted_plane_correction_matrix(max(self.Np_p), ...
                                                p.z ,p.detectors{1}.pixel_size, ...
                                                p.sample_rotation_angles(1),...
                                                p.sample_rotation_angles(2),...
                                                p.sample_rotation_angles(3));
            if self.Np_p(1) ~= self.Np_p(2)
                % quick fix for asymmetric probe dimensions
                blank = true(self.Np_p); 
                blank_ind = find(utils.crop_pad(blank,[ max(self.Np_p), max(self.Np_p)])); 
                deform_mat = deform_mat(blank_ind,blank_ind); 
            end
                
            apply_deform = @(x,D)single(reshape(full(D * double(reshape(x,[],size(x,3)))), size(x)));
            plotting.smart_figure(3423)
            ax(1)=subplot(1,2,1);
            imagesc(log(1+max(self.diffraction, [],3))); axis off image xy ; colormap(plotting.colormaps.franzmap)
            title('Diffraction BEFORE tilted plane correction')
            % apply deformation effects caused by tilted sample , be sure to enforce the mask before
            % interpolation, the hotpixels can spread around after the correction
            self.diffraction = apply_deform(self.diffraction .* ~self.mask, deform_mat');
            self.mask =        apply_deform(self.mask, deform_mat') > 0;
            % plotting.imagesc3D(log(1+max(self.diffraction,[],3))); grid on  
            ax(2)=subplot(1,2,2);
            imagesc(log(1+max(self.diffraction, [],3))); axis off image xy ; colormap(plotting.colormaps.franzmap)
            title('Diffraction AFTER tilted plane correction')
            plotting.suptitle('Effect of tilted ptychography correction')
            linkaxes(ax, 'xy')
            drawnow
            self.diffraction_deform_matrix = deform_mat;
        end


        warning on 

        self.diffraction = ifftshift_2D(self.diffraction);
        self.mask = ifftshift_2D(self.mask);
    end
    
    
    self.filename = [p.detector.data_prefix, p.run_name]; 
    
    % other basic parameters 
    self.path = '';
    if check_option(p, 'prop_regime', 'nearfield')
        self.z_distance = p.z; 
    else
        self.z_distance = inf;  % farfield  
    end
    % multilayer extension 
    if isfield(p, 'delta_z')
        assert(all(isfinite(p.delta_z)), 'Some of the provided layer distanced delta_z is not finite' )
        self.z_distance = [p.delta_z(:)', self.z_distance]; 
    end
    
    self.lambda = p.lambda;
    self.diff_pattern_blur = 0;  % incoherent smoothing 
    self.modes = []; 
    
   
       
    % keep p structure for plotting purposes 
    param.p = p; 

    
    % if requested, shift the average probe to center and shift the object
    % to correspond to the probe shift
    if check_option(p, 'auto_center_probe')
        [x,y] = center(mean(abs(self.probe{1}(:,:,:,1)))); 
        for ii = 1:numel(self.probe)
            self.probe{ii} = imshift_fft(self.probe{ii}, -x,-y);
        end
        for ii = 1:numel(self.object)
            self.object{ii} = imshift_fft(self.object{ii}, -x, -y);
        end
    end
    
    % check if all is ok (remove in future !!!)
    positions = self.probe_positions_0(:,[2,1]); 
    positions = bsxfun(@plus, positions, ceil(self.Np_o/2-self.Np_p/2));    
    positions = round(positions); 
    range = ([min(positions(:,1)), max(positions(:,1))+ Np_p(1), min(positions(:,2)), max(positions(:,2))+ Np_p(2)]);
    if range(1) < 0 || range(2) > self.Np_o(1) || range(3) < 0 || range(4) > self.Np_o(2)
        warning('Object size is too small, not enough space for probes !! ') 
    end

    
end
