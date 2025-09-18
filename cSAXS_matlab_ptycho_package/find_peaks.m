%FIND_PEAKS
% Routine to find diffraction peaks of a grating and compute its angle
% and the distance to the detector
% Be careful to avoid saturation of the center pixel as it will affect
% the calculaton
%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Change parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
clear
scannumber = 96;        % Scan number
gratingperiod = 0.4e-6; % Period of grating
Detector = 'eiger_4';   % 'pilatus_2' or 'eiger_4'

center = [568 810];     %[cy, cx]   Center of pattern
window = [-3:3];        % Window for computation of centroid [-3:3]
searchsteps = [0 35];    % As accurate as possible step search for the next peak
howmany = 20;            % Number of peaks is (howmany*2-1)
removecenter = true;    % If the center is saturated = true to remove the center from the fitting
symetricfitplot = false; % The plot in Fig.2 shows the same range in x and y
thresh = 10;             % Below this value peak is ignored
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%

addpath ~/Data10/cxs_software/base
addpath ./utils
import beamline.*
import io.*

[x y] = meshgrid(window,window);
threshold = round(1.5*sqrt(sum(searchsteps.^2)));     % Used in center computation, if the center between peaks exceeds threshold then it is not used in the computation

scanstring = sprintf('%05d',scannumber);

[data, mask] = io.nexus_read('~/Data10', 'scanNr', scannumber, 'filter', Detector, 'mask', sprintf('~/Data10/cxs_software/base/detector_masks/%s_valid_mask.h5', Detector));
%mask = io.nexus_read('./myFancyMask.h5', 'filter', 'mask');
detpix = data.x_pixel_size*1e-6;
specdata = io.nexus_read('~/Data10', 'scanNr', scannumber, 'filter','spec');%spec_read('~/Data10','ScanNr',scannumber);
lambda = 1.23984197e-9/specdata.mokev;        % Wavelength

%mask = load('~/Data10/cxs_software/ptycho/+detector/+eiger1p5M/eiger1p5m_mask.mat');

figure(1); imagesc(log(double(data.data(:,:,1)))); axis equal xy tight
data2 = sum(data.data,3) .* double(mask);
%data2 = sum(data.data,3); % To show without mask, not recommended for actual processing

figure(20);
imagesc(log10(data2))
axis xy equal tight

%% Save TIFF

clear centers
centers(1,1) = center(1) + sum(y.*data2(center(1)+[window],center(2)+[window]))/sum(data2(center(1)+[window],center(2)+[window]));
centers(1,2) = center(2) + sum(x.*data2(center(1)+[window],center(2)+[window]))/sum(data2(center(1)+[window],center(2)+[window]));

for ii = 2:howmany
    searchstepcounter = 1;
    try
        while data2(round(centers(ii-1,1)+searchstepcounter*searchsteps(1)),round(centers(ii-1,2)+searchstepcounter*searchsteps(2)) ) <= thresh
            searchstepcounter = searchstepcounter+1;
        end
        ceny = centers(ii-1,1)+searchstepcounter*searchsteps(1);
        cenx = centers(ii-1,2)+searchstepcounter*searchsteps(2);
        ROI = data2(round(ceny+window),round(cenx+window));
        centers(ii,1) = round(ceny) + sum(y(:).*ROI(:))/sum(ROI(:));
        centers(ii,2) = round(cenx) + sum(x(:).*ROI(:))/sum(ROI(:));
    catch
        warning(['I reached the end of the detector frame at howmany = ' num2str(ii)])
        howmany = ii-1;
        break
    end
end

tam = size(centers,1);

for jj = 2:howmany
    ii = jj+tam-1;
    searchstepcounter = -1;
    try
        if jj == 2
            while data2(round(centers(1,1)+searchstepcounter*searchsteps(1)),round(centers(1,2)+searchstepcounter*searchsteps(2)) ) <= thresh
                searchstepcounter = searchstepcounter-1;
            end
            ceny = centers(1,1)+searchstepcounter*searchsteps(1);
            cenx = centers(1,2)+searchstepcounter*searchsteps(2);
            ROI = data2(round(ceny+window),round(cenx+window));
            centers(ii,1) = round(ceny) + sum(y(:).*ROI(:))/sum(ROI(:));
            centers(ii,2) = round(cenx) + sum(x(:).*ROI(:))/sum(ROI(:));
        else
            while data2(round(centers(ii-1,1)+searchstepcounter*searchsteps(1)),round(centers(ii-1,2)+searchstepcounter*searchsteps(2)) ) <= thresh
                searchstepcounter = searchstepcounter-1;
            end
            ceny = centers(ii-1,1)+searchstepcounter*searchsteps(1);
            cenx = centers(ii-1,2)+searchstepcounter*searchsteps(2);
            ROI = data2(round(ceny+window),round(cenx+window));
            centers(ii,1) = round(ceny) + sum(y(:).*ROI(:))/sum(ROI(:));
            centers(ii,2) = round(cenx) + sum(x(:).*ROI(:))/sum(ROI(:));
        end
    catch
        warning(['I reached the end of the detector frame at howmany = ' num2str(jj)])
        howmany = jj-1;
        break
    end
end

if removecenter
    whichdisplay = [2:size(centers,1)];
    display('Removing center pixel from the fit')
else
    whichdisplay = [1:size(centers,1)];
end
c = polyfit(centers(whichdisplay,2),centers(whichdisplay,1),1);


figure(1);
set(gcf,'outerPosition',[1   294   712   731])
imagesc(log10(data2))
axis xy equal tight
hold on
plot(centers(whichdisplay,2),centers(whichdisplay,1),'+w')
plot(centers(whichdisplay,2),centers(whichdisplay,1),'or')
hold off
fov = max(max(centers(:,2))-min(centers(:,2)),max(centers(:,1))-min(centers(:,1)))*1.5;
axis([center(2)-fov/2 center(2)+fov/2 center(1)-fov/2 center(1)+fov/2])
%
figure(2)
set(gcf,'outerPosition',[713   521   568   504])
plot(centers(whichdisplay,2),centers(whichdisplay,1),'.r')
hold on
xplot = linspace(min(centers(:,2))*0.999,max(centers(:,2))*1.001,100);
yplot = c(1)*xplot+c(2);%linspace(min(centers(:,1)),max(centers(:,1)),100);
plot(xplot,yplot,'--')
hold off
if symetricfitplot
    axis([center(2)-fov/2 center(2)+fov/2 center(1)-fov/2 center(1)+fov/2])
    axis equal
else
    axis([min(centers(:,2)) max(centers(:,2)) min(centers(:,1)) max(centers(:,1))]);
end
title('Fitting to centers')


alpha = atan(1/c(1))*180/pi;
display(['alpha = ' num2str(alpha) ' degrees'])
%%%
% Find distance to detector
clear deltapoints

if ~removecenter
    for ii = 2:howmany
        deltapoints(ii-1) = sqrt((centers(ii,1)-centers(ii-1,1))^2+(centers(ii,2)-centers(ii-1,2))^2);
    end
    ii = howmany+1;
    deltapoints(ii-1) = sqrt((centers(ii,1)-centers(1,1))^2+(centers(ii,2)-centers(1,2))^2);
    for ii = tam+2:howmany-1+tam
        deltapoints(ii-1) = sqrt((centers(ii,1)-centers(ii-1,1))^2+(centers(ii,2)-centers(ii-1,2))^2);
    end
else
    for ii = 3:howmany
        deltapoints(ii-1) = sqrt((centers(ii,1)-centers(ii-1,1))^2+(centers(ii,2)-centers(ii-1,2))^2);
    end
    for ii = tam+2:howmany-1+tam
        deltapoints(ii-1) = sqrt((centers(ii,1)-centers(ii-1,1))^2+(centers(ii,2)-centers(ii-1,2))^2);
    end
    deltapoints = deltapoints(deltapoints~=0);
end

deltapoints = deltapoints(deltapoints<threshold);

display(['Detector distance ' num2str(mean(deltapoints)*detpix*gratingperiod/lambda) ' m']),

figure(3);
set(gcf,'outerPosition',[713    29   568   492])
plot(deltapoints*detpix*gratingperiod/lambda)
hold on
plot([1:length(deltapoints)],ones(size(deltapoints))*mean(deltapoints)*detpix*gratingperiod/lambda,'--r')
hold off
title(['Detector distance per pixel separation'])

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
