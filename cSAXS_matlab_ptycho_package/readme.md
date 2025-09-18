PTYCHOSHELVES
=======================================================================
version 20180816


Overview
-----------

Ptychoshelves [1] is a general and modular Matlab-based toolkit for ptychography that include an implementation of  maximum likelihood [2,3], difference map [4], and LSQ-ML [5] methods. 

The source code is available at [www.psi.ch/sls/csaxs/software](http://www.psi.ch/sls/csaxs/software)



Requirements
------------

The toolkit was tested with Matlab 2017b with CUDA 8 and Matlab 2018a-b with CUDA 9. Matlab version older than 2017b is not supported. 

Matlab parallel toolbox is required for GPU acceleration. Only Nvidia cards are supported. 

Fast CPU-implementation of ptychography engines is availible only for unix systems with modern CPU architectures supporting AVX2 instruction set. The code was tested only with Intel processors. 


Instalation
------------

Drag & drop both Matlab toolbox packages `cSAXS_matlab_base.mltbx` and `ptychoshelves.mltbx` into the matlab command window. 
Alternativelly, the packages can be installed using  following Matlab command 

`installedToolbox = matlab.addons.toolbox.installToolbox('MyToolbox.mltbx')`

Standard installation directory for Matlab toolboxes is: $HOME/Documents/{Add-Ons,MATLAB}

Troubleshooting 
---------------

The c_solver engine needs several common libraries and it was tested with the following versions: 
gcc/7.3.0
openmpi/3.0.1

Their presence can be checked by following command 

`$ ldd ~/Documents/MATLAB/Add-Ons/Toolboxes/ptychoshelves/code/+engines/+c_solver/ptycho_single_OMP `





Run our example code 
-------------------------

Use command `test_ptychoshelves` to run multiple tests that will check all engines and their methods 

Or edit our template  `edit template_ptycho.m` for your own reconstructions 



References
----------

[1] K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic data, *J. Appl. Cryst.* **53(2)** (2020).

[2] P. Thibault, M. Guizar-Sicairos, "Maximum-likelihood refinement for coherent diffractive imaging." New Journal of Physics **14.6** (2012): 063004.

[3] E. H. R. Tsai,  I. Usov, A. Diaz, A. Menzel, M. Guizar-Sicairos, "X-ray ptychography with extended depth of field." Optics express **24.25** (2016): 29089-29108.

[4] P. Thibault,  et al., *Science* **321**, 7 (2009)

[5] M. Odstrcil, A. Menzel, M. Guizar-Sicairos, *Opt. Express* **26**, 3108-3123 (2018)




>  Academic License Agreement
>================================
>
> Introduction 
> •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
>   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
>   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
>
> Terms and Conditions of the LICENSE
> 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
>       hereinafter set out and until termination of this license as set forth below.
> 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
>       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
>       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
> 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
>       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
>       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
>       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
> 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
>       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
>       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
>       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
>       in the commercial use, application or exploitation of works similar to the PROGRAM.
> 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
>       another computing language:
>       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
>       Scherrer Institut, Switzerland."
>
>    Additionally, any publication using the package, or any translation of the code into another computing language should cite:
> K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
> high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
> and for difference map:
> P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
>   (doi: 10.1126/science.1158573),
> for mixed coherent modes:
> P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
> for LSQ-ML method 
>  M. Odstrcil, A. Menzel, M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express, 3108-3123 (2018)  (2018). (doi: 10.1364/OE.26.003108).
> for OPRP method 
>  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation, 
>    Opt. Express 24, 8360 (2016). (doi: 10.1364/OE.24.008360).
> 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
>       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
> 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
>       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
>       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
>       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
> 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
> 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
>       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
>       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
>       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
> 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
>       the courts of Zürich, Switzerland. 
>
>
>

