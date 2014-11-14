; docformat = 'rst'
;
; NAME:
;       MrSim_tParticle
;
;*****************************************************************************************
;   Copyright (c) 2014, Matthew Argall                                                   ;
;   All rights reserved.                                                                 ;
;                                                                                        ;
;   Redistribution and use in source and binary forms, with or without modification,     ;
;   are permitted provided that the following conditions are met:                        ;
;                                                                                        ;
;       * Redistributions of source code must retain the above copyright notice,         ;
;         this list of conditions and the following disclaimer.                          ;
;       * Redistributions in binary form must reproduce the above copyright notice,      ;
;         this list of conditions and the following disclaimer in the documentation      ;
;         and/or other materials provided with the distribution.                         ;
;       * Neither the name of the <ORGANIZATION> nor the names of its contributors may   ;
;         be used to endorse or promote products derived from this software without      ;
;         specific prior written permission.                                             ;
;                                                                                        ;
;   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY  ;
;   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES ;
;   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT  ;
;   SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,       ;
;   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED ;
;   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR   ;
;   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     ;
;   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN   ;
;   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH  ;
;   DAMAGE.                                                                              ;
;*****************************************************************************************
;
; PURPOSE:
;+
;   Create input files for the "tparticle" simulation.
;
; :Property:
;       THESIM:         in, required, type=string/integer
;                       Name or number of the simulation for which tracing is to be
;                           performed.
;       TIME:           in, required, type=integer
;                       Simulation time index from which to take magnetic and electric
;                           field data.
;       BIN_CENTER:     in, required, type=fltarr(2)/fltarr(3)
;                       Location of the center of the spacial bin from which the
;                           the distribution function will be drawn. This is a 2- (3-)
;                           element array for 2D (3D) simulations
;       HALF_WIDTH:     in, required, type=fltarr(2)/fltarr(3)
;                       Half width of each special dimension of the box from which the
;                           distribution function is drawn.
;       VX_RANGE:       in, required, type=fltarr(2)
;                       Subset of Vx-space to be used for particle tracing.
;       VY_RANGE:       in, required, type=fltarr(2)
;                       Subset of Vy-space to be used for particle tracing.
;       VZ_RANGE:       in, required, type=fltarr(2)
;                       Subset of Vz-space to be used for particle tracing.
;
; :Keywords:
;       DIRECTION:      in, optional, type=int, default=3
;                       Direction which particles will be traced. Options are::
;                           1 - Forward tracing
;                           2 - Backward tracing
;                           3 - Both
;       DT:             in, optional, type=double, default=1/sqrt(2*!pi)
;                       Time interval at which to update particle position and velocity.
;       DX:             in, optional, type=double, default=from info file
;                       X-size of the domain in electron skin depths (de)
;       DY:             in, optional, type=double, default=from info file
;                       Y-size of the domain in electron skin depths (de)
;       DZ:             in, optional, type=double, default=from info file
;                       Z-size of the domain in electron skin depths (de)
;       EMASS:          in, optional, type=double, default=1.0
;                       Mass of an electron.
;       IMASS:          in, optional, type=double, default=from info file
;                       Mass of an ion.
;       N_ELECTRONS:    in, optional, type=double, default=100
;                       Maximum number of electrons to trace.
;       N_IONS:         in, optional, type=double, default=0
;                       Maximum number of ions to trace.
;       NX:             in, optional, type=double, default=from info file
;                       Number of grid cells in the x-dimension.
;       NY:             in, optional, type=double, default=from info file
;                       Number of grid cells in the y-dimension.
;       NZ:             in, optional, type=double, default=from info file
;                       Number of grid cells in the z-dimension.
;       OINITFILE:      in, optional, type=string, default="Init_<sim name>_t<`TIME`>_n<`N_ELECTRONS`>"
;                       Name of the file to which initial positions and velocities is written.
;       OFILENAME:      in, optional, type=string, default="Data_<sim name>_t<`TIME`>_n<`N_ELECTRONS`>"
;                       Name of the output data file. "_forward.dat" or "_backward.dat"
;       OPARAMFILE:     in, optional, type=string, default="Params_<sim name>_t<`TIME`>_n<`N_ELECTRONS`>"
;                       Name of the output file to which all parameters will be written.
;       STEPMAX:        in, optional, type=double, default=20000
;                       Maximum number of iterations to perform.
;       STEPPRINT:      in, optional, type=double, default=20
;                       Number of iterations between saves.
;       XMAX:           in, optional, type=double, default=`DX`/`NX`
;                       Number of (de) per grid cell in the x-direction.
;       YMAX:           in, optional, type=double, default=`DY`/`NY`
;                       Number of (de) per grid cell in the y-direction.
;       ZMAX:           in, optional, type=double, default=`DZ`/`NZ`
;                       Number of (de) per grid cell in the z-direction.
;
; :History:
;   Modification History::
;       2014-11-14  -   Written by Matthew Argall
;-
pro MrSim_Trace_Init, theSim, time, bin_center, half_width, vx_range, vy_range, vz_range, $
DIRECTION=direction, $
DT=dt, $
DX=dx, $
DY=dy, $
DZ=dz, $
EMASS=eMass, $
GDA_DIR=gda_dir, $
IMASS=iMass, $
N_ELECTRONS=n_electrons, $
N_IONS=n_ions, $
NX=nx, $
NY=ny, $
NZ=nz, $
OINITFILE=oInitFile, $
OPARAMFILE=oParamFile, $
OFILENAME=oFilename, $
STEPMAX=stepMax, $
STEPPRINT=stepPrint, $
XMAX=xmax, $
YMAX=ymax, $
ZMAX=zmax
	compile_opt idl2

	;catch errors
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		if n_elements(lun) gt 0 then free_lun, lun
		if obj_valid(oSim) then obj_destroy, oSim
		void = cgErrorMsg()
		return
	endif

	;Check parameters
	if n_params() eq 0 then message, 'Use: MrSim_tParticle_Init, theSim, time, bin_center, half_width, vx_range, vy_range, vz_range'
	if n_params() ne 7 then message, 'Incorrect number of parameters'

	;Defaults
	if n_elements(n_electrons) eq 0 then n_electrons = 100
	if n_elements(n_ions)      eq 0 then n_ions      = 0
	
	;Output files
	cd, CURRENT=dir
	basename = string(FORMAT='(%"_%s_t%i_n%i%s")', theSim, time, n_electrons, '.txt')
	if n_elements(oInitFile)   eq 0 then oInitFile   = filepath('Init'   + basename, ROOT_DIR=dir)
	if n_elements(oParamFile)  eq 0 then oParamFile  = filepath('Params' + basename, ROOT_DIR=dir)
	if n_elements(oFilename)   eq 0 then oFilename   = filepath('Data'   + basename, ROOT_DIR=dir)

;---------------------------------------------------------------------
; Get the Data ///////////////////////////////////////////////////////
;---------------------------------------------------------------------

	;Simulation object
	oSim   = MrSim_Create(theSim, time)

	;Convert y-location to grid cell
	ycell = oSim -> GetCell(bin_center[0], /Y)
	MrSim_Which, theSim, EFILE=eFile, DIST3D=dist3D, DIMENSION=dimension, DIRECTORY=_gda_dir, $
	             YSLICE=ycell, TINDEX=time

	;Read the data
	data = oSim -> ReadElectrons(eFile, $
	                             /VELOCITY, $
	                             DIST3D   = dist3D, $
	                             XRANGE   = bin_center[0] + [-half_width[0], half_width[0]], $
	                             YRANGE   = bin_center[1] + [-half_width[1], half_width[1]], $
	                             ZRANGE   = bin_center[2] + [-half_width[2], half_width[2]], $
	                             VX_RANGE = vx_range, $
	                             VY_RANGE = vy_range, $
	                             VZ_RANGE = vz_range)

	;Ensure a y-position is defined
	dims = size(data, /DIMENSIONS)
	if dims[0] eq 5 then data = [data[0,*], fltarr(1, dims[1]), data[1:4,*]]

	;Reduce the number of electrons
	if n_electrons lt dims[1] then begin
		if n_electrons eq 1 $
			then iElec = 0 $
			else iElec = long(findgen(n_electrons)/(n_electrons-1) * dims[1])
		data  = data[*,iElec]
	endif else if n_electrons gt dims[1] then begin
		n_electrons = dims[1]
	endif

;---------------------------------------------------------------------
; Setup //////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	;Get information about the simulation
	oSim -> GetInfo, MI_ME=_mi_me, NX=_nx, NY=_ny, NZ=_nz, DX_DE=_dx_de, DY_DE=_dy_de, DZ_DE=_dz_de
	obj_destroy, oSim
	
	;Default values
	is3D      = dimension eq '2D' ? 0 : 1
	n_e       = n_electrons
	n_i       = n_ions
	gda_dir   = n_elements(gda_dir)   eq 0 ? _gda_dir          : gda_dir
	dt        = n_elements(dt)        eq 0 ? 1.0 / sqrt(4*!pi) : dt
	eMass     = n_elements(eMass)     eq 0 ? 1.0               : eMass
	iMass     = n_elements(iMass)     eq 0 ? _mi_me            : iMass
	direction = n_elements(direction) eq 0 ? 3                 : direction
	nx        = n_elements(nx)        eq 0 ? _nx               : nx
	ny        = n_elements(ny)        eq 0 ? _ny               : ny
	nz        = n_elements(nz)        eq 0 ? _nz               : nz
	dx        = n_elements(dx)        eq 0 ? _dx_de            : dx
	dy        = n_elements(dy)        eq 0 ? _dy_de            : dy
	dz        = n_elements(dz)        eq 0 ? _dz_de            : dz
	xmax      = n_elements(xmax)      eq 0 ? nx * dx           : xmax
	ymax      = n_elements(ymax)      eq 0 ? ny * dy           : ymax
	zmax      = n_elements(zmax)      eq 0 ? nz * dz           : zmax
	stepMax   = n_elements(stepMax)   eq 0 ? 20000             : stepMax
	stepPrint = n_elements(stepPrint) eq 0 ? 20                : stepPrint

;---------------------------------------------------------------------
; Create Parameter File ///////////////////////////////////////////////
;---------------------------------------------------------------------

	;Open the output parameter file
	openw, lun, oParamFile, /GET_LUN

	;Inputs
	printf, lun, fix(gda_dir,    TYPE=7)    ;String
	printf, lun, fix(time,       TYPE=2)    ;Short
	printf, lun, fix(is3D,       TYPE=2)    ;Short
	printf, lun, fix(direction,  TYPE=2)    ;Short
	printf, lun, fix(stepMax,    TYPE=2)    ;Short
	printf, lun, fix(dt,         TYPE=5)    ;Double

	;Simulation parameters
	printf, lun, fix(n_e,        TYPE=3)    ;Short
	printf, lun, fix(n_i,        TYPE=2)    ;Short
	printf, lun, fix(eMass,      TYPE=5)    ;Double
	printf, lun, fix(iMass,      TYPE=5)    ;Double

	;Simulation Size
	printf, lun, fix(nx,         TYPE=3)    ;Long
	printf, lun, fix(ny,         TYPE=3)    ;Long
	printf, lun, fix(nz,         TYPE=3)    ;Long
	printf, lun, fix(dx,         TYPE=5)    ;Double
	printf, lun, fix(dy,         TYPE=5)    ;Double
	printf, lun, fix(dz,         TYPE=5)    ;Double
	printf, lun, fix(xmax,       TYPE=5)    ;Double
	printf, lun, fix(ymax,       TYPE=5)    ;Double
	printf, lun, fix(zmax,       TYPE=5)    ;Double

	;Output
	printf, lun, fix(stepPrint,  TYPE=2)    ;Short
	printf, lun, fix(oInitFile,  TYPE=7)    ;String
	printf, lun, fix(oFilename,  TYPE=7)    ;String

	;Close the file
	free_lun, lun

;---------------------------------------------------------------------
; Create Initial Conditions File /////////////////////////////////////
;---------------------------------------------------------------------

	;Particle data
	openw,    lun, oInitFile, /GET_LUN
	printf,   lun, data, FORMAT='( (3(f10.4, 3x), 2(f8.4, 3x), f8.4) )'
	free_lun, lun
	
	;Print where files were output to
	message, '', /INFORMATIONAL
	print, '    Parameter file output to:           "' + oParamFile + '".'
	print, '    Initial contiditons file output to: "' + oInitFile  + '".'
	print, '    Tracing data will be output to:     "' + oFilename  + '_<direction>.dat".'
end