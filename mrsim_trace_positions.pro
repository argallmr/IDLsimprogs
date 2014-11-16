; docformat = 'rst'
;
; NAME:
;    MrSim_Trace_Positions
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
;       * Neither the name of the University of New Hampshire nor the names of its       ;
;         contributors may be used to endorse or promote products derived from this      ;
;         software without specific prior written permission.                            ;
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
;   Plot the positions of a set of particles.
;
; :Categories:
;    Bill Daughton, Simulation, tParticle
;
; :Params:
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                               valid simulation object. See MrSim_Create.pro.
;       POSITION:           in, required, type=3xN fltarr
;                           Position of the particle for which to plot the trajectory.
;                               N is the number of particles to be plotted.
;       VELOCITY:           in, optional, type=3xN fltarr
;                           Velocity of the particle for which the positions are plotted.
;                               If provided, trajectory will be color-coated with the
;                               velocity magnitude.
;
; :Keywords:
;       C_NAME:             in, optional, type=string, default=''
;                           Name of the data product whose contours are to be 
;                               overplotted on the image.
;       CURRENT:            in, optional, type=boolean, default=0
;                           If set, the graphics will be added to the current MrWindow
;                               graphics window. The default is to create a new window.
;       NLEVELS:            in, optional, type=int, default=15
;                           The number of contour lines to draw. Used with `C_NAME`.
;       OFILENAME:          out, optional, type=string, default=''
;                           If provided, graphics will be output to a file by this
;                               name. The file extension is used to determine which
;                               type of image to create. "PNG" is the default.
;       SIM_OBJECT:         out, optional, type=object
;                           If `THESIM` is the name or number of a simulation, then
;                               this keyword returns the object reference to the
;                               corresponding simulation object that is created.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
;
; :Returns:
;       TWIN:               MrWindow graphic window containing the requested graphics.
;                               If the image data does not exist, an invalid object
;                               will be returned.
;
; :Author:
;    Matthew Argall::
;    University of New Hampshire
;    Morse Hall Room 113
;    8 College Road
;    Durham, NH 03824
;    matthew.argall@wildcats.unh.edu
;
; :History:
;    Modification History::
;       2014/11/13  -   Written by Matthew Argall
;-
function MrSim_Trace_Positions, theSim, position, velocity, $
BUFFER = buffer, $
C_NAME = c_name, $
CURRENT = current, $
IM_NAME = im_name, $
NLEVELS = nLevels, $
OFILENAME = ofilename, $
SIM_OBJECT = oSim, $
XRANGE = xrange, $
XY_PLANE = xy_plane, $
XZ_PLANE = xz_plane, $
YRANGE = yrange, $
YZ_PLANE = yz_plane, $
ZRANGE = zrange, $
_REF_EXTRA = extra
	compile_opt strictarr

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		if osim_created        && arg_present(oSim)    eq 0 then obj_destroy, oSim
		if obj_valid(colorWin) && keyword_set(current) eq 0 then obj_destroy, colorWin
		void = cgErrorMSG()
		return, obj_new()
	endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
	osim_created = 0B

	;Simulation name or number?
	if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
		oSim = MrSim_Create(theSim, time, yslice, _STRICT_EXTRA=extra)
		if obj_valid(oSim) eq 0 then return, obj_new()
		osim_created = 1B
	
	;Object?
	endif else if MrIsA(theSim, 'OBJREF') then begin
		if obj_isa(theSim, 'MRSIM') eq 0 $
			then message, 'THESIM must be a subclass of the MrSim class.' $
			else oSim = theSim
		
	;Somthing else
	endif else begin
		MrSim_Which
		message, 'THESIM must be a simulation name, number, or object.'
	endelse
	sim_class = obj_class(oSim)

;-------------------------------------------------------
; Defaults /////////////////////////////////////////////
;-------------------------------------------------------
	buffer   = keyword_set(buffer)
	current  = keyword_set(current)
	xy_plane = keyword_set(xy_plane)
	xz_plane = keyword_set(xz_plane)
	yz_plane = keyword_set(yz_plane)
	if n_elements(c_name)    eq 0 then c_name    = ''
	if n_elements(im_name)   eq 0 then im_name   = ''
	if n_elements(ofilename) eq 0 then ofilename = ''
	if n_elements(nlevels)   eq 0 then nlevels   = 15
	if n_elements(xrange)    eq 0 then xrange    = [min(position[0,*], MAX=xMax), xMax]
	if n_elements(yrange)    eq 0 then yrange    = [min(position[1,*], MAX=yMax), yMax]
	if n_elements(zrange)    eq 0 then zrange    = [min(position[2,*], MAX=zMax), zMax]

	;If none were chosen, pick all
	if xy_plane + xz_plane + yz_plane eq 0 then begin
		xy_plane = 1
		xz_plane = 1
		yz_plane = 1
	endif

	;Buffer the output?
	if current eq 0 then $
		if ofilename ne '' then buffer = 1

;-------------------------------------------------------
;Prepare Ranges ////////////////////////////////////////
;-------------------------------------------------------
	;Get the simulation size and time
	oSim -> GetProperty, TIME=time
	oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units

	;Rename items
	;    _name = MrSim_Rename(name, coord_system, MVA_FRAME=mva_frame, /SUBSCRIPT)
	units = MrSim_Rename(units, /SUBSCRIPT)

	;Time interval in wci given?
	;   - Display time in terms of t*wci.
	;   - Display the time index.
	if n_elements(dtxwci) gt 0 $
		then title = 't$\Omega$$\downci$=' + string(time*dtxwci, FORMAT='(f0.1)') $
		else title = 't$\downindex$=' + string(time, FORMAT='(i0)')
	
	;Destroy the object
	if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim

;-------------------------------------------------------
;Make Plots ////////////////////////////////////////////
;-------------------------------------------------------
    ;Convert velocity to |v|
;    if n_elements(velocity) gt 0 then begin
;        vmag = sqrt(total(velocity^2, 1))
;        vmag = bytscl(vmag)
;    endif

	;Create a window
	if current $
		then colorWin = GetMrWindows(/CURRENT) $
		else colorWin = MrWindow(XSIZE=600, YSIZE=500, YGAP=4, BUFFER=buffer, REFRESH=0)

	;Plot in XY-Plane
	if xy_plane then begin
		xyPlot = MrPlot([0], [0], /CURRENT, $
		                /NODATA, $
		                NAME   ='Positions XY-plane', $
		                XRANGE = xrange, $
		                XTITLE = 'X (' + units + ')', $
		                YRANGE = yrange, $
		                YTITLE = 'Y (' + units + ')')

		;Plot the trajectory
		xyPath = MrPlotS(position[0,*], position[1,*], $
	                 COLOR  = vmag, $
	                 NAME   = 'XY Positions', $
	                 NOCLIP = 0, $
	                 PSYM   = 2, $
	                 TARGET = xyPlot)
	endif

	;Plot in XZ-Plane
	if xz_plane then begin
		;Contour plot
		if c_name ne '' then begin
			xzPlot = MrSim_Contour(oSim, c_name, /CURRENT, $
			                       NAME   = 'Positions XZ-plane', $
			                       TITLE  = '', $
			                       XRANGE = xrange, $
			                       XTITLE = 'X (' + units + ')', $
			                       YRANGE = zrange, $
			                       YTITLE = 'Z (' + units + ')')
	
		;Create a set of axes
		endif else begin
			xzPlot = MrPlot([0], [0], /CURRENT, $
			                /NODATA, $
			                NAME   = 'Positions XZ-plane', $
			                XRANGE = xrange, $
			                XTITLE = 'X (' + units + ')', $
			                YRANGE = zrange, $
			                YTITLE = 'Z (' + units + ')')
		endelse

		;Plot the trajectory
		xzPath = MrPlotS(position[0,*], position[2,*], $
		                 COLOR  = vmag, $
		                 NAME   = 'XZ Positions', $
		                 NOCLIP = 0, $
		                 PSYM   = 2, $
		                 TARGET = xzPlot)
	endif

	;Plot in YZ-Plane
	if yz_plane then begin
		;Create a set of axes
		yzPlot = MrPlot([0], [0], /CURRENT, $
		                /NODATA, $
		                NAME   = 'Positions YZ-plane', $
		                XRANGE = yrange, $
		                XTITLE = 'Y (' + units + ')', $
		                YRANGE = zrange, $
		                YTITLE = 'Z (' + units + ')')

		;Plot the trajectory
		yzPath = MrPlotS(position[1,*], position[2,*], $
		                 COLOR  = vmag, $
		                 NAME   = 'YZ Positions', $
		                 NOCLIP = 0, $
		                 PSYM   = 2, $
		                 TARGET = yzPlot)
	endif

	;Put a title on the plot    
	case 1 of
		xy_plane: xyPlot.TITLE = title
		xz_plane: xyPlot.TITLE = title
		yz_plane: xyPlot.TITLE = title
	endcase

	;Refresh and output, if requested.
	if current eq 0 then begin
		colorWin -> Refresh
		if ofilename ne '' then colorWin -> Save, ofilename
	endif

	return, colorWin
end
