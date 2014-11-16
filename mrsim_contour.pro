; docformat = 'rst'
;
; NAME:
;    MrSim_ColorSlab
;
;*****************************************************************************************
;   Copyright (c) 2013, Matthew Argall                                                   ;
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
;   The purpose of this program is to create a color image of a single data product with
;   the option of overlaying contours and vertical or horizontal lines.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                               valid simulation object. See MrSim_Create.pro.
;       NAME:               in, optional, type=string, default='Ay'
;                           The name of the GDA quantity to be contoured.
;       TIME:               in, required, type=int
;                           The simulation time for which a velocity vector field
;                               is to be plotted. Required if `THESIM` is not an
;                               object.
; :Keywords:
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
;       COLORWIN:           MrWindow graphic window containing the requested graphics.
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
;       2014/11/16  -   Written by Matthew Argall
;-
function MrSim_Contour, oSim, name, time, yslice, $
CURRENT = current, $
NLEVELS = nLevels, $
OFILENAME = ofilename, $
OVERPLOT = overplot, $
_REF_EXTRA = extra
	compile_opt strictarr

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		if obj_valid(c_plot) && keyword_set(current) eq 0 then obj_destroy, c_plot.window
		void = cgErrorMSG()
		return, obj_new()
	endif

;-------------------------------------------------------
; Defaults /////////////////////////////////////////////
;-------------------------------------------------------

	;Simulation object
	if size(oSim, /TNAME) ne 'OBJREF' && obj_isa(oSim, 'MRSIM') eq 0 $
		then message, 'OSIM must be an object and a subclass of the MrSim class.'
	
	current = keyword_set(current)
	if n_elements(name)       eq 0 then name      = 'Ay'
	if n_elements(ofilename)  eq 0 then ofilename = ''
	if n_elements(nlevels)    eq 0 then nlevels   = 15

	;Buffer the output?
	if current eq 0 then $
		if ofilename eq '' then buffer = 0 else buffer = 1

;-------------------------------------------------------
;Read Data /////////////////////////////////////////////
;-------------------------------------------------------
    ;Create the object and get the data. Use Ay for contours.
    data = oSim -> getData(name)

    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, XSIM=XSim, YSIM=YSim, ZSIM=ZSim, AXIS_LABELS=axLabls, $
                         COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame, ORIENTATION=orientation
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units

;-------------------------------------------------------
;Prepare Ranges ////////////////////////////////////////
;-------------------------------------------------------
	;Rename items
	_name = MrSim_Rename(name, coord_system, MVA_FRAME=mva_frame, /SUBSCRIPT)
	units = MrSim_Rename(units, /SUBSCRIPT)

	;Time interval in wci given?
	;   - Display time in terms of t*wci.
	;   - Display the time index.
	if n_elements(dtxwci) gt 0 $
		then title = 't$\Omega$$\downci$=' + string(time*dtxwci, FORMAT='(f0.1)') $
		else title = 't$\downindex$=' + string(time, FORMAT='(i0)')

	;If a 3D simulation was given, several orientations are possible
	case orientation of
		'XY': begin
			x = temporary(xSim)
			y = temporary(ySim)
			xtitle = axLabls[0] + ' (' + units + ')'
			ytitle = axLabls[1] + ' (' + units + ')'
		
			;Get ranges
			oSim -> GetProperty, XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange

			title  += '  ' + axLabls[2] + '=' + string(zrange[0], FORMAT='(f0.1)') + units
		endcase

		'XZ': begin
			x = temporary(xSim)
			y = temporary(zSim)
			xtitle = axLabls[0] + ' (' + units + ')'
			ytitle = axLabls[2] + ' (' + units + ')'
		
			;Get ranges
			oSim -> GetProperty, XRANGE=xrange, YRANGE=yr, ZRANGE=yrange
		
			if obj_class(oSim) eq 'MRSIM3D' $
				then title += '  ' + axLabls[1] + '=' + string(yr[0], FORMAT='(f0.1)') + units
		endcase

		'YZ': begin
			x = temporary(ySim)
			y = temporary(zSim)
			xtitle = axLabls[1] + ' (' + units + ')'
			ytitle = axLabls[2] + ' (' + units + ')'
		
			;Get ranges
			oSim -> GetProperty, YRANGE=xrange, ZRANGE=yrange
		endcase

		else: message, 'Orientation unknown: "' + orientation + '".'
	endcase

	;Make sure the data exists. Do this /after/ OSIM has a chance to be destroyed
	;and before the graphic window is created.
	if data eq !Null then return, obj_new()

;-------------------------------------------------------
;Add Contours? /////////////////////////////////////////
;-------------------------------------------------------
	;The X-point is a saddle point in Ay
	;   - Find the minimum along the bottom boundary.
	;   - Find the maximum along a vertical cut that passes through min point.
	if strupcase(name) eq 'AY' then begin
		nx = n_elements(x)
		ny = n_elements(y)

		;Take the derivative
		dx = data[1:nx-1, *] - data[0:nx-2, *]
		dy = data[*, 1:ny-1] - data[*, 0:ny-2]
	
		;Find the minimum of Ay along each row
		;   - Turn the 1D indices into 2D indices
		xMin = min(dx, ixMin, DIMENSION=1, /ABSOLUTE)
		inds = array_indices([nx-1, ny], ixMin, /DIMENSIONS)

		;Find the maximum along the path of minimums just found.
		yMin = min(dy[inds[0,*], inds[1,*]], iyMin, DIMENSION=2, /ABSOLUTE)
		sepAy  = data[inds[0, iyMin], inds[1, iyMin]]

		;Create the contour levels
		levels = [sepAy, cgConLevels(data, NLEVELS=nLevels)]
		levels = levels[sort(levels)]
	endif else begin
		levels = cgConLevels(data, NLEVELS=nLevels)
	endelse

	;Create the contour
	c_plot = MrContour(data, x, y, $
	                   CURRENT  = current, $
	                   C_LABELS = 0, $
	                   LEVELS   = levels, $
	                   NAME     = name + ' Contours', $
	                   OVERPLOT = overplot, $
	                   TITLE    = title, $
	                   XRANGE   = xrange, $
	                   XTITLE   = xtitle, $
	                   YRANGE   = yrange, $
	                   YTITLE   = ytitle, $
	                   _EXTRA   = extra)
;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------

	;Refresh and output, if requested.
	if current eq 0 then begin
		c_plot.window -> Refresh
		if ofilename ne '' then c_plot.window -> Save, ofilename
	endif

	return, c_plot
end
