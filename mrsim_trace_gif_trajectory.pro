; docformat = 'rst'
;
; NAME:
;    MrSim_Trace_GIF_Trajectory
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
;   Create a GIF of traced particle trajectories.
;
; :Categories:
;    Bill Daughton, Simulation, tParticle
;
; :Params:
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                               valid simulation object. See MrSim_Create.pro.
;       FORWARD:            in, optional, type=string, default=''
;                           Name of a file in which forward-traced particle position and
;                               velocity data is stored.
;       BACKWARD:           in, optional, type=string, default=''
;                           Name of a file in which backward-traced particle position and
;                               velocity data is stored.
;
; :Keywords:
;       DELAY_TIME:         in, optional, type=integer, default=50
;                           Time delay between GIF frames, in multiples of 1/100.
;       C_NAME:             in, optional, type=string, default=''
;                           Name of the data product whose contours are to be 
;                               overplotted on the image.
;       IM_NAME:            in, optional, type=string, default=''
;                           Name of a GDA data product to be shown in color underneath
;                               the particle positions.
;       NLEVELS:            in, optional, type=int, default=15
;                           The number of contour lines to draw. Used with `C_NAME`.
;       OFILENAME:          out, optional, type=string, default=''
;                           If provided, graphics will be output to a file by this
;                               name. The file extension is used to determine which
;                               type of image to create. "PNG" is the default.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
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
pro MrSim_Trace_GIF_Trajectory, theSim, forward, backward, $
DELAY_TIME = delay_time, $
C_NAME = c_name, $
IM_NAME = im_name, $
NLEVELS = nLevels, $
OFILENAME = ofilename, $
XRANGE = xrange, $
YRANGE = yrange, $
ZRANGE = zrange, $
_REF_EXTRA = extra
	compile_opt strictarr

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if n_elements(thisDevice) gt 0 then set_plot, thisDevice
		if obj_valid(win) then begin
			win.SaveAs -> GIF_Close
			obj_Destroy, win
		endif
		void = cgErrorMSG()
		return
	endif

	;Defaults
	direction = 0
	if n_elements(backward)   eq 0 then backward   = ''
	if n_elements(delay_time) eq 0 then delay_time = 100
	if n_elements(forward)    eq 0 then forward    = ''
	
	;Output file
	if n_elements(ofilename) eq 0 then begin
		cd, CURRENT=pwd
		ofilename = filepath('MrSim_Trace_Trajectory.gif', ROOT_DIR=pwd)
	endif
	
;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------
	;Backward traced particles
	if file_test(backward) then begin
		direction += 1
		bdata      = MrSim_Trace_Read(backward)
	endif
	
	;Forward traced particles
	if file_test(forward) then begin
		direction += 2
		fdata      = MrSim_Trace_Read(forward)
	endif

	;Backward and forward together.
	case direction of
		0: message, 'A valid file for FORWARD or BACKWARD or both must be supplied.'
		1: data = reverse(temporary(bdata), 2)
		2: data = temporary(fdata)
		3: data = [[reverse((temporary(bdata))[*,1:*,*], 2)], [temporary(fdata)]]
	endcase
	
	;Define ranges
	if n_elements(xrange) eq 0 then xrange = [min(data[0,*,*], MAX=xMax), xMax]
	if n_elements(yrange) eq 0 then yrange = [min(data[1,*,*], MAX=yMax), yMax]
	if n_elements(zrange) eq 0 then zrange = [min(data[2,*,*], MAX=zMax), zMax]
	
;-------------------------------------------------------
; First Frame //////////////////////////////////////////
;-------------------------------------------------------
	;Plot the positions of each particle
	win = MrSim_Trace_Trajectory(theSim, reform(data[0:2,*,0]), $
	                             /BUFFER, $
	                             C_NAME     = c_name, $
	                             IM_NAME    = im_name, $
	                             NLEVELS    = nLevels, $
	                             XRANGE     = xrange, $
	                             YRANGE     = yrange, $
	                             ZRANGE     = zrange, $
	                             _REF_EXTRA = extra)
	
	;Switch to the Z-buffer
	thisDevice = !d.name
	set_plot, 'Z'
	
	;Open the GIF file
	win.SaveAs -> GIF_Open, ofilename

;-------------------------------------------------------
; Other frames /////////////////////////////////////////
;-------------------------------------------------------
	dims       = size(data, /DIMENSIONS)
	nParticles = dims[2]
	title      = win['Trajectory XY-plane'].TITLE

	;Update the positions and write a new frame.
	for i = 0, nParticles - 1 do begin
		if i eq 0 || i mod 5 eq 0 || i eq nParticles-1 then $
			print, FORMAT='(%"Saving frame %i of %i")', i+1, nParticles
	
		;Skip frames that do not have particles in them
		ix = where(data[0,i,*] gt xrange[0] and data[0,i,*] lt xrange[1], nx)
		iy = where(data[1,i,*] gt yrange[0] and data[1,i,*] lt yrange[1], ny)
		iz = where(data[2,i,*] gt zrange[0] and data[2,i,*] lt zrange[1], nz)
		if nx eq 0 && ny eq 0 && nz eq 0 then continue

		;Update particle positions
		win -> Refresh, /DISABLE
		win['Trajectory XY-plane'].TITLE = title + ' n=' + strtrim(i,2)
		win['XY Trajectory'] -> SetProperty, XCOORDS=reform(data[0, *, i]), YCOORDS=reform(data[1, *, i])
		win['XZ Trajectory'] -> SetProperty, XCOORDS=reform(data[0, *, i]), YCOORDS=reform(data[2, *, i])
		win['YZ Trajectory'] -> SetProperty, XCOORDS=reform(data[1, *, i]), YCOORDS=reform(data[2, *, i])
		win -> Refresh

		;Write frame
		win.SaveAs -> GIF_Write, DELAY_TIME=delay_time
	endfor

	;Close the GIF file
	win.SaveAs -> GIF_Close
	
	print, 'GIF saved to "' + oFileName + '".'
	set_plot, thisDevice
	obj_destroy, win
end
