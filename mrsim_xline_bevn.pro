; docformat = 'rst'
;
; NAME:
;    MrSim_ProxX
;
;*****************************************************************************************
;   Copyright (c) 2015, Matthew Argall                                                   ;
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
;       TIME:               in, required, type=int
;                           The simulation time for which a velocity vector field
;                               is to be plotted. Required if `THESIM` is not an
;                               object.
; :Keywords:
;       SIM_OBJECT:         out, optional, type=object
;                           If `THESIM` is the name or number of a simulation, then
;                               this keyword returns the object reference to the
;                               corresponding simulation object that is created.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
;
; :Returns:
;       WIN:                MrWindow graphic window containing the requested graphics.
;                               If the image data does not exist, an invalid object
;                               will be returned.
;
; :Author:
;    Matthew Argall::
;    University of New Hampshire
;    Morse Hall Room 348
;    8 College Road
;    Durham, NH 03824
;    matthew.argall@unh.edu
;
; :History:
;    Modification History::
;       2015/05/25  -   Written by Matthew Argall
;-
function MrSim_XLine_BEVn, theSim, time, $
HORIZONTAL=horizontal, $
SIM_OBJECT=oSim, $
_REF_EXTRA=extra
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win)  then obj_destroy, win
		if obj_valid(win2) then obj_destroy, win2
		if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
	horizontal   = keyword_set(horizontal)
	osim_created = 0B

	;Simulation name or number?
	if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
		oSim = MrSim_Create(theSim, time, yslice, _STRICT_EXTRA=extra)
		if obj_valid(oSim) eq 0 then return, obj_new()
		osim_created = 1B
	
	;Object?
	endif else if MrIsA(theSim, 'OBJREF') then begin
		if obj_isa(theSim, 'MRSIM2') eq 0 $
			then message, 'THESIM must be a subclass of the MrSim class.' $
			else oSim = theSim
		
	;Somthing else
	endif else begin
		MrSim_Which
		message, 'THESIM must be a simulation name, number, or object.'
	endelse
	sim_class = obj_class(oSim)

;-------------------------------------------------------
; Determine XPoint & ReconRate \\\\\\\\\\\\\\\\\\\\\\\\\
;-------------------------------------------------------
	;First, get the time because computing the reconnection rate will change it.
	if n_elements(time) eq 0 then time = oSim.TIME
	
	;Create a window
	win = MrWindow(LAYOUT=[2,1], XGAP=25, XSIZE=900, YSIZE=250, OXMARGIN=[10,15], REFRESH=0)

	;Find the X-Point
	xpt  = MrSim_XPoint(oSim)

;-------------------------------------------------------
; 2D Image of Jey with X-Point Marked \\\\\\\\\\\\\\\\\\
;-------------------------------------------------------
	;Color image of Jey in a new window
	!Null = MrSim_ColorSlab(oSim, 'Jey', /CURRENT, C_NAME='Ay')

	;Update asthetics
	win['Ay Contours']   -> SetProperty, C_THICK=1.0, C_COLOR='Grey'
	win['CB: Color Jey'] -> SetProperty, WIDTH=0.5, XTICKS=2, XTICKFORMAT='(f0.2)'
	
	;Change title location?
	if oSim.coord_sys eq 'MAGNETOPAUSE' then begin
		title                 = win['Color Jey'].TITLE
		win['Color Jey']     -> SetProperty, TITLE=''
		win['CB: Color Jey'] -> SetProperty, TITLE=title
	endif

;-------------------------------------------------------
; Plot 1D Cuts of BL, ni, and EN through X-Point \\\\\\\
;-------------------------------------------------------
	coord = horizontal ? xpt[1] : xpt[0]

	;Plot cuts of BL, n, and EN through the X-line
	lc1 = MrSim_LineCut(oSim, 'BL', coord, /CURRENT, HORIZONTAL=horizontal)
	pos = lc1.position
	lc2 = MrSim_LineCut(oSim, 'ni',  coord, /CURRENT, HORIZONTAL=horizontal)
	lc3 = MrSim_LineCut(oSim, 'EN',  coord, /CURRENT, HORIZONTAL=horizontal)
	lc4 = MrSim_LineCut(oSim, 'UiL', coord, /CURRENT, HORIZONTAL=horizontal)
	
	;Move ni and Ex plots onto Bz to create independent data spaces for them
	lc2 -> SetProperty, COLOR='Blue', POSITION=pos, XSTYLE=5, YSTYLE=5, TITLE=''
	lc3 -> SetProperty, COLOR='Red',  POSITION=pos, XSTYLE=5, YSTYLE=5, TITLE=''
	lc4 -> SetProperty, COLOR='Forest Green', POSITION=pos, XSTYLE=5, YSTYLE=5, TITLE=''
	
	;Remove extra space created by moving graphics
	win -> TrimLayout
	
	;Reorder to put UiL being and BL in front
	lc1 -> Order, /BRING_TO_FRONT
	lc4 -> Order, /SEND_TO_BACK
	
	;Create axes for ni and Ex
	ax1 = MrAxis('Y', NAME='ni',  TARGET=lc2, LOCATION='Right', COLOR='Blue', TICKDIR=1, TITLE='n$\downi$')
	ax2 = MrAxis('Y', NAME='EN',  TARGET=lc3, LOCATION='Right', COLOR='Red',  TICKDIR=1, TITLE='E$\downN$', OFFSET=7)
	ax3 = MrAxis('Y', NAME='UiL', TARGET=lc4, LOCATION='Left',  COLOR='Forest Green',  TICKDIR=0, TITLE='U$\downiL$', OFFSET=7.5)
	
	;Find where Bz = 0
	lc1   -> GetData, x, Bz
	!Null  = min(Bz, iZero, /ABSOLUTE)
	
	;Draw a vertical line throubh BL = 0
	!Null = MrPlotS(x[[iZero, iZero]], lc1.yrange, $
	                COLOR     = 'Black', $
	                LINESTYLE = 'Dash', $
	                NAME      = 'BL=0 Vert', $
	                TARGET    = lc1)
	
	;Draw a horizontal line throubh BL = 0
	!Null = MrPlotS(lc1.xrange, [0, 0], $
	                COLOR     = 'Black', $
	                LINESTYLE = 'Dash', $
	                NAME      = 'BL=0 Horiz', $
	                TARGET    = lc1)
	
	;Draw a horizontal line throubh EN = 0
	!Null = MrPlotS(lc3.xrange, [0, 0], $
	                COLOR     = 'Red', $
	                LINESTYLE = 'Dash', $
	                NAME      = 'EN=0 Horiz', $
	                TARGET    = lc3)

;-------------------------------------------------------
; Return \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-------------------------------------------------------
	;Destroy the simulation object
	if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
	
	win -> Refresh
	return, win
end