; docformat = 'rst'
;
; NAME:
;    MrSim_XPoint
;
; PURPOSE:
;+
;   Find the X-point within a 2D simulation. X-points are located at saddle points of the
;   magnetic flux function Ay. To find an X-point, we sum the x- and z-derivatives of the
;   vector potential and search for the absolute minimum.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                               valid simulation object. See MrSim_Create.pro.
;       XRANGE:             in, required, type=fltarr(2)
;                           X-range in which to search for the X-point.
;       ZRANGE:             in, required, type=fltarr(2)
;                           Z-range in which to search for the X-point.
;       TIME:               in, required, type=int
;                           The simulation time at which to locate the X-point. Required if
;                               `THESIM` is not an object.
; :Keywords:
;       VIEW_WIN:           out, optional, type=object
;                           Set equal to a named variable into which a MrGraphics window
;                               object reference is returned. The window will contain
;                               contours of Ay within the specified range and an "X" at
;                               the location of the X-point.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted by MrSim_Create.pro
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
;   Modification History::
;       2015-03-16  -   Written by Matthew Argall
;-
function mrsim_xpoint, theSim, time, guess, $
SIM_OBJECT=oSim, $
VIEW_WIN=win, $
XRANGE=xrange, $
ZRANGE=zrange, $
_REF_EXTRA=extra

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
		if obj_valid(win) then obj_destroy, win
		void = cgErrorMSG()
		return, !Null
	endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
	osim_created = 0B

	;Simulation name or number?
	if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
		oSim = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, _STRICT_EXTRA=extra)
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
; Find the X-Point /////////////////////////////////////
;-------------------------------------------------------
	;Read the data and X- and Z-coordinates
	Ay = oSim -> GetData('Ay')
	oSim -> GetProperty, XSIM=xSim, ZSIM=zSim

	;Compute the derivative (centered difference)
	;  - dAy = ∂/∂x Ay dx + ∂/∂z Ay dz
	dAy_dx = oSim -> Scalar_Derivative('Ay', /D_DX)
	dAy_dz = oSim -> Scalar_Derivative('Ay', /D_DZ)
	
	;Trim off the edges
	dAy_dx = dAy_dx[1:-2, 1:-2]
	dAy_dz = dAy_dz[1:-2, 1:-2]
	Ay     = Ay[1:-2, 1:-2]
	xSim   = xSim[1:-2]
	zSim   = zSim[1:-2]
	
	;Locate the absolute minimum
;	if n_elements(guess) eq 0 then begin
		void = min(abs(dAy_dx) + abs(dAy_dz), ixpoint)
	
		;Determine the x- and z-index of the minimum
		inds = array_indices(dAy_dx, ixpoint)
		ix   = inds[0]
		iz   = inds[1]
	
;	endif else begin
;		;Initial guess
;		ix = value_locate(xSim, guess[0])
;		iz = value_locate(zSim, guess[1])
;	
;		count     = 0L
;		tolerance = 1e-4
;		max_iter  = 500
;		path = fltarr(2, max_iter)
;		
;		while count lt max_iter && dd gt tolerance do begin
;			void = max(dAy_dx[ix-1:ix+1, iz-1:iz+1], xMax, /ABSOLUTE) 
;			void = max(dAy_dz[ix-1:ix+1, iz-1:iz+1], zMax, /ABSOLUTE)
;			
;			
;			
;			;Magnitude of the gradient
;			mag = sqrt( dAy_dx[ix-1:ix+1, iz-1:iz+1]^2 + dAy_dz[ix-1:ix+1, iz-1:iz+1]^2 )
;			
;			;Move in the direction opposite to the steepest gradient
;			;  - Subtract 2 and take abs() to move to the opposite side of the square.
;			;  - Subtract 1 to put the middle cell at [0,0]
;			void = max(mag, iMax)
;			inds = array_indices([3,3], iMax, /DIMENSIONS)
;			inds = abs(inds - 2) - 1
;			
;			;Move the direction opposite to the steepest gradient.
;			ix += inds[0]
;			iz += inds[1]
;			
;			;Value of the derivative at this point
;			path[0, count] = [xSim[ix], zSim[iz]]
;			dd             = mag[inds[0], inds[1]]
;			count += 1
;		endwhile
;		
;		if count lt max_iter then path = path[*, 0:count-1]
;	endelse

;-------------------------------------------------------
; Mark the X-Point on Ay Contours //////////////////////
;-------------------------------------------------------
	;Pick out the contour levels
	if arg_present(win) then begin
		;Contour levels, including the separatrices
		levels = [Ay[ix, iz], cgConLevels(Ay, NLEVELS=nLevels)]
		levels = levels[sort(levels)]
		
		;Time
		twci  = oSim -> tIndex2txWci(oSim.TIME)
		title = 'Ay t$\Omega$$\downci$=' + strtrim(twci, 2)

		;Plot the contours
		if n_elements(win) gt 0 && obj_valid(win) then begin
			win            -> Refresh, /DISABLE
			win['Ay']      -> SetData, Ay, xSim, zSim
			win['Ay']      -> SetProperty, LEVELS=levels, TITLE=title
			win['X-Point'] -> SetProperty, XCOORDS=xSim[ix], YCOORDS=zSim[iz]
			win            -> Refresh
		endif else begin
			c_Ay = MrContour(Ay, xSim, zSim, $
			                 C_LABELS = 0, $
			                 LEVELS   = levels, $
			                 NAME     = 'Ay', $
			                 XTITLE   = 'x (de)', $
			                 YTITLE   = 'z (de)', $
			                 TITLE    = title)
			                 
			;Mark the X-point
			xpt = MrPlotS(xSim[ix], zSim[iz], /DATA, $
			             LINESTYLE = 'None', $
			             NAME      = 'X-Point', $
			             PSYM      = 7, $
			             SYMCOLOR  = 'Red', $
			             SYMSIZE   = 2, $
			             TARGET    = c_Ay, $
			             THICK     = 2)
		
			;Get the window
			win = c_Ay.window
		endelse
	endif

;-------------------------------------------------------
; Return ///////////////////////////////////////////////
;-------------------------------------------------------
	;Return the location of the X-point
	if osim_created then obj_destroy, oSim
	return, [xSim[ix], zSim[iz]]
end
