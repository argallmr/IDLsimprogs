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
function mrsim_xpoint, theSim, time, $
GUESS=guess, $
MAX_ITER=max_iter, $
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
; Derivatives of Ay ////////////////////////////////////
;-------------------------------------------------------
	;Read the data and X- and Z-coordinates
	Ay = oSim -> GetData('Ay')
	oSim -> GetProperty, XSIM=xSim, ZSIM=zSim

	;Compute the derivative (centered difference): B = Del x A
	;  - Bx   = - ∂Ay / ∂z
	;  - Bz   =   ∂Ay / dx
	;  - Bmag = sqrt( Bx^2 + By^2 )
	Bx = - oSim -> Scalar_Derivative('Ay', /D_DZ)
	Bz =   oSim -> Scalar_Derivative('Ay', /D_DX)
	
	;Trim off the edges
	Bx   = Bx[1:-2, 1:-2]
	Bz   = Bz[1:-2, 1:-2]
	Ay   = Ay[1:-2, 1:-2]
	xSim = xSim[1:-2]
	zSim = zSim[1:-2]

;-------------------------------------------------------
; Absolute Minimum /////////////////////////////////////
;-------------------------------------------------------
	if n_elements(guess) eq 0 then begin
		Bmag = sqrt(Bx^2 + Bz^2)
		!Null = min(Bmag, ixpoint)

		;Determine the x- and z-index of the minimum
		inds = array_indices(Bz, ixpoint)
		ix   = inds[0]
		iz   = inds[1]
		guess = [xSim[ix], zSim[iz]]
	endif

;-------------------------------------------------------
; Wander Toward X-Point ////////////////////////////////
;-------------------------------------------------------
	;Initial guess
	ix = value_locate(xSim, guess[0])
	iz = value_locate(zSim, guess[1])
	Bmag = sqrt(Bx^2 + Bz^2)
	Btemp = !values.f_infinity

	;Loop conditions
	if n_elements(max_iter) eq 0 then max_iter = 500
	count     = 0L
	ixTrans   = -2
	izTrans   = -2
	path      = fltarr(2, max_iter)
	index     = lonarr(2, max_iter)

	;Walk through the data grid
	while (count lt max_iter) do begin
		Btemp = Bmag[ix, iz]
		dx = 5
		dz = 10

		;Ascend in x -- Translate in the direction of steepest gradient
		BxSubArr = Bx[ix-dx:ix+dx, iz-dz:iz+dz]
		Bxmin    = min(BxSubArr, iNew, /ABSOLUTE)
		iNew     = array_indices([2*dx+1, 2*dz+1], iNew, /DIMENSIONS)
		izTrans = iNew[1] - dz

		;Descend in z -- Translate opposite to the steepest gradient
		BzSubArr = Bz[ix-dx:ix+dx, iz-dz:iz+dz]
		Bzmin    = min(BzSubArr, iNew, /ABSOLUTE)
		iNew     = array_indices([2*dx+1, 2*dz+1], iNew, /DIMENSIONS)
		ixTrans  = iNew[0] - dx

		;New  value of the derivative
		ix += ixTrans
		iz += izTrans
		
		;If elements repeat, we will trace the same path over and over.
		;   - This typically happens when we reach the x-point. We step away
		;     then circle back along the same trajectory.
		irepeat = where(ix eq index[0, *])
		if irepeat[0] ne -1 then irepeat = where(iz eq index[1, irepeat])
		if irepeat[0] ne -1 then begin
;				print, 'Repeated elements.'
			ix -= ixTrans
			iz -= izTrans
			break
		endif
		
		;Record the path and indices
		path[*, count]  = [xSim[ix], zSim[iz]]
		index[*, count] = [ix, iz]
		count += 1
	endwhile

	;Trim the path
	if count lt max_iter then path = path[*, 0:count-1]

;-------------------------------------------------------
; Plot Meander Toward X-Point //////////////////////////
;-------------------------------------------------------
	;Pick out the contour levels
	if arg_present(win) then begin
		;Contour levels, including the separatrices
		Ay_levels = [Ay[ix, iz],   cgConLevels(Ay, NLEVELS=nLevels)]
		B_levels  = [Bmag[ix, iz], cgConLevels(Bmag, NLEVELS=nlevels)]
		Bx_levels = [Bx[ix, iz],   cgConLevels(Bx, NLEVELS=nlevels)]
		Bz_levels = [Bz[ix, iz],   cgConLevels(Bz, NLEVELS=nlevels)]
		Ay_levels = Ay_levels[sort(Ay_levels)]
		B_levels  = B_levels[sort(B_levels)]
		Bx_levels = Bx_levels[sort(Bx_levels)]
		Bz_levels = Bz_levels[sort(Bz_levels)]
		
		;Time
		twci  = oSim -> tIndex2txWci(oSim.TIME)
		title = 'Ay t$\Omega$$\downci$=' + strtrim(twci, 2)

		;Plot the contours
		if n_elements(win) gt 0 && obj_valid(win) then begin
			win            -> Refresh, /DISABLE
			
			;Update data and paths
			win['Ay']      -> SetData, Ay, xSim, zSim
			win['B']       -> SetData, Bmag, xSim, zSim
			win['Bx']      -> SetData, Bx, xSim, zSim
			win['Bz']      -> SetData, Bz, xSim, zSim
			win['Ay']      -> SetProperty, LEVELS=Ay_levels, TITLE=title
			win['B']       -> SetProperty, LEVELS=B_levels
			win['Path_Ay'] -> SetProperty, XCOORDS=path[0,*],  YCOORDS=path[1,*]
			win['Path_B']  -> SetProperty, XCOORDS=path[0,*],  YCOORDS=path[1,*]
			win['Path_Bx'] -> SetProperty, XCOORDS=path[0,*],  YCOORDS=path[1,*]
			win['Path_Bz'] -> SetProperty, XCOORDS=path[0,*],  YCOORDS=path[1,*]
			win['X-Point'] -> SetProperty, XCOORDS=path[0,-1], YCOORDS=path[1,-1]
			win            -> Refresh
		endif else begin
			;Draw contours of Ay
			c_Ay = MrContour(Ay, xSim, zSim, $
			                 C_LABELS = 0, $
			                 LAYOUT   = [2,2,1], $
			                 LEVELS   = Ay_levels, $
			                 NAME     = 'Ay', $
			                 XTITLE   = 'x (de)', $
			                 YTITLE   = 'z (de)', $
			                 TITLE    = title)
			c_Ay -> Refresh, /DISABLE
			c_Ay.window -> SetProperty, XSIZE=750, YSIZE=600
			
			c_Bmag = MrContour(Bmag, xSim, zSim, /CURRENT, $
			                   C_LABELS = 0, $
			                   LAYOUT   = [2,2,3], $
			                   LEVELS   = B_levels, $
			                   NAME     = 'B', $
			                   XTITLE   = 'x (de)', $
			                   YTITLE   = 'z (de)', $
			                   TITLE    = '|B|')
			
			c_Bx = MrContour(Bx, xSim, zSim, /CURRENT, $
			                 C_LABELS = 0, $
			                 LAYOUT   = [2,2,2], $
			                 LEVELS   = Bx_levels, $
			                 NAME     = 'Bx', $
			                 XTITLE   = 'x (de)', $
			                 YTITLE   = 'z (de)', $
			                 TITLE    = 'Bx = -dAy/dz')
			
			c_Bz = MrContour(Bz, xSim, zSim, /CURRENT, $
			                   C_LABELS = 0, $
			                   LAYOUT   = [2,2,4], $
			                   LEVELS   = Bz_levels, $
			                   NAME     = 'Bz', $
			                   XTITLE   = 'x (de)', $
			                   YTITLE   = 'z (de)', $
			                   TITLE    = 'Bz = dAy/dx')
			                 
			;Mark the path
			gPath = MrPlotS(path[0,*], path[1,*], /DATA, $
			                LINESTYLE = 'None', $
			                NAME      = 'Path_Ay', $
			                PSYM      = 1, $
			                SYMCOLOR  = 'Blue', $
			                SYMSIZE   = 1, $
			                TARGET    = c_Ay, $
			                THICK     = 2)
			                 
			;Mark the path
			gPath = MrPlotS(path[0,*], path[1,*], /DATA, $
			                LINESTYLE = 'None', $
			                NAME      = 'Path_B', $
			                PSYM      = 1, $
			                SYMCOLOR  = 'Blue', $
			                SYMSIZE   = 1, $
			                TARGET    = c_Bmag, $
			                THICK     = 2)
			                 
			;Mark the path
			gPath = MrPlotS(path[0,*], path[1,*], /DATA, $
			                LINESTYLE = 'None', $
			                NAME      = 'Path_Bx', $
			                PSYM      = 1, $
			                SYMCOLOR  = 'Blue', $
			                SYMSIZE   = 1, $
			                TARGET    = c_Bx, $
			                THICK     = 2)
			                 
			;Mark the path
			gPath = MrPlotS(path[0,*], path[1,*], /DATA, $
			                LINESTYLE = 'None', $
			                NAME      = 'Path_Bz', $
			                PSYM      = 1, $
			                SYMCOLOR  = 'Blue', $
			                SYMSIZE   = 1, $
			                TARGET    = c_Bz, $
			                THICK     = 2)

			;Mark the X-Point
			gXpt = MrPlotS(path[0,-1], path[1,-1], /DATA, $
			               LINESTYLE = 'None', $
			               NAME      = 'X-Point', $
			               PSYM      = 1, $
			               SYMCOLOR  = 'Red', $
			               SYMSIZE   = 1, $
			               TARGET    = c_Ay, $
			               THICK     = 2)
		
			;Get the window
			win = c_Ay.window
			win -> Refresh
		endelse
	endif

;-------------------------------------------------------
; Return ///////////////////////////////////////////////
;-------------------------------------------------------
	;Return the location of the X-point
	if ~arg_present(osim) && osim_created then obj_destroy, oSim
	return, path[*,-1]
end
