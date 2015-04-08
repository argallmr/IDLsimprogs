; docformat = 'rst'
;
; NAME:
;    MrSim_ReconRate
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
function MrSim_ReconRate, theSim, time, $
CURRENT=current, $
GUESS=guess, $
T0=t0, $
VIEW_WIN=win, $
XRANGE=xrange, $
ZRANGE=zrange, $
_REF_EXTRA=extra

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
		void = cgErrorMSG()
		return, !Null
	endif
	
	current = keyword_set(current)
	
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
	osim_created = 0

	;If a simulation object was given, update the time step
	if isa(theSim, 'MrSim2') then begin
		oSim = theSim
		if n_elements(time) gt 0 then if oSim.TIME ne time then oSim.TIME = time
		
	;Create a simulation object (first iteration only)
	endif else if isa(theSim, /NUMBER) || isa(theSim, 'STRING') then begin
		oSim = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, _EXTRA=extra)
		osim_created = 1
	
	;Bad input
	endif else begin
		message, 'THESIM must be a simulation name, number or object reference.'
	endelse
	
;-------------------------------------------------------
; Find the X-Point /////////////////////////////////////
;-------------------------------------------------------
	;Get some properties. Store the initial time.
	oSim -> GetProperty, XRANGE=xr, ZRANGE=zr, TIME=t, NTIMES=nTimes
	if n_elements(t0) eq 0 then t0 = t
	
	;Create an initial guess
	if n_elements(guess) eq 0 then guess = mean([[xr], [zr]], DIMENSION=1)

	;Find the X-point
	xpt = MrSim_XPoint(oSim, GUESS=guess)
	dx  = (max(xr) - min(xr)) / 2
	dz  = (max(zr) - min(zr)) / 2

;-------------------------------------------------------
; Zoom-In on the X-Point ///////////////////////////////
;-------------------------------------------------------
	case oSim.coord_system of
		'MAGNETOPAUSE': begin
			xorder = [ 1, -1]
			zorder = [-1,  1]
		endcase
		'SIMULATION': begin
			xorder = [-1, 1]
			zorder = [-1, 1]
		endcase
		'MAGNETOTAIL': begin
			xorder = [ 1, -1]
			zorder = [-1,  1]
		endcase
		else: message, 'Invalid coordinate system "' + oSim.coord_system + '".'
	endcase
	
	;Zoom in to +/-6 de of the X-point
	oSim -> SetProperty, XRANGE=xpt[0] + dx*xorder, ZRANGE=xpt[1] + dz*zorder

;-------------------------------------------------------
; Calculate the Reconnection Rate //////////////////////
;-------------------------------------------------------
	;Get the electric field, x-, and z-coordinates
	Ey = oSim -> GetData('Ey')
	oSim -> GetProperty, XSIM=xSim, ZSIM=zSim
	
	;Find the region within 1di of the X-point
	roix = xpt[0] + 3.0 * xorder
	roiz = xpt[1] + 3.0 * zorder
	ix   = MrIndexRange(xSim, roix)
	iz   = MrIndexRange(zSim, roiz)
	
	;Take the average electric field over this area
	rate = mean(Ey[ix[0]:ix[1], iz[0]:iz[1]])

;-------------------------------------------------------
; Recurse for All Other Times //////////////////////////
;-------------------------------------------------------
	;Step backward to the beginning of the simulation.
	if (t le t0) && (t gt 0) then begin
		result = MrSim_ReconRate(oSim, t-1, T0=t0, GUESS=xpt)
		rate   = [result, rate]
	endif
	
	;Step forward to the end of the simulation.
	if (t ge t0) && (t lt nTimes-1) then begin
		;Must reset the range the first time
		if t eq t0 then oSim -> SetProperty, XRANGE=xr, ZRANGE=zr
		
		;Go forward in time
		result = MrSim_ReconRate(oSim, t+1, T0=t0, GUESS=xpt)
		rate   = [rate, result]
	endif

;-------------------------------------------------------
; Plot the Reconnection Rate ///////////////////////////
;-------------------------------------------------------
	;Reset to the original configuration
	if t eq t0 then oSim -> SetProperty, TIME=t0, XRANGE=xr, ZRANGE=zr
	return, rate
end
