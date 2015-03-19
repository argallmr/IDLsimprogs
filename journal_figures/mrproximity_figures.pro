; docformat = 'rst'
;
; NAME:
;    MrProximity_Figures
;
; PURPOSE:
;+
;   Create figures for my X-Line proximity paper.
;
; :Categories:
;    Bill Daughton, Simulation
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
;       2015/01/22  -   Written by Matthew Argall
;       2015/03/03  -   Draw distribution boxes on Figure 1. Remove images from Figure 2. - MRA
;-
;*****************************************************************************************
;+
;   Create Figure 1: X-Line proximity for Asymm-Scan/By0, Asymm-Scan/By1, and Asymm-3D.
;-
function MrProxFigs_Figure1, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win1)   then obj_destroy, win1
		if obj_valid(win2)   then obj_destroy, win2
		if obj_valid(win3)   then obj_destroy, win3
		if obj_valid(theSim) then obj_destroy, theSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0'
	cuts      = [36.77, 34.9, 33.0]
	tIndex    = 28
	xrange    = [2,-2]
	zrange    = 36.77 + [-5, 5]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	im_name   = 'Jey'
	c_name    = 'Ay'

;	cut_range = [2,-2]
	xrange = [2,-2]
;	xrange    = [2.0,-1.5]
;	zrange    = 36.77 + [-4.0, 4.0]
	zrange    = 36.77 + [-5, 5]

	win1 = MrSim_XProximity(theSim, cuts, tIndex, $
	                        C_NAME     = c_name, $
	                        COORD_SYS  = coord_sys, $
	                        IM_NAME    = im_name, $
	                        ION_SCALE  = ion_scale, $
	                        MVA_FRAME  = mva_frame, $
	                        SIM_OBJECT = oSim, $
	                        XRANGE     = xrange, $
	                        ZRANGE     = zrange)

;-----------------------------------------------------
; Asymm-Scan/By0: Mark Distribution Functions \\\\\\\\
;-----------------------------------------------------
	;Distribution information (de)
	dist_width  = [0.5, 0.5]
	dist_loc    = [ $
	               ;1st Cut
	               [cuts[0],   8.3], $     ;MSP Inflow
	               [cuts[0],   0.3], $     ;MSP Separatrix
	               [cuts[0],  -0.8], $     ;X-Point
	               [cuts[0],  -1.8], $     ;MSH Separatrix
	               [cuts[0],  -8.3], $     ;MSH Inflow
	               ;2nd Cut
	               [cuts[1],   8.3], $     ;MSP Inflow
	               [cuts[1],   3.3], $     ;MSP Separatrix
	               [cuts[1],  -1.6], $     ;Central Current Sheet
	               [cuts[1],  -7.8], $     ;MSH Separatrix
	               [cuts[1], -10.0], $     ;MSH Inflow
	               ;3rd Cut
	               [cuts[2],  10.0], $     ;MSP Inflow
	               [cuts[2],   6.7], $     ;MSP Separatrix
	               [cuts[2],  -2.8], $     ;Central Current Sheet
	               [cuts[2], -13.1], $     ;MSH Separatrix
	               [cuts[2], -16.7]]       ;MSH Inflow
	
	;Magnetopause?
	if coord_sys eq 'MAGNETOPAUSE' then begin
		dist_loc      = reverse(dist_loc, 1)
		dist_loc[0,*] = -dist_loc[0,*]
	endif
	
	;Ion Scale?
	if ion_scale then begin
		osim -> GetInfo, MI_ME=mi_me
		dist_width    /= sqrt(mi_me)
		dist_loc[0,*] /= sqrt(mi_me)
	endif

	;Retrieve the color image
	img_temp = win1['Color ' + im_name]
	
	;Draw boxes where the distributions are taken
	for i = 0, n_elements(dist_loc[0,*]) - 1 do begin
		;Coordinates of the box corners.
		dx_minus = dist_loc[0,i] - dist_width[0]
		dx_plus  = dist_loc[0,i] + dist_width[0]
		dy_minus = dist_loc[1,i] - dist_width[1]
		dy_plus  = dist_loc[1,i] + dist_width[1]
		xcoords = [dx_minus, dx_plus,  dx_plus, dx_minus, dx_minus]
		ycoords = [dy_minus, dy_minus, dy_plus, dy_plus,  dy_minus]
		
		;Draw the box
		box = MrPlotS(xcoords, ycoords, COLOR='white', NAME='DistBox' + strtrim(i,2), TARGET=img_temp, /DATA)
	endfor

;-----------------------------------------------------
; Asymm-Scan/By1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By1'
	cuts      = [43.5, 43.8, 44.5]
	tIndex    = 30
	xrange    = [1,-1]
	zrange    = [42.0, 45.0]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	im_name   = 'Jey'
	c_name    = 'Ay'
	win2 = MrSim_XProximity(theSim, cuts, tIndex, $
	                        C_NAME    = c_name, $
	                        COORD_SYS = coord_sys, $
	                        IM_NAME   = im_name, $
	                        ION_SCALE = ion_scale, $
	                        MVA_FRAME = mva_frame, $
	                        XRANGE    = xrange, $
	                        ZRANGE    = zrange)

;-----------------------------------------------------
; Asymm-Scan/By1: Mark Distribution Functions \\\\\\\\
;-----------------------------------------------------
	;Distribution information (de)
	dist_width = [0.5, 0.5]
	dist_loc = [ $
	            ;1st Cut
	            [cuts[0],   4.2], $     ;MSP Inflow
	            [cuts[0],   1.8], $     ;MSP Separatrix
	            [cuts[0],   0.3], $     ;X-Point
	            [cuts[0],  -0.4], $     ;MSH Separatrix
	            [cuts[0],  -1.7], $     ;MSH Inflow
	            ;2nd Cut
	            [cuts[1],   5.0], $     ;MSP Inflow
	            [cuts[1],   1.4], $     ;MSP Separatrix
	            [cuts[1],   0.0], $     ;Central Current Sheet
	            [cuts[1],  -1.4], $     ;MSH Separatrix
	            [cuts[1],  -4.2], $     ;MSH Inflow
	            ;3rd Cut
	            [cuts[2],   5.8], $     ;MSP Inflow
	            [cuts[2],   2.9], $     ;MSP Separatrix
	            [cuts[2],  -3.3], $     ;Central Current Sheet
	            [cuts[2],  -4.5], $     ;MSH Separatrix
	            [cuts[2],  -5.8]]       ;MSH Inflow
	
	;Magnetopause?
	if coord_sys eq 'MAGNETOPAUSE' then begin
		dist_loc      = reverse(dist_loc, 1)
		dist_loc[0,*] = -dist_loc[0,*]
	endif
	
	;Ion Scale?
	if ion_scale then begin
		osim -> GetInfo, MI_ME=mi_me
		dist_width    /= sqrt(mi_me)
		dist_loc[0,*] /= sqrt(mi_me)
	endif

	;Retrieve the color image
	img_temp = win2['Color ' + im_name]
	
	;Draw boxes where the distributions are taken
	for i = 0, n_elements(dist_loc[0,*]) - 1 do begin
		;Coordinates of the box corners.
		dx_minus = dist_loc[0,i] - dist_width[0]
		dx_plus  = dist_loc[0,i] + dist_width[0]
		dy_minus = dist_loc[1,i] - dist_width[1]
		dy_plus  = dist_loc[1,i] + dist_width[1]
		xcoords = [dx_minus, dx_plus,  dx_plus, dx_minus, dx_minus]
		ycoords = [dy_minus, dy_minus, dy_plus, dy_plus,  dy_minus]
		
		;Draw the box
		box = MrPlotS(xcoords, ycoords, COLOR='white', NAME='DistBox' + strtrim(i,2), TARGET=img_temp, /DATA)
	endfor

;-----------------------------------------------------
; Asymm-3D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	simName   = 'Asymm-3D'
	cuts      = [45.2, 43.7, 35.2]
	tIndex    = 120816
	xrange    = [5, -4]
	ycell     = 860
	zrange    = 45.3 + [-15, 15]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	im_name   = 'Jey'
	
	;Create the simulation object
	theSim = MrSim_Create(simName, tIndex, $
	                      COORD_SYS = coord_sys, $
	                      ION_SCALE = ion_scale, $
	                      MVA_FRAME = mva_frame, $
	                      XRANGE    = xrange, $
	                      ZRANGE    = zrange)
	ycoord = theSim -> Cell2Coord(ycell, /Y)
	theSim.YRANGE = [ycoord, ycoord]
	
	;Create the graphic
	win3 = MrSim_XProximity(theSim, cuts, IM_NAME = im_name)

	;Destroy the simulation object
	obj_destroy, theSim

;-----------------------------------------------------
; Finish \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
    ;Group the figures
    wins = [win1, win2, win3]
    
    ;Name the figures
    fnames = ['Prox_Asymm-Scan-By0', 'Prox_Asymm-Scan-By1', 'Prox_Asymm-3D']
    
    return, wins
end


;+
;   Create Figure 2: eMap for Asymm-Scan/By0 and Asymm-Scan/By1.
;-
function MrProxFigs_Figure2, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win1)   then obj_destroy, win1
		if obj_valid(win2)   then obj_destroy, win2
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim      = 'Asymm-Scan/By0'
	time        = 28
	coord_sys   = 'Magnetopause'
	ion_scale   = 0
	mva_frame   = 1
	
	;Create the simulation object
	oSim = MrSim_Create(theSim, time, $
	                    COORD_SYS = coord_sys, $
	                    ION_SCALE = ion_scale, $
	                    MVA_FRAME = mva_frame)

	;Describe the distributions
	dist_type   = 'Vpar-Vperp'
	dist_width  = [0.5, 0.5]
	dist_layout = [5, 3]
	v_va        = 1

	;Locations in SIMULATION coordinates.
	cuts        = [367.7, 349.0, 330.0]
	dist_loc    = [ $
	               ;1st Cut
	               [cuts[0],   8.3], $     ;MSP Inflow
	               [cuts[0],   0.3], $     ;MSP Separatrix
	               [cuts[0],  -0.8], $     ;X-Point
	               [cuts[0],  -1.8], $     ;MSH Separatrix
	               [cuts[0],  -8.3], $     ;MSH Inflow
	               ;2nd Cut
	               [cuts[1],   8.3], $     ;MSP Inflow
	               [cuts[1],   3.3], $     ;MSP Separatrix
	               [cuts[1],  -1.6], $     ;Central Current Sheet
	               [cuts[1],  -7.8], $     ;MSH Separatrix
	               [cuts[1], -10.0], $     ;MSH Inflow
	               ;3rd Cut
	               [cuts[2],  10.0], $     ;MSP Inflow
	               [cuts[2],   6.7], $     ;MSP Separatrix
	               [cuts[2],  -2.8], $     ;Central Current Sheet
	               [cuts[2], -13.1], $     ;MSH Separatrix
	               [cuts[2], -16.7]]       ;MSH Inflow

	;Change to MAGNETOPAUSE coordinates
	if coord_sys eq 'MAGNETOPAUSE' then begin
		dist_loc      = reverse(dist_loc, 1)
		dist_loc[0,*] = -dist_loc[0,*]
	endif
	
	;Ion Scale?
	if ion_scale then begin
		oSim -> GetInfo, MI_ME=mi_me
		dist_width    /= sqrt(mi_me)
		dist_loc[0,*] /= sqrt(mi_me)
	endif

	;Create the distribution window.
	win1 = MrSim_eMap(oSim, dist_type, dist_loc, dist_width, $
	                  CIRCLES    = circles, $
	                  HGAP       = hGap, $
	                  NBINS      = nBins, $
	                  LAYOUT     = dist_layout, $
	                  LOCATION   = location, $
	                  POSITIONS  = positions, $
	                  SIM_OBJECT = oSim, $
	                  VGAP       = vGap, $
	                  V_VA       = v_va, $
	                  XSIZE      = xsize, $
	                  YSIZE      = ysize)
	
	;Destroy the simulation object
	obj_destroy, oSim
	
	;Turn refresh off while we update things.
	win1 -> Refresh, /DISABLE
	win1.ysize = 350
	
	;Get the title, then trash all text items
	title = win1['eMap Title'].string
	win1 -> Remove, TYPE='MrText'

	;Set the [xy]tick interval
	win1 -> SetGlobal, XTICKINTERVAL=10, YTICKINTERVAL=5
	
	;Create a new title
	title = MrText(0.5, 0.94, title, ALIGNMENT=0.5, VERTICAL_ALIGNMENT=0.5)
	
;-----------------------------------------------------
; Asymm-Scan/By1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Describe the simulation
	theSim    = 'Asymm-Scan/By1'
	time      = 30
	xrange    = [10, -10]
	zrange    = [420, 450]
	ion_scale = 0
	mva_frame = 1
	coord_sys = 'Magnetopause'
	
	;Create the simulation object
	oSim = MrSim_Create(theSim, time, $
	                    COORD_SYS = coord_sys, $
	                    ION_SCALE = ion_scale, $
	                    MVA_FRAME = mva_frame)
	
	;Distribution information
	v_va        = 1
	dist_type   = 'Vpar-Vperp'
	dist_width  = [0.5, 0.5]
	dist_layout = [5, 3]
	
	;Locations in SIMULATION coordinates.
	cuts     = [434.5, 438.0, 445.0]
	dist_loc = [ $
	            ;1st Cut
	            [cuts[0],   4.2], $     ;MSP Inflow
	            [cuts[0],   1.8], $     ;MSP Separatrix
	            [cuts[0],   0.3], $     ;X-Point
	            [cuts[0],  -0.4], $     ;MSH Separatrix
	            [cuts[0],  -1.7], $     ;MSH Inflow
	            ;2nd Cut
	            [cuts[1],   5.0], $     ;MSP Inflow
	            [cuts[1],   1.4], $     ;MSP Separatrix
	            [cuts[1],   0.0], $     ;Central Current Sheet
	            [cuts[1],  -1.4], $     ;MSH Separatrix
	            [cuts[1],  -4.2], $     ;MSH Inflow
	            ;3rd Cut
	            [cuts[2],   5.8], $     ;MSP Inflow
	            [cuts[2],   2.9], $     ;MSP Separatrix
	            [cuts[2],  -3.3], $     ;Central Current Sheet
	            [cuts[2],  -4.5], $     ;MSH Separatrix
	            [cuts[2],  -5.8]]       ;MSH Inflow

	;Change to MAGNETOPAUSE coordinates
	if coord_sys eq 'MAGNETOPAUSE' then begin
		dist_loc      = reverse(dist_loc, 1)
		dist_loc[0,*] = -dist_loc[0,*]
	endif
	
	;Ion Scale?
	if ion_scale then begin
		oSim -> GetInfo, MI_ME=mi_me
		dist_width    /= sqrt(mi_me)
		dist_loc[0,*] /= sqrt(mi_me)
	endif

	;Create the distribution window.
	win2 = MrSim_eMap(theSim, dist_type, dist_loc, dist_width, time, $
	                  C_NAME     = c_name, $
	                  CIRCLES    = circles, $
	                  HGAP       = hGap, $
	                  IM_NAME    = im_name, $
	                  NBINS      = nBins, $
	                  LAYOUT     = dist_layout, $
	                  LOCATION   = location, $
	                  POSITIONS  = positions, $
	                  SIM_OBJECT = oSim, $
	                  VGAP       = vGap, $
	                  V_VA       = v_va, $
	                  XSIZE      = xsize, $
	                  YSIZE      = ysize, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange, $
	                  ION_SCALE  = ion_scale, $
	                  COORD_SYS  = coord_sys)
	
	;Destroy the simulation object
	obj_destroy, oSim
	
	;Turn refresh off while we update things.
	win2 -> Refresh, /DISABLE
	win2.ysize = 350
	
	;Get the title, then trash all text items
	title = win2['eMap Title'].string
	win2 -> Remove, TYPE='MrText'

	;Set the [xy]tick interval
	win2 -> SetGlobal, XTICKINTERVAL=10, YTICKINTERVAL=5
	
	;Create a new title
	title = MrText(0.5, 0.94, title, ALIGNMENT=0.5, VERTICAL_ALIGNMENT=0.5)

;-----------------------------------------------------
; Finish \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Output file names
	fnames = 'eMap_' + dist_type + ['_Asymm-Scan-By0', '_Asymm-Scan-By1']
	
	win1 -> Refresh
	win2 -> Refresh
	return, [win1, win2]
end


;+
;   Create Figure 2: eMap for Asymm-Scan/By0 and Asymm-Scan/By1.
;-
function MrProxFigs_XPoint, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win)  then obj_destroy, win
		if obj_valid(win2) then obj_destroy, win2
		if obj_valid(oSim) then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim    = 'Asymm-Scan/By0'
	time      = 28
	coord_sys = 'Magnetopause'
	xrange    = [2,-2]
	zrange    = 36.77 + [-5, 5]
	ion_scale = 1
	mva_frame = 1
	
	;Create the simulation object
	oSim = MrSim_Create(theSim, time, $
	                    COORD_SYS = coord_sys, $
	                    ION_SCALE = ion_scale, $
	                    MVA_FRAME = mva_frame, $
	                    XRANGE    = xrange, $
	                    ZRANGE    = zrange)
	
	;Create a window
	win = MrWindow(XSIZE=700, OXMARGIN=[33,15], REFRESH=0)

	;Find the X-Point
	xpt = MrSim_XPoint(oSim)

	;Plot the reconnection rate.
	rate = MrSim_ReconRate(oSim)
	oSim -> SetProperty, TIME=time

	;Convert time to t*wci^-2.
	ti   = lindgen(oSim.nTimes)
	twci = oSim -> tIndex2txwci(ti)
	
	;Plot the reconnection rate
	gPlot = MrPlot(twci, rate, /CURRENT, $
	               TITLE  = 'Reconnection Rate', $
	               XTITLE = 't$\Omega$$\downci$$\up-1$', $
	               YTITLE = 'R')
	
	;Determine the current time
	t = oSim -> tIndex2txwci(oSim.time)

	;Draw a line through the current time
	!Null = MrPlotS([t, t], [min(rate, MAX=rMax), rMax], $
	                COLOR     = 'Forest Green', $
	                LINESTYLE = 'Dash', $
	                TARGET    = gPlot)
	 
	
	
	;Plot cuts of BL, n, and EN through the X-line
	lc1 = MrSim_LineCut(oSim, 'Bz', xpt[1], /CURRENT, /HORIZONTAL)
	pos = lc1.position
	lc2 = MrSim_LineCut(oSim, 'ni', xpt[1], /CURRENT, /HORIZONTAL)
	lc3 = MrSim_LineCut(oSim, 'Ex', xpt[1], /CURRENT, /HORIZONTAL)
	
	;Move ni and Ex plots onto Bz to create independent data spaces for them
	lc1 -> SetProperty, XRANGE=[1, -1]
	lc2 -> SetProperty, XRANGE=[1, -1], COLOR='Blue', POSITION=pos, XSTYLE=5, YSTYLE=5, TITLE=''
	lc3 -> SetProperty, XRANGE=[1, -1], COLOR='Red',  POSITION=pos, XSTYLE=5, YSTYLE=5, TITLE=''
	win -> TrimLayout
	
	;Create axes for ni and Ex
	ax1 = MrAxis('Y', TARGET=lC2, LOCATION='Right', COLOR='Blue', TICKDIR=1, TITLE='n$\downi$')
	ax2 = MrAxis('Y', TARGET=lC3, LOCATION='Right', COLOR='Red',  TICKDIR=1, TITLE='E$\downN$', OFFSET=7)
	
	;Find where Bz = 0
	lc1   -> GetData, x, Bz
	!Null  = min(Bz, iZero, /ABSOLUTE)
	
	;Draw a vertical line throubh Bz = 0
	!Null = MrPlotS(x[[iZero, iZero]], lc1.yrange, $
	                COLOR     = 'Forest Green', $
	                LINESTYLE = 'Dash', $
	                TARGET    = lc1)
	
	;Color image of Jey in a new window
	win2 = MrSim_ColorSlab(oSim, 'Jey', C_NAME='Ay')
	win2 -> Refresh, /DISABLE

	;Put an X at the x-point
	img  = win2['Color Jey']
	gXpt = MrPlotS(xpt[0], xpt[1], $
	               LINESTYLE = 'None', $
	               PSYM      = 7, $
	               SYMCOLOR  = 'Cyan', $
	               SYMSIZE   = 2, $
	               TARGET    = img, $
	               THICK     = 2)
	
	;Update asthetics
	title = img.TITLE
	img                   -> SetProperty, POSITION=[0.1, 0.2, 0.3, 0.77], TITLE='', XRANGE=[4,-2]
	win2['Ay Contours']   -> SetProperty, C_THICK=1.0, C_COLOR='Grey'
	win2['CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', TITLE=title, WIDTH=0.5, XTICKS=2, $
	                                      XTICKFORMAT='(f0.2)'
	
	;Put the contents of the ColorSlab window into the LineCut window.
	allGfx = win2 -> Get(/ALL)
	foreach gfx, allGfx do gfx -> SwitchWindows, win
	obj_destroy, win2
	
	
	;Destroy the simulation object
	obj_destroy, oSim
	
	win -> Refresh
	return, win
end


;+
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function MrProximity_Figures, figure, $
PNG=png, $
EPS=eps, $
PS=ps, $
IM_PNG=im_png, $
SAVE=tf_save
	compile_opt strictarr

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if max(obj_valid(win)) then obj_destroy, win
		void = cgErrorMSG()
		return, obj_new()
	endif

;---------------------------------------------------------------------
; Info ///////////////////////////////////////////////////////////////
;---------------------------------------------------------------------

	;Current list of figures
	list_of_figures = [['Figure 1', 'X-Line proximity for Asymm-Scan/By0, Asymm-Scan/By1, and Asym3D'], $
	                   ['Figure 2', 'eDists at inflow, separatrix, center for Asymm-Scan/By0 and Asymm-Scan/By1'], $
	                   ['X-Point',  'Cuts of BL, n, and EN through the X-point, Jey with X-point marked, ReconRate.']]

	;Print the list of figures?
	if n_elements(figure) eq 0 then begin
		len = strtrim(max(strlen(list_of_figures[0,*])) > 7, 2)
		print, 'FIGURES', 'DESCRIPTION', FORMAT='(a' + len + ', 4x, a11)'
		print, list_of_figures, FORMAT='(a-' + len + ', 4x, a0)'
		return, obj_new()
	endif

;---------------------------------------------------------------------
; Create Figures /////////////////////////////////////////////////////
;---------------------------------------------------------------------

	_figure = strupcase(figure)
	tf_save = keyword_set(tf_save)

	case _figure of
		'FIGURE 1': win = MrProxFigs_Figure1(FNAMES=fnames)
		'FIGURE 2': win = MrProxFigs_Figure2(FNAMES=fnames)
		'X-POINT':  win = MrProxFigs_XPoint(FNAMES=fnames)
		else: message, 'Figure "' + figure + '" not an option.', /INFORMATIONAL
	endcase
    
;---------------------------------------------------------------------
; Save to File? //////////////////////////////////////////////////////
;---------------------------------------------------------------------
	if keyword_set(tf_save) then begin
		;Save as what?
		eps    = keyword_set(eps)
		ps     = keyword_set(ps)
		png    = keyword_set(png)
		im_png = keyword_set(im_png)
		if eps + png + ps + im_png eq 0 then begin
			eps    = 1
			ps     = 1
			png    = 1
			im_png = 1
		endif

		;Number of files and where to save them.
		nWins = n_elements(win)
		froot = '/home/argall/figures/'

		;Single window
		if nWins eq 1 then begin
			;Create the file name
			fname = 'MrProximity_' + idl_validname(figure, /CONVERT_ALL)
			fbase = filepath(fname, ROOT_DIR=froot)
	
			;Save a variety of file types.
			win -> Refresh
			if im_png then win -> Save, fbase + '_im.png'
			if eps    then win -> Save, fbase + '.eps'
			if ps     then win -> Save, fbase + '.ps'
	
			;Take a snapshot
			if png then begin
				win.SAVEAS -> SetProperty, IM_RASTER=0
				win -> Save, fbase + '-ss.png'
				win.SAVEAS -> SetProperty, IM_RASTER=1
			endif
		
		;Multiple windows
		endif else begin
			for i = 0, nWins - 1 do begin
				fname = FilePath('MrProximity_' + fnames[i], ROOT_DIR=froot)
				if im_png then win[i] -> Save, fname + '_im.png'
				if eps    then win[i] -> Save, fname + '.eps'
				if ps     then win[i] -> Save, fname + '.ps'
				
				;Snapshot
				if png then begin
					win[i].SAVEAS -> SetProperty, IM_RASTER=0
					win[i] -> Save, fname + '-ss.png'
					win[i].SAVEAS -> SetProperty, IM_RASTER=1
				endif
			endfor
		endelse
	endif

	return, win
end