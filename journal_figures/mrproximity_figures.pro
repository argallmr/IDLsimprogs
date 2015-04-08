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
function MrProxFigs_OhmsLaw, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win) then obj_destroy, win
		if n_elements(by0_wins)  gt 0 then obj_destroy, by0_wins
		if n_elements(by01_wins) gt 0 then obj_destroy, by01_wins
		if n_elements(by1_wins)  gt 0 then obj_destroy, by1_wins
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0 t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Simulation
	theSim       = 'Asymm-Scan/By0'
	tIndex       = 28
	xrange       = [1,-2]
	zrange       = 36.77 + [-5, 5]
	coord_system = 'Magnetopause'
	mva_frame    = 1
	ion_scale    = 1
	horizontal   = 1
	component    = 'X'
	oSim   = MrSim_Create(theSim, tIndex, $
	                      COORD_SYSTEM = coord_system, $
	                      MVA_FRAME    = mva_frame, $
	                      ION_SCALE    = ion_scale, $
	                      XRANGE       = xrange, $
	                      ZRANGE       = zrange)

	;Do for various times
	times      = [28, 32, 44]
	nTimes     = n_elements(times)
	by0_wins   = objarr(nTImes)
	by0_fnames = strarr(nTimes)
	
	;Step through each time
	for i = 0, nTimes - 1 do begin
		;Set the time
		oSim.TIME = times[i]

		;Find the X-Point
		xpt = MrSim_XPoint(oSim)

		;Ohm's Law
		win = MrSim_OhmsLaw(oSim, component, xpt[1], HORIZONTAL=horizontal)
		win -> Refresh, /DISABLE
	
		;Extract graphics to use as targets
		gET = win['Total EX']
		gEC = win['EX vs. EC']
		gEH = win['EX vs. Hall E']
		gEP = win['EX vs. E Pressure']
		gEI = win['EX vs. E eInert']

		;Draw vertical lines to mark the x-point
		oSim -> GetProperty, XSIM=xSim
		line = MrPlotS( [xpt[0], xpt[0]], gET.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 E', $
		                TARGET    = gET)
		line = MrPlotS( [xpt[0], xpt[0]], gEC.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EC', $
		                TARGET    = gEC)
		line = MrPlotS( [xpt[0], xpt[0]], gEH.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EH', $
		                TARGET    = gEH)
		line = MrPlotS( [xpt[0], xpt[0]], gEP.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EP', $
		                TARGET    = gEP)
		line = MrPlotS( [xpt[0], xpt[0]], gEI.YRANGE, $
		                COLOR     = 'Grey', $
		                NAME      = 'BL=0 EI', $
		                LINESTYLE = '--', $
		                TARGET    = gEI)
		
		;Change the size of the legend
		ohLegend = win["Ohm's Law"]
		ohLegend.VERTICAL_SPACING = 1.0
		for j = 0, 5 do ohLegend -> SetProperty, j, TEXT_SIZE=1.0
		
		;Keep the window
		by0_winS[i]    = win
		by0_fnames[i]  = 'asymm-scan-by0_OhmsLaw_t' + strtrim(times[i], 2) + '_xline'
		win           -> Refresh
	endfor

	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By1 t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Simulation
	theSim       = 'Asymm-Scan/By1'
	tIndex       = 30
	xrange       = [1,-2]
	zrange       = [42.0, 45.0]
	coord_system = 'Magnetopause'
	mva_frame    = 1
	ion_scale    = 1
	horizontal   = 1
	component    = 'X'
	oSim   = MrSim_Create(theSim, tIndex, $
	                      COORD_SYSTEM = coord_system, $
	                      MVA_FRAME    = mva_frame, $
	                      ION_SCALE    = ion_scale, $
	                      XRANGE       = xrange, $
	                      ZRANGE       = zrange)

	;Do for various times
	times      = [30]
	nTimes     = n_elements(times)
	by1_wins   = objarr(nTImes)
	by1_fnames = strarr(nTimes)
	
	;Step through each time
	for i = 0, nTimes - 1 do begin
		;Set the time
		oSim.TIME = times[i]

		;Find the X-Point
		xpt = MrSim_XPoint(oSim)

		;Ohm's Law
		win = MrSim_OhmsLaw(oSim, component, xpt[1], HORIZONTAL=horizontal)
		win -> Refresh, /DISABLE
	
		;Extract graphics to use as targets
		gET = win['Total EX']
		gEC = win['EX vs. EC']
		gEH = win['EX vs. Hall E']
		gEP = win['EX vs. E Pressure']
		gEI = win['EX vs. E eInert']

		;Draw vertical lines to mark the x-point
		oSim -> GetProperty, XSIM=xSim
		line = MrPlotS( [xpt[0], xpt[0]], gET.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 E', $
		                TARGET    = gET)
		line = MrPlotS( [xpt[0], xpt[0]], gEC.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EC', $
		                TARGET    = gEC)
		line = MrPlotS( [xpt[0], xpt[0]], gEH.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EH', $
		                TARGET    = gEH)
		line = MrPlotS( [xpt[0], xpt[0]], gEP.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EP', $
		                TARGET    = gEP)
		line = MrPlotS( [xpt[0], xpt[0]], gEI.YRANGE, $
		                COLOR     = 'Grey', $
		                NAME      = 'BL=0 EI', $
		                LINESTYLE = '--', $
		                TARGET    = gEI)
		
		;Change the size of the legend
		ohLegend = win["Ohm's Law"]
		ohLegend.VERTICAL_SPACING = 1.0
		for j = 0, 5 do ohLegend -> SetProperty, j, TEXT_SIZE=1.0
		
		;Keep the window
		by1_wins[i]    = win
		by1_fnames[i]  = 'asymm-scan-by1_OhmsLaw_t' + strtrim(times[i], 2) + '_xline'
		win           -> Refresh
	endfor
	
	;Destroy the simulation object
	obj_destroy, oSim
	
;-----------------------------------------------------
; Asymm-Scan/By0.1 t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Simulation
	theSim       = 'Asymm-Scan/By0.1'
	tIndex       = 28
	xrange       = [1,-2]
	zrange       = 40.5 + [-2.0, 2.0]
	coord_system = 'Magnetopause'
	mva_frame    = 1
	ion_scale    = 1
	horizontal   = 1
	component    = 'X'
	oSim   = MrSim_Create(theSim, tIndex, $
	                      COORD_SYSTEM = coord_system, $
	                      MVA_FRAME    = mva_frame, $
	                      ION_SCALE    = ion_scale, $
	                      XRANGE       = xrange, $
	                      ZRANGE       = zrange)

	;Do for various times
	times       = [28, 32, 37]
	nTimes      = n_elements(times)
	by01_wins   = objarr(nTImes)
	by01_fnames = strarr(nTimes)
	
	;Step through each time
	for i = 0, nTimes - 1 do begin
		;Set the time
		oSim.TIME = times[i]

		;Find the X-Point
		xpt = MrSim_XPoint(oSim)

		;Ohm's Law
		win = MrSim_OhmsLaw(oSim, component, xpt[1], HORIZONTAL=horizontal)
		win -> Refresh, /DISABLE
	
		;Extract graphics to use as targets
		gET = win['Total EX']
		gEC = win['EX vs. EC']
		gEH = win['EX vs. Hall E']
		gEP = win['EX vs. E Pressure']
		gEI = win['EX vs. E eInert']

		;Draw vertical lines to mark the x-point
		oSim -> GetProperty, XSIM=xSim
		line = MrPlotS( [xpt[0], xpt[0]], gET.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 E', $
		                TARGET    = gET)
		line = MrPlotS( [xpt[0], xpt[0]], gEC.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EC', $
		                TARGET    = gEC)
		line = MrPlotS( [xpt[0], xpt[0]], gEH.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EH', $
		                TARGET    = gEH)
		line = MrPlotS( [xpt[0], xpt[0]], gEP.YRANGE, $
		                COLOR     = 'Grey', $
		                LINESTYLE = '--', $
		                NAME      = 'BL=0 EP', $
		                TARGET    = gEP)
		line = MrPlotS( [xpt[0], xpt[0]], gEI.YRANGE, $
		                COLOR     = 'Grey', $
		                NAME      = 'BL=0 EI', $
		                LINESTYLE = '--', $
		                TARGET    = gEI)
		
		;Change the size of the legend
		ohLegend = win["Ohm's Law"]
		ohLegend.VERTICAL_SPACING = 1.0
		for j = 0, 5 do ohLegend -> SetProperty, j, TEXT_SIZE=1.0
		                
		;Keep the window
		by01_wins[i]   = win
		by01_fnames[i] = 'asymm-scan-by01_OhmsLaw_t' + strtrim(times[i], 2) + '_xline'
		win           -> Refresh
	endfor
	
	;Destroy the simulation object
	obj_destroy, oSim
	
;-----------------------------------------------------
; Asymm-Scan/By0.1 t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Gather windows
	wins = [by0_wins, by01_wins, by1_wins]
	
	;Gather file names
	fnames = [by0_fnames, by01_fnames, by1_fnames]

	return, wins
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
		if n_elements(win1) gt 0 then obj_destroy, win1
		if n_elements(win2) gt 0 then obj_destroy, win2
		if n_elements(win3) gt 0 then obj_destroy, win3
		if obj_valid(win)        then obj_destroy, win
		if obj_valid(oSim)       then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim    = 'Asymm-Scan/By0'
	time      = 28
	coord_sys = 'Magnetopause'
	xrange    = [2,-2]
	zrange    = 36.77 + [-5, 5]
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;How many windows to create
	times  = [28, 32, 44]
	nTimes = n_elements(times)
	win1   = objarr(nTimes)
	fname1 = strarr(nTimes)
	
	;Make the plots
	for i = 0, nTimes - 1 do begin
		;Change the simulation time
		oSim.TIME = times[i]

		;Create the figure
		win = MrSim_ProxX(oSim, times[i])
		win -> Refresh, /DISABLE

		;Zoom in
;		win['Color Jey']  -> SetProperty, XRANGE=[3, -2]
		win['Cut BL']     -> SetProperty, XRANGE=[1, -1]
		win['Cut ni']     -> SetProperty, XRANGE=[1, -1]
		win['Cut EN']     -> SetProperty, XRANGE=[1, -1]
		win['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
		win['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
		win               -> Remove, win['X-Point']

		;Store the window and file anme
		win -> Refresh
		win1[i]   = win
		fname1[i] = 'Asymm-Scan-By0_ProxX_t' + strtrim(oSim.time, 2)
	endfor
	
	;Destroy the object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By1, t=30 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By1'
	time      = 30
	xrange    = [2,-2]
	zrange    = [42.0, 45.0]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;How many windows to create
	times  = [30]
	nTimes = n_elements(times)
	win2   = objarr(nTimes)
	fname2 = strarr(nTimes)
	
	;Make the plots
	for i = 0, nTimes - 1 do begin
		;Change the simulation time
		oSim.TIME = times[i]
	
		;Create the figure
		win = MrSim_ProxX(oSim, times[i])
		win -> Refresh, /DISABLE
		
		;Zoom in
;		win['Color Jey']  -> SetProperty, XRANGE=[3, -2]
		win['Cut BL']     -> SetProperty, XRANGE=[1, -1]
		win['Cut ni']     -> SetProperty, XRANGE=[1, -1]
		win['Cut EN']     -> SetProperty, XRANGE=[1, -1]
		win['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
		win['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
		win               -> Remove, win['X-Point']

		;Refresh
		win -> Refresh
		fname2[i] = 'Asymm-Scan-By1_ProxX_t' + strtrim(oSim.time, 2)
		win2[i]   = win
	endfor
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By0.1, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0.1'
	time      = 28
	xrange    = [2,-2]
	zrange    = 40.5 + [-2.0, 2.0]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;How many windows to create
	times  = [28, 32, 37]
	nTimes = n_elements(times)
	win3   = objarr(nTimes)
	fname3 = strarr(nTimes)
	
	;Make the plots
	for i = 0, nTimes - 1 do begin
		;Change the simulation time
		oSim.TIME = times[i]

		;Create the figure
		win = MrSim_ProxX(oSim, times[i])
		win -> Refresh, /DISABLE
	
		;Zoom in
;		win['Color Jey']  -> SetProperty, XRANGE=[3, -2]
		win['Cut BL']     -> SetProperty, XRANGE=[1, -1]
		win['Cut ni']     -> SetProperty, XRANGE=[1, -1]
		win['Cut EN']     -> SetProperty, XRANGE=[1, -1]
		win['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
		win['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
		win               -> Remove, win['X-Point']

		;Refresh
		win -> Refresh
		win3[i]   = win
		fname3[i] = 'Asymm-Scan-By01_ProxX_t' + strtrim(oSim.time, 2)
	endfor
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Return \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Return the window
	fnames = [fname1, fname2, fname3]
	return, [win1, win2, win3]
end


;+
;   Create Figure 2: eMap for Asymm-Scan/By0 and Asymm-Scan/By1.
;-
function MrProxFigs_TProx, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win1)  then obj_destroy, win1
		if obj_valid(win2)  then obj_destroy, win2
		if obj_valid(win3)  then obj_destroy, win3
		if obj_valid(win4)  then obj_destroy, win4
		if obj_valid(win5)  then obj_destroy, win5
		if obj_valid(win6)  then obj_destroy, win6
		if obj_valid(win7)  then obj_destroy, win7
		if obj_valid(oSim)  then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim    = 'Asymm-Scan/By0'
	time      = 28
	coord_sys = 'Magnetopause'
	xrange    = [2,-2]
	zrange    = 36.77 + [-5, 5]
	ion_scale = 1
	mva_frame = 1
	
	;Create the figure
	win1 = MrSim_ProxT(theSim, time, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  SIM_OBJECT = oSim, $s
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win1 -> Refresh, /DISABLE
	
	;Zoom in
;	win1['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win1['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win1['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win1['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win1['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win1['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win1['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win1['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win1 -> Refresh
	fname1 = 'Asymm-Scan-By0_ProxT_t' + strtrim(time, 2)

;-----------------------------------------------------
; Asymm-Scan/By0, t=32 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Change times
	oSim.TIME = 32
	
	;Create the figure
	win2 = MrSim_ProxT(oSim, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win2 -> Refresh, /DISABLE
	
	;Zoom in
;	win2['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win2['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win2['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win2['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win2['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win2['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win2['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win2['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win2 -> Refresh
	fname2 = 'Asymm-Scan-By0_ProxT_t' + strtrim(oSim.time, 2)

;-----------------------------------------------------
; Asymm-Scan/By0, t=44 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Change times
	oSim.TIME = 44
	
	;Create the figure
	win3 = MrSim_ProxT(oSim, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win3 -> Refresh, /DISABLE
	
	;Zoom in
;	win3['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win3['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win3['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win3['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win3['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win3['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win3['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win3['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win3 -> Refresh
	fname3 = 'Asymm-Scan-By0_ProxT_t' + strtrim(oSim.time, 2)
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By1, t=30 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By1'
	time      = 30
	xrange    = [2,-2]
	zrange    = [42.0, 45.0]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	
	;Create the figure
	win4 = MrSim_ProxT(theSim, time, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  SIM_OBJECT = oSim, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win4 -> Refresh, /DISABLE
	
	;Zoom in
;	win4['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win4['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win4['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win4['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win4['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win4['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win4['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win4['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win4 -> Refresh
	fname4 = 'Asymm-Scan-By1_ProxT_t' + strtrim(oSim.time, 2)
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By0.1, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0.1'
	time      = 28
	xrange    = [2,-2]
	zrange    = 40.5 + [-2.0, 2.0]
	coord_sys = 'Magnetopause'
	ion_scale = 1
	mva_frame = 1
	
	;Create the figure
	win5 = MrSim_ProxT(theSim, time, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  SIM_OBJECT = oSim, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win5 -> Refresh, /DISABLE
	
	;Zoom in
;	win5['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win5['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win5['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win5['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win5['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win5['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win5['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win5['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win5 -> Refresh
	fname5 = 'Asymm-Scan-By01_ProxT_t' + strtrim(oSim.time, 2)

;-----------------------------------------------------
; Asymm-Scan/By0.1, t=32 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Change times
	oSim.TIME = 32
	
	;Create the figure
	win6 = MrSim_ProxT(oSim, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win6 -> Refresh, /DISABLE
	
	;Zoom in
;	win6['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win6['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win6['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win6['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win6['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win6['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win6['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win6['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win6 -> Refresh
	fname6 = 'Asymm-Scan-By01_ProxT_t' + strtrim(oSim.time, 2)

;-----------------------------------------------------
; Asymm-Scan/By0.1, t=37 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Change times
	oSim.TIME = 37
	
	;Create the figure
	win7 = MrSim_ProxT(oSim, $
	                  COORD_SYS  = coord_sys, $
	                  ION_SCALE  = ion_scale, $
	                  MVA_FRAME  = mva_frame, $
	                  XRANGE     = xrange, $
	                  ZRANGE     = zrange)
	win7 -> Refresh, /DISABLE
	
	;Zoom in
;	win7['Color Jey']  -> SetProperty, XRANGE=[3, -2]
	win7['Cut Bz']     -> SetProperty, XRANGE=[1, -1]
	win7['Cut Jey']    -> SetProperty, XRANGE=[1, -1]
	win7['Cut ne']     -> SetProperty, XRANGE=[1, -1]
	win7['Cut Ex']     -> SetProperty, XRANGE=[1, -1]
	win7['Cut Dng_e']  -> SetProperty, XRANGE=[1, -1]
	win7['BL=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	win7['EN=0 Horiz'] -> SetProperty, XCOORDS=[1, -1]
	
	;Refresh
	win7 -> Refresh
	fname7 = 'Asymm-Scan-By01_ProxT_t' + strtrim(oSim.time, 2)
	obj_destroy, oSim

;-----------------------------------------------------
; Return \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Return the window
	fnames = [fname1, fname2, fname3, fname4, fname5, fname6, fname7]
	return, [win1, win2, win3, win4, win5, win6, win7]
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
	                   ['X-Point',  'Cuts of BL, n, and EN through the X-point, Jey with X-point marked, ReconRate.'], $
	                   ['TProx',    'Cuts of BL, Jey, n, EN, and Dng through the X-point, Jey with X-point marked, ReconRate.'], $
	                   ['OhmsLaw',  'Terms in Ohms Law at cuts through the X-line.']]

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
		'OHMSLAW':  win = MrProxFigs_OhmsLaw(FNAMES=fnames)
		'X-POINT':  win = MrProxFigs_XPoint(FNAMES=fnames)
		'TPROX':    win = MrProxFigs_TProx(FNAMES=fnames)
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