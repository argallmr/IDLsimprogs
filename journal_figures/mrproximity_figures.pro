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
;   Asymm-Scan/By0 and Asymm-Scan/By1, Asymm-3D.
;       - Proximity to the X-line
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
;   Asymm-Scan/By0 and Asymm-Scan/By1.
;       - eMaps at points along cuts that show proximity to the X-line
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
;   Asymm-Scan/By0 and Asymm-Scan/By1, Asymm-3D.
;       - Proximity to the X-line
;-
function MrProxFigs_XLineProxBy0, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win) then obj_destroy, win
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
	im_name   = 'Jy'
	c_name    = 'Ay'

	win = MrSim_XProximity(theSim, cuts, tIndex, $
	                       C_NAME     = c_name, $
	                       COORD_SYS  = coord_sys, $
	                       IM_NAME    = im_name, $
	                       ION_SCALE  = ion_scale, $
	                       MVA_FRAME  = mva_frame, $
	                       XRANGE     = xrange, $
	                       ZRANGE     = zrange)
	
	;Label in MVA components
	win['Cut Bz']       -> SetProperty, YTITLE='B$\downL$'
	win['Cut By']       -> SetProperty, YTITLE='B$\downM$'
	win['Cut Uiz']      -> SetProperty, YTITLE='U$\downiL$'
	win['Cut Ex']       -> SetProperty, YTITLE='E$\downN$'
	win['CB: Color Jy'] -> SetProperty, TITLE='J$\downM$'
	win['HLines Jy']    -> SetProperty, THICK=2
	
	fnames = 'proxfigs_by0-mrx'
	return, win
end


;+
;   Plots of By0, By0.1, By1
;       - Jey on left
;       - BL overplotted with EN, ni, UiL on right
;       - Each 1D plot has its own axis
;-
function MrProxFigs_XLine_BEVn, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win)      then obj_destroy, win
		if obj_valid(win_temp) then obj_destroy, win_temp
		if obj_valid(oSim)     then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;-----------------------------------------------------
; Asymm-Scan/By0, t=30 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0'
	time      = 35
	xrange    = 36.77 + [-5, 5]
	zrange    = [-1,2]
	coord_sys = 'Simulation'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;Create the figure
	win = MrSim_XLine_BEVn(oSim)
	win -> Refresh, /DISABLE
	win.ysize = 600
	win.ygap = 4
	
	;Change the title
	title = win['Color JM'].title
	win['Color JM'].title = 'B$\downG$=0 ' + title
	
	;Rename all of the graphics
	graphics = win -> Get(/ALL)
	foreach gfx, graphics do begin
		name = gfx.name
		gfx.name = 'By0 ' + name
	endforeach
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By0.1, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0.1'
	time      = 28
	xrange    = 40.5 + [-2.0, 2.0]
	zrange    = [-1,2]
	coord_sys = 'Simulation'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;Create the figure
	win_temp = MrSim_XLine_BEVn(oSim)
	win_temp -> Refresh, /DISABLE
	
	;Change the title
	title = win_temp['Color JM'].title
	win_temp['Color JM'].title = 'B$\downG$=0.1 ' + title
	
	;Move into row 2
	win_temp['Color JM'] -> SetLayout, [1,2]
	win_temp['Cut BL']    -> SetLayout, [2,2]

	;Rename and move all graphics
	graphics = win_temp -> Get(/ALL)
	foreach gfx, graphics do begin
		;Rename
		name = gfx.name
		gfx.name = 'By0.1 ' + name
		
		;Move
		gfx -> SwitchWindows, win
	endforeach

	;Destroy the old window
	obj_destroy, win_temp
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By1, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By1'
	time      = 28
	xrange    = 42.9 + [-2, 2]
	zrange    = [-1,2]
	coord_sys = 'Simulation'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;Create the figure
	win_temp = MrSim_XLine_BEVn(oSim)
	win_temp -> Refresh, /DISABLE
	
	;Change the title
	title = win_temp['Color JM'].title
	win_temp['Color JM'].title = 'B$\downG$=1 ' + title
	
	;Move into row 3
	win_temp['Color JM'] -> SetLayout, [1,3]
	win_temp['Cut BL']    -> SetLayout, [2,3]
	
	;Rename and move all graphics
	graphics = win_temp -> Get(/ALL)
	foreach gfx, graphics do begin
		;Rename
		name = gfx.name
		gfx.name = 'By1 ' + name
		
		;Move
		gfx -> SwitchWindows, win
	endforeach
	
	;Destroy the old window
	obj_destroy, win_temp
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asym-3D, t=120816 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	win_temp = MrProxFigs_3D_BEVn()
	win_temp -> Refresh, /DISABLE
	
	;Move into row 4
	win_temp['Color JM'] -> SetLayout, [1,4]
	win_temp['Cut BL']    -> SetLayout, [2,4]
	
	;Rename and move all graphics
	graphics = win_temp -> Get(/ALL)
	foreach gfx, graphics do begin
		;Rename
		name = gfx.name
		gfx.name = '3D-By1 ' + name
		
		;Move
		gfx -> SwitchWindows, win
	endforeach
	
	;Destroy the old window
	obj_destroy, win_temp

;-----------------------------------------------------
; General Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	win -> SetGlobal, CHARSIZE=2.0, CHARTHICK=2, XTICKLEN=0.04, YTICKLEN=0.04, $
	                  THICK=2, XTHICK=2, YTHICK=2
	win -> SetProperty, XSIZE=1000, YSIZE=700

	win -> Refresh
	win.ygap = 3.5
	win -> Refresh, /DISABLE
	
	for icol = 1, 2 do begin
		for irow = 1, 3 do begin
			gfx = win -> FindByColRow([icol,irow])
			gfx.xtitle = ''
		endfor
	endfor
	
	sims = ['By0', 'By0.1', 'By1', '3D-By1']
	for i = 0, n_elements(sims)-1 do begin
		pos = win[sims[i] + ' Cut BL'].POSITION
		win[sims[i] + ' Cut ni'].POSITION=pos
		win[sims[i] + ' Cut UiL'].POSITION=pos
		win[sims[i] + ' Cut EN'].POSITION=pos
	endfor

	;Make the sun at the right
	;   - Negate independent data to make earth in -N direction
	;   - Reverse x-range for +N (sun) to the right and -N (earth) to left
	;   - Negate N- and M-components
	graphics = win -> Get(/ALL, ISA='MrPlot')
	foreach gfx, graphics do begin
		gfx -> GetData,  x, y
		
		;Negate N- & M-components
		;   - UiL will have units of 1e-3 (update MrAxis)
		;   - EN will have XTICKFORMAT='(f0.2)' (update MrAxis)
		if stregex(gfx.NAME, 'EN', /BOOLEAN) then begin
			gfx -> SetData, -x, -y
;			gfx.yrange = -reverse(gfx.yrange)
		endif else if stregex(gfx.NAME, 'UiL', /BOOLEAN) then begin
			gfx -> SetData, -x, y*1e2
;			gfx.yrange = gfx.yrange * 1e3
		endif else begin
			gfx -> SetData, -x,  y
		endelse
		
		gfx.TITLE         = ''
		gfx.XRANGE        = -zrange
		gfx.XTICKINTERVAL = 1.0
	endforeach

	;The images are of J_eM, so must be negated
	;   - Instead of multiplying the image by -1, we multiply the colorbar
	;     range by -1. That way the colors in the colorbar
	;     do not have to be reversed, just the axis range, to acheive the
	;     correct color scheme
	graphics = win -> Get(/ALL, ISA='MrImage')
	foreach gfx, graphics do begin
		gfx -> GetData, im, x,  y
		gfx -> SetData, im, x, -y
		
		title = gfx.TITLE
		title = (stregex(gfx.name, '^3D', /BOOLEAN) ? '3D ' : '2D ') + title

		gfx -> SetProperty, YRANGE=[-1, 1], YTICKINTERVAL=1.0, TITLE=title
	endforeach
	
	;Contours are of A_M and must be negated
	graphics = win -> Get(/ALL, ISA='MrContour')
	foreach gfx, graphics do begin
		;Get Data
		gfx  -> GetData,  z, x,  y
		name  = gfx.NAME
		gfx -> SetData, -z, x, -y
		
		;Set contour levels and ranges
		gfx -> GetProperty, LEVELS=levels, YRANGE=yrange, RANGE=range
		gfx -> SetProperty, LEVELS=-reverse(levels), YRANGE=[-1,1], RANGE=-reverse(range)
	endforeach
	
	;Colorbars
	graphics = win -> Get(/ALL, ISA='MrColorbar')
	foreach gfx, graphics do begin
		gfx.title  = 'J$\downM$'
		gfx.range  = -reverse(gfx.range)
		gfx.border = 1
	endforeach
	
	;Lines
	graphics = win -> Get(/ALL, ISA='MrPlotS')
	foreach gfx, graphics do begin
		if ~stregex(gfx.name, 'XPt', /BOOLEAN) then begin
			gfx -> GetProperty, XCOORDS=x, YCOORDS=[-0.25,0.5]
			gfx -> SetProperty, XCOORDS=-reverse(x)
		endif else begin
			gfx -> SetProperty, YCOORDS=[-1,1]
		endelse
	endforeach

	;Axes
	graphics = win -> Get(/ALL, ISA='MrAxis')
	foreach gfx, graphics do begin
		if stregex(gfx.name, 'EN', /BOOLEAN) then begin
			gfx.TICKFORMAT = '(f0.2)'
;			gfx.AXIS_RANGE = -gfx.AXIS_RANGE
		endif else if stregex(gfx.name, 'UiL', /BOOLEAN) then begin
			gfx.title = 'V$\downiL$'
		endif
	endforeach

;-----------------------------------------------------
; Individual Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;EN
	win['By0 Cut EN']    -> SetProperty, YRANGE=[-0.02, 0.03]
	win['By0.1 Cut EN']  -> SetProperty, YRANGE=[-0.02, 0.03]
	win['By1 Cut EN']    -> SetProperty, YRANGE=[-0.02, 0.03]
	win['3D-By1 Cut EN'] -> SetProperty, YRANGE=[-0.02, 0.03]

	;EN Axis
	win['By0 EN']    -> SetProperty, TICKINTERVAL=0.02, TICKLEN=0.04
	win['By0.1 EN']  -> SetProperty, TICKINTERVAL=0.02, TICKLEN=0.04
	win['By1 EN']    -> SetProperty, TICKINTERVAL=0.02, TICKLEN=0.04
	win['3D-By1 EN'] -> SetProperty, TICKINTERVAL=0.02, TICKLEN=0.04
	
	;NI
	win['By0 Cut ni']    -> SetProperty, YRANGE=[0.0, 1.0]
	win['By0.1 Cut ni']  -> SetProperty, YRANGE=[0.0, 1.0]
	win['By1 Cut ni']    -> SetProperty, YRANGE=[0.0, 1.0]
	win['3D-By1 Cut ni'] -> SetProperty, YRANGE=[0.0, 1.0]
	
	;NI Axis
	win['By0 ni']    -> SetProperty, TICKINTERVAL=0.5
	win['By0.1 ni']  -> SetProperty, TICKINTERVAL=0.5
	win['By1 ni']    -> SetProperty, TICKINTERVAL=0.5
	win['3D-By1 ni'] -> SetProperty, TICKINTERVAL=0.5
	
	;UiL
	win['By0 Cut UiL']    -> SetProperty, YRANGE=[-3,0]
	win['By0.1 Cut UiL']  -> SetProperty, YRANGE=[-3,0]
	win['By1 Cut UiL']    -> SetProperty, YRANGE=[-3,0]
	win['3D-By1 Cut UiL'] -> SetProperty, YRANGE=[-3,0]
	
	;UiL Axis
	win['By0 UiL']    -> SetProperty, TICKINTERVAL=1
	win['By0.1 UiL']  -> SetProperty, TICKINTERVAL=1
	win['By1 UiL']    -> SetProperty, TICKINTERVAL=1
	win['3D-By1 UiL'] -> SetProperty, TICKINTERVAL=1
	
	;BL
	win['By0 Cut BL']    -> SetProperty, YRANGE=[-0.25,0.5], YTICKINTERVAL=0.2
	win['By0.1 Cut BL']  -> SetProperty, YRANGE=[-0.25,0.5], YTICKINTERVAL=0.2
	win['By1 Cut BL']    -> SetProperty, YRANGE=[-0.25,0.5], YTICKINTERVAL=0.2
	win['3D-By1 Cut BL'] -> SetProperty, YRANGE=[-0.25,0.5], YTICKINTERVAL=0.2

	;JEY 3D
	win['3D-By1 Color JM']  -> SetProperty, YRANGE=[-5,5], RANGE=[-0.03,0.12], YTICKINTERVAL=2.5
	win['3D-By1 XPt']        -> SetProperty, YCOORDS=[-5,5]
	win['3D-By1 EN=0 Horiz'] -> SetProperty, XCOORDS=[1,-2]
	win['3D-By1 BL=0 Horiz'] -> SetProperty, XCOORDS=[1,-2]

	;Vertical Lines
	win['By0 BL=0 Vert']    -> SetProperty, YCOORDS=win['By0 Cut BL'].YRANGE
	win['By0.1 BL=0 Vert']  -> SetProperty, YCOORDS=win['By0.1 Cut BL'].YRANGE
	win['By1 BL=0 Vert']    -> SetProperty, YCOORDS=win['By1 Cut BL'].YRANGE
	win['3D-By1 BL=0 Vert'] -> SetProperty, YCOORDS=win['3D-By1 Cut BL'].YRANGE

;-----------------------------------------------------
; Return \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	fnames = 'proxfigs_bevn'
	
	;Return the window
	win -> Refresh
	return, win
end


;+
;   Plots of By0, By0.1, By1
;       - EN overplotted with contributing Hall terms.
;-
function MrProxFigs_3D_BEVn, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win)  then obj_destroy, win
		if obj_valid(oSim) then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

;---------------------------------------------------------------------
; 3D Sim, t=120816, y=1400 ///////////////////////////////////////////
;---------------------------------------------------------------------
	theSim       = 'Asymm-3D'
	time         = 120816
	ycell        = 860
	coord_system = 'SIMULATION'
	mva_frame    = 1
	ion_scale    = 1
	
	if coord_system eq 'SIMULATION' then begin
		horizontal = 0
		xrange     = 45.3 + [-15, 15]
		zrange     = [-4, 5]
	endif else begin
		horizontal = 1
		xrange     = [5, -4]
		zrange     = 45.3 + [-15, 15]
	endelse
	xpt = [45.3,0]
	
	;Show a larger domain size to fit in the three islands.
	xrange = 45.3 + [-15,15]
	zrange = [-1,2]
	
	;Create the simulation
	oSim   = MrSim_Create(theSim, time, $
	                      ION_SCALE = ion_scale, $
	                      XRANGE = xrange, $
	                      ZRANGE = zrange, $
	                      COORD_SYSTEM = coord_system, $
	                      MVA_FRAME = mva_frame)
	ycoord = oSim -> Cell2Coord(ycell, /Y)
	oSim.yrange = [ycoord, ycoord]

;---------------------------------------------------------------------
; Create the Plot ////////////////////////////////////////////////////
;---------------------------------------------------------------------

	;Create the color image
	win = mrsim_colorslab(osim, 'JM')
	win -> Refresh, /DISABLE
	win -> SetProperty, LAYOUT = [2,1], XSIZE=1000, OXMARGIN=[8,15], XGAP=29
	coord = horizontal ? xpt[1] : xpt[0]

	;Plot cuts of BL, n, and EN through the X-line
	lc1 = MrSim_LineCut(oSim, 'BL', coord, /CURRENT, HORIZONTAL=horizontal)
	pos = lc1.position
	lc2 = MrSim_LineCut(oSim, 'ni',  coord, /CURRENT, HORIZONTAL=horizontal)
	lc3 = MrSim_LineCut(oSim, 'EN',  coord, /CURRENT, HORIZONTAL=horizontal)
	lc4 = MrSim_LineCut(oSim, 'UiL', coord, /CURRENT, HORIZONTAL=horizontal)
	
	;Move ni and Ex plots onto Bz to create independent data spaces for them
	lc1 -> SetProperty, YSTYLE=9
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

	;Zoom in on BL transition
;	graphics = win -> Get(/ALL, ISA='MrPlot')
;	foreach gfx, graphics do gfx.xrange = [-1,2]

;---------------------------------------------------------------------
; Add Lines //////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	;Find where Bz = 0
	lc1   -> GetData, x, Bz
	!Null  = min(Bz, iZero, /ABSOLUTE)

	;Draw a vertical line throubh BL = 0
	!Null = MrPlotS(x[[iZero, iZero]], lc1.yrange, $
	                /CLIP, $
	                COLOR     = 'Black', $
	                LINESTYLE = 'Dash', $
	                NAME      = 'BL=0 Vert', $
	                TARGET    = lc1)
	
	;Draw a horizontal line throubh BL = 0
	!Null = MrPlotS(lc1.xrange, [0, 0], $
	                /CLIP, $
	                COLOR     = 'Black', $
	                LINESTYLE = 'Dash', $
	                NAME      = 'BL=0 Horiz', $
	                TARGET    = lc1)
	
	;Draw a horizontal line throubh EN = 0
	!Null = MrPlotS(lc3.xrange, [0, 0], $
	                /CLIP, $
	                COLOR     = 'Red', $
	                LINESTYLE = 'Dash', $
	                NAME      = 'EN=0 Horiz', $
	                TARGET    = lc3)

	;Line through the X-point
	!Null = MrPlotS( [xpt[0],xpt[0]], win['Color JM'].yrange, $
	                 COLOR     = 'White', $
	                 LINESTYLE = 'Dash', $
	                 NAME      = 'XPt', $
	                 TARGET    = win['Color JM'])

;---------------------------------------------------------------------
; Make Pretty ////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	
	;Add guidefield strength to title
	title = win['Color JM'].TITLE
	win['Color JM'].TITLE = 'B$\downG$=1 ' + title
	
	;Decrease CB width
	win['CB: Color JM'].width = 0.5

	win -> Refresh
	return, win
end


;+
;   Plots of By0, By0.1, By1
;       - EN overplotted with contributing Hall terms.
;-
function MrProxFigs_ETot_Terms, $
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
; Asymm-Scan/By0, t=30 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0'
	time      = 28
	xrange    = 36.77 + [-5, 5]
	zrange    = [-1,2]
	coord_sys = 'Simulation'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;Figure out where the x-point is
	xpt = MrSim_XPoint(oSim)
	
	;Create the figure
	win = MrSim_OhmsLaw(oSim, 'Z', xpt[0], /TOTAL)
	win -> Refresh, /DISABLE
	
	;Rename all of the graphics
	graphics = win -> Get(/ALL)
	foreach gfx, graphics do begin
		name = gfx.name
		gfx.name = 'By0 ' + name
	endforeach
	
	;Change title
	title = win['By0 Total EZ'].TITLE
	win['By0 Total EZ'].TITLE = 'B$\downG$=0 ' + title
	
	;Draw line at X-point
	gET  = win['By0 Total Ez']
	line = MrPlotS( [xpt[1], xpt[1]], gET.YRANGE, $
	                COLOR     = 'Grey', $
	                LINESTYLE = '--', $
	                NAME      = 'By0 BL=0', $
	                TARGET    = gET)
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By0.1, t=28 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By0.1'
	time      = 28
	xrange    = 40.5 + [-2.0, 2.0]
	zrange    = [-1,2]
	coord_sys = 'Simulation'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;Figure out where the x-point is
	xpt = MrSim_XPoint(oSim)
	
	;Create the figure
	!Null = MrSim_OhmsLaw(oSim, 'Z', xpt[0], /TOTAL, /CURRENT)
	
	;Rename all of the graphics
	graphics = win -> Get(/ALL)
	foreach gfx, graphics do begin
		name = gfx.name
		if ~stregex(name, '^By', /BOOLEAN) then gfx.name = 'By0.1 ' + name
	endforeach
	
	;Change title
	title = win['By0.1 Total EZ'].TITLE
	win['By0.1 Total EZ'].TITLE = 'B$\downG$=0.1 ' + title
	
	;Draw line at X-point
	gET  = win['By0.1 Total Ez']
	line = MrPlotS( [xpt[1], xpt[1]], gET.YRANGE, $
	                COLOR     = 'Grey', $
	                LINESTYLE = '--', $
	                NAME      = 'By0.1 BL=0', $
	                TARGET    = gET)
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-Scan/By1, t=30 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	theSim    = 'Asymm-Scan/By1'
	time      = 28
	xrange    = [42.0, 45.0]
	zrange    = [-1,2]
	coord_sys = 'Simulation'
	ion_scale = 1
	mva_frame = 1
	oSim      = MrSim_Create(theSim, time, $
	                         COORD_SYSTEM = coord_sys, $
	                         MVA_FRAME    = mva_frame, $
	                         ION_SCALE    = ion_scale, $
	                         XRANGE       = xrange, $
	                         ZRANGE       = zrange)
	
	;Figure out where the x-point is
	xpt = MrSim_XPoint(oSim)
	
	;Create the figure
	!Null = MrSim_OhmsLaw(oSim, 'Z', xpt[0], /TOTAL, /CURRENT)
	
	;Rename all of the graphics
	graphics = win -> Get(/ALL)
	foreach gfx, graphics do begin
		name = gfx.name
		if ~stregex(name, '^By', /BOOLEAN) then gfx.name = 'By1 ' + name
	endforeach
	
	;Change title
	title = win['By1 Total EZ'].TITLE
	win['By1 Total EZ'].TITLE = 'B$\downG$=1 ' + title
	
	;Draw line at X-point
	gET  = win['By1 Total Ez']
	line = MrPlotS( [xpt[1], xpt[1]], gET.YRANGE, $
	                COLOR     = 'Grey', $
	                LINESTYLE = '--', $
	                NAME      = 'By1 BL=0', $
	                TARGET    = gET)
	
	;Destroy the simulation object
	obj_destroy, oSim

;-----------------------------------------------------
; Asymm-3D/By1, t=120816 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
;	theSim    = 'Asymm-3D'
;	time      = 120816
;	ycell     = 860
;	coord_sys = 'SIMULATION'
;	ion_scale = 1
;	mva_frame = 1
;	
;	if coord_sys eq 'SIMULATION' then begin
;		horizontal = 0
;		xrange     = 45.3 + [-15, 15]
;		zrange     = [-1,2]
;	endif else begin
;		horizontal = 1
;		xrange     = [2,-1]
;		zrange     = 45.3 + [-15, 15]
;	endelse
;	xpt = [45.3,0]
;	
;	oSim = MrSim_Create(theSim, time, $
;	                    COORD_SYSTEM = coord_sys, $
;	                    MVA_FRAME    = mva_frame, $
;	                    ION_SCALE    = ion_scale, $
;	                    XRANGE       = xrange, $
;	                    ZRANGE       = zrange)
;	ycoord = oSim -> Cell2Coord(ycell, /Y)
;	oSim.yrange = [ycoord, ycoord]
;	
;
;	;Create the figure
;	!Null = MrSim_OhmsLaw(oSim, 'Z', xpt[0], /TOTAL, /CURRENT)
;
;	;Rename all of the graphics
;	graphics = win -> Get(/ALL)
;	foreach gfx, graphics do begin
;		name = gfx.name
;		if ~stregex(name, '^By', /BOOLEAN) then gfx.name = 'D3-By1 ' + name
;	endforeach
;	
;	;Change title
;	title = win['D3-By1 Total EZ'].TITLE
;	win['D3-By1 Total EZ'] -> SetProperty, TITLE='B$\downG$=1 ' + title
;	
;	;Draw line at X-point
;	gET  = win['D3-By1 Total Ez']
;	line = MrPlotS( [0.4967, 0.4967], gET.YRANGE, $
;	                COLOR     = 'Grey', $
;	                LINESTYLE = '--', $
;	                NAME      = '3D-By1 BL=0', $
;	                TARGET    = gET)
;	
;	;Destroy the simulation object
;	obj_destroy, oSim

;-----------------------------------------------------
; Make Pretty \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Make spacing larger
	win -> SetProperty, YGAP=0.0, OXMARGIN=[10,18]

	;Remove and move legends
	win -> Remove, win["By0.1 Ohm's Law"]
	win -> Remove, win["By1 Ohm's Law"]
;	win -> Remove, win["D3-By1 Ohm's Law"]
	win["By0 Ohm's Law"] -> SetProperty, ALIGNMENT='NW', LINESTYLE='None', FILL_COLOR=''

	;global y-range
	yrange = [-0.03, 0.06]

	;Make the sun at the right
	;   - Must negate all N-components to point the correct direction
	graphics = win -> Get(/ALL, ISA='MrPlot')
	foreach gfx, graphics do begin
		gfx -> GetProperty, XRANGE=xrange
		gfx -> GetData,  x,  y
		gfx -> SetData, -x, -y
		gfx -> SetProperty, XRANGE=-xrange, YRANGE=yrange
		
		;Reverse x-range for +N (sun) to the right and -N (earth) to left
;		gfx -> GetProperty, XRANGE=xrange, YRANGE=yrange
;		gfx -> SetProperty, XRANGE=-zrange, YRANGE=-reverse(yrange)
	endforeach

	;Vertical lines
	graphics = win -> Get(/ALL, ISA='MrPlotS')
	foreach gfx, graphics do begin
		gfx -> GetProperty, XCOORDS=xcoords
		gfx -> SetProperty, XCOORDS=-xcoords, YCOORDS=yrange
	endforeach

	;Remove title and ticks from intermediate plots
	win['By0 Total EZ']   -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
	win['By0.1 Total EZ'] -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
;	win['By1 Total EZ']   -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'

;---------------------------------------------------------------------
; Add B_G=X As A Legend //////////////////////////////////////////////
;---------------------------------------------------------------------
	bg = ['By0', 'By0.1', 'By1'];, 'D3-By1']

	for i = 0, n_elements(bg) - 1 do begin
		gfx       = win[bg[i] + ' Total EZ']
		title     = gfx.TITLE
		gfx.TITLE = ''
	
		;Title Segments
		parts = strsplit(title, ' ', /EXTRACT)
		
		;Create a Legend
		!Null = MrLegend(ALIGNMENT    ='NW', $
		                 FILL_COLOR   = '', $
		                 NAME         = bg[i], $
		                 POSITION     = [0,1], $
		                 /RELATIVE, $
		                 LABEL        = parts[0], $
		                 LINESTYLE    = 6, $
		                 SAMPLE_WIDTH = 0.0, $
		                 TARGET       = gfx)
		
		;Draw a line through EN = 0
		line = MrPlotS( gfx.XRANGE, [0,0], $
		                COLOR     = 'Black', $
		                LINESTYLE = '--', $
		                NAME      = bg[i] + ' EN=0', $
		                TARGET    = gfx)
	endfor

;-----------------------------------------------------
; Return \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	win -> SetProperty, XSIZE=450, YSIZE=450
	win -> SetGlobal, YTITLE='E$\downN$'
	fnames = 'proxfigs_E-total'
	
	;Return the window
	win -> Refresh
	return, win
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
		'FIGURE 1':   win = MrProxFigs_Figure1(FNAMES=fnames)
		'FIGURE 2':   win = MrProxFigs_Figure2(FNAMES=fnames)
		'OHMSLAW':    win = MrProxFigs_OhmsLaw(FNAMES=fnames)
		'3D BEVN':    win = MrProxFigs_3D_BEVn(FNAMES=fnames)
		'XLINE BEVN': win = MrProxFigs_XLine_BEVn(FNAMES=fnames)
		'BY0-MRX':    win = MrProxFigs_XLineProxBy0(FNAMES=fnames)
		'ETOTAL':     win = MrProxFigs_ETot_Terms(FNAMES=fnames)
		'X-POINT':    win = MrProxFigs_XPoint(FNAMES=fnames)
		'TPROX':      win = MrProxFigs_TProx(FNAMES=fnames)
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
	
			;Take a snapshot
			if png then begin
				;Set a true type font
				device, SET_FONT='Helvetica Bold', /TT_FONT
				win -> SetGlobal, FONT=1
			
				win.SAVEAS -> SetProperty, IM_RASTER=0
				win -> Save, fbase + '-ss.png'
				win.SAVEAS -> SetProperty, IM_RASTER=1
			endif
	
			;
			; PostScript
			;
			
			;Hardware fonts
			win -> SetGlobal, FONT=0;, CHARSIZE=1.0
	
			;Save a variety of file types.
			win -> Refresh
			if im_png then win -> Save, fbase + '_im.png'
			if eps    then win -> Save, fbase + '.eps'
			if ps     then win -> Save, fbase + '.ps'
		
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