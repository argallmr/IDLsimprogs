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
;-
;*****************************************************************************************
;+
;   Create Figure 1: X-Line proximity for Asymm-Scan/By0, Asymm-Scan/By1, and Asymm-3D.
;-
function ProxFigs_Figure1, $
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
    win1 = MrSim_XProximity(theSim, cuts, tIndex, $
                            C_NAME    = c_name, $
                            COORD_SYS = coord_sys, $
                            IM_NAME   = im_name, $
                            ION_SCALE = ion_scale, $
                            MVA_FRAME = mva_frame, $
                            XRANGE    = xrange, $
                            ZRANGE    = zrange)

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
	                   ['        ', '']]

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
		'FIGURE 1': win = ProxFigs_Figure1(FNAMES=fnames)
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