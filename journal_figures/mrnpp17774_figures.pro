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
;   Create Figure 2: eMap for Asymm-Scan/By0 and Asymm-Scan/By1.
;-
function MrNPP17774_Figure1, $
FNAMES=fnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if obj_valid(win1) then obj_destroy, win1
		if obj_valid(oSim) then obj_destroy, oSim
		void = cgErrorMSG()
		return, obj_new()
	endif

    ;Sym distributions
    edist_sym  = [[847.85, -0.88], $
                  [847.85, -0.025], $
                  [847.85,  0.98]]
    
    ;Asym distributions
    edist_asym = [[367.62, -1.79], $
                  [367.62, -0.908], $
                  [367.62, 1.53]]

;-----------------------------------------------------
; Sim1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim      = 'Sim1'
	time        = 76
	coord_sys   = 'Simulation'
	c_name      = 'Ay'
	xrange      = [820,880]
	zrange      = [-8, 8]

	;Create the simulation object.
	oSim = MrSim_Create(theSim, time, /BINARY, XRANGE=xrange, ZRANGE=zrange)
    
    ;Create the graphics windows
    win1 = MrWindow(OXMARGIN=[10,10], XSIZE=800, YSIZE=500, YGAP=0.5, XGAP=15, REFRESH=0)
    
    ;Plot of ni, Jey, Dng
	!Null = MrSim_ColorSlab(oSim, 'ni',    C_NAME=c_name, /CURRENT)
    !Null = MrSim_ColorSlab(oSim, 'Jey',   C_NAME=c_name, /CURRENT)
    !Null = MrSim_ColorSlab(oSim, 'Dng_e', C_NAME=c_name, /CURRENT)
    
    ;Destroy the simulation object
    obj_destroy, osim
    
    ;Rename graphics
    allGfx = win1 -> Get(/ALL)
    foreach gfx, allGfx do begin
        name = gfx.NAME
        gfx.NAME = 'Sym ' + name
    endforeach

    ;Set properties
    title = win1['Sym Color ni'].TITLE
    win1['Sym Color ni']        -> SetProperty, XTICKFORMAT='(a1)', XTITLE='', TITLE='Symmetric ' + title
    win1['Sym Color Jey']       -> SetProperty, XTICKFORMAT='(a1)', XTITLE='', TITLE=''
    win1['Sym Color Dng_e']     -> SetProperty, TITLE=''
    win1['Sym CB: Color ni']    -> SetProperty, WIDTH=1
    win1['Sym CB: Color Jey']   -> SetProperty, WIDTH=1
    win1['Sym CB: Color Dng_e'] -> SetProperty, WIDTH=1, TITLE='Dng'

;-----------------------------------------------------
; Asymm-Scan/By0 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim      = 'Asymm-Scan/By0'
	time        = 28
	c_name      = 'Ay'
	xrange      = 367.7 + [-20, 20]
	zrange      = [-8, 5]

	;Create the simulation object.
	oSim = MrSim_Create(theSim, time, /BINARY, XRANGE=xrange, ZRANGE=zrange)

	;Plot ni, Jey, Dng.
	!Null = MrSim_ColorSlab(oSim, 'ni',    C_NAME=c_name, /CURRENT)
    !Null = MrSim_ColorSlab(oSim, 'Jey',   C_NAME=c_name, /CURRENT)
    !Null = MrSim_ColorSlab(oSim, 'Dng_e', C_NAME=c_name, /CURRENT)
    
    ;Destroy the simulation object
    obj_destroy, osim

    ;Move graphics to column 2
    win1['Color ni']    -> SetLayout, [2,1]
    win1['Color Jey']   -> SetLayout, [2,2]
    win1['Color Dng_e'] -> SetLayout, [2,3]
    win1 -> TrimLayout
    
    ;Set properties
    title = win1['Color ni'].TITLE
    win1['Color ni']        -> SetProperty, RANGE=[0.56, 0.8], XTICKFORMAT='(a1)', XTITLE='', YTITLE='', TITLE='Asymmetric ' + title
    win1['Color Jey']       -> SetProperty, XTICKFORMAT='(a1)', XTITLE='', YTITLE='', TITLE=''
    win1['Color Dng_e']     -> SetProperty, YTITLE='', TITLE=''
    win1['CB: Color ni']    -> SetProperty, RANGE=[0.56, 0.8], WIDTH=1
    win1['CB: Color Jey']   -> SetProperty, WIDTH=1
    win1['CB: Color Dng_e'] -> SetProperty, WIDTH=1, TITLE='Dng'
	
	win1 -> Refresh
	return, win1
end


;+
;   Create Figure 2: eMap for Asymm-Scan/By0 and Asymm-Scan/By1.
;-
function MrNPP17774_Figure2, $
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

    ;Asym distributions
    edist_asym = [[367.62, -1.79], $
                  [367.62, -0.908], $
                  [367.62, 1.53]]

;-----------------------------------------------------
; Sim1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim      = 'Sim1'
	time        = 76
	coord_sys   = 'Simulation'
	xrange      = [820,880]
	zrange      = [-8, 8]

	;Create the simulation object.
	oSim = MrSim_Create(theSim, time, /BINARY, XRANGE=xrange, ZRANGE=zrange)
    
    ;Describe eMap
    type       = 'Vx-Vy'
    bin_center = [[847.85, -0.88], $
                  [847.85, -0.025], $
                  [847.85,  0.98]]
    half_width = [0.5]
    layout     = [1,3]
    v_va       = 1
    xsize      = 300
    ysize      = 400
    
    ;Create the eMap
    win = MrSim_eMap(oSim, type, bin_center, half_width, time, $
                     NBINS=nBins, $
                     LAYOUT=layout, $
                     V_VA=v_va, $
                     XSIZE=xsize, $
                     YSIZE=ysize)
    
    ;Destroy the simulation object
    obj_destroy, oSim
    
    ;Set window properties
    win -> Refresh, /DISABLE
    win -> SetProperty, XGAP=20, XSIZE=600
    
    ;Rename everything
    allGfx = win -> Get(/ALL)
    foreach gfx, allGfx do begin
        name = gfx.NAME
        gfx.NAME = 'Sym ' + name
    endforeach
    
    ;Set Properties
    win['Sym eMap Title'] -> GetProperty, STRING=title, XLOC=xloc
    title = strmid(title, strpos(title, 't*\Omega'))
    win['Sym eMap Title'] -> SetProperty, STRING='Symmetric ' + title, XLOC=0.25
    
    pos = win['Sym CB: Counts'].POSITION
    win['Sym CB: Counts'] -> SetProperty, WIDTH=1.5, POSITION=[0.42, pos[1], 0.43, pos[3]]
    
;-----------------------------------------------------
; Asymm-Scan/By0 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Descibe the simulation
	theSim      = 'Asymm-Scan/By0'
	time        = 28
	xrange      = 367.7 + [-20, 20]
	zrange      = [-8, 5]

	;Create the simulation object.
	oSim = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange)
    
    ;Describe eMap
    type       = 'Vx-Vy'
    bin_center = [[367.62, -1.79], $
                  [367.62, -0.908], $
                  [367.62, 1.53]]
    half_width = [0.5]
    layout     = [1,3]
    v_va       = 1
    xsize      = 300
    ysize      = 400
    
    ;Create the eMap
    win1 = MrSim_eMap(oSim, type, bin_center, half_width, time, $
                     NBINS=nBins, $
                     LAYOUT=layout, $
                     V_VA=v_va, $
                     XSIZE=xsize, $
                     YSIZE=ysize)
    win1 -> Refresh, /DISABLE
    
    ;Destroy the simulation object
    obj_destroy, oSim
    
    ;Move distributions to the second column, then to the other window
    allIm = win1 -> Get(/ALL, ISA='MrImage')
    foreach gfx, allIm do begin
        layout = gfx.LAYOUT
        colrow = win -> ConvertLocation(layout[2], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
        gfx -> SwitchWindows, win
    endforeach
    
    ;Move all remaining objects into the first window
    allGfx = win1 -> Get(/ALL)
    foreach gfx, allGfx do gfx -> SwitchWindows, win
    
    ;Set Properties
    win['eMap Title']     -> GetProperty, STRING=title, XLOC=xloc
    title = strmid(title, strpos(title, 't*\Omega'))
    win['eMap Title']     -> SetProperty, STRING='Asymmetric ' + title, XLOC=0.75
    win['Sym CB: Counts'] -> SetProperty, WIDTH=1.5
    
    ;Destroy the unused window
    obj_destroy, win1
	
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
function MrNPP17774_Figures, figure, $
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
	list_of_figures = [['Figure 1', 'Symmetric & Asymmetric n, Jey, Ng'], $
	                   ['Figure 2', 'Symmetric & Asymmetric triangular eDist']]

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
		'FIGURE 1': win = MrNPP17774_Figure1(FNAMES=fnames)
		'FIGURE 2': win = MrNPP17774_Figure2(FNAMES=fnames)
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
			fname = 'MrNPP17774_' + idl_validname(figure, /CONVERT_ALL)
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
				fname = FilePath('MrNPP17774_' + fnames[i], ROOT_DIR=froot)
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