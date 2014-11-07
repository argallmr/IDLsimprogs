; docformat = 'rst'
;
; NAME:
;    Figures_XLine_Proximity
;
; PURPOSE:
;+
;   Create figures in preparation of my X-Line proximity paper..
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
;       2014/03/15  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Create Figure 2: 2D simulation.
;-
function Prox_Figure1
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(Fig1)   then obj_destroy, Fig1
        if obj_valid(difwin) then obj_destroy, difwin
        if obj_valid(cwin)   then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 2.0
    col_width = [0.2, 0.3, 0.2, 0.3]
    layout    = [4,5]
    oymargin  = [4,4]
    xgap      = [10, 7, 11]
    xsize     = 1200
    ygap      = 0.5
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim       = 'Asymm-Large-2D'
    time         = 32
    xrange       = [2, -2]
    zrange       = 150.9 + [-3.0, 3.0]
    coord_system = 'Magnetopause'
    ion_scale    = 1
    oSim = MrSim_Create(theSim, time, $
                        ION_SCALE    = ion_scale, $
                        MVA_FRAME    = mva_frame, $
                        XRANGE       = xrange, $
                        ZRANGE       = zrange, $
                        COORD_SYSTEM = coord_system)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create cuts of Bx, By, ni, Uix, Ez
    Fig1 = MrSim_XProximity([150.9, 150.3, 148.4], SIM_OBJECT=oSim)
    Fig1 -> Refresh, /DISABLE

    ;Format the y-axes
    Fig1['Cut Bx']  -> SetProperty, YTICKINTERVAL=0.25
    Fig1['Cut ni']  -> SetProperty, YTICKINTERVAL=0.5,  YTITLE='n$\downi$', YRANGE=[0,1.2]
    Fig1['Cut Uix'] -> SetProperty, YTICKINTERVAL=0.01, YTITLE='U$\downiL$'
    Fig1['Cut Ez']  -> SetProperty, YTICKINTERVAL=0.01

    ;Set the layout
    ;   Move plots to the second column.
    Fig1 -> SetProperty, XSIZE=1200, YSIZE=500
    for i = 0, 4 do Fig1[i] -> SetLayout, [2,i+1]

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Add a 2D color plot
    ;   - Will automatically be added to location [1,1]
    !Null = MrSim_ColorSlab(oSim, 'Jey', C_NAME='Ay', /CURRENT, $
                            HORIZ_LINE=[150.9, 150.3, 148.4])

    ;Format the X-axis
    Fig1['Color Jey']  -> SetProperty, TITLE='', XTICKS=2, XTICKFORMAT='(f4.1)'
    Fig1['HLines Jey'] -> SetProperty, THICK=3, LINESTYLE=2, COLOR = ['Cyan', 'Green', 'Red']
    
    ;Move the colorbar to the bottom of the figure
    Fig1['CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', OFFSET=0.5, WIDTH=1, $
                                          XTICKS=2, FORMAT='(f5.2)'
    
    ;Change the names of all figures be pre-pending a "t=32"
    graphics = Fig1 -> Get(/ALL)
    foreach gfx, graphics do begin
        name = gfx -> GetName()
        gfx.NAME = 't=32 ' + name
    endforeach

    ;Destroy the simulation object
    obj_destroy, oSim
    
;---------------------------------------------------------------------
; 3D Sim, t=108090, y=905 ////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim       = 'Asymm-3D'
    time         = 108090
    xrange       = [5, -4]
    yslice       = 905
    zrange       = 45.3 + [-15.0, 15.0]
    coord_system = 'Magnetopause'
    ion_scale    = 1
    oSim = MrSim_Create(theSim, time, yslice, $
                        ION_SCALE    = ion_scale, $
                        MVA_FRAME    = mva_frame, $
                        XRANGE       = xrange, $
                        ZRANGE       = zrange, $
                        COORD_SYSTEM = coord_system)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    lines  = [47.0, 43.2, 37.0]
    difwin = MrSim_XProximity(lines, SIM_OBJECT=oSim)
    difwin -> Refresh, /DISABLE

    ;Change properties
    difwin['Cut Bx'] -> SetProperty, YTICKINTERVAL=0.3
    difwin['Cut By'].YTICKINTERVAL  = 0.1
    difwin['Cut Ez'].YTICKINTERVAL  = 0.02
    difwin['Cut Uix'].YTICKINTERVAL = 0.02
    
    ;We want to move the plots to window Fig3
    ;   - Resize to be similar to Fig3
    difwin -> SetProperty, XSIZE=xsize, YSIZE=ysize, XGAP=xgap, LAYOUT=layout
    
    ;Move the figures to the 4th column
    graphics = difwin -> Get(/ALL, ISA='MRPLOT', COUNT=nPlots)
    for i = 0, nPlots - 1 do begin
        ;Move to the 4th column
        thisPlot = graphics[i]
        thisPlot -> SetLayout, [4,i+1]
        
        ;Rename by prepending a "y=860"
        name          = thisPlot -> GetName()
        thisPlot.NAME = 'y=905 ' + name
        
        ;Switch windows
        thisPlot -> SwitchWindows, Fig1
    endfor

    ;Move the legend (last remaining item)
    theLegend      = difwin[0]
    name           = theLegend -> GetName()
    theLegend.NAME = 'y=905 ' + name
    theLegend     -> SwitchWindows, Fig1
    
    ;Destroy the empty DIFWIN window
    obj_destroy, difwin

;---------------------------------------------------------------------
; 2D Color Plot  /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Plot Jey
    cwin = MrSim_ColorSlab(oSim, 'Jey', HORIZ_LINE=lines)
    cwin -> Refresh, /DISABLE
    obj_destroy, oSim
    
    ;Move to location [1,2]. Split the title into two lines
    cwin['Color Jey'] -> SetLayout, [1,2]
    cwin['Color Jey'] -> SetProperty, TITLE='', XTICKINTERVAL=3
    
    ;Make the cut locations more visible.
    cwin['HLines Jey'] -> SetProperty, THICK=3, LINESTYLE=2, COLOR=['Cyan', 'Green', 'Red']
    
    ;Rename and move graphics into FIG2
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do begin
        name     = gfx -> GetName()
        gfx.NAME = 'y=905 ' + name
        gfx     -> SwitchWindows, Fig1
    endforeach

    ;Move the colorbar to the bottom of the figure
    Fig1['y=905 CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', OFFSET=0.5, WIDTH=1, $
                                                XTICKS=2, XTICKFORMAT='(f5.2)'

    ;Delete CWIN and refocus FIG2
    cwin -> Cleanup
    Fig1 -> SetCurrent
    
    ;Move the color plot to location [3,1]
    Fig1['y=905 Color Jey'] -> SetLayout, [3,1]

;---------------------------------------------------------------------
; Set Positions //////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Set the column widths
    Fig1 -> SetProperty, COL_WIDTH=col_width, OYMARGIN=oymargin, XGAP=xgap, YGAP=ygap
    Fig1 -> SetGlobal, CHARSIZE=charsize

    ;Set up the positions
    pos = MrLayout(layout, OYMARGIN=oymargin, XGAP=xgap, YGAP=ygap, $
                   COL_WIDTH=col_width, CHARSIZE=charsize)

    ;Change the position of the color image
    top      = pos[1,0] + (pos[3,0] - pos[1,0])*0.25
    position = [pos[0,16], pos[1,16], pos[2,4], top]
    Fig1['t=32 Color Jey'] -> SetProperty, POSITION=position

    ;Reposition to fill the middle three rows.
    top      = pos[1,2] + (pos[3, 2] - pos[1, 2])*0.25
    position = [pos[0,18], pos[1,18], pos[2,2], top]
    Fig1['y=905 Color Jey'].POSITION=position

    ;Adjust the legend
    ;t=32
    oTemp = Fig1['t=32 Legend: Cuts']
    oTemp.LOCATION=4
    oTemp -> GetProperty, TITLES=titles, LOCATION=location
    titles = strsplit(titles, 'L=', /EXTRACT)
    location -= [0, 0.015]
    oTemp -> SetProperty, TITLES=titles, LOCATION=location
    
    ;t=90
    oTemp = Fig1['y=905 Legend: Cuts']
    oTemp.LOCATION=4
    oTemp -> GetProperty, TITLES=titles, LOCATION=location
    titles = strsplit(titles, 'L=', /EXTRACT)
    location -= [0, 0.015]
    oTemp -> SetProperty, TITLES=titles, LOCATION=location
    
    ;Put "(a)" and "(b)" on the figure
    aText = MrText(0.03, 0.9, '(a)', CHARSIZE=2, /CURRENT, /NORMAL)
    bText = MrText(0.50, 0.9, '(b)', CHARSIZE=2, /CURRENT, /NORMAL)

    Fig1 -> Refresh
    return, Fig1
end


;+
;   Simulate an MMS flyby across the X-line.
;-
function Prox_Asymm3D_FlyBy
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(win)   then obj_destroy, win
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    time   = 108090
    ycell  = 905
    xrange = 460 + [-40, 40]
    zrange = [-20, 20]
    root   = '/data2/Asymm-3D/'
    oSim   = MrSim_Create('Asymm-3D', time, XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange, $
                          DIRECTORY=root, INFO_BINARY=root+'/data/info', INFO_ASCII=root+'info', $
                          /BINARY)
    ycoord      = oSim -> GetCoord(ycell, /Y)
    oSim.YRANGE = [ycoord, ycoord]

    r0          = [457, ycoord, -15]
    r1          = [457, ycoord,  15]
    name        = 'Bx'
    nDist       = 5
    dist_type   = 'Vperp1-Vperp2'
    dist_layout = [3,2]
    dist_size   = [1,1,1]
    im_name     = 'Dng_e'
    velocity    = 60.0
    win = MrSim_FlyBy( oSim, name, r0, r1, $
                       /ADD_LEGEND, $
                       IM_NAME='Dng_e', $
                       NDIST=nDist, $
                       DIST_LAYOUT=dist_layout, $
                       DIST_TYPE=dist_type, $
                       DIST_SIZE=dist_size, $
                       DIST_WIN=dist_win, $
                       VELOCITY=velocity)

    return, win
end


;+
;   eMap for the Asymm-3D simulation.
;-
function Prox_Asymm3D_eMap, $
FNAME=fname
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(im_win)   then obj_destroy, im_win
        if obj_valid(tempWin)  then obj_destroy, tempWin
        if max(obj_valid(win)) then obj_destroy, win
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Create the simulation object
    theSim = 'Asymm-3D'
    time   = 108090
    ycell  = 905
    xrange = 458.8 + [-100, 100]
    zrange = [-30, 30]
    root   = '/data2/Asymm-3D/'
    oSim   = MrSim_Create(theSim, time, XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange, $
                          /BINARY)
    
    ;Set the y-range
    ycoord = oSim -> GetCoord(ycell, /Y)
    oSim.yrange = [ycoord, ycoord]

    ;Create the distribution map    
    type       = ['Vx-Vy', 'Vx-Vz', 'Vy-Vz', 'Vpar-Vperp', 'Vpar-Vperp1', 'Vpar-Vperp2', 'Vperp1-Vperp2']
    bin_center = [458.8, ycoord, 2.9]
    half_width = [1,1,1]
    im_name    = 'Jey'
    hgap       = 10
    vgap       = 5
    layout     = [5,5]
    location   = 6
    
    nDist = n_elements(type)
    win   = objarr(nDist+1)
    
    ;Step through each type of distribution.
    for i = 0, n_elements(type)-1 do begin
        ;Create a 2D overview plot for only the first distribution type
        if i eq 0 $
            then imname = im_name $
            else imname = !Null
    
        ;Create the distribution
        tempWin = MrSim_eMap(oSim, type[i], bin_center, half_width, $
                             HGAP     = hgap, $
                             IM_NAME  = imname, $
                             IM_WIN   = im_win, $
                             LAYOUT   = layout, $
                             LOCATION = location, $
                             VGAP     = vgap)
        
        ;Store the window
        if i eq 0 then begin
            win[0] = im_win
            win[1] = tempWin
        endif else begin
            win[i+1] = tempWin
        endelse
    endfor

    ;Create an output file name
    fpart1 = string(FORMAT='(%"%s_t%i_y%i")', theSim, time, ycell)
    fpart2 = string(FORMAT='(%"x%iz%i_w%i_h%i_v%i_lay%ix%i_loc%i")', $
             bin_center[[0,2]], half_width[0], hgap, vgap, layout, location)
    fname = [fpart1 + '_' + im_name + '_'   + fpart2, $
             fpart1 + '_eMap_' + type + '_' + fpart2]

    return, win
end


;+
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function Figures_Proximity, figure, $
SAVE=tf_save
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if max(obj_valid(win)) then obj_destroy, win
        void = cgErrorMSG()
        return, obj_new()
    endif

    _figure = strupcase(figure)
    tf_save = keyword_set(tf_save)

    ;Create the figure    
    case _figure of
        'FIGURE 1':       win = Prox_Figure1()
        'ASYMM-3D EMAP':  win = Prox_Asymm3D_eMap(FNAME=fname)
        'ASYMM-3D FLYBY': win = Prox_Asymm3D_FlyBy()
        else: message, 'Figure "' + figure + '" not an option.', /INFORMATIONAL
    endcase

    ;Save the figure
    if tf_save then begin
        directory = '/home/argall/figures/Asymm-3D/eMap-Scan/'
        if n_elements(fname) eq 0 then fname = _figure
        
        nPlots = n_elements(win)
        
        ;Single plot
        if nPlots eq 1 then begin
            ;Take a snapshot
            saveas = win.SaveAs
            saveas -> SetProperty, IM_RASTER=0
            win -> Save, filepath(fname + '.png', ROOT_DIR=directory)
            
        ;Multiple plots
        endif else begin
            for i = 0, nPlots-1 do begin
                ;Take a snapshot
                saveas = win[i].SaveAs
                saveas -> SetProperty, IM_RASTER=0
                win[i] -> Save, filepath(fname[i] + '.png', ROOT_DIR=directory)
            endfor
        endelse
    endif
    
    return, win
end