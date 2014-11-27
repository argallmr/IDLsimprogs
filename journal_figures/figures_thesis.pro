; docformat = 'rst'
;
; NAME:
;    Figures_Thesis
;
; PURPOSE:
;+
;   Create figures for my thesis.
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
;       2014/11/06  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_Asymm3D_y860_Prox
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    xgap      = 6
    xsize     = 700
    ysize     = 550
    col_width = [0.4, 0.6]

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim       = 'Asymm-3D'
    tIndex       = 120816
    xrange       = [5, -4]
    ycell        = 860
    zrange       = 45.3 + [-15, 15]
    coord_system = 'Magnetopause'
    mva_frame    = 1
    ion_scale    = 1
    im_name      = 'Jey'
    oSim         = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange, $
                                ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                                COORD_SYSTEM=coord_system)
    
    ;Set the y-range
    ycoord      = oSim -> Cell2Coord(ycell, /Y)
    oSim.yrange = [ycoord, ycoord]
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts         = [45.2, 43.7, 35.2]
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, CHARSIZE=charsize, COL_WIDTH=col_widht, LAYOUT=layout, $
                        OYMARGIN=oymargin, XGAP=xgap, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        thisLay = gfx.LAYOUT
        colrow  = win -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    cwin = MrSim_ColorSlab(oSim, im_name, HORIZ_LINE=cuts)
    
    ;Move into the other window
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, win
    obj_destroy, cwin
    
    ;Calculate the positions
    win -> SetCurrent
    pos = MrLayout(layout, COL_WIDTH=col_width, CHARSIZE=charsize, OYMARGIN=oymargin, XGAP=xgap)
    im_pos     = [pos[0,8], pos[1,8], pos[2,0], pos[1,1]]
    im_pos[3] -= 0.06

    ;Adjust properties
    win['Cut Uiz']              -> SetProperty, YTICKINTERVAL=0.02
    win['Cut Ex']               -> SetProperty, YTICKINTERVAL=0.02
    win['Color '     + im_name] -> SetProperty, POSITION=im_pos, TITLE=''
    win['CB: Color ' + im_name] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05
    win['HLines '    + im_name] -> SetProperty, COLOR=['Cyan', 'Green', 'Red']
    
    ;Change character size and refresh
    win -> SetGlobal, CHARSIZE=charsize
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_Asymm3D_y1440_Prox
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    xgap      = 6
    xsize     = 700
    ysize     = 550
    col_width = [0.4, 0.6]

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim       = 'Asymm-3D'
    tIndex       = 120816
    xrange       = [8,-8]
    ycell        = 1440
    zrange       = 44+[-15,15]
    coord_system = 'Magnetopause'
    mva_frame    = 1
    ion_scale    = 1
    im_name      = 'Jey'
    oSim         = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange, $
                                ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                                COORD_SYSTEM=coord_system)
    
    ;Set the y-range
    ycoord      = oSim -> Cell2Coord(ycell, /Y)
    oSim.yrange = [ycoord, ycoord]
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts         = [44.7, 39.0, 32.0]
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, CHARSIZE=charsize, COL_WIDTH=col_widht, LAYOUT=layout, $
                        OYMARGIN=oymargin, XGAP=xgap, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        thisLay = gfx.LAYOUT
        colrow  = win -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    cwin = MrSim_ColorSlab(oSim, im_name, HORIZ_LINE=cuts)
    
    ;Move into the other window
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, win
    obj_destroy, cwin
    
    ;Calculate the positions
    win -> SetCurrent
    pos = MrLayout(layout, COL_WIDTH=col_width, CHARSIZE=charsize, OYMARGIN=oymargin, XGAP=xgap)
    im_pos     = [pos[0,8], pos[1,8], pos[2,0], pos[1,1]]
    im_pos[3] -= 0.06

    ;Adjust properties
    win['Cut By']               -> SetProperty, YTICKINTERVAL=0.1
    win['Color '     + im_name] -> SetProperty, POSITION=im_pos, TITLE=''
    win['CB: Color ' + im_name] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05
    win['HLines '    + im_name] -> SetProperty, COLOR=['Cyan', 'Green', 'Red']
    
    ;Change character size and refresh
    win -> SetGlobal, CHARSIZE=charsize
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmScanBy0_FlyBy, $
FNAMES=fnames
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    oymargin  = [4,4]
    xsize     = 700
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim    = 'Asymm-Scan/By0'
    time      = 28
    xrange    = 367.7 + [-50, 50]
    zrange    = [-20, 20]
    ion_scale = 0
    mva_frame = 0
    coord_sys = 'Simulation'
    im_name   = 'Dng_e'
    oSim      = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, $
                             ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                             COORD_SYSTEM=coord_sys)
    
    ;Get the yrange
    ycoord = oSim -> Cell2Coord(0)
    yrange = [ycoord, ycoord]
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts        = [367.7, 349, 330]
    name        = 'ne'
    c_name      = 'Ay'
    im_name     = 'Dng_e'
    nDist       = 25
    dist_layout = [5,5]
    dist_size   = [0.5,0.5,0.5]
    dist_type   = ['Vx-Vy', 'Vx-Vz', 'Vy-Vz', 'Vpar-Vperp', 'Vpar-Vperp1', 'Vpar-Vperp2', 'Vperp1-Vperp2']
    
    nCuts  = n_elements(cuts)
    nTypes = n_elements(dist_type)
    
    imarr    = objarr(nCuts)
    distarr  = objarr(nCuts*nTypes)
    fnames   = strarr(nCuts*nTypes)
    imfnames = strarr(nCuts)
    
    ;Perform the Fly-By
    for i = 0, nCuts - 1 do begin
        for j = 0, nTypes - 1 do begin
            print, FORMAT='(%"Cut %i, Type %s")', cuts[i], dist_type[j]
        
            index = i*nTypes + j
                
            ;Points on satellite path
            r0 = [cuts[i], zrange[0]]
            r1 = [cuts[i], zrange[1]]
        
            ;Create the distributions
            win = MrSim_MMS_FlyBy(oSim, im_name, r0, r1, $
                                  C_NAME      = c_name, $
                                  IM_NAME     = im_name, $
                                  NDIST       = nDist, $
                                  DIST_LAYOUT = dist_layout, $
                                  DIST_SIZE   = dist_size, $
                                  DIST_WIN    = dist_win, $, $
                                  DIST_TYPE   = dist_type[j])
            
            ;Save the windows
            distarr[index] = dist_win
            fnames[index] = string(FORMAT='(%"Asymm-Scan-By0_MMS-FlyBy_%s_x%i.png")', dist_type[j], cuts[i])
            
            ;Images
            if j eq 0 then begin
                imarr[i] = win
                imfnames[i] = string(FORMAT='(%"Asymm-Scan-By0_MMS-FlyBy_%s_x%i.png")', 'Dng', cuts[i])
            endif else begin
                obj_destroy, win
            endelse
            
            ;Make file names
        endfor
    endfor
    
    ;Combine images and distributions
    distarr = [distarr, imarr]
    fnames  = [fnames, imfnames]
    
    ;Return
    return, distarr
end


;+
;   eMap for the Asymm-3D simulation.
;-
function FigThesis_AsymmScanBy0_OhmsLaw, $
FNAME=fname
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        if obj_valid(win2) then obj_destroy, win2
        if obj_valid(win3) then obj_destroy, win3
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Simulation
    theSim       = 'Asymm-Scan/By0'
    tIndex       = 28
    xrange       = [2,-2]
    zrange       = 36.77 + [-5, 5]
    coord_system = 'Magnetopause'
    mva_frame    = 1
    ion_scale    = 1
    horizontal   = 1
    component    = 'X'
    im_name      = 'Jey' 
    oSim   = MrSim_Create(theSim, tIndex, COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame, $
                          ION_SCALE=ion_scale, XRANGE=xrange, ZRANGE=zrange)

    ;Ohm's Law
    cut       = [36.77, 34.9, 33.0]
    win       = MrSim_OhmsLaw(oSim, component, cut[0], HORIZONTAL=horizontal)
    win2      = MrSim_OhmsLaw(oSim, component, cut[1], HORIZONTAL=horizontal)
    win3      = MrSim_OhmsLaw(oSim, component, cut[2], HORIZONTAL=horizontal)
    
    nWins = 2
    nCols = 3
    nRows = 4
    win  -> Refresh, /DISABLE
    win  -> SetProperty, LAYOUT=[nCols, nRows], XGAP=3, XSIZE=900, OXMARGIN=[10,15]
    win2 -> Refresh, /DISABLE
    win3 -> Refresh, /DISABLE

;-------------------------------------------------------
; Switch Windows ///////////////////////////////////////
;-------------------------------------------------------
    for iWin = 1, nWins do begin
        ;Select the window
        case iWin of
            1: theWin = win2
            2: theWin = win3
        endcase
        theWin -> SetProperty, LAYOUT=[nCols, nRows], XGAP=1.5, XSIZE=900
    
        ;Step through each graphic in the window.
        graphics = theWin -> Get(/ALL)
        foreach gfx, graphics do begin
            ;Throw away the legends
            if obj_isa(gfx, 'MrLegend') then continue
        
            ;Append the cut-location to the name
            name = gfx.name
            gfx.name = strtrim(cut[iWin-1], 2) + ' ' + name

            ;Change columns
            ;   - So that they do not push objects in WIN to next column
            isOPlot = gfx -> GetOverplot()
            if ~isOPlot then begin
                layout = gfx.layout
                colrow = win -> ConvertLocation(layout[2], layout[0:1], /PINDEX, /TO_COLROW)
                gfx -> SetLayout, [iWin+1, colrow[1]]
            endif else begin
                layout = !Null
            endelse
            
            ;Switch Windows
            gfx -> SwitchWindows, win
        endforeach

        ;Destroy the window
        obj_destroy, theWin
    endfor

;-------------------------------------------------------
; Format Annotations ///////////////////////////////////
;-------------------------------------------------------
    ;Step through each row
    for row = 1, nRows do begin
        yrange = [!values.f_infinity, -!values.f_infinity]
    
        ;Step through each column
        for col = 1, nCols do begin
            ;Find the graphic
            gfx     = win -> FindByColRow([col,row])
            isOPlot = gfx -> GetOverplot(TARGET=target)
            if isOPlot then gfx = target
        
            ;Get the YRANGE
            yr        = gfx.yrange
            yrange[0] = yr[0] < yrange[0]
            yrange[1] = yr[1] > yrange[1]
        endfor
    
        ;Step through each column
        for col = 1, 3 do begin
            ;Find the graphic
            gfx     = win -> FindByColRow([col,row])
            isOPlot = gfx -> GetOverplot(TARGET=target)
            if isOPlot then gfx = target
        
            ;Set Properties
            gfx.yrange = yrange
            if col gt 1 then gfx -> SetProperty, YTITLE='', YTICKFORMAT='(a1)'
            if col eq 1 then gfx -> SetProperty, YTITLE='E$\downN$'
            if row eq 1 then gfx -> SetProperty, TITLE='$\Omega$$\downci$$\up-1$=64.0 L=' + string(cut[col-1], FORMAT='(f0.1)') + 'd$\downi$'
        endfor

        ;Set Properties
        if row gt 1     then gfx -> SetProperty, TITLE=''
        if row lt nRows then gfx -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
    endfor

;-------------------------------------------------------
; Move Legends /////////////////////////////////////////
;-------------------------------------------------------

    ;Relocate the legends
    win["Ohm's Law"]               -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['34.9000 Total E' + component]
    win["Ohm's Law: VxB term"]     -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['34.9000 E' + component + ' vs. Ec']
    win["Ohm's Law: JxB term"]     -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['34.9000 E' + component + ' vs. Hall E']
    win["Ohm's Law: div(Pe) term"] -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['34.9000 E' + component + ' vs. E inert']    
    
    ;Overview
    im_win = MrSim_ColorSlab(oSim, im_name, C_NAME='Ay', HORIZ_LINES=cut, $
                             LINE_COLOR=['White', 'Forest Green', 'Red'])

    if obj_valid(oSim) then obj_destroy, oSim
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmScanBy0_Prox_3IM
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    col_width = [0.2, 0.2, 0.2, 0.4]
    layout    = [4,5]
    oymargin  = [4,2]
    xsize     = 1000
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim    = 'Asymm-Scan/By0'
    time      = 28
    xrange    = [2,-2]
    zrange    = 36.77 + [-5, 5]
    ion_scale = 1
    mva_frame = 1
    coord_sys = 'Magnetopause'
    im_name   = ['Jey', 'ni', 'Dng_e']
    oSim      = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, $
                             ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                             COORD_SYSTEM=coord_sys)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts = [36.77, 34.9, 33.0]
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, LAYOUT=layout, CHARSIZE=charsize, COL_WIDTH=col_width, $
                        OYMARGIN=oymargin, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        pos          = gfx.POSITION
        gfx.position = [0.67, pos[1], 0.95, pos[3]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    !Null = MrSim_ColorSlab(oSim, im_name[0], C_NAME='Ay', /CURRENT, HORIZ_LINE=cuts)
    !Null = MrSim_ColorSlab(oSim, im_name[1], C_NAME='Ay', /CURRENT, HORIZ_LINE=cuts)
    !Null = MrSim_ColorSlab(oSim, im_name[2], C_NAME='Ay', /CURRENT, HORIZ_LINE=cuts)

    ;Reposition the images
    win['Color ' + im_name[0]] -> SetProperty, POSITION=[0.070, 0.115, 0.220, 0.85], TITLE=''
    win['Color ' + im_name[1]] -> SetProperty, POSITION=[0.247, 0.115, 0.397, 0.85], TITLE='', YTITLE='', YTICKFORMAT='(a1)', RANGE=[0.3,1.37]
    win['Color ' + im_name[2]] -> SetProperty, POSITION=[0.424, 0.115, 0.574, 0.85], TITLE='', YTITLE='', YTICKFORMAT='(a1)'
    
    ;Adjust color properties
    win['CB: Color ' + im_name[0]] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.1
    win['CB: Color ' + im_name[1]] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.4, RANGE=[0.3,1.37]
    win['CB: Color ' + im_name[2]] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05

    win -> SetGlobal, CHARSIZE=charsize
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmScanBy0_Prox
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 2.0
    col_width = [0.4, 0.6]
    layout    = [2,5]
    oymargin  = [4,2]
    xgap      = 10
    xsize     = 700
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim    = 'Asymm-Scan/By0'
    time      = 28
    xrange    = [2,-2]
    zrange    = 36.77 + [-5, 5]
    ion_scale = 1
    mva_frame = 1
    coord_sys = 'Magnetopause'
    im_name   = 'Jey'
    oSim      = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, $
                             ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                             COORD_SYSTEM=coord_sys)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts = [36.77, 34.9, 33.0]
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, LAYOUT=layout, CHARSIZE=charsize, COL_WIDTH=col_width, $
                        OYMARGIN=oymargin, XGAP=xgap, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        thisLay = gfx.LAYOUT
        colrow  = win -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    cwin = MrSim_ColorSlab(oSim, im_name, C_NAME='Ay', HORIZ_LINE=cuts)
    
    ;Move into the other window
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, win
    obj_destroy, cwin
    
    ;Calculate the positions
    pos = MrLayout(layout, CHARSIZE=charsize, COL_WIDTH=col_width, XGAP=xgap, OYMARGIN=oymargin)
    im_pos     = [pos[0,8], pos[1,8], pos[2,0], pos[3,0]]
    im_pos[3] -= 0.08
    
    ;Adjust properties
    win['Cut Bz']               -> SetProperty, YTICKINTERVAL=0.3
    win['Cut By']               -> SetProperty, YTICKINTERVAL=0.05
    win['Cut ni']               -> SetProperty, YTICKINTERVAL=0.5, YRANGE=[0,1.2]
    win['Cut Uiz']              -> SetProperty, YTICKINTERVAL=0.005
    win['Color '     + im_name] -> SetProperty, POSITION=im_pos, TITLE=''
    win['CB: Color ' + im_name] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05

    win -> SetGlobal, CHARSIZE=charsize
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmScanBy1_Prox
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    oymargin  = [4,4]
    xsize     = 700
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim    = 'Asymm-Scan/By1'
    time      = 30
    xrange    = [1,-1]
    zrange    = [42.0, 45.0]
    ion_scale = 1
    mva_frame = 1
    coord_sys = 'Magnetopause'
    im_name   = 'Dng_e'
    oSim      = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, $
                             ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                             COORD_SYSTEM=coord_sys)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts = [43.5, 43.8, 44.5]
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, LAYOUT=layout, CHARSIZE=charsize, OYMARGIN=oymargin, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        thisLay = gfx.LAYOUT
        colrow  = win -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    cwin = MrSim_ColorSlab(oSim, im_name, C_NAME='Ay', HORIZ_LINE=cuts)
    
    ;Move into the other window
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, win
    obj_destroy, cwin
    
    ;Calculate the positions
    pos = MrLayout(layout, CHARSIZE=charsize, OYMARGIN=oymargin)
    im_pos     = [pos[0,8], pos[1,8], pos[2,0], pos[3,0]]
    im_pos[3] -= 0.05
    
    ;Adjust properties
    win['Cut Bz']               -> SetProperty, YTICKINTERVAL=0.25
    win['Cut Uiz']              -> SetProperty, YTICKINTERVAL=0.01
    win['Color '     + im_name] -> SetProperty, POSITION=im_pos, TITLE=''
    win['CB: Color ' + im_name] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05

    win -> SetGlobal, CHARSIZE=charsize
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmScanBy1_FlyBy, $
FNAMES=fnames
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    oymargin  = [4,4]
    xsize     = 700
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim    = 'Asymm-Scan/By1'
    time      = 30
    xrange    = [420, 450]
    zrange    = [-10, 10]
    ion_scale = 0
    mva_frame = 0
    coord_sys = 'Simulation'
    im_name   = 'Dng_e'
    oSim      = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, $
                             ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                             COORD_SYSTEM=coord_sys)
    
    ;Get the yrange
    ycoord = oSim -> Cell2Coord(0)
    yrange = [ycoord, ycoord]
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts        = [435, 438, 445]
    name        = 'ne'
    c_name      = 'Ay'
    im_name     = 'Dng_e'
    nDist       = 25
    dist_layout = [5,5]
    dist_size   = [0.5,0.5,0.5]
    dist_type   = ['Vx-Vy', 'Vx-Vz', 'Vy-Vz', 'Vpar-Vperp', 'Vpar-Vperp1', 'Vpar-Vperp2', 'Vperp1-Vperp2']
    
    nCuts  = n_elements(cuts)
    nTypes = n_elements(dist_type)
    
    imarr    = objarr(nCuts)
    distarr  = objarr(nCuts*nTypes)
    fnames   = strarr(nCuts*nTypes)
    imfnames = strarr(nCuts)
    
    ;Perform the Fly-By
    for i = 0, nCuts - 1 do begin
        for j = 0, nTypes - 1 do begin
            print, FORMAT='(%"Cut %i, Type %s")', cuts[i], dist_type[j]
        
            index = i*nTypes + j
                
            ;Points on satellite path
            r0 = [cuts[i], zrange[0]]
            r1 = [cuts[i], zrange[1]]
        
            ;Create the distributions
            win = MrSim_MMS_FlyBy(oSim, im_name, r0, r1, $
                                  C_NAME      = c_name, $
                                  IM_NAME     = im_name, $
                                  NDIST       = nDist, $
                                  DIST_LAYOUT = dist_layout, $
                                  DIST_SIZE   = dist_size, $
                                  DIST_WIN    = dist_win, $, $
                                  DIST_TYPE   = dist_type[j])
            
            ;Save the windows
            distarr[index] = dist_win
            fnames[index] = string(FORMAT='(%"Asymm-Scan-By1_MMS-FlyBy_%s_x%i.png")', dist_type[j], cuts[i])
            
            ;Images
            if j eq 0 then begin
                imarr[i] = win
                imfnames[i] = string(FORMAT='(%"Asymm-Scan-By1_MMS-FlyBy_%s_x%i.png")', 'Dng', cuts[i])
            endif else begin
                obj_destroy, win
            endelse
            
            ;Make file names
        endfor
    endfor
    
    ;Combine images and distributions
    distarr = [distarr, imarr]
    fnames  = [fnames, imfnames]
    
    ;Return
    return, distarr
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmLarge2D_FlyBy, $
FNAMES=fnames
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    oymargin  = [4,4]
    xsize     = 700
    ysize     = 550

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim    = 'Asymm-Large-2D-NEW'
    time      = 13
    xrange    = 1486 + [-30, 30]
    zrange    = [-20, 20]
    ion_scale = 0
    mva_frame = 0
    coord_sys = 'Simulation'
    im_name   = 'Dng_e'
    oSim      = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, $
                             ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                             COORD_SYSTEM=coord_sys)
    
    ;Get the yrange
    ycoord = oSim -> Cell2Coord(0)
    yrange = [ycoord, ycoord]
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts        = [1486, 1480, 1462]
    c_name      = 'Ay'
    im_name     = 'Dng_e'
    nDist       = 25
    dist_layout = [5,5]
    dist_size   = [0.5,0.5,0.5]
    dist_type   = ['Vx-Vy', 'Vx-Vz', 'Vy-Vz', 'Vpar-Vperp', 'Vpar-Vperp1', 'Vpar-Vperp2', 'Vperp1-Vperp2']
    
    nCuts  = n_elements(cuts)
    nTypes = n_elements(dist_type)
    
    imarr    = objarr(nCuts)
    distarr  = objarr(nCuts*nTypes)
    fnames   = strarr(nCuts*nTypes)
    imfnames = strarr(nCuts)
    
    ;Perform the Fly-By
    for i = 0, nCuts - 1 do begin
        for j = 0, nTypes - 1 do begin
            print, FORMAT='(%"Cut %i, Type %s")', cuts[i], dist_type[j]
        
            index = i*nTypes + j
                
            ;Points on satellite path
            r0 = [cuts[i], zrange[0]]
            r1 = [cuts[i], zrange[1]]
        
            ;Create the distributions
            win = MrSim_MMS_FlyBy(oSim, im_name, r0, r1, $
                                  C_NAME      = c_name, $
                                  IM_NAME     = im_name, $
                                  NDIST       = nDist, $
                                  DIST_LAYOUT = dist_layout, $
                                  DIST_SIZE   = dist_size, $
                                  DIST_WIN    = dist_win, $, $
                                  DIST_TYPE   = dist_type[j])
            
            ;Save the windows
            distarr[index] = dist_win
            fnames[index] = string(FORMAT='(%"Asymm-2D-Large-NEW_MMS-FlyBy_%s_x%i.png")', dist_type[j], cuts[i])
            
            ;Images
            if j eq 0 then begin
                imarr[i] = win
                imfnames[i] = string(FORMAT='(%"Asymm-Large-2D-NEW_MMS-FlyBy_%s_x%i.png")', 'Dng', cuts[i])
            endif else begin
                obj_destroy, win
            endelse
            
            ;Make file names
        endfor
    endfor
    
    ;Combine images and distributions
    distarr = [distarr, imarr]
    fnames  = [fnames, imfnames]
    
    ;Return
    return, distarr
end


;+
;   eMap for the Asymm-3D simulation.
;-
function FigThesis_AsymmLarge2D_OhmsLaw, $
FNAME=fname
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        if obj_valid(win2) then obj_destroy, win2
        if obj_valid(win3) then obj_destroy, win3
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Simulation
    theSim       = 'Asymm-Large-2D'
    tIndex       = 32
    xrange       = [2.0, -2.0]
    zrange       = 150.9 + [-3.0, 3.0]
    coord_system = 'Magnetopause'
    mva_frame    = 1
    ion_scale    = 1
    horizontal   = 1
    component    = 'X'
    im_name      = 'JxB_x' 
    oSim   = MrSim_Create(theSim, tIndex, COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame, $
                          ION_SCALE=ion_scale, XRANGE=xrange, ZRANGE=zrange)

    ;Ohm's Law
    cut       = [150.9, 150.3, 148.4] 
    win       = MrSim_OhmsLaw(oSim, component, cut[0], HORIZONTAL=horizontal)
    win2      = MrSim_OhmsLaw(oSim, component, cut[1], HORIZONTAL=horizontal)
    win3      = MrSim_OhmsLaw(oSim, component, cut[2], HORIZONTAL=horizontal)
    
    nWins = 2
    nCols = 3
    nRows = 4
    win  -> Refresh, /DISABLE
    win  -> SetProperty, LAYOUT=[nCols, nRows], XGAP=3, XSIZE=900, OXMARGIN=[10,15]
    win2 -> Refresh, /DISABLE
    win3 -> Refresh, /DISABLE

;-------------------------------------------------------
; Switch Windows ///////////////////////////////////////
;-------------------------------------------------------
    for iWin = 1, nWins do begin
        ;Select the window
        case iWin of
            1: theWin = win2
            2: theWin = win3
        endcase
        theWin -> SetProperty, LAYOUT=[nCols, nRows], XGAP=1.5, XSIZE=900
    
        ;Step through each graphic in the window.
        graphics = theWin -> Get(/ALL)
        foreach gfx, graphics do begin
            ;Throw away the legends
            if obj_isa(gfx, 'MrLegend') then continue
        
            ;Append the cut-location to the name
            name = gfx.name
            gfx.name = strtrim(cut[iWin-1], 2) + ' ' + name

            ;Change columns
            ;   - So that they do not push objects in WIN to next column
            isOPlot = gfx -> GetOverplot()
            if ~isOPlot then begin
                layout = gfx.layout
                colrow = win -> ConvertLocation(layout[2], layout[0:1], /PINDEX, /TO_COLROW)
                gfx -> SetLayout, [iWin+1, colrow[1]]
            endif else begin
                layout = !Null
            endelse
            
            ;Switch Windows
            gfx -> SwitchWindows, win
        endforeach

        ;Destroy the window
        obj_destroy, theWin
    endfor

;-------------------------------------------------------
; Format Annotations ///////////////////////////////////
;-------------------------------------------------------
    ;Step through each row
    for row = 1, nRows do begin
        yrange = [!values.f_infinity, -!values.f_infinity]
    
        ;Step through each column
        for col = 1, nCols do begin
            ;Find the graphic
            gfx     = win -> FindByColRow([col,row])
            isOPlot = gfx -> GetOverplot(TARGET=target)
            if isOPlot then gfx = target
        
            ;Get the YRANGE
            yr        = gfx.yrange
            yrange[0] = yr[0] < yrange[0]
            yrange[1] = yr[1] > yrange[1]
        endfor
    
        ;Step through each column
        for col = 1, 3 do begin
            ;Find the graphic
            gfx     = win -> FindByColRow([col,row])
            isOPlot = gfx -> GetOverplot(TARGET=target)
            if isOPlot then gfx = target
        
            ;Set Properties
            gfx.yrange = yrange
            if col gt 1 then gfx -> SetProperty, YTITLE='', YTICKFORMAT='(a1)'
            if col eq 1 then gfx -> SetProperty, YTITLE='E$\downN$'
            if row eq 1 then gfx -> SetProperty, TITLE='$\Omega$$\downci$$\up-1$=64.0 L=' + string(cut[col-1], FORMAT='(f0.1)') + 'd$\downi$'
        endfor

        ;Set Properties
        if row gt 1     then gfx -> SetProperty, TITLE=''
        if row lt nRows then gfx -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
    endfor

;-------------------------------------------------------
; Move Legends /////////////////////////////////////////
;-------------------------------------------------------
    ;Relocate the legends
    win["Ohm's Law"]               -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['150.300 Total E' + component]
    win["Ohm's Law: VxB term"]     -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['150.300 E' + component + ' vs. Ec']
    win["Ohm's Law: JxB term"]     -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['150.300 E' + component + ' vs. Hall E']
    win["Ohm's Law: div(Pe) term"] -> SetProperty, LOCATION=8, VSPACE=2.0, TARGET=win['150.300 E' + component + ' vs. E inert']    
    
    ;Overview
    im_win = MrSim_ColorSlab(oSim, im_name, C_NAME='Ay', HORIZ_LINES=cut, $
                             LINE_COLOR=['White', 'Forest Green', 'Red'])

    if obj_valid(oSim) then obj_destroy, oSim
    win -> Refresh
    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmLarge2D_t32_Prox
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    xgap      = 6
    xsize     = 700
    ysize     = 550
    col_width = [0.4, 0.6]

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim       = 'Asymm-Large-2D'
    tIndex       = 32
    xrange       = [2.0, -2.0]
    zrange       = 150.9 + [-3.0, 3.0]
    coord_system = 'Magnetopause'
    mva_frame    = 1
    ion_scale    = 1
    im_name      = 'Jey'
    oSim         = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange, $
                                ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                                COORD_SYSTEM=coord_system)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts       = [150.9, 150.3, 148.4] 
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, CHARSIZE=charsize, COL_WIDTH=col_widht, LAYOUT=layout, $
                        OYMARGIN=oymargin, XGAP=xgap, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        thisLay = gfx.LAYOUT
        colrow  = win -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    cwin = MrSim_ColorSlab(oSim, im_name, C_NAME='Ay', HORIZ_LINE=cuts)
    
    ;Move into the other window
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, win
    obj_destroy, cwin
    
    ;Calculate the positions
    win -> SetCurrent
    pos = MrLayout(layout, COL_WIDTH=col_width, CHARSIZE=charsize, OYMARGIN=oymargin, XGAP=xgap)
    im_pos     = [pos[0,8], pos[1,8], pos[2,0], pos[1,0]]
    im_pos[3] -= 0.06
    
    ;Adjust properties
    win['Cut Uiz']              -> SetProperty, YTITLE='U$\downiZ$', YTICKINTERVAL=0.01
    win['Cut Ex']               -> SetProperty, YTITLE='E$\downiX$'
    win['Color '     + im_name] -> SetProperty, POSITION=im_pos, TITLE=''
    win['CB: Color ' + im_name] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05
    win['HLines '    + im_name] -> SetProperty, COLOR=['Cyan', 'Forest Green', 'Red']

    ;Change character size and refresh
    win -> SetGlobal, CHARSIZE=charsize
    win -> Refresh

    return, win
end


;+
;   Create Figure 2: 2D simulation.
;-
function FigThesis_AsymmLarge2D_t90_Prox
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(cwin) then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Layout
    charsize  = 1.5
    layout    = [2,5]
    xgap      = 6
    xsize     = 700
    ysize     = 550
    col_width = [0.4, 0.6]

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    theSim       = 'Asymm-Large-2D'
    tIndex       = 90
    xrange       = [3.0, -3.0]
    zrange       = 153.3 + [-7.0, 7.0]
    coord_system = 'Magnetopause'
    mva_frame    = 1
    ion_scale    = 1
    im_name      = 'Jey'
    oSim         = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange, $
                                ION_SCALE=ion_scale, MVA_FRAME=mva_frame, $
                                COORD_SYSTEM=coord_system)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    cuts         = [153.3, 150.8, 148.0]
    
    ;Create cuts of Bx, By, ni, Uix, Ez
    win = MrSim_XProximity(oSim, cuts)
    win -> Refresh, /DISABLE
    
    ;Create a second column
    win -> SetProperty, CHARSIZE=charsize, COL_WIDTH=col_widht, LAYOUT=layout, $
                        OYMARGIN=oymargin, XGAP=xgap, XSIZE=xsize, YSIZE=ysize

    ;Move all graphics into it
    graphics = win -> Get(/ALL, ISA='MrPlot')
    foreach gfx, graphics do begin
        thisLay = gfx.LAYOUT
        colrow  = win -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
        gfx -> SetLayout, [2, colrow[1]]
    endforeach

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create a 2D color plot
    cwin = MrSim_ColorSlab(oSim, im_name, C_NAME='Ay', HORIZ_LINE=cuts)
    
    ;Move into the other window
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, win
    obj_destroy, cwin
    
    ;Calculate the positions
    win -> SetCurrent
    pos = MrLayout(layout, COL_WIDTH=col_width, CHARSIZE=charsize, OYMARGIN=oymargin, XGAP=xgap)
    im_pos     = [pos[0,8], pos[1,8], pos[2,0], pos[1,1]]
    im_pos[3] -= 0.06

    ;Adjust properties
    win['Cut Uiz']              -> SetProperty, YTITLE='U$\downiZ$', YTICKINTERVAL=0.02
    win['Cut Ex']               -> SetProperty, YTITLE='E$\downiX$', YTICKINTERVAL=0.02
    win['Color '     + im_name] -> SetProperty, POSITION=im_pos, TITLE=''
    win['CB: Color ' + im_name] -> SetProperty, CBLOCATION='Top', OFFSET=0.5, WIDTH=1.5, TICKINTERVAL=0.05
    win['HLines '    + im_name] -> SetProperty, COLOR=['Cyan', 'Forest Green', 'Red']
    
    ;Change character size and refresh
    win -> SetGlobal, CHARSIZE=charsize
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
function Figures_Thesis, figure, $
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
    list_of_figures = [['Asymm-3D',       ''], $
                       ['    y860  Prox', ''], $
                       ['    y1440 Prox', ''], $
                       ['Asymm-Scan-By1', ''], $
                       ['    FlyBy',      ''], $
                       ['    Prox',       ''], $
                       ['Asymm-Scan-By0', ''], $
                       ['    FlyBy',      ''], $
                       ['    Prox',       ''], $
                       ['    Ohms Law',   ''], $
                       ['Asymm-Large-2D', ''], $
                       ['    FlyBy',      ''], $
                       ['    t32 Prox',   ''], $
                       ['    t90 Prox',   ''], $
                       ['    Ohms Law',   '']]
    
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
        'ASYMM-3D Y860 PROX':      win = FigThesis_Asymm3D_Y860_Prox()
        'ASYMM-3D Y1440 PROX':     win = FigThesis_Asymm3D_Y1440_Prox()
        'ASYMM-SCAN-BY0 FLYBY':    win = FigThesis_AsymmScanBy0_FlyBy(FNAMES=fnames)
        'ASYMM-SCAN-BY0 OHMS LAW': win = FigThesis_AsymmScanBy0_OhmsLaw()
        'ASYMM-SCAN-BY0 PROX':     win = FigThesis_AsymmScanBy0_Prox()
        'ASYMM-SCAN-BY0 PROX 3IM': win = FigThesis_AsymmScanBy0_Prox_3IM()
        'ASYMM-SCAN-BY1 FLYBY':    win = FigThesis_AsymmScanBy1_FlyBy(FNAMES=fnames)
        'ASYMM-SCAN-BY1 PROX':     win = FigThesis_AsymmScanBy1_Prox()
        'ASYMM-LARGE-2D FLYBY':    win = FigThesis_AsymmLarge2D_FlyBy(FNAMES=fnames)
        'ASYMM-LARGE-2D OHMS LAW': win = FigThesis_AsymmLarge2D_OhmsLaw()
        'ASYMM-LARGE-2D T32 PROX': win = FigThesis_AsymmLarge2D_t32_Prox()
        'ASYMM-LARGE-2D T90 PROX': win = FigThesis_AsymmLarge2D_t90_Prox()
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
            fname = 'MrThesis_' + idl_validname(figure, /CONVERT_ALL)
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
            for i = 0, nWins - 1 do win[i] -> Save, FilePath(fnames[i], ROOT_DIR=froot)
        endelse
    endif
    
    return, win
end