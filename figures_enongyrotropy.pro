;+
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function FEN_eMap, simname, tIndex, yslice, $
LEFT=left, $
EMAP=eMap, $
STRNAME=strname, $
SIM_OBJECT=oSim
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(im_win)    then obj_destroy, im_win
        if max(obj_valid(wins)) then obj_destroy, wins
        void = cgErrorMSG()
        return, !Null
    endif
    
    if keyword_set(left) then xline = 'LEFT' else xline = 'RIGHT'

    ;Electron bins
    bin_loc    = 5
    bin_hsize  = 0.5
    bin_vsize  = 0.5
    bin_layout = [5,5]
    bin_hspace = 0.5
    bin_vspace = 0.5
    nBins      = 75
    
;-------------------------------------------------------
; Data Range ///////////////////////////////////////////
;-------------------------------------------------------
    case simname of
        'ASYMM-SCAN/BY1': begin
            case tIndex of
                30: begin
                    case xline of
                        'LEFT': message, 'No particle data for left x-line yet'
                        'RIGHT': begin
                            v_va       = 1
                            bin_center = [434.6, 0.4]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            xrange     = [420, 450]
                            zrange     = [-10, 10]
                        endcase
                    endcase
                endcase
            endcase
        endcase
        'ASYMM-SCAN/BY0': begin
            case tIndex of
                28: begin
                    v_va       = 1
                    bin_center = [367.6, -0.9]
                    bin_hspace = 0.0
                    bin_vspace = 0.0
                    bin_loc    = 5
                    xrange     = [345, 390]
                    zrange     = [-10, 8]
                endcase
            endcase
        endcase
        'ASYMM-LARGE-2D-NEW': begin
            case tIndex of
                26: begin
                    case xline of
                        'LEFT': begin
                            v_va       = 1
                            bin_center = [1490, 0]
                            xrange     = bin_center[0] + [-20, 20]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            zrange     = [-10, 10]
                        endcase
                        'RIGHT': begin
                            v_va       = 1
                            bin_center = [1782.1, 2.4]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            xrange     = [1770, 1795]
                            zrange     = [-10,10]
                        endcase
                    endcase
                endcase
                36: begin
                    case xline of
                        'LEFT': begin
                            v_va       = 1
                            bin_center = [1517, 2.7]
                            xrange     = bin_center[0] + [-20, 20]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            zrange     = [-10, 10]
                        endcase
                        'RIGHT': begin
                            v_va       = 1
                            bin_center = [2021, 4.6]
                            xrange     = bin_center[0] + [-20, 20]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            zrange     = [-10, 15]
                        endcase
                    endcase
                endcase
            endcase
        endcase
        'ASYMM-3D': begin
            case tIndex of
                108090: begin
                    case yslice of
                        650: begin
                            v_va       = 1
                            bin_center = [471.5, 3.2]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            xrange     = bin_center[0] + [-20,20]
                            zrange     = [-10, 15]
                        endcase
                        905: begin
                            v_va       = 1
                            bin_center = [459, 3.1]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            xrange     = bin_center[0] + [-20,20]
                            zrange     = [-10, 15]
                        endcase
                        1440: begin
                            v_va       = 1
                            bin_center = [420, 7]
                            bin_hspace = 0.0
                            bin_vspace = 0.0
                            bin_loc    = 5
                            xrange     = [400,490]
                            zrange     = [-10, 15]
                        endcase
                    endcase
                endcase
            endcase
        endcase
        'SIM1': begin
            case tIndex of
                68: begin
                    v_va       = 1
                    binary     = 1
                    bin_center = [839, 0.0]
                    bin_loc    = 5
                    bin_layout = [5,5]
                    bin_hspace = 0.0
                    bin_vspace = 0.0
                    xrange     = bin_center[0] + [-25,25]
                    zrange     = [-7, 7]
                endcase
                72: begin
                    v_va       = 1
                    binary     = 1
                    bin_center = [842, 0.0]
                    bin_loc    = 5
                    bin_layout = [5,5]
                    bin_hspace = 0.0
                    bin_vspace = 0.0
                    xrange     = [810, 880]
                    zrange     = [-10, 10]
                endcase
                80: begin
                    v_va       = 1
                    binary     = 1
                    bin_center = [858, 0.0]
                    bin_loc    = 5
                    bin_layout = [5,5]
                    bin_hspace = 0.0
                    bin_vspace = 0.0
                    xrange     = bin_center[0] + [-35,35]
                    zrange     = [-7, 7]
                endcase
                116: begin
                    v_va       = 1
                    binary     = 1
                    bin_center = [900, 0.4]
                    bin_loc    = 5
                    bin_layout = [5,5]
                    bin_hspace = 0.0
                    bin_vspace = 0.0
                    xrange     = bin_center[0] + [-25,25]
                    zrange     = [-10,10]
                endcase
            endcase
        endcase
        else: message, 'Simulation name "' + simname + '" not recognized.'
    endcase
    
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------

    ;Create the simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = MrSim_Create(simname, tIndex, yslice, XRANGE=xrange, ZRANGE=zrange, BINARY=binary)
                                          
    ;Make sure the parameters are correct
    endif else begin
        ;Get the name of the current simulation
        oSim -> GetProperty, SIMNUM=number
        MrSim_Which, number, NAME=name
        
        ;Switch simulations?
        if strupcase(simname) ne strupcase(name) then begin
            obj_destroy, oSim
            oSim = MrSim_Create(simname, tIndex, yslice, XRANGE=xrange, ZRANGE=zrange, BINARY=binary)
        endif else begin
            oSim -> GetProperty, SIMNUM=simnum, TIME=tt, YSLICE=yy, XRANGE=xx, ZRANGE=zz
            if (tIndex ne tt) || (array_equal(xx, xrange) eq 0) || (array_equal(zz, zrange) eq 0) $
                then oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
            if obj_class(oSim) eq 'MRSIM3D' then if yslice ne yy $
                then oSim -> SetProperty, YSLICE=yslice
        endelse
    endelse

;-------------------------------------------------------
;Create eMaps //////////////////////////////////////////
;-------------------------------------------------------
    ;Determine which distribution type was selected
    ;   - 2^0 bit stored as last elements, must reverse
    type = ['x-z', 'x-Vx', 'z-Vz', 'Vx-Vy', 'Vx-Vz', 'Vy-Vz', $
            'Vpar-Vperp', 'Vpar-Vperp1', 'Vpar-Vperp2', 'Vperp1-Vperp2']
    bitSet = cgBitGet(eMap)
    iBits   = where(reverse(bitSet), nBits)
    if nBits eq 0 then message, 'No distributions selected. See keyword list.', /INFORMATIONAL

    ;Allocate space for the plotting windows
    ;   - Include an extra element for the overview window.
    wins    = objarr(nBits+1)
    name    = 'Jey'
    
    ;Step through each distribution type.
    for i = 0, nBits - 1 do begin
        index  = iBits[i]
        if i eq 0 then im_name = name else im_name = ''
        if i eq 0 && n_elements(yslice) eq 0 then c_name  = 'Ay' else c_name  = ''
        
        ;Create the desired distribution functions
        wins[i] = MrSim_eMap(oSim, type[index], bin_center[0], bin_center[1], bin_hsize, bin_vsize, $
                             LAYOUT     = bin_layout, $
                             LOCATION   = bin_loc, $
                             HGAP       = bin_hspace, $
                             VGAP       = bin_vspace, $
                             V_VA       = v_va, $
                             IM_NAME    = im_name, $
                             C_NAME     = c_name, $
                             IM_WIN     = im_win, $
                             NBINS      = nBins)
        
        ;Store the overview plot as well
        if i eq 0 && obj_valid(im_win) then wins[-1] = im_win
    endfor

    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
    
;---------------------------------------------------------------------
; Figure Names ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;time-slice
    strTime = string(FORMAT='(%"t%i_")', tIndex)
    
    ;y-slice, left/right x-line
    if n_elements(yslice) gt 0 then begin
        strYSlice = string(FORMAT='(%"y%04i_")', yslice)
        strLR     = ''
    endif else begin
        strYSlice = ''
        if n_elements(left) gt 0 $
            then strLR = left ? 'left_' : 'right_' $
            else strLR = ''
    endelse
    
    ;Combine components
    strhead = 'eMap_'
    strbutt = string(FORMAT='(%"_%ix%i_%04ic%02i_loc%i_%ix%ide_%03ibins")', $
                     bin_layout, bin_center, bin_loc, 2*bin_hsize, 2*bin_vsize, nBins[0])
    strname = strTime + strYSlice + strLR + strhead + [type[iBits], name] + strbutt
    
    ;Return
    return, wins
end


;+
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function FEN_Overview_t06, $
SIM_OBJECT=oSim, $
STRNAME=strname
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif

    ;Data Range
    time       = 6
    xrange     = [1300, 1700]
    zrange     = [ -10,    5]
    dir        = '/home/argall/simulations/Asymm-Large-2D-NEW/'
    
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
    ;Create the simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = obj_new('MrSim2D', time, XRANGE=xrange, ZRANGE=zrange, $
                                        DIRECTORY=dir+'data', INFO_FILE=dir+'info')
                                          
    ;Make sure the parameters are correct
    endif else begin
        oSim -> GetProperty, TIME=tt, XRANGE=xx, ZRANGE=zz
        if (time ne tt) || (array_equal(xx, xrange) eq 0) || (array_equal(zz, zrange) eq 0) $
            then oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    endelse

;-------------------------------------------------------
; Overview Plot ////////////////////////////////////////
;-------------------------------------------------------
    win = MrSim_MultiSlab(['Dng_e', 'ne', 'Jey', 'A0_e', 'An_e'], C_NAME='Ay', SIM_OBJECT=oSim)

    ;Make pretty
    win -> Refresh, /DISABLE
    win -> SetProperty, XSIZE=550, YSIZE=600

    ;Create a filename
    strname = string(FORMAT='(%"t%02i_overview_%04ix%04i_%02iz%02i")', $
                    time, xrange, abs(zrange))


    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function FEN_Overview_t09, $
SIM_OBJECT=oSim, $
STRNAME=strname
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif

    ;Data Range
    time       = 9
    xrange     = [1300, 1700]
    zrange     = [  -6,    8]
    dir        = '/home/argall/simulations/Asymm-Large-2D-NEW/'
    
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
    ;Create the simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = obj_new('MrSim2D', time, XRANGE=xrange, ZRANGE=zrange, $
                                        DIRECTORY=dir+'data', INFO_FILE=dir+'info')
                                          
    ;Make sure the parameters are correct
    endif else begin
        oSim -> GetProperty, TIME=tt, XRANGE=xx, ZRANGE=zz
        if (time ne tt) || (array_equal(xx, xrange) eq 0) || (array_equal(zz, zrange) eq 0) $
            then oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    endelse

;-------------------------------------------------------
; Overview Plot ////////////////////////////////////////
;-------------------------------------------------------
    win = MrSim_MultiSlab(['Dng_e', 'ne', 'Jey', 'A0_e', 'An_e'], C_NAME='Ay', SIM_OBJECT=oSim)

    ;Make pretty
    win -> Refresh, /DISABLE
    win -> SetProperty, XSIZE=550, YSIZE=600

    ;Create a filename
    strname = string(FORMAT='(%"t%02i_overview_%04ix%04i_%02iz%02i")', $
                    time, xrange, abs(zrange))

    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function FEN_Overview_t13, $
SIM_OBJECT=oSim, $
STRNAME=strname
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif

    ;Data Range
    time       = 13
    xrange     = [1300, 1700]
    zrange     = [ -20,   15]
    dir        = '/home/argall/simulations/Asymm-Large-2D-NEW/'
    
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
    ;Create the simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = obj_new('MrSim2D', time, XRANGE=xrange, ZRANGE=zrange, $
                                        DIRECTORY=dir+'data', INFO_FILE=dir+'info')
                                          
    ;Make sure the parameters are correct
    endif else begin
        oSim -> GetProperty, TIME=tt, XRANGE=xx, ZRANGE=zz
        if (time ne tt) || (array_equal(xx, xrange) eq 0) || (array_equal(zz, zrange) eq 0) $
            then oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    endelse

;-------------------------------------------------------
; Overview Plot ////////////////////////////////////////
;-------------------------------------------------------
    win = MrSim_MultiSlab(['Dng_e', 'ne', 'Jey', 'A0_e', 'An_e'], C_NAME='Ay', SIM_OBJECT=oSim)

    ;Make pretty
    win -> Refresh, /DISABLE
    win -> SetProperty, XSIZE=550, YSIZE=600

    ;Create a filename
    strname = string(FORMAT='(%"t%02i_overview_%04ix%04i_%02iz%02i")', $
                    time, xrange, abs(zrange))

    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function FEN_Overview_t18, $
SIM_OBJECT=oSim, $
STRNAME=strname
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif

    ;Data Range
    time       = 18
    xrange     = [1300, 1850]
    zrange     = [ -60,   45]
    dir        = '/home/argall/simulations/Asymm-Large-2D-NEW/'
    
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
    ;Create the simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = obj_new('MrSim2D', time, XRANGE=xrange, ZRANGE=zrange, $
                                        DIRECTORY=dir+'data', INFO_FILE=dir+'info')
                                          
    ;Make sure the parameters are correct
    endif else begin
        oSim -> GetProperty, TIME=tt, XRANGE=xx, ZRANGE=zz
        if (time ne tt) || (array_equal(xx, xrange) eq 0) || (array_equal(zz, zrange) eq 0) $
            then oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    endelse

;-------------------------------------------------------
; Overview Plot ////////////////////////////////////////
;-------------------------------------------------------
    win = MrSim_MultiSlab(['Dng_e', 'ne', 'Jey', 'A0_e', 'An_e'], C_NAME='Ay', SIM_OBJECT=oSim)

    ;Make pretty
    win -> Refresh, /DISABLE
    win -> SetProperty, XSIZE=550, YSIZE=600

    ;Create a filename
    strname = string(FORMAT='(%"t%02i_overview_%04ix%04i_%02iz%02i")', $
                    time, xrange, abs(zrange))

    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function FEN_Figure1_LiJen, $
SIM_OBJECT=oSim
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)     then obj_destroy, win
        if obj_valid(uexwin)  then obj_destroy, uexWin
        if obj_valid(winVxVy) then obj_destroy, winVxVy
        if obj_valid(winVxVz) then obj_destroy, winVxVz
        if obj_valid(winVyVz) then obj_destroy, winVyVz
        void = cgErrorMSG()
        return, !Null
    endif
    
;-------------------------------------------------------
; Create a Window & Simulation Object //////////////////
;-------------------------------------------------------
    oxmargin = [10, 15]
    xsize    = 800
    ysize    = 590
    win = MrWindow(OXMARGIN=oxmargin, XSIZE=xsize, YSIZE=ysize, NAME='Fig1: LiJen', REFRESH=0)
    
    oSim = MrSim_Create('Sim1', 76, $
                        /BINARY, $
                        XRANGE=[800, 900], $
                        ZRANGE=[-20,20])
                        
;-------------------------------------------------------
; Line Cut /////////////////////////////////////////////
;-------------------------------------------------------
    !Null = MrSim_LineCut(oSim, 'Uex', 0, /HORIZONTAL, /CURRENT)
                        
;-------------------------------------------------------
; Generate an eMap /////////////////////////////////////
;-------------------------------------------------------
    nBins      = 75
    bin_center = [847.8, 0]
    bin_layout = [5, 1]
    bin_hsize  = 0.5
    bin_vsize  = 0.5
    bin_hspace = 5
    location   = 4
    
    ;Vx-Vy
    winVxVy = MrSim_eMap(oSim, 'Vx-Vy', bin_center[0], bin_center[1], bin_hsize, bin_vsize, $
                         C_NAME   = 'Ay', $
                         HGAP     = bin_hspace, $
                         IM_NAME  = 'Uex', $
                         IM_WIN   = uexWin, $
                         LAYOUT   = bin_layout, $
                         LOCATION = location, $
                         NBINS    = nBins)
    
    ;Vx-Vz
    winVxVz = MrSim_eMap(oSim, 'Vy-Vz', bin_center[0], bin_center[1], bin_hsize, bin_vsize, $
                         HGAP     = bin_hspace, $
                         LAYOUT   = bin_layout, $
                         LOCATION = location, $
                         NBINS    = nBins)
    
    ;Vy-Vz
    winVyVz = MrSim_eMap(oSim, 'Vy-Vz', bin_center[0], bin_center[1], bin_hsize, bin_vsize, $
                         HGAP     = bin_hspace, $
                         LAYOUT   = bin_layout, $
                         LOCATION = location, $
                         NBINS    = nBins)
    
    ;Turn refresh off
    winVxVy -> Refresh, /DISABLE
    winVxVz -> Refresh, /DISABLE
    winVyVz -> Refresh, /DISABLE
                        
;-------------------------------------------------------
; Move Color and Line Plot into Main Window ////////////
;-------------------------------------------------------
    colorPos = [0.25, 0.77, 0.75, 0.95]
    linePos  = [0.25, 0.70, 0.75, 0.77]

    ;Move Color Uix
    tempGfx = uexWin -> Get(/ALL)
    foreach gfx, tempGfx do gfx -> SwitchWindows, win
    
    ;Destroy the empty window
    obj_destroy, uexWin

    ;Reposition the color and line plot
    win['Color Uex'] -> SetLayout, POSITION=colorPos
    win['Cut Uex']   -> SetLayout, POSITION=linePos
                        
;-------------------------------------------------------
; Move Distributions into Main Window //////////////////
;-------------------------------------------------------
    ;Width and height of the distributions
    ;   - Make them square
    height = 0.17
    width  = height * xsize / ysize
    
    ;Positions of the distributions
    x0 = 0.1 + width*indgen(5)
    x1 = x0  + width
    y1 = 0.6 - height*indgen(3)
    y0 = y1  - height

    ;Move VxVy
    iDist = 1
    tempGfx = winVxVy -> Get(/ALL)
    foreach gfx, tempGfx do begin
        gfx -> SwitchWindows, win
        if obj_class(gfx) eq 'MRIMAGE' then begin
            gfx -> SetLayout, POSITION=[x0[i], y0[0], x1[i], y1[0]]
            i += 1
        endif
    endforeach
    obj_destroy, winVxVy

    ;Move VxVz
    iDist = 1
    tempGfx = winVxVz -> Get(/ALL)
    foreach gfx, tempGfx do begin
        gfx -> SwitchWindows, win
        if obj_class(gfx) eq 'MRIMAGE' then begin
            gfx -> SetLayout, POSITION=[x0[i], y0[1], x1[i], y1[1]]
            i += 1
        endif
    endforeach
    obj_destroy, winVxVz

    ;Move VxVz
    iDist = 1
    tempGfx = winVyVz -> Get(/ALL)
    foreach gfx, tempGfx do begin
        gfx -> SwitchWindows, win
        if obj_class(gfx) eq 'MRIMAGE' then begin
            gfx -> SetLayout, POSITION=[x0[i], y0[2], x1[i], y1[2]]
            i += 1
        endif
    endforeach
    obj_destroy, winVyVz

;-------------------------------------------------------
; Create the Image /////////////////////////////////////
;-------------------------------------------------------

    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function FEN_Figure1_Dng, $
SIM_OBJECT=oSim
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif
    
;-------------------------------------------------------
; Create a Window //////////////////////////////////////
;-------------------------------------------------------
    win = MrWindow(OXMARGIN=[10,15], YGAP=4, YSIZE=590, NAME='Fig1: Dng Overviews', REFRESH=0)
    
;-------------------------------------------------------
; Sim1 /////////////////////////////////////////////////
;-------------------------------------------------------
    ;Data Range
    simname = 'Sim1'
    time    = 72
    xrange  = 842 + [-30, 30]
    zrange  = [-10, 10]
    
    ;Simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = MrSim_Create(simname, time, /BINARY, XRANGE=xrange, ZRANGE=zrange)
    endif else begin
        oSim -> SetSim, simname, /BINARY
        oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    endelse
    
    ;Plot
    win = MrSim_ColorSlab(oSim, 'Dng_e', C_NAME='Ay', /CURRENT)

    ;SetProperties
    title = win['Color Dng_e'].TITLE
    win['Color Dng_e']     -> SetProperty, NAME='Sim1 Dng', XTITLE='', TITLE='Symmetric B$\downy$=0 ' + title
    win['CB: Color Dng_e'] -> SetProperty, NAME='Sim1 CB'
    win['Ay Contours']     -> SetProperty, NAME='Sim1 Ay Contours'
    
;-------------------------------------------------------
; Asymm-Scan/By0 ///////////////////////////////////////
;-------------------------------------------------------

    ;Data Range
    simname = 'Asymm-Scan/By0'
    time    = 28
    xrange  = 367 + [-30, 30]
    zrange  = [-10, 10]
    
    ;Change simulations
    oSim -> SetSim, simname
    oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    
    ;Plot
    win = MrSim_ColorSlab(oSim, 'Dng_e', C_NAME='Ay', /CURRENT)
    
    ;SetProperties
    title = win['Color Dng_e'].TITLE
    win['Color Dng_e']     -> SetProperty, NAME='By0 Dng', XTITLE='', TITLE='Asymmetric B$\downy$=0 ' + title
    win['CB: Color Dng_e'] -> SetProperty, NAME='By0 CB'
    win['Ay Contours']     -> SetProperty, NAME='By0 Ay Contours'

    
;-------------------------------------------------------
; Asymm-Scan/By1 ///////////////////////////////////////
;-------------------------------------------------------

    ;Data Range
    simname = 'Asymm-Scan/By1'
    time    = 30
    xrange  = 434 + [-30,30]
    zrange  = [-10, 10]
    
    ;Change simulations
    oSim -> SetSim, simname
    oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    
    ;Plot
    win = MrSim_ColorSlab(oSim, 'Dng_e', C_NAME='Ay', /CURRENT)
    
    ;SetProperties
    title = win['Color Dng_e'].TITLE
    win['Color Dng_e']     -> SetProperty, NAME='By1 Dng', TITLE='Asymmetric B$\downy$=1 ' + title
    win['CB: Color Dng_e'] -> SetProperty, NAME='By1 CB'
    win['Ay Contours']     -> SetProperty, NAME='By1 Ay Contours'

;-------------------------------------------------------
; Create the Image /////////////////////////////////////
;-------------------------------------------------------

    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function FEN_Sym_Overview_t72, $
SIM_OBJECT=oSim
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif

    ;Data Range
    time       = 72
    xrange     = 844 + [-40, 40]
    zrange     = [-20, 20]
    dir        = '/data1/sim1/'
    
;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
    ;Create the simulation object
    if obj_valid(oSim) eq 0 then begin
        oSim = obj_new('MrSim2D', time, XRANGE=xrange, ZRANGE=zrange, $
                                        DIRECTORY=dir+'data', /BINARY)
                                          
    ;Make sure the parameters are correct
    endif else begin
        oSim -> GetProperty, TIME=tt, XRANGE=xx, ZRANGE=zz
        if (time ne tt) || (array_equal(xx, xrange) eq 0) || (array_equal(zz, zrange) eq 0) $
            then oSim -> SetProperty, TIME=time, XRANGE=xrange, ZRANGE=zrange
    endelse

;-------------------------------------------------------
; Overview Plot ////////////////////////////////////////
;-------------------------------------------------------
    win = MrSim_MultiSlab(['Dng_e', 'ne', 'Jey', 'A0_e', 'An_e'], C_NAME='Ay', SIM_OBJECT=oSim)

    ;Make pretty
    win -> Refresh, /DISABLE
    win -> SetProperty, XSIZE=550, YSIZE=600


    ;Destroy the simulation object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
function Figures_eNongyrotropy, figure, $
SAVE=tf_save, $
SHOW_FIGS=show_figs, $
SIM_OBJECT=oSim, $
X_Z=x_z, $
X_VX=x_vx, $
Z_VZ=z_vz, $
VX_VZ=vx_vz, $
VX_VY=vx_vy, $
VY_VZ=vy_vz, $
VPERP_VPAR=vpar_vperp, $
VPAR_VPERP1=vpar_vperp1, $
VPAR_VPERP2=vpar_vperp2, $
VPERP1_VPERP2=vperp1_vperp2

    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if max(obj_valid(win)) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif
    
;---------------------------------------------------------------------
; Info ///////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    
    ;Current list of figures
    list_of_figures = [['Figures eNongyrotropy', 'Figures used in my study of nongyrotropy in asymmetric guide field case.'], $
                       ['Figure1 Dng          ', ''], $
                       ['Figure1 LiJen        ', ''], $
                       ['Asymm-Scan/By0       ', 'Asymm-Scan/By0 simulation figures'], $
                       ['    t28 eMap         ', ''], $
                       ['Asymm-Scan/By1       ', 'Asymm-Scan/By1 simulation figures'], $
                       ['    t30 eMap         ', ''], $, $
                       ['3D                   ', 'Asymm-3D simulation figures'], $
                       ['    t108090 y650 eMap  ', ''], $
                       ['    t108090 y905 eMap  ', ''], $
                       ['    t108090 y1440 eMap ', ''], $
                       ['Large                ', 'Asymm-Large-2D'], $
                       ['Asymm-Large-2D-NEW   ', 'Asymm-Large-2D-NEW'], $
                       ['    t06 Overview     ', ''], $
                       ['    t09 Overview     ', ''], $
                       ['    t13 Overview     ', ''], $
                       ['    t18 Overview     ', ''], $
                       ['    t26 Left eMap    ', ''], $
                       ['    t26 Right eMap   ', ''], $
                       ['    t36 Left eMap    ', ''], $
                       ['    t36 Right eMap   ', ''], $
                       ['Sim1                 ', 'Sim1 symmetric simulation figures'], $
                       ['    t68 eMap         ', ''], $
                       ['    t72 eMap         ', ''], $
                       ['    t80 eMap         ', ''], $
                       ['   t116 eMap         ', '']]
    
    ;Print the list of figures?
    if keyword_set(show_figs) then begin
        len = strtrim(max(strlen(list_of_figures[0,*])), 2)
        print, 'FIGURES', 'DESCRIPTION', FORMAT='(a' + len + ', 4x, a11)'
        print, list_of_figures, FORMAT='(a-' + len + ', 4x, a0)'
        return, obj_new()
    endif
    
;---------------------------------------------------------------------
; Create Figure //////////////////////////////////////////////////////
;---------------------------------------------------------------------
    eMap    = 0L
    eMap    =     keyword_set(x_z) + $
              2 * keyword_set(x_vx) + $
              4 * keyword_set(z_vz) + $
              8 * keyword_set(vx_vy) + $
             16 * keyword_set(vx_vz) + $
             32 * keyword_set(vy_vz) + $
             64 * keyword_set(vpar_vperp) + $
            128 * keyword_set(vpar_vperp1) + $
            256 * keyword_set(vpar_vperp2) + $
            512 * keyword_set(vperp1_vperp2) 

    ;Defaults
    tf_save = keyword_set(tf_save)
    _figure = strupcase(figure)
    
    
;-------------------------------------------------------
; eMap? ////////////////////////////////////////////////
;-------------------------------------------------------
    iseMap = stregex(_figure, 'EMAP', /BOOLEAN)
    if iseMap then begin
        parts = stregex(_figure, '([A-Z0-9/-]+)( T[0-9]+)( Y[0-9]+)? (LEFT|RIGHT)?[ ]?EMAP', /SUBEXP, /EXTRACT)
        if parts[0] eq '' then message, 'Regex mistake.'
        
        simname = parts[1]
        tIndex  = long(strmid(parts[2], 2))
        if parts[3] ne '' then yslice = long(strmid(parts[3], 2))
        if parts[4] ne '' then left   = parts[4] eq 'LEFT'
        
        win = FEN_eMap(simname, tIndex, yslice, LEFT=left, SIM_OBJECT=oSim, STRNAME=strname, EMAP=eMap)

;-------------------------------------------------------
; Other Figure /////////////////////////////////////////
;-------------------------------------------------------
    endif else begin

        ;Create the figure    
        case _figure of
            'FIGURE1 DNG':   win = FEN_Figure1_Dng()
            'FIGURE1 LIJEN': win = FEN_Figure1_LiJen()
            'ASYMM-SCAN/BY1 T30 DNG': win = FEN_AsymmScan_By1_t30_Dng(SIM_OBJECT=oSim)
            'ASYMM-2D-LARGE-NEW T06 OVERVIEW': win = FEN_Overview_t06(SIM_OBJECT=oSim, STRNAME=strname)
            'ASYMM-2D-LARGE-NEW T09 OVERVIEW': win = FEN_Overview_t09(SIM_OBJECT=oSim, STRNAME=strname)
            'ASYMM-2D-LARGE-NEW T13 OVERVIEW': win = FEN_Overview_t13(SIM_OBJECT=oSim, STRNAME=strname)
            'ASYMM-2D-LARGE-NEW T18 OVERVIEW': win = FEN_Overview_t18(SIM_OBJECT=oSim, STRNAME=strname)
            'SIM1 OVERVIEW T72': win = FEN_Overview_Sym_t72_xline(SIM_OBJECT=oSim)
            else: message, 'Figure "' + figure + '" not an option.'
        endcase
    endelse

    ;Save the plots?
    if tf_save then begin
        froot = '/home/argall/figures/'
        fbase = strlowcase(strjoin(strsplit(simname, '/', /EXTRACT), '-')) + '_' + strname
        fname = filepath(fbase, ROOT_DIR=froot)

        nWins = n_elements(strname)
        if nWins eq 1 then begin
            ;Do not use ImageMagick
            win.saveas -> SetProperty, IM_RASTER=0
            win -> Save, fname + '.png'
        endif else begin
            ;Save each plot        
            for i = 0, nWins - 1 do begin
                ;Do not use ImageMagick
                win[i].saveas -> SetProperty, IM_RASTER=0
                win[i] -> Save, fname[i] + '.png'
            endfor
        endelse
    endif
    
    if arg_present(oSim) eq 0 && obj_valid(oSim) then obj_destroy, oSim
    return, win
end