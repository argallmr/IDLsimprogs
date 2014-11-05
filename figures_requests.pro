;+
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function Requests_Sim1_ExB_Xline, simname, tIndex, yslice, $
LEFT=left, $
EMAP=eMap, $
STRNAME=strname, $
SIM_OBJECT=oSim
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        void = cgErrorMSG()
        return, !Null
    endif
    
    theSim = 'Sim1'
    time   = 76
    xrange = 850 + [-50, 50]
    zrange = [-20, 20]
    xline  = [847.7, 0]
    
    ;Create the simulation object.
    oSim = MrSim_Create(theSim, time, XRANGE=xrange, ZRANGE=zrange, /BINARY)
    
    ;Create an overview plot
    win = MrSim_ColorSlab(oSim, 'Jey', C_NAME='Ay')
    win -> Refresh, /DISABLE
    
    ;Change window size
    win.YSIZE = 600
    win.YGAP  = 0
    
    ;Create the line plots
    !Null = MrSim_LineCut(oSim, 'ExB_mag', xline[0], /CURRENT, VCUT_RANGE=[-5, 5])
    !Null = MrSim_LineCut(oSim, 'ExB_x', xline[0], /CURRENT,   VCUT_RANGE=[-5, 5])
    !Null = MrSim_LineCut(oSim, 'ExB_y', xline[0], /CURRENT,   VCUT_RANGE=[-5, 5])
    !Null = MrSim_LineCut(oSim, 'ExB_z', xline[0], /CURRENT,   VCUT_RANGE=[-5, 5])
    
    ;Adjust
    win['Color Jey']   -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
    win['Cut ExB_mag'] -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
    win['Cut ExB_x']   -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
    win['Cut ExB_y']   -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
    win['Cut ExB_z']   -> SetProperty, TITLE=''
    
    ;Draw a vertical line through the x-line
    !Null = MrPlotS([xline[0], xline[0]], [-5,5], TARGET=win['Color Jey'], COLOR='White')
    
    ;Return
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
function Figures_Requests, figure, $
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
    
    _figure = strupcase(figure)
    tf_save = keyword_set(tf_save)
    
;---------------------------------------------------------------------
; Info ///////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    
    ;Current list of figures
    list_of_figures = [['SIM1 EXB XLINE', 'Vertical cut of ExB through the X-line at twci=19 for the Sim1 run'], $
                       ['Figure1 Dng          ', '']]
    
    ;Print the list of figures?
    if keyword_set(show_figs) then begin
        len = strtrim(max(strlen(list_of_figures[0,*])), 2)
        print, 'FIGURES', 'DESCRIPTION', FORMAT='(a' + len + ', 4x, a11)'
        print, list_of_figures, FORMAT='(a-' + len + ', 4x, a0)'
        return, obj_new()
    endif
    
;-------------------------------------------------------
; Create Figures ///////////////////////////////////////
;-------------------------------------------------------

    ;Create the figure    
    case _figure of
        'SIM1 EXB XLINE': win = Requests_Sim1_ExB_Xline()
        else: message, 'Figure "' + figure + '" not an option.'
    endcase

;-------------------------------------------------------
; Save /////////////////////////////////////////////////
;-------------------------------------------------------
    if tf_save then begin
;        simname = oSim.SIMNAME
;        if n_elements(strname) eq 0 then strname = strjoin(strsplit(strlowcase(_figure), ' ', /EXTRACT), '_')
;        if strpos(simname, '/') ne -1 $
;            then fbase = strlowcase(strjoin(strsplit(simname, '/', /EXTRACT), '-')) + '_' + strname $
;            else fbase = strlowcase(simname) + '_' + strname
        
        froot = '/home/argall/figures/'
        fname = filepath(strname, ROOT_DIR=froot)

        nWins = n_elements(strname)
        if nWins eq 1 then begin
            ;Do not use ImageMagick
            win.saveas -> SetProperty, IM_RASTER=0
            win -> Save, fname
        endif else begin
            ;Save each plot        
            for i = 0, nWins - 1 do begin
                ;Do not use ImageMagick
                win[i].saveas -> SetProperty, IM_RASTER=0
                win[i] -> Save, fname[i]
            endfor
        endelse
    endif
    
    if arg_present(oSim) eq 0 && obj_valid(oSim) then obj_destroy, oSim
    return, win
end