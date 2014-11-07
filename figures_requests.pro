;+
;   Make an array of vx-vy dist.'s that roughly covers the entire outflow jet from
;   T=29 for the Bg=3% run
;-
function FigRequests_By03_Jet_eMap, $
FNAME=fname
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim) then obj_destroy, oSim
        if max(obj_valid(win)) then obj_destroy, win
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    ;Simulation
    theSim = 'data-by0.03-NEW'
    tIndex = 29
    xrange = [740, 850]
    zrange = [-35, 35]
    oSim   = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange, /BINARY)
    
    ;Ohm's Law
;    xrange = 771 + [-70, 200]
;    centers = [[744.185181, 0.261692], $
;               [747.360107, 0.261692], $
;               [753.392517, 0.261692], $
;               [757.837402, 0.523380], $
;               [763.234863, 0.261692], $
;               [767.679749, 0.654222], $
;               [772.124695, 0.785069], $
;               [776.252075, 0.915911], $
;               [780.379517, 1.046758], $
;               [782.919495, 1.439293], $
;               [786.729431, 2.093513], $
;               [792.444275, 2.747737], $
;               [799.111694, 2.878583], $
;               [806.096558, 2.616895], $
;               [815.938904, 2.093513], $
;               [825.781189, 1.046758], $
;               [827.051147, 1.308446], $
;               [840.385925, 0.261692], $
;               [863.880493, -2.616894], $
;               [878.167725, -4.841252], $
;               [888.962524, -6.934766], $
;               [896.899902, -8.243213], $
;               [907.059692, -9.028281], $
;               [919.124512, -7.981524], $
;               [928.331848, -6.411386]]
    centers = [[744.649536, 0.332720], $
               [746.938171, 0.332720], $
               [749.226807, 0.332720], $
               [750.657227, 0.332720], $
               [752.373657, 0.332720], $
               [755.234497, 0.000001], $
               [757.809204, 0.332720], $
               [760.383911, 0.332720], $
               [763.244690, 0.332720], $
               [765.533325, 0.332720], $
               [767.535889, 0.332720], $
               [769.538452, 0.332720], $
               [772.113159, 0.332720], $
               [773.829651, 0.998158], $
               [778.120850, 1.330877], $
               [782.698120, 1.330877], $
               [784.414551, 1.996315], $
               [789.850098, 2.661753], $
               [794.999512, 1.996315], $
               [800.148926, 1.996315], $
               [806.156616, 1.996315], $
               [810.161743, 1.996315], $
               [813.594666, 1.996315], $
               [817.599792, 1.663596], $
               [823.035278, 1.663596]]
    half_width = fltarr(size(centers, /DIMENSIONS)) + 1
    im_name    = 'Uex'
    layout     = [5,5]
    
    ;Overview
    win = MrSim_eMap(oSim, 'Vx-Vy', centers, half_width, $
                     LAYOUT  = layout, $
                     IM_NAME = im_name, $
                     IM_WIN  = im_win, $
                     C_NAME  = 'Ay')
                    

    if obj_valid(oSim) then obj_destroy, oSim
    return, win
end


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
        if obj_valid(win) then obj_destroy, win
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
    !Null = MrSim_LineCut(oSim, 'vExB_mag', xline[0], /CURRENT, VCUT_RANGE=[-5, 5])
    !Null = MrSim_LineCut(oSim, 'vExB_x', xline[0], /CURRENT,   VCUT_RANGE=[-5, 5])
    !Null = MrSim_LineCut(oSim, 'vExB_y', xline[0], /CURRENT,   VCUT_RANGE=[-5, 5])
    !Null = MrSim_LineCut(oSim, 'vExB_z', xline[0], /CURRENT,   VCUT_RANGE=[-5, 5])
    
    ;Adjust
    win['Color Jey']   -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
    win['Cut vExB_mag'] -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)', YTITLE='|v$\downExB$|'
    win['Cut vExB_x']   -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)', YTITLE='v$\downExB,x$'
    win['Cut vExB_y']   -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)', YTITLE='v$\downExB,y$'
    win['Cut vExB_z']   -> SetProperty, TITLE='', YTITLE='v$\downExB,z$'
    
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
    list_of_figures = [['SIM1 EXB XLINE',      'Vertical cut of ExB through the X-line at twci=19 for the Sim1 run'], $
                       ['By0.03-NEW Jet eMap', 'Make an array of vx-vy dist.s that roughly covers the entire outflow jet from T=29 for the Bg=3% run']]
    
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
        'BY0.03-NEW JET EMAP': win = FigRequests_By03_Jet_eMap()
        'SIM1 EXB XLINE':      win = Requests_Sim1_ExB_Xline()
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