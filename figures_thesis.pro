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
;   eMap for the Asymm-3D simulation.
;-
function FigThesis_AsymmLarge2D_OhmsLaw, $
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
    theSim = 'Asymm-Large-2D'
    tIndex = 32
    xrange = 1509 + [-30, 30]
    zrange = [-20, 20]
    oSim   = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange)
    
    ;Ohm's Law
    component = 'X'
    cut       = 1509
    win       = MrSim_OhmsLaw(oSim, component, cut)
    
    ;Overview
    im_win = MrSim_ColorSlab(oSim, 'Jey', C_NAME='Ay', VERT_LINES=cut, LINE_COLOR='White')

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
function Figures_Thesis, figure, $
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
    list_of_figures = [['Asymm-Large-2D', ''], $
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
        'ASYMM-LARGE-2D OHMS LAW': win = FigThesis_AsymmLarge2D_OhmsLaw()
        else: message, 'Figure "' + figure + '" not an option.', /INFORMATIONAL
    endcase

;---------------------------------------------------------------------
; Save Figures ///////////////////////////////////////////////////////
;---------------------------------------------------------------------  
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