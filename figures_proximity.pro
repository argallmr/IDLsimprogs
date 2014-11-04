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
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function Figures_Proximity, figure
    compile_opt strictarr
    on_error, 2

    theFig = strtrim(figure, 2)

    ;Create the figure    
    case figure of
        '1': win = Prox_Figure1()
        else: message, 'Figure "' + figure + '" not an option.', /INFORMATIONAL
    endcase
    
    return, win
end