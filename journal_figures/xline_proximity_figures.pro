;+
;   2D color plot of Jey surrounding the X-line near the peak reconnection rate.
;-
function XLP_Figure2_Scaled
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(Fig2)   then obj_destroy, Fig2
        if obj_valid(difwin) then obj_destroy, difwin
        if obj_valid(cwin)   then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    

;---------------------------------------------------------------------
; 2D Sim, t=32 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    dir  = '/home/argall/Work/AssymectricSim/asymm-2D-large/'
    time = 32
    oSim = obj_new('MrSim2D', time, ION_SCALE=1, XRANGE=[2,-2], ZRANGE=150.9+[-3.0,3.0], $
                              DIRECTORY=dir+'data', INFO_FILE=dir+'info', $
                              COORD_SYSTEM='MAGNETOPAUSE', /MVA_FRAME)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create cuts of Bx, By, ni, Uix, Ez
    Fig2 = MrSim_XProximity([150.9, 150.3, 148.4], SIM_OBJECT=oSim)
    Fig2 -> Refresh, /DISABLE

    ;Bx
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = Fig2['Cut Bx']
    oTemp -> GetData, Bx
    Bx *= 10.0
    oTemp -> SetData, temporary(Bx)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*10.0, YTICKFORMAT='(f4.1)', YTITLE='B$\downL$!Cx10$\up-1$', YTICKINTERVAL=2.5

    ;By
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = Fig2['Cut By']
    oTemp -> GetData, By
    By *= 10.0
    oTemp -> SetData, temporary(By)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*10.0, YTICKFORMAT='(f4.1)', YTITLE='B$\downM$!Cx10$\up-1$'

    ;ni
    Fig2['Cut ni'] -> SetProperty, YTITLE='n$\downi$', YRANGE=[0,1.2], YTICKINTERVAL = 0.5

    ;Uix
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = Fig2['Cut Uix']
    oTemp -> GetData, Uix
    Uix *= 100.0
    oTemp -> SetData, temporary(Uix)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*100.0, YTICKFORMAT='(f4.1)', YTITLE='U$\downiL$!Cx10$\up-2$', YTICKINTERVAL=1

    ;Ez
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = Fig2['Cut Ez']
    oTemp -> GetData, Ez
    Ez *= 100.0
    oTemp -> SetData, temporary(Ez)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*100.0, YTICKFORMAT='(f4.1)', YTITLE='E$\downN$!Cx10$\up-2$', YTICKINTERVAL=1

    ;Set the layout
    ;   Move plots to the second column.
    Fig2 -> SetProperty, XSIZE=1000, YSIZE=500
    for i = 0, 4 do Fig2[i] -> SetLayout, [2,i+1]

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Add a 2D color plot
    ;   - Will automatically be added to location [1,1]
    !Null = MrSim_ColorSlab('Jey', C_NAME='Ay', /CURRENT, $
                            HORIZ_LINE=[150.9, 150.3, 148.4], SIM_OBJECT=oSim)

    ;Move the colorbar to the bottom of the figure
    Fig2['CB: Color Jey'] -> SetProperty, NAME='t=32 CB: Color Jey', $
                                          CBLOCATION='BOTTOM', OFFSET=3.5, WIDTH=1, $
                                          XTICKS=1, XTICKFORMAT='(f5.2)'
    
    ;Change the names of all figures be pre-pending a "t=32"
    graphics = Fig2 -> Get(/ALL)
    foreach gfx, graphics do begin
        name = gfx -> GetName()
        gfx.NAME = 't=32 ' + name
    endforeach

;---------------------------------------------------------------------
; 2D Sim, t=90 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Change to t=90
    ;   - Re-center around the X-line
    oSim -> SetProperty, TIME=90, XRANGE=[3,-3], ZRANGE=153.3+[-7.0,7.0]

;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Line plots of cuts surrounding the X-line
    difwin = MrSim_XProximity([153.3, 150.8, 148.0], SIM_OBJECT=oSim)
    difwin -> Refresh, /DISABLE

    ;Bx
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = difwin['Cut Bx']
    oTemp -> GetData, Bx
    Bx *= 10.0
    oTemp -> SetData, temporary(Bx)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*10.0, YTICKFORMAT='(f4.1)', YTITLE='B$\downL$!Cx10$\up-1$', $
                          YTICKINTERVAL=2.5

    ;By
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = difwin['Cut By']
    oTemp -> SetProperty, YTICKFORMAT='(f4.1)', YTITLE='B$\downM$'

    ;ni
    difwin['Cut ni'] -> SetProperty, YTITLE='n$\downi$', YRANGE=[0,1], YTICKINTERVAL = 0.5

    ;Uix
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = difwin['Cut Uix']
    oTemp -> GetData, Uix
    Uix *= 10.0
    oTemp -> SetData, temporary(Uix)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*10.0, YTICKFORMAT='(f4.1)', YTITLE='U$\downiL$!Cx10$\up-1$', $
                          YTICKINTERVAL=0.2

    ;Ez
    ;   - Rescale the data
    ;   - Change Properties
    oTemp = difwin['Cut Ez']
    oTemp -> GetData, Ez
    Ez *= 10.0
    oTemp -> SetData, temporary(Ez)
    oTemp -> GetProperty, YRANGE=yrange
    oTemp -> SetProperty, YRANGE=yrange*10.0, YTICKFORMAT='(f4.1)', YTITLE='E$\downN$!Cx10$\up-1$', $
                          YTICKINTERVAL=0.2
    
    ;We want to move the plots to window Fig2
    ;   - Resize to be similar to Fig2
    difwin -> SetProperty, XSIZE=1000, YSIZE=500, XGAP=4, LAYOUT=[4,5]
    
    ;Move the figures to the 4th column
    graphics = difwin -> Get(/ALL, ISA='MRPLOT', COUNT=nPlots)
    for i = 0, nPlots - 1 do begin
        ;Move to the 4th column
        thisPlot = graphics[i]
        thisPlot -> SetLayout, [4,i+1]
        
        ;Rename by prepending a "t=90"
        name          = thisPlot -> GetName()
        thisPlot.NAME = 't=90 ' + name
        
        ;Switch windows
        thisPlot -> SwitchWindows, Fig2
    endfor

    ;Move the legend (last remaining item)
    theLegend      = difwin[0]
    name           = theLegend -> GetName()
    theLegend.NAME = 't=90 ' + name
    theLegend     -> SwitchWindows, Fig2
    
    ;Destroy the empty DIFWIN window
    obj_destroy, difwin
    
    ;Refocus on Fig2
    Fig2 -> SetCurrent

;---------------------------------------------------------------------
; 2D Color Plot  /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Plot Jey
    cwin = MrSim_ColorSlab('Jey', C_NAME='Ay', $
                           HORIZ_LINE=[153.3, 150.8, 148.0], SIM_OBJECT=oSim)
    cwin -> Refresh, /DISABLE
    
    ;Move to location [1,2]
    cwin['Color Jey'] -> SetLayout, [1,2]

    ;Rename and move graphics into FIG2
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do begin
        name     = gfx -> GetName()
        gfx.NAME = 't=90 ' + name
        gfx     -> SwitchWindows, Fig2
    endforeach

    ;Move the colorbar to the bottom of the figure
    Fig2['t=90 CB: Color Jey'] -> SetProperty, CBLOCATION='BOTTOM', OFFSET=3.5, WIDTH=1, $
                                               XTICKS=1, XTICKFORMAT='(f5.2)'

    ;Delete CWIN and refocus FIG2
    cwin -> Cleanup
    Fig2 -> SetCurrent
    
    ;Move the color plot to location [3,1]
    Fig2['t=90 Color Jey'] -> SetLayout, [3,1]
    
;---------------------------------------------------------------------
; Set Positions //////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Set the column widths
    fig2 -> SetProperty, XGAP=[11, 8, 11], COL_WIDTH=[0.15, 0.35, 0.15, 0.35]
    fig2 -> SetGlobal, CHARSIZE=2
    
    ;Set up the positions
    pos = MrLayout([4,5], WDIMS=[1100,500], XGAP=[11, 4, 11], YGAP=0.5, $
                   COL_WIDTH=[0.15,0.35,0.15,0.35], CHARSIZE=2)

    ;Change the position of the color image
    position = [pos[0,12], pos[1,12], pos[2,4], pos[3,4]]
    Fig2['t=32 Color Jey'] -> SetProperty, POSITION=position
    
    ;Reposition to fill the middle three rows.
    position = [pos[[0,1],14], pos[[2,3],6]]
    Fig2['t=90 Color Jey'].POSITION=position
    
    Fig2 -> Refresh
    return, Fig2
end


;+
;   Create Figure 2: 2D simulation.
;-
function XLP_Figure2
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(Fig2)   then obj_destroy, Fig2
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
    dir  = '/home/argall/Work/AssymectricSim/asymm-2D-large/'
    time = 32
    oSim = obj_new('MrSim2D', time, ION_SCALE=1, XRANGE=[2,-2], ZRANGE=150.9+[-3.0,3.0], $
                              DIRECTORY=dir+'data', INFO_FILE=dir+'info', $
                              COORD_SYSTEM='MAGNETOPAUSE', /MVA_FRAME)
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Create cuts of Bx, By, ni, Uix, Ez
    Fig2 = MrSim_XProximity([150.9, 150.3, 148.4], SIM_OBJECT=oSim)
    Fig2 -> Refresh, /DISABLE

    ;Format the y-axes
    Fig2['Cut Bx']  -> SetProperty, YTICKINTERVAL=0.25
    Fig2['Cut ni']  -> SetProperty, YTICKINTERVAL=0.5,  YTITLE='n$\downi$', YRANGE=[0,1.2]
    Fig2['Cut Uix'] -> SetProperty, YTICKINTERVAL=0.01, YTITLE='U$\downiL$'
    Fig2['Cut Ez']  -> SetProperty, YTICKINTERVAL=0.01

    ;Set the layout
    ;   Move plots to the second column.
    Fig2 -> SetProperty, XSIZE=1200, YSIZE=500
    for i = 0, 4 do Fig2[i] -> SetLayout, [2,i+1]

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Add a 2D color plot
    ;   - Will automatically be added to location [1,1]
    !Null = MrSim_ColorSlab('Jey', C_NAME='Ay', /CURRENT, $
                            HORIZ_LINE=[150.9, 150.3, 148.4], SIM_OBJECT=oSim)

    ;Format the X-axis
    Fig2['Color Jey']  -> SetProperty, TITLE='', XTICKS=2, XTICKFORMAT='(f4.1)'
    Fig2['HLines Jey'] -> SetProperty, THICK=3, LINESTYLE=2, COLOR = ['Cyan', 'Green', 'Red']
    
    ;Move the colorbar to the bottom of the figure
    Fig2['CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', OFFSET=0.5, WIDTH=1, $
                                          XTICKS=2, FORMAT='(f5.2)'
    
    ;Change the names of all figures be pre-pending a "t=32"
    graphics = Fig2 -> Get(/ALL)
    foreach gfx, graphics do begin
        name = gfx -> GetName()
        gfx.NAME = 't=32 ' + name
    endforeach

;---------------------------------------------------------------------
; 2D Sim, t=90 ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Change to t=90
    ;   - Re-center around the X-line
    oSim -> SetProperty, TIME=90, XRANGE=[3,-3], ZRANGE=153.3+[-7.0,7.0]

;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Line plots of cuts surrounding the X-line
    difwin = MrSim_XProximity([153.3, 150.8, 148.0], SIM_OBJECT=oSim)
    difwin -> Refresh, /DISABLE

    difwin['Cut Bx']  -> SetProperty, YTICKINTERVAL=0.25
    difwin['Cut ni']  -> SetProperty, YTICKINTERVAL=0.4, YTITLE='n$\downi$', YRANGE=[0,0.85]
    difwin['Cut Uix'] -> SetProperty, YTICKINTERVAL=0.02
    difwin['Cut Ez']  -> SetProperty, YTICKINTERVAL=0.02
    
    ;We want to move the plots to window Fig2
    ;   - Resize to be similar to Fig2
    difwin -> SetProperty, LAYOUT=layout, XSIZE=xsuze, YSIZE=ysize, XGAP=xgap, YGAP=ygap
    
    ;Move the figures to the 4th column
    graphics = difwin -> Get(/ALL, ISA='MRPLOT', COUNT=nPlots)
    for i = 0, nPlots - 1 do begin
        ;Move to the 4th column
        thisPlot = graphics[i]
        thisPlot -> SetLayout, [4,i+1]
        
        ;Rename by prepending a "t=90"
        name          = thisPlot -> GetName()
        thisPlot.NAME = 't=90 ' + name
        
        ;Switch windows
        thisPlot -> SwitchWindows, Fig2
    endfor

    ;Move the legend (last remaining item)
    theLegend      = difwin[0]
    name           = theLegend -> GetName()
    theLegend.NAME = 't=90 ' + name
    theLegend     -> SwitchWindows, Fig2
    
    ;Destroy the empty DIFWIN window
    obj_destroy, difwin
    
    ;Refocus on Fig2
    Fig2 -> SetCurrent

;---------------------------------------------------------------------
; 2D Color Plot  /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Plot Jey
    cwin = MrSim_ColorSlab('Jey', C_NAME='Ay', $
                           HORIZ_LINE=[153.3, 150.8, 148.0], SIM_OBJECT=oSim)
    obj_destroy, oSim
    
    ;Turn refresh off
    cwin -> Refresh, /DISABLE
    
    ;Move to location [1,2]
    cwin['HLines Jey'] -> SetProperty, THICK=3, LINESTYLE=2, COLOR = ['Cyan', 'Green', 'Red']
    oTemp = cwin['Color Jey']
    oTemp -> SetLayout, [1,2]
    oTemp -> SetProperty, XTICKS=2, TITLE=''

    ;Rename and move graphics into FIG2
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do begin
        name     = gfx -> GetName()
        gfx.NAME = 't=90 ' + name
        gfx     -> SwitchWindows, Fig2
    endforeach

    ;Move the colorbar to the bottom of the figure
    Fig2['t=90 CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', OFFSET=0.5, WIDTH=1, $
                                               XTICKS=2, FORMAT='(f5.2)'

    ;Delete CWIN and refocus FIG2
    cwin -> Cleanup
    Fig2 -> SetCurrent
    
    ;Move the color plot to location [3,1]
    Fig2['t=90 Color Jey'] -> SetLayout, [3,1]
    
;---------------------------------------------------------------------
; Set Positions //////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Set the column widths
    Fig2 -> SetProperty, COL_WIDTH=col_width, OYMARGIN=oymargin, XGAP=xgap, YGAP=ygap
    Fig2 -> SetGlobal, CHARSIZE=charsize

    ;Set up the positions
    pos = MrLayout(layout, OYMARGIN=oymargin, XGAP=xgap, YGAP=ygap, $
                   COL_WIDTH=col_width, CHARSIZE=charsize)

    ;Change the position of the color image
    top      = pos[1,0] + (pos[3,0] - pos[1,0])*0.25
    position = [pos[0,16], pos[1,16], pos[2,4], top]
    Fig2['t=32 Color Jey'] -> SetProperty, POSITION=position

    ;Reposition to fill the middle three rows.
    top      = pos[1,2] + (pos[3, 2] - pos[1, 2])*0.25
    position = [pos[0,18], pos[1,18], pos[2,2], top]
    Fig2['t=90 Color Jey'].POSITION=position

    ;Adjust the legend
    ;t=32
    oTemp = Fig2['t=32 Legend: Cuts']
    oTemp.LOCATION=4
    oTemp -> GetProperty, TITLES=titles, LOCATION=location
    titles = strsplit(titles, 'L=', /EXTRACT)
    location -= [0, 0.015]
    oTemp -> SetProperty, TITLES=titles, LOCATION=location
    
    ;t=90
    oTemp = Fig2['t=90 Legend: Cuts']
    oTemp.LOCATION=4
    oTemp -> GetProperty, TITLES=titles, LOCATION=location
    titles = strsplit(titles, 'L=', /EXTRACT)
    location -= [0, 0.015]
    oTemp -> SetProperty, TITLES=titles, LOCATION=location
    
    ;Put "(a)" and "(b)" on the figure
    aText = MrText(0.03, 0.9, '(a)', CHARSIZE=2, /CURRENT, /NORMAL)
    bText = MrText(0.50, 0.9, '(b)', CHARSIZE=2, /CURRENT, /NORMAL)

    Fig2 -> Refresh
    return, Fig2
end


;+
;   Create Figure 3: 3D simulation.
;-
function XLP_Figure3
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(Fig3)   then obj_destroy, Fig3
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
; 3D Sim, t=120816, y=1400 ///////////////////////////////////////////
;---------------------------------------------------------------------
    dir    = '/home/argall/Work/AssymectricSim/asymm-3D/'
    time   = 120816
    yslice = 860
    oSim   = obj_new('MrSim3D', time, yslice, /ION_SCALE, XRANGE=[5,-4], ZRANGE=45.3+[-15,15], $
                                              DIRECTORY=dir+'data', INFO_FILE=dir+'info', $
                                              COORD_SYSTEM='MAGNETOPAUSE', /MVA_FRAME)

;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    Fig3 = MrSim_XProximity([45.2, 43.7, 35.2], SIM_OBJECT=oSim)
    Fig3 -> Refresh, /DISABLE
    
    ;Change properties
    Fig3['Cut Bx'] -> SetProperty, YTICKINTERVAL=0.3
    Fig3['Cut Uix'].YTICKINTERVAL=0.02
    Fig3['Cut Ez'].YTICKINTERVAL=0.02

    ;Set the layout
    ;   Move plots to the second column.
    Fig3-> SetProperty, XSIZE=xsize, YSIZE=ysize
    for i = 0, 4 do Fig3[i] -> SetLayout, [2,i+1]

;---------------------------------------------------------------------
; Jey Color Plot /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Aadd a 2D color plot
    ;   - Automatically added to location [1,1]
    cwin = MrSim_ColorSlab('Jey', /CURRENT, HORIZ_LINE=[45.2, 43.7, 35.2], SIM_OBJECT=oSim)
    
    ;Split the title into two lines. Change the tick interval.
    cwin['Color Jey'] -> SetProperty, TITLE='', XTICKINTERVAL=3, XTICKS=2
    cwin['HLines Jey'] -> SetProperty, THICK=3, LINESTYLE=2, COLOR=['Cyan', 'Green', 'Red']
    
    ;Move the colorbar to the bottom of the plot. Adjust properties.
    cwin['CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', WIDTH=1, OFFSET=0.5, $
                                          XTICKS=2, FORMAT='(f5.2)'
    
    ;Change the names of all figures by pre-pending a "y=1400"
    graphics = Fig3 -> Get(/ALL)
    foreach gfx, graphics do begin
        name = gfx -> GetName()
        gfx.NAME = 'y=860 ' + name
    endforeach

;---------------------------------------------------------------------
; 3D Sim, t=120816, y=860 ///////////////////////////////////////////
;---------------------------------------------------------------------
    ;Change to y=860
    ;   - Re-center around the x-line
    oSim -> SetProperty, YSLICE=1400, XRANGE=[8,-6], ZRANGE=44+[-15,15]
    
;---------------------------------------------------------------------
; Cuts within the Exhaust ////////////////////////////////////////////
;---------------------------------------------------------------------
    difwin = MrSim_XProximity([44.7, 39.0, 32.0], SIM_OBJECT=oSim)
    difwin -> Refresh, /DISABLE

    ;Change properties
    difwin['Cut Bx'] -> SetProperty, YTICKINTERVAL=0.3
    difwin['Cut By'].YTICKINTERVAL=0.1
    
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
        thisPlot.NAME = 'y=1400 ' + name
        
        ;Switch windows
        thisPlot -> SwitchWindows, Fig3
    endfor

    ;Move the legend (last remaining item)
    theLegend      = difwin[0]
    name           = theLegend -> GetName()
    theLegend.NAME = 'y=1400 ' + name
    theLegend     -> SwitchWindows, Fig3
    
    ;Destroy the empty DIFWIN window
    obj_destroy, difwin
    
    ;Refocus on Fig3
    Fig3 -> SetCurrent

;---------------------------------------------------------------------
; 2D Color Plot  /////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Plot Jey
    cwin = MrSim_ColorSlab('Jey', HORIZ_LINE=[44.7, 39.0, 32.0], SIM_OBJECT=oSim)
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
        gfx.NAME = 'y=1400 ' + name
        gfx     -> SwitchWindows, Fig3
    endforeach

    ;Move the colorbar to the bottom of the figure
    Fig3['y=1400 CB: Color Jey'] -> SetProperty, CBLOCATION='TOP', OFFSET=0.5, WIDTH=1, $
                                               XTICKS=2, XTICKFORMAT='(f5.2)'

    ;Delete CWIN and refocus FIG2
    cwin -> Cleanup
    Fig3 -> SetCurrent
    
    ;Move the color plot to location [3,1]
    Fig3['y=1400 Color Jey'] -> SetLayout, [3,1]
    
;---------------------------------------------------------------------
; Set Positions //////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Set the column widths
    Fig3 -> SetProperty, COL_WIDTH=col_width, OYMARGIN=oymargin, XGAP=xgap, YGAP=ygap
    Fig3 -> SetGlobal, CHARSIZE=charsize

    ;Set up the positions
    Fig3 -> SetCurrent
    pos = MrLayout(layout, OYMARGIN=oymargin, XGAP=xgap, YGAP=yap, $
                   COL_WIDTH=col_width, CHARSIZE=charsize)

    ;Change the position of the color image
    top      = Fig3['y=860 Cut Bx'].position
    top      = top[1] + (top[3] - top[1]) * 0.25
    position = [pos[0,16], pos[1,16], pos[2,0], top]
    Fig3['y=860 Color Jey'] -> SetProperty, POSITION=position
    
    ;Reposition to fill the middle three rows.
    top      = Fig3['y=1400 Cut Bx'].position
    top      = top[1] + (top[3] - top[1]) * 0.25
    position = [pos[0,18], pos[1,18], pos[2,2], top]
    Fig3['y=1400 Color Jey'].POSITION=position

    ;Adjust the legend
    ;y=1400
    oTemp = Fig3['y=860 Legend: Cuts']
    oTemp.LOCATION=4
    oTemp -> GetProperty, TITLES=titles, LOCATION=location
    titles = strsplit(titles, 'L=', /EXTRACT)
    location -= [0,0.015]
    oTemp -> SetProperty, TITLES=titles, LOCATION=location
    
    ;y=860
    oTemp = Fig3['y=1400 Legend: Cuts']
    oTemp.LOCATION=4
    oTemp -> GetProperty, TITLES=titles, LOCATION=location
    titles = strsplit(titles, 'L=', /EXTRACT)
    location -= [0, 0.015]
    oTemp -> SetProperty, TITLES=titles, LOCATION=location
    
    ;Put "(a)" and "(b)" on the figure
    aText = MrText(0.05, 0.92, '(a)', CHARSIZE=2, /CURRENT, /NORMAL)
    bText = MrText(0.52, 0.92, '(b)', CHARSIZE=2, /CURRENT, /NORMAL)
    
    Fig3 -> Refresh
    return, Fig3
end


function XLP_Figure1_GEM
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(Fig3)   then obj_destroy, Fig3
        if obj_valid(difwin) then obj_destroy, difwin
        if obj_valid(cwin)   then obj_destroy, cwin
        void = cgErrorMSG()
        return, obj_new()
    endif
    
    
    charsize  = 2.0
    charthick = 2.0
    col_width = [0.3, 0.7]
    layout    = [2,5]
    thick     = 2
    xgaps     = 8
    xsize     = 800
    ysize     = 600
;---------------------------------------------------------------------
; t=32 ///////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    dir  = '/home/argall/Work/AssymectricSim/asymm-2D-large/'
    time = 32
    oSim = obj_new('MrSim2D', time, ION_SCALE=1, XRANGE=[2,-2], ZRANGE=150.9+[-3.0,3.0], $
                              DIRECTORY=dir+'data', INFO_FILE=dir+'info', $
                              COORD_SYSTEM='MAGNETOPAUSE', /MVA_FRAME)
                              
;---------------------------------------------------------------------
; Jey ////////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    cwin = MrSim_ColorSlab('Jey', C_NAME='Ay', HORIZ_LINE=[150.9, 150.3, 148.4], SIM_OBJECT=oSim)
    cwin -> SetProperty, XSIZE=xsize, YSIZE=ysize, OYMARGIN=[3,6]

    ;Set properties 
    cwin['Color Jey'] -> SetProperty, XTICKLEN=-0.02, YTICKLEN=-0.02
    cwin['HLines Jey'] -> SetProperty, COLOR=['Cyan', 'Green', 'Red'], THICK=thick
    cwin['CB: Color Jey'] -> SetProperty, WIDTH=1, CBLOCATION='BOTTOM', OFFSET=4, $
                                          XTICKS=1, TEXTTHICK=thick
                              
;---------------------------------------------------------------------
; Cuts Across Exhaust ////////////////////////////////////////////////
;---------------------------------------------------------------------
    Fig2 = MrSim_XProximity([150.9, 150.3, 148.4], SIM_OBJECT=oSim)
    Fig2 -> Refresh, /DISABLE
    
    ;Resize the window
    Fig2 -> SetProperty, COL_WIDTH=col_width, LAYOUT=layout, XSIZE=xsize, YSIZE=ysize
    
    ;Move plots to column 2
    graphics = Fig2 -> Get(/ALL, ISA='MRPLOT', COUNT=nPlots)
    for i = 0, nPlots-1 do graphics[i] -> SetLayout, [2,i+1]
    
    ;Set properties
    Fig2['Cut Bx'] -> SetProperty,  YTICKINTERVAL=0.3
    Fig2['Cut ni'] -> SetProperty,  YTICKINTERVAL=0.5, YRANGE=[0,1.2]
    Fig2['Cut Uix'] -> SetProperty, YTICKINTERVAL=0.01
    Fig2['Cut Ez'] -> SetProperty,  YTICKINTERVAL=0.01
    
    ;Change the position of the legend
    oTemp = Fig2['Legend: Cuts']
    loc   = oTemp.LOCATION
    location = loc + [0, -0.02]
    oTemp -> SetProperty, LOCATION=location, CHARSIZE=charsize

;---------------------------------------------------------------------
; Combine ////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Move the 2D color plot of Jey into the window with the line cuts.
    graphics = cwin -> Get(/ALL)
    foreach gfx, graphics do gfx -> SwitchWindows, Fig2
    cwin -> Cleanup
    
    ;Set Global Properties
    Fig2 -> SetGlobal,  CHARSIZE=charsize, CHARTHICK=charthick
    
    ;Reposition
    pos = MrLayout(layout, WDIMS=[xsize, ysize], COL_WIDTH=col_width)
    bottom = (pos[3,6] + pos[1,6]) / 2.0
    top    = (pos[3,0] + pos[1,0]) / 2.0
    position = [pos[0,6], bottom, pos[2,2], top]
    Fig2['Color Jey'].POSITION=position
    

    return, win
end




;+
;   Create the desired figure.
;
; Params:
;       FIGURE:         in, optional, type=string
;                       Figure number of the figure to be created.
;-
function XLine_Proximity_Figures, figure
    compile_opt strictarr
    on_error, 2

    ;Create the figure    
    case figure of
        '2': win = XLP_Figure2()
        '3': win = XLP_Figure3()
        else: message, 'Figure "' + figure + '" not an option.', /INFORMATIONAL
    endcase
    
    return, win
end