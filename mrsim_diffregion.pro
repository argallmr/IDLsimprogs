; docformat = 'rst'
;
; NAME:
;    MrSim_DiffRegion
;
; PURPOSE:
;+
;    The purpose of this program is to make plots of::
;       NAME        DESCRIPTION
;       Bx          Reconnecting component of the magnetic field
;       By          Out-of-Plane component of the magnetic field (Hall field)
;       ni          Ion density
;       Uix         Ion outflow velocity
;       Ez          Normal component of the electric field
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Examples:
;    IDL> density_asym
;
; :Params:
;       CUTS:                   in, required, type=int/intarr
;                               The X-location where vertical cuts are to be taken.
;       TIME:                   in, required, type=int
;                               The simulation time for which at which to create the plot.
;                                   Igored if `SIM_OBJECT` is given.
;
; :Keywords:
;       C_NAME:                 in, optional, type=string, default=''
;                               The name of a quantity to be plotted in color.
;       CUT_ZRANGE:             in, optional, type=intarr(2), default=`ZRANGE`
;                               The vertical range over which the cut will be taken. Units
;                                   are those specified by `SIM_OBJECT`.
;       NLEVELS:                in, optional, type=int, default=15
;                               The number of contour lines to draw
;       OFILENAME:              out, optional, type=string, default=''
;                               If provided, graphics will be output to a file by this
;                                   name. The file extension is used to determine which
;                                   type of image to create. "PNG" is the default.
;       SIM_OBJECT:             out, optional, type=object
;                               The "Sim_Object" (or subclass) reference containing the 
;                                   simulation data and domain information used in
;                                   creating the plot. If not provided, one will be
;                                   created according to the `SIM3D` keyword.
;       SIM3D:                  in, optional, type=boolean, default=0
;                               If set, and `SIM_OBJECT` is not given, then a 3D
;                                   simulation object will be created. The default is
;                                   to make 2D simulation object.
;       _REF_EXTRA:             in, optional, type=structure
;                               Any keyword accepted Sim_Object::Init(), or any of its
;                                   subclasses is also accepted here. Ignored of
;                                   `SIM_OBJECT` is present.
;
; :Uses:
;   Uses the following external programs::
;       error_message.pro (Coyote Graphics)
;       sim_color_plot.pro
;       MrPlot__define.pro
;       getIndexRange.pro
;       MrReadSim__define.pro
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
;       2014/07/25  -   Written by Matthew Argall
;-
function MrSim_DiffRegion, cuts, time, $
CUT_RANGE = cut_range, $
C_NAME = c_name, $
OFILENAME = ofilename, $
NLEVELS = nlevels, $
SIM_OBJECT = oSim, $
SIM3D = Sim3D, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(difWin) then obj_destroy, difWin
        void = cgErrorMSG()
        return, obj_new()
    endif

;-------------------------------------------------------
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------

    ;Set defaults
    nCuts = n_elements(cuts)
    Sim3D = keyword_set(Sim3D)
    if n_elements(c_name) eq 0 then c_name = ''
    if n_elements(ofilename) eq 0 then ofilename = ''
    if n_elements(nlevels) eq 0 then nlevels = 15

    ;Create a simulation object
    if n_elements(oSim) eq 0 then begin
        if Sim3D then sim_class = 'MrSim3D' else sim_class = 'MrReadSim'
        oSim = obj_new(sim_class, time, DIRECTORY=directory, NSMOOTH=nsmooth, $
                                  XRANGE=xrange, ZRANGE=zrange)
        if obj_valid(oSim) eq 0 then return, obj_new()
    endif else begin
        if obj_valid(oSim) eq 0 || obj_isa(oSim, 'MrSim') eq 0 then $
            message, 'SIM_OBJECT must valid SIM_OBJECT object.'
    endelse

    ;Get/Make the window
    if ofilename eq '' then buffer = 0 else buffer = 1
    difWin = MrWindow(XSIZE=500, YSIZE=700, BUFFER=buffer, REFRESH=0, YGAP=0.5, OXMARGIN=[10,11])
    
;-------------------------------------------------------
; Create Plot Annotations //////////////////////////////
;-------------------------------------------------------
    oSim -> GetProperty, TIME=time, YSLICE=yslice, COORD_SYSTEM=coord_system, AXIS_LABELS=axLabls
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    units = MrSim_Rename(units, /SUBSCRIPT)
    
    ;HORIZONTAL CUTS
    if coord_system eq 'MAGNETOPAUSE' then begin
        horizontal  = 1
        horiz_lines = cuts
        legendLabl  = axLabls[2]
        if n_elements(cut_range) gt 0 then hRange = cut_range
        
        ;Location of the horizontal cuts
        if obj_isa(oSim, 'MRSIM3D') then begin
            if nCuts gt 1 $
                then cutLoc = string(FORMAT='(%"' + axLabls[1] + '=%0.1f%s")', yslice, units) $
                else cutLoc = string(FORMAT='(%"[' + strjoin(axLabls[1:2], ',') + ']=[%0.1f, %0.1f]%s")', cuts, yslice, units)
        endif else begin
            if nCuts gt 1 $
                then cutLoc = '' $
                else cutLoc = axLabls[2] + '=' + string(cuts, FORMAT='(i0.1)') + units
        endelse
        
    ;VERTICAL CUTS
    endif else begin
        horizontal = 0
        vert_lines = cuts
        legendLabl = axLabls[0]
        if n_elements(cut_range) gt 0 then vRange = cut_range
        
        ;Location of the vertical cuts
        if obj_isa(oSim, 'MRSIM3D') then begin
            if nCuts gt 1 $
                then cutLoc = string(FORMAT='(%"' + axLabls[1] + '=%0.1f%s")', yslice, units) $
                else cutLoc = string(FORMAT='(%"[' + strjoin(axLabls[0:1], ',') + ']=[%0.1f, %0.1f]%s")', cuts, yslice, units)
        endif else begin
            if nCuts gt 1 $
                then cutLoc = '' $
                else cutLoc = axLabls[0] + '=' + string(cuts, FORMAT='(i0.1)') + units
        endelse
    endelse

;-------------------------------------------------------
; Create a 2D Color Plot? //////////////////////////////
;-------------------------------------------------------
    if c_name ne '' then begin
        !Null = MrSim_ColorSlab(c_name, time, /CURRENT, $
                                NLEVELS = nlevels, $
                                VERT_LINES = vert_lines, $
                                HORIZ_LINES = horiz_lines, $
                                SIM_OBJECT = oSim)

        ;Change some properties
        difWin['Color ' + c_name] -> SetProperty, XTICKFORMAT='(a1)', XTITLE=''
    endif

;-------------------------------------------------------
; Create Line Cut Plot /////////////////////////////////
;-------------------------------------------------------
    ;Create a title
    if n_elements(dtxwci) eq 0 $
        then title = 't$\downindex$=' + string(time, FORMAT='(f0.1)') + '  ' + cutLoc $
        else title = 't$\Omega$$\downci$=' + string(time*dtxwci, FORMAT='(f0.1)') + '  ' + cutLoc

    ;Plot the individual cuts.
    !Null = MrSim_LineCut('Bx',  cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('By',  cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('ni',  cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('ne',  cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('Uix', cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('Uex', cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('Ui_perp', cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('Ue_perp', cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('ExB_mag', cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)
    !Null = MrSim_LineCut('Ez',  cuts, /CURRENT, HORIZONTAL=horizontal, SIM_OBJECT=oSim, HCUT_RANGE=hRange, VCUT_RANGE=vRange)

    ;Destroy the object.
    if arg_present(oSim) eq 0 then obj_destroy, oSim

;-------------------------------------------------------
; Organize Plots ///////////////////////////////////////
;-------------------------------------------------------
    ;Overplot ne onto ni
    ne_plot = difWin['Cut ne']
    ni_plot = difWin['Cut ni']
    ne_plot -> Overplot, ni_plot
    ne_plot.COLOR = 'Red'
    ni_plot.COLOR = 'Blue'
    
    ;Create a legend
    !Null = MrLegend(TARGET=ni_plot, /CURRENT, LOCATION=8, NAME='Legend: Density', $
                     TITLES=['n$\downi$', 'n$\downe$'], TCOLORS=['Blue', 'Red'], LENGTH=0)
    
    ;Overplot Uex onto Uix
    ;   - Flux will be overplotted later. Create legend then.
    Uix_plot = difWin['Cut Uix']
    Uex_plot = difWin['Cut Uex']
    Uex_plot -> Overplot, Uix_plot
    Uex_plot.COLOR  = 'Red'
    Uix_plot.COLOR  = 'Blue'
    
    ;Overplot Ue_perp and v_ExB onto Ui_perp
    Ue_perp_plot = difWin['Cut Ue_perp']
    Ui_perp_plot = difWin['Cut Ui_perp']
    ExB_plot     = difWin['Cut ExB_mag']
    Ue_perp_plot -> Overplot, Ui_perp_plot
    ExB_plot     -> Overplot, Ui_perp_plot
    Ue_perp_plot.COLOR = 'Red'
    Ue_perp_plot -> GetData, Ue_perp
    Ui_perp_plot -> SetProperty, COLOR='Blue', YRANGE=[0, max(temporary(Ue_perp))]
    
    ;Create a legend
    !Null = MrLegend(TARGET=Ui_perp_plot, /CURRENT, LOCATION=8, NAME='Legend: NonIdeal', $
                     TITLES=['U$\downExB$', 'U$\downix$', 'U$\downex$'], $
                     TCOLORS=['Black', 'Blue', 'Red'], LENGTH=0)
    
    ;Fill the holes
    difWin -> FillHoles, /TRIMLAYOUT
    

;-------------------------------------------------------
; Particle Flux ////////////////////////////////////////
;-------------------------------------------------------
    ni_plot -> GetProperty, XRANGE=xrange, ZRANGE=zrange

    ;Calculate flux
    ;   - First, gather required data
    ni_plot  -> GetData, indep, ni_data
    ne_plot  -> GetData, ne_data
    Uix_plot -> GetData, Uix_data
    Uex_plot -> GetData, Uex_data
    iFlux = temporary(ni_data) * temporary(Uix_data)
    eFlux = temporary(ne_data) * Uex_data
    
    ;Ion Flux
    ;   - Make sure to use the ranges from the ni_plot.
    ;       + For MAGNETOPAUSE orientations, the zrange is [max, min]
    iFlux_plot = MrPlot(indep, iFlux, /CURRENT, $
                        NAME='Ion Flux', $
                        OVERPLOT=Uix_plot, $
                        COLOR='Blue', $
                        LINESTYLE=2)

    ;Electron Flux
    eFlux_plot = MrPlot(indep, eFlux, /CURRENT, $
                        NAME='e Flux', $
                        OVERPLOT=Uix_plot, $
                        LINESTYLE=2, $
                        COLOR='Red')
    
    ;Set the data range
    Uix_plot.YRANGE = [min(temporary(Uex_data), MAX=maxUe), maxUe]
    
    ;Create a legend
    !Null = MrLegend(TARGET=iFlux_plot, /CURRENT, LOCATION=8, NAME='Legend: Flux', $
                     TITLES=['Uix', 'Uex', 'Ion Flux', 'e- Flux'], $
                     TCOLORS=['Blue', 'Red', 'Blue', 'Red'], COLORS=['Blue', 'Red', 'Blue', 'Red'], LINESTYLE=[0,0,2,2])

;-------------------------------------------------------
; Prepare Output ///////////////////////////////////////
;-------------------------------------------------------
    ;Get all of the plots
    allGr = difWin -> Get(/ALL, ISA=['MRIMAGE', 'MRPLOT'])
    foreach gfx, allGr do gfx -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
    allGr[0].TITLE=title
    allGr[-1] -> SetProperty, XTICKFORMAT='', XTITLE=axLabls[0] + ' (' + units + ')'
    
    ;Change the character size
    difWin -> SetGlobal, CHARSIZE=1.8
    
    ;Set title and margin
    if c_name ne '' then begin
        difWin['Color ' + c_name].TITLE=title
        difWin['CB: Color ' + c_name].WIDTH=1.5
        difWin.OXMARGIN=[10,11]
    endif
    
    ;Refresh and output, if requested.
    difWin -> Refresh
    if ofilename ne '' then difWin -> Save, ofilename
    
    return, difWin
end
