; docformat = 'rst'
;
; NAME:
;    SIM_DIFFREGION
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
;       04/10/2013  -   Written by Matthew Argall
;       05/06/2013  -   Updated to use MrReadSim__define.pro and getIndexRange.pro. Added
;                           the XSIZE and YSIZE keywords. - MRA
;       2014/01/22  -   Made the color plot optional via the COLOR keyword. - MRA
;       2014/01/29  -   The COLOR keyword accepts the name of the quantity to be plotted
;                            in color. Added the SIM_OBJECT and T_WCI keywords. - MRA
;       2014/01/31  -   Added the _REF_EXTRA keyword and removed others that it replaced.
;                           Renamed the COLOR keyword to C_NAME. Added SIM3D. - MRA
;       2014/03/05  -   Parameter CUT renamed to CUTS and can now be an array. Make use
;                           of the MrSim_LineCut routine. Added y-slice to title for 3D
;                           simulations. - MRA
;       2014/07/25  -   Renamed from Sim_DiffRegion.pro to MrSim_XProximity.pro - MRA
;-
function MrSim_XProximity, theSim, cuts, time, $
C_NAME = c_name, $
CUT_RANGE = cut_range, $
IM_NAME = im_name, $
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
        if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(drWin) then obj_destroy, drWin
        void = cgErrorMSG()
        return, obj_new()
    endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
    osim_created = 0B
    
    ;Simulation name or number?
    if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
        oSim = MrSim_Create(theSim, time, yslice, _STRICT_EXTRA=extra)
        if obj_valid(oSim) eq 0 then return, obj_new()
        osim_created = 1B
        
    ;Object?
    endif else if MrIsA(theSim, 'OBJREF') then begin
        if obj_isa(theSim, 'MRSIM2') eq 0 $
            then message, 'THESIM must be a subclass of the MrSim class.' $
            else oSim = theSim
            
    ;Somthing else
    endif else begin
        MrSim_Which
        message, 'THESIM must be a simulation name, number, or object.'
    endelse
    sim_class = obj_class(oSim)

;-------------------------------------------------------
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------

    ;Set defaults
    nCuts = n_elements(cuts)
    Sim3D = keyword_set(Sim3D)
    if n_elements(im_name)   eq 0 then im_name   = ''
    if n_elements(ofilename) eq 0 then ofilename = ''
    if n_elements(nlevels)   eq 0 then nlevels   = 15

    ;Get/Make the window
    if ofilename eq '' then buffer = 0 else buffer = 1
    drWin = MrWindow(XSIZE=500, YSIZE=700, BUFFER=buffer, REFRESH=0, YGAP=0.5, OXMARGIN=[10,3])
    
;-------------------------------------------------------
; Create Plot Annotations //////////////////////////////
;-------------------------------------------------------
    oSim -> GetProperty, TIME=time, COORD_SYSTEM=coord_system, AXIS_LABELS=axLabls, $
                         MVA_FRAME=mva_frame, YRANGE=yrange
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    units = MrSim_GDA_Format(units)
    
    ;HORIZONTAL CUTS
    if coord_system eq 'MAGNETOPAUSE' then begin
        names       = ['Bz', 'By', 'ni', 'Uiz', 'Ex']
        horizontal  = 1
        horiz_lines = cuts
        legendLabl  = axLabls[2]
        if n_elements(cut_range) gt 0 then hRange = cut_range
        
        ;Location of the horizontal cuts
        if oSim.dimension eq '3D' then begin
            oSim -> GetProperty, YRANGE=yrange
            if nCuts gt 1 $
                then cutLoc = string(FORMAT='(%"' + axLabls[1] + '=%0.1f%s")', yrange[0], units) $
                else cutLoc = string(FORMAT='(%"[' + strjoin(axLabls[1:2], ',') + ']=[%0.1f, %0.1f]%s")', cuts, yrange[0], units)
        endif else begin
            if nCuts gt 1 $
                then cutLoc = '' $
                else cutLoc = axLabls[2] + '=' + string(cuts, FORMAT='(i0.1)') + units
        endelse
        
    ;VERTICAL CUTS
    endif else begin
        names       = ['Bx', 'By', 'ni', 'Uix', 'Ez']
        horizontal = 0
        vert_lines = cuts
        legendLabl = axLabls[0]
        if n_elements(cut_range) gt 0 then vRange = cut_range
        
        ;Location of the vertical cuts
        if oSim.dimension eq '3D' then begin
            if nCuts gt 1 $
                then cutLoc = string(FORMAT='(%"' + axLabls[1] + '=%0.1f%s")', yrange[0], units) $
                else cutLoc = string(FORMAT='(%"[' + strjoin(axLabls[0:1], ',') + ']=[%0.1f, %0.1f]%s")', cuts, yrange[0], units)
        endif else begin
            if nCuts gt 1 $
                then cutLoc = '' $
                else cutLoc = axLabls[0] + '=' + string(cuts, FORMAT='(i0.1)') + units
        endelse
    endelse

;-------------------------------------------------------
;Make a Plot For Each X-Line ///////////////////////////
;-------------------------------------------------------
    ;Create a title
    if n_elements(dtxwci) eq 0 $
        then title = 't$\downindex$=' + string(time, FORMAT='(f0.1)') + '  ' + cutLoc $
        else title = 't$\Omega$$\downci$=' + string(time*dtxwci, FORMAT='(f0.1)') + '  ' + cutLoc

    ;Plot the individual cuts.
    !Null = MrSim_LineCut(oSim, names[0], cuts, /CURRENT, HORIZONTAL=horizontal, CUT_RANGE=cut_range)
    !Null = MrSim_LineCut(oSim, names[1], cuts, /CURRENT, HORIZONTAL=horizontal, CUT_RANGE=cut_range)
    !Null = MrSim_LineCut(oSim, names[2], cuts, /CURRENT, HORIZONTAL=horizontal, CUT_RANGE=cut_range)
    !Null = MrSim_LineCut(oSim, names[3], cuts, /CURRENT, HORIZONTAL=horizontal, CUT_RANGE=cut_range)
    !Null = MrSim_LineCut(oSim, names[4], cuts, /CURRENT, HORIZONTAL=horizontal, CUT_RANGE=cut_range)

    ;Do not outline or fill the legend
    drWin['LineCut: ' + names[0]] -> SetProperty, FILL_COLOR='', LINESTYLE='None'
    drWin['LineCut: ' + names[1]] -> SetProperty, FILL_COLOR='', LINESTYLE='None'
    drWin['LineCut: ' + names[2]] -> SetProperty, FILL_COLOR='', LINESTYLE='None'
    drWin['LineCut: ' + names[3]] -> SetProperty, FILL_COLOR='', LINESTYLE='None'
    drWin['LineCut: ' + names[4]] -> SetProperty, FILL_COLOR='', LINESTYLE='None'

;-------------------------------------------------------
; Create a 2D Color Plot? //////////////////////////////
;-------------------------------------------------------
    if im_name ne '' then begin
        oSim -> GetProperty, COORD_SYSTEM=coord_sys
        
    ;-------------------------------------------------------
    ; Vertical Orientation /////////////////////////////////
    ;-------------------------------------------------------
        if coord_sys eq 'MAGNETOPAUSE' then begin
            ;Create a second column
            drWin -> SetProperty, LAYOUT=[2,5], COL_WIDTH=[0.4, 0.6], XGAP=10, XSIZE=650, YSIZE=550
            
            ;Move plots into the second column
            allP = drWin -> Get(/ALL, ISA='MRPLOT')
            foreach gfx, allP do begin
                thisLay = gfx.layout
                colrow  = drWin -> ConvertLocation(thisLay[2], thisLay[0:1], /PINDEX, /TO_COLROW)
                gfx    -> SetLayout, [2, colrow[1]]
            endforeach
            
            ;Create the 2D color image
            !Null = MrSim_ColorSlab(oSim, im_name, /CURRENT, $
                                    C_NAME      = c_name, $
                                    NLEVELS     = nlevels, $
                                    VERT_LINES  = vert_lines, $
                                    HORIZ_LINES = horiz_lines)
            
            ;Relocate it to fill the first column
            position = [0.12, 0.13, 0.41, 0.82]
            drWin['Color '     + im_name] -> SetProperty, POSITION=position, TITLE=''
            drWin['CB: Color ' + im_name] -> SetProperty, CBLOCATION='TOP', OFFSET=0.5, WIDTH=1.0, XTICKS=2
        
    ;-------------------------------------------------------
    ; Horizontal Orientation ///////////////////////////////
    ;-------------------------------------------------------
        endif else begin
            
            ;Create the 2D color image
            !Null = MrSim_ColorSlab(oSim, im_name, /CURRENT, $
                                    C_NAME      = c_name, $
                                    NLEVELS     = nlevels, $
                                    VERT_LINES  = vert_lines, $
                                    HORIZ_LINES = horiz_lines)
            
            ;Move it to location [1,1]
            drWin['Color ' + im_name] -> SetLayout, [1,1]
            drWin['Color '     + im_name].TITLE=title
            drWin['CB: Color ' + im_name].WIDTH=1.5
            drWin -> SetProperty, XGAP=[5, 0.5, 0.5, 0.5, 0.5]
        endelse
        
    ;-------------------------------------------------------
    ; ColorSlab Properties /////////////////////////////////
    ;-------------------------------------------------------
    endif

    ;Destroy the object.
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim

;-------------------------------------------------------
;Prepare output ////////////////////////////////////////
;-------------------------------------------------------
    ;Get all of the plots
    drWin['Cut ' + names[0]] -> SetProperty, XTITLE='', XTICKFORMAT='(a1)', TITLE=title
    drWin['Cut ' + names[1]] -> SetProperty, XTITLE='', XTICKFORMAT='(a1)', TITLE=''
    drWin['Cut ' + names[2]] -> SetProperty, XTITLE='', XTICKFORMAT='(a1)', TITLE=''
    drWin['Cut ' + names[3]] -> SetProperty, XTITLE='', XTICKFORMAT='(a1)', TITLE=''
    drWin['Cut ' + names[4]] -> SetProperty, XTITLE=axLabls[0] + ' (' + units + ')', TITLE=''
    
    ;Change the character size
    drWin -> SetGlobal, CHARSIZE=1.8
    
    ;Refresh and output, if requested.
    drWin -> Refresh
    if ofilename ne '' then drWin -> Save, ofilename
    
    return, drWin
end
