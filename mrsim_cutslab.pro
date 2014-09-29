; docformat = 'rst'
;
; NAME:
;    MrSim_CutSlab
;
; PURPOSE:
;+
;   The purpose of this program is to make two columns of plots. The right column will
;   be color plots of all quantities specified by `NAMES` and the left column will be
;   line plots of those quantities taken at `CUT`.::
;       NAME                    DESCRIPTION
;       "Color: [name]"     -   2D image of the quantity specified by "name"
;       "Cut: [name] x=#.#" -   1D vertical cut through the 2D image. #.# locates the cut.
;       "[name] Contours    -   Contours overlayed on the image.
;       "CB: Color [name]   -   Colorbar associated with the 2D image.
;       "VLine [name] #     -   Vertical line overplotted on the image. # = cut number
;       "HLines [name]      -   Horiztonal line overplotted on the image.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       NAMES:              in, required, type=string
;                           The name of the quantitis to be plotted in color.
;       CUT:                in, required, type=int
;                           The x-location (in de) at which the vertical cut will
;                               be taken.
;       TIME:               in, required, type=int
;                           The simulation time for which at which to create the plot.
;                               Igored if `SIM_OBJECT` is given.
;
; :Keywords:
;       C_NAME:             in, optional, type=string, default=''
;                           Name of the data product whose contours are to be 
;                               overplotted on the image.
;       FRACTION:           in, optional, type=float
;                           Fraction of the vectors to be drawn. Used with `VX_NAME`
;                               and `VY_NAME`.
;       HCUT_RANGE:         in, optional, type=fltarr(2)
;                           The horizontal range, in data coordinates, over which to cut.
;                               "Horizontal" is defined by the `HORIZONTAL` keyword. Data
;                               coordinates are determined by the "ion_scale" property
;                               within `SIM_OBJECT`. The default is to take the
;                               appropriate simulation range.
;       HORIZONTAL:         in, optional, type=boolean, default=0
;                           If set, a horizontal cut will be taken. The default is
;                               to take a vertical cut. For an "XY" orientation, "X" is
;                               horizontal, "Y" is vertical. Similar for "XZ", etc.
;                               The orientationis taken from the `SIM_OBJECT` property.
;       NLEVELS:            in, optional, type=int, default=15
;                           The number of contour lines to draw
;       OFILENAME:          in, optional, type=string, default=''
;                           If provided, graphics will be output to a file by this
;                               name. The file extension is used to determine which
;                               type of image to create. "PNG" is the default.
;       SIM_OBJECT:         out, optional, type=object
;                           The "Sim_Object" (or subclass) reference containing the 
;                               simulation data and domain information used in
;                               creating the plot. If not provided, one will be
;                               created according to the `SIM3D` keyword.
;       SIM3D:              in, optional, type=boolean, default=0
;                           If set, and `SIM_OBJECT` is not given, then a 3D
;                               simulation object will be created. The default is
;                               to make 2D simulation object.
;       VCUT_RANGE:         in, optional, type=fltarr(2)
;                           The verical range, in data coordinates, over which to cut.
;                               "Vertical" is defined by the `HORIZONTAL` keyword. Data
;                               coordinates are determined by the "ion_scale" property
;                               within `SIM_OBJECT`. The default is to take the
;                               appropriate simulation range.
;       VX_NAME:            in, optional, type=string, default=''
;                           Name of the x-component of a data product whose vector
;                               field is to be overplotted on the image. Must be used
;                               with `VY_NAME`.
;       VY_NAME:            in, optional, type=string, default=''
;                           Name of the y-component of a data product whose vector
;                               field is to be overplotted on the image. Must be used
;                               with `VX_NAME`.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted Sim_Object::Init(), or any of its
;                               subclasses is also accepted here. Ignored if
;                               `SIM_OBJECT` is present.
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       MrSim_LineCut.pro
;       MrSim_ColorSlab.pro
;       MrSim2D__Define.pro
;       MrSim3D__Define.pro
;       MrWindow.pro
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
;       2014/02/07  -   Written by Matthew Argall
;       2014/03/05  -   Parameter CUT renamed to CUTS and can now be an array. - MRA
;       2014/03/27  -   Added the VX_NAME and VY_NAME keywords. - MRA
;-
function MrSim_CutSlab, names, cuts, time, $
C_NAME = c_name, $
FRACTION = fraction, $
HCUT_RANGE = hcut_range, $
HORIZONTAL = horizontal, $
NLEVELS = nlevels, $
OFILENAME = ofilename, $
SIM_OBJECT = oSim, $
SIM3D = Sim3D, $
VCUT_RANGE = vcut_range, $
VX_NAME = vx_name, $
VY_NAME = vy_name, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(drWin) then obj_destroy, drcWin
        void = cgErrorMSG()
        return, obj_new()
    endif

;-------------------------------------------------------
;Defaults //////////////////////////////////////////////
;-------------------------------------------------------
    ;Set defaults
    horizontal = keyword_set(horizontal)
    sim_class  = keyword_set(sim3d) ? 'MrSim3D' : 'MrSim2D'
    if n_elements(ofilename) eq 0  then ofilename = ''
    if ofilename             eq '' then buffer    = 0 else buffer = 1

    ;Create a simulation object
    if n_elements(oSim) eq 0 then begin
        oSim = obj_new(sim_class, time, _STRICT_EXTRA=extra)
        if obj_valid(oSim) eq 0 then return, obj_new()
    endif else begin
        if obj_isa(oSim, 'MRSIM') eq 0 then $
            message, 'SIM_OBJECT must be a subclass of "MrSim"'
    endelse

    nCuts  = n_elements(cuts)
    nNames = n_elements(names)

;-------------------------------------------------------
;Graphic Annotations ///////////////////////////////////
;-------------------------------------------------------
    oSim -> GetProperty, ORIENTATION=orientation, TIME=time, COORD_SYSTEM=coord_system
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    
    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then t_title = 't$\Omega$$\downc,i$$\up-1$=' + string(time*dtxwci, FORMAT='(f0.1)') $
        else t_title = 't$\downindex$=' + string(time, FORMAT='(i0)')
    
    ;Labels, etc.
    if horizontal eq 1 then begin
        xaxis = strmid(strlowcase(orientation), 0, 1)
        cutaxis = strmid(strlowcase(orientation), 1, 1)
        horiz_lines = cuts
    endif else begin
        xaxis = strmid(strlowcase(orientation), 1, 1)
        cutaxis = strmid(strlowcase(orientation), 0, 1)
        vert_lines = cuts
    endelse

    ;Title for where cut was taken.
    if nCuts gt 1 $
        then d_title = cutaxis + '= [' + strjoin(string(cuts, FORMAT='(f0.1)'), ', ') + ']' $
        else d_title = cutaxis + '=' + string(cuts, FORMAT='(f0.1)') + units
    if obj_class(oSim) eq 'MRSIM3D' $
        then d_title += '  y=' + string(oSim.yslice, FORMAT='(f0.1)') + units

;-------------------------------------------------------
;Create Graphics ///////////////////////////////////////
;-------------------------------------------------------
    ;Set the layout.
    ;   Magnetopause plots have color plots on top and line plots on the bottom.
    ;   Other orienatations have color plots on the left and line plots on the right.
    if coord_system eq 'MAGNETOPAUSE' $
        then layout = [nNames,2] $
        else layout = [2,nNames]

    ;Make the window
    drcWin = MrWindow(XSIZE=750, YSIZE=400, BUFFER=buffer, REFRESH=0, YGAP=1, XGAP=8, $
                      OXMARGIN=[10,12], LAYOUT=layout)
                      
    ;Make the graphics
    iGr = intarr(nNames)
    count = 0
    for i = 0, nNames - 1 do begin
        
        ;Make a line plot of the data
        lcWin = MrSim_LineCut(names[i], cuts, /CURRENT, $
                              HCUT_RANGE=hcut_range, $
                              HORIZONTAL=horizontal, $
                              SIM_OBJECT=oSim, $
                              VCUT_RANGE=vcut_range)

        ;Skip to the next iteration if the graphic could not be made
        if obj_valid(lcWin) eq 0 then continue
        iGr[count] = i
        count += 1

        ;Create an overview image for DATA.
        !Null = MrSim_ColorSlab(names[i], /CURRENT, $
                                C_NAME = c_name, $
                                FRACTION = fraction, $
                                HORIZ_LINES = horiz_lines, $
                                NLEVELS = nlevels, $
                                SIM_OBJECT = oSim, $
                                VERT_LINES = vert_lines, $
                                VX_NAME = vx_name, $
                                VY_NAME = vy_name)

        ;Move the colorbar to the top
        if coord_system eq 'MAGNETOPAUSE' $
            then drcWin['CB: Color ' + names[i]] -> SetProperty, CBLOCATION='BOTTOM', OFFSET=4, WIDTH=1
    endfor

    ;Destroy the object.
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
    ;Get rid of the holes.
    iGr = iGr[0:count-1]
    if count ne nNames then drcWin -> FillHoles, /TRIMLAYOUT

;-------------------------------------------------------
;Prep for Output ///////////////////////////////////////
;-------------------------------------------------------

    ;Turn off some properties.
    graphics = drcWin -> Get(/ALL, ISA=['MRPLOT', 'MRIMAGE'])
    foreach gfx, graphics do gfx -> SetProperty, XTICKFORMAT='(a1)', XTITLE='', TITLE=''

    ;Put the X-Axis back on the bottom graphics. Put titles back on the top graphics
    drcWin['Cut '   + names[iGr[0]]].TITLE=d_title
    drcWin['Color ' + names[iGr[0]]].TITLE=t_title
    drcWin['Cut '   + names[iGr[-1]]] -> SetProperty, XTICKFORMAT='', XTITLE=xaxis   + ' (' + units + ')'
    drcWin['Color ' + names[iGr[-1]]] -> SetProperty, XTICKFORMAT='', XTITLE=cutaxis + ' (' + units + ')'

    ;Refresh and output, if requested.
    if coord_system eq 'MAGNETOPAUSE' $
        then drcWin -> SetProperty, XSIZE=((400 > nNames*125) < 1000), YSIZE=500, OYMARGIN=[10,3] $
        else drcWin.YSIZE = nNames*125 < 690
        

    drcWin -> Refresh
    if ofilename ne '' then begin
        drcWin -> Save, ofilename
        obj_destroy, drcWin
        drcWin = obj_new()
    endif
    
    return, drcWin
end
