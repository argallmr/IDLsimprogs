; docformat = 'rst'
;
; NAME:
;    MrSim_eFlyBy
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
;       2015/01/23  -   Written by Matthew Argall
;-
function MrSim_eFlyBy, theSim, type, bin_center, half_width, time, $
C_NAME=c_name, $
CIRCLES=circles, $
HGAP=hGap, $
IM_NAME=im_name, $
NBINS=nBins, $
LAYOUT=layout, $
LOCATION=location, $
POSITIONS=positions, $
SIM_OBJECT=oSim, $
VGAP=vGap, $
V_VA=v_va, $
XSIZE=xsize, $
YSIZE=ysize, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim)   then obj_destroy, oSim
        if obj_valid(im_win) then obj_destroy, im_win
        if obj_valid(fbWin)  then obj_destroy, fbWin
        void = cgErrorMSG()
        return, obj_new()
    endif

    ;Set defaults
    if n_elements(ofilename) eq 0 then ofilename = ''
    
;-------------------------------------------------------
; Create the eMap //////////////////////////////////////
;-------------------------------------------------------

    ;Create the distribution functions and the image.    
    fbWin = MrSim_eMap(theSim, type, bin_center, half_width, time, $
                       C_NAME     = c_name, $
                       CIRCLES    = circles, $
                       HGAP       = hGap, $
                       IM_NAME    = im_name, $
                       IM_WIN     = im_win, $
                       NBINS      = nBins, $
                       LAYOUT     = layout, $
                       LOCATION   = location, $
                       POSITIONS  = positions, $
                       SIM_OBJECT = oSim, $
                       VGAP       = vGap, $
                       V_VA       = v_va, $
                       XSIZE      = xsize, $
                       YSIZE      = ysize, $
                       _EXTRA     = extra)

    ;Disable refresh
    fbWin  -> Refresh, /DISABLE
    im_win -> Refresh, /DISABLE
    
;-------------------------------------------------------
; Move the Image into the Window ///////////////////////
;-------------------------------------------------------
    
    ;Was an image created?
    if obj_valid(im_win) then begin
        ;Get the image
        gfx = im_win -> Get(ISA='MrImage')
        
        ;Make the eMap window the current window
        fbWin -> SetCurrent
       
        ;Make room for the image
        oSim -> GetProperty, COORD_SYSTEM=coord_sys
        if strupcase(coord_sys) eq 'MAGNETOPAUSE' then begin
            fbWin.OXMARGIN = [20, 15]
            gfx.POSITION   = [0.15, 0.1, 0.4, 0.77]
        endif else begin
            fbWin.OYMARGIN = [4, 18]
            gfx.POSITION   = [0.2, 0.7, 0.8, 0.95]
        endelse

        ;Move everything into the eMap window
        gfx = im_win -> Get(/ALL)
        foreach graph, gfx do graph -> SwitchWindows, fbWin
        obj_destroy, im_win
    endif
    
    ;Adjust the color bar and title position
    fbWin['CB: Counts'].position = [0.92, 0.1, 0.94, 0.5]
    fbWin['eMap Title'] -> SetProperty, YLOC=0.575
    
    ;Destroy the simulation
    obj_destroy, oSim
    
    ;Refresh and output, if requested.
    fbWin -> Refresh
    if ofilename ne '' then drWin -> Save, ofilename
    
    return, fbWin
end
