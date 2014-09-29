; docformat = 'rst'
;
; NAME:
;    VECTOR_PLOTS
;
; PURPOSE:
;+
;   The purpose of this program is to create vector field from physical quantities.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       TIME:                   in, required, type=int
;                               The simulation time for which a velocity vector field
;                                   is to be plotted.
;       NAME:                   in, optional, type=string, default='By'
;                               The name of the quantity to be plotted in color underneath
;                                   the velocity vector field. Type "IDL> void = mr_readsim()"
;                                   for a list of all available names.
; :Keywords:
;       DIRECTORY:              in, optional, type=string, default=pwd
;                               The directory in which to find all of the simulation data
;       FRACTION:               in, optional, type=float, default=0.02
;                               The fraction of the total number of data points to be
;                                   turned into velocity vectors.
;       LENGTH:                 in, optional, type=float, default=0.08
;                               The maximum vectorl ength relative to the plot data window
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
;       2014/03/28  -   Updated to use newer library routines and renamed from
;                           vector_plots to MrSim_VectorField. - MRA
;-
function MrSim_VectorFields, vx_name, vy_name, time, $
CURRENT = current, $
FRACTION = fraction, $
LENGTH = length, $
OFILENAME = ofilename, $
ORDERED = ordered, $
OVERPLOT = overplot, $
SIM3D = sim3D, $
SIM_OBJECT = oSim, $
_REF_EXTRA = extra
    compile_opt idl2
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim) && arg_present(oSim) eq 0 then obj_destroy, oSim
        void = cgErrorMSG()
        return, !Null
    endif

;-------------------------------------------------------
;Check Inputs //////////////////////////////////////////
;-------------------------------------------------------

    ;Fraction of vectors to be plotted
    ordered = n_elements(ordered) eq 0 ? 1 : keyword_set(ordered)
    if n_elements(fraction) eq 0 then fraction = 0.005
    if n_elements(oFilename) eq 0 then oFilename = ''
    
    ;Create a simulation object
    if n_elements(oSim) gt 0 then begin
        if obj_valid(oSim) eq 0 || obj_isa(oSim, 'MRSIM') eq 0 then $
            message, 'SIM_OBJECT must be valid and a subclass of "MrSim"'
        sim_class = obj_class(oSim)
    endif else begin
        Sim3D   = keyword_set(Sim3D)
        if Sim3D then sim_class = 'MRSIM3D' else sim_class = 'MRSIM2D'
        oSim = obj_new(sim_class, time, _STRICT_EXTRA=extra)
    endelse
    
;-------------------------------------------------------
;Read Data /////////////////////////////////////////////
;-------------------------------------------------------
    ;Read in the velocities
    vecx = oSim -> GetData(vx_name)
    vecy = oSim -> GetData(vy_name)

    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, XSIM=XSim, ZSIM=ZSim
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
        
    ;Destroy the object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
    ;Make sure the data exists. Do this /after/ OSIM has a chance to be destroyed
    ;and before the graphic window is created.
    if vecx eq !Null || vecy eq !Null then return, obj_new()

;-------------------------------------------------------
;Prepare Ranges ////////////////////////////////////////
;-------------------------------------------------------
    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't=' + string(time*dtxwci, FORMAT='(f0.1)') + ' $\Omega$$\downc,i$$\up-1$' $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')

    if sim_class eq 'MRSIM3D' then begin
        oSim -> GetProperty, YSLICE=yslice
        title += '  y=' + string(yslice, FORMAT='(f0.1)') + units
    endif

;-------------------------------------------------------
;Overview //////////////////////////////////////////////
;-------------------------------------------------------

    ;Create a Zoomed-In Plot of By
    vField = MrVector(vecx, vecy, XSim, ZSim, $
                      CURRENT=current, OVERPLOT=overplot, $
                      NAME=vx_name + '-' + vy_name, TITLE=title, $
                      FRACTION=fraction, LENGTH=length, ORDERED=ordered, $
                      XTITLE='x (' + units + ')', YTITLE='z + (' + units + ')')

;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------    
    ;Refresh and output, if requested.
    if current eq 0 then begin
        vField -> Refresh
        if ofilename ne '' then colorWin -> Save, ofilename
    endif

end
