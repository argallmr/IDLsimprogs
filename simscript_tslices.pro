; docformat = 'rst'
;
; NAME:
;    SimScript_CutSlab
;
; PURPOSE:
;+
;   The purpose of this program is to serve as a script for creating and saving
;   graphics from MrSim_MultiSlab.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       MrSim_CutSlab.pro
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
;       2014/02/19  -   Written by Matthew Argall
;-
pro simScript_tSlices
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim) then obj_destroy, oSim
        if obj_valid(csWin) then obj_destroy, csWin
        void = cgErrorMSG()
        return
    endif

    ;Inputs
    simname  = 'Asymm-Large-2D'
    output_dir = '/home/argall/figures/tslices/tslices_Asym-Large-2D_nongyro-fixed/'
    delta_t    = 2
    names      = ['Dng_e', 'ne', 'Jey', 'An_e', 'A0_e']
    c_name     = 'Ay'
    ion_scale  = 0
    xrange     = [1300,1900]
    zrange     = [-80,50]
    
    ranges = [[0.00, 0.16], $
              [0.00, 3.00], $
              [0.00, 0.16], $
              [0.00, 5.00], $
              [0.00, 0.42]]

;---------------------------------------------------------------------
;DO NOT EDIT BELOW HERE!! ////////////////////////////////////////////
;---------------------------------------------------------------------

    ;Create the simulation object
    oSim = MrSim_Create(simname, time, ION_SCALE=ion_scale, $
                        XRANGE=xrange, ZRANGE=zrange)
    
    ;How big is a "de"?
    nt = 65
    oSim -> GetInfo, MI_ME=mi_me, UNITS=units
    if ion_scale then delta = 1/sqrt(mi_me) else delta = 1
    
    ;Create a base for the output file name.
    oFileBase = output_dir + 'sim2D_Asymm-Large_tslices_nongyro-fixed_t'

    ;Create and save all of the plots
    for i = 0, nt - 1, delta_t do begin
        oSim.time = i
        oFileName = oFileBase + string(i, FORMAT='(i03)') + '.png'
        csWin = MrSim_MultiSlab(names, C_NAME=c_name, RANGES=ranges, $
                                SIM_OBJECT=oSim, OFILENAME=oFileName)
    endfor
    
    ;Destroy the simulation object.
    obj_destroy, oSim
end
