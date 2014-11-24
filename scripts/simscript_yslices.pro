; docformat = 'rst'
;
; NAME:
;    SimScript_CutSlab
;
; PURPOSE:
;+
;   The purpose of this program is to serve as a script for creating and saving
;   graphics from MrSim_CutSlab.
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
;       2014/02/16  -   Written by Matthew Argall
;-
pro simScript_ySlices
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim)  then obj_destroy, oSim
        if obj_valid(csWin) then obj_destroy, csWin
        void = cgErrorMSG()
        return
    endif

    ;Inputs
;    directory  = '/home/argall/Work/AssymectricSim/asymm-3D/'
    directory  = '/home/daughton/Asymm-3D/'
    data_dir   = directory + 'data/'
    ascii_info = directory + 'info'
    output_dir = '/home/argall/Desktop/3DPlots/yslices_t108090/'
    time       = 108090
    delta_y    = 5
    names      = 'Jey'
    ion_scale  = 1
    xrange     = 44+[-20, 25]
    zrange     = [-8,8]

;---------------------------------------------------------------------
;DO NOT EDIT BELOW HERE!! ////////////////////////////////////////////
;---------------------------------------------------------------------

    ;Create the simulation object
    oSim = obj_new('MrSim3D', time, 0, ION_SCALE=ion_scale, $
                   XRANGE=xrange, ZRANGE=zrange, $
                   DIRECTORY=data_dir, INFO_FILE=ascii_info)
    
    ;How big is a "de"?
    oSim -> GetInfo, MI_ME=mi_me, NY=ny, UNITS=units
    if ion_scale then delta = 1/sqrt(mi_me) else delta = 1
    
    ;Create a base for the output file name.
    oFileBase = output_dir + 'sim3D_yslices_t' + string(time, FORMAT='(i06)') + '_y'

    ;Create and save all of the plots
    for i = 0, ny - 1, delta_y do begin
        ;Set the y-slice and create a file name
        oSim.yslice = i
        oFileName = oFileBase + string(i, FORMAT='(i04)') + '.png'
        
        ;Create the plot
        csWin = MrSim_MultiSlab(names, SIM_OBJECT=oSim, OFILENAME=oFileName)
        obj_destroy, csWin
    endfor
    
    ;Destroy the simulation object.
    obj_destroy, oSim
end
