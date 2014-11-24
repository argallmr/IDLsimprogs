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
;       2014/02/12  -   Written by Matthew Argall
;-
pro simScript_CutSlab
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
    directory = '/home/argall/Work/AssymectricSim/asymm-3D/'
    data_dir = directory + 'data/'
    ascii_info = directory + 'info'
    output_dir = '/home/argall/Desktop/3DPlots/t103286_y1108/'
    names = ['Bx', 'By', 'ni', 'Uix', 'Ez', 'Jey']
    time = 103286
    c_name = ''
    ion_scale = 1
    yslice = 1108
    Sim3D = 1
    xrange = [35, 65]
    zrange = [-3, 4]
    
;    directory = '/home/argall/Work/AssymectricSim/asymm-2D-large/'
;    data_dir = directory + 'data/'
;    ascii_info = directory + 'info'
;    output_dir = '/home/argall/Desktop/2DPlots/t90/'
;    names = ['Bx', 'By', 'ni', 'Uix', 'Ez', 'Jey']
;    time = 90
;    c_name = 'Ay'
;    ion_scale = 1
;    yslice = 0
;    Sim3D = 0
;    xrange = [131, 270]
;    zrange = [-13, 9]

;---------------------------------------------------------------------
;DO NOT EDIT BELOW HERE!! ////////////////////////////////////////////
;---------------------------------------------------------------------

    ;Create the simulation object
    class = Sim3D eq 1 ? 'MrSim3D' : 'MrSim2D'
    oSim = obj_new(class, time, yslice, ION_SCALE=ion_scale, $
                   XRANGE=xrange, ZRANGE=zrange, $
                   DIRECTORY=data_dir, INFO_FILE=ascii_info)
    
    ;How big is a "de"?
    oSim -> GetInfo, MI_ME=mi_me, UNITS=units
    if ion_scale then delta = 1/sqrt(mi_me) else delta = 1

    ;Where will the cuts be taken?
    vCuts = linspace(xrange[0], xrange[1]-.01, 2*delta, /INTERVAL)
    
    ;Create a base for the output file name.
    slice = Sim3D ? '_y' + string(yslice, FORMAT='(i04)') : ''
    oFileBase = output_dir + 'sim_cutslab_t' + string(time, FORMAT='(i0)') + '_' + $
                units + slice + '_x'

    ;Create and save all of the plots
    for i = 0, n_elements(vCuts) - 1 do begin
        x = ion_scale ? vCuts[i] * sqrt(mi_me) : vCuts[i]
        oFileName = oFileBase + string(x, FORMAT='(i0)') + '.png'
    
        csWin = MrSim_CutSlab(names, vCuts[i], SIM3D=Sim3D, SIM_OBJECT=oSim, $
                                               C_NAME=c_name, OFILENAME=oFileName)
    endfor
    
    ;Destroy the simulation object.
    obj_destroy, oSim
end
