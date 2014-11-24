; docformat = 'rst'
;
; NAME:
;    SimScript_CutSlab
;
; PURPOSE:
;+
;   Create a set of standard plots at different locations across the reconnection exhaust.
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
pro SimScript_Standard
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim) then obj_destroy, oSim
        if obj_valid(win)  then obj_destroy, win
        void = cgErrorMSG()
        return
    endif

    ;Inputs
    directory  = '/home/argall/Work/AssymectricSim/asymm-2D-large/'
    data_dir   = directory + 'data/'
    ascii_info = directory + 'info'
    output_dir = '/home/argall/Desktop/paper_xline/'
    cuts       = [150.9, 150.3, 148.4]
    time       = 32
    Sim3D      = 0
    c_name     = 'By'
    ion_scale  = 1
    mva_frame  = 1
    coord_sys  = 'Magnetopause'     ;{'Magnetopause' | 'Magnetotail' | 'Simulation'}
    xrange     = [4, -3]
    zrange     = 150.9+[-3.0,3.0]
    yslice     = 1108
    
;---------------------------------------------------------------------
;DO NOT EDIT BELOW HERE!! ////////////////////////////////////////////
;---------------------------------------------------------------------

    ;Create the simulation object
    class = Sim3D eq 1 ? 'MrSim3D' : 'MrSim2D'
    oSim = obj_new(class, time, yslice, ION_SCALE=ion_scale, $
                   MVA_FRAME=mva_frame, COORD_SYSTEM=coord_sys, $
                   XRANGE=xrange, ZRANGE=zrange, $
                   DIRECTORY=data_dir, INFO_FILE=ascii_info)
        
    ;Horizontal or vertical lines?
    case strupcase(coord_sys) of
        'MAGNETOPAUSE': hLines = cuts
        'MAGNETOTAIL':  vLines = cuts
        'SIMULATION':   vLines = cuts
        else: ;Do nothing
    endcase
    
    ;Get the mass ratio
    ;   - Build filename in units of de to avoid decimal places.
    if ion_scale $
        then oSim -> GetInfo, MI_ME=mi_me $
        else mi_me = 1.0
    mi_me = sqrt(mi_me)
    
    ;Output Filename
    if Sim3D $
        then fileBase = 'sim3D_t' + strtrim(time, 2) + '_y' + strtrim(yslice, 2) + '_' $
        else fileBase = 'sim2D_t' + strtrim(time, 2) + '_'
    fileTail = string(FORMAT='(%"_%ix%i_%iz%i")', abs(xrange*mi_me), abs(zrange*mi_me))

    ;Create and save all of the plots
    for i = 0, n_elements(cuts) - 1 do begin
        fileCut = '_c' + strtrim(fix(cuts[i]*mi_me), 2) + '.png'
        
        ;Color plot with lines
        oFile = filepath(fileBase + 'ColorSlab' + fileTail + fileCut, ROOT_DIR=output_dir)
        win   = MrSim_ColorSlab(c_name, SIM_OBJ=oSim, HORIZ_LINES=hLines, VERT_LINES=vLines, $
                                OFILENAME=oFile)
        obj_destroy, win
        
        ;Diffusion Indicators
        oFile = filepath(fileBase + 'DiffRegion' + fileTail + fileCut, ROOT_DIR=output_dir)
        win   = MrSim_DiffRegion(cuts[i], SIM_OBJ=oSim, OFILENAME=oFile)
        obj_destroy, win
        
        ;X-Line Proximity Indicators
        oFile = filepath(fileBase + 'XProx' + fileTail + fileCut, ROOT_DIR=output_dir)
        win   = MrSim_XProximity(cuts[i], SIM_OBJ=oSim, OFILENAME=oFile)
        obj_destroy, win
    endfor
    
    ;Destroy the simulation object.
    obj_destroy, oSim
end
