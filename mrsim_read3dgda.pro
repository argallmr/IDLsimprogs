; docformat = 'rst'
;
;+
;   Read 3D .gda field and moment files produced by VPIC at LANL.
;
;*****************************************************************************************
;   Copyright (c) 2014, Matthew Argall                                                   ;
;   All rights reserved.                                                                 ;
;                                                                                        ;
;   Redistribution and use in source and binary forms, with or without modification,     ;
;   are permitted provided that the following conditions are met:                        ;
;                                                                                        ;
;       * Redistributions of source code must retain the above copyright notice,         ;
;         this list of conditions and the following disclaimer.                          ;
;       * Redistributions in binary form must reproduce the above copyright notice,      ;
;         this list of conditions and the following disclaimer in the documentation      ;
;         and/or other materials provided with the distribution.                         ;
;       * Neither the name of the University of New Hampshire nor the names of its       ;
;         contributors may be used to endorse or promote products derived from this      ;
;         software without specific prior written permission.                            ;
;                                                                                        ;
;   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY  ;
;   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES ;
;   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT  ;
;   SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,       ;
;   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED ;
;   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR   ;
;   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     ;
;   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN   ;
;   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH  ;
;   DAMAGE.                                                                              ;
;*****************************************************************************************
;
; :Categories:
;   Pic Simulation, Data Reader, Bill Daughton
;
; :Author:
;   Matthew Argall::
;       University of New Hampshire
;       Morse Hall, Room 113
;       8 College Rd.
;       Durham, NH, 03824
;       matthew.argall@wildcats.unh.edu
;
; :Copyright:
;       Copyright 2014 by the University of New Hampshire
;
; :History:
;   Modification History::
;       2014/10/11  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Read 3D gda files.
;
; :Params:
;       THESIM:         in, required,  type=string/integer/object
;                       The name, number, or object of the simulation for which the data
;                           is to be read. See MrSim_Which and MrSim_Create. If an object
;                           is given, all defaults are taken from the properties of the
;                           object.
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;       TINDEX:         in, optional, type=integer, default=0
;                       Time-index at which data is to be read. For 3D simulations, this
;                           is the indicated in the file name.
;
; :Keywords:
;       COORD_SYSTEM:   in, optional, type=string, default='SIMULATION'
;                       Coordinate system in which data is to be returned. Choices are::
;                           'SIMULATION'
;                           'MAGNETOPAUSE'
;                           'MAGNETOTAIL'
;       DIRECTORY:      in, optional, type=string
;                       Directory in which to find .gda data. By default, the directory
;                           indicated by MrSim_Which is used.
;       NSMOOTH:        in, optional, type=integer, default=0
;                       Number of points over which the data should be smoothed.
;       ORIENTATION:    in, optional, type=string, default='XZ'
;                       Orientation of the 2D plane to be read. Options are::
;                           'XY'
;                           'XZ'
;                           'YZ'
;       SIM_OBJECT:     out, optional, type=object
;                       If `THESIM` is not an object, one will be created Set SIM_OBJECT
;                           equal to a named variable into which the simulation object
;                           may be returned.
;       XRANGE:         in, optional, type=dblarr(2), default=entire x-range
;                       X-range over which data should be read. If `YZ_PLANE` is set,
;                           XRANGE[0] indicates the location in X at which the 2D slice
;                           will be taken. By default, the entire simulation domain in 
;                           the x-dimension is read.
;       YRANGE:         in, optional, type=dblarr(2), default=entire y-range
;                       Y-range over which data should be read. If `XZ_PLANE` is set,
;                           YRANGE[0] indicates the location in Y at which the 2D slice
;                           will be taken. By default, the entire simulation domain in 
;                           the y-dimension is read.
;       ZRANGE:         in, optional, type=dblarr(2), default=entire z-range
;                       Z-range over which data should be read. If `XY_PLANE` is set,
;                           ZRANGE[0] indicates the location in Z at which the 2D slice
;                           will be taken. By default, the entire simulation domain in 
;                           the z-dimension is read.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrSim_Create is also accepted. Ingored
;                           if `THESIM` is an object.
;-
function MrSim_Read3Dgda, theSim, name, tIndex, $
COORD_SYSTEM=coord_system, $
DIRECTORY=directory, $
NSMOOTH=nSmooth, $
ORIENTATION=orientation, $
SIM_OBJECT=oSim, $
XRANGE=xrange, $
YRANGE=yrange, $
ZRANGE=zrange, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
        void = cgErrorMsg()
        return, !Null
    endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
    osim_created = 0B
    
    ;Simulation name or number?
    if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
        oSim = MrSim_Create(theSim)
        if obj_valid(oSim) eq 0 then return, obj_new()
        osim_created = 1B
        
    ;Object?
    endif else if MrIsA(theSim, 'OBJREF') then begin
        if obj_isa(theSim, 'MRSIM') eq 0 $
            then message, 'THESIM must be a subclass of the MrSim class.' $
            else oSim = theSim
            
    ;Unknown
    endif else begin
        MrSim_Which
        message, 'THESIM must be a simulation name, number, or object.'
    endelse
    sim_class = obj_class(oSim)

;-------------------------------------------------------
; Read the Data ////////////////////////////////////////
;-------------------------------------------------------
    data = oSim -> ReadMoment( name, tIndex, $
                               COORD_SYSTEM = coord_system, $
                               DIRECTORY    = directory, $
                               NSMOOTH      = nSmooth, $
                               ORIENTATION  = orientation, $
                               XRANGE       = xrange, $
                               YRANGE       = yrange, $
                               ZRANGE       = zrange $
                             )

    ;Destroy the simulation object if it is not being returned
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    
    return, data
end
