; docformat = 'rst'
;
; NAME:
;       MrSim_tParticle
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
;       * Neither the name of the <ORGANIZATION> nor the names of its contributors may   ;
;         be used to endorse or promote products derived from this software without      ;
;         specific prior written permission.                                             ;
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
; PURPOSE:
;+
;   Create input files for the "tparticle" simulation.
;
; :History:
;   Modification History::
;       2014-11-11  -   Written by Matthew Argall
;-
pro MrSim_tParticle
    compile_opt idl2
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        if obj_valid(oSim) then obj_destroy, oSim
        void = cgErrorMsg()
        return
    endif
    
    ;Files to which information will be written
    save_dir  = '~/programs/tracer/inputs/'
    tfile     = filepath('timestep.txt',       ROOT_DIR=save_dir)
    ofile     = filepath('ofname.txt',         ROOT_DIR=save_dir)
    nfile     = filepath('n_electrons.txt',    ROOT_DIR=save_dir)
    ifile     = filepath('initialization.txt', ROOT_DIR=save_dir)
    dfile     = filepath('direct.txt',         ROOT_DIR=save_dir)
    nElec     = 1
    direction = 3
    ofile_suffix = ''
    
    ;Simulation object
    theSim = 'Asymm-3D'
    tIndex = 108090
    xrange = 457.0 + [-40, 40]
    yrange = [102.2, 102.2]
    zrange = [-15, 15]
    oSim   = MrSim_Create(theSim, time, XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    
    ;Read particles
    x        = 457.0
    y        = 102.25
    z        = 1.86
    dx       = 1.0
    dy       = 1.0
    dz       = 1.0
    vx_range = [-0.1, 0.1]
    vy_range = [-0.1, 0.1]
    vz_range = [-0.1, 0.1]

    ;Convert y-location to grid cell
    ycell = oSim -> GetCell(y, /Y)
    MrSim_Which, theSim, EFILE=eFile, DIST3D=dist3D, YSLICE=ycell, TINDEX=tIndex

    ;Read the data
    data = oSim -> ReadElectrons(eFile, $
                                 DIST3D=dist3D, $
                                 /VELOCITY, $
                                 XRANGE=x + [-dx, dx], $
                                 YRANGE=y + [-dy, dy], $
                                 ZRANGE=z + [-dz, dz], $
                                 VX_RANGE=vx_range, $
                                 VY_RANGE=vy_range, $
                                 VZ_RANGE=vz_range)
    obj_destroy, oSim
    
    ;Ensure a y-position is defined
    dims = size(data, /DIMENSIONS)
    if dims[0] eq 5 then data = [data[0,*], fltarr(1, dims[2]), data[1:4,*]]

    ;Reduce the number of electrons
    if nElec lt dims[1] then begin
        if nElec eq 1 $
            then iElec = 0 $
            else iElec = long(findgen(nElec)/(nElec-1) * dims[1])
        data  = data[*,iElec]
    endif else if nElec gt dims[1] then begin
        nElec = dims[1]
    endif

    ;Number of electrons found.
    message, '', /INFORMATIONAL
    print, FORMAT='(%"  # Electrons Found: %i")', dims[1]
    print, FORMAT='(%"  # Electrons Used:  %i")', nElec
    
    ;Time Index
    openw,    lun, tfile, /GET_LUN
    printf,   lun, tIndex, FORMAT='(i0)'
    free_lun, lun
    
    ;Electron counts
    openw,    lun, nfile, /GET_LUN
    printf,   lun, nElec, FORMAT='(i0)'
    free_lun, lun
    
    ;Direction
    openw,    lun, dfile, /GET_LUN
    printf,   lun, direction, FORMAT='(i0)'
    free_lun, lun
    
    ;Output filename
    ofname = string(FORMAT='(%"%s_t%i_n%i")', theSim, tIndex, nElec) + ofile_suffix
    openw,    lun, ofile, /GET_LUN
    printf,   lun, ofname, FORMAT='(a0)'
    free_lun, lun
    
    ;Particle data
    openw,    lun, ifile, /GET_LUN
    printf,   lun, data, FORMAT='( (3(f10.4, 3x), 2(f8.4, 3x), f8.4) )'
    free_lun, lun
    
    ;Indicate where files were saved.
    print, FORMAT='(%"  Files written to:  %s")', save_dir
end