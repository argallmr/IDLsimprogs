; docformat = 'rst'
;
; NAME:
;    MrSim_Which
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
; PURPOSE:
;+
;   The purpose of this program is to provide information about the available simulations.
;   If information is available for a particular simulation parameter, !Null is returned.
;
; :Categories:
;    Bill Daughton, Simulation
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
;       2014/09/12  -   Written by Matthew Argall
;       2014/10/01  -   Replaced property-specific helper functions with simulation-
;                           specific helper functions. - MRA
;       2014/10/02  -   Forgot the EREGIONS keyword. - MRA
;-
;*****************************************************************************************
;+
;   Return the electron count factor.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time shown in `EFILE`
;                                 YSLICE: intarr(N), $          y-slice in the .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_mime1836_by00, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
compile_opt strictarr
on_error, 2

    directory = '/data2/daughton/mime1836/by0.0/data/'
    dtxwci    = 2.0
    if arg_present(ascii_version) then ascii_version = 2
    if arg_present(fMap_dir)      then fMap_dir      = '/home/guanlai/pic/by00/'
    if arg_present(eCountFactor)  then eCountFactor  = 2

    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)

    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex*1.e5/6.38), 2)

        ;Filenames
        case tIndex of
            50: eFile = '/data2/daughton/mime1836/by0.0/electrons-783450.bin'
            60: eFile = '/data2/daughton/mime1836/by0.0/electrons-940140.bin'
            else: message, 'No electron file for time index ' + strtrim(tIndex) + '.'
        endcase
    endif

    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = { tIndex: [50, 60], $
                     xrange: [[300, 600], $
                              [300, 600]], $
                     zrange: [[-20,  20], $
                              [-20,  20]] $
                   }
    endif
end


;+
;   Return the electron count factor.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time shown in `EFILE`
;                                 YSLICE: intarr(N), $          y-slice in the .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_mime1836_by05, tIndex, $
    ASCII_INFO=ascii_info, $
    ASCII_VERSION=ascii_version, $
    BINARY_INFO=binary_info, $
    DTXWCI=dtxwci, $
    ECOUNTFACTOR=eCountFactor, $
    EFILE=eFile, $
    EREGIONS=eRegions, $
    FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2
    
    directory = '/data2/daughton/mime1836/by0.05/data/'
    dtxwci = 2.0
    if arg_present(ascii_version) then ascii_version = 2
    if arg_present(fMap_dir)     then fMap_dir     = '/home/guanlai/pic/by05/'
    if arg_present(eCountFactor) then eCountFactor = 2
    
    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif
    
    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex*1.e5/6.38), 2)
    
        ;Filenames
        eFile = '/data2/daughton/mime1836/by0.05/electrons-' + twci + '.bin'
    endif
    
    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = { tIndex: [50, 60], $
                     xrange: [[300, 600], $
                              [300, 600]], $
                     zrange: [[-20,  20], $
                              [-20,  20]] $
                   }
    endif
end


;+
;   Return the electron count factor.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time shown in `EFILE`
;                                 YSLICE: intarr(N), $          y-slice in the .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_mime1836_by10, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2
    
    directory = '/data2/daughton/mime1836/by0.1/data/'
    dtxwci    = 2.0
    if arg_present(ascii_version) then ascii_version = 2
    if arg_present(fMap_dir)      then fMap_dir      = '/home/guanlai/pic/by10/'
    if arg_present(eCountFactor)  then eCountFactor  = 2
    
    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif
    
    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex*1.e5/6.38), 2)
    
        ;Filenames
        eFile = '/data2/daughton/mime1836/by0.1/electrons-' + twci + '.bin'
    endif
    
    ;Regions with electrons
    if arg_present(eRegions) then begin
    eRegions = { tIndex: [50, 60], $
                 xrange: [[300, 600], $
                          [300, 600]], $
                 zrange: [[-20,  20], $
                          [-20,  20]] $
               }
    endif
end


;+
;   Return the electron count factor.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time shown in `EFILE`
;                                 YSLICE: intarr(N), $          y-slice in the .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_mime1836_by40, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2
    
    directory = '/data2/daughton/mime1836/by0.4/data/'
    dtxwci = 2.0
    if arg_present(ascii_version) then ascii_version = 2
    if arg_present(fMap_dir)      then fMap_dir      = '/home/guanlai/pic/by40/'
    if arg_present(eCountFactor)  then eCountFactor  = 2
    
    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif
    
    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex*1.e5/6.38), 2)
    
        ;Filenames
        eFile = '/data2/daughton/mime1836/by0.4/electrons-' + twci + '.bin'
    endif
    
    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = { tIndex: [50, 60], $
                     xrange: [[100, 700], $
                              [300, 600]], $
                     zrange: [[-50,  50], $
                              [-20,  20]] $
        }
    endif
end


;+
;   Return the electron count factor.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time shown in `EFILE`
;                                 YSLICE: intarr(N), $          y-slice in the .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_Asymm3D, tIndex, yslice, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2

    ;DIRECTORY and DTXWCI are needed for ASCII_INFO, BINARY_INFO, and EFILE
    directory = '/data2/Asymm-3D/data/'
    dtxwci    = 8.324609e-04 ;From info file - no t-slices in 3D .gda files
    if arg_present(ascii_version) then ascii_version = 1
    if arg_present(fMap_dir)      then fMap_dir      = '/home/argall/simulations/Asymm-3D/'
    if arg_present(eCountFactor)  then eCountFactor  = 2.0

    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        ;   - y-slice
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex * dtxwci), 2)
        yStr = strtrim(yslice, 2)
        
        ;File name
        eFile = '/data2/Asymm-3D/data/electrons-y' + yStr + '.bin'
    endif
    
    ;Regions with electrons
    ;   - Have moment data for   t=108090
    ;   - Have particle data for t=108108
    if arg_present(eRegions) then begin
        eRegions = { tIndex: [108090, 108090, 108090], $
                     yslice: [650, 905, 1440], $
                     xrange: [[250, 700], $
                              [250, 700], $
                              [250, 700]], $
                     zrange: [[-80,  80], $
                              [-80,  80], $
                              [-80,  80]] $
                   }
    endif
end


;+
;   Configuration program for the Asymm-Large-2D simulation.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time-index into .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_AsymmLarge2D, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2

    directory = '/data2/Asymm-Large-2D/data/'
    dtxwci    = 2.0
    if arg_present(ascii_version) then ascii_version = 1
    if arg_present(fMap_dir)      then fMap_dir      = !Null
    if arg_present(eCountFactor)  then eCountFactor  = !Null

    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        message, 'No particle data available for t*wci = ' + twci + '.', /INFORMATIONAL
        eFile = !Null
    endif
    
    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = !Null
    endif
end


;+
;   Configuration program for the Asymm-Large-2D-NEW simulation.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time-index into .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_AsymmLarge2D_NEW, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2

    ;DIRECTORY and DTXWCI are needed for ASCII_INFO, BINARY_INFO, and EFILE
    directory = '/data2/Asymm-Large-2D-NEW/data/'
    dtxwci    = 5.0
    if arg_present(ascii_version) then ascii_version = 1
    if arg_present(fMap_dir)      then fMap_dir      = '/home/argall/simulations/Asymm-Large-2D-NEW/'

    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex * dtxwci), 2)

        ;Filenames
        eFile = '/data2/Asymm-Large-2D-NEW/electrons-t' + twci + '.bin'
    endif
    
    ;Electron count factor
    if arg_present(eCountFactor) then begin
        case tIndex of
             6: eCountFactor = 2
             9: eCountFactor = 2
            13: eCountFactor = 2
            18: eCountFactor = 2
            26: eCountFactor = 1
            36: eCountFactor = 1
            else: message, 'No particle data for sim # ' + strtrim(simnum,2) + $
                           ' at tIndex ' + strtrim(tIndex,2), /INFORMATIONAL
        endcase
    endif
    
    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = { tIndex: [6, 9, 13, 18, 26, 36], $
                     xrange: [[1410, 1700], $
                              [1250, 1920], $
                              [ 970, 2230], $
                              [ 740, 2620], $
                              [1100, 2100], $
                              [1100, 2400]], $
                     zrange: [[ -60,   60], $
                              [-100,   80], $
                              [-150,  100], $
                              [-200,  150], $
                              [-100,   70], $
                              [-140,   80]] $
                   }
    endif
end


;+
;   Configuration program for the Asymm-Scan/By0 simulation.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time-index into .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_AsymmScan_By0, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2

    ;DIRECTORY and DTXWCI are needed for ASCII_INFO, BINARY_INFO, and EFILE
    directory = '/data2/Asymm-Scan/By0/data/'
    dtxwci    = 2.0
    if arg_present(ascii_version) then ascii_version = 1
    if arg_present(fMap_dir)      then fMap_dir      = '/home/argall/simulations/Asymm-Scan/By0/'
    if arg_present(eCountFactor)  then eCountFactor  = 1L

    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex * dtxwci), 2)

        ;Filenames
        case tIndex of
            28: eFile = '/data2/Asymm-Scan/By0/electrons-twci' + twci + '-it103604.bin'
            else: message, 'No particle data available for t*wci = ' + twci + '.'
        endcase
    endif
    
    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = !Null ;Unknown
    endif
end


;+
;   Configuration program for the Asymm-Scan/By1 simulation.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time-index into .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_AsymmScan_By1, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2

    ;DIRECTORY and DTXWCI are needed for ASCII_INFO, BINARY_INFO, and EFILE
    directory = '/data2/Asymm-Scan/By1/data/'
    dtxwci    = 2.0
    if arg_present(ascii_version) then ascii_version = 1
    if arg_present(fMap_dir)      then fMap_dir      = '/home/argall/simulations/Asymm-Scan/By1/'
    if arg_present(eCountFactor)  then eCountFactor  = 1L

    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex * dtxwci), 2)
        
        ;File names
        case tIndex of
            30: eFile = '/data2/Asymm-Scan/By1/electrons-twci' + twci + '-it70200.bin'
            else: begin
                eFile = ''
                message, 'No particle data available for t*wci = ' + twci + '.', /INFORMATIONAL
            endelse
        endcase
    endif

    ;Regions with electrons
    if arg_present(eRegions) then begin
        eRegions = { tIndex: 30, $
                     xrange: [420, 450], $
                     zrange: [-10,  10]  $
                   }
    endif
end


;+
;   Configuration program for the Sim1 simulation.
;
; :Params:
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;
; :Keywords:
;       ASCII_VERSION:      in, optional, type=integer, default=1
;                           Version of the ASCII info file.
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
;       DTXWCI:             out, optional, type=integer
;                           Time between time-slice of of the ".gda" field and moment
;                               files. Multiply `TINDEX` by DTXWCI to get unitless
;                               simulation time t*wci. Note that DTXWCI is different from
;                               the "dtxwci" in the ASCII info file (data is not saved
;                               at every iteration of the simulation).
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       EREGIONS:           out, optional, type=structure
;                           Regions and times for which we have particle data::
;                               { TINDEX: intarr(N), $          time-index into .gda file
;                                 XRANGE: intarr(2xN), $        [xmin, xmax]
;                                 ZRANGE: intarr(2xN) $         [zmin, zmax]
;                               }
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;-
pro MrSim_Which_Sim1, tIndex, $
ASCII_INFO=ascii_info, $
ASCII_VERSION=ascii_version, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir
    compile_opt strictarr
    on_error, 2

    ;DIRECTORY and DTXWCI are needed for ASCII_INFO, BINARY_INFO, and EFILE
    directory = '/data1/sim1/data/'
    dtxwci    = 0.25
    if arg_present(ascii_version) then ascii_version = 1
    if arg_present(fMap_dir)      then fMap_dir     = '/home/argall/simulations/Sim1/'
    if arg_present(eCountFactor)  then eCountFactor = 2L
    if arg_present(dtxwci)        then dtxwci       = 0.25
    
    ;Binary Info File
    ;   - Stored in the data directory
    if arg_present(binary_info)  then binary_info = filepath('info', ROOT_DIR=directory)
    
    ;Ascii Info File
    ;   - Stored one directory up from the data directory
    if arg_present(ascii_info) then begin
        cd, CURRENT=pwd
        cd, directory
        cd, '..'
        cd, CURRENT=ascii_dir
        cd, pwd
        ascii_info = filepath('info', ROOT_DIR=ascii_dir)
    endif

    ;Electron files
    if arg_present(eFile) then begin
        ;Time information
        ;   - Time index indicating time-slices in the .gda field and moment files.
        ;   - t*wci corresponding to those time-slices
        tStr = strtrim(tIndex, 2)
        twci = strtrim(long(tIndex * dtxwci), 2)

        ;Filenames
        eFile = '/data1/sim1/electrons-'  + strtrim(long(tIndex)*1828L, 2) + '.bin'
    endif
    
    ;Regions with electron data
    ;   - Use MrSim_Create_fMap with the VERBOSE keyword to find out.
    if arg_present(eRegions) then begin
        eRegions = { tIndex: [52, 68, 72, 76, 80, 84, 88, 92, 104, 116], $
                     xrange: [[200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400], $
                              [200, 1400]], $
                     zrange: [[-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50], $
                              [-50,   50]] $
                   }
    endif
end


;+
;   Return simulation-specific information.
;
; :Params:
;       THISSIM:            in, optional, type=string/integer
;                           Either the name or the number of the simulation for which
;                               information is desired. If not provided, a list of
;                               simulations will be printed to the command window.
;
; :Keywords:
;       ASYMMETRIC:         out, optional, type=boolean
;                           Returns 1 if the simulation is asymmetic and 0 if it is
;                               symmetric.
;       DIMENSION:          out, optional, type=string
;                           Returns the dimension of the simulation. Possible values are::
;                               '2D'
;                               '3D'
;       DIRECTORY:          out, optional, type=string
;                           Returns the directory in which the simulation field and moment
;                               data can be found.
;       NAME:               out, optional, type=string
;                           Returns the name of the simulation.
;       SIMNUM:             out, optional, type=integer
;                           Returns the number of the simulation.
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;       _REF_EXTRA:         in, optional, type=any
;                           The following output keywords are also accepted::
;                               ASCII_INFO   - Ascii info file
;                               BINARY_INFO  - Binary info file
;                               DTXWCI       - Time interval between time-slices in .gda files
;                               ECOUNTFACTOR - Electron count multiplication factor.
;                               EFILE        - Electron particle file.
;                               EREGIONS     - Regions and time slices that have electron data.
;                               FMAP_DIR     - Directory containing fMaps.
;-
pro MrSim_Which, thisSim, $
ASYMMETRIC=asymmetric, $
DIMENSION=dimension, $
DIRECTORY=directory, $
NAME=name, $
NUMBER=number, $
TINDEX=tIndex, $
YSLICE=yslice, $
_REF_EXTRA=extra
  compile_opt strictarr
  on_error, 2
    
    ;Information about each simulation.
    simInfo = [['#',  'D',   'Asymmetric', 'Name',                'Directory'], $
               ['1',  '2D',     'yes',     'Asymm-Large-2D',      '/data2/Asymm-Large-2D/data/'], $
               ['2',  '2D',     'yes',     'Asymm-Large-2D-NEW',  '/data2/Asymm-Large-2D-NEW/data/'], $
               ['3',  '3D',     'yes',     'Asymm-3D',            '/data2/Asymm-3D/data/'], $
               ['4',  '2D',     'yes',     'Asymm-2D',            '/data2/Asymm-2D/data/'], $
               ['5',  '2D',     'no',      'sim1',                '/data1/sim1/data/'], $
               ['6',  '2D',     'no',      'mime400',             '/data1/sim1/mime400/data/'], $
               ['7',  '2D',     'no',      'mime400-b',           '/data1/sim1/mime400-b/data/'], $
               ['8',  '2D',     'no',      'data-by0.03-OLD',     '/data1/sim1/data-by0.03-OLD/data/'], $
               ['9',  '2D',     'no',      'data-by0.03-NEW',     '/data1/sim1/data-by0.03-NEW/data/'], $
               ['10', '2D',     'yes',     'Asymm-Scan/By1',      '/data2/Asymm-Scan/By1/data/'], $
               ['11', '2D',     'yes',     'Asymm-Scan/By0',      '/data2/Asymm-Scan/By0/data/'],$
               ['12', '2D',     'no',      'mime1836/by00',       '/data2/daughton/mime1836/by0.0/data/'],$
               ['13', '2D',     'no',      'mime1836/by05',       '/data2/daughton/mime1836/by0.05/data/'],$
               ['14', '2D',     'no',      'mime1836/by10',       '/data2/daughton/mime1836/by0.1/data/'],$
               ['15', '2D',     'no',      'mime1836/by40',       '/data2/daughton/mime1836/by0.4/data/']]
    
    ;Print the info if no input was given.
    if n_elements(thisSim) eq 0 then begin
        nameLen = strtrim(max(strlen(simInfo[3,*])), 2)
        print, simInfo[*,0],   FORMAT='(a2, 4x, a1, 5x, a10, 3x, a-' + nameLen + ', 3x, a0)'
        print, simInfo[*,1:*], FORMAT='(a2, 4x, a2, 7x, a3, 7x, a-' + nameLen + ', 3x, a0)'
        return
    endif
    
    ;Which simulation was chosen?
    ;   - Search by name or number.
    tname = size(thisSim, /TNAME)
    case 1 of
        MrIsA(thisSim, 'STRING'): index = where(strupcase(thisSim) eq strupcase(simInfo[3,1:*]), count)
        MrIsA(thisSim, /INTEGER): index = where(thisSim eq fix(simInfo[0,1:*]), count)
        else: message, 'THISSIM must be a simulation name or number.'
    endcase
    if count ne 1 then message, 'Simulation not found: ' + string(thisSim) + '.'
    index = index[0] + 1
    
    ;Extract the information.
    number     = fix(simInfo[0,index])
    dimension  = simInfo[1,index]
    asymmetric = simInfo[2,index] eq 'yes' ? 1 : 0
    name       = simInfo[3,index]
    directory  = simInfo[4,index]
    
    ;Other simulation-specific parameters.
    if n_elements(extra) gt 0 then begin
        case number of
            1: MrSim_Which_AsymmLarge2D,     tIndex,         _STRICT_EXTRA=extra
            2: MrSim_Which_AsymmLarge2D_NEW, tIndex,         _STRICT_EXTRA=extra
            3: MrSim_Which_Asymm3D,          tIndex, yslice, _STRICT_EXTRA=extra
            4: message, 'Information not available for "' + simname + '".'
            5: MrSim_Which_Sim1,             tIndex,         _STRICT_EXTRA=extra
            6: message, 'Information not available for "' + simname + '".'
            7: message, 'Information not available for "' + simname + '".'
            8: message, 'Information not available for "' + simname + '".'
            9: message, 'Information not available for "' + simname + '".'
            10: MrSim_Which_AsymmScan_By1,    tIndex,        _STRICT_EXTRA=extra
            11: MrSim_Which_AsymmScan_By0,    tIndex,        _STRICT_EXTRA=extra
            12: MrSim_Which_mime1836_by00,    tIndex,        _STRICT_EXTRA=extra
            13: MrSim_Which_mime1836_by05,    tIndex,        _STRICT_EXTRA=extra
            14: MrSim_Which_mime1836_by10,    tIndex,        _STRICT_EXTRA=extra
            15: MrSim_Which_mime1836_by40,    tIndex,        _STRICT_EXTRA=extra
            else: ;Do nothing
        endcase
    endif
end
