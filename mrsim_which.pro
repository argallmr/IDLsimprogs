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
;       2013/09/12  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Return the ASCII Info file.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       ASCII_INFO:         ASCII info file with information about the particular
;                               simulation.
;-
function MrSim_Which_ASCII_Info, simnum
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(pwd) gt 0 then cd, pwd
        void = cgErrorMSG(/QUIET)
        return, !Null
    endif
        
    ;Get the data directory
    MrSim_Which, simnum, DIRECTORY=directory

    ;Multiplier of electron counts
    case simnum of
        ;The ASCII info file is stored one level up from the data directory
        else: begin
            cd, CURRENT=pwd
            cd, directory
            cd, '..'
            cd, CURRENT=ascii_dir
            cd, pwd
            ascii_info = filepath('info', ROOT_DIR=ascii_dir)
        endcase
    endcase
    
    return, ascii_info
end

;+
;   Return the binary Info file.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       BINARY_INFO:        Binary info file with information about the particular
;                               simulation.
;-
function MrSim_Which_Binary_Info, simnum
    compile_opt strictarr
    on_error, 2
        
    ;Get the data directory
    MrSim_Which, simnum, DIRECTORY=directory
    
    ;Multiplier of electron counts
    case simnum of
        ;The binary info file is stored in the data directory
        else: binary_info = filepath('info', ROOT_DIR=directory)
    endcase
    
    return, binary_info
end


;+
;   Return the time-step at which filed and moment data was saved.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       DTXWCI:             Time step at which simulation data was saved.
;-
function MrSim_Which_dtxwci, simnum
    compile_opt strictarr
    on_error, 2

    ;Hard code -- the one simulation that we have has a subset of time indices.
    case simnum of
         1:   dtxwci = 2.0
         2:   dtxwci = 5.0
         3:   dtxwci = 8.324609e-04     ;From info file.
         5:   dtxwci = 0.25
        10:   dtxwci = 2.0
        11:   dtxwci = 2.0
        else: dtxwci = !Null
    endcase
    
    return, dtxwci
end


;+
;   Return the electron count factor.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       ECOUNTFACTOR:       Multiplication factor of electron counts. Typically, only a
;                               fraction of the number of electrons are saved.
;-
function MrSim_Which_eCountFactor, simnum, tIndex
    compile_opt strictarr
    on_error, 2
    
    ;Multiplier of electron counts
    case simnum of
         2: begin
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
         endcase
         3:   eCountFactor = 2
         5:   eCountFactor = 2
        10:   eCountFactor = 1
        11:   eCountFactor = 1
        else: eCountFactor = !Null
    endcase
    
    return, eCountFactor
end


;+
;   Return the file name of the binary file containing electron data.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;
; :Returns:
;       EFILE:              Name of the file in which electron data is saved.
;-
function MrSim_Which_eFile, simnum, tIndex, yslice, $
DTXWCI=dtxwci
    compile_opt strictarr
    on_error, 2

    ;Get the interval between time slices.
    if n_elements(dtxwci) eq 0 then dtxwci = MrSim_Which_dtxwci(simnum)

    ;Time information
    ;   - Time index
    ;   - t*wci derived from the saved time-step
    ;   - t*wci derived from the simulation time-step
    tStr = strtrim(tIndex, 2)
    twci = strtrim(long(tIndex * dtxwci), 2)
    if n_elements(yslice) ne 0 then yStr = strtrim(yslice, 2)

    ;Guess a file name.
    case simnum of
         2:   eFile = '/data2/Asymm-Large-2D-NEW/'    + 'electrons-t' + twci + '.bin'
         3:   eFile = '/data2/Asymm-3D/data/'         + 'electrons-y' + yStr + '.bin'
         5:   eFile = '/data1/sim1/'                  + 'electrons-'  + strtrim(long(tIndex*1828), 2) + '.bin'
         6:   eFile = '' ;'/data1/sim1/mime400/particles/' + 'electron' + + '.dat'
         9:   eFile = '/data1/sim1/particles-by0.03/' + 'electrons-t' + dtxwci + '.bin'
        10: begin
            case tIndex of
                30: eFile = '/data2/Asymm-Scan/By1/electrons-twci' + twci + '-it70200.bin'
                else: message, 'No particle data available for t*wci = ' + twci + '.'
            endcase
        endcase
        11: begin
            case tIndex of
                28: eFile = '/data2/Asymm-Scan/By0/electrons-twci' + twci + '-it103604.bin'
                else: message, 'No particle data available for t*wci = ' + twci + '.'
            endcase
        endcase
        else: eFile = !Null
    endcase
    
    return, eFile
end


;+
;   Return the regions for which we have electron data.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       EREGIONS:           A structure containing the time index, yslice, xrange, and
;                               zrange over which we have particle data.
;-
function MrSim_Which_eRegions, simnum
    compile_opt strictarr
    on_error, 2
    
    ;Regions for which we have particle data
    case simnum of
         2: eRegions = { tIndex: [30, 45, 65, 90, 130, 180], $
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
         3: eRegions = { tIndex: [120816, 120816], $
                         yslice: [860, 1400], $
                         xrange: [[285, 585], $
                                  [290, 590]], $
                         zrange: [[-40,  50], $
                                  [-60,  80]] $
                       }
         5: eRegions = { tIndex: [52, 68, 72, 76, 80, 84, 88, 92, 104, 116], $
                         xrange: [[  0,    0], $
                                  [200, 1400], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0]], $
                         zrange: [[  0,    0], $
                                  [-50,   50], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0], $
                                  [  0,    0]] $
                       }
        10: eRegions = { tIndex: 30, $
                         xrange: [420, 450], $
                         zrange: [-10,  10]  $
                       }
        11: Message, 'EREGION unknown for simulation #11.', /INFORMATIONAL
        else: eRegions = !Null
    endcase
    
    return, eRegions
end


;+
;   Return the electron count factor.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       FMAP_DIR:           Directory in which to find fMaps.
;-
function MrSim_Which_fMapDir, simnum
    compile_opt strictarr
    on_error, 2

    ;fMap directory
    case simnum of
         2:   fMap_dir = '/home/argall/simulations/Asymm-Large-2D-NEW/'
         3:   fMap_dir = '/home/argall/simulations/Asymm-3D/'
         5:   fMap_dir = '/home/argall/simulations/Sim1/'
        10:   fMap_dir = '/home/argall/simulations/Asymm-Scan/By1/'
        11:   fMap_dir = '/home/argall/simulations/Asymm-Scan/By0/'
        else: fMap_dir = !Null
    endcase
    
    return, fMap_dir
end


;+
;   Return the electron count factor.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       FMAP_DIR:           Directory in which to find fMaps.
;-
function MrSim_Which_AsymmScan_By0, simnum, tIndex, $
ASCII_INFO=ascii_info, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
FMAP_DIR=fMap_dir, $
    compile_opt strictarr
    on_error, 2

    if arg_present(directory)    then directory    = '/data2/Asymm-Scan/By0/data/'
    if arg_present(fMap_dir)     then fMap_dir     = '/home/argall/simulations/Asymm-Scan/By0/'
    if arg_present(eCountFactor) then eCountFactor = 1L
    if arg_present(dtxwci)       then dtxwci       = 2.0

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
;   Return the electron count factor.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       FMAP_DIR:           Directory in which to find fMaps.
;-
function MrSim_Which_AsymmScan_By1, simnum, tIndex, $
ASCII_INFO=ascii_info, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
FMAP_DIR=fMap_dir, $
    compile_opt strictarr
    on_error, 2

    if arg_present(directory)    then directory    = '/data2/Asymm-Scan/By1/data/'
    if arg_present(fMap_dir)     then fMap_dir     = '/home/argall/simulations/Asymm-Scan/By1/'
    if arg_present(eCountFactor) then eCountFactor = 1L
    if arg_present(dtxwci)       then dtxwci       = 2.0

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
;   Return the electron count factor.
;
; :Params:
;       SIMNUM:             in, optional, type=string/integer
;                           Simulation number for which dtxwci is desired.
;
; :Returns:
;       FMAP_DIR:           Directory in which to find fMaps.
;-
function MrSim_Which_Sim1, simnum, tIndex, $
ASCII_INFO=ascii_info, $
BINARY_INFO=binary_info, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
FMAP_DIR=fMap_dir, $
    compile_opt strictarr
    on_error, 2

    if arg_present(directory)    then directory    = '/data2/sim1/'
    if arg_present(fMap_dir)     then fMap_dir     = '/home/argall/simulations/Sim1/'
    if arg_present(eCountFactor) then eCountFactor = 2L
    if arg_present(eRegions)     then eRegions     = !Null
    if arg_present(dtxwci)       then dtxwci       = 0.25

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
                     xrange: [[  0,    0], $
                              [200, 1400], $
                              [  0,    0], $
                              [  0,    0], $
                              [200, 1400], $
                              [  0,    0], $
                              [  0,    0], $
                              [  0,    0], $
                              [  0,    0], $
                              [  0,    0], $
                              [200, 1400]], $
                     zrange: [[  0,    0], $
                              [-50,   50], $
                              [  0,    0], $
                              [  0,    0], $
                              [-50,   50], $
                              [  0,    0], $
                              [  0,    0], $
                              [  0,    0], $
                              [  0,    0], $
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
;       ASCII_INFO:         out, optional, type=string
;                           ASCII info file containing human readable information about
;                               the simlation.
;       BINARY_INFO:        out, optional, type=string
;                           Binary info file containing information about the simulation.
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
;       DTXWCI:             out, optional, type=integer
;                           Time step at which simulation data is saved. This is different
;                               from the time step at which the simulation is evolved.
;                               For the latter, see the info file. Multiply `TINDEX` by
;                               DTXWCI to get unitless simulation time t*wci.
;       ECOUNTFACTOR:       out, optional, type=integer, default=1
;                           Factor by which electrons need to be multiplied. Often, only
;                               every other particle is saved.
;       EFILE:              out, optional, type=string
;                           File in which electron data is saved. Requires `TINDEX` and,
;                               for 3D simulations, `YSLICE`.
;       FMAP_DIR:           out, optional, type=string
;                           Directory in which fMaps are saved.
;       NAME:               out, optional, type=string
;                           Returns the name of the simulation.
;       SIMNUM:             out, optional, type=integer
;                           Returns the number of the simulation.
;       TINDEX:             in, optional, type=long
;                           Simulation time index.
;       YSLICE:             in, optional, type=long
;                           Y-slice within a 3D simulation.
;-
pro MrSim_Which, thisSim, $
ASYMMETRIC=asymmetric, $
ASCII_INFO=ascii_info, $
BINARY_INFO=binary_info, $
DIMENSION=dimension, $
DIRECTORY=directory, $
DTXWCI=dtxwci, $
ECOUNTFACTOR=eCountFactor, $
EFILE=eFile, $
EREGIONS=eRegions, $
FMAP_DIR=fMap_dir, $
NAME=name, $
NUMBER=number, $
TINDEX=tIndex, $
YSLICE=yslice
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
               ['11', '2D',     'yes',     'Asymm-Scan/By0',      '/data2/Asymm-Scan/By0/data/']]

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
    if count ne 1 then message, 'Simulation not found: ' + string(thisSim)
    index = index[0] + 1
    
    ;Extract the information.
    number     = fix(simInfo[0,index])
    dimension  = simInfo[1,index]
    asymmetric = simInfo[2,index] eq 'yes' ? 1 : 0
    name       = simInfo[3,index]
    directory  = simInfo[4,index]

    ;Other simulation-specific parameters.
    if arg_present(ascii_info)   then ascii_info   = MrSim_Which_ASCII_Info(number)
    if arg_present(binary_info)  then binary_info  = MrSim_Which_Binary_Info(number)
    if arg_present(dtxwci)       then dtxwci       = MrSim_Which_dtxwci(number)
    if arg_present(eCountFactor) then eCountFactor = MrSim_Which_eCountFactor(number, tIndex)
    if arg_present(eFile)        then eFile        = MrSim_Which_eFile(number, tIndex, yslice)
    if arg_present(eRegions)     then eRegions     = MrSim_Which_eRegions(number)
    if arg_present(fmap_dir)     then fmap_dir     = MrSim_Which_fMapDir(number)
end
