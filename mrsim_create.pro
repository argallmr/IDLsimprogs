; docformat = 'rst'
;
; NAME:
;    MrSim_Create
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
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       THESIM:         in, optional, type=string/integer
;                       Name or index of the simulation to be created. If not present,
;                           a list of all available simulations will be printed to the
;                           display window.
;       TIME:           in, optional, type=long, default=0
;                       Time index at which simulation data will be read.
;       YSLICE:         in, optional, type=long, default=0
;                       The Y-index at which to obtain data from the XZ plane. Ignored
;                           in 2D simulations.
;
; :Keywords:
;       ASCII_VERSION:  in, optional, type=integer, default=1
;                       Version of the ASCII info file. Ignored if `BINDARY`=1.
;                           See MrSim_Which.pro.
;       ASCII_INFO:     in, optional, type=string, default=`DIRECTORY`/../info
;                       File name of the ASCII info file containing information abou the
;                           simulation to be used.
;       AXIS_LABELS:    in, optional, type=strarr(3), default="['x', 'y', 'z']"
;                       Labels for the axes.
;       BINARY:         in, optional, type=boolean, default=0
;                       If set, `INFO_FILE` points to the binary info file.
;       BINARY_INFO:    in, optional, type=string, default=`DIRECTORY`/info
;                       File name of the binary info file containing information abou the
;                           simulation to be used.
;       COORD_SYSTEM:   in, optional, type=string, default='SIMULATION'
;                       Coordinate system in which to display the data. Options are::
;                           'SIMULATION'
;                           'MAGNETOPAUSE'
;                           'MAGNETOTAIL'
;       DIRECTORY:      in, optional, type=string, default=pwd
;                       Directory in which to find the ".gda" data.
;       INFO_FILE:      in, optional, type=string, default=`DIRECTORY`/../info`
;                       The ASCII info file containing information about the simulation
;                           setup. If `BINARY` is set, the default file will be
;                           `DIRECTORY`/info.
;       ION_SCALE:      in, optional, type=boolean, default=0
;                       Construct the simulation domain in units of "di" instead of "de"
;                           (the ion and electron skin depth, respectively).
;       MVA_FRAME:      in, optional, type=boolean, default=0
;                       If `COORD_SYSTEM` is "Magnetopause" or "Magnetotail", then this
;                           this keyword will automatically pick the proper `AXIS_LABELS`.
;       NSMOOTH:        in, optional, type=integer, default=3
;                       Number of points to smooth data after they is read.
;       ORIENTATION:    in, optional, type=string, default='XZ'
;                       The 2D plane in which data will be viewed. Also sets the direction
;                           of "vertical" and "horizontal" when making 1D cuts. Possible
;                           choices are "XY" and "XZ". This keyword has no meaning for
;                           2D simulations.
;       SIM_INFO:       in, optional, type=boolean, default=0
;                       If set, information about each simulation will be printed to the
;                           command window. Initialization will be aborted.
;       XRANGE:         in, optional, type=dblarr(2), default=width of domain
;                       X-range over which data will be read and stored, in units defined
;                           by `ION_SCALE`.
;       YRANGE:         in, optional, type=dblarr(2), default=depth of domain
;                       Y-range over which data will be read and stored, in units defined
;                           by `ION_SCALE`. This keyword is ignored in 2D simulations.
;       ZRANGE:         in, optional, type=dblarr(2), default=height of domain
;                       Z-range over which data will be read and stored, in units defined
;                           by `ION_SCALE`.
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
;       2014/10/03  -   Added the ASCII_VERSION keyword. - MRA
;       2014/10/11  -   Forgot the DIRECTORY keyword. Fixed. - MRA
;       2014/11/04  -   Added BINARY_INFO and ASCII_INFO keywords. Depricated the
;                           INFO_FILE and BINARY keywords. - MRA
;-
function MrSim_Create, thisSim, time, yslice, $
ASCII_VERSION = ascii_version, $
AXIS_LABELS = axis_labels, $
BINARY = binary, $
COORD_SYSTEM = coord_system, $
DIRECTORY = directory, $
INFO_ASCII = info_ascii, $
INFO_BINARY = info_binary, $
ION_SCALE = ion_scale, $
MVA_FRAME = mva_frame, $
NSMOOTH = nsmooth, $
ORIENTATION = orientation, $
SIM_INFO = sim_info, $
XRANGE = xrange, $
YRANGE = yrange, $
ZRANGE = zrange
    compile_opt strictarr
    on_error, 2
    
    ;Which Simulation was chosen?
    MrSim_Which, thisSim, NAME=simname, DIMENSION=dimension
    if n_elements(simname) eq 0 then return, !Null

    ;Which simulation object should be used?
    class = dimension eq '2D' ? 'MrSim2D' : 'MrSim3D'

    ;Create the simulation object
    sim_object = obj_new(class, simname, time, yslice, $
                         ASCII_VERSION = ascii_version, $
                         AXIS_LABELS   = axis_labels, $
                         BINARY        = binary, $
                         COORD_SYSTEM  = coord_system, $
                         DIRECTORY     = directory, $
                         INFO_ASCII    = info_ascii, $
                         INFO_BINARY   = info_binary, $
                         ION_SCALE     = ion_scale, $
                         MVA_FRAME     = mva_frame, $
                         NSMOOTH       = nsmooth, $
                         ORIENTATION   = orientation, $
                         SIM_INFO      = sim_info, $
                         XRANGE        = xrange, $
                         YRANGE        = yrange, $
                         ZRANGE        = zrange)

    return, sim_object
end