; docformat = 'rst'
;
; NAME:
;       MrSim_Trace_Init_Script
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
;   Main-level program for initializing particle tracing.
;
; :History:
;   Modification History::
;       2014-11-13  -   Written by Matthew Argall
;-
;*****************************************************************************************

;Set Inputs
theSim     = 'Asymm-Scan/By0'
time       = 28
bin_center = [367.0,  0.0, 0.0]
half_width = [  1.0,  1.0, 1.0]
vx_range   = [ -0.1,  0.1]
vy_range   = [  0.3,  0.4]
vz_range   = [ -0.2, -0.1]
stepPrint  = 20
stepMax    = 20000
oParamFile = '/home/argall/programs/tracer/inputs/Params_Asymm-Scan-By0_moon.txt'
oInitFile  = '/home/argall/programs/tracer/inputs/Init_Asymm-Scan-By0_moon.txt'
oFilename  = '/home/argall/programs/tracer/savedat/Data_Asymm-Scan-By0_moon'

;Create the initialization file
MrSim_Trace_Init, theSim, time, bin_center, half_width, vx_range, vy_range, vz_range, $
                  OPARAMFILE=oParamFile, OINITFILE=oInitFile, OFILENAME=oFileName, $
                  STEPPRINT=stepPrint, STEPMAX=stepMax

end