; docformat = 'rst'
;
; NAME:
;    MrSim_ProbeSpec
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
;   Read data from a simulation probe.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       FILENAME:       in, optional, type=string, default=''
;                       Name of the probe data file to read. If no file is given, or the
;                           file does not exist, a dialog box will open asking you to
;                           choose a file.
;
; :Returns:
;       DATA:           out, required, type=structure
;                       A structure of probe data having fields EX, EY, EZ,
;                           BX, BY, BY, and NE.
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       MrReadAscii.pro
;       MrPSD.pro
;       MrWindow.pro
;       MrPlot.pro
;       MrImage.pro
;       MrColorbar.pro
;       MrLegend.pro
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
;       2014/10/19  -   Written by Matthew Argall
;-
function MrSim_ProbeGrid
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(win) then obj_destroy, win
        if obj_valid(oFD) then obj_destroy, oFD
        void = cgErrorMSG()
        return, !Null
    endif
    
    point     = [46, 102, 1]
    directory = '/data2/Asymm-3D/fd/'
    
    ;Create Field Data object
    oFD = MrFD(DIRECTORY=directory)
    
    ;Get the probe location
    theProbe = oFD -> GetProbe(point[0], point[1], point[2])
    
    ;Read data from the probe
    Bx = oFD -> ReadProbe('Bx', theProbe[0], theProbe[2], theProbe[1], /TIME)
    By = oFD -> ReadProbe('By', theProbe[0], theProbe[2], theProbe[1], /TIME)
    Bz = oFD -> ReadProbe('Bz', theProbe[0], theProbe[2], theProbe[1], /TIME)
    Ex = oFD -> ReadProbe('Ex', theProbe[0], theProbe[2], theProbe[1], /TIME)
    Ey = oFD -> ReadProbe('Ey', theProbe[0], theProbe[2], theProbe[1], /TIME)
    Ez = oFD -> ReadProbe('Ez', theProbe[0], theProbe[2], theProbe[1], /TIME)
    obj_destroy, oFD
    
    ;Combine the data
    t = findgen(n_elements(Bx))
    B = transpose([[temporary(Bx)], [temporary(By)], [temporary(Bz)]])
    E = transpose([[temporary(Ex)], [temporary(Ey)], [temporary(Ez)]])

    ;Plot the data
    Bplot = MrPlot(t, B, $
                   DIMENSION   = 2, $
                   NAME        = 'Bxyz', $
                   TITLE       = string(FORMAT='(%"Probe (%i, %i, %i)")', theProbe), $
                   XTICKFORMAT = '(a1)', $
                   YTITLE      = 'B')
    Bplot -> Refresh, /DISABLE
                   
    Eplot = MrPlot(t, E, $
                   /CURRENT, $
                   DIMENSION     = 2, $
                   NAME          = 'Exyz', $
                   XTITLE        = 'Time', $
                   YTICKINTERVAL = 1.0, $
                   YTITLE        = 'E')
                   
    Blegend = MrLegend(/BOX, $
                       BX_COLOR = 'White', $
                       CHARSIZE = 2.0, $
                       LENGTH   = 0, $
                       LOCATION = 8, $
                       TARGET   = Bplot, $
                       TCOLORS  = ['Blue', 'Forest Green', 'Red'], $
                       TITLES   = ['Bx', 'By', 'Bz'])
                   
    Elegend = MrLegend(/BOX, $
                       BX_COLOR = 'White', $
                       CHARSIZE = 2.0, $
                       LENGTH   = 0, $
                       LOCATION = 8, $
                       TARGET   = Eplot, $
                       TCOLORS  = ['Blue', 'Forest Green', 'Red'], $
                       TITLES   = ['Ex', 'Ey', 'Ez'])

    ;Set window properties
    win = Bplot.window
    win -> SetProperty, YGAP=0, XSIZE=500, OXMARGIN=[10, 7]
    win -> SetGlobal, CHARSIZE=2.0
    win -> Refresh
    return, win
end
