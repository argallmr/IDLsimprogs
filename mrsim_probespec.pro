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
function MrSim_ProbeSpec, filename
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        if obj_valid(win) then obj_destroy, win
        void = cgErrorMSG()
        return, !Null
    endif
    
    fname1 = '/home/argall/Work/probe1.dat'
    fname2 = '/home/argall/Work/probe2.dat'
    fname3 = '/home/argall/Work/probe3.dat'
    
    ;Ensure the file exists
    if file_test(fname1) eq 0 then begin
        filename = dialog_pickfile(/OPEN, TITLE='Choose Data Probe File')
        if filename eq '' then return, !Null
    endif
    
    ;Read the data
    probe1 = MrRead_Ascii( fname1, $
                           DATA_START   = 1, $
                           GROUPS       = [0, 1, 1, 2, 2, 2, 3, 3, 3, 4], $
                           COLUMN_NAMES = ['Time', 'XZ', 'XZ', 'E', 'E', 'E', 'B', 'B', 'B', 'n_e'] $
                         )
    probe2 = MrRead_Ascii( fname2, $
                           DATA_START   = 1, $
                           GROUPS       = [0, 1, 1, 2, 2, 2, 3, 3, 3, 4], $
                           COLUMN_NAMES = ['Time', 'XZ', 'XZ', 'E', 'E', 'E', 'B', 'B', 'B', 'n_e'] $
                         )
    probe3 = MrRead_Ascii( fname3, $
                           DATA_START   = 1, $
                           GROUPS       = [0, 1, 1, 2, 2, 2, 3, 3, 3, 4], $
                           COLUMN_NAMES = ['Time', 'XZ', 'XZ', 'E', 'E', 'E', 'B', 'B', 'B', 'n_e'] $
                         )
    
    ;Ranges
    nRange = [ min( [probe1.n_e, probe2.n_e, probe3.n_e], MAX=nMax ), nMax ]
    bRange = [ min( [probe1.B,   probe2.B,   probe3.B],   MAX=bMax ), bMax ]
    eRange = [ min( [probe1.E,   probe2.E,   probe3.E],   MAX=eMax ), eMax ]
    
;    ;Determine sample rate.
;    nPts = n_elements(data.Time)
;    dt   = mean(probe1.Time[1:nPts-1] - probe1.Time[0:nPts-2])
;
;    ;Create the spectrograms
;    B1_psd = MrPSD(probe1.B, nPts, dt, DIMENSION=2, T0=t0, FREQUENCIES=f, TIME=t)
;    B2_psd = MrPSD(probe2.B, nPts, dt, DIMENSION=2, T0=t0)
;    B3_psd = MrPSD(probe3.B, nPts, dt, DIMENSION=2, T0=t0)
;    E1_psd = MrPSD(probe1.E, nPts, dt, DIMENSION=2, T0=t0)
;    E2_psd = MrPSD(probe2.E, nPts, dt, DIMENSION=2, T0=t0)
;    E3_psd = MrPSD(probe3.E, nPts, dt, DIMENSION=2, T0=t0)
;    f      = f * 0.1

    ;Create a window
    win = MrWindow(LAYOUT=[3,3], OXMARGIN=[10,5], REFRESH=0, XGAP=0.5, XSIZE=800, YGAP=0, YSIZE=600)
    
;-----------------------------------------------------
; B Field \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
;    pbRange = [min(B_psd, MAX=pMax), pMax]
;    peRange = [min(E_psd, MAX=eMax), eMax]
;
;    positions = [[0.23, 0.81, 0.92, 0.94], $
;                 [0.23, 0.67, 0.92, 0.81], $
;                 [0.23, 0.53, 0.92, 0.67], $
;                 [0.23, 0.26, 0.92, 0.40], $
;                 [0.23, 0.12, 0.92, 0.26]]

    ;t-n
    n1plot = MrPlot(probe1.Time, probe1.n_e, $
                    /CURRENT, $
                    NAME          = 'Probe1: ne', $
                    TITLE         = 'Probe1', $
                    XTICKFORMAT   = '(a1)', $
                    YRANGE        = nRange, $
                    YTICKINTERVAL = 0.2, $
                    YTITLE        = 'n$\downe$')

    ;t-n
    n2plot = MrPlot(probe2.Time, probe2.n_e, $
                    /CURRENT, $
                    NAME          = 'Probe2: ne', $
                    TITLE         = 'Probe2', $
                    XTICKFORMAT   = '(a1)', $
                    YRANGE        = nRange, $
                    YTICKFORMAT   = '(a1)')

    ;t-n
    n3plot = MrPlot(probe3.Time, probe3.n_e, $
                    /CURRENT, $
                    NAME          = 'Probe3: ne', $
                    TITLE         = 'Probe3', $
                    XTICKFORMAT   = '(a1)', $
                    YRANGE        = nRange, $
                    YTICKFORMAT   = '(a1)')

    ;t-B
    B1plot = MrPlot(probe1.Time, probe1.B, $
                    /CURRENT, $
                    NAME        = 'Probe1: Bxyz', $
                    DIMENSION   = 2, $
                    XTICKFORMAT = '(a1)', $
                    YRANGE      = bRange, $
                    YTITLE      = 'B')

    ;t-B
    B2plot = MrPlot(probe2.Time, probe2.B, $
                    /CURRENT, $
                    NAME        = 'Probe2: Bxyz', $
                    DIMENSION   = 2, $
                    XTICKFORMAT = '(a1)', $
                    YRANGE      = bRange, $
                    YTICKFORMAT = '(a1)')

    ;t-B
    B3plot = MrPlot(probe3.Time, probe3.B, $
                    /CURRENT, $
                    NAME        = 'Probe3: Bxyz', $
                    DIMENSION   = 2, $
                    XTICKFORMAT = '(a1)', $
                    YRANGE      = bRange, $
                    YTICKFORMAT = '(a1)')

    ;t-E
    E1plot = MrPlot(probe1.Time, probe1.E, $
                    /CURRENT, $
                    NAME          = 'Probe1: Exyz', $
                    DIMENSION     = 2, $
                    XTITLE        = 'Time (0.02$\Omega$$\downi$$\up-1$)', $
                    YRANGE        = eRange, $
                    YTICKINTERVAL = 0.01, $
                    YTITLE        = 'E')

    ;t-E
    E2plot = MrPlot(probe2.Time, probe2.E, $
                    /CURRENT, $
                    NAME        = 'Probe2: Exyz', $
                    DIMENSION   = 2, $
                    XTITLE      = 'Time (0.02$\Omega$$\downi$$\up-1$)', $
                    YRANGE      = eRange, $
                    YTICKFORMAT = '(a1)')

    ;t-E
    E3plot = MrPlot(probe3.Time, probe3.E, $
                    /CURRENT, $
                    NAME        = 'Probe3: Exyz', $
                    DIMENSION     = 2, $
                    XTITLE        = 'Time (0.02$\Omega$$\downi$$\up-1$)', $
                    YRANGE        = eRange, $
                    YTICKFORMAT   = '(a1)')

    
    ;E Legend
    Elegend = MrLegend(/BOX, $
                       BX_COLOR = 'White', $
                       CHARSIZE = 2.0, $
                       LENGTH   = 0, $
                       LOCATION = 8, $
                       TARGET   = E3plot, $
                       TCOLORS  = ['Blue', 'Forest Green', 'Red'], $
                       TITLES   = ['Ex', 'Ey', 'Ez'])
    
    ;B Legend
    Blegend = MrLegend(/BOX, $
                       BX_COLOR = 'White', $
                       CHARSIZE = 2.0, $
                       LENGTH   = 0, $
                       LOCATION = 8, $
                       TARGET   = B3plot, $
                       TCOLORS  = ['Blue', 'Forest Green', 'Red'], $
                       TITLES   = ['Bx', 'By', 'Bz'])
    
    Elegend -> Order, /SEND_TO_BACK
    Blegend -> Order, /SEND_TO_BACK
;-----------------------------------------------------
; E Field \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

    win -> SetGlobal, CHARSIZE=2.0, CHARTHICK=2
    win -> Refresh
    return, win
end
