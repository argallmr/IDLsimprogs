; docformat = 'rst'
;
; NAME:
;    MrSim_MultiSlab
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
;   The purpose of this program is to create a color image of a single data product with
;   the option of overlaying contours and vertical or horizontal lines.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                           valid simulation object. See MrSim_Create.pro.
;       NAMES:              in, required, type=string
;                           The name of the vector quantity to be plotted in color.
;       TIME:               in, required, type=int
;                           The simulation time for which a velocity vector field
;                               is to be plotted.
; :Keywords:
;       C_NAME:             in, optional, type=string, default=''
;                           Name of the data product whose contours are to be 
;                               overplotted on the image.
;       NLEVELS:            in, optional, type=int, default=15
;                           The number of contour lines to draw. Used with `C_NAME`.
;       OFILENAME:          out, optional, type=string, default=''
;                           If provided, graphics will be output to a file by this
;                               name. The file extension is used to determine which
;                               type of image to create. "PNG" is the default.
;       RANGES:             in, optional, type=fltarr(2)/fltarr(2\,N)
;                           Data range of the quantities to be displayed. If a single
;                                   [min,max] pair is given, it will be applied to all
;                               quantities. Otherwise, there must be N pairs, where N
;                               is the number of `NAMES` given.
;       SIM_OBJECT:         out, optional, type=object
;                           If `THESIM` is the name or number of a simulation, then
;                               this keyword returns the object reference to the
;                               corresponding simulation object that is created.
;       XSIZE:              in, optional, type=long, default=500
;                           Horizontal size of the display window, in pixels.
;       YSIZE:              in, optional, type=long, default=(nSlabs+1)*100 < 690
;                           Vertical size of the display window, in pixels.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
;
; :Returns:
;       msWin:              MrWindow graphic window containing the requested graphics.
;                               If the image data does not exist, an invalid object
;                               will be returned.
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       cgConLevels.pro
;       GetMrWindows.pro
;       MrImage.pro
;       MrPlot.pro
;       MrColorbar.pro
;       MrContour.pro
;       MrSim2D__Define.pro
;       MrSim3D__Define.pro
;       MrWindow.pro
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
;       2014/03/18  -   Written by Matthew Argall
;       2014/09/11  -   Added the RANGES keyword. - MRA
;       2014/09/14  -   Added the XSIZE and YSIZE keyword. - MRA
;       2014/09/30  -   Removed the SIM3D keyword and added the THESIM parameter.
;                           Repurposed the SIM_OBJECT keyword. - MRA
;-
function MrSim_MultiSlab, theSim, names, time, $
C_NAME = c_name, $
NLEVELS = nLevels, $
OFILENAME = ofilename, $
RANGES = ranges, $
SIM_OBJECT = oSim, $
SIM3D = Sim3D, $
XSIZE=xsize, $
YSIZE=ysize, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(colorWin) then obj_destroy, colorWin
        void = cgErrorMSG()
        return, obj_new
    endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
    osim_created = 0B
    
    ;Simulation name or number?
    if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
        oSim = MrSim_Create(theSim, time, _STRICT_EXTRA=extra)
        if obj_valid(oSim) eq 0 then return, obj_new()
        osim_created = 1B
        
    ;Object?
    endif else if MrIsA(theSim, 'OBJREF') then begin
        if obj_isa(theSim, 'MRSIM') eq 0 $
            then message, 'THESIM must be a subclass of the MrSim class.' $
            else oSim = theSim
            
    ;Somthing else
    endif else begin
        MrSim_Which
        message, 'THESIM must be a simulation name, number, or object.'
    endelse
    sim_class = obj_class(oSim)

;-------------------------------------------------------
; Defaults /////////////////////////////////////////////
;-------------------------------------------------------
    ;Set defaults
    ;   -YSIZE is determined later and is based on the number of slabs.
    Sim3D   = keyword_set(Sim3D)
    if n_elements(c_name)     eq 0 then c_name    = ''
    if n_elements(ofilename)  eq 0 then ofilename = ''
    if n_elements(xsize)      eq 0 then xsize     = 500
    
    ;Buffer the output?
    if ofilename eq '' then buffer = 0 else buffer = 1
    
    ;Number of plots to make
    nNames  = n_elements(names)
    nRanges = n_elements(ranges)/2
    
;-------------------------------------------------------
; Color Images /////////////////////////////////////////
;-------------------------------------------------------
    ;Create a window
    msWin = MrWindow(LAYOUT=[1,nNames], XSIZE=_xsize, YSIZE=ysize, YGAP=0.5, $
                     OXMARGIN=[10,14], BUFFER=buffer, REFRESH=0)

    ;Step through each name
    for i = 0, nNames - 1 do begin
        ;Set the range?
        if nRanges gt 0 then range = nRanges eq 1 ? ranges : ranges[*,i]

        ;Create the image.
        !Null = MrSim_ColorSlab(oSim, names[i], /CURRENT, $
                                RANGE      = range, $
                                C_NAME     = c_name, $
                                NLEVELS    = nLevels)
    endfor
    
    ;Destroy the simulation object
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim

;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------

    ;Prettify by setting axis properties
    allIm = msWin -> Get(/ALL, ISA=['MRIMAGE'], COUNT=count)
    if nNames gt 2 then $
        foreach gfx, allIm[1:nNames-2] do gfx -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
        
    if nNames gt 1 then begin
        allIm[0] -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
        allIm[-1].TITLE=''
    endif
    
    ;Resize the window
    if n_elements(_ysize) eq 0 then msWin.YSIZE = 100*(count+1) < 690

    ;Refresh and output, if requested.
    msWin -> Refresh
    if ofilename ne '' then msWin -> Save, ofilename
    
    return, msWin
end
