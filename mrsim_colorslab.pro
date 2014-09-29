; docformat = 'rst'
;
; NAME:
;    MrSim_ColorSlab
;
;*****************************************************************************************
;   Copyright (c) 2013, Matthew Argall                                                   ;
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
;       NAME:                   in, required, type=string
;                               The name of the vector quantity to be plotted in color.
;       TIME:                   in, required, type=int
;                               The simulation time for which a velocity vector field
;                                   is to be plotted.
; :Keywords:
;       C_NAME:                 in, optional, type=string, default=''
;                               Name of the data product whose contours are to be 
;                                   overplotted on the image.
;       CURRENT:                in, optional, type=boolean, default=0
;                               If set, the graphics will be added to the current MrWindow
;                                   graphics window. The default is to create a new window.
;       FRACTION:               in, optional, type=float
;                               Fraction of the vectors to be drawn. Used with `VX_NAME`
;                                   and `VY_NAME`.
;       HORIZ_LINES:            in, optional, type=flt/fltarr
;                               Locations, in data coordinates, along the y-axis of the
;                                   image, at which horizontal lines are to be drawn.
;       LINE_COLOR:             in, optional, type=strarr, default="['Blue', 'Forest Green', 'Red', 'Magenta', 'Orange']"
;                               Colors of the lines to be drawn on the image.
;       NLEVELS:                in, optional, type=int, default=15
;                               The number of contour lines to draw. Used with `C_NAME`.
;       OFILENAME:              out, optional, type=string, default=''
;                               If provided, graphics will be output to a file by this
;                                   name. The file extension is used to determine which
;                                   type of image to create. "PNG" is the default.
;       RANGE:                  in, optional, type=fltarr(2), default=[min\,max]
;                               Data range of the quantity to be displayed.
;       SIM_OBJECT:             out, optional, type=object
;                               The "Sim_Object" (or subclass) reference containing the 
;                                   simulation data and domain information used in
;                                   creating the plot. If not provided, one will be
;                                   created according to the `SIM3D` keyword.
;       SIM3D:                  in, optional, type=boolean, default=0
;                               If set, and `SIM_OBJECT` is not given, then a 3D
;                                   simulation object will be created. The default is
;                                   to make 2D simulation object.
;       VX_NAME:                in, optional, type=string, default=''
;                               Name of the x-component of a data product whose vector
;                                   field is to be overplotted on the image. Must be used
;                                   with `VY_NAME`.
;       VY_NAME:                in, optional, type=string, default=''
;                               Name of the y-component of a data product whose vector
;                                   field is to be overplotted on the image. Must be used
;                                   with `VX_NAME`.
;       VERT_LINES:             in, optional, type=flt/fltarr
;                               Locations, in data coordinates, along the x-axis of the
;                                   image, at which vertical lines are to be drawn.
;       _REF_EXTRA:             in, optional, type=structure
;                               Any keyword accepted Sim_Object::Init(), or any of its
;                                   subclasses is also accepted here. Ignored of
;                                   `SIM_OBJECT` is present.
;
; :Returns:
;       COLORWIN:               MrWindow graphic window containing the requested graphics.
;                                   If the image data does not exist, an invalid object
;                                   will be returned.
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
;       MrVector.pro
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
;       2013/04/30  -   Written by Matthew Argall
;       2014/01/22  -   Contours no longer have labels. Read the Info file to convert
;                           numbers to physical values. Added the ION_SCALE and INFO
;                           keywords. - MRA
;       2014/01/26  -   Added the SIM_OBJECT and T_WCI keywords. Removed the ION_SCALE
;                           and INFO keywords.
;       2014/01/30  -   Added the C_NAME keyword. - MRA
;       2014/01/31  -   Added the 3DSIM keyword. Removed keywords accepted by SIM_OBJECT. - MRA
;       2014/02/07  -   Renamed from sim_color_plot.pro to MrSim_ColorSlab.pro. Added
;                           the V_NAME keyword. Removed the T_WCI keyword. - MRA
;       2014/02/09  -   If the image data does not exist, an invalid object is returned. - MRA
;       2014/03/05  -   Added the LINE_COLOR keyword. - MRA
;       2014/03/26  -   Added ability to overplot vector fields. Changed V_NAME to VX_NAME
;                           and VY_NAME. Added the FRACTION keyword. - MRA
;       2014/09/11  -   Added the RANGE keyword. - MRA
;-
function MrSim_ColorSlab, name, time, $
C_NAME = c_name, $
CURRENT = current, $
FRACTION = fraction, $
HORIZ_LINES = horiz_lines, $
LINE_COLOR = line_color, $
NLEVELS = nLevels, $
OFILENAME = ofilename, $
RANGE = range, $
SIM_OBJECT = oSim, $
SIM3D = Sim3D, $
VX_NAME = vx_name, $
VY_NAME = vy_name, $
VERT_LINES = vert_lines, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim)     && arg_present(oSim)    eq 0 then obj_destroy, oSim
        if obj_valid(colorWin) && keyword_set(current) eq 0 then obj_destroy, colorWin
        void = cgErrorMSG()
        return, obj_new
    endif

;-------------------------------------------------------
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------

    ;Set defaults
    current = keyword_set(current)
    
    if n_elements(c_name)     eq 0 then c_name    = ''
    if n_elements(fraction)   eq 0 then fraction  = 0.005
    if n_elements(vx_name)    eq 0 then vx_name   = ''
    if n_elements(vy_name)    eq 0 then vy_name   = ''
    if n_elements(ofilename)  eq 0 then ofilename = ''
    if n_elements(nlevels)    eq 0 then nlevels   = 15
    if vx_name ne '' xor vy_name ne '' then $
        message, 'VX_NAME and VY_NAME must be used together.'
    
    nVLines = n_elements(vert_lines)
    nHLines = n_elements(horiz_lines)
    
    ;Create a simulation object
    if n_elements(oSim) gt 0 then begin
        if obj_valid(oSim) eq 0 || obj_isa(oSim, 'MRSIM') eq 0 then $
            message, 'SIM_OBJECT must be valid and a subclass of "MrSim"'
        sim_class = obj_class(oSim)
    endif else begin
        Sim3D   = keyword_set(Sim3D)
        if Sim3D then sim_class = 'MRSIM3D' else sim_class = 'MRSIM2D'
        oSim = obj_new(sim_class, time, _STRICT_EXTRA=extra)
    endelse
    
    ;Buffer the output?
    if current eq 0 then $
        if ofilename eq '' then buffer = 0 else buffer = 1
    
;-------------------------------------------------------
;Read Data /////////////////////////////////////////////
;-------------------------------------------------------
    ;Create the object and get the data. Use Ay for contours.
    data = oSim -> getData(name)
    if c_name ne '' then c_data = oSim -> GetData(c_name)
    if vx_name ne '' and vy_name ne '' then begin
        vx_data = oSim -> GetData(vx_name)
        vy_data = oSim -> GetData(vy_name)
    endif

    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, XSIM=XSim, ZSIM=ZSim, AXIS_LABELS=axLabls, $
                         COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
        
    ;Destroy the object
    if arg_present(oSim) eq 0 then obj_destroy, oSim
    
    ;Make sure the data exists. Do this /after/ OSIM has a chance to be destroyed
    ;and before the graphic window is created.
    if data eq !Null then return, obj_new()

;-------------------------------------------------------
;Prepare Ranges ////////////////////////////////////////
;-------------------------------------------------------
    ;Rename items
    _name = MrSim_Rename(name, coord_system, MVA_FRAME=mva_frame, /SUBSCRIPT)
    units = MrSim_Rename(units, /SUBSCRIPT)
    
    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't$\Omega$$\downci$=' + string(time*dtxwci, FORMAT='(f0.1)') $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')

    if sim_class eq 'MRSIM3D' then begin
        oSim -> GetProperty, YSLICE=yslice
        title += '  ' + axLabls[1] + '=' + string(yslice, FORMAT='(f0.1)') + units
    endif

;-------------------------------------------------------
;Make Plots ////////////////////////////////////////////
;-------------------------------------------------------

    ;Create a window
    if current $
        then colorWin = GetMrWindows(/CURRENT) $
        else colorWin = MrWindow(XSIZE=600, YSIZE=300, BUFFER=buffer, OXMARGIN=[10,14], REFRESH=0)

    ;Create a Color Plot of the named variable
    colorIm = MrImage(data, XSim, ZSim, $
                      TITLE=title, RANGE=range, $
                      XTITLE=axLabls[0] + ' (' + units + ')', XRANGE=[XSim[0],XSim[-1]], XSTYLE=1, $
                      YTITLE=axLabls[2] + ' (' + units + ')', YRANGE=[ZSim[0],ZSim[-1]], YSTYLE=1, $
                      /AXES, /SCALE, CTINDEX=13, NAME='Color ' + name, /CURRENT)
                                
    colorCB = MrColorbar(TARGET=colorIm, NAME='CB: Color ' + name, TITLE=_name)

    ;Bind the colorbar to the image.
    colorWin -> Bind, [colorCB, colorIm], /CAXIS
    
;-------------------------------------------------------
;Add Contours? /////////////////////////////////////////
;-------------------------------------------------------
    if c_name ne '' then begin
        ;Create the contour plot
        levels = cgConLevels(c_data, NLEVELS=nLevels)
        colorCo = MrContour(c_data, XSim, ZSim, $
                            XRANGE=[XSim[0],XSim[-1]], $
                            YRANGE=[ZSim[0],ZSim[-1]], $
                            LEVELS=levels, NAME=c_name + ' Contours', C_LABELS=0, $
                            OVERPLOT=colorIm, /CURRENT)
                            
        ;Bind the contours to the image
        colorWin -> Bind, [colorIm, colorCo], /XAXIS
    endif
    
;-------------------------------------------------------
;Add Vector Field? /////////////////////////////////////
;-------------------------------------------------------
    if vx_name ne '' then begin
        colorVec = MrVector(vx_data, vy_data, XSim, ZSim, COLOR='Black', $
                            XRANGE=[XSim[0],XSim[-1]], FRACTION=fraction, $
                            YRANGE=[ZSim[0],ZSim[-1]], /ORDERED, $
                            NAME=vx_name + '-' + vy_name + ' Vectors', $
                            OVERPLOT=colorIm, /CURRENT)
    endif
                            
;-------------------------------------------------------
;Add Lines? ////////////////////////////////////////////
;-------------------------------------------------------
    nLineColors = n_elements(line_color)
    if nLineColors eq 0 then begin
        if nVLines gt 1 || nHLines gt 1 $
            then line_color = ['Blue', 'Forest Green', 'Red', 'Magenta', 'Orange'] $
            else line_color = 'Black'
        nLineColors = n_elements(line_color)
    endif

    ;Vertical Lines
    if nVLines gt 0 then begin
        for i = 0, nVLines - 1 do begin
            vlPlot = MrPlot([vert_lines[i], vert_lines[i]], [ZSim[0],ZSim[-1]], $
                            /CURRENT, OVERPLOT=colorIm, COLOR=line_color[i mod nLineColors], $
                            NAME='VLine ' + name + ' ' + strtrim(i, 2))
        endfor
    endif

    ;Horizontal Lines
    if nHLines gt 0 then begin
        hlPlot = MrPlot([XSim[0],XSim[-1]], rebin(reform(horiz_lines, 1, nHLines), 2, nHLines), $
                        /CURRENT, OVERPLOT=colorIm, COLOR=line_color[0:nHLines-1], $
                        DIMENSION=1, NAME='HLines ' + name)
    endif
;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------
    
    ;Refresh and output, if requested.
    if current eq 0 then begin
        colorWin -> Refresh
        if ofilename ne '' then colorWin -> Save, ofilename
    endif

    return, colorWin
end
