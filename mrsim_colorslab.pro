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
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                               valid simulation object. See MrSim_Create.pro.
;       NAME:               in, required, type=string
;                           The name of the vector quantity to be plotted in color.
;       TIME:               in, required, type=int
;                           The simulation time for which a velocity vector field
;                               is to be plotted. Required if `THESIM` is not an
;                               object.
; :Keywords:
;       C_NAME:             in, optional, type=string, default=''
;                           Name of the data product whose contours are to be 
;                               overplotted on the image.
;       CURRENT:            in, optional, type=boolean, default=0
;                           If set, the graphics will be added to the current MrWindow
;                               graphics window. The default is to create a new window.
;       FRACTION:           in, optional, type=float
;                           Fraction of the vectors to be drawn. Used with `VX_NAME`
;                               and `VY_NAME`.
;       HORIZ_LINES:        in, optional, type=flt/fltarr
;                           Locations, in data coordinates, along the y-axis of the
;                               image, at which horizontal lines are to be drawn.
;       LINE_COLOR:         in, optional, type=strarr, default="['Blue', 'Forest Green', 'Red', 'Magenta', 'Orange']"
;                           Colors of the lines to be drawn on the image.
;       NLEVELS:            in, optional, type=int, default=15
;                           The number of contour lines to draw. Used with `C_NAME`.
;       OFILENAME:          out, optional, type=string, default=''
;                           If provided, graphics will be output to a file by this
;                               name. The file extension is used to determine which
;                               type of image to create. "PNG" is the default.
;       RANGE:              in, optional, type=fltarr(2), default=[min\,max]
;                           Data range of the quantity to be displayed.
;       SIM_OBJECT:         out, optional, type=object
;                           If `THESIM` is the name or number of a simulation, then
;                               this keyword returns the object reference to the
;                               corresponding simulation object that is created.
;       VX_NAME:            in, optional, type=string, default=''
;                           Name of the x-component of a data product whose vector
;                               field is to be overplotted on the image. Must be used
;                               with `VY_NAME`.
;       VY_NAME:            in, optional, type=string, default=''
;                           Name of the y-component of a data product whose vector
;                               field is to be overplotted on the image. Must be used
;                               with `VX_NAME`.
;       VERT_LINES:         in, optional, type=flt/fltarr
;                           Locations, in data coordinates, along the x-axis of the
;                               image, at which vertical lines are to be drawn.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
;
; :Returns:
;       COLORWIN:           MrWindow graphic window containing the requested graphics.
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
;       2014/09/30  -   Removed the SIM3D keyword and added the THESIM parameter.
;                           Repurposed the SIM_OBJECT keyword. - MRA
;       2014/10/02  -   If C_NAME='Ay', add a contour passing through the X-line. - MRA
;-
function MrSim_ColorSlab, theSim, name, time, yslice, $
C_NAME = c_name, $
CURRENT = current, $
FRACTION = fraction, $
HORIZ_LINES = horiz_lines, $
LINE_COLOR = line_color, $
NLEVELS = nLevels, $
OFILENAME = ofilename, $
RANGE = range, $
SIM_OBJECT = oSim, $
VX_NAME = vx_name, $
VY_NAME = vy_name, $
VERT_LINES = vert_lines, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if osim_created        && arg_present(oSim)    eq 0 then obj_destroy, oSim
        if obj_valid(colorWin) && keyword_set(current) eq 0 then obj_destroy, colorWin
        void = cgErrorMSG()
        return, obj_new()
    endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
    osim_created = 0B
    
    ;Simulation name or number?
    if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
        oSim = MrSim_Create(theSim, time, yslice, _STRICT_EXTRA=extra)
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
    oSim -> GetProperty, TIME=time, XSIM=XSim, YSIM=YSim, ZSIM=ZSim, AXIS_LABELS=axLabls, $
                         COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame, $
                         ORIENTATION=orientation
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    
    ;Make sure the data exists. Do this /after/ OSIM has a chance to be destroyed
    ;and before the graphic window is created.
    if data eq !Null then return, obj_new()

;-------------------------------------------------------
;Prepare Ranges ////////////////////////////////////////
;-------------------------------------------------------
    ;Rename items
    _name = MrSim_Rename(name, coord_system, MVA_FRAME=mva_frame, /SUBSCRIPT)
    units = MrSim_Rename(units, /SUBSCRIPT)
    
    ;Time interval in wci given?
    ;   - Display time in terms of t*wci.
    ;   - Display the time index.
    if n_elements(dtxwci) gt 0 $
        then title = 't$\Omega$$\downci$=' + string(time*dtxwci, FORMAT='(f0.1)') $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')

    ;If a 3D simulation was given, several orientations are possible
    case orientation of
        'XY': begin
            x = temporary(xSim)
            y = temporary(ySim)
            xtitle  = axLabls[0] + ' (' + units + ')'
            ytitle  = axLabls[1] + ' (' + units + ')'
            
            oSim -> GetProperty, ZRANGE=zrange
            title  += '  ' + axLabls[2] + '=' + string(zrange[0], FORMAT='(f0.1)') + units
        endcase
    
        'XZ': begin
            x = temporary(xSim)
            y = temporary(zSim)
            xtitle = axLabls[0] + ' (' + units + ')'
            ytitle = axLabls[2] + ' (' + units + ')'
            
            if sim_class eq 'MRSIM3D' then begin
                oSim -> GetProperty, YRANGE=yrange
                title += '  ' + axLabls[1] + '=' + string(yrange[0], FORMAT='(f0.1)') + units
            endif
        endcase
    
        'YZ': begin
            x = temporary(ySim)
            y = temporary(zSim)
            xtitle = axLabls[1] + ' (' + units + ')'
            ytitle = axLabls[2] + ' (' + units + ')'
        endcase
    
        else: message, 'Orientation unknown: "' + orientation + '".'
    endcase
        
    ;Destroy the object
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim

;-------------------------------------------------------
;Make Plots ////////////////////////////////////////////
;-------------------------------------------------------

    ;Create a window
    if current $
        then colorWin = GetMrWindows(/CURRENT) $
        else colorWin = MrWindow(XSIZE=600, YSIZE=300, BUFFER=buffer, OXMARGIN=[10,14], REFRESH=0)

    ;Create a Color Plot of the named variable
    colorIm = MrImage(data, x, y, $
                      /AXES, /CURRENT, /SCALE, $
                      CTINDEX = 13, $
                      NAME    = 'Color ' + name, $
                      RANGE   = range, $
                      TITLE   = title, $
                      XTITLE  = xtitle, $
                      YTITLE  = ytitle)

    ;Create a color bar
    colorCB = MrColorbar(TARGET=colorIm, NAME='CB: Color ' + name, TITLE=_name)

    ;Bind the colorbar to the image.
    colorWin -> Bind, [colorCB, colorIm], /CAXIS
    
;-------------------------------------------------------
;Add Contours? /////////////////////////////////////////
;-------------------------------------------------------
    if c_name ne '' then begin
        ;The X-point is a saddle point in Ay
        ;   - Find the minimum along the bottom boundary.
        ;   - Find the maximum along a vertical cut that passes through min point.
        if strupcase(c_name) eq 'AY' then begin
            nx = n_elements(x)
            ny = n_elements(y)

            ;Take the derivative
            dx = c_data[1:nx-1, *] - c_data[0:nx-2, *]
            dz = c_data[*, 1:nz-1] - c_data[*, 0:nz-2]
            
            ;Find the minimum of Ay along each row
            ;   - Turn the 1D indices into 2D indices
            xMin = min(dx, ixMin, DIMENSION=1, /ABSOLUTE)
            inds = array_indices([nx-1, nz], ixMin, /DIMENSIONS)
            
            ;Find the maximum along the path of minimums just found.
            zMin = min(dz[inds[0,*], inds[1,*]], izMin, DIMENSION=2, /ABSOLUTE)
            sepAy  = c_data[inds[0, izMin], inds[1, izMin]]

            ;Create the contour levels
            levels = [sepAy, cgConLevels(c_data, NLEVELS=nLevels)]
            levels = levels[sort(levels)]
        endif else begin
            levels = cgConLevels(c_data, NLEVELS=nLevels)
        endelse
        
        ;Create the contour
        colorCo = MrContour(c_data, x, y, $
                            /CURRENT, $
                            C_LABELS = 0, $
                            LEVELS   = levels, $
                            NAME     = c_name + ' Contours', $
                            OVERPLOT = colorIm)
                            
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
