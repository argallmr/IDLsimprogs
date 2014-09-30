; docformat = 'rst'
;
; NAME:
;    MrSim_LineCut
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
;   The purpose of this program is to create a line plot at a particular location. Data
;   is taken from a 1D cut either horizontally or vertically through a 2D image
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
;       CUTS:               in, required, type=long/lonarr
;                           Locations at which cuts are to be made. If an array is
;                               provided, the trace of each cut will be overplotted into
;                               the same set of axes.
;       TIME:               in, required, type=int
;                           The simulation time at which cuts are to be taken.
;
; :Keywords:
;       CURRENT:            in, optional, type=boolean, default=0
;                           If set, the graphics will be added to the current MrWindow
;                               graphics window. The default is to create a new window.
;       HCUT_RANGE:         in, optional, type=fltarr(2)
;                           The horizontal range, in data coordinates, over which to cut.
;                               "Horizontal" is defined by the `HORIZONTAL` keyword. Data
;                               coordinates are determined by the "ion_scale" property
;                               within `SIM_OBJECT`. The default is to take the
;                               appropriate simulation range.
;       HORIZONTAL:         in, optional, type=boolean, default=0
;                           If set, a horizontal cut will be taken. The default is
;                               to take a vertical cut. For an "XY" orientation, "X" is
;                               horizontal, "Y" is vertical. Similar for "XZ", etc.
;                               The orientationis taken from the `SIM_OBJECT` property.
;       OFILENAME:          out, optional, type=string, default=''
;                           If provided, graphics will be output to a file by this
;                               name. The file extension is used to determine which
;                               type of image to create. "PNG" is the default.
;       OVERPLOT:           in, optional, type=boolean/object
;                           Set either to 1 (one) or a MrGraphics object. In the former
;                               case, line cuts will be place over the selected or highest
;                               ordered graphic. In the latter case, the MrGraphics object
;                               on which to overplot the line cuts.
;       SIM_OBJECT:         out, optional, type=object
;                           If `THESIM` is the name or number of a simulation, then
;                               this keyword returns the object reference to the
;                               corresponding simulation object that is created.
;       VCUT_RANGE:         in, optional, type=fltarr(2)
;                           The verical range, in data coordinates, over which to cut.
;                               "Vertical" is defined by the `HORIZONTAL` keyword. Data
;                               coordinates are determined by the "ion_scale" property
;                               within `SIM_OBJECT`. The default is to take the
;                               appropriate simulation range.
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
;
; :Returns:
;       LCWIN:              A MrWindow graphic window containing the line cut. If the
;                               requested data product does not exist, a null object will
;                               be returned.
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       GetMrWindows.pro
;       MrPlot.pro
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
;       2014/02/07  -   Written by Matthew Argall
;       2014/03/05  -   The parameter CUT was renamed to CUTS and can now be an array. - MRA
;       2014/03/28  -   Added the COLOR and OVERPLOT keyword. Return the plot's object
;                           reference, not the window's. - MRA
;       2014/04/16  -   If one cut is taken, add the location to the title. If more than
;                           one, add a legend indicating where the cuts were taken. - MRA
;       2014/09/30  -   Removed the SIM3D keyword and added the THESIM parameter.
;                           Repurposed the SIM_OBJECT keyword. - MRA
;-
function MrSim_LineCut, theSim, name, cuts, time, $
ADD_LEGEND = add_legend, $
CURRENT = current, $
COLOR = color, $
HCUT_RANGE = hcut_range, $
HORIZONTAL = horizontal, $
OFILENAME = ofilename, $
OVERPLOT = overplot, $
SIM_OBJECT = oSim, $
VCUT_RANGE = vcut_range, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
        if current eq 0 && obj_valid(lcWin) then obj_destroy, lcWin
        return, obj_new()
    endif

;-------------------------------------------------------
; Check Simulat/////////////////////////////////////////
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
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------

    ;Set defaults
    add_legend = keyword_set(add_legend)
    current    = keyword_set(current)
    horizontal = keyword_set(horizontal)
    Sim3D      = keyword_set(Sim3D)
    if n_elements(ofilename) eq 0 then ofilename = ''
    nCuts = n_elements(cuts)
    
;-------------------------------------------------------
;Read Data /////////////////////////////////////////////
;-------------------------------------------------------
    ;Create the object and get the data. Use Ay for contours.
    data = oSim -> LineCuts(name, cuts, pos, HORIZONTAL=horizontal, $
                            HCUT_RANGE=hcut_range, VCUT_RANGE=vcut_range)

    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, ORIENTATION=orientation, YSLICE=yslice, $
                         COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
        
    ;Destroy the object
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    
    ;Make sure the data exists. Do this /after/ OSIM has a chance to be destroyed
    ;and before the graphic window is created.
    if data eq !Null then return, obj_new()
    
    ;Rename the data products
    _name = MrSim_Rename(name, coord_system, MVA_FRAME=mva_frame, /SUBSCRIPT)
    units = MrSim_Rename(units, /SUBSCRIPT)

;-------------------------------------------------------
;Prepare Ranges ////////////////////////////////////////
;-------------------------------------------------------
    xaxis = strlowcase(strmid(orientation, 0, 1))
    yaxis = strlowcase(strmid(orientation, 1, 1))
    
    ;Dimension to plot
    dimension = nCuts eq 1 ? 0 : 1

    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't$\Omega$$\downci$$\up-1$=' + string(time*dtxwci, FORMAT='(f0.1)') $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')

    ;Where was the cut taken?
    if horizontal eq 1 then begin
        xaxis = strlowcase(strmid(orientation, 0, 1))
        caxis = strlowcase(strmid(orientation, 1, 1))
        pname = 'Cut ' + name
    endif else begin
        xaxis = strlowcase(strmid(orientation, 1, 1))
        caxis = strlowcase(strmid(orientation, 0, 1))
        pname = 'Cut ' + name
    endelse
    
    ;Add the location to the title
    if n_elements(cuts) eq 1 then begin
        if sim_class eq 'MRSIM2D' $
            then title += '  '  + caxis +   '='  + string(cuts, FORMAT='(f0.1)')  + units $
            else title += '  (' + caxis + 'y)=(' + strjoin(string(cuts, yslice, FORMAT='(f0.1)'), ',') + ')' + units
    endif else begin
        if sim_class eq 'SIM3D' $
            then title += '  ' + caxis + 'y=' + string(yslice, FORMAT='(f0.1)') + units
    endelse

;-------------------------------------------------------
;Make Plots ////////////////////////////////////////////
;-------------------------------------------------------
    if horizontal $
        then xrange = hcut_range $
        else xrange = vcut_range

    ;Create a Line Plot of the named variable
    lineCut = MrPlot(pos, data, CURRENT=current, $
                     COLOR=color, $
                     DIMENSION=dimension, $
                     OVERPLOT=overplot, $
                     NAME=pname, $
                     TITLE=title, $
                     XTITLE=xaxis + ' (' + units + ')', $
                     XRANGE=xrange, $
                     YRANGE=yrange, $
                     YTITLE=_name)
    
;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------
    ;Add a legend if more than one cut was taken
    if add_legend && nCuts gt 1 then begin
        lcLegend = MrLegend(TARGET=lineCut, LOCATION=7, COLOR=color, $
                            TITLE=string(cuts, FORMAT='(f0.1)'), LENGTH=0)
    endif
    
    ;Refresh and output, if requested.
    if ofilename ne '' then lineCut -> Output, ofilename
    
    return, lineCut
end
