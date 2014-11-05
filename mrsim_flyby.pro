; docformat = 'rst'
;
; NAME:
;    MrSim_FlyBy
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
;       TIME:               in, required, type=int
;                           The simulation time at which cuts are to be taken.
;
; :Keywords:
;       CURRENT:            in, optional, type=boolean, default=0
;                           If set, the graphics will be added to the current MrWindow
;                               graphics window. The default is to create a new window.
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
;       _REF_EXTRA:         in, optional, type=structure
;                           Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                               is an object.
;
; :Returns:
;       fbWin:              A MrWindow graphic window containing the line cut. If the
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
;       2014/10/18  -   Written by Matthew Argall
;-
function MrSim_FlyBy, theSim, name, r0, r1, time, $
ADD_LEGEND = add_legend, $
C_NAME = c_name, $
CURRENT = current, $
COLOR = color, $
DIST_LAYOUT = dist_layout, $
DIST_SIZE = dist_size, $
DIST_TYPE = dist_type, $
DIST_WIN = distWin, $
IM_NAME = im_name, $
NDIST = nDist, $
OFILENAME = ofilename, $
OVERPLOT = overplot, $
SIM_OBJECT = oSim, $
VELOCITY   = velocity, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
        if current eq 0 && obj_valid(fbWin) then obj_destroy, fbWin
        if obj_valid(distWin) then obj_destroy, distWin
        void = cgErrorMSG()
        return, obj_new()
    endif

;-------------------------------------------------------
; Simulation Object ////////////////////////////////////
;-------------------------------------------------------
    osim_created = 0B
    current      = keyword_set(current)
    
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
; Defaults ////////////////////////////////////////////
;-------------------------------------------------------

    ;Set defaults
    add_legend = keyword_set(add_legend)
    if n_elements(im_name)   eq 0 then im_name   = ''
    if n_elements(nDist)     eq 0 then nDist     = 0
    if n_elements(ofilename) eq 0 then ofilename = ''
    
;-------------------------------------------------------
; Data, Metadata, Window ///////////////////////////////
;-------------------------------------------------------

    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, ORIENTATION=orientation, YRANGE=yrange, $
                         COORD_SYSTEM=coord_system, MVA_FRAME=mva_frame
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    
    ;Name
    ;   - Apply subscripts and labels
    _name = MrSim_Rename(name, coord_system, MVA_FRAME=mva_frame, /SUBSCRIPT)
    units = MrSim_Rename(units, /SUBSCRIPT)

    ;Time
    ;   - Inverse gyro-time or time index?
    if n_elements(dtxwci) gt 0 $
        then title = 't$\Omega$$\downci$$\up-1$=' + string(time*dtxwci, FORMAT='(f0.1)') $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')

    ;Graphics window
    fbWin = MrWindow(NAME='FlyBy', XSIZE=500, YSIZE=500, REFRESH=0)

;-------------------------------------------------------
; 2D Overview //////////////////////////////////////////
;-------------------------------------------------------
    if im_name ne '' then begin
        ;Make room for the colorbar
        fbWin.OXMARGIN = [10, 15]

        ;Overview plot
        !Null   = MrSim_ColorSlab(oSim, im_name, /CURRENT, C_NAME=c_name)
        colorIm = fbWin['Color ' + im_name]
    endif

;-------------------------------------------------------
; FlyBy ///////////////////////////////////////////////
;-------------------------------------------------------
    ;Get data.
    cells = oSim -> GridLine(r0, r1, COORDS=coords)
    data  = oSim -> ReadGDA_Cells(cells, name)
    
    ;Clean up
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    if data eq !Null then message, 'NAME must be the name of a GDA file.'

    ;Distance along the line
    nPts           = n_elements(coords[0,*])
    dl             = fltarr(3, nPts)
    dl[0,1:nPts-1] = coords[0,1:nPts-1] - coords[0,0]
    dl[1,1:nPts-1] = coords[1,1:nPts-1] - coords[1,0]
    dl[2,1:nPts-1] = coords[2,1:nPts-1] - coords[2,0]
    dl             = sqrt(total(dl^2, 1))

    ;Convert distance to time?
    if n_elements(velocity) gt 0 then begin
        x      = temporary(dl) / velocity
        xtitle = 'Time (s)'
    endif else begin
        x = temporary(dl)
        xtitle = 'Distance (' + units + ')'
    endelse

    ;Create a Line Plot of the named variable
    lineCut = MrPlot(x, data, $
                     /CURRENT, $
                     COLOR     = color, $
                     DIMENSION = dimension, $
                     OVERPLOT  = overplot, $
                     NAME      = 'Cut ' + name, $
                     TITLE     = title, $
                     XTITLE    = xtitle, $
                     YTITLE    = _name)
    
    ;Draw a line where the cut was taken
    if obj_valid(colorIm) then begin
        path = MrPlotS([r0[0], r1[0]], [r0[2], r1[2]], $
                       /DATA, $
                       COLOR  = 'White', $
                       NAME   = 'Trajectory', $
                       TARGET = colorIm)
    endif

;-------------------------------------------------------
; Make Distributions ///////////////////////////////////
;-------------------------------------------------------
    if nDist gt 0 then begin
        ;Default distribution types
        if n_elements(dist_layout) eq 0 then begin
            nCols = ceil(sqrt(nDist))
            nRows = ceil(float(nDist) / nCols)
            layout = [nCols, nRows]
        endif
        if n_elements(dist_size) eq 0 then dist_size = [1,1,1]
        if n_elements(dist_type) eq 0 then dist_type = 'Vpar_Vperp'
        
        ;Make everything 3D
        p0 = n_elements(r0)        eq 3 ? r0        : [r0[0],        0, r0[1]]
        p1 = n_elements(r1)        eq 3 ? r1        : [r1[0],        0, r1[1]]
        sz = n_elements(dist_size) eq 3 ? dist_size : [dist_size[0], 0, dist_size[1]]
        
        ;Locations of distributions
        bin_centers      = fltarr(3, nDist)
        half_width       = fltarr(3, nDist)
        bin_centers[0,*] = linspace(p0[0], p1[0], nDist)
        bin_centers[1,*] = make_array(nDist, TYPE=6, VALUE=p0[1])
        bin_centers[2,*] = linspace(p0[2], p1[2], nDist)
        half_width[0,*]  = sz[0]
        half_width[1,*]  = sz[1]
        half_width[2,*]  = sz[2]

        ;Create the distributions
        distWin = MrSim_eMap(oSim, dist_type, bin_centers, half_width, $
                             LAYOUT    = dist_layout, $
                             POSITIONS = positions)

        ;Draw the location of the distributions on the overview plot
        if obj_valid(colorIm) then begin
            ;Draw boxes where the distribution functions are being taken.
            for i = 0, nDist - 1 do begin
                theBin = positions[*,i]
                xpoly  = theBin[[0,3,3,0,0]]
                ypoly  = theBin[[2,2,5,5,2]]
                !Null  = MrPlotS(xpoly, ypoly, $
                                 /DATA, $
                                 COLOR  = 'White', $
                                 NAME   = 'Bin ' + strtrim(i, 2), $
                                 TARGET = colorIm, $
                                 THICK  = 2.0)
            endfor
        endif
    endif
        
;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------
    ;Add a legend if more than one cut was taken
    if add_legend then begin
        !Null = MrLegend(TARGET   = lineCut, $
                         /BOX, $
                         BX_COLOR = 'White', $
                         LENGTH   = 0, $
                         LOCATION = 7, $
                         COLOR    = color, $
                         TITLE    = [string(FORMAT='(%"r0 = (%0.1f, %0.1f, %0.1f)")', r0), $
                                     string(FORMAT='(%"r1 = (%0.1f, %0.1f, %0.1f)")', r1)])
    endif
    
    ;Refresh and output, if requested.
    if ofilename ne '' then lineCut -> Output, ofilename
    
    fbWin -> Refresh
    return, fbWin
end
