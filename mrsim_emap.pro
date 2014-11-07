; docformat = 'rst'
;
; NAME:
;    MrSim_eMap
;
; PURPOSE:
;+
;   The purpose of this program is to create an array of electron distributions, a.k.a.
;   an "eMap".
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       THESIM:         in, required, type=string/integer/object
;                       The name or number of the simulation to be used, or a
;                           valid simulation object. See MrSim_Create.pro.
;       TYPE:           in, optional, type=string, default='Vx-Vz'
;                       The type of distribution function to make. Choices are::
;                           'x-z'
;                           'x-Vx'
;                           'z-Vz'
;                           'Vx-Vz'
;                           'Vx-Vy'
;                           'Vy-Vz'
;                           'Vpar-Vperp'
;                           'Vpar-Vperp1'
;                           'Vpar-Vperp2'
;                           'Vperp1-Vperp2'
;       X:              in, required, type=float/fltarr(N)
;                       X-location(s) of the center of the spacial bin from which electron
;                           distribution information is collected. If both `X` and `Z` are
;                           vectors, they must be the same length.
;       Z:              in, required, type=float/fltarr(N)
;                       Z-location(s) of the center of the spacial bin from which electron
;                           distribution information is collected. If both `X` and `Z` are
;                           vectors, they must be the same length.
;       DX:             in, optional, type=float/fltarr(N), default=1
;                       Half-width of the bin(s).
;       DZ:             in, optional, type=float/fltarr(N), default=1
;                       Half-height of the bin(s).
;
; :Keywords:
;       C_NAME:         in, optional, type=string, default=''
;                       Name of the data product whose contours are to be 
;                           overplotted on the image.
;       CIRCLES:        in, optional, type=boolean, default=0
;                       If set, concentric cirlces will be drawn at v=0.25 and v=0.5.
;       HGAP:           in, optional, type=float, default=0.0
;                       Horizontal space between adjacent distributions. Ignored of either
;                           `X` or `Z` are vectors.
;       IM_NAME:        in, optional, type=string, default=''
;                       Name of a data product for which a 2D overview plot is to be
;                           made. White square boxes will be drawn on the image to
;                           indicate from where the distributions were taken.
;       IM_WIN:         out, optional, type=object
;                       If `IM_NAME` is provided, then use this keyword to obtain
;                           the object reference of the window in which the image
;                           was created.
;       NBINS:          in, optional, type=long/lonarr(2), default=75
;                       Number of bins into which a distribution will be subdivided.
;                           If a scalar is given, the horizontal and vertical
;                           components will be split into the same number of bins.
;                           Provide a 2-element vector to have unequal number of bins.
;       LAYOUT:         in, optional, type=intarr(2), default=[1,1]
;                       Specifies the number of columns and rows: [nCols, nRows]; in
;                           the eMap. There will be nCols*nRows number of distributions
;                           total. If `X` and `Z` are scalars, a regular grid will be
;                           formed. If either or both are vectors, nCols*nRows must be
;                           at least as large as the length of the vector. The layout is
;                           filled from left-to-right then top-to-bottom.
;       LOCATION:       in, optional, type=integer, default=5
;                       Specifies which bin within the eMap is centered on `X` and `Y`.
;                           Options are::
;                                   1: Upper-Left
;                                   2: Upper-Middle
;                                   3: Upper-Right
;                                   4: Center-Left
;                                   5: Center-Middle
;                                   6: Center-Right
;                                   7: Bottom-Left
;                                   8: Bottom-Middle
;                                   9: Bottom-Right
;       POSITIONS:      out, optional, type=4xN float
;                       A named variable to receive the positions of each bin, in data
;                           coordinates. "Position" has the IDL definition.
;       SIM_OBJECT:     out, optional, type=object
;                       If `THESIM` is the name or number of a simulation, then
;                           this keyword returns the object reference to the
;                           corresponding simulation object that is created.
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, velocity will be normalized to v_A instead of c::
;                               c / vA = sqrt(mi_me) * f_pe / f_ce.
;                           so multiplying v/c by the above equation does the trick.
;       VGAP:           in, optional, type=float, default=0.0
;                       Vertical space between adjacent distributions. Ignored of either
;                           `X` or `Z` are vectors.
;       _REF_EXTRA:     in, optional, type=structure
;                       Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                           is an object.
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       MrSim_LineCut.pro
;       MrSim_ColorSlab.pro
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
;       2014/09/01  -   Written by Matthew Argall
;       2014/09/05  -   Added the IM_WIN keyword. - MRA
;       2014/09/30  -   Removed the SIM3D keyword and added the THESIM parameter.
;                           Repurposed the SIM_OBJECT keyword. Added the CIRCLES and V_VA
;                           keyword. - MRA
;       2014/10/13  -   Added the XSIZE and YSIZE keywords. X, Y, DX, DY can be vectors. - MRA.
;       2014/10/18  -   Added the POSITIONS keyword. - MRA
;       2014/10/30  -   Changed input parameters X, Z, DX, and, DZ to BIN_CENTER and
;                           HALF_WIDTH. - MRA
;-
function MrSim_eMap, theSim, type, bin_center, half_width, $
C_NAME=c_name, $
CIRCLES=circles, $
HGAP=hGap, $
IM_NAME=im_name, $
IM_WIN=im_win, $
NBINS=nBins, $
LAYOUT=layout, $
LOCATION=location, $
POSITIONS=positions, $
SIM_OBJECT=oSim, $
VGAP=vGap, $
V_VA=v_va, $
XSIZE=xsize, $
YSIZE=ysize, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
        if obj_valid(im_win)  then obj_destroy, im_win
        if obj_valid(win2)    then obj_destroy, win2
        if obj_valid(imgTemp) then imgTemp -> Close
        void = cgErrorMSG()
        return, obj_new()
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
    nDims = size(bin_center, /N_DIMENSIONS)
    dims  = size(bin_center, /DIMENSIONS)
    nDist = nDims eq 1 ? dims[0] : dims[1]

    if n_elements(nBins)    eq 0 then nBins    = 75
    if n_elements(layout)   eq 0 then layout   = [1,nDist]
    if n_elements(location) eq 0 then location = 5
    if n_elements(c_name)   eq 0 then c_name   = ''
    if n_elements(im_name)  eq 0 then im_name  = ''
    if n_elements(vGap)     eq 0 then vGap     = 0
    if n_elements(hGap)     eq 0 then hGap     = 0
    if n_elements(xsize)    eq 0 then xsize    = 800
    if n_elements(ysize)    eq 0 then ysize    = 600
    
    ;Make sure the distribution is large enough
    if nDist gt 1 then begin
        if nDist gt layout[0]*layout[1] then $
            message, string(FORMAT='(%"LAYOUT is not large enough to hold %i distributions.")', nDist)
    endif

;-------------------------------------------------------
; Determine Location of Dist Fns ///////////////////////
;-------------------------------------------------------
    if nDims eq 1 then begin
    ;-------------------------------------------------------
    ; Create Array of Distribution Functions ///////////////
    ;-------------------------------------------------------
        ;Center of eMap and half-widths of distributions.
        if dims[0] eq 2 then begin
            center = [bin_center[0], 0, bin_center[1]]
            width  = [half_width[0], 0, half_width[1]]
        endif else begin
            center = bin_center
            width  = half_width
        endelse
        
        ;Position of eMap bin indicated by LOCATION
        ;   - [x0, y0, z0, x1, y1, z1]
        positions = [center[0] - width[0], $
                     center[1] - width[1], $
                     center[2] - width[2], $
                     center[0] + width[0], $
                     center[1] + width[1], $
                     center[2] + width[2]]
        
        ;Bin-center offsets from the center of the bin defined by LOCATION
        ;   - Dimensions are [position, col, row]
        ;   - Value is +/-N
        ;       o (+/-) indicates left/right (up/down) offset
        ;       o Offset by N bins
        hoffset = indgen(1, layout[0])
        voffset = indgen(1, 1, layout[1])

        ;XOFFSET
        if layout[0] gt 1 then begin
            case location of
                1: ;Do nothing                      ;LEFT
                2: hoffset -= layout[0] / 2         ;CENTER
                3: hoffset -= layout[0] - 1         ;RIGHT
                4: ;Do nothing                      ;LEFT
                5: hoffset -= layout[0] / 2         ;CENTER
                6: hoffset -= layout[0] - 1         ;RIGHT
                7: ;Do nothing                      ;LEFT
                8: hoffset -= layout[0] / 2         ;CENTER
                9: hoffset -= layout[0] - 1         ;RIGHT
                else: message, 'LOCATION must be an integer between 1 and 9.'
            endcase
        endif
    
        ;If layout[2] = 2
        ;   - VOFFSET will have only one dimension.
        ;   - The call to Reverse below will fail.
        if layout[1] gt 1 then begin
            case location of
                1: voffset  = -voffset                              ;TOP
                2: voffset  = -voffset                              ;TOP
                3: voffset  = -voffset                              ;TOP
                4: voffset  = reverse(voffset - layout[1]/2, 3)     ;MIDDLE
                5: voffset  = reverse(voffset - layout[1]/2, 3)     ;MIDDLE
                6: voffset  = reverse(voffset - layout[1]/2, 3)     ;MIDDLE
                7: ;Do nothing                                      ;BOTTOM
                8: ;Do nothing                                      ;BOTTOM
                9: ;Do nothing                                      ;BOTTOM
                else: message, 'LOCATION must be an integer between 1 and 9.'
            endcase
        endif

        ;Include bin widths and gaps between bins
        hoffset = (2*width[0] + hGap) * rebin(hoffset, 2, layout[0], layout[1])
        voffset = (2*width[2] + vGap) * rebin(voffset, 2, layout[0], layout[1])
        
        ;Data locations of distribution functions
        ;   - Saved as [position, index]
        ;   - Position is the IDL definition of the word
        ;   - Index starts with 0 in the upper-right corner, then moves right then down
        nDist                 = product(layout)
        positions             = rebin(positions, 6, layout[0], layout[1])
        positions[[0,3],*,*] += hoffset
        positions[[2,5],*,*] += voffset
        positions             = reform(positions, 6, nDist)
        hoffset               = reform(hoffset,   2, nDist)
        voffset               = reform(voffset,   2, nDist)

        ;Bin centers
        centers = fltarr(3, nDist)
        centers[0,*] = center[0] + (temporary(hoffset))[0,*]
        centers[1,*] = center[1]
        centers[2,*] = center[2] + (temporary(voffset))[0,*]
        center = !Null
        
        ;Half widths
        widths = rebin(temporary(width), 3, nDist)
        
    endif else begin
    ;-------------------------------------------------------
    ; Individual Locations Given ///////////////////////////
    ;-------------------------------------------------------
        ;Determine bin centers and half-widths.
        ;   - Make 3D if 2D.
        if dims[0] eq 2 then begin
            centers = [bin_center[0,*], fltarr(1, dims[1]), bin_center[1,*]]
            widths  = [half_width[0,*], fltarr(1, dims[1]), half_width[1,*]]
        endif else begin
            centers = [bin_center[0,*], bin_center[1,*], bin_center[2,*]]
            widths  = [half_width[0,*], half_width[1,*], half_width[2,*]]
        endelse
        
        ;Get positions: [x0, y0, z0, x1, y1, z0]
        positions = [centers[0,*] - widths[0,*], $
                     centers[1,*] - widths[1,*], $
                     centers[2,*] - widths[2,*], $
                     centers[0,*] + widths[0,*], $
                     centers[1,*] + widths[1,*], $
                     centers[2,*] + widths[2,*]]
    endelse

;-------------------------------------------------------
; Overview Plot ////////////////////////////////////////
;-------------------------------------------------------
    if im_name ne '' then begin
        ;Create the color plot
        im_win  = MrSim_ColorSlab(oSim, im_name, C_NAME=c_name)
        im_win -> Refresh, /DISABLE
        
        ;Draw boxes where the distribution functions are being taken.
        theIm = im_win['Color ' + im_name]
        for i = 0, nDist - 1 do begin
            theBin = positions[*,i]
            xpoly  = theBin[[0,3,3,0,0]]
            ypoly  = theBin[[2,2,5,5,2]]
            !Null  = MrPlotS(xpoly, ypoly, TARGET=theIm, NAME='Bin ' + strtrim(i, 2), $
                             /DATA, COLOR='White', THICK=2.0)
        endfor
    endif

;-------------------------------------------------------
; Setup Window & Positions /////////////////////////////
;-------------------------------------------------------
    oxmargin = [10,12]
    xgap     = 0
    ygap     = 0
    charsize = 1.5
    win2 = MrWindow(LAYOUT=layout, CHARSIZE=charsize, OXMARGIN=oxmargin, NAME='Dist Win', $
                    XGAP=xgap, YGAP=ygap, XSIZE=xsize, YSIZE=ysize, REFRESH=0)

;-------------------------------------------------------
; Distribution Functions ///////////////////////////////
;-------------------------------------------------------
    ;Velocity space?
    space = strsplit(type, '-', /EXTRACT)
    xaxis = space[0]
    yaxis = space[1]
    tf_xvel = strmid(xaxis, 0, 1) eq 'V'
    tf_yvel = strmid(xaxis, 0, 1) eq 'V'

    ;Keep track of ranges to make all distributions uniform
    xrange = [!values.f_infinity, -!values.f_infinity]
    yrange = [!values.f_infinity, -!values.f_infinity]
    cmax = 0D
    
    ;Step through each distribution.
    for i = 0, nDist - 1 do begin
        ;Create the distribution function
        ;   - Remove the colorbar
        imgTemp  = MrSim_eDist(oSim, type, centers[*,i], widths[*,i], $
                               NBINS=nBins, CIRCLES=circles, /CURRENT, V_VA=v_va)
        win2    -> Remove, win2['CB: eDist'], /DESTROY

        ;Number the distribution
        distNo               = string(i, FORMAT='(i02)')
        imgTemp.NAME         = 'eDist ' + type + ' ' + distNo
        win2['BinX'].NAME    = 'BinX '  + type + ' ' + distNo
        win2['BinZ'].NAME    = 'BinZ '  + type + ' ' + distNo
;        win2['Circles'].NAME = 'Circles ' + distNo

        ;Get maximum data range
        ;   - [xy]ranges are only applicable to velocity-space
        ;   - [xy]-space is handled below
        imgTemp -> GetProperty, XRANGE=xr, YRANGE=yr, RANGE=r
        cmax     = cmax  > max(r)
        xRange[0] <= xr[0]
        xRange[1] >= xr[1]
        yRange[0] <= yr[0]
        yRange[1] >= yr[1]
    endfor

;-------------------------------------------------------
; Make Pretty //////////////////////////////////////////
;-------------------------------------------------------
    if tf_xvel eq 0 then xrange = x + [-dx, dx]
    if tf_yvel eq 0 then yrange = y + [-dy, dy]
    if v_va    eq 0 then begin
        xtickinterval = 0.5
        ytickinterval = 0.5
    endif

    ;Determine the aspect ratio of the plots and create positions
    aspect = 1.0;(yrange[1] - yrange[0]) / (xrange[1] - xrange[0])
    pos  = MrLayout(layout, CHARSIZE=charsize, OXMARGIN=oxmargin, ASPECT=aspect, $
                    XGAP=xgap, YGAP=ygap)

    ;Step through each distribution
    allIm = win2 -> Get(/ALL, ISA='MrImage')
    foreach imgTemp, allIm do begin
        ;Determine location in layout.
        name = imgTemp.NAME
        i = stregex(name, 'eDist .+ ([0-9]+)', /SUBEXP, /EXTRACT)
        i = fix(i[1])
        
        ;Change the position, axis ranges, and tickmarks
        if tf_xvel then imgTemp.XRANGE = xrange
        if tf_yvel then imgTemp.YRANGE = yrange
        imgTemp -> SetProperty, POSITION=pos[*,i], RANGE=[0, cmax], $
                                XTICKINTERVAL=xtickinterval, XMINOR=5, YTICKINTERVAL=ytickinterval, YMINOR=5

        ;All Except Left Column
        ;   - Do not label y-axis
        if i mod layout[0] ne 0 $
            then imgTemp -> SetProperty, YTITLE='', YTICKFORMAT='(a1)'

        ;All Except Bottom Row
        ;   - Do not label x-axis
        if (layout[1] ne 1) && (i/layout[0] ne layout[1]-1) $
            then imgTemp -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
    endforeach

    if layout[0] gt 3 then begin
        charsize = 1.0
        allText  = win2 -> Get(/ALL, ISA='MrText')
        foreach txt, allText do txt.CHARSIZE=charsize
    endif

;-------------------------------------------------------
; Annotate /////////////////////////////////////////////
;-------------------------------------------------------
    sim_class = obj_class(oSim)
    oSim -> GetProperty, TIME=time, YRANGE=yrange
    oSim -> GetInfo, DTXWCI=dtxwci
    
    ;Create a title
    title  = 'eMap ' + type + ' t*$\Omega$$\downci$='
    title += string(dtxwci*time, FORMAT='(f0.1)')
    if sim_class eq 'MRSIM3D' then title += ' y=' + string(yrange[0], FORMAT='(f0.1)') + 'de'
    
    ;Add the title
    !Null = MrText(0.5, 0.965, title, ALIGNMENT=0.5, /CURRENT, $
                   NAME='eMap Title', /NORMAL, CHARSIZE=2)
    
    ;Create a colorbar
    !Null = MrColorbar(CTINDEX=13, RANGE=[0, cmax], POSITION=[0.9, 0.25, 0.92, 0.75], $
                       TITLE='Counts', /CURRENT, NAME='CB: Counts', $
                       /VERTICAL, TLOCATION='Right')

    win2   -> Refresh
    if obj_valid(im_win) then im_win -> Refresh
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    return, win2
end
