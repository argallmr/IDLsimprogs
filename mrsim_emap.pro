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
;       BIN_CENTER:     in, required, type=float(2)/float(3)
;                       Coordinate of the distribution that defines the `LOCATION` of
;                           the eMap. If the particle distribution data file has two
;                           (three) spacial dimensions, then BIN_CENTER has 2 (3) elements.
;                           If an array of bin center locations are given, they define
;                           the center of each distribution in the eMap.
;       HALF_WIDTH:     in, required, type=float(2)/float(3)
;                       Half-width of each dimension of the distribution that defines the
;                           `LOCATION` of the eMap. If the particle distribution data file
;                           has two (three) spacial dimensions, then HALF_WIDTH has 2 (3)
;                           elements. All distributions in the eMap will have the same
;                           half-width.
;
; :Keywords:
;       C_NAME:         in, optional, type=string, default=''
;                       Name of the data product whose contours are to be 
;                           overplotted on the image.
;       CIRCLES:        in, optional, type=boolean, default=0
;                       If set, concentric cirlces will be drawn at v=0.25 and v=0.5.
;       HGAP:           in, optional, type=float, default=0.0
;                       Horizontal space between adjacent distributions. Ignored more than
;                           one `BIN_CENTER` is given. "Horizontal" is taken to mean the
;                           first spacial dimension.
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
;                       Vertical space between adjacent distributions. Ignored if more
;                           than one `BIN_CENTER` is given. "Vertical" is taken to mean
;                           the last spacial dimension.
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
;       2014/11/21  -   Program is more general to simulation dimensions.
;       2015/01/23  -   Set the ASPECT property instead of calculating positions. MrText
;                           uses TeXToIDL instead of cgCheckForSymbols, so update was
;                           required. - MRA.
;       2015/03/09  -   Added the LOG and RANGE keywords. - MRA
;-
function MrSim_eMap, theSim, type, bin_center, half_width, time, $
C_NAME=c_name, $
CIRCLES=circles, $
HGAP=hGap, $
IM_NAME=im_name, $
IM_WIN=im_win, $
NBINS=nBins, $
LAYOUT=layout, $
LOCATION=location, $
LOG=log, $
POSITIONS=positions, $
RANGE=range_in, $
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
        if obj_isa(theSim, 'MRSIM2') eq 0 $
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
    ;Number of distributions to make & number of spacial dimensions (2 or 3) in
    ;the particle distribution file.
    nDims = size(bin_center, /N_DIMENSIONS)
    dims  = size(bin_center, /DIMENSIONS)
    is2D  = dims[0] eq 2
    nDist = nDims eq 1 ? 1 : dims[1]

    log = keyword_set(log)
    if n_elements(nBins)    eq 0 then nBins    = 75
    if n_elements(layout)   eq 0 then layout   = [ceil(sqrt(nDist)), ceil(float(nDist) / ceil(sqrt(nDist)))]
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
    endif else begin
        nDist = product(layout)
    endelse
    
;-------------------------------------------------------
; Position of Distribution Functions ///////////////////
;-------------------------------------------------------
    ;
    ;Position of the distribution indicated by LOCATION
    ;
    
    ;Particle data has 2 spacial dimensions
    ;   - [x0, z0, x1, z1]
    if is2D then begin
        ;Indexing is different for 2D and 3D
        ;   - Number of spacial dimensions
        ;   - Indices of the horizontal coordinates
        ;   - Indices of the vertical coordinates
        nD = 2
        ih = [0,2]
        iv = [1,3]
    
        ;Half widths are the same for all distributions.
        widths = rebin(half_width, nD, nDist)
        
        ;Position of bin(s)
        position = [bin_center[0,*] - widths[0,*], $
                    bin_center[1,*] - widths[1,*], $
                    bin_center[0,*] + widths[0,*], $
                    bin_center[1,*] + widths[1,*]]
        
    ;Particle data has 3 spacial dimensions
    ;   - [x0, y0, z0, x1, y1, z1]
    endif else begin
        ;Indexing is different for 2D and 3D
        ;   - Number of spacial dimensions
        ;   - Indices of the horizontal coordinates
        ;   - Indices of the vertical coordinates
        nD = 3
        ih = [0,3]
        iv = [2,5]
    
        ;Half widths are the same for all distributions.
        widths = rebin(half_width, nD, nDist)
        
        ;Position of bin(s)
        position = [bin_center[0,*] - widths[0,*], $
                    bin_center[1,*] - widths[1,*], $
                    bin_center[2,*] - widths[2,*], $
                    bin_center[0,*] + widths[0,*], $
                    bin_center[1,*] + widths[1,*], $
                    bin_center[2,*] + widths[2,*]]
    endelse

;-------------------------------------------------------
; Determine Location of Dist Fns ///////////////////////
;-------------------------------------------------------
    if nDims eq 1 then begin
        ;Bin-center offsets from the center of the bin defined by LOCATION
        ;   - Dimensions are [position, col, row]
        ;   - Value is +/-N
        ;       o (+/-) indicates left/right (up/down) offset
        ;       o Offset by N bins
        hoffset = indgen(1, layout[0])
        voffset = indgen(1, 1, layout[1])

        ;Horizontal Offset
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
    
        ;Vertical Offset
        ;   - If layout[2] = 2
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
        hoffset = (2*half_width[0] + hGap) * rebin(hoffset, 2, layout[0], layout[1])
        voffset = (2*half_width[1] + vGap) * rebin(voffset, 2, layout[0], layout[1])

        ;Data locations of distribution functions
        ;   - Saved as [position, index]
        ;   - Position is the IDL definition of the word
        ;   - Index starts with 0 in the upper-right corner, then moves right then down
        nDist              = product(layout)
        positions          = rebin(position,   2*nD, layout[0], layout[1])
        positions[ih,*,*] += hoffset
        positions[iv,*,*] += voffset
        positions          = reform(positions, 2*nD, nDist)
        hoffset            = reform(hoffset,   2,    nDist)
        voffset            = reform(voffset,   2,    nDist)

        ;Bin centers
        ;   - Each distribution is offset from the one indicated by LOCATION
        centers = fltarr(nD, nDist)
        centers[0,*] = bin_center[0] + (temporary(hoffset))[0,*]
        if is2D then begin
            centers[1,*] = bin_center[1] + (temporary(voffset))[0,*]
        endif else begin
            centers[1,*] = bin_center[1]
            centers[2,*] = bin_center[2] + (temporary(voffset))[0,*]
        endelse
    endif else begin
        centers   = bin_center
        positions = temporary(position)
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
            xpoly  = theBin[ih[[0,1,1,0,0]]]
            ypoly  = theBin[iv[[0,0,1,1,0]]]
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
    win2 = MrWindow(LAYOUT=layout, CHARSIZE=charsize, OXMARGIN=oxmargin, NAME='eMap Win', $
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
                               /CURRENT, $
                               CIRCLES = circles, $
                               LOG     = log, $
                               NBINS   = nBins, $
                               RANGE   = range, $
                               V_VA    = v_va)

        ;Number the distribution
        distNo            = string(i, FORMAT='(i02)')
        imgTemp.NAME      = 'eDist ' + type + ' ' + distNo
        win2['BinX'].NAME = 'BinX '  + type + ' ' + distNo
        win2['BinZ'].NAME = 'BinZ '  + type + ' ' + distNo
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
    
    ;Remove all colorbars
    win2 -> Remove, TYPE='weColorbar'
    if n_elements(range_in) eq 0 $
        then range = log ? [1, cmax] : [0, cmax] $
        else range = range_in

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
    aspect = (yrange[1] - yrange[0]) / (xrange[1] - xrange[0])
    win2 -> SetProperty, ASPECT=aspect

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
        imgTemp -> SetProperty, RANGE=range, $
                                XTICKINTERVAL=xtickinterval, XMINOR=5, $
                                YTICKINTERVAL=ytickinterval, YMINOR=5

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
    title  = 'eMap ' + type + ' t*\Omega_{ci}='
    title += string(dtxwci*time, FORMAT='(f0.1)')
    if sim_class eq 'MRSIM3D' then title += ' y=' + string(yrange[0], FORMAT='(f0.1)') + 'de'
    
    ;Add the title
    !Null = MrText(0.5, 0.965, title, ALIGNMENT=0.5, /CURRENT, $
                   NAME='eMap Title', /NORMAL, CHARSIZE=2)
    
    ;Create a colorbar
    !Null = MrColorbar(/CURRENT, $
                       CTINDEX   = 13, $
                       YLOG      = log, $
                       NAME      = 'CB: Counts', $
                       POSITION  = [0.9, 0.25, 0.92, 0.75], $
                       RANGE     = range, $
                       TITLE     = 'Counts', $
                       TLOCATION = 'Right', $
                       /VERTICAL)

    win2   -> Refresh
    if obj_valid(im_win) then im_win -> Refresh
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    return, win2
end
