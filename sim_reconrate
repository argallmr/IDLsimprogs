; docformat = 'rst'
;
; NAME:
;    HALL_ASYM
;
; PURPOSE:
;+
;   The purpose of this program is to create a color image of a single data product. By
;   using the `MyPlot` keyword, many plots may be combined together.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       TIME:                   in, required, type=int
;                               The simulation time for which a velocity vector field
;                                   is to be plotted.
;       NAME:                   in, required, type=string
;                               The name of the vector quantity to be plotted in color.
; :Keywords:
;       DIRECTORY:              in, optional, type=string, default=pwd
;                               The directory in which to find all of the simulation data
;       MYPLOT:                 in, optional, type=object
;                               A "MrPlot" object to which the data will be added and
;                                   displayed.
;       NLEVELS:                in, optional, type=int, default=15
;                               The number of contour lines to draw
;       NSMOOTH:                in, optional, type=int, default=5
;                               The smoothing width to apply to the simulation data.
;       OFILENAME:              out, optional, type=string, default=''
;                               If provided, graphics will be output to a file by this
;                                   name. The file extension is used to determine which
;                                   type of image to create. "PNG" is the default.
;       OSIM:                   out, optional, type=object
;                               The MrReadSim object reference containing the simulation
;                                   data and domain information used in creating the plot.
;       XRANGE:                 in, optional, type=int, default=[1400, 1700]
;                               The x-range (in "de" - electron skin depth) of the
;                                   simulation domain to be plotted.
;       ZRANGE:                 in, optional, type=int, default=[-20, 20]
;                               The z-range (in "de" - electron skin depth) of the
;                                   simulation domain to be plotted.
;       _REF_EXTRA:             in, optional, type=structure
;                               Any keyword accepted by MrPlot__define.pro
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
;        04/30/2013 - Written by Matthew Argall
;
;-
function sim_reconrate, tRange, $
 CURRENT = current, $
 DIAGNOSTIC = diagnostic, $
 DIRECTORY = directory, $
 OFILENAME = ofilename, $
 NSMOOTH = nsmooth, $
 QUIET = quiet, $
 XRANGE = xrange, $
 ZRANGE = zrange, $
_REF_EXTRA = extra
    compile_opt idl2
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(oSim) then obj_destroy, oSim
        if obj_valid(colorWin) then obj_destroy, colorWin
        if n_elements(tempWin) gt 0 then if windowavailable(tempWin) then wdelete, tempWin
        void = error_message()
        return, obj_new
    endif

;-------------------------------------------------------
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------

    ;Set defaults
    current = keyword_set(current)
    diagnostic = keyword_set(diagnostic)
    quiet = keyword_set(quiet)
    if n_elements(ofilename) eq 0 then ofilename = ''
    if n_elements(xrange) eq 0 then xrange = [1470, 1640]
    if n_elements(zrange) eq 0 then zrange = [-20, 20]
    if n_elements(nsmooth) eq 0 then nsmooth = 5

    ;Set the simulation data directory
    if n_elements(directory) eq 0 then void = cgRootName(DIRECTORY=directory)
    
    ;Buffer the output?
    if current eq 0 then $
        if ofilename eq '' then buffer = 0 else buffer = 1
    
    if diagnostic $
        then tempWin = MrGetWindow(/FREE) $
        else tempWin = MrGetWindow(/FREE, /PIXMAP)
        
;-------------------------------------------------------
;Read Data ////////// //////////////////////////////////
;-------------------------------------------------------
    ;Create the object and get the data. Use Ay for contours.
    oSim = obj_new('MrReadSim', DIRECTORY=directory)
    oSim -> GetProperty, XSIM=XSim, ZSIM=ZSim
    
    ;Make sure the values are inside of [XZ]RANGE
    ixrange = getIndexRange(XSim, xrange)
    izrange = getIndexRange(ZSim, zrange)
    x = XSim[ixrange[0]:ixrange[1]]
    z = ZSim[izrange[0]:izrange[1]]

    ;Find the minimum value of ZSim. The separatrices will cross here.
    void = min(abs(z), izmin)

    ;Step through each time
    ER = fltarr(tRange[1]-tRange[0]+1)
    AAy = fltarr(tRange[1]-tRange[0]+1)
    xpt = fltarr(2,tRange[1]-tRange[0]+1)
    position = MrLayout([1,1,1], CHARSIZE=1.5, XMARGIN=[10,15])
    for time = tRange[0], tRange[1] do begin
        ;Get some data
        Ay = oSim -> getData('Ay', time, NSMOOTH=nsmooth, XRANGE=xrange, ZRANGE=zrange)
        Ey = oSim -> getData('Ey', time, NSMOOTH=nsmooth, XRANGE=xrange, ZRANGE=zrange)

        ;Find the minimum value of Ay at Z~0 This will be the separatrix.
        level = min(Ay[*,izmin])

        ;Get the path information of a contour of Ay along the separatrix. Make sure
        ;the contour is open so that the upper and lower separatrix can be directly
        ;compared.
        cgContour, Ay, x, z, POSITION=position, $
                   XRANGE=[x[0],x[-1]], XSTYLE=1, YRANGE=[z[0],z[-1]], YSTYLE=1, $
                   LEVELS=level, CLOSED=0, PATH_INFO=path_info, PATH_XY=path_xy, $
                   /PATH_DATA_COORDS
    
        ;If PATH_INFO does not have two elements, then the separatrix is on a closed
        ;contour. Reduce the XRANGE until the path is open.
        if n_elements(path_info) ne 2 then begin
            if quiet eq 0 $
                then message, 'Separatrix at t=' + strtrim(time,2) + ' is closed. ' + $
                              'Reducing XRANGE by 20%.', /INFORMATIONAL
            
            ;Remove 20% from each side.
            xlen = (xrange[1] - xrange[0]) * 0.1
            xrange += [xlen, -xlen]
            ixrange = getIndexRange(XSim, xrange)
            x = XSim[ixrange[0]:ixrange[1]]
            time -= 1
            continue
        endif
        
        ;Get the index range of the separatrices within PATH_XY
        isep1 = [path_info[0].offset, path_info[0].offset + path_info[0].n - 1]
        isep2 = [path_info[1].offset, path_info[1].offset + path_info[1].n - 1]
        
        ;Interpolate so that the separatrices have values at the same points.
        ;   The number of points in each path with be reduced to the smallest number.
        if path_info[0].n gt path_info[0].n then begin
            px = path_xy[0,isep1[0]:isep1[1]]
            p1z = path_xy[1,isep1[0]:isep1[1]]
            p2z = interpol(path_xy[1,isep2[0]:isep2[1]], path_xy[0,isep2[0]:isep2[1]], px)
        endif else begin
            px = path_xy[0,isep2[0]:isep2[1]]
            p2z = path_xy[1,isep2[0]:isep2[1]]
            p1z = interpol(path_xy[1,isep1[0]:isep1[1]], path_xy[0,isep1[0]:isep1[1]], px)
        endelse

        ;Find the minimum distance between separatrices
        void = min(abs(p1z - p2z), imin)
        
        ;Find where in XSim and ZSim the value of PATH_XY[*,ixpoint] occurs.
        ;   For Z, take the midpoint between the two separatrices.
        ;   Also, Z values are positive and negative.
        void = min(abs(z - p2z[imin]), iz2)
        void = min(abs(z - p1z[imin]), iz1)
        void = min(abs(z - mean([z[iz1], z[iz2]])), iz)
        void = min(abs(x - px[imin]), ix)

        ;Find the average Ey within a 3x3de box around the x-point.
        ER[time-tRange[0]]  = mean(Ey[ix-1:ix+1, iz-2:iz+2])
        AAy[time-tRange[0]] = mean(Ay[ix-1:ix+1, iz-2:iz+2])
        xpt[*,time-tRange[0]] = [x[ix], z[iz]]

        ;Make a plot of the separatrices and the box of integration to make sure
        ;things are being done properly.
        if diagnostic then begin        
            ;Plot Ay
            cgImage, Ey, /AXES, /SCALE, CTINDEX=13, POSITION=position, $
                     XTITLE='X (di)', YTITLE='Z (di)', $
                     TITLE='Ey t=' + string(time/2.0, FORMAT='(f0.1)') + ' $\Omega$$\downi$$\up-1$', $
                     XRANGE=[x[0],x[-1]], YRANGE=[z[0],z[-1]]
            
            ;Draw the contours
            levels = cgConLevels(Ay)
            cgContour, Ay, x, z, C_LABELS=0, LEVELS=levels, /OVERPLOT

            ;Overplot the separatrices onto the image.
            cgPlot, px, p1z, COLOR='White', /OVERPLOT
            cgPlot, px, p2z, COLOR='White', /OVERPLOT
            
            ;Draw an "X" at the X-line and a 1di box around where Ey is being integrated
            cgText, x[ix], z[iz], 'X', COLOR='Black', ALIGNMENT=0.5
            cgPlotS, x[[ix-10,ix+10,ix+10,ix-10,ix-10]], z[[iz-10,iz-10,iz+10,iz+10,iz-10]], $
                     COLOR='black'
        endif
        
        ;If the x-range of the contour is less than 50% of the total x-range,
        ;increase the z-range by 50%.
        if abs(max(px) - min(px)) lt 0.5*(xrange[1] - xrange[0]) then begin
            if quiet eq 0 then $
                message, 'Separatrix at t=' + strtrim(time,2) + ' is too short.' + $
                         'Increasing [X,Z]RANGE by [10,50]%', /INFORMATIONAL

            ;Increase the z-range by 50%
            zheight = 0.25*(zrange[1] - zrange[0])
            zrange += [-zheight, zheight]
            izrange = getIndexRange(ZSim, zrange)
            z = ZSim[izrange[0]:izrange[1]]
            void = min(abs(z), izmin)
            
            ;Increase the x-range by 10% from each side.
            xlen = (xrange[1] - xrange[0]) * 0.05
            xrange += [xlen, -xlen]
            ixrange = getIndexRange(XSim, xrange)
            x = XSim[ixrange[0]:ixrange[1]]
        endif
        
        ;Re-center the viewing box every 10 steps
        if (time mod 10) eq 0 then begin
            if quiet eq 0 then $
                message, string('("Time step %i of %i")', time, tRange[1]-tRange[0])
            
            ;Re-center around the x-point
            xlen = xrange[1] - xrange[0]
            xrange = x[ix] + [-xlen, xlen]
            ixrange = getIndexRange(XSim, xrange)
            x = XSim[ixrange[0]:ixrange[1]]
        endif
    endfor
    
    ;Delete the temporary window.
    wdelete, tempWin
    
    ;Create the time array. Time is 1/2 of the index (not sure why)
    time = linspace(tRange[0], tRange[1], 1, /INTERVAL) / 2.0
    
    ;Calulate the derivative of Ey
    dAy = shift(AAy, -1) - AAy
    dAy = dAy[0:-2]
    
    ;Create a window
    if current $
        then rrateWin = GetMrWindows(/CURRENT) $
        else rrateWin = MrWindow(XSIZE=550, YSIZE=400, LAYOUT=[2,2], BUFFER=buffer, XMARGIN=[10,14], REFRESH=0)
    
    ;Ey
    Ey_plot = MrPlot(time, ER, LAYOUT=[2,2,1,1], /CURRENT, $
                     TITLE='Reconnection Electric Field', YTITLE='<Ey>', $
                     XTITLE='Time $\Omega$$\downion$$\up-1$', NAME='ER')
                     
    ;dAy/dt
    dAy_plot = MrPlot(time, dAy, LAYOUT=[2,2,2,1], /CURRENT, $
                     TITLE='d<Ay>/dt', YTITLE='d<Ay>/dt', NAME='dAy', $
                     XTITLE='Time $\Omega$$\downion$$\up-1$')
    
    ;Location of the X-Point
    xpt_plot = MrPlot(xpt[0,*]/10.0, xpt[1,*]/10.0, LAYOUT=[2,2,1,2], /CURRENT, $
                      TITLE='Location of X-Point', XTITLE='X (di)', YTITLE='Z (di)', $
                      PSYM=-1, SYMCOLOR='Blue', NAME='Xpt')
    
    
    ;Refresh and output, if requested.
    if current eq 0 then begin
        rrateWin -> Refresh
        if ofilename ne '' then rrateWin -> Output, ofilename
    endif
    
    return, rrateWin
end
