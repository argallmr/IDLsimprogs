;+
;   Determine the direction parallel to the magnetic field. The average magnetic
;   field throughout the area defined by XRANGE and ZRANGE is used.
;
; :Params:
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;       SIM_OBJECT:     in, required, type=object
;                       A MrSim subclassed object used to read simulation data.
;
; :Returns:
;       B_HAT:          A unit vector pointing in the direction of the magnetic field.
;-
function MrSim_eDist_B_hat, xrange, zrange, sim_object
    compile_opt strictarr
    on_error, 2
    
    ;Average magnetic field within the bin
    Bx    = sim_object -> GetMoment('Bx', xrange, zrange)
    By    = sim_object -> GetMoment('By', xrange, zrange)
    Bz    = sim_object -> GetMoment('Bz', xrange, zrange)
    B_avg = [mean(Bx), mean(By), mean(Bz)]
    B_hat = B_avg / sqrt(total(B_avg^2))
    undefine, Bx, By, Bz
        
    ;Return the unit vector
    return, B_hat
end


;+
;   Determine the Perp-1 direction (same as ExB-direction). The average magnetic and
;   electric fields throughout the area defined by XRANGE and ZRANGE is used.
;
; :Params:
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;       SIM_OBJECT:     in, required, type=object
;                       A MrSim subclassed object used to read simulation data.
;
; :Keywords:
;       B_HAT:          out, optional, type=float(3)
;                       A unit vector pointing in the direction of the magnetic field.
;       E_HAT:          out, optional, type=float(3)
;                       A unit vector pointing in the direction of the electric field.
;
; :Returns:
;       P1_HAT:         A unit vector pointing in the ExB-direction, perpendicular to the
;                           magnetic field.
;-
function MrSim_eDist_Perp1_hat, xrange, zrange, sim_object, $
B_HAT=b_hat, $
E_HAT=e_hat
    compile_opt strictarr
    on_error, 2
    
    ;Average magnetic field within the bin
    Bx    = sim_object -> GetMoment('Bx', xrange, zrange)
    By    = sim_object -> GetMoment('By', xrange, zrange)
    Bz    = sim_object -> GetMoment('Bz', xrange, zrange)
    B_avg = [mean(Bx), mean(By), mean(Bz)]
    B_hat = B_avg / sqrt(total(B_avg^2))
    undefine, Bx, By, Bz
    
    ;Average electric field within the bin
    Ex    = sim_object -> GetMoment('Ex', xrange, zrange)
    Ey    = sim_object -> GetMoment('Ey', xrange, zrange)
    Ez    = sim_object -> GetMoment('Ez', xrange, zrange)
    E_avg = [mean(Ex), mean(Ey), mean(Ez)]
    E_hat = E_avg / sqrt(total(E_avg^2))
    undefine, Ex, Ey, Ez
    
    ;perp1 -- ExB direction
    perp1  = [ E_hat[1]*B_hat[2] - E_hat[2]*B_hat[1], $
               E_hat[2]*B_hat[0] - E_hat[0]*B_hat[2], $
               E_hat[0]*B_hat[1] - E_hat[1]*B_hat[0] ]
    p1_hat = perp1 / sqrt(total(perp1^2))
    
    ;Return the unit vector
    return, p1_hat
end


;+
;   Determine the Perp-2 direction (same as the Bx(ExB)-direction). The average magnetic
;   and electric fields throughout the area defined by XRANGE and ZRANGE is used.
;
; :Params:
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;       SIM_OBJECT:     in, required, type=object
;                       A MrSim subclassed object used to read simulation data.
;
; :Keywords:
;       B_HAT:          out, optional, type=float(3)
;                       A unit vector pointing in the direction of the magnetic field.
;       E_HAT:          out, optional, type=float(3)
;                       A unit vector pointing in the direction of the electric field.
;       P1_HAT:         A unit vector pointing in the ExB-direction, perpendicular to the
;                           magnetic field.
;
; :Returns:
;       P2_HAT:         A unit vector pointing in the Bx(ExB)-direction, perpendicular to
;                           the magnetic field and the ExB-drift direction.
;-
function MrSim_eDist_Perp2_hat, xrange, zrange, sim_object, $
B_HAT=b_hat, $
E_HAT=E_hat, $
P1_HAT=p1_hat
    compile_opt strictarr
    on_error, 2
    
    ;Average magnetic field within the bin
    Bx    = sim_object -> GetMoment('Bx', xrange, zrange)
    By    = sim_object -> GetMoment('By', xrange, zrange)
    Bz    = sim_object -> GetMoment('Bz', xrange, zrange)
    B_avg = [mean(Bx), mean(By), mean(Bz)]
    B_hat = B_avg / sqrt(total(B_avg^2))
    undefine, Bx, By, Bz
    
    ;Average electric field within the bin
    Ex    = sim_object -> GetMoment('Ex', xrange, zrange)
    Ey    = sim_object -> GetMoment('Ey', xrange, zrange)
    Ez    = sim_object -> GetMoment('Ez', xrange, zrange)
    E_avg = [mean(Ex), mean(Ey), mean(Ez)]
    E_hat = E_avg / sqrt(total(E_avg^2))
    undefine, Ex, Ey, Ez
    
    ;perp1 -- ExB direction
    perp1  = [ E_hat[1]*B_hat[2] - E_hat[2]*B_hat[1], $
               E_hat[2]*B_hat[0] - E_hat[0]*B_hat[2], $
               E_hat[0]*B_hat[1] - E_hat[1]*B_hat[0] ]
    p1_hat = perp1 / sqrt(total(perp1^2))
    
    ;Perp2 -- Bx(ExB) direction
    perp2  = [ B_hat[1]*p1_hat[2] - B_hat[2]*p1_hat[1], $
               B_hat[2]*p1_hat[0] - B_hat[0]*p1_hat[2], $
               B_hat[0]*p1_hat[1] - B_hat[1]*p1_hat[0] ]
    p2_hat = perp2 / sqrt(total(perp2^2))
    
    ;Return the unit vector
    return, p2_hat
end


;+
;   Create a distribution function in x-Vx space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_x_z, e_data, xrange, zrange, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Histogram in velocity space
    e_dist = hist_nd(e_data[[0,1], *], NBINS=nBins)
    
    ;Determine the data range and dimensions
    dims   = size(e_dist, /DIMENSIONS)
    
    ;Create the axis locations
    xloc    = linspace(xrange[0], xrange[1], dims[0])
    zloc    = linspace(zRange[0], zRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'x (de)'
    ytitle  = 'z (de)'

    ;Create the distribution function
    img    = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                     NAME='eDist: x-z', XRANGE=xrange, YRANGE=yrange)
    
    return, img
end


;+
;   Create a distribution function in x-Vx space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_x_Vx, e_data, xrange, zrange, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Histogram in velocity space
    e_dist = hist_nd(e_data[[0,2], *], NBINS=nBins)
    
    ;Determine the data range and dimensions
    dims    = size(e_dist, /DIMENSIONS)
    vxRange = [min(e_data[4,*], MAX=vxMax), vxMax]
    
    ;Create the axis locations
    xloc    = linspace(xrange[0], xrange[1], dims[0])
    yloc    = linspace(vxRange[0], vxRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'x (de)'
    ytitle  = 'V$\downx$'
    vmax    = max(abs(vxRange))
    yrange  = [-vmax, vmax]
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                     NAME='eDist: x-Vx', XRANGE=xrange, YRANGE=yrange)
    
    return, img
end


;+
;   Create a distribution function in z-Vz space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_z_Vz, e_data, xrange, zrange, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Histogram in velocity space
    e_dist = hist_nd(e_data[[1,4], *], NBINS=nBins)
    
    ;Determine the data range and dimensions
    dims    = size(e_dist, /DIMENSIONS)
    vzRange = [min(e_data[4,*], MAX=vzMax), vzMax]
    
    ;Create the spacial grid
    xloc    = linspace(zrange[0], zrange[1], dims[0])
    yloc    = linspace(vzRange[0], vzRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'z (de)'
    ytitle  = 'V$\downz$'
    vmax    = max(abs(vzRange))
    yrange  = [-vmax, vmax]
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                     NAME='eDist: z-Vz', XRANGE=xrange, YRANGE=yrange)
    
    return, img
end


;+
;   Create a distribution function in Vx-Vz space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vx_Vz, e_data, xrange, zrange, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif

    ;Histogram in velocity space
    e_dist  = hist_nd(e_data[[2,4],*], NBINS=nBins)
    
    ;Determine the data range and dimensions
    dims    = size(e_dist, /DIMENSIONS)
    vxRange = [min(e_data[2,*], MAX=xmax), xmax]
    vzRange = [min(e_data[4,*], MAX=zmax), zmax]
    
    ;Create the spacial grid
    xloc    = linspace(vxRange[0], vxRange[1], dims[0])
    zloc    = linspace(vzRange[0], vzRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downx$'
    ytitle  = 'V$\downz$'
    vmax    = max(abs([vxRange, vzRange]))
    xrange  = [-vmax, vmax]
    yrange  = xrange
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                     NAME='eDist: Vx-Vz', POSITION=pos, XRANGE=xrange, YRANGE=yrange)
    
    return, img
end


;+
;   Create a distribution function in Vx-Vy space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vx_Vy, e_data, xrange, zrange, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif

    ;Histogram in velocity space
    e_dist  = hist_nd(e_data[[2,3],*], NBINS=nBins)
    
    ;Determine the data range and dimensions
    dims    = size(e_dist, /DIMENSIONS)
    vxRange = [min(e_data[2,*], MAX=xmax), xmax]
    vyRange = [min(e_data[3,*], MAX=ymax), ymax]
    
    ;Create the spacial grid
    xloc    = linspace(vxRange[0], vxRange[1], dims[0])
    yloc    = linspace(vyRange[0], vyRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downx$'
    ytitle  = 'V$\downy$'
    vmax    = max(abs([vxRange, vyRange]))
    xrange  = [-vmax, vmax]
    yrange  = xrange
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                     NAME='eDist: Vx-Vy', POSITION=pos, XRANGE=xrange, YRANGE=yrange)
    
    return, img
end


;+
;   Create a distribution function in Vx-Vy space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vy_Vz, e_data, xrange, zrange, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif

    ;Histogram in velocity space
    e_dist  = hist_nd(e_data[[3,4],*], NBINS=nBins)
    
    ;Determine the data range and dimensions
    dims    = size(e_dist, /DIMENSIONS)
    vyRange = [min(e_data[3,*], MAX=vyMax), vyMax]
    vzRange = [min(e_data[4,*], MAX=vzMax), vzMax]
    
    ;Create the spacial grid
    xloc    = linspace(vyRange[0], vyRange[1], dims[0])
    yloc    = linspace(vzRange[0], vzRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downy$'
    ytitle  = 'V$\downz$'
    vmax    = max(abs([vyRange, vzRange]))
    xrange  = [-vmax, vmax]
    yrange  = xrange
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                     NAME='eDist: Vy-Vz', POSITION=pos, XRANGE=xrange, YRANGE=yrange)
    
    return, img
end



;+
;   Create a distribution function in Vx-Vz space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vpar_Vperp, e_data, xrange, zrange, sim_object, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        if (size(e_data, /DIMENSIONS))[0] ne 5 then e_data = transpose(e_data)
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Get the average magnetic field in the bin
    B_hat = MrSim_eDist_B_hat(xrange, zrange, sim_object)
    
    ;Magnitude, parallel, and perpendicular components
    e_data = transpose(e_data)
    v_mag  = sqrt(total(e_data[*,[2,3,4]]^2, 2))
    v_par  = e_data[*,2] * B_hat[0] + $
             e_data[*,3] * B_hat[1] + $
             e_data[*,4] * B_hat[2]
    v_perp = sqrt(v_mag^2 - v_par^2)
    
    ;Histogram the data
    v_par_perp = transpose([[temporary(v_par)], [temporary(v_perp)]])
    e_dist     = hist_nd(v_par_perp, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax       = max(abs(v_par_perp))
    vParRange  = [-vMax, vMax]
    vPerpRange = [0, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vParRange[0],  vParRange[1],  dims[0])
    zloc    = linspace(vPerpRange[0], vPerpRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downpar$'
    ytitle  = 'V$\downperp$'

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                  NAME='eDist: Vpar-Vperp', POSITION=pos, XRANGE=vParRange, YRANGE=vPerpRange)
    
    return, img
end


;+
;   Create a distribution function in Vx-Vz space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vpar_Vperp1, e_data, xrange, zrange, sim_object, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        if (size(e_data, /DIMENSIONS))[0] ne 5 then e_data = transpose(e_data)
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Get the parallel and perp1-directions
    p1_hat = MrSim_eDist_Perp1_hat(xrange, zrange, sim_object, B_HAT=B_hat)
    
    ;Perpendicular directions
    e_data  = transpose(e_data)
    v_par   = e_data[*,2] * B_hat[0] + $
              e_data[*,3] * B_hat[1] + $
              e_data[*,4] * B_hat[2]
    v_perp1 = e_data[*,2] * p1_hat[0] + $
              e_data[*,3] * p1_hat[1] + $
              e_data[*,4] * p1_hat[2]
    
    ;Histogram the data
    v_par_perp = transpose([[temporary(v_par)], [temporary(v_perp1)]])
    e_dist     = hist_nd(v_par_perp, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax        = max(abs(v_par_perp))
    vParRange   = [-vMax, vMax]
    vPerp1Range = [-vMax, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vParRange[0],   vParRange[1],   dims[0])
    zloc    = linspace(vPerp1Range[0], vPerp1Range[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\down||$'
    ytitle  = 'V$\downPerp1$'

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                  NAME='eDist: Vpar-Vperp1', POSITION=pos, XRANGE=vParRange, YRANGE=vPerp1Range)
    
    return, img
end


;+
;   Create a distribution function in Vx-Vz space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vpar_Vperp2, e_data, xrange, zrange, sim_object, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        if (size(e_data, /DIMENSIONS))[0] ne 5 then e_data = transpose(e_data)
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Get the parallel and perp1-directions
    p2_hat = MrSim_eDist_Perp2_hat(xrange, zrange, sim_object, B_HAT=B_hat)
    
    ;Perpendicular directions
    e_data  = transpose(e_data)
    v_par   = e_data[*,2] * B_hat[0] + $
              e_data[*,3] * B_hat[1] + $
              e_data[*,4] * B_hat[2]
    v_perp2 = e_data[*,2] * p2_hat[0] + $
              e_data[*,3] * p2_hat[1] + $
              e_data[*,4] * p2_hat[2]
    
    ;Histogram the data
    v_par_perp = transpose([[temporary(v_par)], [temporary(v_perp2)]])
    e_dist     = hist_nd(v_par_perp, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax        = max(abs(v_par_perp))
    vParRange   = [-vMax, vMax]
    vPerp2Range = [-vMax, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vParRange[0],   vParRange[1],   dims[0])
    zloc    = linspace(vPerp2Range[0], vPerp2Range[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\down||$'
    ytitle  = 'V$\downPerp2$'

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                  NAME='eDist: Vpar-Vperp2', POSITION=pos, XRANGE=vParRange, YRANGE=vPerp2Range)
    
    return, img
end


;+
;   Create a distribution function in Vx-Vz space.
;
; :Params:
;       E_DATA:         in, required, type=5xN float
;                       Electron distribution data.
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;
; :Keywords:
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vperp1_Vperp2, e_data, xrange, zrange, sim_object, $
NBINS=nBins
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        if (size(e_data, /DIMENSIONS))[0] ne 5 then e_data = transpose(e_data)
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Get the perp1 and perp2-directions
    p2_hat = MrSim_eDist_Perp2_hat(xrange, zrange, sim_object, P1_HAT=p1_hat)
    
    ;Perpendicular directions
    e_data = transpose(e_data)
    v_p1  = e_data[*,2] * p1_hat[0] + $
            e_data[*,3] * p1_hat[1] + $
            e_data[*,4] * p1_hat[2]
    v_p2  = e_data[*,2] * p2_hat[0] + $
            e_data[*,3] * p2_hat[1] + $
            e_data[*,4] * p2_hat[2]
    
    ;Histogram the data
    v_perp_perp = transpose([[temporary(v_p1)], [temporary(v_p2)]])
    e_dist     = hist_nd(v_perp_perp, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax        = max(abs(v_perp_perp))
    vPerp1Range = [-vMax, vMax]
    vPerp2Range = [-vMax, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vPerp1Range[0], vPerp1Range[1], dims[0])
    zloc    = linspace(vPerp2Range[0], vPerp2Range[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downperp1$'
    ytitle  = 'V$\downperp2$'

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, $
                  NAME='eDist: Vperp1-Vperp2', POSITION=pos, XRANGE=vPerp1Range, YRANGE=vPerp2Range)
    
    return, img
end


;+
;   Create an electron distribution function.
;
; :Params:
;       TYPE:           in, optional, type=string, default='Vx-Vz'
;                       The type of distribution function to make. Choices are::
;                           'x-Vx'
;                           'z-Vz'
;                           'Vx-Vz'
;                           'Vx-Vy'
;                           'Vy-Vz'
;                           'Vpar-Vperp'
;                           'Vpar-Vperp1'
;                           'Vpar-Vperp2'
;                           'Vperp1-Vperp2'
;       X:              in, required, type=float
;                       X-location of the center of the spacial bin from which electron
;                           distribution information is collected.
;       Z:              in, required, type=float
;                       Z-location of the center of the spacial bin from which electron
;                           distribution information is collected.
;       DX:             in, optional, type=float, default=1
;                       Half-width of the bin.
;       DZ:             in, optional, type=float, default=1
;                       Half-height of the bin.
;
; :Keywords:
;       CURRENT:        in, optional, type=boolean, default=0
;                       If set, the distribution function will be placed in the current
;                           MrWindow graphics window.
;       ENERGY:         in, optional, type=boolean, default=0
;                       If set, the distribution function will be plotted in energy space.
;                           The default is velocity space.
;       FILENAME:       in, optional, type=string, default=simulation specific
;                       Name of the binary file containing electron distribution functions
;                           data.
;       MOMENTUM:       in, optional, type=boolean, default=0
;                       If set, the distribution function will be plotted in momentum
;                           space. The default is velocity space.
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;       SIM3D:          in, optional, type=boolean, default=0
;                       If `SIM_OBJECT` is not given, then set this keyword to indicate
;                           that a 3D simulation is to be used. The default is to use
;                           a 2D simulation.
;       SIM_OBJECT:     in, optional, type=object
;                       A MrSim object or one of its subclasses that contains information
;                           about the simulation to be used. Not specified, one will
;                           be created. See `SIM3D`.
;       _REF_EXTRA:     in, optional, type=any
;                       If `SIM_OBJECT` is not given, then use any keyword available to
;                           MrSim and its subclasses to initialize one.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist, type, x, z, dx, dz, $
CURRENT=current, $
ENERGY=energy, $
FILENAME=filename, $
MOMENTUM=momentum, $
NBINS=nBins, $
SIM3D=sim3d, $
SIM_OBJECT=oSim, $
_REF_EXTRA=ref_extra
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(win)  && keyword_set(current) eq 0 then obj_destroy, win
        if obj_valid(oSim) && arg_present(oSim)    eq 0 then obj_destroy, oSim
        void = cgErrorMsg()
        return, obj_new()
    endif

    ;Defaults
    _type    = strupcase(type)
    current  = keyword_set(current)
    energy   = keyword_set(energy)
    momentum = keyword_set(momentum)
    if n_elements(type)     eq 0 then type     = 'Vx-Vz'
    if n_elements(dx)       eq 0 then dx       = 1.5
    if n_elements(dz)       eq 0 then dz       = 1.5
    if n_elements(nBins)    eq 0 then nBins    = 75
    if n_elements(filename) eq 0 then filename = ''
    
    ;Dependencies
    if energy + momentum gt 1 then $
        message, 'Keywords ENERGY and MOMENTUM are mutually exclusive.'
    
    ;Create the spacial range of the distribution function
    xrange = x + [-dx, dx]
    zrange = z + [-dz, dz]

    ;Create a simulation object
    if n_elements(oSim) eq 0 then begin
        if keyword_set(sim3d) then sim_class = 'MrSim3D' else sim_class = 'MrSim2D'
        oSim = obj_new(sim_class, time, _STRICT_EXTRA=extra)
        if obj_valid(oSim) eq 0 then return, obj_new()
    endif else begin
        if obj_isa(oSim, 'MRSIM') eq 0 then $
            message, 'SIM_OBJECT must be a subclass of "MrSim"'
    endelse
    
;-------------------------------------------------------
;Select Data ///////////////////////////////////////////
;-------------------------------------------------------

    ;Pick out the data
    e_data = oSim -> GetElectrons(xrange, zrange, /VELOCITY, FILENAME=filename)
    if n_elements(e_data) eq 0 then return, !Null
    
    ;Count factor
    oSim -> GetInfo, ECOUNTFACTOR=eCountFactor
    
;-------------------------------------------------------
; Create the Distrubtion ///////////////////////////////
;-------------------------------------------------------
    ;Get a window
    refresh_in = 1
    if current then begin
        win = GetMrWindows(/CURRENT)
        refresh_in = win.REFRESH
        win -> Refresh, /DISABLE
    endif else begin
        win = MrWindow(NAME='eDist', XSIZE=500, YSIZE=400, OXMARGIN=[10,15], REFRESH=0)
    endelse

    ;Make the distribution
    case _type of
        'X-Z':           img = MrSim_eDist_x_z(e_data,   xrange, zrange, NBINS=nBins)
        'X-VX':          img = MrSim_eDist_x_Vx(e_data,  xrange, zrange, NBINS=nBins)
        'Z-VZ':          img = MrSim_eDist_z_Vz(e_data,  xrange, zrange, NBINS=nBins)
        'VX-VZ':         img = MrSim_eDist_Vx_Vz(e_data, xrange, zrange, NBINS=nBins)
        'VX-VY':         img = MrSim_eDist_Vx_Vy(e_data, xrange, zrange, NBINS=nBins)
        'VY-VZ':         img = MrSim_eDist_Vy_Vz(e_data, xrange, zrange, NBINS=nBins)
        'VPAR-VPERP':    img = MrSim_eDist_Vpar_Vperp(e_data,    xrange, zrange, oSim, NBINS=nBins)
        'VPAR-VPERP1':   img = MrSim_eDist_Vpar_Vperp1(e_data,   xrange, zrange, oSim, NBINS=nBins)
        'VPAR-VPERP2':   img = MrSim_eDist_Vpar_Vperp2(e_data,   xrange, zrange, oSim, NBINS=nBins)
        'VPERP1-VPERP2': img = MrSim_eDist_Vperp1_Vperp2(e_data, xrange, zrange, oSim, NBINS=nBins)
        else:         message, 'Distribution type "' + type + '" not recognized.'
    endcase
    e_data = !Null
    
    ;Correct the counts
    img -> GetData, counts
    img -> SetData, temporary(counts)*eCountFactor

;-------------------------------------------------------
;Plot the Distribution /////////////////////////////////
;-------------------------------------------------------
    ;Create a title indicating where the distribution function was taken
    ;   x +/- dx
    ;   z +/- dz
    title = 'Electron Distribution Function'
    
    ;Add circles
    if _type ne 'VPAR-VPERP' then begin
        !Null = MrCircle([0.25, 0.5], 0, 0, LINESTYLE=0, TARGET=img, NAME='Circles', $
                         /DATA, COLOR='Magenta', FILL_BACKGROUND=0)
    endif
    
    ;Add text annotations
    xText = string(FORMAT='(%"X=%0.1f$\\+-$%0.1f")', x, dx)
    zText = string(FORMAT='(%"Z=%0.1f$\\+-$%0.1f")', z, dz)
    !Null = MrText(0.04, 0.87, xText, /RELATIVE, TARGET=img, NAME='BinX')
    !Null = MrText(0.04, 0.05, zText, /RELATIVE, TARGET=img, NAME='BinZ')
    
    ;Colorbar
    cb  = MrColorbar(NAME='CB: eDist', TARGET=img, /CURRENT, TITLE='Counts')
    
    ;Return
    win -> Refresh, DISABLE=~refresh_in
    return, img
end