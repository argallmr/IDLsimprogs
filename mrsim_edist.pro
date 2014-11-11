; docformat = 'rst'
;
; NAME:
;    MrSim_eDist
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
;   The purpose of this program is to create an electron distribution function from
;   simulation particle data.
;
; :Categories:
;    Bill Daughton, Simulation
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
;       2014/09/27  -   Written by Matthew Argall
;       2014/09/30  -   Removed the SIM3D keyword and added the THESIM parameter.
;                           Repurposed the SIM_OBJECT keyword. Added the V_VA and CIRCLES
;                           keywords. - MRA
;       2014/10/14  -   Electrons are now read from the file. With fmaps, the method is
;                           faster than the previous method. - MRA
;       2014/10/30  -   Changed input parameters X and Z to BIN_CENTER and HALF_WIDTH.
;                           Added YRANGE parameter to helper functions. - MRA
;-
;*****************************************************************************************
;+
;   Determine the direction parallel to the magnetic field. The average magnetic
;   field throughout the area defined by XRANGE and ZRANGE is used.
;
; :Params:
;       XRANGE:         in, required, type=fltarr(2)
;                       x-range over which particle counts will be collected.
;       YRANGE:         in, required, type=fltarr(2)
;                       y-range over which particle counts will be collected.
;       ZRANGE:         in, required, type=fltarr(2)
;                       z-range over which particle counts will be collected.
;       SIM_OBJECT:     in, required, type=object
;                       A MrSim subclassed object used to read simulation data.
;
; :Returns:
;       B_HAT:          A unit vector pointing in the direction of the magnetic field.
;-
function MrSim_eDist_B_hat, xrange, yrange, zrange, sim_object
    compile_opt strictarr
    on_error, 2
    
    ;Average magnetic field within the bin
    Bx    = sim_object -> ReadGDA('Bx', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    By    = sim_object -> ReadGDA('By', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Bz    = sim_object -> ReadGDA('Bz', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
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
;       YRANGE:         in, required, type=fltarr(2)
;                       y-range over which particle counts will be collected.
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
function MrSim_eDist_Perp1_hat, xrange, yrange, zrange, sim_object, $
B_HAT=b_hat, $
E_HAT=e_hat
    compile_opt strictarr
    on_error, 2
    
    ;Average magnetic field within the bin
    Bx    = sim_object -> ReadGDA('Bx', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    By    = sim_object -> ReadGDA('By', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Bz    = sim_object -> ReadGDA('Bz', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    B_avg = [mean(Bx), mean(By), mean(Bz)]
    B_hat = B_avg / sqrt(total(B_avg^2))
    undefine, Bx, By, Bz
    
    ;Average electric field within the bin
    Ex    = sim_object -> ReadGDA('Ex', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Ey    = sim_object -> ReadGDA('Ey', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Ez    = sim_object -> ReadGDA('Ez', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
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
;       YRANGE:         in, required, type=fltarr(2)
;                       y-range over which particle counts will be collected.
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
function MrSim_eDist_Perp2_hat, xrange, yrange, zrange, sim_object, $
B_HAT=b_hat, $
E_HAT=E_hat, $
P1_HAT=p1_hat
    compile_opt strictarr
    on_error, 2
    
    ;Average magnetic field within the bin
    Bx    = sim_object -> ReadGDA('Bx', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    By    = sim_object -> ReadGDA('By', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Bz    = sim_object -> ReadGDA('Bz', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    B_avg = [mean(Bx), mean(By), mean(Bz)]
    B_hat = B_avg / sqrt(total(B_avg^2))
    undefine, Bx, By, Bz
    
    ;Average electric field within the bin
    Ex    = sim_object -> ReadGDA('Ex', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Ey    = sim_object -> ReadGDA('Ey', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    Ez    = sim_object -> ReadGDA('Ez', XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
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
function MrSim_eDist_x_z, e_data, xrange, yrange, zrange, $
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
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    if dims[0] eq 5 $
        then e_dist = hist_nd(e_data[[0,1], *], NBINS=nBins) $
        else e_dist = hist_nd(e_data[[0,2], *], NBINS=nBins)
    
    ;Create the axis locations
    dims = size(e_dist, /DIMENSIONS)
    xloc = linspace(xrange[0], xrange[1], dims[0])
    zloc = linspace(zRange[0], zRange[1], dims[1])
    
    ;Plotting parameters
    xtitle = 'x (de)'
    ytitle = 'z (de)'

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_x_Vx, e_data, xrange, yrange, zrange, $
NBINS=nBins, $
V_VA=v_va
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)
    
    ;Histogram in velocity space
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    if dims[0] eq 5 then begin
        e_dist = hist_nd(e_data[[0,3], *], NBINS=nBins)
        vxRange = [min(e_data[3,*], MAX=vxMax), vxMax]
    endif else begin
        e_dist = hist_nd(e_data[[0,3], *], NBINS=nBins)
        vxRange = [min(e_data[3,*], MAX=vxMax), vxMax]
    endelse
    
    ;Create the axis locations
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(xrange[0], xrange[1], dims[0])
    yloc    = linspace(vxRange[0], vxRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'x (de)'
    ytitle  = 'V$\downx$' + (v_va ? '/V$\downA$' : '/c')
    vmax    = max(abs(vxRange))
    yrange  = [-vmax, vmax]
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_z_Vz, e_data, xrange, yrange, zrange, $
NBINS=nBins, $
V_VA=v_va
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)
    
    ;Histogram in velocity space
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    if dims[0] eq 5 then begin
        e_dist = hist_nd(e_data[[1,4], *], NBINS=nBins)
        vzRange = [min(e_data[4,*], MAX=vzMax), vzMax]
    endif else begin
        e_dist = hist_nd(e_data[[1,5], *], NBINS=nBins)
        vzRange = [min(e_data[5,*], MAX=vzMax), vzMax]
    endelse
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(zrange[0], zrange[1], dims[0])
    yloc    = linspace(vzRange[0], vzRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'z (de)'
    ytitle  = 'V$\downz$' + (v_va ? '/V$\downA$' : '/c')
    vmax    = max(abs(vzRange))
    yrange  = [-vmax, vmax]
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vx_Vz, e_data, xrange, yrange, zrange, $
NBINS=nBins, $
V_VA=v_va
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)

    ;Histogram in velocity space
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    if dims[0] eq 5 then begin
        e_dist  = hist_nd(e_data[[2,4],*], NBINS=nBins)
        vxRange = [min(e_data[2,*], MAX=xmax), xmax]
        vzRange = [min(e_data[4,*], MAX=zmax), zmax]
    endif else begin
        e_dist  = hist_nd(e_data[[3,5],*], NBINS=nBins)
        vxRange = [min(e_data[3,*], MAX=xmax), xmax]
        vzRange = [min(e_data[5,*], MAX=zmax), zmax]
    endelse
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vxRange[0], vxRange[1], dims[0])
    zloc    = linspace(vzRange[0], vzRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downx$' + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downz$' + (v_va ? '/V$\downA$' : '/c')
    vmax    = max(abs([vxRange, vzRange]))
    xrange  = [-vmax, vmax]
    yrange  = xrange
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vx_Vy, e_data, xrange, yrange, zrange, $
NBINS=nBins, $
V_VA=v_va
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)

    ;Histogram in velocity space
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    if dims[0] eq 5 then begin
        e_dist  = hist_nd(e_data[[2,3],*], NBINS=nBins)
        vxRange = [min(e_data[2,*], MAX=xmax), xmax]
        vyRange = [min(e_data[3,*], MAX=ymax), ymax]
    endif else begin
        e_dist  = hist_nd(e_data[[3,4],*], NBINS=nBins)
        vxRange = [min(e_data[3,*], MAX=xmax), xmax]
        vyRange = [min(e_data[4,*], MAX=ymax), ymax]
    endelse
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vxRange[0], vxRange[1], dims[0])
    yloc    = linspace(vyRange[0], vyRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downx$' + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downy$' + (v_va ? '/V$\downA$' : '/c')
    vmax    = max(abs([vxRange, vyRange]))
    xrange  = [-vmax, vmax]
    yrange  = xrange
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vy_Vz, e_data, xrange, yrange, zrange, $
NBINS=nBins, $
V_VA=v_va
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(img) then img -> Close
        void = cgErrorMsg()
        return, obj_new()
    endif
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)

    ;Histogram in velocity space
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    if dims[0] eq 5 then begin
        e_dist  = hist_nd(e_data[[3,4],*], NBINS=nBins)
        vyRange = [min(e_data[3,*], MAX=vyMax), vyMax]
        vzRange = [min(e_data[4,*], MAX=vzMax), vzMax]
    endif else begin
        e_dist  = hist_nd(e_data[[4,5],*], NBINS=nBins)
        vyRange = [min(e_data[4,*], MAX=vyMax), vyMax]
        vzRange = [min(e_data[5,*], MAX=vzMax), vzMax]
    endelse
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vyRange[0], vyRange[1], dims[0])
    yloc    = linspace(vzRange[0], vzRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downy$' + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downz$' + (v_va ? '/V$\downA$' : '/c')
    vmax    = max(abs([vyRange, vzRange]))
    xrange  = [-vmax, vmax]
    yrange  = xrange
    
    ;Create the distribution function
    img    = MrImage(e_dist, xloc, yloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                     XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vpar_Vperp, e_data, xrange, yrange, zrange, sim_object, $
NBINS=nBins, $
V_VA=v_va
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
    
    ;Inputs
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)
    
    ;Get the average magnetic field in the bin
    B_hat = MrSim_eDist_B_hat(xrange, yrange, zrange, sim_object)
    
    ;Magnitude, parallel, and perpendicular components
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    e_data = transpose(e_data)
    if dims[0] eq 5 then begin
        v_mag = sqrt(total(e_data[*,[2,3,4]]^2, 2))
        v_par  = e_data[*,2] * B_hat[0] + $
                 e_data[*,3] * B_hat[1] + $
                 e_data[*,4] * B_hat[2]
    endif else begin
        v_mag = sqrt(total(e_data[*,[3,4,5]]^2, 2))
        v_par  = e_data[*,3] * B_hat[0] + $
                 e_data[*,4] * B_hat[1] + $
                 e_data[*,5] * B_hat[2]
    endelse
    v_perp = sqrt(v_mag^2 - v_par^2)
    
    ;Histogram the data
    v_par_perp = transpose([[temporary(v_par)], [temporary(v_perp)]])
    e_dist     = hist_nd(v_par_perp, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax       = max(abs(temporary(v_par_perp)))
    vParRange  = [-vMax, vMax]
    vPerpRange = [0, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vParRange[0],  vParRange[1],  dims[0])
    zloc    = linspace(vPerpRange[0], vPerpRange[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\downpar$'  + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downperp$' + (v_va ? '/V$\downA$' : '/c')

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vpar_Vperp1, e_data, xrange, yrange, zrange, sim_object, $
NBINS=nBins, $
V_VA=v_va
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
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)
    
    ;Get the parallel and perp1-directions
    p1_hat = MrSim_eDist_Perp1_hat(xrange, yrange, zrange, sim_object, B_HAT=B_hat)
    
    ;Perpendicular directions
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    e_data  = transpose(e_data)
    if dims[0] eq 5 then begin
        v_par   = e_data[*,2] * B_hat[0] + $
                  e_data[*,3] * B_hat[1] + $
                  e_data[*,4] * B_hat[2]
        v_perp1 = e_data[*,2] * p1_hat[0] + $
                  e_data[*,3] * p1_hat[1] + $
                  e_data[*,4] * p1_hat[2]
    endif else begin
        v_par   = e_data[*,3] * B_hat[0] + $
                  e_data[*,4] * B_hat[1] + $
                  e_data[*,5] * B_hat[2]
        v_perp1 = e_data[*,3] * p1_hat[0] + $
                  e_data[*,4] * p1_hat[1] + $
                  e_data[*,5] * p1_hat[2]
    endelse
    
    ;Histogram the data
    v_par_perp1 = transpose([[temporary(v_par)], [temporary(v_perp1)]])
    e_dist      = hist_nd(v_par_perp1, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax        = max(abs(temporary(v_par_perp1)))
    vParRange   = [-vMax, vMax]
    vPerp1Range = [-vMax, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vParRange[0],   vParRange[1],   dims[0])
    zloc    = linspace(vPerp1Range[0], vPerp1Range[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\down||$'    + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downPerp1$' + (v_va ? '/V$\downA$' : '/c')

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vpar_Vperp2, e_data, xrange, yrange, zrange, sim_object, $
NBINS=nBins, $
V_VA=v_va
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
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)
    
    ;Get the parallel and perp1-directions
    p2_hat = MrSim_eDist_Perp2_hat(xrange, yrange, zrange, sim_object, B_HAT=B_hat)
    
    ;Perpendicular directions
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    e_data  = transpose(e_data)
    if dims[0] eq 5 then begin
        v_par   = e_data[*,2] * B_hat[0] + $
                  e_data[*,3] * B_hat[1] + $
                  e_data[*,4] * B_hat[2]
        v_perp2 = e_data[*,2] * p2_hat[0] + $
                  e_data[*,3] * p2_hat[1] + $
                  e_data[*,4] * p2_hat[2]
    endif else begin
        v_par   = e_data[*,3] * B_hat[0] + $
                  e_data[*,4] * B_hat[1] + $
                  e_data[*,5] * B_hat[2]
        v_perp2 = e_data[*,3] * p2_hat[0] + $
                  e_data[*,4] * p2_hat[1] + $
                  e_data[*,5] * p2_hat[2]
    endelse
    
    ;Histogram the data
    v_par_perp2 = transpose([[temporary(v_par)], [temporary(v_perp2)]])
    e_dist      = hist_nd(v_par_perp2, NBINS=nBins)
    
    ;Determine the data range and dimensions
    vMax        = max(abs(temporary(v_par_perp2)))
    vParRange   = [-vMax, vMax]
    vPerp2Range = [-vMax, vMax]
    
    ;Create the spacial grid
    dims    = size(e_dist, /DIMENSIONS)
    xloc    = linspace(vParRange[0],   vParRange[1],   dims[0])
    zloc    = linspace(vPerp2Range[0], vPerp2Range[1], dims[1])
    
    ;Plotting parameters
    xtitle  = 'V$\down||$'    + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downPerp2$' + (v_va ? '/V$\downA$' : '/c')

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
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
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, then `E_DATA` is in units of v/v_A instead of v/c (i.e.
;                           normalized to the alvfen speed instead of the speed of light).
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist_Vperp1_Vperp2, e_data, xrange, yrange, zrange, sim_object, $
NBINS=nBins, $
V_VA=v_va
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
    
    ;Units
    dims = size(e_data, /DIMENSIONS)
    v_va = keyword_set(v_va)
    
    ;Get the perp1 and perp2-directions
    p2_hat = MrSim_eDist_Perp2_hat(xrange, yrange, zrange, sim_object, P1_HAT=p1_hat)
    
    ;Perpendicular directions
    ;   - [x, z, ux, uy, uz]
    ;   - [x, y, z, ux, uy, uz]
    e_data = transpose(e_data)
    if dims[0] eq 5 then begin
        v_p1  = e_data[*,2] * p1_hat[0] + $
                e_data[*,3] * p1_hat[1] + $
                e_data[*,4] * p1_hat[2]
        v_p2  = e_data[*,2] * p2_hat[0] + $
                e_data[*,3] * p2_hat[1] + $
                e_data[*,4] * p2_hat[2]
    endif else begin
        v_p1  = e_data[*,3] * p1_hat[0] + $
                e_data[*,4] * p1_hat[1] + $
                e_data[*,5] * p1_hat[2]
        v_p2  = e_data[*,3] * p2_hat[0] + $
                e_data[*,4] * p2_hat[1] + $
                e_data[*,5] * p2_hat[2]
    endelse
    
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
    xtitle  = 'V$\downperp1$' + (v_va ? '/V$\downA$' : '/c')
    ytitle  = 'V$\downperp2$' + (v_va ? '/V$\downA$' : '/c')

    ;Create the distribution function
    img = MrImage(e_dist, xloc, zloc, /CURRENT, /SCALE, /AXES, CTINDEX=13, $
                  XTITLE=xtitle, YTITLE=ytitle, TITLE=title, MISSING_VALUE=0, MISSING_COLOR='White', $
                  NAME='eDist: Vperp1-Vperp2', POSITION=pos, XRANGE=vPerp1Range, YRANGE=vPerp2Range)
    
    return, img
end


;+
;   Create an electron distribution function.
;
; :Params:
;       THESIM:         in, required, type=string/integer/object
;                       The name or number of the simulation to be used, or a
;                           valid simulation object. See MrSim_Create.pro.
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
;       BIN_CENTER:     in, required, type=fltarr(2)/fltarr(3)
;                       An array specifying the coordinates of the distribution function
;                           center. For 2D simulations, this is a 2-element array, and for
;                           3D, a 3-element array is required.
;       HALF_WIDTH:     in, required, type=fltarr(2)/fltarr(3)
;                       Half-widge of each dimension of the bin.
;
; :Keywords:
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, velocity will be normalized to v_A instead of c::
;                               c / vA = sqrt(mi_me) * f_pe / f_ce.
;                           so multiplying v/c by c from the above equation does the trick.
;       CIRCLES:        in, optional, type=boolean, default=0
;                       If set, concentric cirlces will be drawn at v=0.25 and v=0.5.
;       CURRENT:        in, optional, type=boolean, default=0
;                       If set, the distribution function will be placed in the current
;                           MrWindow graphics window.
;       ENERGY:         in, optional, type=boolean, default=0
;                       If set, the distribution function will be plotted in energy space.
;                           The default is velocity space.
;       FILENAME:       in, optional, type=string, default=simulation specific
;                       Name of the binary file containing electron distribution functions
;                           data. The default comes from MrSim_Which.pro
;       FMAP_DIR:       in, optional, type=string
;                       Directory in which to find fMap save files. The default is
;                           determined by MrSim_Which.
;       MOMENTUM:       in, optional, type=boolean, default=0
;                       If set, the distribution function will be plotted in momentum
;                           space. The default is velocity space.
;       NBINS:          in, optional, type=integer/intarr(2), default=256
;                       Number of bins to use when grouping counts. E.g., for Vx-Vz
;                           distrubution, the ranges of Vx and Vz will be broken into
;                           this many bins.
;       SIM_OBJECT:     out, optional, type=object
;                       If `THESIM` is the name or number of a simulation, then
;                           this keyword returns the object reference to the
;                           corresponding simulation object that is created.
;       V_VA:           in, optional, type=boolean, default=0
;                       If set, velocity will be normalized to v_A instead of c::
;                               c / vA = sqrt(mi_me) * f_pe / f_ce.
;                           so multiplying v/c by the above equation does the trick.
;       _REF_EXTRA:     in, optional, type=structure
;                       Any keyword accepted MrSim_Create.pro. Ignored if `THESIM`
;                           is an object.
;
; :Returns:
;       IMG:            A MrImage object refernce to the plot of the distribution function.
;-
function MrSim_eDist, theSim, type, bin_center, half_width, $
CIRCLES=circles, $
CURRENT=current, $
ENERGY=energy, $
FMAP_DIR=fmap_dir, $
FILENAME=filename, $
MOMENTUM=momentum, $
NBINS=nBins, $
SIM_OBJECT=oSim, $
V_VA=v_va, $
_REF_EXTRA=ref_extra
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if obj_valid(win) && keyword_set(current) eq 0 then obj_destroy, win
        if osim_created   && arg_present(oSim)    eq 0 then obj_destroy, oSim
        void = cgErrorMsg()
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
    _type    = strupcase(type)
    circles  = keyword_set(circles)
    current  = keyword_set(current)
    energy   = keyword_set(energy)
    momentum = keyword_set(momentum)
    v_va     = keyword_set(v_va)
    if n_elements(type)     eq 0 then type     = 'Vx-Vz'
    if n_elements(dx)       eq 0 then dx       = 1.5
    if n_elements(dz)       eq 0 then dz       = 1.5
    if n_elements(nBins)    eq 0 then nBins    = 75
    if n_elements(filename) eq 0 then filename = ''
    
    ;Dependencies
    if energy + momentum gt 1 then $
        message, 'Keywords ENERGY and MOMENTUM are mutually exclusive.'
    
    ;Location and size of the distribution function
    center = n_elements(bin_center) eq 2 ? [bin_center[0], 0, bin_center[1]] : bin_center
    width  = n_elements(half_width) eq 0 ? [half_width[0], 0, half_width[1]] : half_width
    
    ;Volume of distribution function.
    xrange = center[0] + [-width[0], width[0]]
    yrange = center[1] + [-width[1], width[1]]
    zrange = center[2] + [-width[2], width[2]]
    
;-------------------------------------------------------
;Select Data ///////////////////////////////////////////
;-------------------------------------------------------
    ;Pick out the data
    e_data = oSim -> ReadElectrons(filename, FMAP_DIR=fmap_dir, /VELOCITY, $
                                             XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange)
    if n_elements(e_data) eq 0 then return, !Null
    
    ;Convert to units of v/vA from v/c?
    if v_va then begin
        oSim -> GetInfo, MI_ME=mi_me, WPE_WCE=wpe_wce

        if n_elements(mi_me) eq 0 then begin
            message, 'mi_me unknown. Make sure ASCII info file exists. Setting V_VA=0', /INFORMATIONAL
            v_va = 0

            ;Special case for Sim1
;            mi_me   = 200.0
;            wpe_wce = 2.0
;            e_data[2:4,*] *= sqrt(mi_me) * wpe_wce
        endif else begin
            e_data[2:4,*] *= sqrt(mi_me) * wpe_wce
        endelse
    endif
    
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
        'X-Z':           img = MrSim_eDist_x_z(e_data,   xrange, yrange, zrange, NBINS=nBins)
        'X-VX':          img = MrSim_eDist_x_Vx(e_data,  xrange, yrange, zrange, NBINS=nBins, V_VA=v_va)
        'Z-VZ':          img = MrSim_eDist_z_Vz(e_data,  xrange, yrange, zrange, NBINS=nBins, V_VA=v_va)
        'VX-VZ':         img = MrSim_eDist_Vx_Vz(e_data, xrange, yrange, zrange, NBINS=nBins, V_VA=v_va)
        'VX-VY':         img = MrSim_eDist_Vx_Vy(e_data, xrange, yrange, zrange, NBINS=nBins, V_VA=v_va)
        'VY-VZ':         img = MrSim_eDist_Vy_Vz(e_data, xrange, yrange, zrange, NBINS=nBins, V_VA=v_va)
        'VPAR-VPERP':    img = MrSim_eDist_Vpar_Vperp(e_data,    xrange, yrange, zrange, oSim, NBINS=nBins, V_VA=v_va)
        'VPAR-VPERP1':   img = MrSim_eDist_Vpar_Vperp1(e_data,   xrange, yrange, zrange, oSim, NBINS=nBins, V_VA=v_va)
        'VPAR-VPERP2':   img = MrSim_eDist_Vpar_Vperp2(e_data,   xrange, yrange, zrange, oSim, NBINS=nBins, V_VA=v_va)
        'VPERP1-VPERP2': img = MrSim_eDist_Vperp1_Vperp2(e_data, xrange, yrange, zrange, oSim, NBINS=nBins, V_VA=v_va)
        else:            message, 'Distribution type "' + type + '" not recognized.'
    endcase
    e_data = !Null
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    
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
    if circles && _type ne 'VPAR-VPERP' then begin
        !Null = MrCircle([0.25, 0.5], 0, 0, LINESTYLE=0, TARGET=img, NAME='Circles', $
                         /DATA, COLOR='Magenta', FILL_BACKGROUND=0)
    endif
    
    ;Add text annotations
    xText = string(FORMAT='(%"X=%0.1f$\\+-$%0.1f")', center[0], width[0])
    zText = string(FORMAT='(%"Z=%0.1f$\\+-$%0.1f")', center[2], width[2])
    !Null = MrText(0.04, 0.87, xText, /RELATIVE, TARGET=img, NAME='BinX')
    !Null = MrText(0.04, 0.05, zText, /RELATIVE, TARGET=img, NAME='BinZ')
    
    ;Colorbar
    cb  = MrColorbar(NAME='CB: eDist', TARGET=img, /CURRENT, TITLE='Counts')
    
    ;Return
    win -> Refresh, DISABLE=~refresh_in
    return, img
end