; docformat = 'rst'
;
;+
;   The purpose of this program is to provide a base class for 2D and 3D simulation
;   data provided by Bill Daughton from Los Alamos National Laboratory.
;
;   MrSim is meant to be subclassed by a class that over-rides the ReadData method.
;   At the end of the new ReadData method, ::SetData should be called.
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
; :Categories:
;   Pic Simulation, Data Reader, Bill Daughton
;
;   NOTES:
;       SimNum
;           Simulation-specific parameters will be set according to the SIMNUM property.
;           Use the SIM_INFO keyword to get a description of each simulation.
;
; :Author:
;   Matthew Argall::
;       University of New Hampshire
;       Morse Hall, Room 113
;       8 College Rd.
;       Durham, NH, 03824
;       matthew.argall@wildcats.unh.edu
;
; :Copyright:
;       Copyright 2013 by the University of New Hampshire
;
; :History:
;   Modification History::
;
;       2014/01/30  -   Written by Matthew Argall
;       2014/02/07  -   Added the YRange, YSim, and Orientation properties. Added the
;                           XYCut and XZCut methods, which are called by LineCut.
;                           Transferred other data product methods here from MrSim2D. - MRA
;       2014/02/09  -   A missing data product will issue an informational message and
;                           return !Null instead of issuing an error. This is true for
;                           GetData, LineCuts, XYCuts, and XZCuts. - MRA
;       2014/02/16  -   Added the YSLICE property. - MRA
;       2014/03/28  -   Added divPe_y, UexB_[xyz], V[xyz], VxB_[xyz], JxB_[xyz],
;                           JixB_[xyz], and JexB_[xyz] methods. Renamed *_inPlane methods
;                           to *_xz. - MRA
;       2014/03/29  -   Taking centered differences to obtain derivatives reduces the
;                           effective simulation domain by 2 gridpoints in every direction.
;                           This is now accounted for in the XZCuts method. - MRA
;       2014/04/08  -   Added the derivative keywords to GetData. Added the divPi_x,
;                           divPi_z, D_DX, and D_DZ methods. - MRA
;       2014/05/01  -   Renamed the method E_dot_J to E_par. - MRA
;       2014/07/25  -   Added the Ui_par, Ui_perp, Ue_par, and Ue_perp data products.
;                           Removed the Free_Data method. - MRA
;       2014/07/27  -   Added support for the electron parallel and perpendicular
;                           pressures, anisotropy, agyrotropy, and nongyrotropy. - MRA
;       2014/09/15  -   Added capacity for calculating total pressure. - MRA
;       2014/09/30  -   Added the SetSim method. - MRA.
;       2014/10/03  -   Added the ASCII_VERSION keyword. - MRA
;       2014/10/14  -   Added the GetSlice keyword. Y=0 now occurs at the center grid-cell. - MRA
;-
;*****************************************************************************************
;+
;   This method initializes the MrSim class.
;
; :Params:
;       TIME:           in, optional, type=long, default=0
;                       Time index at which simulation data will be read.
;       YSLICE:         in, optional, type=long, default=0
;                       The Y-index at which to obtain data from the XZ plane. Ignored
;                           in 2D simulations.
;
; :Keywords:
;       ASCII_VERSION:  in, optional, type=integer, default=1
;                       Version of the info file to read. Ignored if `BINARY`=1.
;                           See MrSim_Which.pro.
;       AXIS_LABELS:    in, optional, type=strarr(3), default="['x', 'y', 'z']"
;                       Labels for the axes.
;       BINARY:         in, optional, type=boolean, default=0
;                       If set, use the binary info file instead of the ascii info file.
;       COORD_SYSTEM:   in, optional, type=string, default='SIMULATION'
;                       Coordinate system in which to display the data. Options are::
;                           'SIMULATION'
;                           'MAGNETOPAUSE'
;                           'MAGNETOTAIL'
;       DIRECTORY:      in, optional, type=string, default=pwd
;                       Directory in which to find the ".gda" data.
;       INFO_ASCII:     in, optional, type=string, default=`DIRECTORY`/../info
;                       The ASCII info file containing information about the simulation
;                           setup.
;       INFO_BINARY:    in, optional, type=string, default=`DIRECTORY`/info
;                       The binary info file containing information about the simulation
;                           setup.
;       ION_SCALE:      in, optional, type=boolean, default=0
;                       Construct the simulation domain in units of "di" instead of "de"
;                           (the ion and electron skin depth, respectively).
;       MVA_FRAME:      in, optional, type=boolean, default=0
;                       If `COORD_SYSTEM` is "Magnetopause" or "Magnetotail", then this
;                           this keyword will automatically pick the proper `AXIS_LABELS`.
;       NSMOOTH:        in, optional, type=integer, default=3
;                       Number of points to smooth data after they is read.
;       ORIENTATION:    in, optional, type=string, default='XZ'
;                       The 2D plane in which data will be viewed. Also sets the direction
;                           of "vertical" and "horizontal" when making 1D cuts. Possible
;                           choices are "XY" and "XZ". This keyword has no meaning for
;                           2D simulations.
;       SIM_INFO:       in, optional, type=boolean, default=0
;                       If set, information about each simulation will be printed to the
;                           command window. Initialization will be aborted.
;       XRANGE:         in, optional, type=dblarr(2), default=width of domain
;                       X-range over which data will be read and stored, in units defined
;                           by `ION_SCALE`.
;       YRANGE:         in, optional, type=dblarr(2), default=depth of domain
;                       Y-range over which data will be read and stored, in units defined
;                           by `ION_SCALE`. This keyword is ignored in 2D simulations.
;       ZRANGE:         in, optional, type=dblarr(2), default=height of domain
;                       Z-range over which data will be read and stored, in units defined
;                           by `ION_SCALE`.
;-
function MrSim::INIT, theSim, time, yslice, $
ASCII_VERSION = ascii_version, $
AXIS_LABELS = axis_labels, $
BINARY = binary, $
COORD_SYSTEM = coord_system, $
DIRECTORY = directory, $
INFO_ASCII = info_ascii, $
INFO_BINARY = binary_info, $
ION_SCALE = ion_scale, $
MVA_FRAME = mva_frame, $
NSMOOTH = nsmooth, $
ORIENTATION = orientation, $
SIM_INFO = sim_info, $
XRANGE = xrange, $
YRANGE = yrange, $
ZRANGE = zrange
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, 0
    endif
    

;-------------------------------------------------------
; Simulation Number ////////////////////////////////////
;-------------------------------------------------------
    
    ;List the simulations?
    if keyword_set(sim_info) then begin
        MrSim_Which
        return, 0
    endif      

    ;Get information about the simulation
    MrSim_Which, theSim, NAME=simname, NUMBER=simnum, DIRECTORY=_dir, $
                 INFO_ASCII=_ascii_info, INFO_BINARY=_bin_info, $
                 ASCII_VERSION=vASCII
    
    ;Defaults
    binary = keyword_set(binary)
    if n_elements(directory)     eq 0 then directory     = _dir
    if n_elements(ascii_version) eq 0 then ascii_version = vASCII
    if n_elements(info_ascii)    eq 0 then info_ascii    = _ascii_info
    if n_elements(info_binary)   eq 0 then info_binary   = _bin_info

    ;Properties
    self.simnum    = simnum
    self.simname   = simname
    self.directory = directory
    
    ;Use binary file?
    if ~binary && ~file_test(info_ascii) then begin
        message, 'ASCII info file not found. Switching to binary info file.', /INFORMATIONAL
        binary = 1
    endif
;-------------------------------------------------------
;Read the Info File ////////////////////////////////////
;-------------------------------------------------------
    
    ;Make sure the file exists
    self.info = ptr_new(/ALLOCATE_HEAP)

    ;Binary info file
    if binary then begin
        if keyword_set(ion_scale) then $
            message, 'ION_SCALE is only possible with ASCII info file. ' + $
                     'Setting ION_SCALE=0', /INFORMATIONAL
        ion_scale = 0
        self -> ReadInfo_Binary, info_binary
        
    ;Ascii info file
    endif else begin
        self -> ReadInfo_Ascii, info_ascii, VERSION=ascii_version
    endelse

;-------------------------------------------------------
;Set Defaults //////////////////////////////////////////
;-------------------------------------------------------
    ;Allocate Heap
    self.XSim = ptr_new(/ALLOCATE_HEAP)
    self.YSim = ptr_new(/ALLOCATE_HEAP)
    self.ZSim = ptr_new(/ALLOCATE_HEAP)

    ;Defaults
    ion_scale    = keyword_set(ion_scale)
    mva_frame    = keyword_set(mva_frame)
    coord_system = n_elements(coord_system) eq 0 ? 'SIMULATION'    : strupcase(coord_system)
    _axis_labels = n_elements(axis_labels)  eq 0 ? ['x', 'y', 'z'] : axis_labels
    _orientation = n_elements(orientation)  eq 0 ? 'XZ'            : strupcase(orientation)
    
    ;Convert to units of "di"
    self -> GetInfo, LX_DE=xsize, LY_DE=ysize, LZ_DE=zsize, MI_ME=mi_me
    if ion_scale then begin
        xsize /= sqrt(mi_me)
        ysize /= sqrt(mi_me)
        zsize /= sqrt(mi_me)
    endif
    
    if n_elements(nsmooth) eq 0 then nsmooth = 3
    if n_elements(time)    eq 0 then time    = 0
    if n_elements(xrange)  eq 0 then xrange  = [0, xsize]
    if n_elements(yrange)  eq 0 then yrange  = [-ysize/2.0, ysize/2.0]
    if n_elements(yslice)  eq 0 then yslice  = 0L
    if n_elements(zrange)  eq 0 then zrange  = [-zsize/2.0, zsize/2.0]
    
    if max(coord_system eq ['SIMULATION', 'MAGNETOPAUSE', 'MAGNETOTAIL']) eq 0 then $
        message, 'Coordinate system "' + coord_system + '" not recognized.'
    if max(_orientation eq ['XY', 'XZ', 'YZ']) eq 0 then $
        message, 'Orienatation "' + _orientation + '" not recognized.'
    
    ;MVA Frame?
    if mva_frame then begin
        case coord_system of
            'MAGNETOPAUSE': _axis_labels = ['N', 'M', 'L']      ;N is along X-GSE
            'MAGNETOTAIL':  _axis_labels = ['L', 'M', 'N']      ;N is along Z-GSE
            else:           ;Do nothing
        endcase
    endif

    ;Set Properties
    self.axis_labels  = _axis_labels
    self.coord_system = coord_system
    self.directory    = directory
    self.ion_scale    = ion_scale
    self.mva_frame    = mva_frame
    self.nsmooth      = nsmooth
    self.orientation  = _orientation
    self.time         = time
    self.xrange       = xrange
    self.yrange       = yrange
    self.yslice       = yslice
    self.zrange       = zrange

    ;Set the domain coordinates
    self -> SetScale, ION_SCALE=ion_scale
    self -> MakeSimDomain

;-------------------------------------------------------
;Allocate Heap /////////////////////////////////////////
;-------------------------------------------------------

    ;Data Products
    self.electrons = ptr_new(/ALLOCATE_HEAP)
    self.Ay        = ptr_new(/ALLOCATE_HEAP)
    self.Bx        = ptr_new(/ALLOCATE_HEAP)
    self.By        = ptr_new(/ALLOCATE_HEAP)
    self.Bz        = ptr_new(/ALLOCATE_HEAP)
    self.Ex        = ptr_new(/ALLOCATE_HEAP)
    self.Ey        = ptr_new(/ALLOCATE_HEAP)
    self.Ez        = ptr_new(/ALLOCATE_HEAP)
    self.n_e       = ptr_new(/ALLOCATE_HEAP)
    self.n_i       = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xx     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xy     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xz     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yx     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yy     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yz     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zx     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zy     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zz     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xx     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xy     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xz     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yx     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yy     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yz     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_zx     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_zy     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_zz     = ptr_new(/ALLOCATE_HEAP)
    self.Uex       = ptr_new(/ALLOCATE_HEAP)
    self.Uey       = ptr_new(/ALLOCATE_HEAP)
    self.Uez       = ptr_new(/ALLOCATE_HEAP)
    self.Uix       = ptr_new(/ALLOCATE_HEAP)
    self.Uiy       = ptr_new(/ALLOCATE_HEAP)
    self.Uiz       = ptr_new(/ALLOCATE_HEAP)

    return, 1
end


;+
;   Clean up after the object is destroyed.
;-
pro MrSim::CLEANUP
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Free all of the pointers
    ptr_free, self.info
    ptr_free, self.XSim
    ptr_free, self.ZSim
    
    ;Data products
    ptr_free, self.Ay
    ptr_free, self.Bx
    ptr_free, self.By
    ptr_free, self.Bz
    ptr_free, self.Ex
    ptr_free, self.Ey
    ptr_free, self.Ez
    ptr_free, self.n_e
    ptr_free, self.n_i
    ptr_free, self.Pe_xx
    ptr_free, self.Pe_xy
    ptr_free, self.Pe_xz
    ptr_free, self.Pe_yx
    ptr_free, self.Pe_yy
    ptr_free, self.Pe_yz
    ptr_free, self.Pe_zx
    ptr_free, self.Pe_zy
    ptr_free, self.Pe_zz
    ptr_free, self.Pi_xx
    ptr_free, self.Pi_xy
    ptr_free, self.Pi_xz
    ptr_free, self.Pi_yx
    ptr_free, self.Pi_yy
    ptr_free, self.Pi_yz
    ptr_free, self.Pi_zx
    ptr_free, self.Pi_zy
    ptr_free, self.Pi_zz
    ptr_free, self.Uex
    ptr_free, self.Uey
    ptr_free, self.Uez
    ptr_free, self.Uix
    ptr_free, self.Uiy
    ptr_free, self.Uiz
end


;+
;   The purpose of this method is to release memory taken up by simulation data.
;
; :Private:
;
; :Params:
;       DATA_PRODUCT:           in, optional, type=string
;                               The name of the data product whose data is to be released.
;                                   If not given, the data for all data products is freed.
;-
pro MrSim::Clear_Data, data_product
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Free all of the data
    if n_params() eq 0 then begin
        data_product = ['Ay', 'Bx', 'By', 'Bz', 'E-', 'Ex', 'Ey', 'Ez', 'ne', 'ni', $
                        'Pe-xx', 'Pe-xy', 'Pe-xz', 'Pe-yx', 'Pe-yy', 'Pe-yz', $
                        'Pe-zx', 'Pe-zy', 'Pe-zz', 'Pi-xx', 'Pi-xy', 'Pi-xz', $
                        'Pi-yx', 'Pi-yy', 'Pi-yz', 'Pi-zx', 'Pi-zy', 'Pi-zz', $
                        'Uex', 'Uey', 'Uez', 'Uix', 'Uiy', 'Uiz']
    endif
    
    ;Step through all of the data products given and free the data
    foreach name, data_product do begin
        case name of
            'Ay':    *self.Ay        = !Null
            'Bx':    *self.Bx        = !Null
            'By':    *self.By        = !Null
            'Bz':    *self.Bz        = !Null
            'E-':    *self.electrons = !Null
            'Ex':    *self.Ex        = !Null
            'Ey':    *self.Ey        = !Null
            'Ez':    *self.Ez        = !Null
            'ne':    *self.n_e       = !Null
            'ni':    *self.n_i       = !Null
            'Pe-xx': *self.Pe_xx     = !Null
            'Pe-xy': *self.Pe_xy     = !Null
            'Pe-xz': *self.Pe_xz     = !Null
            'Pe-yx': *self.Pe_yx     = !Null
            'Pe-yy': *self.Pe_yy     = !Null
            'Pe-yz': *self.Pe_yz     = !Null
            'Pe-zx': *self.Pe_zx     = !Null
            'Pe-zy': *self.Pe_zy     = !Null
            'Pe-zz': *self.Pe_zz     = !Null
            'Pi-xx': *self.Pi_xx     = !Null
            'Pi-xy': *self.Pi_xy     = !Null
            'Pi-xz': *self.Pi_xz     = !Null
            'Pi-yx': *self.Pi_yx     = !Null
            'Pi-yy': *self.Pi_yy     = !Null
            'Pi-yz': *self.Pi_yz     = !Null
            'Pi-zx': *self.Pi_zx     = !Null
            'Pi-zy': *self.Pi_zy     = !Null
            'Pi-zz': *self.Pi_zz     = !Null
            'Uex':   *self.Uex       = !Null
            'Uey':   *self.Uey       = !Null
            'Uez':   *self.Uez       = !Null
            'Uix':   *self.Uix       = !Null
            'Uiy':   *self.Uiy       = !Null
            'Uiz':   *self.Uiz       = !Null
            else: ;Do nothing
        endcase
    endforeach
end


;+
;   Compute the component of a vector parallel to the magnetic field.
;
; :Params:
;       DATA:               in, required, type=NxM float
;                           The data of which the derivative will be taken.
;
; :Keywords:
;       OVERWRITE:          in, optional, type=boolean, default=0
;                           If set, the derivative will overwrite `DATA` and avoids
;                               having an extra copy in memory.
;
; :Returns:
;       DERIVATIVE:         The derivative of `DATA` with respect to X
;-
function MrSim::Parallel, vx, vy, vz, $
FIELD=field, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    ;Defaults
    magnitude = keyword_set(magnitude)
    x         = keyword_set(x)
    y         = keyword_set(y)
    z         = keyword_set(z)
    if magnitude + x + y + z gt 1 then message, 'MAGNITUDE, X, Y, and Z are mutually exclusive.'
    
    _field = n_elements(field) eq 0 ? 'B' : strupcase(field)
    if MrIsMember(['B', 'E'], _field) eq 0 $
        then message, 'FIELD must be either "B" or "E".'
    
    ;Get Data
    x   = self -> GetData(_field + 'x')
    y   = self -> GetData(_field + 'y')
    z   = self -> GetData(_field + 'z')
    mag = sqrt(x^2 + y^2 + z^2)
        
    ;Magnitude of the parallel component
    par = (vx*x + vy*y + vz*z) / temporary(mag)
    
    ;Components
    ;   par_x = par * x_hat
    ;   par_y = par * y_hat
    ;   par_z = par * z_hat
    case 1 of
        magnitude: ;Do nothing. PAR had already been calculated.
        x:         par = par * x / temporary(mag)
        y:         par = par * y / temporary(mag)
        z:         par = par * z / temporary(mag)
        else:      par = [[[par * x / mag], $
                           [par * y / mag], $
                           [par * z / mag]]]
    endcase
    
    ;Return
    return, par
end


;+
;   Compute the component of a vector parallel to the magnetic field.
;
; :Params:
;       DATA:               in, required, type=NxM float
;                           The data of which the derivative will be taken.
;
; :Keywords:
;       OVERWRITE:          in, optional, type=boolean, default=0
;                           If set, the derivative will overwrite `DATA` and avoids
;                               having an extra copy in memory.
;
; :Returns:
;       DERIVATIVE:         The derivative of `DATA` with respect to X
;-
function MrSim::Perpendicular, vx, vy, vz, $
FIELD=field, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    magnitude = keyword_set(magnitude)
    x         = keyword_set(x)
    y         = keyword_set(y)
    z         = keyword_set(z)
    if magnitude + x + y + z gt 1 then message, 'MAGNITUDE, X, Y, and Z are mutually exclusive.'
    
    ;Compute the parallel component
    par_xyz = self -> Parallel(vx, vy, vz, FIELD=field)
    
    ;Compute the perpendicular components
    perp_x = vx - par_xyz[*,*,0]
    perp_y = vy - par_xyz[*,*,1]
    perp_z = vz - par_xyz[*,*,2]
    
    ;Compute the perpendicular component
    case 1 of
        magnitude: perp = sqrt(perp_x^2 + perp_y^2 + perp_z^2)
        x:         perp = temporary(perp_x)
        y:         perp = temporary(perp_y)
        z:         perp = temporary(perp_z)
        else:      perp = [[[temporary(perp_x)], $
                            [temporary(perp_y)], $
                            [temporary(perp_z)]]]
    endcase
    
    return, perp 
end


;+
;   The purpose of this program is to compute the x-derivative of a data array.
;
; :Params:
;       DATA:               in, required, type=NxM float
;                           The data of which the derivative will be taken.
;
; :Keywords:
;       OVERWRITE:          in, optional, type=boolean, default=0
;                           If set, the derivative will overwrite `DATA` and avoids
;                               having an extra copy in memory.
;
; :Returns:
;       DERIVATIVE:         The derivative of `DATA` with respect to X
;-
function MrSim::D_DX, data, $
OVERWRITE=overwrite
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    overwrite = keyword_set(overwrite)

    ;Get the x-size of the simulation in electron skin depths
    self -> GetInfo, DX_DE=dx_de
    dims = size(data, /DIMENSIONS)

    ;Overwrite?
    if overwrite then begin
        ;Take the central difference
        data[1:dims[0]-2,*] = (data[2:dims[0]-1,*] - data[0:dims[0]-3,*]) / (2.0 * dx_de)
        data[0,*]           = 0
        data[dims[0]-1,*]   = 0
        
        return, data
    endif
    
    ;Make a copy?
    derivative = fltarr(dims)
    derivative[1:dims[0]-2,*] = (data[2:dims[0]-1,*] - data[0:dims[0]-3,*]) / (2.0 * dx_de)
    
    return, derivative
end


;+
;   The purpose of this program is to compute the Z-derivative of a data array.
;
; :Params:
;       DATA:               in, required, type=NxM float
;                           The data of which the derivative will be taken.
;
; :Keywords:
;       OVERWRITE:          in, optional, type=boolean, default=0
;                           If set, the derivative will overwrite `DATA` and avoids
;                               having an extra copy in memory.
;
; :Returns:
;       DERIVATIVE:         The derivative of `DATA` with respect to Z
;-
function MrSim::D_DZ, data, $
OVERWRITE=overwrite
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    overwrite = keyword_set(overwrite)

    ;Get the x-size of the simulation in electron skin depths
    self -> GetInfo, DZ_DE=dz_de
    dims = size(data, /DIMENSIONS)

    ;Overwrite?
    if overwrite then begin
        ;Take the central difference
        data[*,1:dims[1]-2] = (data[*,2:dims[1]-1] - data[*,0:dims[1]-3]) / (2.0 * dz_de)
        data[*,0]           = 0
        data[*,dims[1]-1]   = 0
        
        return, data
    endif
    
    ;Make a copy?
    derivative = fltarr(dims)
    derivative[*,1:dims[1]-2] = (data[*,2:dims[1]-1] - data[*,0:dims[1]-3]) / (2.0 * dx_de)
    
    return, derivative
end


;+
;   The purpose of this program is to read electron data from a ".bin" file produced by 
;   one of Bill Daughton's simulation runs.
;
; :Private:
;
; :Params:
;       FILENAME:           in, required, type=string
;                           Name of the file containing particle data.
;
; :Keywords:
;       ENERGY:             in, optional, type=boolean, default=0
;                           If set, momentum will be converted to energy.
;       FMAP_DIR:           in, optional, type=string, default=pwd
;                           Directory in which to find an fMap. See MrSim_Create_fMap.pro.
;       N_RECS_PER_CHUNK:   in, optional, type=ulong64, default=1000000ULL
;                           Number of records per data chunk.
;       REC_SAMPLE:         in, optional, type=any, default=fltarr(5)
;                           An example of how records are to be read. The default assumes
;                               one record consists of five 4-byte (16-bit) floating point
;                               data values of the form [x, z, ux, uy, uz], where x and
;                               z are two spatial locations and u[xyz] are the three
;                               components of the momentum.
;       VELOCITY:           in, optional, type=boolean, default=0
;                           If set, momentum will be converted to velocity.
;       VERBOSE:            in, optional, type=boolean, default=0
;                           If set, information about particle will be printed data to
;                               the command window.
;       XRANGE:             in, required, type=fltarr(2)
;                           X-range (in de) over which particle data is to be kept.
;       ZRANGE:             in, required, type=fltarr(2)
;                           Z-range (in de) over which particle data is to be kept.
;
; :Returns:
;       DATA:               Electron particle data.
;-
function MrSim::ReadElectrons, filename, $
ENERGY=energy, $
FMAP_DIR=fmap_dir, $
N_RECS_PER_CHUNK=n_recs_per_chunk, $
REC_SAMPLE=rec_sample, $
VELOCITY=velocity, $
VERBOSE=verbose, $
XRANGE=xrange, $
ZRANGE=zrange
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(pwd) gt 0 then cd, pwd
        void = cgErrorMsg()
        return, !Null
    endif
    
    ;Get defaults
    MrSim_Which, self.simname, FMAP_DIR=fmap, EFILE=eFile, TINDEX=self.time
    
    ;Set defaults
    if n_elements(filename) eq 0 then filename = eFile
    if n_elements(fmap_dir) eq 0 then fmap_dir = fmap
    if n_elements(xrange)   eq 0 then xrange   = self.xrange
    if n_elements(zrange)   eq 0 then zrange   = self.zrange
    
    ;Read the simulation info file
    data = MrSim_ReadParticles(filename, xrange, zrange, $
                               ENERGY           = energy, $
                               FMAP_DIR         = fMap_dir, $
                               N_RECS_PER_CHUNK = n_rec_per_chunk, $
                               REC_SAMPLE       = rec_sample, $
                               VELOCITY         = velocity, $
                               VERBOSE          = verbose)
    
    return, data
end


;+
;   The purpose of this method is to read the ASCII "info" file relating to Bill
;   Daughton's simulations.
;
; :Private:
;
; :Params:
;       FILENAME:           in, required, type=string
;                           Name of the "info" file to be read.
;
; :Keywords:
;       ASCII_VERSION:  in, optional, type=integer, default=1
;                       Version of the ASCII info file. See MrSim_Which.pro.
;-
pro MrSim::ReadInfo_Ascii, filename, $
VERSION=version
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        void = temporary(*self.info)
        void = cgErrorMsg()
        return
    endif
    
    ;Read the simulation info file
    *self.info = MrSim_ReadInfoAscii(filename, VERSION=version)
end


;+
;   The purpose of this program is to read parameters out of the "info" file relating
;   to the domain size in pixels and in physical units, "de" -- electron skin depth
;
; :Obsolete:
;
; :Private:
;
; :Params:
;       STATUS:                 out, optional, type=int
;                               A status indicator telling if the file was read::
;                                   0   -   File was read without error.
;                                   1   -   An error occured while reading the file.
;-
pro MrSim::ReadInfo_Binary, info_file, status
    compile_opt strictarr

    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        void = cgErrorMSG()
        return
    endif

    status = 0
    if n_elements(info_file) eq 0 $
        then info_file = filepath('info', ROOT_DIR=self.directory)

;---------------------------------------------------------------------
;Check Inputs ////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    if file_test(info_file) eq 0 then begin
        status = 1
        message, 'Cannot find binary info file: "' + info_file + '"'
    endif

    ;Set data type for variables to be read.
    nx = 0L     ;X-Size of simulation in pixels
    ny = 0L     ;Y-Size of simulation in pixels
    nz = 0L     ;Z-Size of simulation in pixels

;---------------------------------------------------------------------
;Try Little Endian ///////////////////////////////////////////////////
;---------------------------------------------------------------------

    ;Try little endian, then switch to big endian if this fails
    little = 0
    on_ioerror, switch_endian
    openr, lun, info_file, /GET_LUN, /F77_UNFORMATTED, /SWAP_IF_BIG_ENDIAN

    readu, lun, nx, ny, nz

    little=1

;---------------------------------------------------------------------
;Try Big Endian //////////////////////////////////////////////////////
;---------------------------------------------------------------------

    switch_endian: if (not little) then begin
        message, " ** Little endian failed --> Switch to big endian", /INFORMATIONAL
        close, lun
        free_lun, lun
        openr, unit, info_file, /GET_LUN, /F77_UNFORMATTED, /SWAP_IF_LITTLE_ENDIAN
        readu, unit, nx, ny, nz
    endif

    ; Read the problem desciption
    on_ioerror, halt_error2
    readu, lun, xmax, ymax, zmax

;    readu,unit,mime,wpewce,rhoi,teti,toutput

    close, lun
    free_lun, lun

    ;Set object properties
    info = { nx: nx, $
             ny: ny, $
             nz: nz, $
             lx_de: xmax, $
             ly_de: ymax, $
             lz_de: zmax $
           }

    ;If the info has already been read, simply assign the proper values. Since
    ;the ASCII info file has more information, we do not want to overwrite it.
    if n_elements(*self.info) gt 0 $
        then struct_assign, info, *self.info, /NOZERO $
        else *self.info = info

    self.little = little

    return

    halt_error2:
        status = 1
        message, 'Info file cannot be read: "' + info_file + '"'
end


;+
;   The purpose of this program is to read data from a ".gda" file produced by 
;   one of Bill Daughton's simulation runs.
;
; :Private:
;
; :Params:
;       FILENAME:           in, required, type=string
;                           Name of the "info" file to be read.
;-
pro MrSim::ReadMoment, filename
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(pwd) gt 0 then cd, pwd
        void = cgErrorMsg()
        return
    endif

    ;Read the simulation info file
    message, 'The ReadMoment method must be overridden by a superclass.'
end


;+
;   The purpose of this program is to return a subset of electrons within the current
;   dataspace.
;
;   The ::GetData method will read, store, and return electron particle data within
;   the dataspace defined by the XRANGE and ZRANGE properties. This method will return
;   a subset of the electrons within this dataspace and avoid, when possible, reading
;   data from files.
;
; :Params:
;       XRANGE:             in, required, type=string, default=
;                           X-range in which electrons are to be found.
;       ZRANGE:             in, required, type=string, default=
;                           Z-range in which electrons are to be found.
;
; :Keywords:
;       ENERGY:             in, optional, type=boolean, default=0
;                           If set, momentum will be converted to energy.
;       FILENAME:           in, optional, type=string
;                           Name of the binary file containing particle data. Used only
;                               if electrons have not already been read.
;       VELOCITY:           in, optional, type=boolean, default=0
;                           If set, momentum will be converted to velocity.
;
; :Returns:
;       E_DATA:             Position and momentum [x, z, ux, uy, uz] for all electrons
;                               within the specified region.
;-
function MrSim::GetElectrons, xrange, zrange, $
ENERGY=energy, $
FILENAME=filename, $
VELOCITY=velocity
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    energy   = keyword_set(energy)
    velocity = keyword_set(velocity)
    if energy + velocity gt 1 then $
        message, 'Keywords ENERGY and VELOCITY are mutually exclutive.'
    if n_elements(filename) eq 0 then filename = ''

    ;Make sure electrons are saved.
    if n_elements(*self.electrons) eq 0 then begin
        message, 'Reading electron particle data...', /INFORMATIONAL
        self -> ReadElectrons, filename, ENERGY=energy, VELOCITY=velocity
        if n_elements(*self.electrons) eq 0 then return, !Null
    endif
    
    ;Pick out the data
    iDist = where((*self.electrons)[0,*] ge xrange[0] and $
                  (*self.electrons)[0,*] le xrange[1] and $
                  (*self.electrons)[1,*] ge zrange[0] and $
                  (*self.electrons)[1,*] le zrange[1], nDist)
    
    ;Ensure electrons were found.
    if nDist eq 0 $
        then message, 'No particles in the spacial range given.' $
        else e_data = (*self.electrons)[*, temporary(iDist)]
    
    return, e_data
end


;+
;   The purpose of this program is to return a subset of data from the current data space.
;
;   The ::GetData method will read, store, and return data within the dataspace defined
;   by the XRANGE and ZRANGE properties. This method will search for a subset of data
;   within this dataspace and will avoid reading data from the GDA files, if possible.
;
; :Params:
;       XRANGE:             in, required, type=string, default=
;                           X-range in which electrons are to be found.
;       ZRANGE:             in, required, type=string, default=
;                           Z-range in which electrons are to be found.
;
; :Keywords:
;       ENERGY:             in, optional, type=boolean, default=0
;                           If set, momentum will be converted to energy.
;       FILENAME:           in, optional, type=string
;                           Name of the binary file containing particle data. Used only
;                               if electrons have not already been read.
;       VELOCITY:           in, optional, type=boolean, default=0
;                           If set, momentum will be converted to velocity.
;
; :Returns:
;       E_DATA:             Position and momentum [x, z, ux, uy, uz] for all electrons
;                               within the specified region.
;-
function MrSim::GetMoment, name, xrange, zrange, $
XLOC=xloc, $
ZLOC=zloc, $
FILENAME=filename
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    if n_elements(filename) eq 0 then filename = ''

    ;Get the data
    data = self -> GetData(name)

    ;Find the index range
    ;   - Using the GetProperty method will select the current subset of the simulation
    ;       coordinates.
    self -> GetProperty, XSIM=xSim, ZSIM=zSim
    ix = MrIndexRange(xSim, xrange, STRIDE=xStride)
    iz = MrIndexRange(zSim, zrange, STRIDE=zStride)

    ;Select the subset of data to be returned
    data = data[ix[0]:ix[1]:xStride, iz[0]:iz[1]:zStride]
    if arg_present(xloc) then xloc = xSim[ix[0]:ix[1]:xStride]
    if arg_present(zloc) then zloc = zSim[iz[0]:iz[1]:zStride]

    return, data
end


;+
;   The purpose of this program is to read data from a ".gda" file produced by 
;   one of Bill Daughton's simulation runs.
;
; :Params:
;       DATA_PRODUCT:       in, required, type=string, default=
;                           The name of the data product to be read. For a list of
;                               available data product, call mr_readSIM without any
;                               arguments.
;
; :Keywords:
;       DX:                 in, optional, type=boolean, default=0
;                           If set, the derivative of `DATA_PRODUCT` with respect to X
;                               will be taken.
;       DZ:                 in, optional, type=boolean, default=0
;                           If set, the derivative of `DATA_PRODUCT` with respect to Z
;                               will be taken.
;
; :Returns:
;       DATA:               The requested data. If the data product does not exist,
;                               then !Null will be returned.
;-
function MrSim::GetData, name, $
DX=dx, $
DZ=dz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    ;If no parameters were given, print a list of available data products.
    if n_params() eq 0 then begin
        MrSim -> listProducts
        return, !Null
    endif
    
    ;Defaults
    dx = keyword_set(dx)
    dz = keyword_set(dz)
    if dx + dz gt 1 then message, 'DX and DZ are mutually exclusive.'

    _name = strupcase(name)
;-------------------------------------------------------
;Had Data been Read? ///////////////////////////////////
;-------------------------------------------------------

    ;If the data has not been read, then try to read it from a file.
    if _name eq 'E-' then $
        if self -> HasData(_name) eq 0 then self -> ReadElectrons

    ;Check to see if the data has already been read first.
    case strupcase(_name) of
        'AY':     data = self -> GetGDA('Ay')
        'BX':     data = self -> GetGDA('Bx')
        'BY':     data = self -> GetGDA('By')
        'BZ':     data = self -> GetGDA('Bz')
        'E-':     *self.electrons
        'EX':     data = self -> GetGDA('Ex')
        'EY':     data = self -> GetGDA('Ey')
        'EZ':     data = self -> GetGDA('Ez')
        'NE':     data = self -> GetGDA('ne')
        'NI':     data = self -> GetGDA('ni')
        'PE-XX':  data = self -> GetGDA('Pe-xx')
        'PE-XY':  data = self -> GetGDA('Pe-xy')
        'PE-XZ':  data = self -> GetGDA('Pe-xz')
        'PE-YX':  data = self -> GetGDA('Pe-xy')
        'PE-YY':  data = self -> GetGDA('Pe-yy')
        'PE-YZ':  data = self -> GetGDA('Pe-yz')
        'PE-ZX':  data = self -> GetGDA('Pe-xz')
        'PE-ZY':  data = self -> GetGDA('Pe-yz')
        'PE-ZZ':  data = self -> GetGDA('Pe-zz')
        'PI-XX':  data = self -> GetGDA('Pi-xx')
        'PI-XY':  data = self -> GetGDA('Pi-xy')
        'PI-XZ':  data = self -> GetGDA('Pi-xz')
        'PI-YX':  data = self -> GetGDA('Pi-xy')
        'PI-YY':  data = self -> GetGDA('Pi-yy')
        'PI-YZ':  data = self -> GetGDA('Pi-yz')
        'PI-ZX':  data = self -> GetGDA('Pi-xz')
        'PI-ZY':  data = self -> GetGDA('Pi-yz')
        'PI-ZZ':  data = self -> GetGDA('Pi-zz')
        'UEX':    data = self -> GetGDA('Uex')
        'UEY':    data = self -> GetGDA('Uey')
        'UEZ':    data = self -> GetGDA('Uez')
        'UIX':    data = self -> GetGDA('Uix')
        'UIY':    data = self -> GetGDA('Uiy')
        'UIZ':    data = self -> GetGDA('Uiz')
        
        ;Custom Data Products
        'A0_E':       data = self -> A0_e()
        'AN_E':       data = self -> An_e()
        'BETA_E':     data = self -> Beta_e()
        'BETA_I':     data = self -> Beta_i()
        'BETA':       data = self -> Beta_p()
        'BMAG':       data = self -> Bmag()
        'DIVPE_MAG':  data = self -> divPe_mag()
        'DIVPE_X':    data = self -> divPe_x()
        'DIVPE_Z':    data = self -> divPe_z()
        'DIVPI_X':    data = self -> divPi_x()
        'DIVPI_Z':    data = self -> divPi_z()
        'DNG_E':      data = self -> Dng_e()
        'E_DOT_J':    data = self -> E_dot_J()
        'EJPAR':      data = self -> EJpar()
        'E_PAR':      data = self -> E_par()
        'EMAG':       data = self -> Emag()
        'EPAR':       data = self -> Epar()
        'EPERP':      data = self -> Eperp()
        'EXJX':       data = self -> ExJx()
        'EYBY':       data = self -> EyBy()
        'EYJY':       data = self -> EyJy()
        'EZJZ':       data = self -> EzJz()
        'FE_MAG':     data = self -> Fe_mag()
        'FEX':        data = self -> Fex()
        'FEY':        data = self -> Fey()
        'FEZ':        data = self -> Fez()
        'FI_MAG':     data = self -> Fi_mag()
        'FIX':        data = self -> Fix()
        'FIY':        data = self -> Fiy()
        'FIZ':        data = self -> Fiz()
        'JE_XZ':      data = self -> Je_xz()
        'JEX':        data = self -> Jex()
        'JEY':        data = self -> Jey()
        'JEZ':        data = self -> Jez()
        'JEXB_X':     data = self -> JexB_x()
        'JEXB_Y':     data = self -> JexB_y()
        'JEXB_Z':     data = self -> JexB_z()
        'JIX':        data = self -> Jix()
        'JIY':        data = self -> Jiy()
        'JIZ':        data = self -> Jiz()
        'JPAR':       data = self -> Jpar()
        'JXB_X':      data = self -> JxB_x()
        'JXB_Y':      data = self -> JxB_y()
        'JXB_Z':      data = self -> JxB_z()
        'JX':         data = self -> Jx()
        'JY':         data = self -> Jy()
        'JZ':         data = self -> Jz()
        'LAMBDA_E':   data = self -> lambda_e()
        'LAMBDA_I':   data = self -> lambda_i()
        'NE_NI':      data = self -> ne_ni()
        'NE_RATIO':   data = self -> ne_ratio()
        'NI_RATIO':   data = self -> ni_ratio()
        'P_TOT':      data = self -> P_tot()
        'PE_PAR':     data = self -> Pe_par()
        'PE_PERP':    data = self -> Pe_perp()
        'PB':         data = self -> Pb()
        'PBE':        data = self -> Pbe()
        'PBI':        data = self -> Pbi()
        'PE':         data = self -> Pe()
        'PI':         data = self -> Pi()
        'UE_MAG':     data = self -> Ue_mag()
        'UE_PAR':     data = self -> Ue_par()
        'UE_PERP':    data = self -> Ue_perp()
        'UEXB_X':     data = self -> UexB_x()
        'UEXB_Y':     data = self -> UexB_y()
        'UEXB_Z':     data = self -> UexB_z()
        'UI_MAG':     data = self -> Ui_mag()
        'UI_PAR':     data = self -> Ui_par()
        'UI_PERP':    data = self -> Ui_perp()
        'UIXB_MAG':   data = self -> UixB_mag()
        'UIXB_X':     data = self -> UixB_x()
        'UIXB_Y':     data = self -> UixB_y()
        'UIXB_Z':     data = self -> UixB_z()
        'VEXB_MAG':   data = self -> vExB_mag()
        'VEXB_X':     data = self -> vExB_x()
        'VEXB_Y':     data = self -> vExB_y()
        'VEXB_Z':     data = self -> vExB_z()
        'VXB_X':      data = self -> VxB_x()
        'VXB_Y':      data = self -> VxB_y()
        'VXB_Z':      data = self -> VxB_z()
        'VX':         data = self -> Vx()
        'VY':         data = self -> Vy()
        'VZ':         data = self -> Vz()
        'V_AE':       data = self -> v_Ae()
        'V_A':        data = self -> v_Ae()
        'W_PE':       data = self -> w_pe()
        'W_PI':       data = self -> w_pi()
        'W_CE':       data = self -> w_ce()
        'W_CI':       data = self -> w_ci()
        else: message, 'Data product not available: "' + data_product + '".'
    endcase
    
    ;Take the derivative?
    if dx then data = self -> D_DX(data, /OVERWRITE)
    if dz then data = self -> D_DZ(data, /OVERWRITE)
    
    return, data
end


;+
;   The purpose of this program is to read data from a ".gda" file produced by 
;   one of Bill Daughton's simulation runs.
;
; :Params:
;       DATA_PRODUCT:       in, required, type=string, default=
;                           The name of the data product to be read. For a list of
;                               available data product, call mr_readSIM without any
;                               arguments.
;
; :Keywords:
;       DX:                 in, optional, type=boolean, default=0
;                           If set, the derivative of `DATA_PRODUCT` with respect to X
;                               will be taken.
;       DZ:                 in, optional, type=boolean, default=0
;                           If set, the derivative of `DATA_PRODUCT` with respect to Z
;                               will be taken.
;
; :Returns:
;       DATA:               The requested data. If the data product does not exist,
;                               then !Null will be returned.
;-
function MrSim::GetGDA, name
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    ;Convert name from current coordinate system to simulation coordinates
    simName = MrSim_GDA_Rename(name, self.coord_system)

    ;Read the data
    if self -> HasData(simName) eq 0 then self -> ReadData, simName
    
    ;Return the data
    case strupcase(simName) of
        'AY':    return, *self.Ay
        'BX':    return, *self.Bx
        'BY':    return, *self.By
        'BZ':    return, *self.Bz
        'EX':    return, *self.Ex
        'EY':    return, *self.Ey
        'EZ':    return, *self.Ez
        'NE':    return, *self.n_e
        'NI':    return, *self.n_i
        'PE-XX': return, *self.Pe_xx
        'PE-XY': return, *self.Pe_xy
        'PE-XZ': return, *self.Pe_xz
        'PE-YX': return, *self.Pe_xy
        'PE-YY': return, *self.Pe_yy
        'PE-YZ': return, *self.Pe_yz
        'PE-ZX': return, *self.Pe_xz
        'PE-ZY': return, *self.Pe_yz
        'PE-ZZ': return, *self.Pe_zz
        'PI-XX': return, *self.Pi_xx
        'PI-XY': return, *self.Pi_xy
        'PI-XZ': return, *self.Pi_xz
        'PI-YX': return, *self.Pi_xy
        'PI-YY': return, *self.Pi_yy
        'PI-YZ': return, *self.Pi_yz
        'PI-ZX': return, *self.Pi_xz
        'PI-ZY': return, *self.Pi_yz
        'PI-ZZ': return, *self.Pi_zz
        'UEX':   return, *self.Uex
        'UEY':   return, *self.Uey
        'UEZ':   return, *self.Uez
        'UIX':   return, *self.Uix
        'UIY':   return, *self.Uiy
        'UIZ':   return, *self.Uiz
        else: message, 'Not a valid GDA product: "' + name + '".'
    endcase
end


;+
;   The purpose of this method is to retrieve data from the info file.
;
; :Keywords:
;       ECOUNTFACTOR:       out, optional, type=long
;                           Bill Daughton does not save every particle, only every-other.
;                               When dealing with electron counts, you must multiply by
;                               this factor to obtain the actual total counts.
;       L_DI:               out, optional, type=double
;                           Harris sheet thickness over ion inertial length
;       L_DE:               out, optional, type=double
;                           Number of electron inertial lengths per cell
;       TI_TE:              out, optional, type=double
;                           Ion to electron temperature ratio
;       TBI_TBE:            out, optional, type=double
;                           Ion to electron background temperature ratio
;       TBE_TE:             out, optional, type=double
;                           Background to Harris electron temperature ratio
;       WPE_WCE:            out, optional, type=double
;                           Electron plasma to cyclotron frequency ratio
;       MI_ME:              out, optional, type=double
;                           Ion to electron mass ratio
;       THETA:              out, optional, type=double
;                           Angle between background field and Bx
;       TAUI:               out, optional, type=double
;                           Time to run simulation, in wci's
;       NUM_STEP:           out, optional, type=double
;                           ???
;       LX_DE:              out, optional, type=double
;                           X-size of simulation, in electron skin depths
;       LY_DE:              out, optional, type=double
;                           Y-size of simulation, in electron skin depths
;       LZ_DE:              out, optional, type=double
;                           Z-size of simulation, in electron skin depths
;       LX_DI:              out, optional, type=double
;                           X-size of simulation, in ion skin depths
;       LY_DI:              out, optional, type=double
;                           Y-size of simulation, in ion skin depths
;       LZ_DI:              out, optional, type=double
;                           Z-size of simulation, in ion skin depths
;       NX:                 out, optional, type=double
;                           Number if cells in the x-direction
;       NY:                 out, optional, type=double
;                           Number of cells in the y-direction
;       NZ:                 out, optional, type=double
;                           Number of cells in the z-direction
;-
pro MrSim::GetInfo, $
;Derive Info
UNITS=units, $
ECOUNTFACTOR=eCountFactor, $

;From INFO file
B0=b0, $
B2_B1=B2_B1, $
COURANT=courant, $
DAMP=damp, $
DI=di, $
DTXWCE=dtxwce, $
DTXWCI=dtxwci, $
DTXWPE=dtxwpe, $
DX_RHOI=dx_rhoi, $
DX_RHOE=dx_rhoe, $
DX_DBYE=dx_dbye, $
DX_DE=dx_de, $
DY_DE=dy_de, $
DZ_DE=dz_de, $
E_INTERVAL=E_interval, $
E_WEIGHT_SHEET=e_weight_sheet, $
E_WEIGHT_BCKGRND=e_weight_bckgrnd, $
FUNCT_DIAGNOSTICS_INTERVAL=funct_diagnostics_interval, $
ION_WEIGHT_SHEET=ion_weight_sheet, $
ION_WEIGHT_BCKGRND=ion_weight_bckgrnd, $
L_DEBYE=L_debye, $
L_DE=L_de, $
L_DI=L_di, $
LX_DE=Lx_de, $
LX_DI=Lx_di, $
LY_DE=Ly_de, $
LY_DI=Ly_di, $
LZ_DE=Lz_de, $
LZ_DI=Lz_di, $
MI_ME=mi_me, $
N_E=N_e, $
N_TOTAL=N_total, $
N0=n0, $
N2_N1=n2_n1, $
NE_SHEET=Ne_sheet, $
NE_BACK=Ne_back, $
NFAC=nfac, $
NPROC=nproc, $
NPPC=nppc, $
NUM_STEP=num_step, $
NX=nx, $
NY=ny, $
NZ=nz, $
TAUI=taui, $
TBI_TI=Tbi_Ti, $
TBE_TE=Tbe_Te, $
THETA=theta, $
TI_TE=Ti_Te, $
V_A=v_A, $
VTHI_C=vthi_c, $
VTHE_C=vthe_c, $
VDRI_C=vdri_c, $
VDRE_C=vdre_c, $
WPE_WCE=wpe_wce
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Derived INFO
    if arg_present(units)        then if self.ion_scale then units = 'di' else units = 'de'
    if arg_present(eCountFactor) then MrSim_Which, self.simnum, TINDEX=self.time, ECOUNTFACTOR=eCountFactor
    
    ;Make sure the info file has been read
    if n_elements(*self.info) eq 0 then return
    
    ;
    ; Binary Info file
    ;   - The parameters are also in the ASCII info file
    ;
    if arg_present(nx)         then nx         = (*self.info).nx
    if arg_present(ny)         then ny         = (*self.info).ny
    if arg_present(nz)         then nz         = (*self.info).nz
    if arg_present(Lx_de)      then Lx_de      = (*self.info).Lx_de
    if arg_present(Ly_de)      then Ly_de      = (*self.info).Ly_de
    if arg_present(Lz_de)      then Lz_de      = (*self.info).Lz_de
    
    ;
    ; ASCII Info file
    ;   - Make sure it has been read.
    ;
    if has_tag(*self.info, 'MI_ME') eq 0 then return
    if arg_present(L_di)       then L_di       = (*self.info).L_di
    if arg_present(L_de)       then L_de       = (*self.info).L_de
    if arg_present(Ti_Te)      then Ti_Te      = (*self.info).Ti_Te
    if arg_present(Tbi_Ti)     then Tbi_Ti     = (*self.info).Tbi_Ti
    if arg_present(Tbe_Te)     then Tbe_Te     = (*self.info).Tbe_Te
    if arg_present(wpe_wce)    then wpe_wce    = (*self.info).wpe_wce
    if arg_present(mi_me)      then mi_me      = (*self.info).mi_me
    if arg_present(theta)      then theta      = (*self.info).theta
    if arg_present(taui)       then taui       = (*self.info).taui
    if arg_present(num_step)   then num_step   = (*self.info).num_step
    if arg_present(Lx_di)      then Lx_di      = (*self.info).Lx_di
    if arg_present(Ly_di)      then Ly_di      = (*self.info).Ly_di
    if arg_present(Lz_di)      then Lz_di      = (*self.info).Lz_di
    if arg_present(damp)       then damp       = (*self.info).damp
    if arg_present(courant)    then courant    = (*self.info).courant
    if arg_present(nproc)      then nproc      = (*self.info).nproc
    if arg_present(nppc)       then nppc       = (*self.info).nppc
    if arg_present(b0)         then b0         = (*self.info).b0
    if arg_present(v_A)        then v_A        = (*self.info).v_A
    if arg_present(di)         then di         = (*self.info).di
    if arg_present(N_e)        then N_e        = (*self.info).N_e
    if arg_present(Ne_sheet)   then Ne_sheet   = (*self.info).Ne_sheet
    if arg_present(Ne_back)    then Ne_back    = (*self.info).Ne_back
    if arg_present(N_total)    then N_total    = (*self.info).N_total
    if arg_present(dtxwpe)     then dtxwpe     = (*self.info).dtxwpe
    if arg_present(dtxwce)     then dtxwce     = (*self.info).dtxwce
    if arg_present(dtxwci)     then dtxwci     = (*self.info).dtxwci
    if arg_present(E_interval) then E_interval = (*self.info).E_interval
    if arg_present(dx_de)      then dx_de      = (*self.info).dx_de
    if arg_present(dy_de)      then dy_de      = (*self.info).dy_de
    if arg_present(dz_de)      then dz_de      = (*self.info).dz_de
    if arg_present(L_debye)    then L_debye    = (*self.info).L_debye
    if arg_present(dx_rhoi)    then dx_rhoi    = (*self.info).dx_rhoi
    if arg_present(dx_rhoe)    then dx_rhoe    = (*self.info).dx_rhoe
    if arg_present(dx_debye)   then dx_debye   = (*self.info).dx_debye
    if arg_present(n0)         then n0         = (*self.info).n0
    if arg_present(vthi_c)     then vthi_c     = (*self.info).vthi_c
    if arg_present(vthe_c)     then vthe_c     = (*self.info).vthe_c
    if arg_present(vdri_c)     then vdri_c     = (*self.info).vdri_c
    if arg_present(vdre_c)     then vdre_c     = (*self.info).vdre_c
    if arg_present(nfac)       then nfac       = (*self.info).nfac
    if arg_present(e_weight_sheet)     then e_weight_sheet     = (*self.info).e_weight_sheet
    if arg_present(e_weight_bckgrnd)   then e_weight_bckgrnd   = (*self.info).e_weight_bckgrnd
    if arg_present(ion_weight_sheet)   then ion_weight_sheet   = (*self.info).ion_weight_sheet
    if arg_present(ion_weight_bckgrnd) then ion_weight_bckgrnd = (*self.info).ion_weight_bckgrnd
    if arg_present(funct_diagnostics_interval) then funct_diagnostics_interval = (*self.info).funct_diagnostics_interval
end


;+
;   The purpose of this method is to retrieve the specified object properties.
;
; :Keywords:
;       DIRECTORY:          out, optional, type=string
;                           The directory in which search for simulation data.
;-
pro MrSim::GetProperty, $
AXIS_LABELS = axis_labels, $
COORD_SYSTEM = coord_system, $
DIRECTORY = directory, $
ION_SCALE = ion_scale, $
MVA_FRAME = mva_frame, $
ORIENTATION = orientation, $
SIMNAME = simname, $
SIMNUM = simnum, $
TIME = time, $

;SIMULATION DOMAIN
IXRANGE = ixrange, $
IYRANGE = iyrange, $
IZRANGE = izrange, $
XRANGE = xrange, $
XSIM = XSim, $
YRANGE = yrange, $
YSLICE = yslice, $
YSIM = YSim, $
ZRANGE = zrange, $
ZSIM = ZSim
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    if arg_present(axis_labels)  then axis_labels  = self.axis_labels
    if arg_present(coord_system) then coord_system = self.coord_system
    if arg_present(directory)    then directory    = self.directory
    if arg_present(ion_scale)    then ion_scale    = self.ion_scale
    if arg_present(ixrange)      then ixrange      = MrIndexRange(*self.XSim, self.xrange)
    if arg_present(iyrange)      then iyrange      = MrIndexRange(*self.YSim, self.yrange)
    if arg_present(izrange)      then izrange      = MrIndexRange(*self.ZSim, self.zrange)
    if arg_present(mva_frame)    then mva_frame    = self.mva_frame
    if arg_present(orientation)  then orientation  = self.orientation
    if arg_present(simname)      then simname      = self.simname
    if arg_present(simnum)       then simnum       = self.simnum
    if arg_present(time)         then time         = self.time
    if arg_present(xrange)       then xrange       = self.xrange
    if arg_present(yrange)       then yrange       = self.yrange
    if arg_present(yslice)       then yslice       = (*self.YSim)[self.yslice]
    if arg_present(zrange)       then zrange       = self.zrange

    ;Simulation Domain
    if arg_present(XSim) && n_elements(*self.XSim) gt 0 then begin
        ix = MrIndexRange(*self.XSim, self.xrange)
        XSim = (*self.XSim)[ix[0]:ix[1]]
    endif
    
    ;2D simulations have only 1 grid cell in Y
    if arg_present(YSim) then begin
        if n_elements(*self.YSim) gt 1 then begin
            iy = MrIndexRange(*self.YSim, self.yrange)
            YSim = (*self.YSim)[iy[0]:iy[1]]
        endif else begin
            YSim = (*self.YSim)
        endelse
    endif

    if arg_present(ZSim) && n_elements(*self.ZSim) gt 0 then begin
        iz = MrIndexRange(*self.ZSim, self.zrange)
        ZSim = (*self.ZSim)[iz[0]:iz[1]]
    endif
end


;+
;   Given a particular slice/index, determine the location in units of de (di if
;   the ION_SCALE property is set).
;
; :Params:
;       CELL:           in, required, type=integer
;                       Index (or grid-cell) at which the location in "de" is to be
;                           returned.
;
; :Keywords:
;       X:              in, optional, type=boolean, default=0
;                       If set, `CELL` reference to a grid-cell along the x-dimension.
;                           This is the default if no keywords are set.
;       Y:              in, optional, type=boolean, default=0
;                       If set, `CELL` reference to a grid-cell along the y-dimension
;       Z:              in, optional, type=boolean, default=0
;                       If set, `CELL` reference to a grid-cell along the z-dimension
;
; :Returns:
;       COORD:          Location in de (di if the ION_SCALE property is set) of the
;                           given grid cell along the chosen dimension.
;-
function MrSim::GetCoord, cell, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    on_error, 2
    
    ;Defaults
    x = keyword_set(x)
    y = keyword_set(y)
    z = keyword_set(z)
    if x + y + z eq 0 then x = 1
    if x + y + z gt 1 then message, 'X, Y, and Z are mutually exclusive.'
    
    ;Get the location
    case 1 of
        x: coord = (*self.xSim)[cell]
        y: coord = (*self.ySim)[cell]
        z: coord = (*self.zSim)[cell]
    endcase
    
    return, coord
end


;+
;   Given a particular location units of de (di of the ION_SCALE property is set) determine
;   the equivalent grid cell. For a grid cell spanning [min, max] de, the provided
;   location is rounded down and matched to min.
;
; :Params:
;       COORD:          in, required, type=integer
;                       Location at which the grid cell is to be determined. Units are
;                           de unless ION_SCALE is specified.
;
; :Keywords:
;       X:              in, optional, type=boolean, default=0
;                       If set, `COORD` is an x-coordinate.
;                           This is the default if no keywords are set.
;       Y:              in, optional, type=boolean, default=0
;                       If set, `COORD` is an y-coordinate
;       Z:              in, optional, type=boolean, default=0
;                       If set, `COORD` is an z-coordinate
;
; :Returns:
;       CELL:           Grid cell/index corresponding to the given location along the
;                           desired dimension.
;-
function MrSim::GetCell, coord, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    on_error, 2
    
    ;Defaults
    x = keyword_set(x)
    y = keyword_set(y)
    z = keyword_set(z)
    if x + y + z eq 0 then x = 1
    if x + y + z gt 1 then message, 'X, Y, and Z are mutually exclusive.'
    
    ;Get the location
    case 1 of
        x: cell = value_locate(*self.xSim, coord)
        y: cell = n_elements(*self.ySim) eq 1 ? 0 : value_locate(*self.ySim, coord)
        z: cell = value_locate(*self.zSim, coord)
    endcase
    
    return, cell
end


;+
;   Convert time, normalized to the ion gyroperiod, to time indices into the .gda data
;   files.
;
; :Params:
;       TXWCI:          in, required, type=integer
;                       Time, normalized by the ion cyclotron frequency.
;
; :Returns:
;       TINDEX:         Time index into the .gda files.
;-
function MrSim::GetTIndex, txwci
    compile_opt strictarr
    on_error, 2
    
    ;Get the interval between saves.
    self -> GetInfo, DTXWCI=dtxwci
    if n_elements(dtxwci) eq 0 then $
        message, 'Cannot convert: dt*wci unknown. Read ASCII info file and see MrSim_Which.'

    ;Create a long integer of the time index.
    tIndex = fix(txwci / dtxwci, TYPE=3)
    
    return, tIndex
end


;+
;   Convert time indices into the .gda data files to time, normalized to the ion
;   gyroperiod.
;
; :Params:
;       TINDEX:         in, required, type=integer
;                       Time index into the .gda files.
;
; :Returns:
;       TXWCI:          Time, normalized by the ion cyclotron frequency.
;-
function MrSim::GetTxWci, tIndex
    compile_opt strictarr
    on_error, 2
    
    ;Get the interval between saves.
    self -> GetInfo, DTXWCI=dtxwci
    if n_elements(dtxwci) eq 0 then $
        message, 'Cannot convert: dt*wci unknown. Read ASCII info file and see MrSim_Which.'

    ;Create a long integer of the time index.
    txwci = tIndex * dtxwci
    
    return, txwci
end


;+
;   The purpose of this program is to read data from a ".gda" file produced by 
;   one of Bill Daughton's simulation runs.
;
; :Private:
;
; :Params:
;       DATA_PRODUCT:           in, required, type=string, default=
;                               The name of the data product to be read. For a list of
;                                   available data product, call mr_readSIM without any
;                                   arguments.
;
; :Returns:
;       TF_HAS:                  Returns 1 (true) if the requested data product has
;                                   already been read, and 0 (false) otherwise. If the
;                                   data product does not exist as a '.gda' data file,
;                                   then -1 is returned.
;-
function MrSim::HasData, data_product
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

;-------------------------------------------------------
;Had Data been Read? ///////////////////////////////////
;-------------------------------------------------------
    tf_has = 0B
    
    ;Check to see if the data has already been read first.
    case strupcase(data_product) of
        'AY':    if n_elements(*self.Ay)        gt 0 then tf_has = 1B
        'BX':    if n_elements(*self.Bx)        gt 0 then tf_has = 1B
        'BY':    if n_elements(*self.By)        gt 0 then tf_has = 1B
        'BZ':    if n_elements(*self.Bz)        gt 0 then tf_has = 1B
        'E-':    if n_elements(*self.electrons) gt 0 then tf_has = 1B
        'EX':    if n_elements(*self.Ex)        gt 0 then tf_has = 1B
        'EY':    if n_elements(*self.Ey)        gt 0 then tf_has = 1B
        'EZ':    if n_elements(*self.Ez)        gt 0 then tf_has = 1B
        'NE':    if n_elements(*self.n_e)       gt 0 then tf_has = 1B
        'NI':    if n_elements(*self.n_i)       gt 0 then tf_has = 1B
        'PE-XX': if n_elements(*self.Pe_xx)     gt 0 then tf_has = 1B
        'PE-XY': if n_elements(*self.Pe_xy)     gt 0 then tf_has = 1B
        'PE-XZ': if n_elements(*self.Pe_xz)     gt 0 then tf_has = 1B
        'PE-YX': if n_elements(*self.Pe_xy)     gt 0 then tf_has = 1B
        'PE-YY': if n_elements(*self.Pe_yy)     gt 0 then tf_has = 1B
        'PE-YZ': if n_elements(*self.Pe_yz)     gt 0 then tf_has = 1B
        'PE-ZX': if n_elements(*self.Pe_xz)     gt 0 then tf_has = 1B
        'PE-ZY': if n_elements(*self.Pe_yz)     gt 0 then tf_has = 1B
        'PE-ZZ': if n_elements(*self.Pe_zz)     gt 0 then tf_has = 1B
        'PI-XX': if n_elements(*self.Pi_xx)     gt 0 then tf_has = 1B
        'PI-XY': if n_elements(*self.Pi_xy)     gt 0 then tf_has = 1B
        'PI-XZ': if n_elements(*self.Pi_xz)     gt 0 then tf_has = 1B
        'PI-YX': if n_elements(*self.Pi_xy)     gt 0 then tf_has = 1B
        'PI-YY': if n_elements(*self.Pi_yy)     gt 0 then tf_has = 1B
        'PI-YZ': if n_elements(*self.Pi_yz)     gt 0 then tf_has = 1B
        'PI-ZX': if n_elements(*self.Pi_xz)     gt 0 then tf_has = 1B
        'PI-ZY': if n_elements(*self.Pi_yz)     gt 0 then tf_has = 1B
        'PI-ZZ': if n_elements(*self.Pi_zz)     gt 0 then tf_has = 1B
        'UEX':   if n_elements(*self.Uex)       gt 0 then tf_has = 1B
        'UEY':   if n_elements(*self.Uey)       gt 0 then tf_has = 1B
        'UEZ':   if n_elements(*self.Uez)       gt 0 then tf_has = 1B
        'UIX':   if n_elements(*self.Uix)       gt 0 then tf_has = 1B
        'UIY':   if n_elements(*self.Uiy)       gt 0 then tf_has = 1B
        'UIZ':   if n_elements(*self.Uiz)       gt 0 then tf_has = 1B
        else: tf_has = -1
    endcase
    
    return, tf_has
end


;+
;   The purpose of this program is to recreate the simulation domain so that a map
;   can be made between physical units and pixel locations.
;
; :Private:
;-
pro MrSim::MakeSimDomain
    compile_opt strictarr

    ;Get the simulations size in pixels and in electron skin depths.
    self -> GetInfo, LX_DE=xsize, LY_DE=ysize, LZ_DE=zsize, NX=nx, NY=ny, NZ=nz, MI_ME=mi_me

    ;Convert to units of "di"
    if self.ion_scale then begin
        xsize /= sqrt(mi_me)
        ysize /= sqrt(mi_me)
        zsize /= sqrt(mi_me)
    endif

    ;Set the simulation domain    
    case strupcase(self.coord_system) of
        'SIMULATION': begin
            *self.XSim = linspace(0, xsize, nx)
            *self.YSim = linspace(0, ysize, ny) - ysize/2.0
            *self.ZSim = linspace(0, zsize, nz) - zsize/2.0
        endcase
        
        'MAGNETOPAUSE': begin
            *self.ZSim =   linspace(0, xsize, nx)
            *self.YSim =   linspace(0, ysize, ny) - ysize/2.0
            *self.XSim = -(linspace(0, zsize, nz) - zsize/2.0)
        endcase
        
        'MAGNETOTAIL': begin
            *self.XSim = -linspace(0, xsize, nx)
            *self.YSim =  linspace(0, ysize, ny) - ysize/2.0
            *self.ZSim =  linspace(0, zsize, nz) - zsize/2.0
        endcase
        
        else: message, 'Coordinate system unknown: "' + self.coord_system + '".'
    endcase
end


;+
;   The purpose of this program is to read data from a ".gda" file. It must be
;   over-ridden and, at the end, must store the data via the SetData method.
;
; :Private:
;
; :Params:
;       NAME:                   in, required, type=string, default=
;                               The name of the ".gda" data file (without the ".gda"
;                                   file extension). GDA files are typically named after
;                                   the parameter whose data they contain.
;-
pro MrSim::ReadData, name
    compile_opt strictarr

    ;Catch any errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif

    message, 'This method must be over-ridden by a subclass object.'
end


;+
;   The purpose of this method is to store the data in the object.
;
; :Private:
;
; :Keywods:
;       NAME:           in, required, type=string, default=
;                       The name of the ".gda" data file (without the ".gda"
;                           file extension). GDA files are typically named after
;                           the parameter whose data they contain.
;       DATA:           in, required, type=any
;                       The data to be stored.
;-
pro MrSim::SetData, name, data
    compile_opt strictarr

    ;Catch any errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Store the data in the proper location
    case strupcase(name) of
        'AY':    *self.Ay = data
        'BX':    *self.Bx = data
        'BY':    *self.By = data
        'BZ':    *self.Bz = data
        'EX':    *self.Ex = data
        'EY':    *self.Ey = data
        'EZ':    *self.Ez = data
        'NE':    *self.n_e = data
        'NI':    *self.n_i = data
        'PE-XX': *self.Pe_xx = data
        'PE-XY': *self.Pe_xy = data
        'PE-XZ': *self.Pe_xz = data
        'PE-YX': *self.Pe_xy = data
        'PE-YY': *self.Pe_yy = data
        'PE-YZ': *self.Pe_yz = data
        'PE-ZX': *self.Pe_xz = data
        'PE-ZY': *self.Pe_yz = data
        'PE-ZZ': *self.Pe_zz = data
        'PI-XX': *self.Pi_xx = data
        'PI-XY': *self.Pi_xy = data
        'PI-XZ': *self.Pi_xz = data
        'PI-YX': *self.Pi_xy = data
        'PI-YY': *self.Pi_yy = data
        'PI-YZ': *self.Pi_yz = data
        'PI-ZX': *self.Pi_xz = data
        'PI-ZY': *self.Pi_yz = data
        'PI-ZZ': *self.Pi_zz = data
        'UEX':   *self.Uex = data
        'UEY':   *self.Uey = data
        'UEZ':   *self.Uez = data
        'UIX':   *self.Uix = data
        'UIY':   *self.Uiy = data
        'UIZ':   *self.Uiz = data
        else: message, 'Data "' + name + '" is cannot be set.'
    endcase
end


;+
;   The purpose of this method is to set object properties.
;
;   NOTE:
;       More than likely, you will have to set the TIME, [XYZ]Range and possibly the
;       YSLICE properties after switching simulations.
;
; :Params:
;       SIM_ID:         in, optional, type=string/integer
;                       Either the name or the number of the simulation. If not given, a
;                           list of simulations is printed to the command window.
;
; :Keywords:
;       BINARY:         in, optional, type=boolean, default=0
;                       If set, `INFO_FILE` points to the binary info file.
;       DIRECTORY:      in, optional, type=string, default=pwd
;                       Directory in which to find the ".gda" data.
;       INFO_FILE:      in, optional, type=string, default=`DIRECTORY`/../info`
;                       The ASCII info file containing information about the simulation
;                           setup. If `BINARY` is set, the default file will be
;                           `DIRECTORY`/info.
;-
pro MrSim::SetSim, sim_id, $
BINARY=binary, $
DIRECTORY=directory, $
INFO_FILE=info_file
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;If no simulation was given, print info to command window.
    if n_elements(sim_id) eq 0 then begin
        MrSim_Which
        return
    endif
    
    ;Get the simulation name, number, etc.
    MrSim_Which, sim_id, NAME=name, NUMBER=number, $
                         ASCII_INFO=ascii_info, BINARY_INFO=binary_info, DIRECTORY=dir
    
    ;Is the simulation different?
    if number eq self.simnum then return
    
    ;Defaults
    binary = keyword_set(binary)
    if n_elements(directory) eq 0 then directory = dir
    if n_elements(info_file) eq 0 then begin
        if binary $ 
            then info_file = binary_info $
            else info_file = ascii_info
    endif
    
    ;Make sure the info file exists.
    if file_test(info_file) eq 0 then $
        message, 'Cannot find info file: "' + info_file + '"'
    
    ;Clear all data
    self -> Clear_Data
    void = temporary(*self.info)
    void = temporary(*self.xSim)
    void = temporary(*self.ySim)
    void = temporary(*self.zSim)
    void = -1
    
    ;Binary info file
    if binary then begin
        if keyword_set(ion_scale) then $
            message, 'ION_SCALE is only possible with ASCII info file. ' + $
                     'Setting ION_SCALE=0', /INFORMATIONAL
        ion_scale = 0
        self -> ReadInfo_Binary, info_file
        
    ;Ascii info file
    endif else begin
        self -> ReadInfo_Ascii, info_file
    endelse
    
    ;Create the simulation domain
    self -> MakeSimDomain
    
    ;Properties
    self.simname   = name
    self.simnum    = number
    self.directory = directory
end


;+
;   The purpose of this method is to set object properties.
;
; :Keywords:
;       ION_SCALE:      in, optional, type=boolean, default=0
;                       If set, the unit of measurement will be the ion skin depth (di).
;                           The default is to use the electron skin depth (de).
;-
pro MrSim::SetScale, $
ION_SCALE=ion_scale
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ion_scale = keyword_set(ion_scale)
    if ion_scale eq self.ion_scale then return
    
    ;Set the scale size and get the mass ratio    
    self.ion_scale = ion_scale
    self -> GetInfo, MI_ME=mi_me
    if n_elements(mi_me) eq 0 then begin
        message, 'MI_ME not defined. Cannot change scale.', /INFORMATIONAL
        self.ion_scale = 0
        return
    endif

    ;Convert to "di" units
    if ion_scale then begin
        self -> MakeSimDomain
        self.xrange /= sqrt(mi_me)
        self.yrange /= sqrt(mi_me)
        self.zrange /= sqrt(mi_me)
        
    ;Convert to "de" units
    endif else begin
        self -> MakeSimDomain
        self.xrange *= sqrt(mi_me)
        self.yrange *= sqrt(mi_me)
        self.zrange *= sqrt(mi_me)
    endelse
end


;+
;   The purpose of this method is to set object properties.
;
; :Params:
;       DIRECTORY:              out, optional, type=string
;                               The directory in which search for simulation data.
;       COORD_SYSTEM:           in, optional, type=string
;                               Coordinate system in which to display the data. Options are::
;                                   'SIMULATION'
;                                   'MAGNETOPAUSE'
;                                   'MAGNETOTAIL'
;       XRANGE:                 in, optional, type=fltarr(2)
;                               The x-range, in physical units specified by the `ION_SCALE`,
;                                   over which to retrieve data from the simulation files.
;                                   When the x-range is changed, all internal data is purged.
;       XRANGE:                 in, optional, type=fltarr(2)
;                               The y-range, in physical units specified by the `ION_SCALE`,
;                                   over which to retrieve data from the simulation files.
;                                   When the y-range is changed, all internal data is purged.
;       ZRANGE:                 in, optional, type=fltarr(2)
;                               The z-range, in physical units specified by the `ION_SCALE`,
;                                   over which to retrieve data from the simulation files.
;                                   When the z-range is changed, all internal data is purged.
;-
pro MrSim::SetProperty, $
AXIS_LABELS = axis_labels, $
DIRECTORY = directory, $
COORD_SYSTEM = coord_system, $
MVA_FRAME = mva_frame, $
NSMOOTH = nsmooth, $
ORIENTATION = orientation, $
TIME = time, $
XRANGE = xrange, $
YRANGE = yrange, $
YSLICE = yslice, $
ZRANGE = zrange
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    if n_elements(axis_labels) gt 0 then self.axis_labels = axis_labels
    if n_elements(directory)   gt 0 then self.directory   = directory
    if n_elements(nsmooth)     gt 0 then self.nsmooth     = nsmooth
    if n_elements(time)        gt 0 then self.time        = time
    if n_elements(xrange)      gt 0 then self.xrange      = xrange
    if n_elements(yrange)      gt 0 then self.yrange      = yrange
    if n_elements(yslice)      gt 0 then self.yslice      = yslice
    if n_elements(zrange)      gt 0 then self.zrange      = zrange

    ;ORIENTATION
    if n_elements(orientation) gt 0 then begin
        _orientation = strupcase(orientation)
        planes = ['XY', 'XZ', 'YZ']
        if max(_orientation eq planes) eq 0 $
            then message, 'Orientation "' + orientation + '" not recognized.', /INFORMATIONAL
        
        ;Set the orientation    
        self.orientation = orientation
    endif

    ;COORDINATE SYSTEM
    if n_elements(coord_system) gt 0 then begin
        _coord_system = strupcase(coord_system)
        old_sys       = strupcase(self.coord_system)

        ;Translate range from current system to new system
        case old_sys of
            'SIMULATION': begin
                case _coord_system of
                    'SIMULATION':   ;Do nothing
                    'MAGNETOPAUSE': begin
                        xrange = reverse(self.zrange)
                        zrange = self.xrange
                    endcase
                    'MAGNETOTAIL': begin
                        xrange = -self.xrange
                        zrange =  self.zrange
                    endcase
                endcase
            endcase
            
            'MAGNETOPAUSE': begin
                case _coord_system of
                    'SIMULATION': begin
                        xrange = self.zrange
                        zrange = reverse(self.xrange)
                    endcase
                    'MAGNETOPAUSE': ;Do nothing
                    'MAGNETOTAIL': begin
                        xrange = -self.zrange
                        zrange = reverse(self.xrange)
                    endcase
                endcase
            endcase
            
            'MAGNETOTAIL': begin
                case _coord_system of
                    'SIMULATION': begin
                        xrange = -self.xrange
                        zrange =  self.zrange
                    endcase
                    'MAGNETOPAUSE': begin
                        xrange = reverse(self.zrange)
                        zrange = -self.xrange
                    endcase
                    'MAGNETOTAIL': ;Do nothing
                endcase
            endcase
            
            else: message, 'Coordinate system not recognized:"' + coord_system + '".'
        endcase
        
        ;Set object properties
        self.xrange       = xrange
        self.zrange       = zrange
        self.coord_system = coord_system
        
        ;Remake the simulation domain
        self -> MakeSimDomain
    endif
    
    ;MVA_FRAME
    if n_elements(mva_frame) gt 0 then begin
        self.mva_frame = keyword_set(mva_frame)
        
        ;Pick the axis labels
        case self.coord_system of
            'SIMULATION':   self.axis_labels = self.mva_frame ? ['x', 'y', 'z'] : ['x', 'y', 'z']
            'MAGNETOPAUSE': self.axis_labels = self.mva_frame ? ['L', 'M', 'N'] : ['z', 'y', 'x']
            'MAGNETOTAIL':  self.axis_labels = self.mva_frame ? ['N', 'M', 'L'] : ['x', 'y', 'z']
        endcase
    endif
    
    ;Clear all of the data
    ;   - Data is read for a particular time, [xyz]-range, orientation, etc.
    ;   - Data needs to be cleared if any of this changes (so that it is read with the
    ;       new parameters).
    self -> Clear_Data
end


;+
;   The purpose of this program is to make horizontal and vertical cuts through the data.
;   'XZ' and 'XY' cuts are determined by the value of the ORIENTATION property.
;
; :Params:
;       DATA_PRODUCT:       in, required, type=string
;                       The name of the data product to be read. For a list of
;                           available data product, call mr_readSIM without any
;                           arguments.
;       LOCATIONS:          in, required, type=lonarr
;                       The locations along the X-axis where vertical cuts are to
;                           be taken. Locations are povided in units of "di", the
;                           simulation lenght-scale. If `NCUTS` is provided, then
;                           this parameter outputs the locations of `NCUT` number
;                           of cuts equally spaced throughout `CUT_XRANGE`.
;       POS:                out, optional, type=fltarr
;                       The positions along the cut-axis corresponding to
;                               [XYZ]SIM[`IRANGE`[0]:`IRANGE`[1]].
; :Keywords:
;       HCUT_RANGE:     in, optional, type=fltarr(2)
;                       The horizontal range, in data coordinates, over which to cut.
;                           "Horizontal" is defined by the `HORIZONTAL` keyword. Data
;                           coordinates are determined by the "ion_scale" property.
;                           The default is to take the appropriate simulation range.
;       HORIZONTAL:     in, optional, type=boolean, default=0
;                       If set, a horizontal cut will be taken. The default is
;                           to take a vertical cut. For an "XY" orientation, "X" is
;                           horizontal, "Y" is vertical. Similar for "XZ", etc.
;       ICUTS:              out, optional, type=intarr
;                       The index values within [XYZ]Sim that define the cut-axis.
;       IRANGE:         out, optional, type=intarr
;                       The index values within [XYZ]Sim that span `HCUT_RANGE` or
;                           `VCUT_RANGE`.
;       NCUTS:          in, optional, type=int
;                       If provided, the [HV]CUT_RANGE perpendicutlar to the direction
;                           indicated by `HORIZONTAL` will be divided into NCUT evenly
;                           spaced locations. The `LOCATIONS` keyword will be ignored
;                           and these locations will be used for cuts.
;       VCUT_RANGE:     in, optional, type=fltarr(2)
;                       The verical range, in data coordinates, over which to cut.
;                           "Vertical" is defined by the `HORIZONTAL` keyword. Data
;                           coordinates are determined by the "ion_scale" property.
;                           The default is to take the appropriate simulation range.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the GetData method is also accepted
;                           via keyword inheritance.
;
; :Returns:
;       CUT:            The data along the vertical or horizontal cut. If the data
;                           requested does not exist, !Null will be returned.
;-
function MrSim::LineCuts, data_product, locations, pos, $
HCUT_RANGE=hcut_range, $
HORIZONTAL=horizontal, $
ICUTS=icuts, $
IRANGE=iRange, $
NCUTS=ncuts, $
VCUT_RANGE=vcut_range, $
_REF_EXTRA=extra
    compile_opt strictarr

    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, -1
    endif

    case self.orientation of
        'XY': cut = self -> XYCuts(data_product, locations, pos, NCUTS=nCuts, $
                                   CUT_XRANGE=hcut_range, CUT_YRANGE=vcut_range, $
                                   HORIZONTAL=horizontal, ICUTS=iCuts, IRANGE=iRange, $
                                   _EXTRA=extra)
                
        'XZ': cut = self -> XZCuts(data_product, locations, pos, NCUTS=nCuts, $
                                   CUT_XRANGE=hcut_range, CUT_ZRANGE=vcut_range, $
                                   HORIZONTAL=horizontal, ICUTS=iCuts, IRANGE=irange, $
                                   _EXTRA=extra)
    endcase
    
    return, cut
end


;+
;   The purpose of this program is to make vertical cuts through image data.
;
;   Note that if `CUT_[XZ]RANGE` fall outside of the current data range, they will
;   be set equal to the end points of the data range: [XZ]Range.
;
; :Private:
;
; :Params:
;       DATA_PRODUCT:       in, required, type=string
;                           The name of the data product to be read. For a list of
;                               available data product, call mr_readSIM without any
;                               arguments.
;       LOCATIONS:          in, required, type=lonarr
;                           The locations along the X-axis where vertical cuts are to
;                               be taken. Locations are povided in units of "di", the
;                               simulation lenght-scale. If `NCUTS` is provided, then
;                               this parameter outputs the locations of `NCUT` number
;                               of cuts equally spaced throughout `CUT_XRANGE`.
;       POS:                out, optional, type=fltarr
;                           The positions along the cut-axis corresponding to
;                               [XYZ]SIM[`IRANGE`[0]:`IRANGE`[1]].
; :Keywords:
;       CUT_XRANGE:         in, optional, type=fltarr(2), default=`XRANGE`
;                           The x-range over which to take vertical cuts. Use this in
;                               combination with `NCUTS`.
;       CUT_YRANGE:         in, optional, type=fltarr(2), default=`YRANGE`
;                           The z-range over which the vertical cut extents.
;       HORIZONTAL:         in, optional, type=boolean, default=0
;                           If set, a horizontal cut will be taken. The default is
;                               to take a vertical cut.
;       ICUTS:              out, optional, type=intarr
;                           The index values within XSim at which the vertical cuts
;                               were taken.
;       IRANGE:             out, optional, type=intarr
;                           The index values within [XYZ]Sim that span `CUT_XRANGE` or
;                               `CUT_YRANGE`.
;       NCUTS:              in, optional, type=int
;                           If provided, this specifies the number of equally spaced
;                               cuts to take between CUT_XRANGE[0] and CUT_XRANGE[1].
;                               When this keyword is used, `LOCATIONS` is ignored.
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by the GetData method is also accepted
;                               via keyword inheritance.
;
; :Returns:
;   CUT:                    The data along the vertical or horizontal cut. If the data
;                               product does not exist or an error occurs, !NULL will
;                               be returned.
;-
function MrSim::XYCuts, data_product, locations, pos, $
 CUT_XRANGE = cut_xrange, $
 CUT_YRANGE = cut_yrange, $
 HORIZONTAL = horizontal, $
 ICUTS = iCuts, $
 IRANGE = iRange, $
 NCUTS = ncuts, $
_REF_EXTRA = extra
    compile_opt strictarr

    ;Make sure the data suites our needs.
    if self.orientation ne 'XY' then $
        message, 'The orientation is "' + self.orientation + '". Cannot take XYCut.'

;-------------------------------------------------------
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------
    ;By using the GetProperty method, we limit ourselves to the
    ;range from which the data was taken.
    self -> GetProperty, XSIM=XSim, YSIM=YSim

    ;Set defaults
    horizontal = keyword_set(horizontal)
    if n_elements(cut_xrange) eq 0 then self -> GetProperty, XRANGE=cut_xrange
    if n_elements(cut_yrange) eq 0 then self -> GetProperty, YRANGE=cut_yrange
    if n_elements(ncuts) gt 0 then begin
        if horizontal $
            then locations = linspace(cut_yrange[0], cut_yrange[1], ncuts) $
            else locations = linspace(cut_xrange[0], cut_xrange[1], ncuts)
    endif

    ;How many cuts are being taken?
    ncuts = n_elements(locations)
    
;-------------------------------------------------------
;Read Data /////////////////////////////////////////////
;-------------------------------------------------------

    ;Get the data
    data = self -> GetData(data_product, _STRICT_EXTRA=extra)
    if data eq !Null then return, !Null

;-------------------------------------------------------
;Vertical Cut //////////////////////////////////////////
;-------------------------------------------------------
    if horizontal eq 0 then begin
        ;Get the index range over which the vertcal cuts span
        iRange = MrIndexRange(YSim, cut_yrange)
        pos = YSim[iRange[0]:iRange[1]]

        ;X-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(XSim, locations)
        void = MrIsMember(XSim[iCuts], locations, COMPLEMENT=bumpThese)
        if n_elements(bumpThese) ne 0 then iCuts[bumpThese] += 1
    
        ;data as a function of z along the vertical line.
        ;   Transpose to be consistent with vertical cuts: [pos, locations]
        cut = transpose(data[icuts, iRange[0]:iRange[1]])

;-------------------------------------------------------
;Horizontal Cut ////////////////////////////////////////
;-------------------------------------------------------
    endif else begin
        ;Get the index range over which the vertcal cuts span
        iRange = MrIndexRange(XSim, cut_xrange)
        pos = XSim[iRange[0]:iRange[1]]

        ;Z-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(YSim, locations)
        void = MrIsMember(YSim[iCuts], locations, COMPLEMENT=bumpThese)
        if n_elements(bumpThese) ne 0 then icut_pts[bumpThese] += 1

        ;data as a function of x along the horizontal line.
        cut = data[iRange[0]:iRange[1], iCuts]
    endelse
    
    return, cut
end


;+
;   The purpose of this program is to make vertical cuts through image data.
;
;   Note that if `CUT_[XZ]RANGE` fall outside of the current data range, they will
;   be set equal to the end points of the data range: [XZ]Range.
;
; :Private:
;
; :Params:
;       DATA_PRODUCT:   in, required, type=string
;                       The name of the data product to be read. For a list of
;                           available data product, call mr_readSIM without any
;                           arguments.
;       LOCATIONS:      in, required, type=lonarr
;                       The locations along the X-axis where vertical cuts are to
;                           be taken. Locations are povided in units of "di", the
;                           simulation lenght-scale. If `NCUTS` is provided, then
;                           this parameter outputs the locations of `NCUT` number
;                           of cuts equally spaced throughout `CUT_XRANGE`.
;       POS:            out, optional, type=fltarr
;                       The positions along the cut-axis corresponding to
;                                   [XYZ]SIM[`IRANGE`[0]:`IRANGE`[1]].
; :Keywords:
;       CUT_XRANGE:     in, optional, type=fltarr(2), default=`XRANGE`
;                       The x-range over which to take vertical cuts. Use this in
;                           combination with `NCUTS`.
;       CUT_ZRANGE:     in, optional, type=fltarr(2), default=`ZRANGE`
;                       The z-range over which the vertical cut extents.
;       HORIZONTAL:     in, optional, type=boolean, default=0
;                       If set, a horizontal cut will be taken. The default is
;                           to take a vertical cut.
;       ICUTS:          out, optional, type=intarr
;                       The index values within XSim at which the vertical cuts
;                           were taken.
;       IRANGE:         out, optional, type=intarr
;                       The index values within [XYZ]Sim that span `CUT_XRANGE` or
;                                   `CUT_YRANGE`.
;       NCUTS:          in, optional, type=int
;                       If provided, this specifies the number of equally spaced
;                           cuts to take between CUT_XRANGE[0] and CUT_XRANGE[1].
;                           When this keyword is used, `LOCATIONS` is ignored.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the GetData method is also accepted
;                           via keyword inheritance.
;
; :Returns:
;       CUT:            The data along the vertical or horizontal cut. If the data
;                           product does not exist or an error occurs, !NULL will
;                           be returned.
;-
function MrSim::XZCuts, data_product, locations, pos, $
 CUT_XRANGE = cut_xrange, $
 CUT_ZRANGE = cut_zrange, $
 HORIZONTAL = horizontal, $
 ICUTS = iCuts, $
 IRANGE = iRange, $
 NCUTS = ncuts, $
_REF_EXTRA = extra
    compile_opt strictarr

    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    ;Make sure the data suites our needs.
    if self.orientation ne 'XZ' then $
        message, 'The orientation is "' + self.orientation + '". Cannot take XZCut.'
    
;-------------------------------------------------------
;Define Ranges /////////////////////////////////////////
;-------------------------------------------------------
    ;By using the GetProperty method, we limit ourselves to the
    ;range from which the data was taken.
    self -> GetProperty, XSIM=XSim, ZSIM=ZSim

    ;Set defaults
    horizontal = keyword_set(horizontal)
    if n_elements(cut_xrange) eq 0 then self -> GetProperty, XRANGE=cut_xrange
    if n_elements(cut_zrange) eq 0 then self -> GetProperty, ZRANGE=cut_zrange
    if n_elements(ncuts) gt 0 then locations = linspace(cut_xrange[0], cut_xrange[1], ncuts)

    ;How many cuts are being taken?
    ncuts = n_elements(locations)
    
;-------------------------------------------------------
;Read Data /////////////////////////////////////////////
;-------------------------------------------------------

    ;Get the data
    data = self -> GetData(data_product, _STRICT_EXTRA=extra)
    if data eq !Null then return, !Null

;-------------------------------------------------------
;Horizontal Cut ////////////////////////////////////////
;-------------------------------------------------------
    if horizontal then begin
        ;Get the index range over which the vertcal cuts span
        iRange = MrIndexRange(XSim, cut_xrange, STRIDE=stride)
        pos = XSim[iRange[0]:iRange[1]:stride]

        ;Z-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(ZSim, locations)
        void = MrIsMember(ZSim[iCuts], locations, COMPLEMENT=bumpThese)
        if n_elements(bumpThese) ne 0 then icuts[bumpThese] += 1

        ;data as a function of x along the horizontal line.
        cut = data[iRange[0]:iRange[1]:stride, iCuts]

;-------------------------------------------------------
;Vertical Cut //////////////////////////////////////////
;-------------------------------------------------------
    endif else begin
        ;Get the index range over which the vertcal cuts span
        iRange = MrIndexRange(ZSim, cut_zrange, STRIDE=stride)
        pos = ZSim[iRange[0]:iRange[1]:stride]

        ;X-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(XSim, locations)
        void = MrIsMember(XSim[iCuts], locations, COMPLEMENT=bumpThese)
        if n_elements(bumpThese) ne 0 then icuts[bumpThese] += 1
    
        ;data as a function of z along the vertical line.
        ;   Transpose to be consistent with vertical cuts: [pos, locations]
        cut = transpose(data[icuts, iRange[0]:iRange[1]:stride])
    endelse
    
    return, cut
end


;+
;   The purpose of this program is to calculate the electron plasma beta::
;       \beta_{e} = \frac{Tr(P_{e})} {B^{2} / {2 \mu_{0}}
;
; :Private:
;
; :Returns:
;       BETA_E:                 The electron beta.
;-
function MrSim::Beta_e
    compile_opt strictarr, hidden
    on_error, 2
	
	;Read the electric and magnetic field data.
    Pe_xx = self -> getData('Pe-xx')
    Pe_yy = self -> getData('Pe-yy')
    Pe_zz = self -> getData('Pe-zz')
    Bmag = self -> getData('Bmag')
	
	;Calculate the dot product between E and B. Divide by 3 for average.
    Beta_e = (Pe_xx + Pe_yy + Pe_zz)  / (6.0 * Bmag^2)
    
    return, Beta_e
end


;+
;   The purpose of this program is to calculate the ion plasma beta::
;       \beta_{i} = \frac{Tr(P_{i})} {B^{2} / {2 \mu_{0}}
;
; :Private:
;
; :Returns:
;       BETA_I:                 The ion beta.
;-
function MrSim::Beta_i
    compile_opt strictarr, hidden
    on_error, 2
	
	;Read the electric and magnetic field data.
    Pi_xx = self -> getData('Pi-xx')
    Pi_yy = self -> getData('Pi-yy')
    Pi_zz = self -> getData('Pi-zz')
    Bmag = self -> getData('Bmag')
	
	;Calculate the dot product between E and B. Divide by 3 for average.
    Beta_i = (Pi_xx + Pi_yy + Pi_zz)  / (6.0 * Bmag^2)
    
    return, Beta_i
end


;+
;   The purpose of this program is to calculate the plasma beta::
;       \beta_{i} = \frac{Tr(P_{i})} {B^{2} / {2 \mu_{0}}
;
; :Private:
;
; :Returns:
;       BETA_P:                 The overall plasma beta.
;-
function MrSim::Beta_p
    compile_opt strictarr, hidden
    on_error, 2
	
	;Read the electric and magnetic field data.
    Pe_xx = self -> getData('Pe-xx')
    Pe_yy = self -> getData('Pe-yy')
    Pe_zz = self -> getData('Pe-zz')
    Pi_xx = self -> getData('Pi-xx')
    Pi_yy = self -> getData('Pi-yy')
    Pi_zz = self -> getData('Pi-zz')
    Bmag = self -> getData('Bmag')
	
	;Calculate the dot product between E and B. Divide by 6 for average.
    Beta_p = (Pe_xx + Pe_yy + Pe_zz + Pi_xx + Pi_yy + Pi_zz)  / (12.0 * Bmag^2)
    
    return, Beta_p
end


;+
;   The purpose of this program is to calculate the magnitude of the Magnetic Field::
;       \left| B \right| = \sqrt{ \vec{B} \cdot \vec{B} }
;                        = \sqrt{ B_{x}^{2} + B_{y}^{2} + B_{z}^{2} }
;
; :Private:
;
; :Returns:
;       B_MAG:                  The magnitude of the magnetic field.
;-
function MrSim::Bmag
	compile_opt strictarr, hidden
	on_error, 2
	
	;Read the electric field data.
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of E
    B_mag = sqrt(Bx^2 + By^2 + Bz^2)
    
    return, B_mag
end


;+
;   The purpose of this program is to calculate magnitude of the divergence of the
;   electron pressure tensor.
;       \left| \nabla \cdot P_{e} \right| = \sqrt{ \left[ \nabla \cdot (P_{e})_{x} \right]^{2} + 
;                                                  \left[ \nabla \cdot (P_{e})_{y} \right]^{2} + 
;                                                  \left[ \nabla \cdot (P_{e})_{z} \right]^{2} }
;
; :Private:
;
; :Returns:
;       divPe_mag:                The magnitude of the divergence of the electron pressure
;                                   tensor.
;-
function MrSim::divPe_mag
	compile_opt strictarr, hidden

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the x-components of the Electric Field and the Current Density
    divPe_x = self -> getData('divPe_x')
    divPe_z = self -> getData('divPe_z')
    divPe_y = 0     ;d/dy = 0
    
    ;Sum the terms, dividing by electron density
    divPe_mag = sqrt(divPe_x^2 + divPe_z^2)
    
    return, divPe_x
end


;+
;   The purpose of this program is to calculate x-component of the divergence of the
;   electron pressure tensor.
;       (\nabla \cdot P_{e})_{x} = \frac{1} {n_{e}}
;                                  \left( \frac{\partial (P_{e})_{xx}} {\partial x} +
;                                         \frac{\partial (P_{e})_{yx}} {\partial x} +
;                                         \frac{\partial (P_{e})_{zx}} {\partial x} \right)
;
; :Private:
;
; :Returns:
;       divPe_x:                The x-component of the divergence of the electron pressure
;                                   tensor.
;-
function MrSim::divPe_x
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
    Pe_xx = self -> getData('Pe-xx')
    Pe_xy = self -> getData('Pe-xy')
    Pe_xz = self -> getData('Pe-xz')
    n_e   = self -> getData('ne')
    dims = size(Pe_xx, /DIMENSIONS)
    
    ;Get the x-size of a grid cell in electron skin depth
    self -> GetInfo, DX_DE=dx_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    dPe_xx = fltarr(dims)
    dPe_xy = fltarr(dims)
    dPe_xz = fltarr(dims)
    dPe_x  = fltarr(dims)
    
    ;Compute the central difference
    dPe_xx[1:dims[0]-2,*] = (Pe_xx[2:dims[0]-1,*] - Pe_xx[0:dims[0]-3,*]) / (2.0 * dx_de)
    dPe_xy[1:dims[0]-2,*] = (Pe_xy[2:dims[0]-1,*] - Pe_xy[0:dims[0]-3,*]) / (2.0 * dx_de)
    dPe_xz[1:dims[0]-2,*] = (Pe_xz[2:dims[0]-1,*] - Pe_xz[0:dims[0]-3,*]) / (2.0 * dx_de)
    
    ;Sum the terms, dividing by electron density
    divPe_x = (dPe_xx + dPe_xy + dPe_xz) / n_e
    
    return, divPe_x
end


;+
;   The purpose of this program is to calculate z-component of the divergence of the
;   electron pressure tensor.
;       (\nabla \cdot P_{e})_{z} = \frac{1} {n_{e}}
;                                  \left( \frac{\partial (P_{e})_{xz}} {\partial z} +
;                                         \frac{\partial (P_{e})_{yz}} {\partial z} +
;                                         \frac{\partial (P_{e})_{zz}} {\partial z} \right)
;
; :Private:
;
; :Params:
;       TIME:                   in, required, type=long
;                               The time-step at which the data is desired.
;
; :Keywords:
;       XRANGE:                 in, optional, type=fltarr(2), default=[0,xmax]
;                               The xrange over which data is to be returned.
;       ZRANGE:                 in, optional, type=fltarr(2), default=[-zmax/2\, zmax/2]
;                               The zrange over which data is to be returned.
;
; :Returns:
;       divPe_z:                The z-component of the divergence of the electron pressure
;                                   tensor.
;-
function MrSim::divPe_z
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor.
	;Pe-zx and Pe-zy are unavailable -- assume the tensor is symmetric.
    Pe_zx = self -> getData('Pe-xz')
    Pe_zy = self -> getData('Pe-yz')
    Pe_zz = self -> getData('Pe-zz')
    n_e = self -> getData('ne')
    
    dims = size(Pe_zz, /DIMENSIONS)
    
    ;Cell size in de
    self -> GetInfo, DZ_DE=dz_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    dPe_zx  = fltarr(dims)
    dPe_zy  = fltarr(dims)
    dPe_zz  = fltarr(dims)
    divPe_z = fltarr(dims)

    ;Compute the central difference
    dPe_zx[*,1:dims[1]-2] = (Pe_zx[*,2:dims[1]-1] - Pe_zx[*,0:dims[1]-3]) / (2.0 * dz_de)
    dPe_zy[*,1:dims[1]-2] = (Pe_zy[*,2:dims[1]-1] - Pe_zy[*,0:dims[1]-3]) / (2.0 * dz_de)
    dPe_zz[*,1:dims[1]-2] = (Pe_zz[*,2:dims[1]-1] - Pe_zz[*,0:dims[1]-3]) / (2.0 * dz_de)
    
    ;Sum the terms
    divPe_z = (dPe_zx + dPe_zy + dPe_zz) / n_e
    
    return, divPe_z
end


;+
;   The purpose of this program is to calculate x-component of the divergence of the
;   electron pressure tensor.
;       (\nabla \cdot P_{i})_{x} = \frac{1} {n_{i}}
;                                  \left( \frac{\partial (P_{i})_{xx}} {\partial x} +
;                                         \frac{\partial (P_{i})_{yx}} {\partial x} +
;                                         \frac{\partial (P_{i})_{zx}} {\partial x} \right)
;
; :Private:
;
; :Returns:
;       divPe_x:                The x-component of the divergence of the electron pressure
;                                   tensor.
;-
function MrSim::divPi_x
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
    Pi_xx = self -> getData('Pi-xx')
    Pi_xy = self -> getData('Pi-xy')
    Pi_xz = self -> getData('Pi-xz')
    n_i = self -> getData('ni')
    
    dims = size(Pi_xx, /DIMENSIONS)
    
    ;Get the x-size of a grid cell in electron skin depth
    self -> GetInfo, DX_DE=dx_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    dPi_xx = fltarr(dims)
    dPi_xy = fltarr(dims)
    dPi_xz = fltarr(dims)
    dPi_x  = fltarr(dims)
    
    ;Compute the central difference
    dPi_xx[1:dims[0]-2,*] = (Pi_xx[2:dims[0]-1,*] - Pi_xx[0:dims[0]-3,*]) / (2.0 * dx_de)
    dPi_xy[1:dims[0]-2,*] = (Pi_xy[2:dims[0]-1,*] - Pi_xy[0:dims[0]-3,*]) / (2.0 * dx_de)
    dPi_xz[1:dims[0]-2,*] = (Pi_xz[2:dims[0]-1,*] - Pi_xz[0:dims[0]-3,*]) / (2.0 * dx_de)
    
    ;Sum the terms, dividing by electron density
    divPi_x = (dPi_xx + dPi_yx + dPi_zx) / n_i
    
    return, divPi_x
end


;+
;   The purpose of this program is to calculate z-component of the divergence of the
;   electron pressure tensor.
;       (\nabla \cdot P_{i})_{z} = \frac{1} {n_{i}}
;                                  \left( \frac{\partial (P_{i})_{xz}} {\partial z} +
;                                         \frac{\partial (P_{i})_{yz}} {\partial z} +
;                                         \frac{\partial (P_{i})_{zz}} {\partial z} \right)
;
; :Private:
;
; :Params:
;       TIME:                   in, required, type=long
;                               The time-step at which the data is desired.
;
; :Keywords:
;       XRANGE:                 in, optional, type=fltarr(2), default=[0,xmax]
;                               The xrange over which data is to be returned.
;       ZRANGE:                 in, optional, type=fltarr(2), default=[-zmax/2\, zmax/2]
;                               The zrange over which data is to be returned.
;
; :Returns:
;       divPi_z:                The z-component of the divergence of the electron pressure
;                                   tensor.
;-
function MrSim::divPi_z
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor.
	;Pe-zx and Pe-zy are unavailable -- assume the tensor is symmetric.
    Pi_zx = self -> getData('Pi-xz')
    Pi_zy = self -> getData('Pi-yz')
    Pi_zz = self -> getData('Pi-zz')
    n_i   = self -> getData('ni')
    
    dims = size(Pi_zz, /DIMENSIONS)
    
    ;Cell size in de
    self -> GetInfo, DZ_DE=dz_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    dPi_zx  = fltarr(dims)
    dPi_zy  = fltarr(dims)
    dPi_zz  = fltarr(dims)
    divPi_z = fltarr(dims)

    ;Compute the central difference
    dPi_zx[*,1:dims[1]-2] = (Pi_zx[*,2:dims[1]-1] - Pi_zx[*,0:dims[1]-3]) / (2.0 * dz_de)
    dPi_zy[*,1:dims[1]-2] = (Pi_zy[*,2:dims[1]-1] - Pi_zy[*,0:dims[1]-3]) / (2.0 * dz_de)
    dPi_zz[*,1:dims[1]-2] = (Pi_zz[*,2:dims[1]-1] - Pi_zz[*,0:dims[1]-3]) / (2.0 * dz_de)
    
    ;Sum the terms
    divPi_z = (dPi_zx + dPi_zy + dPi_zz) / n_i
    
    return, divPi_z
end


;+
;   The purpose of this program is to calculate the dot product of the Electric and
;   Magnetic Fields::
;       E_par = \vec{E} \cdot \vec{B} = Ex Bx + Ey By + Ez Bz / |B|
;
; :Private:
;
; :Returns:
;       E_PARA:                 The electric field strength in the parallel-to-B
;                                   direction.
;-
function MrSim::E_par
    compile_opt strictarr, hidden
    on_error, 2
	
	;Read the electric and magnetic field data.
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
	
	;Calculate the dot product between E and B.
    E_par = Ex*Bx + Ey*By + Ez*Bz / sqrt(Bx^2 + By^2 + Bz^2)
    
    return, E_par
end


;+
;   The purpose of this program is to calculate the dot product between the electric
;   field and the current density::
;       \vec{E} \cdot \vec{J} = E_{x}*J_{x} + E_{y}*J_{y} + E_{z}*J_{z}
;
; :Private:
;
; :Returns:
;       E_dot_J:                The electric field dotted with the current density.
;-
function MrSim::E_dot_J
    compile_opt strictarr, hidden
    on_error, 2

    ;Get the electric field and current density data
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Jx = self -> getData('Jx')
    Jy = self -> getData('Jy')
    Jz = self -> getData('Jz')

    ;Calculate the dot product
    E_dot_J = Ex*Jx + Ey*Jy + Ez*Jz
    
    return, E_dot_J
end


;+
;   The purpose of this program is to calculate the energy conversion term parallel to
;   the magnetic field::
;       E_{\parallel} J_{\parallel} = ( \vec{E} \cdot \vec{B} ) ( \vec{J} \cdot \vec{B} )
;
; :Private:
;
; :Returns:
;       E_PAR:                  The electric field strength in the parallel-to-B
;                                   direction.
;-
function MrSim::EJpar
	compile_opt strictarr, hidden
	on_error, 2
	
	;Read the electric and magnetic field data.
    Epar = self -> getData('Epar')
    Jpar = self -> getData('Jpar')
    
    ;Dot the electric field with B_had
    EJpar = Epar * Jpar
    
    return, EJpar
end


;+
;   The purpose of this program is to calculate the magnitude of the Electric Field::
;       \left| E \right| = \sqrt{ \vec{E} \cdot \vec{E} }
;                        = \sqrt{ E_{x}^{2} + E_{y}^{2} + E_{z}^{2} }
;
; :Private:
;
; :Returns:
;       E_MAG:                  The magnitude of the electric field.
;-
function MrSim::Emag
	compile_opt strictarr, hidden
	
	;Read the electric field data.
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    
    ;Calculate the magnitude of E
    E_mag = sqrt(Ex^2 + Ey^2 + Ez^2)
    
    return, E_mag
end


;+
;   The purpose of this program is to calculate the Electric Field parallel to the
;   magnetic field::
;       B_hat = [Bx, By, Bz] / |B|
;       E_par = E dot B_hat = Ex*Bx_hat + Ey*By_hat + Ez*Bz_hat
;
; :Private:
;
; :Returns:
;       E_PAR:                  The electric field strength in the parallel-to-B
;                                   direction.
;-
function MrSim::Epar
	compile_opt strictarr, hidden
	
	;Read the electric and magnetic field data.
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate a unit vector pointing in the direction of the magnetic field
    Bx_unit = Bx/SQRT(Bx^2+By^2+Bz^2)
    By_unit = By/SQRT(Bx^2+By^2+Bz^2)
    Bz_unit = Bz/SQRT(Bx^2+By^2+Bz^2)
    
    ;Dot the electric field with B_had
    E_par = Ex*Bx_unit + Ey*By_unit + Ez*Bz_unit
    
    return, E_par
end


;+
;   The purpose of this program is to calculate the Electric Field perpendicular to the
;   magnetic field::
;       E_{\bot} = \left| E \right| - E_{\parallel}
;
; :Private:
;
; :Returns:
;       E_PERP:                 The electric field strength in the perpendicular-to-B
;                                   direction.
;-
function MrSim::Eperp
	compile_opt strictarr, hidden
	
	;Read the electric and magnetic field data.
    E_mag = self -> getData('Emag')
    E_par = self -> getData('Epar')
    
    ;Dot the electric field with B_had
    E_perp = E_mag - E_par
    
    return, E_perp
end




;+
;   The purpose of this program is to calculate the z-component of the cross product 
;   between the Electric and Magnetic Fields::
;       |B| = sqrt( Bx^2 + By^2 + Bz^2 )
;       ExB_x = (Ey*Bz - Ez*By) / |B|^2
;       ExB_y = (Ez*Bx - Ex*Bz) / |B|^2
;       ExB_z = (Ex*By - Ey*Bx) / |B|^2
;       ExB = ExB_x + ExB_y + ExB_z
;       |ExB| = sqrt(ExB dot ExB) = sqrt( ExB_x^2 + ExB_y^2 + ExB_z^2 )
;
; :Private:
;
; :Returns:
;       ExB_z:                  The z-component of the cross product between the Electric
;                                   and Magnetic Fields.
;-
function MrSim::v_ExB, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	
	;Determine what type of plot should be made
	magnitude = keyword_set(magnitude)
	x         = keyword_set(x)
	y         = keyword_set(y)
	z         = keyword_set(z)
	if (x + y + z + magnitude) gt 1 then message, 'MAGNITUDE, X, Y, and Z are mutually exclusive.'
	
	;Get the Electric and Magnetic Field data
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Compute the magnitude of B
    B_mag_squared = Bx^2 + By^2 + Bz^2
    
    ;Calculate the magnitude of ExB
    if x || magnitude then ExB_x   = (Ey*Bz - Ez*By) / B_mag_squared
    if y || magnitude then ExB_y   = (Ez*Bx - Ex*Bz) / B_mag_squared
    if z || magnitude then ExB_z   = (Ex*By - Ey*Bx) / B_mag_squared
    if magnitude      then ExB_mag = sqrt(ExB_x^2 + ExB_y^2 + ExB_z^2)
    
    ;Return the correct quantity
    case 1 of
        magnitude: return, ExB_mag
        x:         return, ExB_x
        y:         return, ExB_y
        z:         return, ExB_z
        else:      ;Not possible
    endcase
end


;+
;   The purpose of this program is to calculate the z-component of the cross product 
;   between the Electric and Magnetic Fields::
;       |B| = sqrt( Bx^2 + By^2 + Bz^2 )
;       ExB_x = (Ey*Bz - Ez*By) / |B|^2
;       ExB_y = (Ez*Bx - Ex*Bz) / |B|^2
;       ExB_z = (Ex*By - Ey*Bx) / |B|^2
;       ExB = ExB_x + ExB_y + ExB_z
;       |ExB| = sqrt(ExB dot ExB) = sqrt( ExB_x^2 + ExB_y^2 + ExB_z^2 )
;
; :Private:
;
; :Returns:
;       ExB_z:                  The z-component of the cross product between the Electric
;                                   and Magnetic Fields.
;-
function MrSim::vExB_mag
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    B_mag_squared = Bx^2 + By^2 + Bz^2
    ExB_x  = (Ey*Bz - Ez*By) / B_mag_squared
    ExB_y  = (Ez*Bx - Ex*Bz) / B_mag_squared
    ExB_z  = (Ex*By - Ey*Bx) / B_mag_squared
    ExB_mag = sqrt(ExB_x^2 + ExB_y^2 + ExB_z^2)
    
    return, ExB_mag
end


;+
;   The purpose of this program is to calculate the x-component of the cross product 
;   between the Electric and Magnetic Fields::
;       |B| = sqrt( Bx^2 + By^2 + Bz^2 )
;       ExB_x = (Ey*Bz - Ez*By) / |B|^2
;       ExB_y = (Ez*Bx - Ex*Bz) / |B|^2
;       ExB_z = (Ex*By - Ey*Bx) / |B|^2
;       ExB = ExB_x + ExB_y + ExB_z
;
; :Private:
;
; :Returns:
;       ExB_x:                  The x-component of the cross product between the Electric
;                                   and Magnetic Fields.
;-
function MrSim::vExB_x
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the x-component of ExB
    B_mag_squared = Bx^2 + By^2 + Bz^2
    ExB_x  = (Ey*Bz - Ez*By) / B_mag_squared
    
    return, ExB_x
end


;+
;   The purpose of this program is to calculate the x-component of the cross product 
;   between the Electric and Magnetic Fields::
;       |B| = sqrt( Bx^2 + By^2 + Bz^2 )
;       ExB_x = (Ey*Bz - Ez*By) / |B|^2
;       ExB_y = (Ez*Bx - Ex*Bz) / |B|^2
;       ExB_z = (Ex*By - Ey*Bx) / |B|^2
;       ExB = ExB_x + ExB_y + ExB_z
;
; :Private:
;
; :Returns:
;       ExB_y:                  The y-component of the cross product between the Electric
;                                   and Magnetic Fields.
;-
function MrSim::vExB_y
	compile_opt strictarr, hidden
	on_error, 2
	
	;Get the Electric and Magnetic Field data
    Ex = self -> getData('Ex')
    Ez = self -> getData('Ez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the y-component of ExB
    B_mag_squared = Bx^2 + By^2 + Bz^2
    ExB_y  = (Ez*Bx - Ex*Bz) / B_mag_squared
    
    return, ExB_y
end


;+
;   The purpose of this program is to calculate the z-component of the cross product 
;   between the Electric and Magnetic Fields::
;       |B| = sqrt( Bx^2 + By^2 + Bz^2 )
;       ExB_x = (Ey*Bz - Ez*By) / |B|^2
;       ExB_y = (Ez*Bx - Ex*Bz) / |B|^2
;       ExB_z = (Ex*By - Ey*Bx) / |B|^2
;       ExB = ExB_x + ExB_y + ExB_z
;
; :Private:
;
; :Returns:
;       ExB_z:                  The z-component of the cross product between the Electric
;                                   and Magnetic Fields.
;-
function MrSim::vExB_z
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the y-component of ExB
    B_mag_squared = Bx^2 + By^2 + Bz^2
    ExB_z  = (Ex*By - Ey*Bx) / B_mag_squared
    
    return, ExB_z
end


;+
;   The purpose of this program is to calculate product of the x-components of the
;   Electric Field and the Current Density.
;       ExJx = Ex*Bx
;
; :Private:
;
; :Returns:
;       ExJx:                   The product of the x-components of the Electric Field
;                                   and the Current Density
;-
function MrSim::ExJx
	compile_opt strictarr, hidden
	
	;Get the x-components of the Electric Field and the Current Density
    Ex = self -> getData('Ex')
    Jx = self -> getData('Jx')
    
    ;Calculate the product
    ExJx = Ex * Jx
    
    return, ExJx
end


;+
;   The purpose of this program is to calculate product of the y-components of the
;   Electric Field and the Current Density::
;       EyJy = Ey*Jy
;
; :Private:
;
; :Returns:
;       EyJy:                   The product of the y-components of the Electric Field
;                                   and the Current Density
;-
function MrSim::EyJy
	compile_opt strictarr, hidden
	
	;Get the y-component of the Electric Field and the Current Density
    Ey = self -> getData('Ey')
    Jy = self -> getData('Jy')

    ;Calculate the product
    EyJy = Ey * Jy
    
    return, EyJy
end


;+
;   The purpose of this program is to calculate product of the y-components of the
;   Electric and Magnetic Fields::
;       EyBy = Ey*By
;
; :Private:
;
; :Returns:
;       EyBy:                   The product of the y-components of the Electric and
;                                   Magnetic Fields.
;-
function MrSim::EyBy
	compile_opt strictarr, hidden
	
	;Get the y-components of the Electric and Magnetic Fields
    Ey = self -> getData('Ey')
    By = self -> getData('By')
    
    ;Calculate the product
    EyBy = Ey * By
    
    return, EyBy
end


;+
;   The purpose of this program is to calculate product of the z-components of the
;   Electric Field and the Current Density::
;       EzJz = Ez * Jz
;
; :Private:
;
; :Returns:
;       EzJz:                   The product of the z-components of the Electric Field
;                                   and the Current Density
;-
function MrSim::EzJz
	compile_opt strictarr, hidden
	
	;Get the z-component of the Electric Field and the Current Density
    Ez = self -> getData('Ez')
    Jz = self -> getData('Jz')

    ;Calculate the product
    EzJz = Ez * Jz

    return, EzJz
end


;+
;   The purpose of this program is to calculate the magnitude of the lorentz force acting
;   on electrons::
;       F = E + (Ue x B))
;       Fx = Ex + (Uey*Bz - Uez*By)
;       Fy = Ey + (Uez*Bx - Uex*Bz)
;       Fz = Ez + (Uex*By - Uey*Bx)
;       |F| = sqrt(Fx^2 + Fy^2 + Fz^2)
;
; :Private:
;
; :Returns:
;       Fe_mag:                  The magnitude of the lorentz force.
;-
function MrSim::Fe_mag
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Fields and the Electron Velocty
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of the Lorentz Force
    Fex = Ex + Uey*Bz - Uez*By
    Fey = Ey + Uez*Bx - Uex*Bz
    Fez = Ez + Uex*By - Uey*Bx
    Fe_mag = sqrt(Fex^2 + Fey^2 + Fez^2)
    
    return, Fe_mag
end


;+
;   The purpose of this program is to calculate the x-component of the lorentz force
;   acting on electrons::
;       F = E + (Ue x B))
;       Fx = Ex + (Uey*Bz - Uez*By)
;
; :Private:
;
; :Returns:
;       Fx:                     The x-components of the lorentz force.
;-
function MrSim::Fex
	compile_opt strictarr, hidden
	
	;Get the data
    Ex = self -> getData('Ex')
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the Force
    Fex = Ex + Uey*Bz - Uez*By
    
    return, Fex
end


;+
;   The purpose of this program is to calculate the y-component of the lorentz force
;   acting on electrons::
;       F = E + (Ue x B))
;       Fy = Ey + (Uez*Bx - Uex*Bz)
;
; :Private:
;
; :Returns:
;       Fy:                     The y-components of the lorentz force.
;-
function MrSim::Fey
	compile_opt strictarr, hidden

    ;Get the relevant Electric and Magnetic Field, and Electron Velocity Data
    Ey = self -> getData('Ey')
    Uex = self -> getData('Uex')
    Uez = self -> getData('Uez')
    Bx = self -> getData('Bx')
    Bz = self -> getData('Bz')
    
    ;Calculate the y-component of the Lorentz Force
    Fey = Ey + Uez*Bx - Uex*Bz
    
    return, Fey
end


;+
;   The purpose of this program is to calculate the z-component of the lorentz force
;   acting on electrons::
;       F = E + (Ue x B))
;       Fz = Ez + (Uex*By - Uey*Bx)
;
; :Private:
;
; :Returns:
;       Fz:                     The z-components of the lorentz force.
;-
function MrSim::Fez
	compile_opt strictarr, hidden

    ;Get the relevant Electric and Magnetic Field, and Electron Velocity Data
    Ez = self -> getData('Ez')
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Bx = self -> getData('Bx')
    By = self -> getData('By')

    ;Calculate the y-component of the Lorentz Force
    Fez = Ez + Uex*By - Uey*Bx
    
    return, Fez
end



;+
;   The purpose of this program is to calculate the magnitude of the lorentz force acting
;   on ions::
;       F = E + (Ui x B))
;       Fx = Ex + (Uiy*Bz - Uiz*By)
;       Fy = Ey + (Uiz*Bx - Uix*Bz)
;       Fz = Ez + (Uix*By - Uiy*Bx)
;       |F| = sqrt(Fx^2 + Fy^2 + Fz^2)
;
; :Private:
;
; :Returns:
;       F_mag:                  The magnitude of the lorentz force acting on ions.
;-
function MrSim::Fi_mag
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Fields and the Ion Velocty
    Ex = self -> getData('Ex')
    Ey = self -> getData('Ey')
    Ez = self -> getData('Ez')
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')

    ;Calculate the magnitude of the Lorentz Force acting on ions
    Fx = Ex + Uiy*Bz - Uiz*By
    Fiy = Ey + Uiz*Bx - Uix*Bz
    Fiz = Ez + Uix*By - Uiy*Bx
    Fi_mag = sqrt(Fx^2 + Fiy^2 + Fiz^2)
    
    return, Fi_mag
end


;+
;   The purpose of this program is to calculate the x-component of the lorentz force
;   acting on electrons::
;       F = E + (Ui x B))
;       Fx = Ex + (Uiy*Bz - Uiz*By)
;
; :Private:
;
; :Returns:
;       Fix:                    The x-components of the lorentz force.
;-
function MrSim::Fix
	compile_opt strictarr, hidden
	on_error, 2
	
    ;Get the relevant Electric and Magnetic Fields and the Ion Velocty components.
    Ex = self -> getData('Ex')
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    By = self -> getData('By')
    Bz = self -> getData('Bz')

    ;Calculate the y-component of the Lorentz Force acting on ions
    Fix = Ex + (Uiy*Bz - Uiz*By)
    
    return, Fix
end


;+
;   The purpose of this program is to calculate the y-component of the lorentz force
;   acting on electrons::
;       F = E + (Ue x B))
;       Fy = Ey + (Uiz*Bx - Uix*Bz)
;
; :Private:
;
; :Returns:
;       Fiy:                    The y-components of the lorentz force.
;-
function MrSim::Fiy
	compile_opt strictarr, hidden
	
    ;Get the relevant Electric and Magnetic Fields and the Ion Velocty components.
    Ey = self -> getData('Ey')
    Uix = self -> getData('Uix')
    Uiz = self -> getData('Uiz')
    Bx = self -> getData('Bx')
    Bz = self -> getData('Bz')

    ;Calculate the y-component of the Lorentz Force acting on ions
    Fiy = Ey + (Uiz*Bx - Uix*Bz)
    
    return, Fiy
end


;+
;   The purpose of this program is to calculate the z-component of the lorentz force
;   acting on electrons::
;       F = E + (Ue x B))
;       Fz = Ez + (Uix*By - Uiy*Bx)
;
; :Private:
;
; :Returns:
;       Fiz:                    The z-components of the lorentz force.
;-
function MrSim::Fiz
	compile_opt strictarr, hidden
	
    ;Get the relevant Electric and Magnetic Fields and the Ion Velocty components.
    Ez = self -> getData('Ez')
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Bx = self -> getData('Bx')
    By = self -> getData('By')

    ;Calculate the y-component of the Lorentz Force acting on ions
    Fiz = Ez + (Uix*By - Uiy*Bx)
    
    return, Fiz
end


;+
;   The purpose of this program is to calculate the gradient of the scalar pressure.
;       \nabla P_{e} = \frac{1} {n_{e}}
;                                  \left( \frac{\partial P_{e}} {\partial x} +
;                                         \frac{\partial P_{e}} {\partial x} +
;                                         \frac{\partial P_{e}} {\partial x} \right)
;
; :Private:
;
; :Returns:
;       gradPe:         The x-component of the divergence of the electron pressure
;                           tensor.
;-
function MrSim::gradP
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
	;Pe-yx, Pe-zx are unavailable, so I am assuming the matrix is symmetric
    Pe  = self -> Read('Pe')
    n_e = self -> getData('ne')
    
    dims = size(Pe_xx, /DIMENSIONS)
    
    ;Make X (1D) the same size as Pe-xx (2D)
    self -> GetProperty, XSIM=xsim
    X = rebin(*self.XSim, dims[0], dims[1])

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    delta_Pe_xx = fltarr(dims)
    delta_Pe_yx = fltarr(dims)
    delta_Pe_zx = fltarr(dims)
    divPe_x = fltarr(dims)
    
    ;Compute the derivative
    delta_Pe_xx[1:*,*] = (Pe_xx[1:dims[0]-1,*] - Pe_xx[0:dims[0]-2,*]) / (X[1:dims[0]-1,*] - X[0:dims[0]-2,*])
    delta_Pe_yx[1:*,*] = (Pe_yx[1:dims[0]-1,*] - Pe_yx[0:dims[0]-2,*]) / (X[1:dims[0]-1,*] - X[0:dims[0]-2,*])
    delta_Pe_zx[1:*,*] = (Pe_zx[1:dims[0]-1,*] - Pe_zx[0:dims[0]-2,*]) / (X[1:dims[0]-1,*] - X[0:dims[0]-2,*])
    
    ;Sum the terms, dividing by electron density
    divPe_x = (delta_Pe_xx + delta_Pe_yx + delta_Pe_zx) / n_e
    
    return, divPe_x
end


;+
;   The purpose of this program is to calculate the z-component of the Electron
;   Current Density::
;       (J_{e})_{xz} = -1.0 * n_{e} * \sqrt{ U_{ex}^{2} + U_{ez}^{2} }
;
; :Private:
;
; :Returns:
;       Je_xz:          The in-plane Electron Current Density.
;-
function MrSim::Je_xz
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density  
    Uex = self -> getData('Uex')
    Uez = self -> getData('Uez')
    n_e = self -> getData('ne')
    
    ;Calculate the in-plane Current Density
    Ue_xz = sqrt(Uex^2 + Uez^2)
    Je_xz = -1.0 * n_e * Ue_xz
    
    return, Je_xz
end


;+
;   The purpose of this program is to calculate the x-component of the Electron
;   Current Density::
;       Jex = -1.0 * ne * Uex
;
; :Private:
;
; :Returns:
;       Jex:                    The x-component of the Electron Current Density.
;-
function MrSim::Jex
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density 
    Uex = self -> getData('Uex')
    n_e = self -> getData('ne')
    
    ;Calculate the Current Density
    Jex = -1.0 * Uex * n_e
    
    return, Jex
end


;+
;   The purpose of this program is to calculate the y-component of the Electron
;   Current Density::
;       Jey = -1.0 * ne * Uey
;
; :Private:
;
; :Returns:
;       Jey:                    The y-component of the Electron Current Density.
;-
function MrSim::Jey
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density 
    Uey = self -> getData('Uey')
    n_e = self -> getData('ne')
    if Uey eq !Null || n_e eq !Null then return, !Null
    
    ;Calculate the Current Density
    Jey = -1.0 * Uey * n_e
    
    return, Jey
end


;+
;   The purpose of this program is to calculate the z-component of the Electron
;   Current Density::
;       Jez = -1.0 * ne * Uez
;
; :Private:
;
; :Returns:
;       Jez:            The z-component of the Electron Current Density.
;-
function MrSim::Jez
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density 
    Uez = self -> getData('Uez')
    n_e = self -> getData('ne')
    
    ;Calculate the Current Density
    Jez = -1.0 * Uez * n_e
    
    return, Jez
end


;+
;   The purpose of this program is to calculate the x-component of the JxB portion
;   of the Generalized Ohm's Law::
;       (J_{e} \times B)_{x} = \frac{1} {n_{e}} (J_{e,y} B_{z} - J_{e,z} B_{y})
;
; :Private:
;
; :Returns:
;       JexB_x:             The x-component of the JxB force.
;-
function MrSim::JexB_x
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
    n_e = self -> getData('ne')
    Jey = self -> getData('Jey')
    Jez = self -> getData('Jez')
    By  = self -> getData('By')
    Bz  = self -> getData('Bz')
    
    ;Calculate the Current Density
    JexB_x = (Jey*Bz - Jez*By) / n_e
    
    return, JexB_x
end


;+
;   The purpose of this program is to calculate the y-component of the JxB force::
;       (J_{e} \times B)_{y} = \frac{1} {e n_{e}} (J_{e,x} B_{z} - J_{e,z} B_{x})
;
; :Private:
;
; :Returns:
;       JexB_y:             The x-component of the JxB force.
;-
function MrSim::JexB_y
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
    n_e = self -> getData('ne')
    Jex = self -> getData('Jex')
    Jez = self -> getData('Jez')
    Bx  = self -> getData('Bx')
    Bz  = self -> getData('Bz')
    
    ;Calculate the Current Density
    JexB_y = (Jex*Bz - Jez*Bx) / n_e
    
    return, JexB_y
end


;+
;   The purpose of this program is to calculate the z-component of the JxB force::
;       (J_{e} \times B)_{z} = \frac{1} {n_{e}} (J_{e,x} B_{y} - J_{e,y} B_{x})
;
; :Private:
;
; :Returns:
;       JexB_z:             The z-component of the JxB force.
;-
function MrSim::JexB_z
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
    n_e = self -> getData('ne')
    Jex = self -> getData('Jex')
    Jey = self -> getData('Jey')
    Bx  = self -> getData('Bx')
    By  = self -> getData('By')
    
    ;Calculate the Current Density
    JexB_z = -(Jex*By - Jey*Bx) / n_e
    
    return, JexB_z
end


;+
;   The purpose of this program is to calculate the z-component of the Electron
;   Current Density::
;       (J_{i})_{xz} = n_{i} * \sqrt{ U_{ix}^{2} + U_{iz}^{2} }
;
; :Private:
;
; :Returns:
;       Ji_xz:          The in-plane Electron Current Density.
;-
function MrSim::Ji_xz
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density  
    Uix = self -> getData('Uix')
    Uiz = self -> getData('Uiz')
    n_i = self -> getData('ni')
    
    ;Calculate the in-plane Current Density
    Ui_xz=sqrt(Uix^2 + Uiz^2)
    Ji_xz = n_i * Ui_xz
    
    return, Ji_xz
end


;+
;   The purpose of this program is to calculate the x-component of the Ion
;   Current Density::
;       J_{ix} = n_{i} * U_{ix}
;
; :Private:
;
; :Returns:
;       Jix:            The x-component of the Ion Current Density.
;-
function MrSim::Jix
	compile_opt strictarr, hidden
	
	;Get the Ion Velocity and Number Density  
    Uix = self -> getData('Uix')
    ni = self -> getData('ni')
    
    ;Calculate the z-component of the Ion Current Density
    Jix = ni * Uix
    
    return, Jix
end


;+
;   The purpose of this program is to calculate the y-component of the Ion
;   Current Density::
;       J_{iy} = n_{i} * U_{iy}
;
; :Private:
;
; :Returns:
;       Jiy:                    The y-component of the Ion Current Density.
;-
function MrSim::Jiy
	compile_opt strictarr, hidden
	
	;Get the Ion Velocity and Number Density  
    Uiy = self -> getData('Uiy')
    ni = self -> getData('ni')
    
    ;Calculate the y-component of the Ion Current Density
    Jiy = ni * Uiy
    
    return, Jiy
end


;+
;   The purpose of this program is to calculate the z-component of the Ion
;   Current Density::
;       J_{iz} = n_{i} * U_{iz}
;
; :Private:
;
; :Returns:
;       Jiz:                    The z-component of the Ion Current Density.
;-
function MrSim::Jiz
	compile_opt strictarr, hidden
	
	;Get the Ion Velocity and Number Density  
    Uiz = self -> getData('Uiz')
    ni = self -> getData('ni')
    
    ;Calculate the z-component of the Ion Current Density
    Jiz = ni * Uiz
    
    return, Jiz
end


;+
;   The purpose of this program is to calculate the x-component of the JxB portion
;   of the Generalized Ohm's Law::
;       (J_{i} \times B)_{x} = \frac{1} {n_{i}} (J_{i,y} B_{z} - J_{i,z} B_{y})
;
; :Private:
;
; :Returns:
;       JixB_x:             The x-component of the JxB force.
;-
function MrSim::JixB_x
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
    n_i = self -> getData('ni')
    Jiy = self -> getData('Jiy')
    Jiz = self -> getData('Jiz')
    By  = self -> getData('By')
    Bz  = self -> getData('Bz')
    
    ;Calculate the Current Density
    JixB_x = (Jiy*Bz - Jiz*By) / n_i
    
    return, JixB_x
end


;+
;   The purpose of this program is to calculate the y-component of the JxB force::
;       (J_{i} \times B)_{y} = \frac{1} {e n_{i}} (J_{i,x} B_{z} - J_{i,z} B_{x})
;
; :Private:
;
; :Returns:
;       JixB_y:             The x-component of the JxB force.
;-
function MrSim::JixB_y
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
    n_i = self -> getData('ni')
    Jix = self -> getData('Jix')
    Jiz = self -> getData('Jiz')
    Bx  = self -> getData('Bx')
    Bz  = self -> getData('Bz')
    
    ;Calculate the Current Density
    JixB_y = (Jix*Bz - Jiz*Bx) / n_i
    
    return, JixB_y
end


;+
;   The purpose of this program is to calculate the z-component of the JxB force::
;       (J_{i} \times B)_{z} = \frac{1} {n_{i}} (J_{i,x} B_{y} - J_{i,y} B_{x})
;
; :Private:
;
; :Returns:
;       JixB_z:             The z-component of the JxB force.
;-
function MrSim::JixB_z
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
    n_i = self -> getData('ne')
    Jix = self -> getData('Jix')
    Jiy = self -> getData('Jiy')
    Bx  = self -> getData('Bx')
    By  = self -> getData('By')
    
    ;Calculate the Current Density
    JixB_z = -(Jix*By - Jiy*Bx) / n_i
    
    return, JixB_z
end


;+
;   The purpose of this program is to calculate the Current Density parallel to the
;   magnetic field::
;       B_hat = [Bx, By, Bz] / |B|
;       J_para = J dot B_hat = Jx*Bx_hat + Jy*By_hat + Jz*Bz_hat
;
; :Private:
;
; :Returns:
;       J_PARA:                 The electric field strength in the parallel-to-B
;                                   direction.
;-
function MrSim::Jpar
	compile_opt strictarr, hidden
	
	;Read the electric and magnetic field data.
    Jx = self -> getData('Jx')
    Jy = self -> getData('Jy')
    Jz = self -> getData('Jz')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate a unit vector pointing in the direction of the magnetic field
    Bx_unit = Bx/SQRT(Bx^2+By^2+Bz^2)
    By_unit = By/SQRT(Bx^2+By^2+Bz^2)
    Bz_unit = Bz/SQRT(Bx^2+By^2+Bz^2)
    
    ;Dot the Current Density with B_had
    J_para = Jx*Bx_unit + Jy*By_unit + Jz*Bz_unit
    
    return, J_para
end


;+
;   The purpose of this program is to calculate the x-component of the JxB portion
;   of the Generalized Ohm's Law::
;       (\vec{J} \times \vec{B})_{x} = \frac{1} {n_{i} + n_{e}} (J_{y} B_{z} - J_{z} B_{y})
;
; :Private:
;
; :Returns:
;       JxB_x:             The x-component of the JxB force.
;-
function MrSim::JxB_x
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
	n_i = self -> getData('ni')
	n_e = self -> getData('ne')
    Jy = self -> getData('Jy')
    Jz = self -> getData('Jz')
    By  = self -> getData('By')
    Bz  = self -> getData('Bz')
    
    ;Calculate the Current Density
    JxB_x = (Jy*Bz - Jz*By) / (n_i + n_e)
    
    return, JxB_x
end


;+
;   The purpose of this program is to calculate the y-component of the JxB force::
;       (\vec{J} \times \vec{B})_{y} = \frac{1} {n_{i} + n_{e}} (J_{z} B_{x} - J_{x} B_{z})
;
; :Private:
;
; :Returns:
;       JxB_y:             The x-component of the JxB force.
;-
function MrSim::JxB_y
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
	n_i = self -> getData('ni')
	n_e = self -> getData('ne')
    Jx = self -> getData('Jx')
    Jz = self -> getData('Jz')
    Bx  = self -> getData('Bx')
    Bz  = self -> getData('Bz')
    
    ;Calculate the Current Density
    JxB_y = (Jx*Bz - Jz*Bx) / (n_i + n_e)
    
    return, JxB_y
end


;+
;   The purpose of this program is to calculate the z-component of the JxB force::
;       (\vec{J} \times \vec{B})_{z} = \frac{1} {n_{i} + n_{e}} (J_{x} B_{y} - J_{y} B_{x})
;
; :Private:
;
; :Returns:
;       JxB_z:             The z-component of the JxB force.
;-
function MrSim::JxB_z
	compile_opt strictarr, hidden
	
	;Get the Electron Velocity and Number Density
	n_i = self -> getData('ni')
	n_e = self -> getData('ne')
    Jx = self -> getData('Jx')
    Jy = self -> getData('Jy')
    Bx  = self -> getData('Bx')
    By  = self -> getData('By')
    
    ;Calculate the Current Density
    JxB_z = (Jx*By - Jy*Bx) / (n_i + n_e)
    
    return, JxB_z
end


;+
;   The purpose of this program is to calculate the x-component of the current density::
;       J_{x} = ( -n_{e} * U_{ex} ) + ( n_{i} * U_{ix} )
;
; :Private:
;
; :Returns:
;       Jx:                     The electric field dotted with the current density.
;-
function MrSim::Jx
    compile_opt strictarr, hidden

    ;Get the electric field and current density data
    n_i = self -> getData('ni')
    n_e = self -> getData('ne')
    Uix = self -> getData('Uix')
    Uex = self -> getData('Uex')

    ;Calculate the dot product
    Jx = (-n_e * Uex) + (n_i * Uix)
    
    return, Jx
end


;+
;   The purpose of this program is to calculate the y-component of the current density::
;       J_{y} = ( -n_{e} * U_{ey} ) + ( n_{i} * U_{iy} )
;
; :Private:
;
; :Returns:
;       Jy:                     The electric field dotted with the current density.
;-
function MrSim::Jy
    compile_opt strictarr, hidden

    ;Get the electric field and current density data
    n_i = self -> getData('ni')
    n_e = self -> getData('ne')
    Uiy = self -> getData('Uiy')
    Uey = self -> getData('Uey')

    ;Calculate the dot product
    Jy = (-n_e * Uey) + (n_i * Uiy)
    
    return, Jy
end


;+
;   The purpose of this program is to calculate the z-component of the current density::
;       J_{z} = ( -n_{e} * U_{ez} ) + ( n_{i} * U_{iz} )
;
; :Private:
;
; :Returns:
;       Jz:                     The electric field dotted with the current density.
;-
function MrSim::Jz
    compile_opt strictarr, hidden

    ;Get the electric field and current densitz data
    n_i = self -> getData('ni')
    n_e = self -> getData('ne')
    Uiz = self -> getData('Uiz')
    Uez = self -> getData('Uez')

    ;Calculate the dot product
    Jz = (-n_e * Uez) + (n_i * Uiz)
    
    return, Jz
end


;+
;   The purpose of this program is to calculate the electron inertial length::
;       \lambda_{e} = c / \omega_{p,e}
;
; :Private:
;
; :Returns:
;       LAMBDA_e:               Electron inertial length (skin depth).
;-
function MrSim::lambda_e
    compile_opt strictarr, hidden

    ;Get the electric field and current densitz data
    c = constants('c')
    w_pe = self -> getData('w_pe')

    ;Calculate the dot product
    lambda_e = 1.0 / w_pe
    
    return, lambda_e
end


;+
;   The purpose of this program is to calculate the electron inertial length::
;       \lambda_{i} = c / \omega_{p,e}
;
; :Private:
;
; :Returns:
;       LAMBDA_i:               Ion inertial length (skin depth).
;-
function MrSim::lambda_i
    compile_opt strictarr, hidden
    on_error, 2

    ;Get the electric field and current densitz data
    c = constants('c')
    w_pe = self -> getData('w_pi')

    ;Calculate the dot product
    lambda_i = c / w_pi
    
    return, lambda_i
end


;+
;   Calculate the total pressure::
;
;       P_{tot} = \frac{|B|^2} {2} + \frac{1}{3} (P_{e,xx}^{2} + P_{e,yy}^{2} + P_{e,zz}^{2}) + 
;                                    \frac{1}{3} (P_{i,xx}^{2} + P_{i,yy}^{2} + P_{i,zz}^{2})
;-
function MrSim::P_tot
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe-xx')
    Pe_yy = self -> GetData('Pe-yy')
    Pe_zz = self -> GetData('Pe-zz')
    Pi_xx = self -> GetData('Pi-xx')
    Pi_yy = self -> GetData('Pi-yy')
    Pi_zz = self -> GetData('Pi-zz')
    
    ;Magnitude and direction of magnetic field
    Pb = (Bx^2 + By^2 + Bz^2) / 2
    Pe = (Pe_xx + Pe_yy + Pe_zz) / 3.0
    Pi = (Pi_xx + Pi_yy + Pi_zz) / 3.0
    
    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    P_tot = Pb + Pe + Pi
    
    return, P_tot
end

;+
;   Calculate the total pressure::
;
;       P_{tot} = \frac{|B|^2} {2} 
;-
function MrSim::Pb
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    
    ;Magnitude and direction of magnetic field
    P_b = (Bx^2 + By^2 + Bz^2) / 2
    
    return, P_b
end


;+
;   Calculate the total pressure::
;
;       P_{tot} = \frac{|B|^2} {2} + \frac{1}{3} (P_{e,xx}^{2} + P_{e,yy}^{2} + P_{e,zz}^{2})
;-
function MrSim::Pbe
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe-xx')
    Pe_yy = self -> GetData('Pe-yy')
    Pe_zz = self -> GetData('Pe-zz')
    
    ;Magnitude and direction of magnetic field
    Pb = (Bx^2 + By^2 + Bz^2) / 2
    Pe = (Pe_xx + Pe_yy + Pe_zz) / 3.0
    
    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    P_be = Pb + Pe
    
    return, P_be
end


;+
;   Calculate the total pressure::
;
;       P_{tot} = \frac{|B|^2} {2} + \frac{1}{3} (P_{i,xx}^{2} + P_{i,yy}^{2} + P_{i,zz}^{2})
;-
function MrSim::Pbi
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    Pi_xx = self -> GetData('Pi-xx')
    Pi_yy = self -> GetData('Pi-yy')
    Pi_zz = self -> GetData('Pi-zz')
    
    ;Magnitude and direction of magnetic field
    Pb = (Bx^2 + By^2 + Bz^2) / 2
    Pi = (Pi_xx + Pi_yy + Pi_zz) / 3.0
    
    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    P_bi = Pb + Pi
    
    return, P_bi
end


;+
;   The purpose of this program is to calculate the scalar pressure.
;       P_{e} = 1/3 (P_{e,xx} + P_{e,yy} + P_{e,zz})
;
; :Private:
;
; :Returns:
;       Pe:         The x-component of the divergence of the electron pressure
;                       tensor.
;-
function MrSim::Pe
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the diagonal components of the pressure tensor
    Pe_xx = self -> getData('Pe-xx')
    Pe_yy = self -> getData('Pe-yy')
    Pe_zz = self -> getData('Pe-zz')

    ;Compute the scalar pressure
    Pe = (Pe_xx + Pe_yy + Pe_zz) / 3.0
        
    return, Pe
end


;+
;   The purpose of this program is to calculate the scalar pressure.
;       P_{i} = 1/3 (P_{i,xx} + P_{i,yy} + P_{i,zz})
;
; :Private:
;
; :Returns:
;       Pi:         The x-component of the divergence of the electron pressure
;                       tensor.
;-
function MrSim::Pi
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the diagonal components of the pressure tensor
    Pi_xx = self -> getData('Pi-xx')
    Pi_yy = self -> getData('Pi-yy')
    Pi_zz = self -> getData('Pi-zz')

    ;Compute the scalar pressure
    Pi = (Pi_xx + Pi_yy + Pi_zz) / 3.0
        
    return, Pi
end



;+
;   Calculate the parallel electron pressure::
;
;       P_{\parallel} = \integral (\vec{v_{e}} \cdot \vec{B}) (\vec{U_{e}} \cdot \vec{B}) d^{3}v
;-
function MrSim::Pe_par
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe-xx')
    Pe_xy = self -> GetData('Pe-xy')
    Pe_xz = self -> GetData('Pe-xz')
    Pe_yy = self -> GetData('Pe-yy')
    Pe_yz = self -> GetData('Pe-yz')
    Pe_zz = self -> GetData('Pe-zz')
    
    ;Magnitude and direction of magnetic field
    Bmag  = sqrt(Bx^2 + By^2 + Bz^2)
    bx_hat = temporary(Bx) / Bmag
    by_hat = temporary(By) / Bmag
    bz_hat = temporary(Bz) / temporary(Bmag)
    
    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    P_par = Pe_xx * bx_hat^2 + 2.0 * Pe_xy * bx_hat * by_hat + 2.0 * Pe_xz * bx_hat * bz_hat + $
            Pe_yy * by_hat^2 + 2.0 * Pe_yz * by_hat * bz_hat + $
            Pe_zz * bz_hat^2
    
    return, P_par
end


;+
;   Calculate the perpendicular electron pressure, assuming the plasma is nearly
;   gyrotropic::
;
;       Tr(P) = P_{\parallel} + 2 P_{\bot}
;
; :Hidden:
;-
function MrSim::Pe_perp
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Pe_xx = self -> GetData('Pe-xx')
    Pe_yy = self -> GetData('Pe-yy')
    Pe_zz = self -> GetData('Pe-zz')
    P_par = self -> GetData('Pe_par')
    
    ;Perpendicular pressure
    ;   - Assuming the pressure tensor is mostly gyrotropic, P = P_par + 2 P_perp
    P_perp = (Pe_xx + Pe_yy + Pe_zz - P_par) / 2.0
    
    return, P_perp
end


;+
;   Compute the electron anisotropy, valid where the plasma is mostly isotropic::
;
;       An_{e} = \frac{ P_{e \bot} } { P_{e \parallel} }
;-
function MrSim::An_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    P_par  = self -> GetData('Pe_par')
    P_perp = self -> GetData('Pe_perp')
    
    ;Compute the anisotropy
    e_An = P_par / P_perp
    
    return, e_An
end


;+
;   Compute the electron Agyrotropy
;
;       An_{e} = \frac{ P_{e \bot} } { P_{e \parallel} }
;
;   Reference::
;       Scudder, J., and W. Daughton (2008), “Illuminating” electron diffusion regions of
;           collisionless magnetic reconnection using electron agyrotropy, 
;           J. Geophys. Res., 113, A06222, doi:10.1029/2008JA013035.
;-
function MrSim::A0_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe-xx')
    Pe_xy = self -> GetData('Pe-xy')
    Pe_xz = self -> GetData('Pe-xz')
    Pe_yy = self -> GetData('Pe-yy')
    Pe_yz = self -> GetData('Pe-yz')
    Pe_zz = self -> GetData('Pe-zz')
    
    ;Magnitude and direction of magnetic field
    Bmag  = sqrt(Bx^2 + By^2 + Bz^2)
    bx = temporary(Bx) / Bmag
    by = temporary(By) / Bmag
    bz = temporary(Bz) / temporary(Bmag)
    
    ;Weighted average of dispersions of velocity vectors perpendicular to the magnetic
    ;field in the electron's zero moment frame.
    Nxx =  by*by*Pe_zz - 2.0*by*bz*Pe_yz + bz*bz*Pe_yy
    Nxy = -by*bx*Pe_zz +     by*bz*Pe_xz + bz*bx*Pe_yz - bz*bz*Pe_xy
    Nxz =  by*bx*Pe_yz -     by*by*Pe_xz - bz*bx*Pe_yy + bz*by*Pe_xy
    Nyy =  bx*bx*Pe_zz - 2.0*bx*bz*Pe_xz + bz*bz*Pe_xx
    Nyz = -bx*bx*Pe_yz +     bx*by*Pe_xz + bz*bx*Pe_xy - bz*by*Pe_xx
    Nzz =  bx*bx*Pe_yy - 2.0*bx*by*Pe_xy + by*by*Pe_xx

    ;Perpendicular directions
    alpha_e = Nxx + Nyy + Nzz
    beta_e  = -(Nxy^2 + Nxz^2 + Nyz^2 - Nxx*Nyy - Nxx*Nzz - Nyy*Nzz)

    ; Return agyrotropy data
    A0 = 2.0 * sqrt(alpha_e^2 - 4.0*beta_e) / alpha_e
    return, A0
end


;+
;   Compute the degree of nongyrotropy
;
;       An_{e} = \frac{ P_{e \bot} } { P_{e \parallel} }
;
;   Reference::
;       Aunai, N., Hesse, M. and Kuznetsova, M., Electron nongyrotropy in the context of
;           collisionless magnetic reconnection, Phys. Plasmas 20, 092903 (2013)
;           http://dx.doi.org/10.1063/1.4820953
;-
function MrSim::Dng_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx    = self -> GetData('Bx')
    By    = self -> GetData('By')
    Bz    = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe-xx')
    Pe_xy = self -> GetData('Pe-xy')
    Pe_xz = self -> GetData('Pe-xz')
    Pe_yy = self -> GetData('Pe-yy')
    Pe_yz = self -> GetData('Pe-yz')
    Pe_zz = self -> GetData('Pe-zz')
    
    ;Magnitude and direction of magnetic field
    Bmag  = sqrt(Bx^2 + By^2 + Bz^2)
    bx_hat = temporary(Bx) / Bmag
    by_hat = temporary(By) / Bmag
    bz_hat = temporary(Bz) / temporary(Bmag)
    
    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    P_par = Pe_xx * bx_hat^2 + 2.0 * Pe_xy * bx_hat * by_hat + 2.0 * Pe_xz * bx_hat * bz_hat + $
            Pe_yy * by_hat^2 + 2.0 * Pe_yz * by_hat * bz_hat + $
            Pe_zz * bz_hat^2

    ;Perpendicular pressure
    ;   - Assuming the pressure tensor is mostly gyrotropic, P = P_par + 2 P_perp
    P_perp = (Pe_xx + Pe_yy + Pe_zz - P_par) / 2.0
    
    ; Compute G in the x-y-z basis, using the definition in equation (1) of Aunai et al., 2013, POP:
    ; G = pper*I + (ppar - pper)*bb, where I is the 3x3 identity matrix, and bb is the outer product
    ; of the normalized magnetic field vector b. In the x-y-z basis, b = [bx, by, bz], and G is not 
    ; gyrotropic in general, while in the frame parallel to the magnetic field, b is simply [0, 0, 1],
    ; so G would be gyrotropic in this frame.
    Gxx = P_perp + (P_par - P_perp) * bx_hat^2
    Gyy = P_perp + (P_par - P_perp) * by_hat^2
    Gzz = P_perp + (P_par - P_perp) * bz_hat^2
    Gxy = (P_par - P_perp) * bx_hat * by_hat 
    Gxz = (P_par - P_perp) * bx_hat * bz_hat
    Gyz = (P_par - P_perp) * by_hat * bz_hat
    
    ;Free memory
    P_perp = !Null
    P_par  = !Null
    bx_hat = !Null
    by_hat = !Null
    bz_hat = !Null
    
    ; Now compute N, which is just P - G.
    Nxx = Pe_xx - Gxx
    Nyy = Pe_yy - Gyy
    Nzz = Pe_zz - Gzz
    Nxy = Pe_xy - Gxy
    Nxz = Pe_xz - Gxz
    Nyz = Pe_yz - Gyz
    
    ;Free memory
    Gxx   = !Null
    Gxy   = !Null
    Gxz   = !Null
    Gyy   = !Null
    Gyz   = !Null
    Gzz   = !Null
    Pe_xy = !Null
    Pe_xz = !Null
    Pe_yz = !Null

    ; Next, compute Dng (the measure of nongyrotropy), a scalar for each grid point.
    Dng = sqrt(Nxx^2 + Nyy^2 + Nzz^2 + 2*Nxy^2 + 2*Nxz^2 + 2*Nyz^2) / (Pe_xx + Pe_yy + Pe_zz)
    return, Dng
end



;+
;   The purpose of this program is to calculate the magnitude squared of the electron
;   velocity::
;       \left| U_{e} \right| = \sqrt{ U_{e} \cdot U_{e} } = \sqrt{ U_{ex}^{2} + U_{ey}^{2} + U_{ez}^{2} }
;       \vec{U}_{e} = \sqrt{ \vec{U}_{e} \cdot \vec{U}_{e} } = \sqrt{ U_{ex}^{2} + U_{ey}^{2} + U_{ez}^{2} }
;
; :Private:
;
; :Returns:
;       Ue_mag:         The magnitude of the electron velocity.
;-
function MrSim::Ue_mag
	compile_opt strictarr, hidden
	
	;Get the electron velocity
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    
    ;Calculate |Ue| = sqrt(Ue dot Ue)
    Ue_mag = sqrt(Uex^2 + Uey^2 + Uez^2)
    
    return, Ue_mag
end



;+
;   Calculate the magnitude of the parallel electron velocity::
;       U_{e \parallel} = \left| \frac{ \vec{U}_{e} \cdot \vec{B} } { |B| } \right|
;
; :Private:
;
; :Returns:
;       Ue_par:         Electron velocity parallel to the magnetic field.
;-
function MrSim::Ue_par
	compile_opt strictarr, hidden
	
	;Get the electron velocity
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    
    ;Get the magnetic field
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the direction of B (B_hat)
    B  = sqrt(Bx^2 + By^2 + Bz^2)
    Bx /= B
    By /= B
    Bz /= B
    
    ;Dot Ui into B_hat
    Ue_par = Uex*Bx + Uey*By + Uez*Bz
    return, Ue_par
end



;+
;   Calculate the magnitude of the parallel electron velocity::
;       U_{e \perp} = | \vec{U}_{e} - \vec{U}_{e \parallel} |
;
; :Private:
;
; :Returns:
;       Ue_perp:         Electron velocity perpendicular to the magnetic field.
;-
function MrSim::Ue_perp
	compile_opt strictarr, hidden
	
	;Get the electron velocity
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    
    ;Get the magnetic field
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Use only the direction of B (B_hat)
    B  = sqrt(Bx^2 + By^2 + Bz^2)
    Bx /= B
    By /= B
    Bz /= temporary(B)
    
    ;Calculate the 3-components of the parallel velocity
    Ue_par   = Uex*Bx + Uey*By + Uez*Bz
    Ue_par_x = Ue_par * temporary(Bx)
    Ue_par_y = Ue_par * temporary(By)
    Ue_Par_z = temporary(Ue_par) * temporary(Bz)
    
    ;Compute the 3-components of the perpendicular velocity
    Ue_perp_x = temporary(Uex) - temporary(Ue_par_x)
    Ue_perp_y = temporary(Uey) - temporary(Ue_par_y)
    Ue_perp_z = temporary(Uez) - temporary(Ue_par_z)
    
    ;Magnitude
    Ue_perp = sqrt(Ue_perp_x^2 + Ue_perp_y^2 + Ue_perp_z^2)
    return, Ue_perp
end


;+
;   The purpose of this program is to calculate the magnitude of the cross product 
;   between the bulk ion velocity and Magnetic Field::
;       UexB_x = Uey*Bz - Uez*By
;       UexB_y = Uez*Bx - Uex*Bz
;       UexB_z = Uex*By - Uey*Bx
;       UexB = UexB_x + UexB_y + UexB_z
;       |UexB| = sqrt(UexB dot UexB) = sqrt( UexB_x^2 + UexB_y^2 + UexB_z^2 )
;
; :Private:
;
; :Returns:
;       UexB_mag:       The magnitude of the cross product between the electron bulk
;                           velocity and the magnetic field.
;-
function MrSim::UexB_mag
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    UexB_x  = Uey*Bz - Uez*By
    UexB_y  = Uez*Bx - Uex*Bz
    UexB_z  = Uex*By - Uey*Bx
    UexB_mag = sqrt(UexB_x^2 + UexB_y^2 + UexB_z^2)
    
    return, UexB_mag
end


;+
;   The purpose of this program is to calculate the x-component of the cross product 
;   between the bulk electron velocity and Magnetic Field::
;       UexB_x = Uey*Bz - Uez*By
;
; :Private:
;
; :Returns:
;       UexB_x:         The x-component of the cross product between the bulk electron
;                           velocity and the magnetic field.
;-
function MrSim::UexB_x
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uey = self -> getData('Uey')
    Uez = self -> getData('Uez')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    UexB_x  = Uey*Bz - Uez*By
    
    return, UexB_x
end


;+
;   The purpose of this program is to calculate the y-component of the cross product 
;   between the bulk electron velocity and Magnetic Field::
;       UexB_y = Uez*Bx - Uex*Bz
;
; :Private:
;
; :Returns:
;       UexB_y:         The y-component of the cross product between the bulk electron
;                           velocity and the magnetic field.
;-
function MrSim::UexB_y
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uex = self -> getData('Uex')
    Uez = self -> getData('Uez')
    Bx = self -> getData('Bx')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    UexB_y  = Uez*Bx - Uex*Bz
    
    return, UexB_y
end


;+
;   The purpose of this program is to calculate the z-component of the cross product 
;   between the bulk electron velocity and Magnetic Field::
;       UexB_z = Uex*By - Ey*Bx
;
; :Private:
;
; :Returns:
;       UexB_z:         The z-component of the cross product between the bulk electron
;                           velocity and the magnetic field.
;-
function MrSim::UexB_z
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uex = self -> getData('Uex')
    Uey = self -> getData('Uey')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    
    ;Calculate the magnitude of ExB
    UexB_z  = Uex*By - Uey*Bx
    
    return, UexB_z
end


;+
;   The purpose of this program is to calculate the magnitude squared of the electron
;   velocity::
;       \left| U_{i} \right| = \sqrt{ U_{i} \cdot U_{i} } = \sqrt{ U_{i,x}^{2} + U_{i,y}^{2} + U_{i,z}^{2} }
;       \vec{U}_{i} = \sqrt{ \vec{U}_{i} \cdot \vec{U}_{i} } = \sqrt{ U_{i,x}^{2} + U_{i,y}^{2} + U_{i,z}^{2} }
;
; :Private:
;
; :Returns:
;       Ui_mag:         The magnitude of the ion velocity.
;-
function MrSim::Ui_mag
	compile_opt strictarr, hidden
	
	;Get the electron velocity
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    
    ;Calculate |Ue| = sqrt(Ue dot Ue)
    Ui_mag = sqrt(Uix^2 + Uiy^2 + Uiz^2)
    
    return, Ui_mag
end


;+
;   Calculate the magnitude of the parallel ion velocity::
;       U_{i \parallel} = \left| \frac{ \vec{U} \cdot \vec{B} } { |B| } \right|
;
; :Private:
;
; :Returns:
;       Ui_mag:         Ion velocity parallel to the magnetic field.
;-
function MrSim::Ui_par
	compile_opt strictarr, hidden
	
	;Get the electron velocity
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    
    ;Get the magnetic field
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the direction of B (B_hat)
    B  = sqrt(Bx^2 + By^2 + Bz^2)
    Bx /= B
    By /= B
    Bz /= B
    
    ;Dot Ui into B_hat
    Ui_par = Uix*Bx + Uiy*By + Uiz*Bz
    return, Ui_par
end


;+
;   Calculate the magnitude of the perpendicular ion velocity::
;       U_{i \perp} = | \vec{U}_{i} - \vec{U}_{i\parallel} |
;
; :Private:
;
; :Returns:
;       Ui_perp:         Ion velocity perpendicular to the magnetic field.
;-
function MrSim::Ui_perp
	compile_opt strictarr, hidden
	
	;Get the electron velocity
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    
    ;Get the magnetic field
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Use only the direction of B (B_hat)
    B  = sqrt(Bx^2 + By^2 + Bz^2)
    Bx /= B
    By /= B
    Bz /= temporary(B)
    
    ;Calculate the 3-components of the parallel velocity
    Ui_par = Uix*Bx + Uiy*By + Uiz*Bz
    Ui_par_x = Ui_par * temporary(Bx)
    Ui_par_y = Ui_par * temporary(By)
    Ui_Par_z = temporary(Ui_par) * temporary(Bz)
    
    ;Compute the 3-components of the perpendicular velocity
    Ui_perp_x = temporary(Uix) - temporary(Ui_par_x)
    Ui_perp_y = temporary(Uiy) - temporary(Ui_par_y)
    Ui_perp_z = temporary(Uiz) - temporary(Ui_par_z)
    
    ;Magnitude
    Ui_perp = sqrt(Ui_perp_x^2 + Ui_perp_y^2 + Ui_perp_z^2)
    return, Ui_perp
end


;+
;   The purpose of this program is to calculate the magnitude of the cross product 
;   between the bulk ion velocity and Magnetic Field::
;       UixB_x = Uiy*Bz - Uiz*By
;       UixB_y = Uiz*Bx - Uix*Bz
;       UixB_z = Uix*By - Uiy*Bx
;       UixB = UixB_x + UixB_y + UixB_z
;       |UixB| = sqrt(UixB dot UixB) = sqrt( UixB_x^2 + UixB_y^2 + UixB_z^2 )
;
; :Private:
;
; :Returns:
;       UixB_mag:       The magnitude of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::UixB_mag
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    UixB_x  = Uiy*Bz - Uiz*By
    UixB_y  = Uiz*Bx - Uix*Bz
    UixB_z  = Uix*By - Uiy*Bx
    UixB_mag = sqrt(UixB_x^2 + UixB_y^2 + UixB_z^2)
    
    return, UixB_mag
end


;+
;   The purpose of this program is to calculate the x-component of the cross product 
;   between the bulk ion velocity and Magnetic Field::
;       UixB_x = Uiy*Bz - Uiz*By
;
; :Private:
;
; :Returns:
;       UixB_x:         The x-component of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::UixB_x
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uiy = self -> getData('Uiy')
    Uiz = self -> getData('Uiz')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    UixB_x  = Uiy*Bz - Uiz*By
    
    return, UixB_x
end


;+
;   The purpose of this program is to calculate the y-component of the cross product 
;   between the bulk ion velocity and Magnetic Field::
;       UixB_y = Uiz*Bx - Uix*Bz
;
; :Private:
;
; :Returns:
;       UixB_mag:       The y-component of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::UixB_y
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uix = self -> getData('Uix')
    Uiz = self -> getData('Uiz')
    Bx = self -> getData('Bx')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of ExB
    UixB_y  = Uiz*Bx - Uix*Bz
    
    return, UixB_y
end


;+
;   The purpose of this program is to calculate the z-component of the cross product 
;   between the bulk ion velocity and Magnetic Field::
;       UixB_z = Uix*By - Ey*Bx
;
; :Private:
;
; :Returns:
;       UixB_mag:       The magnitude of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::UixB_z
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
    Uix = self -> getData('Uix')
    Uiy = self -> getData('Uiy')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    
    ;Calculate the magnitude of ExB
    UixB_z  = Uix*By - Uiy*Bx
    
    return, UixB_z
end


;+
;   The purpose of this program is to calculate the x-component of the cross product 
;   between the center-of-mass velocity and the magnetic field::
;       (\vec{V} \times \vec{B})_x = V_{y} B_{z} - V_{z} B_{y}
;
;   where (j = x, y, z)::
;       V_{j} = \frac{m_{i} U_{i,j} + m_{e} U_{e,j}} {m_{i} + m_{e}}
;
; :Private:
;
; :Returns:
;       VxB_x:          The x-component of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::VxB_x
	compile_opt strictarr, hidden
	
	;Get the MHD velocity and magnetic field data
    Vy = self -> getData('Vy')
    Vz = self -> getData('Vz')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of VxB
    VxB_x  = Vy*Bz - Vz*By
    
    return, VxB_x
end


;+
;   The purpose of this program is to calculate the y-component of the cross product 
;   between the center-of-mass velocity and the magnetic field::
;       (\vec{V} \times \vec{B})_y = V_{x} B_{z} - V_{z} B_{x}
;
;   where (j = x, y, z)::
;       V_{j} = \frac{m_{i} U_{i,j} + m_{e} U_{e,j}} {m_{i} + m_{e}}
;
; :Private:
;
; :Returns:
;       VxB_y:          The x-component of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::VxB_y
	compile_opt strictarr, hidden
	
	;Get the MHD velocity and magnetic field data
    Vx = self -> getData('Vx')
    Vz = self -> getData('Vz')
    Bx = self -> getData('Bx')
    Bz = self -> getData('Bz')
    
    ;Calculate the magnitude of VxB
    VxB_y  = Vx*Bz - Vz*Bx
    
    return, VxB_y
end


;+
;   The purpose of this program is to calculate the z-component of the cross product 
;   between the center-of-mass velocity and the magnetic field::
;       (\vec{V} \times \vec{B})_z = V_{x} B_{y} - V_{y} B_{x}
;
;   where (j = x, y, z)::
;       V_{j} = \frac{m_{i} U_{i,j} + m_{e} U_{e,j}} {m_{i} + m_{e}}
;
; :Private:
;
; :Returns:
;       VxB_z:          The z-component of the cross product between the Electric
;                           and Magnetic Fields.
;-
function MrSim::VxB_z
	compile_opt strictarr, hidden
	
	;Get the MHD velocity and magnetic field data
    Vx = self -> getData('Vx')
    Vy = self -> getData('Vy')
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    
    ;Calculate the magnitude of VxB
    VxB_z  = Vx*By - Vy*Bx
    
    return, VxB_z
end


;+
;   The purpose of this program is to calculate the x-component of the MHD, single-
;   fluid, center-of-mass velocity and the magnetic field::
;       V_{x} = \frac{m_{i} U_{i,x} + m_{e} U_{e,x}} {m_{i} + m_{e}}
;
; :Private:
;
; :Returns:
;       Vx:         The x-component of the cross product between the Electric
;                       and Magnetic Fields.
;-
function MrSim::Vx
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
	self -> GetInfo, MI_ME=mi_me
	mi = mi_me
	me = 1.0
    Uix = self -> getData('Uix')
    Uex = self -> getData('Uex')
    
    ;Calculate the magnitude of ExB
    Vx = (mi*Uix + me*Uex) / (mi + me)
    
    return, Vx
end


;+
;   The purpose of this program is to calculate the x-component of the MHD, single-
;   fluid, center-of-mass velocity and the magnetic field::
;       V_{y} = \frac{m_{i} U_{i,y} + m_{e} U_{e,y}} {m_{i} + m_{e}}
;
; :Private:
;
; :Returns:
;       Vy:         The y-component of the cross product between the Electric
;                       and Magnetic Fields.
;-
function MrSim::Vy
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
	self -> GetInfo, MI_ME=mi_me
	mi = mi_me
	me = 1.0
    Uiy = self -> getData('Uiy')
    Uey = self -> getData('Uey')
    
    ;Calculate the magnitude of ExB
    Vy = (mi*Uiy + me*Uey) / (mi + me)
    
    return, Vy
end


;+
;   The purpose of this program is to calculate the x-component of the MHD, single-
;   fluid, center-of-mass velocity and the magnetic field::
;       V_{z} = \frac{m_{i} U_{i,z} + m_{e} U_{e,z}} {m_{i} + m_{e}}
;
; :Private:
;
; :Returns:
;       Vz:         The z-component of the cross product between the Electric
;                       and Magnetic Fields.
;-
function MrSim::Vz
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
	self -> GetInfo, MI_ME=mi_me
	mi = mi_me
	me = 1.0
    Uiz = self -> getData('Uiz')
    Uez = self -> getData('Uez')
    
    ;Calculate the magnitude of ExB
    Vz = (mi*Uiz + me*Uez) / (mi + me)
    
    return, Vz
end


;+
;   The purpose of this program is to calculate the electron Alfven Speed::
;       V_{ae} = \left| B \right| / \sqrt{ne}
;
; :Private:
;
; :Returns:
;       Vae:                    The Alfven Speed
;-
function MrSim::v_Ae
	compile_opt strictarr, hidden
	
	;Get the Magnetic Field and Electron Density 
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    n_e = self -> getData('ne')
    
    ;Calculate the Alfven Speed
    b_mag = sqrt(bx^2 + by^2 + bz^2)
    Vae = b_mag / sqrt(n_e)
    
    return, Vae
end


;+
;   The purpose of this program is to calculate the ion Alfven Speed::
;       V_{A} = \left| B \right| / \sqrt{n e}
;
; :Private:
;
; :Returns:
;       Vae:                    The Alfven Speed
;-
function MrSim::v_A
	compile_opt strictarr, hidden
	
	;Get the Magnetic Field and Electron Density 
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    n_i = self -> getData('ni')
    
    ;Calculate the Alfven Speed
    b_mag = sqrt(bx^2 + by^2 + bz^2)
    Vae = b_mag / sqrt(n_i)
    
    return, Vae
end


;+
;   The purpose of this program is to calculate the electron plasma frequency::
;       w_{p,e} = \sqrt{ \frac{n_{0} e^{2}} {\epsilon_{0} m} }
;
; :Private:
;
; :Returns:
;       w_pe:                   The electron plasma frequency
;-
function MrSim::w_pe
	compile_opt strictarr, hidden
	
	;Get the Magnetic Field and Electron Density
    n_e = self -> getData('ne')
    
    ;Calculate the Alfven Speed
    w_pe = sqrt(n_e)
    
    return, w_pe
end


;+
;   The purpose of this program is to calculate the ion plasma frequency::
;       w_{p,i} = \sqrt{ \frac{n_{0} e^{2}} {\epsilon_{0} m} }
;
; :Private:
;
; :Returns:
;       w_pi:                   The ion plasma frequency
;-
function MrSim::w_pi
	compile_opt strictarr, hidden
	
	;Get the Magnetic Field and Electron Density
    n_i = self -> getData('ni')
    
    ;Calculate the Alfven Speed
    w_pi = sqrt(n_i)
    
    return, w_pi
end


;+
;   The purpose of this program is to calculate the electron cyclotron frequency::
;       w_{c,e} = \frac{e B} {m}
;
; :Private:
;
; :Returns:
;       w_ce:                   The ion cyclotron frequency
;-
function MrSim::w_ci
	compile_opt strictarr, hidden
	
	;Get the Magnetic Field and Electron Density
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the Alfven Speed
    w_ce = sqrt(bx^2 + by^2 + bz^2)
    
    return, w_ce
end


;+
;   The purpose of this program is to calculate the ion cyclotron frequency::
;       w_{c,i} = \frac{e B} {m}
;
; :Private:
;
; :Returns:
;       w_ci:                   The ion cyclotron frequency
;-
function MrSim::w_ci
	compile_opt strictarr, hidden
	
	;Get the Magnetic Field and Electron Density
    Bx = self -> getData('Bx')
    By = self -> getData('By')
    Bz = self -> getData('Bz')
    
    ;Calculate the Alfven Speed
    w_ci = sqrt(bx^2 + by^2 + bz^2)
    
    return, w_ci
end


;+
;   This definition statement for the MrSim class.
;
; :Params:
;       CLASS           out, optional, type=named structure
;                       The class definition structure.
;
; :FIELDS:
;       AXIS_LABELS:    Coordinate labels for the axes (e.g. 'x', 'y', 'z')
;       COORD_SYSTEM:   Coordinate system in which to setup the simulation domain.
;       DIRECTORY:      Directory in which the data files (.gda) are located
;       INFO:           Information about the simulation setup
;       ION_SCALE:      Sets that the simulation domain to "di"
;       LITTLE:         Are binary files little endian?
;       NSMOOTH:        Data will be smoothed by this amount
;       NLOW_GRID:      Create the simulation grid using density from the low-density side?
;       ORIENTATION:    2D orientation at which data will be displayed.
;       TIME:           Time index at which data will be read
;       XRANGE:         X-range of the simulation domain to be read
;       XSIM:           X-coordinates of the simulation domain
;       ZRANGE:         Z-range of the simulation domain to be read
;       ZSIM:           Z-coordinates of the simulation domain
;       AY:             Data for the flux function
;       Bx:             Data for the x-component of the magnetic field
;       Bx:             Data for the y-component of the magnetic field
;       Bx:             Data for the z-component of the magnetic field
;       Ex:             Data for the x-component of the electric field
;       Ex:             Data for the y-component of the electric field
;       Ex:             Data for the z-component of the electric field
;       n_e:            Data for the electron density
;       n_i:            Data for the ion density
;       Pe_xx:          Data for the [x,x]-component of the electron pressure tensor
;       Pe_xy:          Data for the [x,y]-component of the electron pressure tensor
;       Pe_xz:          Data for the [x,z]-component of the electron pressure tensor
;       Pe_yx:          Data for the [y,x]-component of the electron pressure tensor
;       Pe_yy:          Data for the [y,y]-component of the electron pressure tensor
;       Pe_yz:          Data for the [y,z]-component of the electron pressure tensor
;       Pe_zx:          Data for the [z,x]-component of the electron pressure tensor
;       Pe_zy:          Data for the [z,y]-component of the electron pressure tensor
;       Pe_zz:          Data for the [z,z]-component of the electron pressure tensor
;       Pi_xx:          Data for the [x,x]-component of the ion pressure tensor
;       Pi_xy:          Data for the [x,y]-component of the ion pressure tensor
;       Pi_xz:          Data for the [x,z]-component of the ion pressure tensor
;       Pi_yx:          Data for the [y,x]-component of the ion pressure tensor
;       Pi_yy:          Data for the [y,y]-component of the ion pressure tensor
;       Pi_yz:          Data for the [y,z]-component of the ion pressure tensor
;       Pi_zx:          Data for the [z,x]-component of the ion pressure tensor
;       Pi_zy:          Data for the [z,y]-component of the ion pressure tensor
;       Pi_zz:          Data for the [z,z]-component of the ion pressure tensor
;       Uex:            Data for the x-component of the electron bulk field
;       Uey:            Data for the y-component of the electron bulk field
;       Uez:            Data for the z-component of the electron bulk field
;       Uix:            Data for the x-component of the ion bulk field
;       Uiy:            Data for the y-component of the ion bulk field
;       Uiz:            Data for the z-component of the ion bulk field
;-
pro MrSim__DEFINE, class
    compile_opt strictarr
    
    class = { MrSim, $
              inherits IDL_Object, $
             
              ;Object Properties
              axis_labels:  ['', '', ''], $
              coord_system: '', $
              directory:    '', $
              info:         ptr_new(), $
              ion_scale:    0B, $
              little:       0B, $
              mva_frame:    0B, $
              nSmooth:      0L, $
              orientation:  '', $
              simnum:       0,  $
              simname:      '', $
              time:         0L, $
              
              ;Sim Domain
              XRANGE: [0D, 0D], $
              XSim:   ptr_new(), $
              YSlice: 0L, $
              YRANGE: [0D, 0D], $
              YSim:   ptr_new(), $
              ZRANGE: [0D, 0D], $
              ZSim:   ptr_new(), $
             
              ;Data
              electrons: ptr_new(), $
              Ay:        ptr_new(), $
              Bx:        ptr_new(), $
              By:        ptr_new(), $
              Bz:        ptr_new(), $
              Ex:        ptr_new(), $
              Ey:        ptr_new(), $
              Ez:        ptr_new(), $
              n_e:       ptr_new(), $
              n_i:       ptr_new(), $
              Pe_xx:     ptr_new(), $
              Pe_xy:     ptr_new(), $
              Pe_xz:     ptr_new(), $
              Pe_yx:     ptr_new(), $
              Pe_yy:     ptr_new(), $
              Pe_yz:     ptr_new(), $
              Pe_zx:     ptr_new(), $
              Pe_zy:     ptr_new(), $
              Pe_zz:     ptr_new(), $
              Pi_xx:     ptr_new(), $
              Pi_xy:     ptr_new(), $
              Pi_xz:     ptr_new(), $
              Pi_yx:     ptr_new(), $
              Pi_yy:     ptr_new(), $
              Pi_yz:     ptr_new(), $
              Pi_zx:     ptr_new(), $
              Pi_zy:     ptr_new(), $
              Pi_zz:     ptr_new(), $
              Uex:       ptr_new(), $
              Uey:       ptr_new(), $
              Uez:       ptr_new(), $
              Uix:       ptr_new(), $
              Uiy:       ptr_new(), $
              Uiz:       ptr_new() $
            }
end