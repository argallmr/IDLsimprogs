; docformat = 'rst'
;
;+
;   The purpose of this program is to provide a base class for 2D and 3D simulation
;   data provided by Bill Daughton from Los Alamos National Laboratory.
;
;   MrSim2 is meant to be subclassed by a class that over-rides the ReadData method.
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
;       2014/09/06  -   Written by Matthew Argall
;       2014/11/20  -   Get1DCut returns number of points along cut. - MRA
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
;                       If set, `INFO_FILE` points to the binary info file.
;       COORD_SYSTEM:   in, optional, type=string, default='SIMULATION'
;                       Coordinate system in which to display the data. Options are::
;                           'SIMULATION'
;                           'MAGNETOPAUSE'
;                           'MAGNETOTAIL'
;       DIRECTORY:      in, optional, type=string, default=pwd
;                       Directory in which to find the ".gda" data.
;       INFO_FILE:      in, optional, type=string, default=`DIRECTORY`/../info`
;                       The ASCII info file containing information about the simulation
;                           setup. If `BINARY` is set, the default file will be
;                           `DIRECTORY`/info.
;       INFO_VERSION:   in, optional, type=integer, default=1
;                       Version of the info file to read. Ignored if `BINARY`=1.
;                           See MrSim_Which.pro.
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
function MrSim2::INIT, theSim, time, yslice, $
ASCII_VERSION = ascii_version, $
AXIS_LABELS = axis_labels, $
BINARY = binary, $
COORD_SYSTEM = coord_system, $
DIRECTORY = directory, $
INFO_FILE = info_file, $
INFO_VERSION = info_version, $
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
    MrSim_Which, theSim, NAME=simname, NUMBER=simnum, $
                 DIRECTORY=_dir, INFO_ASCII=_ascii_info, INFO_BINARY=_bin_info, $
                 ASCII_VERSION=vASCII, DIMENSION=dimension, ASYMMETRIC=asymmetric
    
    ;Defaults
    binary = keyword_set(binary)
    if n_elements(directory)     eq 0 then directory     = _dir
    if n_elements(ascii_version) eq 0 then ascii_version = vASCII
    if n_elements(info_ascii)    eq 0 then info_ascii    = _ascii_info
    if n_elements(info_file)     eq 0 then info_binary   = _bin_info
    
    ;Properties
    self.symmetric = ~asymmetric
    self.directory = directory
    self.dimension = dimension
    self.simnum    = simnum
    self.simname   = simname
;-------------------------------------------------------
;Read the Info File ////////////////////////////////////
;-------------------------------------------------------
    
    ;Make sure the file exists
    self.info = ptr_new(/ALLOCATE_HEAP)
    if file_test(info_ascii) eq 0 then $
        message, 'Cannot find ASCII info file: "' + info_ascii + '"'
    
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
    
    ;
    ;From Bill Daughton:
    ;   The 3D data is so big, it is difficult to look at and/or move.
    ;   To make things easier, we usually apply a smoothing operator and then 
    ;   downsample by a factor of 2x in each direction - which gives a 
    ;   reduction of 8x in the file size. This is why nx, ny, nz are 2x 
    ;   smaller in the data I sent. For example, in 1D the smoothing 
    ;   operator would be
    ;
    ;   B_j  =  0.25 *[ B(j-1) + 2*B(j) + B(j+1) ]
    ;
    ;   but instead of keeping every point, I would skip every other one for
    ;   the new "smoothed" array. 
    ;
    
    ;This information is captured in the binary info file, but not the ASCII file.
    if self.dimension eq '3D' && binary eq 0 $
        then self -> ReadInfo_Binary, _bin_info

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
    
    ;Coordinate System
    if max(coord_system eq ['SIMULATION', 'MAGNETOPAUSE', 'MAGNETOTAIL']) eq 0 then $
        message, 'Coordinate system "' + coord_system + '" not recognized.'
    if max(_orientation eq ['XY', 'XZ', 'YZ']) eq 0 then $
        message, 'Orienatation "' + _orientation + '" not recognized.'

    ;Axis labels
    if n_elements(axis_labels) eq 0 then begin
        ;MVA Frame?
        if mva_frame then begin
            case coord_system of
                'SIMULATION':   axis_labels = ['L', 'M', 'N']
                'MAGNETOPAUSE': axis_labels = ['N', 'M', 'L']      ;N is along X-GSE
                'MAGNETOTAIL':  axis_labels = ['L', 'M', 'N']      ;N is along Z-GSE
                else:           ;Do nothing
            endcase
        endif else begin
            axis_labels = ['x', 'y', 'z']
        endelse
    endif

    ;Set Properties
    self.axis_labels  = axis_labels
    self.coord_system = coord_system
    self.directory    = directory
    self.ion_scale    = ion_scale
    self.mva_frame    = mva_frame
    self.nsmooth      = nsmooth
    self.orientation  = _orientation
    self.time         = time
    self.xrange       = xrange
    self.yrange       = yrange
    self.zrange       = zrange

    ;Set the domain coordinates
    self -> SetScale, ION_SCALE=ion_scale
    self -> MakeSimDomain

;-------------------------------------------------------
; Data Properties //////////////////////////////////////
;-------------------------------------------------------
    ;Data Products
    success = self -> MrSim2_Data::Init()
    if ~success then return, 0

;-------------------------------------------------------
; Print Info about the Simulation //////////////////////
;-------------------------------------------------------

    ;Get the information
    self -> GetInfo, NX=nx, NY=ny, NZ=nz, $
                     LX_DE=lx_de, LY_DE=ly_de, LZ_DE=lz_de, $
                     MI_ME=mi_me, DTXWCI=dtxwci
    
    ;Check the ::GetProperty method.
    message, '', /INFORMATIONAL
    print, FORMAT='(%"    DTxWCI          = %0.2f  **hard-coded**")', dtxwci

    ;Determine the number of time steps
    Bx_info = file_info(filepath('Bx.gda', ROOT_DIR=self.directory))
    if Bx_info.exists then begin
        ;Size of one time slice
        ;   - 4-byte floats, plus a leading and trailing character
        rec_size = 4L * (nx * ny * nz + 2L)
        nTimes   = (Bx_info.size / rec_size) - 1
        print, FORMAT='(%"    # t-slices      = %i")', nTimes
    endif
    
    ;Simulation info
    print, FORMAT='(%"    Simulation Size = %i x %i x %i")', nx, ny, nz
    print, FORMAT='(%"    Sim Size (de)   = %i x %i x %i")', lx_de, ly_de, lz_de
    if n_elements(mi_me) gt 0 then print, FORMAT='(%"    m_i / m_e       = %f")', mi_me

    return, 1
end


;+
;   Clean up after the object is destroyed.
;-
pro MrSim2::CLEANUP
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
    ptr_free, self.YSim
    ptr_free, self.ZSim
    
    ;Superclasses
    self -> MrSim2_Data::Cleanup
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
function MrSim2::Cell2Coord, cell, $
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
function MrSim2::Coord2Cell, coord, $
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
;       NAME:           in, required, type=string
;                       Name of a data quantity from which a 1D cut is extracted.
;       R0:             in, required, type=float(2)/float(3)
;                       End point of the line that defines the 1D cut. If 2 elements are
;                           given, it is assumed that the data has coordinates [x,y,z],
;                           but that the y-dimension does not vary.
;       R1:             in, required, type=float(2)/float(3)
;                       End point of the line that defines the 1D cut. If 2 elements are
;                           given, it is assumed that the data has coordinates [x,y,z],
;                           but that the y-dimension does not vary.
;
; :Keywords:
;       COORDS:         out, optional, type=2xN or 3xN float
;                       Coordinates of each point connection `R0` to `R1`. If the
;                           Orientation property is 2D, COORDS will be reduced to the two
;                           relevant dimensions.
;       COUNT:          out, optional, type=long
;                       Number of data points returned.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the ::GetData method.
;
; :Returns:
;       DATA:           The requested data.
;-
function MrSim2::Get1DCut, name, r0, r1, $
COORDS=coords, $
COUNT=count, $
_REF_EXTRA=extra
    compile_opt strictarr
    on_error, 2
    
    ;Make sure three coordinates were given.
    _r0 = n_elements(r0) eq 2 ? [r0[0], 0, r0[1]] : r0
    _r1 = n_elements(r1) eq 2 ? [r1[0], 0, r1[1]] : r1
    
    ;Get the coordinates of the line
    !Null = self -> GridLine(_r0, _r1, COORDS=coords, COUNT=count)
    
    ;Locate the coordinates within the simulation domain
    self -> GetProperty, XSIM=xSim, YSIM=ySim, ZSIM=zSim
    ix = value_locate(xSim, reform(coords[0,*])) > 0
    iz = value_locate(zSim, reform(coords[2,*])) > 0
    iy = n_elements(ySim) eq 1 ? replicate(0, count) : value_locate(ySim, reform(coords[1,*])) > 0
    
    ;Get the desired data product
    data = self -> GetData(name, _STRICT_EXTRA=extra)
    
    ;Extract the subregion
    case self.orientation of
        'XY':  data = data[ix, iy]
        'XZ':  data = data[ix, iz]
        'YZ':  data = data[iy, iz]
        'XYZ': data = data[ix, iy, iz]
    endcase

    ;Extract the coordinates
    if arg_present(coords) then begin
        case self.orientation of
            'XY':  coords = transpose([[xSim[ix]], [ySim[iy]]])
            'XZ':  coords = transpose([[xSim[ix]], [zSim[iz]]])
            'YZ':  coords = transpose([[ySim[ix]], [zSim[iz]]])
            'XYZ': coords = transpose([[xSim[ix]], [ySim[iy]], [zSim[iz]]])
        endcase
    endif

    ;Return the 1D cut
    return, data
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
pro MrSim2::GetInfo, $
;Derive Info
UNITS=units, $
ECOUNTFACTOR=eCountFactor, $
DTXWCI=dtxwci, $
INFO_DTXWCI=dtxwci_info, $

;From INFO file
B0=b0, $
B2_B1=B2_B1, $
COURANT=courant, $
DAMP=damp, $
DI=di, $
DTXWCE=dtxwce, $
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
    ;   - See MrSim_Which
    if arg_present(units)        then units = self.ion_scale ? 'di' : 'de'
    if arg_present(eCountFactor) then MrSim_Which, self.simname, TINDEX=self.time, ECOUNTFACTOR=eCountFactor
    if arg_present(dtxwci)       then MrSim_Which, self.simname, TINDEX=self.time, DTXWCI=dtxwci
    
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
    
    ; 3D-specific
    if self.dimension eq '3D' then begin
        ;Data has been downsampled (see comments in the INIT method). Calculate the modified
        ;size of each cell in de (the true size is stored in (*self.info).dx_de).
        if arg_present(dx_de) then dx_de = (*self.info).lx_de / (*self.info).nx
        if arg_present(dy_de) then dy_de = (*self.info).ly_de / (*self.info).ny
        if arg_present(dz_de) then dz_de = (*self.info).lz_de / (*self.info).nz
    endif else begin
        if arg_present(dx_de) then dx_de = (*self.info).dx_de
        if arg_present(dy_de) then dy_de = (*self.info).dy_de
        if arg_present(dz_de) then dz_de = (*self.info).dz_de
    endelse
    
    ;
    ; ASCII Info file
    ;   - Make sure it has been read.
    ;
    if has_tag(*self.info, 'MI_ME') eq 0 then return
    
    if arg_present(L_di)               then L_di       = (*self.info).L_di
    if arg_present(L_de)               then L_de       = (*self.info).L_de
    if arg_present(Ti_Te)              then Ti_Te      = (*self.info).Ti_Te
    if arg_present(Tbi_Ti)             then Tbi_Ti     = (*self.info).Tbi_Ti
    if arg_present(Tbe_Te)             then Tbe_Te     = (*self.info).Tbe_Te
    if arg_present(wpe_wce)            then wpe_wce    = (*self.info).wpe_wce
    if arg_present(mi_me)              then mi_me      = (*self.info).mi_me
    if arg_present(theta)              then theta      = (*self.info).theta
    if arg_present(taui)               then taui       = (*self.info).taui
    if arg_present(num_step)           then num_step   = (*self.info).num_step
    if arg_present(Lx_di)              then Lx_di      = (*self.info).Lx_di
    if arg_present(Ly_di)              then Ly_di      = (*self.info).Ly_di
    if arg_present(Lz_di)              then Lz_di      = (*self.info).Lz_di
    if arg_present(damp)               then damp       = (*self.info).damp
    if arg_present(courant)            then courant    = (*self.info).courant
    if arg_present(nproc)              then nproc      = (*self.info).nproc
    if arg_present(nppc)               then nppc       = (*self.info).nppc
    if arg_present(b0)                 then b0         = (*self.info).b0
    if arg_present(v_A)                then v_A        = (*self.info).v_A
    if arg_present(di)                 then di         = (*self.info).di
    if arg_present(N_e)                then N_e        = (*self.info).N_e
    if arg_present(Ne_sheet)           then Ne_sheet   = (*self.info).Ne_sheet
    if arg_present(Ne_back)            then Ne_back    = (*self.info).Ne_back
    if arg_present(N_total)            then N_total    = (*self.info).N_total
    if arg_present(dtxwpe)             then dtxwpe     = (*self.info).dtxwpe
    if arg_present(dtxwce)             then dtxwce     = (*self.info).dtxwce
    if arg_present(dtxwci_info)        then dtxwci     = (*self.info).dtxwci
    if arg_present(E_interval)         then E_interval = (*self.info).E_interval
    if arg_present(L_debye)            then L_debye    = (*self.info).L_debye
    if arg_present(dx_rhoi)            then dx_rhoi    = (*self.info).dx_rhoi
    if arg_present(dx_rhoe)            then dx_rhoe    = (*self.info).dx_rhoe
    if arg_present(dx_debye)           then dx_debye   = (*self.info).dx_debye
    if arg_present(n0)                 then n0         = (*self.info).n0
    if arg_present(vthi_c)             then vthi_c     = (*self.info).vthi_c
    if arg_present(vthe_c)             then vthe_c     = (*self.info).vthe_c
    if arg_present(vdri_c)             then vdri_c     = (*self.info).vdri_c
    if arg_present(vdre_c)             then vdre_c     = (*self.info).vdre_c
    if arg_present(nfac)               then nfac       = (*self.info).nfac
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
pro MrSim2::GetProperty, $
AXIS_LABELS = axis_labels, $
COORD_SYSTEM = coord_system, $
DIMENSION = dimension, $
DIRECTORY = directory, $
ION_SCALE = ion_scale, $
MVA_FRAME = mva_frame, $
ORIENTATION = orientation, $
SIMNAME = simname, $
SIMNUM = simnum, $
SYMMETRIC = symmetric, $
TIME = time, $

;SIMULATION DOMAIN
IXRANGE = ixrange, $
IYRANGE = iyrange, $
IZRANGE = izrange, $
XRANGE = xrange, $
XSIM = XSim, $
YRANGE = yrange, $
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
    if arg_present(dimension)    then dimension    = self.dimension
    if arg_present(directory)    then directory    = self.directory
    if arg_present(ion_scale)    then ion_scale    = self.ion_scale
    if arg_present(ixrange)      then ixrange      = getIndexRange(*self.XSim, self.xrange)
    if arg_present(iyrange)      then iyrange      = getIndexRange(*self.YSim, self.yrange)
    if arg_present(izrange)      then izrange      = getIndexRange(*self.ZSim, self.zrange)
    if arg_present(mva_frame)    then mva_frame    = self.mva_frame
    if arg_present(orientation)  then orientation  = self.orientation
    if arg_present(simname)      then simname      = self.simname
    if arg_present(simnum)       then simnum       = self.simnum
    if arg_present(symmetric)    then symmetric    = self.symmetric
    if arg_present(time)         then time         = self.time
    if arg_present(xrange)       then xrange       = self.xrange
    if arg_present(yrange)       then yrange       = self.yrange
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
;   The purpose of this program is to read data from a ".gda" file produced by 
;   one of Bill Daughton's simulation runs.
;
; :Params:
;       DATA_PRODUCT:           in, required, type=string, default=
;                               The name of the data product to be read. For a list of
;                                   available data product, call mr_readSIM without any
;                                   arguments.
;
; :Returns:
;       DATA:                   The requested data. If the data product does not exist,
;                                   then !Null will be returned.
;-
function MrSim2::GridLine, r0, r1, $
COUNT=count, $
COORDS=coords
	compile_opt strictarr

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		void = cgErrorMSG()
		return, !Null
	endif

	;Determine the size of a grid cell in DE
	self -> GetInfo, DX_DE=dx_de, DY_DE=dy_de, DZ_DE=dz_de, MI_ME=mi_me
	if self.ion_scale then begin
		dx_de /= mi_me
		dy_de /= mi_me
		dz_de /= mi_me
	endif

	;Number of grid cells between r0 and r1
	dr    = sqrt( total( (r1 - r0)^2 ) )
	count = fix(dr / dx_de, TYPE=3)

	;Coordinates of each point on the line.
	coords = MrLine3D(r0, r1, NPOINTS=count)

	;Find the cells
	cells = coords
	cells[0,*] = self -> Coord2Cell(coords[0,*], /X)
	cells[1,*] = self -> Coord2Cell(coords[1,*], /Y)
	cells[2,*] = self -> Coord2Cell(coords[2,*], /Z)

	return, cells
end


;+
;   Create a series of point along a line in 3D space, given to points that define the line.
;
; :Params:
;       R0:             in, required, type=fltarr(3)
;                       A point in cartesian space through which the plane should pass.
;       R1:             in, required, type=fltarr(3)
;                       A point in cartesian space through which the plane should pass.
;       R2:             in, required, type=fltarr(3)
;                       A point in cartesian space through which the plane should pass.
;
; :Keywords:
;       NX:             in, optional, type=integer, default=100
;                       In the frame of the plane, the number of points along the
;                           horizontal dimension.
;       NY:             in, optional, type=integer, default=100
;                       In the frame of the plane, the number of points along the
;                           vetical dimension.
;
; :Returns:
;       PLANE:          out, required, type=flaot(3\,`NX`*`NY`)
;                       Points on the plane containing `R0`, `R1`, and `R2`.
;-
function MrSim2::GridPlane, r0, r1, r2, $
NX=nx, $
NY=ny
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    ;Determine the size of a grid cell in DE
    self -> GetInfo, DX_DE=dx_de, DY_DE=dy_de, DZ_DE=dz_de
    
    ;Number of grid cells between r0 and r1
    dr = sqrt( total( (r1 - r0)^2 ) )
    n  = dr * dx_de
    
    ;Coordinates of each point on the line.
    coords = MrPlane3D(r0, r1, r2, NX=nx, NY=ny)
    
    ;Find the cells
    cells = coords
    cells[0,*] = self -> GetCell(coords[0,*], /X)
    cells[1,*] = self -> GetCell(coords[1,*], /Y)
    cells[2,*] = self -> GetCell(coords[2,*], /Z)

    return, cells
end


;+
;   The purpose of this program is to recreate the simulation domain so that a map
;   can be made between physical units and pixel locations.
;
; :Private:
;-
pro MrSim2::MakeSimDomain
	compile_opt strictarr
	on_error, 2

	;Get the simulations size in pixels and in electron skin depths.
	self -> GetInfo, LX_DE=xsize, LY_DE=ysize, LZ_DE=zsize, NX=nx, NY=ny, NZ=nz, MI_ME=mi_me

	;Convert to units of "di"
	if self.ion_scale then begin
		xsize /= sqrt(mi_me)
		ysize /= sqrt(mi_me)
		zsize /= sqrt(mi_me)
	endif

	;Set the simulation domain
	;   - See comments in ::ReadGDA_FilePos for effects of COORD_SYSTEM property.
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
;   The purpose of this method is to read the ASCII "info" file relating to Bill
;   Daughton's simulations.
;
; :Private:
;
; :Params:
;       FILENAME:           in, optional, type=string, default='electrons-y' [y-slice] + '.gda'
;                           Name of the "info" file to be read.
;-
function MrSim2::ReadElectrons, filename, $
DIST3D=dist3D, $
ENERGY=energy, $
FMAP_DIR=fmap_dir, $
MOMENTUM=momentum, $
VERBOSE=verbose, $
VX1_RANGE=vx1_range, $
VX2_RANGE=vx2_range, $
VX3_RANGE=vx3_range, $
X1_RANGE=x1_range, $
X2_RANGE=x2_range, $
X3_RANGE=x3_range
	compile_opt strictarr

	;catch errors
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		void = cgErrorMsg()
		return, !Null
	endif

	;Create defaults
	if n_elements(fMap_dir) eq 0 then MrSim_Which, self.simname, FMAP_DIR=fMap_dir
	if n_elements(filename) eq 0 || filename eq '' then begin
		;Get the y-grid-cell
		yslice = self -> Coord2Cell(self.yrange[0], /Y)
	
		;Get the file name
		MrSim_Which, self.simname, EFILE=filename, DIST3D=dist3D, $
		             TINDEX=self.time, YSLICE=yslice
	endif

	;Set defaults
	dist3D      =  keyword_set(dist3D)
	par_perp    =  keyword_set(par_perp)
	perp1_perp2 =  keyword_set(perp1_perp2)
	velocity    = ~keyword_set(momentum)

	x1r = n_elements(x1_range) gt 0 ? x1_range : self.xrange
	if dist3D then begin
		x2r = n_elements(x2_range) gt 0 ? x2_range : self.yrange
		x3r = n_elements(x3_range) gt 0 ? x3_range : self.zrange
	endif else begin
		x2r = n_elements(x2_range) gt 0 ? x2_range : self.zrange
	endelse

	;Convert ranges from COORD_SYS to SIMULATION coordinates
	;   - ORIENTATION does not matter because the full distribution is returned.
	case self.coord_system of
		'SIMULATION': ;Do nothing
		'MAGNETOPAUSE': begin
			;Convert spatial range
			x1r_temp  = x1r
			if dist3D then begin
				x1r = x3r
				x3r = temporary(x1r_temp)
			endif else begin
				x1r = x2r
				x2r = -reverse(temporary(x1r_temp))
			endelse
			if dist3D then x3r = -reverse(temporary(x1r))
		endcase
		'MAGNETOTAIL': x1r = -x1r
	endcase

	;Cannot be used together.
	if par_perp + perp1_perp2 gt 1 then $
		message, 'PAR_PERP and PERP1_PERP2 are mutually exclusive.'

;-------------------------------------------------------
; Read in Particle Data ////////////////////////////////
;-------------------------------------------------------
	if dist3D then begin
		data = MrSim_ReadParticles(filename, x1r, x2r, x3r, $
		                           FMAP_DIR  = fMap_dir, $
		                           UX1_RANGE = vx1_range, $
		                           UX2_RANGE = vx2_range, $
		                           UX3_RANGE = vx3_range, $
		                           VELOCITY  = velocity, $
		                           VERBOSE   = verbose)
		
		ix = [0,3]
		iz = [2,5]
	endif else begin
		data = MrSim_ReadParticles(filename, x1r, x2r, $
		                           FMAP_DIR  = fMap_dir, $
		                           UX1_RANGE = vx1_range, $
		                           UX2_RANGE = vx2_range, $
		                           UX3_RANGE = vx3_range, $
		                           VELOCITY  = velocity, $
		                           VERBOSE   = verbose)
		ix = [0,2]
		iz = [1,4]
	endelse

;-------------------------------------------------------
; Read in Particle Data ////////////////////////////////
;-------------------------------------------------------
	;Convert from SIMULATION coordinates to COORD_SYS coordinates.
	case self.coord_system of
		'SIMULATION': ;Do nothing
		'MAGNETOPAUSE': begin
			temp_data = data
			data[ix,*] = -temp_data[iz,*]
			data[iz,*] = (temporary(temp_data))[ix,*]
		endcase
		'MAGNETOTAIL': begin
			data[ix,*] = -data[ix,*]
		endcase
	endcase

;-------------------------------------------------------
; Parallel & Perpendicular /////////////////////////////
;-------------------------------------------------------
	if par_perp then begin
		;Indices at which momentum/velocity is stored
		if dist3D then iu = [2,3,4] else iu = [3,4,5]

		;B unit vector
		B_hat = self -> Unit_Vector('B', TIME=time, XRANGE=xrange, ZRANGE=zrange)
	
		;Get the parallel and perpendicular components
		v_mag_sqr = total(data[iu,*]^2, 1)
		v_par     = data[iu[0],*] * B_hat[0] + data[iu[1],*] * B_hat[1] + data[iu[2],*] * B_hat[3]
		v_perp    = sqrt(temporary(v_mag_sqr) - v_par^2)

		;Save the original data?
		if arg_present(original) then original = data

		;Fit into the data
		data = data[0:iu[1],*]
		data[iu[0],*] = temporary(v_par)
		data[iu[1],*] = temporary(v_perp)
	
		;Calculate the temperature of the distribution
		if arg_present(temperature) && velocity then begin
			nParticles     = n_elements(data[0,*])
			temperature    = fltarr(2)
			temperature[0] =       total(data[iu[0],*]^2) / nParticles
			temperature[1] = 0.5 * total(data[iu[1],*]^2) / nParticles
		endif
	endif

;-------------------------------------------------------
; Perp1 & Perp2 ////////////////////////////////////////
;-------------------------------------------------------
	if perp1_perp2 then begin
		;Indices at which momentum/velocity is stored
		if dist3D then iu = [2,3,4] else iu = [3,4,5]
	
		;B unit vector
		p2_hat = self -> Unit_Perp2(B_HAT=B_hat, P1_HAT=p1_hat, TIME=time, XRANGE=xrange, ZRANGE=zrange)
	
		;Get the parallel and perpendicular components
		v_par     = data[iu[0],*] *  B_hat[0] + data[iu[1],*] *  B_hat[1] + data[iu[2],*] *  B_hat[3]
		v_perp1   = data[iu[0],*] * p1_hat[0] + data[iu[1],*] * p1_hat[1] + data[iu[2],*] * p1_hat[3]
		v_perp2   = data[iu[0],*] * p2_hat[0] + data[iu[1],*] * p2_hat[1] + data[iu[2],*] * p2_hat[3]

		;Save the original data?
		if arg_present(original) then original = data

		;Fit into the data
		data[iu[0],*] = temporary(v_par)
		data[iu[1],*] = temporary(v_perp1)
		data[iu[2],*] = temporary(v_perp2)
	
		;Calculate the temperature of the distribution
		if arg_present(temperature) && velocity then begin
			nParticles     = n_elements(data[0,*])
			temperature    = fltarr(3)
			temperature[0] = total(data[iu[0],*]^2) / nParticles
			temperature[1] = total(data[iu[1],*]^2) / nParticles
			temperature[2] = total(data[iu[2],*]^2) / nParticles
		endif
	endif

	return, data
end



;+
;   Given coordinates in a recognized coordinate system, and assuming a data file
;   organized as a [nx, ny, nz] array in simulation coordinates, return the indices
;   into the data file that correspond to the given coordinates.
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
function MrSim2::ReadGDA_FilePos, coords, $
COORD_SYSTEM=coord_system, $
RANGES=ranges
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    ranges  = keyword_set(ranges)
    _system = n_elements(coord_system) eq 0 ? self.coord_system : strupcase(coord_system)

;-------------------------------------------------------
; Grid Points in Current System ////////////////////////
;-------------------------------------------------------
    ;Cell location in the current coordinate system.
    if ranges then begin
        ix = MrIndexRange(*self.xSim, self.xrange)
        iy = MrIndexRange(*self.ySim, self.yrange)
        iz = MrIndexRange(*self.zSim, self.zrange)
    endif else begin
        ix = value_locate(*self.xSim, coords[0,*]) > 0
        iy = value_locate(*self.ySim, coords[1,*]) > 0
        iz = value_locate(*self.zSim, coords[2,*]) > 0 
    endelse

;-------------------------------------------------------
; From Simulation //////////////////////////////////////
;-------------------------------------------------------
    ;We must convert from the current coordinate system back
    ;to simulation coordinates. Coordinate transformations from simulation to
    ;whatever are described below.
    case _system of
        'SIMULATION': ;Do nothing
    
        ; The entire transformation results in
        ;   x -> -z
        ;   y -> +y
        ;   z -> +x
        ; and takes place in in three stages. The first is here.
        ;
        ;
        ; Fields point along Z (N-S).
        ;   - Change order in which data is read from file (::ReadGDA_FilePos).
        ;
        ;           MSP                         MSP
        ; +z  *--------------|       ++x *--------------|
        ;     |              |  ==>      |              |
        ; -z  |--------------+        +x |--------------+
        ;     +x    MSH    ++x           -z     MSH    +z
        ;
        ; Rotate so fields point vertically (with MSH-MSP field pointing S-N).
        ;   - Transpose data (::CoordTransform)
        ;
        ;                            ++z +--------|
        ;                                |        |
        ;           MSP                  |        |
        ; +x  *--------------|         M |        | M
        ;     |              |  ===>   S |        | S
        ; -x  |--------------+         H |        | P
        ;     +z    MSH    ++z           |        |
        ;                                |        |
        ;                             +z |--------*
        ;                               -x       +x
        ;
        ; Make X point toward the sun (to the left)
        ;   - Reverse x-axis (::MakeSimDomain)
        ;
        ;  ++z +--------|        ++z +--------|
        ;      |        |            |        |
        ;      |        |            |        |
        ;    M |        | M        M |        | M
        ;    S |        | S  ===>  S |        | S
        ;    H |        | P        H |        | P
        ;      |        |            |        |
        ;      |        |            |        |
        ;   +z |--------*         +z |--------*
        ;     -x       +x           +x       -x
        ;
        'MAGNETOPAUSE': begin
            ;Interchange x and z
            iTemp = ix
            ix    = iz
            iz    = temporary(iTemp)
        endcase
    
        ;
        ; The entire transformation results in
        ;   x -> -x
        ; and takes place in one stage, none of which occur here.
        ;
        ; Reverse the direction of vectors pointing in the x-direction
        ;   - ::CoordTransform
        ;
        ;          North                       North
        ; +z  *--------------|        +z *--------------|
        ;     |              |  ==>      |              |
        ; -z  |--------------+        -z |--------------+
        ;     +x   South    ++x         -x     South   --x
        ;
        'MAGNETOTAIL': ;Do nothing
    
        else: message, 'Coordinate system not recognized: "' + _system + '".'
    endcase
    
;-------------------------------------------------------
; From Simulation //////////////////////////////////////
;-------------------------------------------------------
    ;Sort ranges
    if ranges then begin
        if ix[0] gt ix[1] then ix = ix[[1,0]]
        if iy[0] gt iy[1] then iy = iy[[1,0]]
        if iz[0] gt iz[1] then iz = iz[[1,0]]
    endif
    
    return, transpose([[ix], [iy], [iz]])
end


;+
;   Transform data from simulation coordinates to another coordinate system.
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
pro MrSim2::ReadGDA_TransformData, name, data, $
COORD_SYSTEM=coord_system, $
ORIENTATION=orientation
	compile_opt strictarr

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /cancel
		void = cgErrorMSG()
		return
	endif

	_orient = n_elements(orientation)  eq 0 ? self.orientation  : strupcase(orientation)
	_system = n_elements(coord_system) eq 0 ? self.coord_system : strupcase(coord_system)

;-------------------------------------------------------
; Change from SIMULATION Coordinates? //////////////////
;-------------------------------------------------------
	;See ::ReadGDA_FilePos for detailed comments about
	;transformations.

	case strupcase(_system) of
		'SIMULATION': ;Do nothing

		'MAGNETOPAUSE': begin
			case _orient of
				'XY': ;Do nothing
				'XZ': data = transpose(data)
				else: message, 'Orienatation "' + _orient + '" not recognized.'
			endcase

			;Now reverse the z-axes
			;   Single negate Bz, Uiz, Pe-xz, Pe-yz, etc.
			case 1 of
				stregex(name, '[A-Z][a-z]*z$',   /BOOLEAN): data = -data
				stregex(name, '-(z[yx]|[xy]z)$', /BOOLEAN): data = -data
				else: ;Do nothing
			endcase
		endcase

		'MAGNETOTAIL': begin
			case _orient of
				'XY': ;Do nothing
				'XZ': ;Do nothing
				else: message, 'Orienatation "' + _orient + '" not recognized.'
			endcase

			;Now reverse the x-axes
			;   Single negate Bx, Uix, Pe-xz, Pe-xy, etc.
			case 1 of
				stregex(name, '[A-Z][a-z]*(x)$', /BOOLEAN): data = -data
				stregex(name, '-(y|z)[x]$',      /BOOLEAN): data = -data
				stregex(name, '-[x](y|z)$',      /BOOLEAN): data = -data
				else: ;Do nothing
			endcase
		endcase
	endcase
end


;+
;   The purpose of this method is to read data. If the data file does not exist, an
;   informational message will be printed and no data will be read.
;
; :Params:
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;
;-
function MrSim2::ReadGDA_Cells, cells, name, tIndex, $
COORD_SYSTEM=coord_system, $
COORDINATES=coordinates, $
DIRECTORY=directory, $
NSMOOTH=nSmooth, $
ORIENTATION=orientation
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMsg()
        return, !Null
    endif

;-------------------------------------------------------
; Defaults /////////////////////////////////////////////
;-------------------------------------------------------
    if n_elements(directory)   eq 0 then directory   = self.directory
    if n_elements(nSmooth)     eq 0 then nSmooth     = self.nSmooth
    if n_elements(orientation) eq 0 then orientation = self.orientation
    if n_elements(tindex)      eq 0 then tindex      = self.time
    coord_system = n_elements(coord_system) eq 0 ? self.coord_system : strupcase(coord_system)

    ;Create file name and test for the file
    ;   - If file does not exist, create dialog to pick the file.
    ;   - If dialog is cancelled, return nothing.
    filename = filepath(name + '_' + string(tindex, FORMAT='(i0)') + '.gda', ROOT_DIR=directory)
    if file_test(filename) eq 0 then begin
        if filename ne '' then message, 'Choose a file. Given file not found: "' + filename + '".', /INFORMATIONAL
        filename = cgPickFile(TITLE='Choose .gda Field or Moment File.', /READ)
        if filename eq '' then return, !Null
    endif

;-------------------------------------------------------
; Read the Data ////////////////////////////////////////
;-------------------------------------------------------
    ;File size
    self -> GetInfo, NX=nx, NY=ny, NZ=nz

    ;Allocate memory
    _cells = ulong64(cells)
    nCells = n_elements(cells[0,*])
    data   = fltarr(nCells)
    record = 0.0

    ;Open the file
    openr, lun, filename, /GET_LUN

    ;Step through each record.
    for i = 0ULL, nCells-1 do begin
        ;Get the cell locations.
        ix = _cells[0,i]
        iy = _cells[1,i]
        iz = _cells[2,i]
        
        ;Jump to the desired cell
        offset = 4ULL * (ix + nx*iy + nx*ny*iz)
        
        ;Read and store the record.
        point_lun, lun, offset
        readu, lun, record
        data[i] = record
    endfor
    
    ;Free the file
    free_lun, lun

;-------------------------------------------------------
; Cleanup //////////////////////////////////////////////
;-------------------------------------------------------
    ;Return the coordinates at which the grid-cells are located
    if arg_present(coordinates) then begin
        coordinates = cells
        coordinates[0,*] = (*self.xSim)[cells[0,*]]
        coordinates[1,*] = (*self.ySim)[cells[1,*]]
        coordinates[2,*] = (*self.zSim)[cells[2,*]]
    endif

    ;Smooth data 
    if nSmooth gt 0 then data = smooth(data, nSmooth, /EDGE_TRUNCATE)  

    ;Store the data so that it does not have to be read again.
    return, data
end


;+
;   Read a 2D slice of data in the XY plane. Data is assumed to be stored in as
;   floating point values ordered in the same as an [nx, ny, xz] array.
;
; :Private:
;
; :Params:
;       FILE:               in, required, type=string
;                           Name of the file to be read.
;       IXRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the x-dimension.
;       IYRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the y-dimension.
;       IZ                  in, required, type=intarr(2)
;                           Grid cell along the z-direction at which to read the 2D slice.
;
; :Keywords:
;       NX:                 in, optional, type=integer
;                           Number of total grid cells along the x-direction.
;       NY:                 in, optional, type=integer
;                           Number of total grid cells along the y-direction.
;       NZ:                 in, optional, type=integer
;                           Number of total grid cells along the z-direction.
;
; :Returns:
;       DATA:               Data extracted from the file.
;-
function MrSim2::ReadGDA_XY, file, ixrange, iyrange, iz, $
NX=nx, $
NY=ny, $
NZ=nz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get/Set defaults
    self -> GetInfo, NX=_nx, NY=_ny, NZ=_nz
    if n_elements(nx) eq 0 then nx = _nx
    if n_elements(ny) eq 0 then ny = _ny
    if n_elements(nz) eq 0 then nz = _nz

    ;Minimum and maximum grid cells to read
    ixMin = min(ixrange, MAX=ixMax)
    iyMin = min(iyrange, MAX=iyMax)
    
    ;Open File and read data
    openr, lun, file, /GET_LUN
    
    ;Declare array for data
    data = fltarr(ixMax-ixMin+1, iyMax-iyMin+1)
    temp = fltarr(ixMax-ixMin+1)

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------
    
    ;How to read:
    ;   - Z is fixed
    ;   - Step through each row (y)
    ;   - Read all columns (x) in each row
    ix      = ulong64(ixMin)
    zOffset = iz*nx*ny*4ULL
    
    ;Read the data
    for iy = ulong64(iyMin), iyMax do begin
        point_lun, lun, zOffset + 4ULL*(ix + iy*nx)
        readu, lun, temp
        data[*,iy] = temp
    endfor
    
    ;Free the file
    free_lun, lun
    
    return, data
end


;+
;   Read a 2D slice of data in the XZ plane. Data is assumed to be stored in as
;   floating point values ordered in the same as an [nx, ny, xz] array.
;
; :Private:
;
; :Params:
;       FILE:               in, required, type=string
;                           Name of the file to be read.
;       IXRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the x-dimension.
;       IZRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the z-dimension.
;       IY                  in, required, type=intarr(2)
;                           Grid cell along the y-direction at which to read the 2D slice.
;
; :Keywords:
;       NX:                 in, optional, type=integer
;                           Number of total grid cells along the x-direction.
;       NY:                 in, optional, type=integer
;                           Number of total grid cells along the y-direction.
;       NZ:                 in, optional, type=integer
;                           Number of total grid cells along the z-direction.
;
; :Returns:
;       DATA:               Data extracted from the file.
;-
function MrSim2::ReadGDA_XZ, file, ixrange, izrange, iy, $
NX=nx, $
NY=ny, $
NZ=nz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get/Set defaults
    self -> GetInfo, NX=_nx, NY=_ny, NZ=_nz
    if n_elements(nx) eq 0 then nx = _nx
    if n_elements(ny) eq 0 then ny = _ny
    if n_elements(nz) eq 0 then nz = _nz

    ;Minimum and maximum grid cells to read
    ixMin = min(ixrange, MAX=ixMax)
    izMin = min(izrange, MAX=izMax)

    ;Open File and read data
    openr, lun, file, /GET_LUN
    
    ;Declare output array for data
    data = fltarr(ixMax-ixMin+1, izMax-izMin+1)

    ;Read this many elements from a row.
    temp = fltarr(ixMax-ixMin+1)

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------

    ;Offset in time
    ;   - Each time slice consists of an [nx,ny,nz] data array.
    ;   - 2D files contain many time slices and one y-slice.
    ;   - 3D files contain a single time, but many y-slices.
    if self.dimension eq '2D' $
        then tOffset = (ulong64(nx)*ulong64(nz) + 2ULL) * self.time * 4ULL $
        else tOffset = 0ULL

    ;How to read:
    ;   - Row (y) is fixed is fixed
    ;   - Jump from z to z
    ;   - Read all columns (x) at each z
    ix = ulong64(ixMin)
    for iz = ulong64(izMin), izMax do begin
        ;Offset to first record.
        ;   ix       = column offset
        ;   iy*ny    = row offset
        ;   iz*nx*ny = 3rd dimension offset
        ;   4        = Bytes per element
        offset = 4ULL * (ix + iy*nx + iz*nx*ny) + tOffset
        
        ;Read and save the records.
        point_lun, lun, offset
        readu, lun, temp
        data[*,iz-izMin] = temp
    endfor
    
    ;Free the file
    free_lun, lun

    return, data
end


;+
;   Read a 2D slice of data in the YZ plane. Data is assumed to be stored in as
;   floating point values ordered in the same as an [nx, ny, xz] array.
;
; :Private:
;
; :Params:
;       FILE:               in, required, type=string
;                           Name of the file to be read.
;       IXRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the x-dimension.
;       IYRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the y-dimension.
;       IZ                  in, required, type=intarr(2)
;                           Grid cell along the z-direction at which to read the 2D slice.
;
; :Keywords:
;       NX:                 in, optional, type=integer
;                           Number of total grid cells along the x-direction.
;       NY:                 in, optional, type=integer
;                           Number of total grid cells along the y-direction.
;       NZ:                 in, optional, type=integer
;                           Number of total grid cells along the z-direction.
;
; :Returns:
;       DATA:               Data extracted from the file.
;-
function MrSim2::ReadGDA_YZ, file, iyrange, izrange, ix, $
NX=nx, $
NY=ny, $
NZ=nz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get/Set defaults
    self -> GetInfo, NX=_nx, NY=_ny, NZ=_nz
    if n_elements(nx) eq 0 then nx = _nx
    if n_elements(ny) eq 0 then ny = _ny
    if n_elements(nz) eq 0 then nz = _nz

    ;Minimum and maximum grid cells to read
    iyMin = min(iyrange, MAX=iyMax)
    izMin = min(izrange, MAX=izMax)

    ;Open File and read data
    openr, lun, file, /GET_LUN
    
    ;Declare array for data
    data = fltarr(iyMax-iyMin+1, izMax-izMin+1)
    temp = fltarr(1)

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------
    
    ;How to read
    ;   - X is fixed
    ;   - Jump to desired row (y)
    ;   - Jump to desired depth (z)
    ;   - y and z are not adjacent, so each element must be read
    for iz = ulong64(izMin), izMax do begin
        for iy = ulong64(iyMin), iyMax do begin
            ;Offset to first record.
            ;   ix       = column offset
            ;   iy*ny    = row offset
            ;   iz*nx*ny = 3rd dimension offset
            ;   4        = Bytes per element
            offset = 4ULL * (ix + iy*nx + iz*nx*ny)
        
            ;Read and save the records.
            point_lun, lun, offset
            readu, lun, temp
            data[iy-iyMin, iz-izMin] = temp
        endfor
    endfor
    
    ;Free the file
    free_lun, lun
    
    return, data
end


;+
;   Read a 3D slice of data. Data is assumed to be stored in as floating point values
;   ordered in the same as an [nx, ny, xz] array.
;
; :Private:
;
; :Params:
;       FILE:               in, required, type=string
;                           Name of the file to be read.
;       IXRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the x-dimension.
;       IYRANGE             in, required, type=intarr(2)
;                           Range of grid cells to be read along the y-dimension.
;       IZ                  in, required, type=intarr(2)
;                           Grid cell along the z-direction at which to read the 2D slice.
;
; :Keywords:
;       NX:                 in, optional, type=integer
;                           Number of total grid cells along the x-direction.
;       NY:                 in, optional, type=integer
;                           Number of total grid cells along the y-direction.
;       NZ:                 in, optional, type=integer
;                           Number of total grid cells along the z-direction.
;
; :Returns:
;       DATA:               Data extracted from the file.
;-
function MrSim2::ReadGDA_XYZ, file, ixrange, iyrange, izrange, $
NX=nx, $
NY=ny, $
NZ=nz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get/Set defaults
    self -> GetInfo, NX=_nx, NY=_ny, NZ=_nz
    if n_elements(nx) eq 0 then nx = _nx
    if n_elements(ny) eq 0 then ny = _ny
    if n_elements(nz) eq 0 then nz = _nz

    ;Minimum and maximum grid cells to read
    ixMin = min(ixrange, MAX=ixMax)
    iyMin = min(iyrange, MAX=iyMax)
    izMin = min(izrange, MAX=izMax)
    
    ;Open File and read data
    openr, lun, file, /GET_LUN
    
    ;Declare array for data
    data = fltarr(ixMax-ixMin+1, iyMax-iyMin+1, izMax-izMin+1)
    temp = fltarr(ixMax-ixMin+1)

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------
    
    ;How to read:
    ;   - Advance to the first desired column
    ;   - Read all desired columns(x)
    ;   - Advance row (y) and repeat.
    ;   - When a full xy-plane is read, advance depth (z), repeat.
    ix = ulong64(ixMin)
    
    ;Read the data
    for iz = ulong64(izMin), izMax do begin
        for iy = ulong64(iyMin), iyMax do begin
            point_lun, lun, 4ULL*(ix + iy*nx + iz*nx*ny)
            readu, lun, temp
            data[*, iy, iz] = temp
        endfor
    endfor
    
    ;Free the file
    free_lun, lun
    
    return, data
end


;+
;   The purpose of this method is to read data. If the data file does not exist, an
;   informational message will be printed and no data will be read.
;
; :Params:
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;
;-
function MrSim2::ReadGDA, name, tIndex, $
COORD_SYSTEM=coord_system, $
DIRECTORY=directory, $
NSMOOTH=nSmooth, $
ORIENTATION=orientation, $
SIM_OBJECT=oSim, $
XRANGE=xrange, $
YRANGE=yrange, $
ZRANGE=zrange
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMsg()
        return, !Null
    endif

;-------------------------------------------------------
; Defaults /////////////////////////////////////////////
;-------------------------------------------------------
    ;Change Pe_xx to Pe-xx
    _name = strpos(name, '_') eq -1 ? name : strjoin(strsplit(name, '_', /EXTRACT), '-')

    ;Defaults
    if n_elements(directory)   eq 0 then directory   = self.directory
    if n_elements(nSmooth)     eq 0 then nSmooth     = self.nSmooth
    if n_elements(orientation) eq 0 then orientation = self.orientation
    if n_elements(tindex)      eq 0 then tindex      = self.time
    if n_elements(xrange)      eq 0 then xrange      = self.xrange
    if n_elements(yrange)      eq 0 then yrange      = self.yrange
    if n_elements(zrange)      eq 0 then zrange      = self.zrange
    coord_system = n_elements(coord_system) eq 0 ? self.coord_system : strupcase(coord_system)

    ;Create file name.
    if self.dimension eq '2D' $
        then file = filepath(_name + '.gda', ROOT_DIR=directory) $
        else file = filepath(_name + '_' + string(tindex, FORMAT='(i0)') + '.gda', ROOT_DIR=directory)
        
    ;Ensure file exists.
    if file_test(file) eq 0 then begin
        if file ne '' then message, 'Choose a file. Given file not found: "' + file + '".', /INFORMATIONAL
        file = cgPickFile(TITLE='Choose .gda Field or Moment File.', /READ)
        if file eq '' then return, !Null
    endif

;-------------------------------------------------------
; Read the Data ////////////////////////////////////////
;-------------------------------------------------------

    ;Data in the gda files are stored in simulation coordinates
    ;   - Translate from COORD_SYSTEM to SIMULATION cell locations.
    cells = self -> ReadGDA_FilePos(/RANGES, COORD_SYSTEM=coord_system)

    case orientation of
        'XY':  data = self -> ReadGDA_XY( file, cells[0,*], cells[1,*], cells[2,0])
        'XZ':  data = self -> ReadGDA_XZ( file, cells[0,*], cells[2,*], cells[1,0])
        'YZ':  data = self -> ReadGDA_YZ( file, cells[1,*], cells[2,*], cells[0,0])
        'XYZ': data = self -> ReadGDA_XYZ(file, cells[1,*], cells[2,*], cells[3,*])
        else: message, 'Orientation "' + orientation + '" not recognized.'
    endcase

    ;Manipulate data from SIMULATION to COORD_SYSTEM coordinates.
    self -> ReadGDA_TransformData, _name, data

;-------------------------------------------------------
; Cleanup //////////////////////////////////////////////
;-------------------------------------------------------
    ;Smooth data 
    if nSmooth gt 0 then data = smooth(data, nSmooth, /EDGE_TRUNCATE)  

    ;Store the data so that it does not have to be read again.
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
pro MrSim2::ReadInfo_Ascii, filename, $
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
pro MrSim2::ReadInfo_Binary, info_file, status
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
             dx_de: xmax / nx, $
             dy_de: ymax / ny, $
             dz_de: zmax / nz, $
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
pro MrSim2::SetSim, sim_id, $
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
pro MrSim2::SetScale, $
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
pro MrSim2::SetProperty, $
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
                    'SIMULATION': begin
                        xrange = self.xrange
                        zrange = self.zrange
                    endcase
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
                    'MAGNETOPAUSE': begin
                        xrange = self.xrange
                        zrange = self.zrange
                    endcase
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
                    'MAGNETOTAIL': begin
                        xrange = self.xrange
                        zrange = self.zrange
                    endcase
                endcase
            endcase
            
            else: message, 'Coordinate system not recognized:"' + coord_system + '".'
        endcase
        
        ;Set object properties
        self.xrange       = xrange
        self.zrange       = zrange
        self.coord_system = _coord_system
        
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
function MrSim2::TxWci2tIndex, txwci
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
function MrSim2::tIndex2TxWci, tIndex
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
;   This definition statement for the MrSim2 class.
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
;-
pro MrSim2__DEFINE, class
    compile_opt strictarr
    
    class = { MrSim2, $
              inherits IDL_Object, $
              inherits MrSim2_Data, $
             
              ;Object Properties
              axis_labels:  ['', '', ''], $
              dimension:    '', $
              directory:    '', $
              info:         ptr_new(), $
              ion_scale:    0B, $
              little:       0B, $
              mva_frame:    0B, $
              nSmooth:      0L, $
              orientation:  '', $
              simnum:       0, $
              simname:      '', $
              symmetric:    0B, $
              time:         0L, $
              coord_system: '', $
              
              ;Sim Domain
              XRANGE: [0D, 0D], $
              XSim:   ptr_new(), $
              YRANGE: [0D, 0D], $
              YSim:   ptr_new(), $
              ZRANGE: [0D, 0D], $
              ZSim:   ptr_new() $
            }
end