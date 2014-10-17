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
;-
;*****************************************************************************************
;+
;   This method initializes the MrSim2 class.
;
; :Params:
;       TIME:           in, optional, type=long, default=0
;                       Time index at which simulation data will be read.
;       YSLICE:         in, optional, type=long, default=0
;                       The Y-index at which to obtain data from the XZ plane. Ignored
;                           in 2D simulations.
;
; :Keywords:
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
AXIS_LABELS = axis_labels, $
BINARY = binary, $
COORD_SYSTEM = coord_system, $
DIRECTORY = directory, $
INFO_FILE = info_file, $
ION_SCALE = ion_scale, $
MVA_FRAME = mva_frame, $
NSMOOTH = nsmooth, $
ORIENTATION = orientation, $
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
;Read the Info File ////////////////////////////////////
;-------------------------------------------------------
    ;Set the data directory
    binary = keyword_set(binary)
    if n_elements(directory) eq 0 then void = cgRootName(DIRECTORY=directory)
    self.directory = directory
    
    ;Guess where the info file is.
    if n_elements(info_file) eq 0 then begin
        if binary then begin
            info_file = filepath('info', ROOT_DIR=directory)
        endif else begin
            info_dir  = strsplit(directory, path_sep(), /EXTRACT, COUNT=nDir)
            info_dir  = strjoin(info_dir[0:nDir-2], path_sep())
            info_file = filepath('info', ROOT_DIR=info_dir)
        endelse
    endif

    ;Make sure the file exists
    self.info = ptr_new(/ALLOCATE_HEAP)
    if file_test(info_file) eq 0 then $
        message, 'Cannot find ASCII info file: "' + info_file + '"'
        
    if binary $
        then self -> ReadBinaryInfo, info_file $
        else self -> ReadAsciiInfo, info_file

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
    if n_elements(yrange)  eq 0 then yrange  = [0, ysize]
    if n_elements(yslice)  eq 0 then yslice  = 0L
    if n_elements(zrange)  eq 0 then zrange  = [-zsize/2.0, zsize/2.0]
    
    if max(coord_system eq ['SIMULATION', 'MAGNETOPAUSE', 'MAGNETOTAIL']) eq 0 then $
        message, 'Coordinate system "' + coord_system + '" not recognized.'
    if max(_orientation eq ['XY', 'XZ']) eq 0 then $
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
    self.directory    = directory
    self.coord_system = coord_system
    self.ion_scale    = ion_scale
    self.orientation  = _orientation
    self.mva_frame    = mva_frame
    self.nsmooth      = nsmooth
    self.time         = time
    self.xrange       = xrange
    self.yrange       = yrange
    self.yslice       = yslice
    self.zrange       = zrange

    ;Set the domain coordinates
    self -> SetScale, ION_SCALE=ion_scale
    self -> MakeSimDomain

;-------------------------------------------------------
; Superclass ///////////////////////////////////////////
;-------------------------------------------------------

    ;Superclass
    success = MrSim2_Data::Init()
    if success eq 0 then return, 0

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
    ptr_free, self.ZSim
    
    ;Data products
    self -> MrSim2_Data::Cleanup
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
pro MrSim2::Clear_Data, data_product
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
        data_product = ['Ay', 'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'ne', 'ni', $
                        'Pe_xx', 'Pe_xy', 'Pe_xz', 'Pe_yx', 'Pe_yy', 'Pe_yz', $
                        'Pe_zx', 'Pe_zy', 'Pe_zz', 'Pi_xx', 'Pi_xy', 'Pi_xz', $
                        'Pi_yx', 'Pi_yy', 'Pi_yz', 'Pi_zx', 'Pi_zy', 'Pi_zz', $
                        'Uex', 'Uey', 'Uez', 'Uix', 'Uiy', 'Uiz']
    endif
    
    ;Step through all of the data products given and free the data
    foreach name, data_product do begin
        case name of
            'Ay':    *self.Ay    = !Null
            'Bx':    *self.Bx    = !Null
            'By':    *self.By    = !Null
            'Bz':    *self.Bz    = !Null
            'Ex':    *self.Ex    = !Null
            'Ey':    *self.Ey    = !Null
            'Ez':    *self.Ez    = !Null
            'ne':    *self.n_e   = !Null
            'ni':    *self.n_i   = !Null
            'Pe-xx': *self.Pe_xx = !Null
            'Pe-xy': *self.Pe_xy = !Null
            'Pe-xz': *self.Pe_xz = !Null
            'Pe-yx': *self.Pe_yx = !Null
            'Pe-yy': *self.Pe_yy = !Null
            'Pe-yz': *self.Pe_yz = !Null
            'Pe-zx': *self.Pe_zx = !Null
            'Pe-zy': *self.Pe_zy = !Null
            'Pe-zz': *self.Pe_zz = !Null
            'Pi-xx': *self.Pi_xx = !Null
            'Pi-xy': *self.Pi_xy = !Null
            'Pi-xz': *self.Pi_xz = !Null
            'Pi-yx': *self.Pi_yx = !Null
            'Pi-yy': *self.Pi_yy = !Null
            'Pi-yz': *self.Pi_yz = !Null
            'Pi-zx': *self.Pi_zx = !Null
            'Pi-zy': *self.Pi_zy = !Null
            'Pi-zz': *self.Pi_zz = !Null
            'Uex':   *self.Uex   = !Null
            'Uey':   *self.Uey   = !Null
            'Uez':   *self.Uez   = !Null
            'Uix':   *self.Uix   = !Null
            'Uiy':   *self.Uiy   = !Null
            'Uiz':   *self.Uiz   = !Null
            else: ;Do nothing
        endcase
    endforeach
end


;+
;   The purpose of this method is to retrieve data from the info file.
;
; :Keywords:
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
    if arg_present(units) then if self.ion_scale then units = 'di' else units = 'de'
    
    ;Density asymmetry: n_MSP / n_MSH
    if arg_present(n1_n2) then begin
        message, 'A hard-coded value of n2/n1 = 0.125 has been set.', /INFORMATIONAL
        ni_ne = 0.125
    endif
    
    ;Magnetic field asymmetry B_MSP / B_MSH
    if arg_present(B1_B2) then begin
        message, 'B2/B1 = 1.37 has been hard-coded.', /INFORMATIONAL
        B2_B1 = 1.27
    endif
    
    ;Data from the INFO file.
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
    if arg_present(Lx_de)      then Lx_de      = (*self.info).Lx_de
    if arg_present(Ly_de)      then Ly_de      = (*self.info).Ly_de
    if arg_present(Lz_de)      then Lz_de      = (*self.info).Lz_de
    if arg_present(Lx_di)      then Lx_di      = (*self.info).Lx_di
    if arg_present(Ly_di)      then Ly_di      = (*self.info).Ly_di
    if arg_present(Lz_di)      then Lz_di      = (*self.info).Lz_di
    if arg_present(nx)         then nx         = (*self.info).nx
    if arg_present(ny)         then ny         = (*self.info).ny
    if arg_present(nz)         then nz         = (*self.info).nz
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
pro MrSim2::GetProperty, $
AXIS_LABELS = axis_labels, $
COORD_SYSTEM = coord_system, $
DIRECTORY = directory, $
ION_SCALE = ion_scale, $
MVA_FRAME = mva_frame, $
ORIENTATION = orientation, $
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
    if arg_present(ixrange)      then ixrange      = getIndexRange(*self.XSim, self.xrange)
    if arg_present(iyrange)      then iyrange      = getIndexRange(*self.YSim, self.yrange)
    if arg_present(izrange)      then izrange      = getIndexRange(*self.ZSim, self.zrange)
    if arg_present(mva_frame)    then mva_frame    = self.mva_frame
    if arg_present(orientation)  then orientation  = self.orientation
    if arg_present(time)         then time         = self.time
    if arg_present(xrange)       then xrange       = self.xrange
    if arg_present(yrange)       then yrange       = self.yrange
    if arg_present(yslice)       then yslice       = (*self.YSim)[self.yslice]
    if arg_present(zrange)       then zrange       = self.zrange

    ;Simulation Domain
    if arg_present(XSim) && n_elements(*self.XSim) gt 0 then begin
        ix = getIndexRange(*self.XSim, self.xrange)
        XSim = (*self.XSim)[ix[0]:ix[1]]
    endif
    
    if arg_present(YSim) && n_elements(*self.YSim) gt 0 then begin
        iy = getIndexRange(*self.YSim, self.yrange)
        YSim = (*self.YSim)[iy[0]:iy[1]]
    endif

    if arg_present(ZSim) && n_elements(*self.ZSim) gt 0 then begin
        iz = getIndexRange(*self.ZSim, self.zrange)
        ZSim = (*self.ZSim)[iz[0]:iz[1]]
    endif
end


;+
;   The purpose of this program is to recreate the simulation domain so that a map
;   can be made between physical units and pixel locations.
;
; :Private:
;-
pro MrSim2::MakeSimDomain
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
    case self.coord_system of
        'SIMULATION': begin
            *self.XSim = linspace(0, xsize, nx)
            *self.YSim = linspace(0, ysize, ny)
            *self.ZSim = linspace(0, zsize, nz) - zsize/2.0
        endcase
        
        'MAGNETOPAUSE': begin
            *self.ZSim =   linspace(0, xsize, nx)
            *self.YSim =   linspace(0, ysize, ny)
            *self.XSim = -(linspace(0, zsize, nz) - zsize/2.0)
        endcase
        
        'MAGNETOTAIL': begin
            *self.XSim = -linspace(0, xsize, nx)
            *self.YSim =  linspace(0, ysize, ny)
            *self.ZSim =  linspace(0, zsize, nz) - zsize/2.0
        endcase
        
        else: ;Nothing
    endcase
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
;-
pro MrSim2::ReadAsciiInfo, filename
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        
        ;close the file and free its logical unit number
        if n_elements(lun) ne 0 then close, lun
        
        void = cgErrorMsg()
        return
    endif
    
    ;Read the simulation info file
    info = MrSim2_ReadInfoAscii(filename)

    ;Store the information
    *self.info = info
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
pro MrSim2::ReadBinaryInfo, info_file, status
    compile_opt strictarr

    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif

    status = 0
    if n_elements(info_file) eq 0 then info_file = self.directory + path_sep() + 'info'

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
    
    ;Convert to "di" units
    if ion_scale then begin
        self -> MakeSimDomain
        self.xrange /= mi_me
        self.yrange /= mi_me
        self.zrange /= mi_me
        
    ;Convert to "de" units
    endif else begin
        self -> MakeSimDomain
        self.xrange *= mi_me
        self.yrange *= mi_me
        self.zrange *= mi_me
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

    ;The following properties require the data to be cleared
    nTime   = n_elements(time)
    nXRange = n_elements(xrange)
    nYRange = n_elements(yrange)
    nYSlice = n_elements(yslice)
    nZRange = n_elements(zrange)
    nOrient = n_elements(orientation)
    nCoord  = n_elements(coord_system)
    nMVA    = n_elements(mva_frame)
    if nTime + nXRange = nYRange + nYSlice + nZRange + nOrient + nCoord + nMVA gt 0 $
        then self -> Clear_Data
        
    if nTime   gt 0 then self.time   = time
    if nXRange gt 0 then self.xrange = xrange
    if nYRange gt 0 then self.yrange = yrange
    if nYSlice gt 0 then self.yslice = yslice
    if NZRange gt 0 then self.zrange = zrange
    
    ;ORIENTATION
    if nOrient gt 0 then begin
        _orientation = strupcase(orientation)
        planes = ['XY', 'XZ']
        if max(_orientation eq planes) $
            then self.orientation = orientation $
            else message, 'Orientation "' + orientation + '" not recognized.', /INFORMATIONAL
    endif
    
    ;COORDINATE SYSTEM
    if nCoord gt 0 then begin
        _coord_system = strupcase(coord_system)
        systems = ['SIMULATION', 'MAGNETPAUSE', 'MAGNETOTAIL']
        if max(_coord_system eq systems) $
            then self.coord_system = _coord_system $
            else message, 'Coordinate system "' + coord_system + '" not recognized.', /INFORMATIONAL
        
        ;Remake the simulation domain
        self -> MakeSimDomain
    endif
    
    ;MVA_FRAME
    if nMVA gt 0 then begin
        self.mva_frame = keyword_set(mva_frame)
        
        ;Pick the axis labels
        case self.coord_system of
            'SIMULATION':   self.axis_labels = self.mva_frame ? ['x', 'y', 'z'] : ['x', 'y', 'z']
            'MAGNETOPAUSE': self.axis_labels = self.mva_frame ? ['L', 'M', 'N'] : ['z', 'y', 'x']
            'MAGNETOTAIL':  self.axis_labels = self.mva_frame ? ['N', 'M', 'L'] : ['x', 'y', 'z']
        endcase
    endif
    
    if clear_flag then self -> Clear_Data
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
function MrSim2::LineCuts, data_product, locations, pos, $
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
function MrSim2::XYCuts, data_product, locations, pos, $
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
        iRange = getIndexRange(YSim, cut_yrange)
        pos = YSim[iRange[0]:iRange[1]]

        ;X-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(XSim, locations)
        void = ismember(XSim[iCuts], locations, NONMEMBER_INDS=bumpThese)
        if n_elements(bumpThese) ne 0 then iCuts[bumpThese] += 1
    
        ;data as a function of z along the vertical line.
        ;   Transpose to be consistent with vertical cuts: [pos, locations]
        cut = transpose(data[icuts, iRange[0]:iRange[1]])

;-------------------------------------------------------
;Horizontal Cut ////////////////////////////////////////
;-------------------------------------------------------
    endif else begin
        ;Get the index range over which the vertcal cuts span
        iRange = getIndexRange(XSim, cut_xrange)
        pos = XSim[iRange[0]:iRange[1]]

        ;Z-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(YSim, locations)
        void = ismember(YSim[iCuts], locations, NONMEMBER_INDS=bumpThese)
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
function MrSim2::XZCuts, data_product, locations, pos, $
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
        iRange = getIndexRange(XSim, cut_xrange, STRIDE=stride)
        pos = XSim[iRange[0]:iRange[1]:stride]

        ;Z-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(ZSim, locations)
        void = ismember(ZSim[iCuts], locations, NONMEMBER_INDS=bumpThese)
        if n_elements(bumpThese) ne 0 then icuts[bumpThese] += 1

        ;data as a function of x along the horizontal line.
        cut = data[iRange[0]:iRange[1]:stride, iCuts]

;-------------------------------------------------------
;Vertical Cut //////////////////////////////////////////
;-------------------------------------------------------
    endif else begin
        ;Get the index range over which the vertcal cuts span
        iRange = getIndexRange(ZSim, cut_zrange, STRIDE=stride)
        pos = ZSim[iRange[0]:iRange[1]:stride]

        ;X-locations of the subset of vertical cuts to be displayed. If matches are
        ;not exact, round up instead of down.
        iCuts = value_locate(XSim, locations)
        void = ismember(XSim[iCuts], locations, NONMEMBER_INDS=bumpThese)
        if n_elements(bumpThese) ne 0 then icuts[bumpThese] += 1
    
        ;data as a function of z along the vertical line.
        ;   Transpose to be consistent with vertical cuts: [pos, locations]
        cut = transpose(data[icuts, iRange[0]:iRange[1]:stride])
    endelse
    
    return, cut
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
              directory:    '', $
              info:         ptr_new(), $
              ion_scale:    0B, $
              little:       0B, $
              mva_frame:    0B, $
              nSmooth:      0L, $
              orientation:  '', $
              time:         0L, $
              coord_system: '', $
              
              ;Sim Domain
              XRANGE: [0D, 0D], $
              XSim:   ptr_new(), $
              YSlice: 0L, $
              YRANGE: [0D, 0D], $
              YSim:   ptr_new(), $
              ZRANGE: [0D, 0D], $
              ZSim:   ptr_new(), $
            }
end