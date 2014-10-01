
; docformat = 'rst'
;
; NAME:
;       MrSim2D__Define
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
;   The purpose of this program is to read data out of '.gda' files created by one
;   of Bill Daughton's simulations.
;
;   NOTES:
;       dt*wci
;           Bill Daughton does not save every time slice. Therefore, the ti*dt*wci does
;           not produce the correct simulation time (ti is the time index and dt*wci is
;           a parameter from the info file). The parameter dt*wci is hard-coded in the
;           ::GetInfo method.
;
;       Electron Particle Data
;           Is typically stored one level up from the field and moment data (for which
;           the location is given by the DIRECTORY property). The file name takes the
;           form "electrons-t[t*wci].bin". If no file name is given to the ::ReadElectrons
;           method, this information will be used to construct a file name.
;
;           Also, Bill only saves every other particle. Therefore, anything computations
;           that lead to particle counts (e.g. distribution functions) should be
;           multiplied by 2.
;
;       MrSim_Which
;           This program is used to obtain all simulation-specific information.
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
;       02/14/2013  -   Modified from Li-Jen Chen's LJ_readInfo by Matthew Argall
;       04/29/2013  -   Converted to an object - MRA
;       04/30/2013  -   Added ::Fix, ::Fiz, ::Emag, ::Eperp, ::Jx, ::Jy, ::Jz, - MRA
;       05/01/2013  -   Added ::divPe_x, ::divPe_z, ::divPe_mag, ::Beta_e, ::Beta_i, 
;                           ::Beta, ::Bmag, ::EJpar. - MRA
;       2014/01/26  -   Removed parameters and keywords from secondary data product
;                           methods. Time and [XZ]Range can only be changed via the
;                           SetProperty method to better control data storage. For a
;                           given time and [XZ]Range, data is stored. Subsequent requests
;                           do not read data from the file. Data is cleared if time or
;                           [xz]range change. Removed Has_Data and Pull_Data methods.
;                           Added the ReadAsciiInfo method. - MRA
;       2014/01/27  -   No longer read the binary info file. Only the ascii one. - MRA
;       2014/01/29  -   Now is a subclasses MrSim. Many methods were moved to the
;                           MrSim class. - MRA.
;       2014/02/07  -   Renamed to MrSim2D__Define. Moved data product methods to the
;                           superclass. - MRA
;       2014/08/13  -   ReadData method was not orienting/selecting the data properly in
;                           MAGNETOPAUSE mode. Fixed. - MRA.
;-
;*****************************************************************************************
;+
;   This method initializes the MrSim2D class.
;-
function MrSim2D::INIT, theSim, time, yslice, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, 0
    endif

    ;Superclasses
    success = self -> MrSim::Init(theSim, time, _STRICT_EXTRA=extra)
    if success eq 0 then message, 'Unable to initialize MrSim superclass'

    ;Print information about the simulation
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
;   The purpose of this program is print a list of data products that can be made
;   using MrSim2D__define.pro
;-
pro MrSim2D::ListProducts
    compile_opt strictarr
    on_error, 2
    
    ;Which data product was requested?
    print, '-----------------------------------------------------------------------'
    print, '/////////////////////////////// LEGEND ////////////////////////////////'
    print, '-----------------------------------------------------------------------'
    print, 'A -- Magnetic Vector Potential'
    print, 'E -- Electric Field'
    print, 'B -- Magnetic Field'
    print, 'P -- Pressure'
    print, 'J -- Current Density'
    print, 'U -- Particle Velocity'
    print, 'n -- Number Density'
    print, 'F -- The Lorentz Force = E + (U cross B)'
    print, 'e -- electrons'
    print, 'i -- ions'
    print, '[x, y, z] -- spacial coordinates'
    print, 'dot -- the dot product'
    print, 'cross -- the cross product'
    print, ''
    print, '-----------------------------------------------------------------------'
    print, '//////////////////////////// DATA PRODUCTS ////////////////////////////'
    print, '-----------------------------------------------------------------------'
    print, 'Ay          -   Ay'
    print, 'Beta_e      -   Tr(Pe) / (|B|^2 / 2 mu0) -- electron plasma beta'
    print, 'Beta_i      -   Tr(Pi) / (|B|^2 / 2 mu0) -- ion plasma beta'
    print, 'Beta_p      -   Tr(Pi + Pe) / (|B|^2 / 2 mu0) -- plasma beta'
    print, 'Bmag        -   |B|'
    print, 'Bx          -   Bx'
    print, 'By          -   By'
    print, 'Bz          -   Bz'
    print, 'divPe_mag   -   |div Pe|'
    print, 'divPe_x     -   (div Pe)_x'
    print, 'divPe_z     -   (div Pe)_z'
    print, 'E_dot_B     -   E dot B'
    print, 'E_dot_J     -   E dot J'
    print, 'EJpar       -   (E dot B/|B|) * (J dot B/|B|)'
    print, 'Emag        -   |E|'
    print, 'Epar        -   E dot B-hat'
    print, 'ExB_mag     -   |E cross B|'
    print, 'ExB_x       -   E cross B, x-component'
    print, 'ExB_y       -   E cross B, y-component'
    print, 'ExB_z       -   E cross B, y-component'
    print, 'Ex          -   Ex'
    print, 'ExJx        -   Ex * Jx'
    print, 'Ey          -   Ey'
    print, 'EyBy        -   Ey * By'
    print, 'EyJy        -   Ey * Jy'
    print, 'Ez          -   Ez'
    print, 'EzJz        -   Ez * Jz'
    print, 'Fe_mag      -   |Fe|'
    print, 'Fex         -   Fex'
    print, 'Fey         -   Fey'
    print, 'Fez         -   Fez'
    print, 'Fi_mag      -   |Fi|'
    print, 'Fix         -   Fix'
    print, 'Fiy         -   Fiy'
    print, 'Fiz         -   Fiz'
    print, 'Je_inPlane  -   -1.0 * ne * (Uex^2 + Uez^2) -- In-plane current density'
    print, 'Jex         -   -1.0 * ne * Uex'
    print, 'Jey         -   -1.0 * ne * Uey'
    print, 'Jez         -   -1.0 * ne * Uez'
    print, 'Jix         -   ni * Uix'
    print, 'Jiy         -   ni * Uiy'
    print, 'Jiz         -   ni * Uiz'
    print, 'Jpar        -   J dot B/|B|'
    print, 'Jx          -   (-ne * Uex) + (ni * Uix)'
    print, 'Jy          -   (-ne * Uey) + (ni * Uiy)'
    print, 'Jz          -   (-ne * Uez) + (ni * Uiz)'
    print, 'ne_ni       -   ne / ni'
    print, 'ne_ratio    -   (ne - ni) / ne'
    print, 'ne          -   ne'
    print, 'ni_ratio    -   (ne - ni) / ni'
    print, 'ni          -   ni'
    print, 'Pe-xx       -   Pe-xx'
    print, 'Pe-xy       -   Pe-xy'
    print, 'Pe-xz       -   Pe-xz'
    print, 'Pe-yy       -   Pe-yy'
    print, 'Pe-yz       -   Pe-yz'
    print, 'Pe-zz       -   Pe-zz'
    print, 'Pi-xx       -   Pi-xx'
    print, 'Pi-xy       -   Pi-xy'
    print, 'Pi-xz       -   Pi-xz'
    print, 'Pi-yy       -   Pi-yy'
    print, 'Pi-yz       -   Pi-yz'
    print, 'Pi-zz       -   Pi-zz'
    print, 'Ue_inPlane  -   sqrt(Uex^2 + Uey^2) -- in-plane electron velocity'
    print, 'Ue_mag      -   Ue dot Ue'
    print, 'Uex         -   Uex'
    print, 'Uey         -   Uey'
    print, 'Uez         -   Uez'
    print, 'Uix         -   Uix'
    print, 'Uiy         -   Uiy'
    print, 'Uiz         -   Uiz'
    print, 'Vae         -   |B| / sqrt(ne) -- Alfven Velocity'
end


;+
;   The purpose of this method is to retrieve data from the info file.
;
; :Keywords:
;       DTXWCI:             out, optional, type=double
;                           Time step multiplied by the ion gyrofrequency. This is the
;                               time interval at which data is saved.
;       INFO_DTXWCI:        out, optional, type=double
;                           Actual simulation time step multiplied by the ion gyrofrequency.
;       ECOUNTFACTOR:       out, optional, type=long
;                           Bill Daughton does not save every particle, only every-other.
;                               When dealing with electron counts, you must multiply by
;                               this factor to obtain the actual total counts.
;       _REF_EXTRA:         out, optional, type=any
;                           Any keyword accepted by MrSim::GetProperty is also
;                               accepted for keyword inheritance.
;-
pro MrSim2D::GetInfo, $
ECOUNTFACTOR=eCountFactor, $
DTXWCI=dtxwci, $
INFO_DTXWCI=dtxwci_info, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif

    ;Hard code -- the one simulation that we have has a subset of time indices.
    get_dtxwci = arg_present(dtxwci)
    get_eCF    = arg_present(eCountFactor)
    if get_dtxwci + get_eCF gt 0 then begin
        MrSim_Which, self.simnum, TINDEX=self.time, ECOUNTFACTOR=eCountFactor, DTXWCI=dtxwci
        
        ;DTXWCI
        if get_dtxwci && n_elements(dtxwci) eq 0 then begin
            message, 'dt*wci unknown for simulatio #' + strtrim(self.simnum, 2) + '. ' + $
                     'Setting dt*wci = 2', /INFORMATIONAL
            dtxwci = 2
        endif
        
        ;ECOUNTFACTOR
        if get_eCF && n_elements(eCountFactor) eq 0 then begin
            message, 'No explicit eCountFactor set. Setting eCountFactor = 1.', /INFORMATIONAL
            eCountFactor = 1
        endif
    endif
    
    if n_elements(extra) gt 0 $
        then self -> MrSim::GetInfo, DTXWCI=dtxwci_info, _STRICT_EXTRA=extra
end


;+
;   The purpose of this method is to read the ASCII "info" file relating to Bill
;   Daughton's simulations.
;
; :Private:
;
; :Params:
;       FILENAME:           in, optional, type=string, default='electrons-t' [time index] + '.gda'
;                           Name of the "info" file to be read.
;-
pro MrSim2D::ReadElectrons, filename
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(pwd) gt 0 then cd, pwd
        void = cgErrorMsg()
        return
    endif

;---------------------------------------------------------------------
; Create a File Name /////////////////////////////////////////////////
;---------------------------------------------------------------------
    if n_elements(filename) eq 0 || filename eq '' then begin
        MrSim_Which, self.simnum, EFILE=filename, FMAP_DIR=fMap_dir, TINDEX=self.time
        if n_elements(filename) eq !Null then filename = ''
        if n_elements(fMap_dir) eq !Null then fMap_dir = ''

        ;Pick a file if the guessed file is not valid.
        if file_test(filename) eq 0 then filename = ''
        if filename eq '' then filename = dialog_pickfile(/READ, TITLE='Pick electron file.')
        if fmap_dir eq '' then fmap_dir = dialog_pickfile(/READ, TITLE='Pick fMap directory.', /DIRECTORY)
    
    ;Ensure the file exists.
    endif else begin
        if file_test(filename) eq 0 then $
            message, 'File does not exist: "' + filename + '".'
    endelse

;---------------------------------------------------------------------
; Read Electron Data /////////////////////////////////////////////////
;---------------------------------------------------------------------
    *self.electrons = MrSim_ReadParticles(filename, self.xrange, self.zrange, $
                                          FMAP_DIR=fMap_dir, /VERBOSE)
end


;+
;   The purpose of this program is to read data from a ".gda" file.
;
; :Private:
;
; :Params:
;       NAME:                   in, required, type=string, default=
;                               The name of the ".gda" data file (without the ".gda"
;                                   file extension). GDA files are typically named after
;                                   the parameter whose data they contain.
;-
pro MrSim2D::ReadData, name
    compile_opt strictarr

    ;Catch any errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif

    ;Read the data
    data = self -> ReadMoment(name)
    
    ;Store the data so the file does not have to be accessed again.
    self -> SetData, name, temporary(data)
end


;+
;   The purpose of this program is to read data from a ".gda" file.
;
; :Private:
;
; :Params:
;       NAME:                   in, required, type=string, default=
;                               The name of the ".gda" data file (without the ".gda"
;                                   file extension). GDA files are typically named after
;                                   the parameter whose data they contain.
;-
function MrSim2D::ReadMoment, name, tindex, xrange, zrange, $
NSMOOTH=nSmooth
    compile_opt strictarr

    ;Catch any errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif

    ;Defaults
    if n_elements(nSmooth) eq 0 then nSmooth = self.nSmooth
    if n_elements(tindex)  eq 0 then tindex  = self.time
    if n_elements(xrange)  eq 0 then xrange  = self.xrange
    if n_elements(zrange)  eq 0 then zrange  = self.zrange

    ;Check if NAME has an underscore. Replace it with a hyphen (i.e. "Pe_xx" -> "Pe-xx")
    _name = stregex(name, '_', /BOOLEAN) ? strjoin(strsplit(name, '_'), '-') : name
    
    ;Check if the file exists. If not, throw an error.
    f_name = filepath(_name + '.gda', ROOT_DIR=self.directory)
    if file_test(f_name) eq 0 then begin
        f_name = cgPickFile(TITLE='Choose .gda Field or Moment File.', /READ)
        if f_name eq '' then return, !Null
    endif

;-------------------------------------------------------
; Coordinate System -> Simulation Coordinates //////////
;-------------------------------------------------------
    ;Get the index range of data to read
    ixrange = getIndexRange(*self.XSim, xrange, STRIDE=xstride)
    izrange = getIndexRange(*self.ZSim, zrange, STRIDE=zstride)
    
    ;Subset of data to read
    ;   - Data is still in SIMULATION coordinates.
    ;   - Ranges are in COORD_SYSTEM coordinates.
    ;   - Must traslate from COORD_SYSTEM to SIMULATION.
    case self.coord_system of
        'SIMULATION': ;Do nothing
        
        ;Interchange indices
        ;   x -> -z (negation occurs below)
        ;   y -> +y
        ;   z -> +x
        'MAGNETOPAUSE': begin
            ;Interchange x and z
            iTemp   = ixrange
            ixrange = izrange
            izrange = temporary(iTemp)
        endcase
        
        ;Interchange indices
        ;   x -> -x (negation occurs below)
        'MAGNETOTAIL': ;Do nothing
        
        else: ;Nothing
    endcase

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------

    ;Get the system size
    self -> GetInfo, NX=nx, NZ=nz
    
    if xstride gt 0 then begin
        ixMin = ixrange[0]
        ixMax = ixrange[1]
    endif else begin
        ixMin = ixrange[1]
        ixMax = ixrange[0]
    endelse
        
    if zstride gt 0 then begin
        izMin = izrange[0]
        izMax = izrange[1]
    endif else begin
        izMin = izrange[1]
        izMax = izrange[0]
    endelse
    
    ;Allocate memory to output array
    data = fltarr(ixMax-ixMin+1, izMax-izMin+1)
    temp = fltarr(ixMax-ixMin+1)

    ;Open the data file.
    openr, lun, f_name, /SWAP_IF_BIG_ENDIAN, /GET_LUN
    
    ;Jump to the proper y-slice
    ;   - each data chunk is nx*nz large and followed by ti, t*wci
    ;   - Jump to the proper time-slice
    ;   - Jump to the proper x- and z-index within that time-slice
    tOffset = (ulong64(nx*nz) + 2ULL) * tindex * 4ULL
    ix      = ulong64(ixMin)
    for iz = izMin, izMax do begin
        offset = tOffset + 4ULL * (ix + iz*nx)
        point_lun, lun, offset
        readu, lun, temp
        data[*,iz-izMin] = temp
    endfor
    
    ;Free the file
    free_lun, lun

;-------------------------------------------------------
; Change from SIMULATION Coordinates? //////////////////
;-------------------------------------------------------
    ;Change coordinate Systems?
    case self.coord_system of
        'SIMULATION': ;Do nothing
        
        'MAGNETOPAUSE': begin
            ;Re-orient the axes
            ;   - x and z have been switch from SIMULATION coordinates arleady (above).
            ;   - IDL draws the bottom row of the image first: img[*,0]
            ;       * When transpose, the bottom row becomes the left-most column.
            ;   - Ordering bottom axis as [+x, -x] requires two steps.
            ;       * Negate all data products with a single "z" in them (below).
            ;       * Reverse the x-axis itself (::MakeSimDomain method).
            ;
            ;                            ++z +--------|
            ;                                |        |
            ;           MSP                  |        |
            ; +x  *--------------|         M |        | M
            ;     |              |   ===>  S |        | S
            ; -x  |--------------+         H |        | P
            ;     +z    MSH    ++z           |        |
            ;                                |        |
            ;                             +z |--------*
            ;                               -x       +x
            ;
            case self.orientation of
                'XY': ;Do nothing
                'XZ': data = transpose(data)
                else: message, 'Orienatation "' + self.orientation + '" not recognized.'
            endcase

            ;Now reverse the z-axes
            ;   Single negate Bz, Uiz, Pe-xz, Pe-yz, etc.
            case 1 of
                stregex(_name, '[A-Z][a-z]*z$',   /BOOLEAN): data = -data
                stregex(_name, '-(z[yx]|[xy]z)$', /BOOLEAN): data = -data
                else: ;Do nothing
            endcase
        end
        
        'MAGNETOTAIL': begin
            ;Re-orient the axes
            ;   - x and z are the same orientation as simulation coordinates
            ;   - Ordering bottom axis as [+x, -x] requires two steps.
            ;       * Negate all data products with a single "x" in them (below).
            ;       * Reverse the x-axis itself (::MakeSimDomain method).
            case self.orientation of
                'XY': ;Do nothing
                'XZ': ;Do nothing
                else: message, 'Orienatation "' + self.orientation + '" not recognized.'
            endcase

            ;Now reverse the x-axes
            ;   single negate Pe-xy, Pe-xz, etc.
            ;   double negate Pe-xx
            case 1 of
                stregex(_name, '[A-Z][a-z]*(x)$', /BOOLEAN): data = -data
                stregex(_name, '-(y|z)[x]$',      /BOOLEAN): data = -data
                stregex(_name, '-[x](y|z)$',      /BOOLEAN): data = -data
                else: ;Do nothing
            endcase
        end
        
        else: message, 'Coordinate system "' + self.coord_system + '" not recognized.'
    endcase
    
    ;Smooth the data
    if nSmooth gt 0 then data = smooth(data, nSmooth, /EDGE_TRUNCATE)
    
    return, data
end


;+
;   The purpose of this program is to read data from a ".gda" file.
;
; :Private:
;
; :Params:
;       NAME:                   in, required, type=string, default=
;                               The name of the ".gda" data file (without the ".gda"
;                                   file extension). GDA files are typically named after
;                                   the parameter whose data they contain.
;-
function MrSim2D::ReadMoment_v1, name, tindex, xrange, zrange, $
NSMOOTH=nSmooth
    compile_opt strictarr

    ;Catch any errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif

    ;Defaults
    if n_elements(nSmooth) eq 0 then nSmooth = self.nSmooth
    if n_elements(tindex)  eq 0 then tindex  = self.time
    if n_elements(xrange)  eq 0 then xrange  = self.xrange
    if n_elements(zrange)  eq 0 then zrange  = self.zrange

    ;Check if NAME has an underscore. Replace it with a hyphen (i.e. "Pe_xx" -> "Pe-xx")
    _name = stregex(name, '_', /BOOLEAN) ? strjoin(strsplit(name, '_'), '-') : name
    
    ;Check if the file exists. If not, throw an error.
    f_name = filepath(_name + '.gda', ROOT_DIR=self.directory)
    if file_test(f_name) eq 0 then begin
        f_name = cgPickFile(TITLE='Choose .gda Field or Moment File.', /READ)
        if f_name eq '' then return, !Null
    endif

;-------------------------------------------------------
; Coordinate System -> Simulation Coordinates //////////
;-------------------------------------------------------
    ;Get the index range of data to read
    ixrange = getIndexRange(*self.XSim, xrange)
    izrange = getIndexRange(*self.ZSim, zrange)
    
    ;Subset of data to read
    ;   - Data is still in SIMULATION coordinates.
    ;   - Ranges are in COORD_SYSTEM coordinates.
    ;   - Must traslate from COORD_SYSTEM to SIMULATION.
    case self.coord_system of
        'SIMULATION': ;Do nothing
        
        ;Interchange indices
        ;   x -> -z (negation occurs below)
        ;   y -> +y
        ;   z -> +x
        'MAGNETOPAUSE': begin
            ;Interchange x and z
            iTemp   = ixrange
            ixrange = izrange
            izrange = temporary(iTemp)
        endcase
        
        ;Interchange indices
        ;   x -> -x (negation occurs below)
        'MAGNETOTAIL': ;Do nothing
        
        else: ;Nothing
    endcase

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------

    ;Get the system size
    self -> GetInfo, NX=nx, NZ=nz

    ;Create an associated variable for the data file. Each data image is [nx] x [nz]
    ;points large and there is 1 time-stamp associated with each image.
    struct = { data: fltarr(nx, nz), $
               time: 0.0, $
               it:   500000 $
             }

    ;Open the data file.
    openr, lun, f_name, /SWAP_IF_BIG_ENDIAN, /GET_LUN
    
    ;FIELD will be an array of structures with one [nx] x [ny] image and one time-stamp
    ;per element. STRUCT will be the image at the desired time step of the simulation.
    field  = assoc(lun, struct)
    struct = field[tindex]
    free_lun, lun

    ;Trim the data
    data   = struct.data[ixrange[0]:ixrange[1], izrange[0]:izrange[1]]
    struct = !Null

;-------------------------------------------------------
; Change from SIMULATION Coordinates? //////////////////
;-------------------------------------------------------
    ;Change coordinate Systems?
    case self.coord_system of
        'SIMULATION': ;Do nothing
        
        'MAGNETOPAUSE': begin
            ;Re-orient the axes
            ;   - x and z have been switch from SIMULATION coordinates arleady (above).
            ;   - IDL draws the bottom row of the image first: img[*,0]
            ;       * When transpose, the bottom row becomes the left-most column.
            ;   - Ordering bottom axis as [+x, -x] requires two steps.
            ;       * Negate all data products with a single "z" in them (below).
            ;       * Reverse the x-axis itself (::MakeSimDomain method).
            ;
            ;                            ++z +--------|
            ;                                |        |
            ;           MSP                  |        |
            ; +x  *--------------|         M |        | M
            ;     |              |   ===>  S |        | S
            ; -x  |--------------+         H |        | P
            ;     +z    MSH    ++z           |        |
            ;                                |        |
            ;                             +z |--------*
            ;                               -x       +x
            ;
            case self.orientation of
                'XY': ;Do nothing
                'XZ': data = transpose(data)
                else: message, 'Orienatation "' + self.orientation + '" not recognized.'
            endcase

            ;Now reverse the z-axes
            ;   Single negate Bz, Uiz, Pe-xz, Pe-yz, etc.
            case 1 of
                stregex(_name, '[A-Z][a-z]*z$',   /BOOLEAN): data = -data
                stregex(_name, '-(z[yx]|[xy]z)$', /BOOLEAN): data = -data
                else: ;Do nothing
            endcase
        end
        
        'MAGNETOTAIL': begin
            ;Re-orient the axes
            ;   - x and z are the same orientation as simulation coordinates
            ;   - Ordering bottom axis as [+x, -x] requires two steps.
            ;       * Negate all data products with a single "x" in them (below).
            ;       * Reverse the x-axis itself (::MakeSimDomain method).
            case self.orientation of
                'XY': ;Do nothing
                'XZ': ;Do nothing
                else: message, 'Orienatation "' + self.orientation + '" not recognized.'
            endcase

            ;Now reverse the x-axes
            ;   single negate Pe-xy, Pe-xz, etc.
            ;   double negate Pe-xx
            case 1 of
                stregex(_name, '[A-Z][a-z]*(x)$', /BOOLEAN): data = -data
                stregex(_name, '-(y|z)[x]$',      /BOOLEAN): data = -data
                stregex(_name, '-[x](y|z)$',      /BOOLEAN): data = -data
                else: ;Do nothing
            endcase
        end
        
        else: message, 'Coordinate system "' + self.coord_system + '" not recognized.'
    endcase
    
    ;Smooth the data
    if nSmooth gt 0 then data = smooth(data, nSmooth, /EDGE_TRUNCATE)
    
    return, data
end


;+
;   The purpose of this method is to set object properties.
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
pro MrSim2D::SetSim, sim_id, $
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
    
    ;Make sure we are switching to a 2D simulation.
    MrSim_Which, sim_id, DIMENSION=dimension
    if dimension ne '2D' then message, 'Cannot change dimensions from 2D to ' + dimension + '.'
    
    ;Call the superclass
    self -> MrSim::SetSim, sim_id, BINARY=binary, DIRECTORY=directory, INFO_FILE=info_file
end


;+
;   Clean up after the object is destroyed.
;-
pro MrSim2D::cleanup
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    self -> MrSim::CleanUp
end


;+
;   This definition statement for the MrSim2D class.
;-
pro MrSim2D__define
    compile_opt strictarr
    
    class = { MrSim2D, $
              inherits MrSim $
            }
end