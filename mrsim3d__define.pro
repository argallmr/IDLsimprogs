; docformat = 'rst'
;
;+
;   The purpose of this program is to provide an interface for the Asymmetric 3D magnetic
;   reconnection simulations provided by Bill Daughton from Los Alamos
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
;       2014/01/30  -   Written by Matthew Argall
;       2014/02/07  -   Contents of the ::MakeSimDomain method were moved to the
;                           superclass, along with the YRANGE, YSim and Orientation
;                           properties. - MRA
;       2014/02/16  -   Moved the YSLICE property into the superclass. - MRA
;       2014/03/29  -   Added the divPe_y, divPi_y, gradPe, gradPi, Pe, and Pi methods. - MRA
;       2014/04/08  -   Added the D_DY and divPi_y methods. - MRA
;-
;*****************************************************************************************
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
function MrSim3D::D_DY, data_product, $
OVERWRITE=overwrite
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Check if the data requested is a standard data product
    if self -> HasData(data_product) eq -1 then $
        message, 'DATA_PRODUCT must be a standard data product with a binary .gda file.'
    
    ;Check the orientation
    if self.orientation ne 'XZ' then $
        message, 'WARNING: Y-derivatives are only performed in the XZ-plane.' + $
                 ' The current orientation is "' + self.orientation + '".'
    
    ;Defaults
    overwrite = keyword_set(overwrite)

    ;Get the data
    data = self -> Read3D(data_product)
    
    ;Get the x-size of the simulation in electron skin depths
    self -> GetInfo, DY_DE=dy_de
    dims = size(data, /DIMENSIONS)

    ;Overwrite?
    if overwrite then begin
        ;Take the central difference
        data = reform(data[*,2,*] - data[*,0,*]) / (2.0 * dy_de)
        return, data
    endif
    
    ;Make a copy?
    derivative = fltarr(dims[[0,2]])
    derivative = reform(data[*,2,*] - data[*,0,*]) / (2.0 * dy_de)
    
    return, derivative
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
function MrSim3D::GetData, data_product, $
DX=dx, $
DY=DY, $
DZ=dz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif

    dx = keyword_set(dx)
    dy = keyword_set(dy)
    dz = keyword_Set(dz)
    if dx + dy + dz gt 1 then message, 'DX, DY, and DZ are mutually exclusive.'

;-------------------------------------------------------
; Y-Derivative? ////////////////////////////////////////
;-------------------------------------------------------
    if dy then begin
        data = self -> D_DY(data_product)
        return, data
    endif

;-------------------------------------------------------
; Check 3D Data Products Before Fowarding //////////////
;-------------------------------------------------------
    case strupcase(data_product) of
        ;3D Data Products
        'GRADPE':   data = self -> gradPe()
        'GRADPE_X': data = self -> gradPe_x()
        'GRADPE_Y': data = self -> gradPe_y()
        'GRADPE_Z': data = self -> gradPe_z()
        'GRADPI':   data = self -> gradPi()
        'DIVPE_Y':  data = self -> divPe_y()
        'DIVPI_Y':  data = self -> divPi_y()
        'PE':       data = self -> Pe()
        'PI':       data = self -> Pi()
        else:       data = self -> MrSim::GetData(data_product)
    endcase

    ;Take the derivative?    
    if dx then data = self -> D_DX(data, /OVERWRITE)
    if dz then data = self -> D_DZ(data, /OVERWRITE)
    
    return, data
end


;+
;   The purpose of this method is to retrieve data from the info file.
;
; :Keywords:
;-
pro MrSim3D::GetInfo, $
DX_DE=dx_de, $
DY_DE=dy_de, $
DZ_DE=dz_de, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Data has been downsampled (see comments in the INIT method). Calculate the modified
    ;size of each cell in de (the true size is stored in (*self.info).dx_de).
    if arg_present(dx_de) then dx_de = (*self.info).lx_de / (*self.info).nx
    if arg_present(dy_de) then dy_de = (*self.info).ly_de / (*self.info).ny
    if arg_present(dz_de) then dz_de = (*self.info).lz_de / (*self.info).nz

    ;More Info
    if n_elements(extra) gt 0 then self -> MrSim::GetInfo, _STRICT_EXTRA=extra
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
pro MrSim3D::ReadElectrons, filename
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        void = cgErrorMsg()
        return
    endif

;---------------------------------------------------------------------
; Create a File Name /////////////////////////////////////////////////
;---------------------------------------------------------------------
    if n_elements(filename) eq 0 || filename eq '' then begin
        MrSim_Which, self.simnum, EFILE=filename, FMAP_DIR=fMap_dir, TINDEX=self.time, YSLICE=self.yslice
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
    ;Read the simulation info file
    *self.electrons = MrSim_ReadParticles(filename, self.xrange, self.zrange, FMAP_DIR=fmap_dir, /VERBOSE)
end


;+
;   The purpose of this method is to read data. If the data file does not exist, an
;   informational message will be printed and no data will be read.
;
; :Params:
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;-
function MrSim3D::ReadMoment, name, tIndex, xrange, yslice, zrange
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMsg()
        return, !Null
    endif

    ;Defaults
    if n_elements(nSmooth) eq 0 then nSmooth = self.nSmooth
    if n_elements(tindex)  eq 0 then tindex  = self.time
    if n_elements(xrange)  eq 0 then xrange  = self.xrange
    if n_elements(yslice)  eq 0 then yslice  = self.yslice
    if n_elements(zrange)  eq 0 then zrange  = self.zrange

    ;Create file name and test for the file
    file = filepath(name + '_' + string(tindex, FORMAT='(i0)') + '.gda', ROOT_DIR=self.directory)
    if file_test(file) eq 0 then begin
        file = cgPickFile(TITLE='Choose .gda Field or Moment File.', /READ)
        if file eq '' then return, !Null
    endif

;-------------------------------------------------------
; Coordinate System -> Simulation Coordinates //////////
;-------------------------------------------------------
    ;Subset of data to read
    ;   - Data is still in SIMULATION coordinates.
    ;   - Must recalculate the indices in simulation coordinates.
    case self.coord_system of
        'SIMULATION':   self -> GetProperty, IXRANGE=ixrange, IYRANGE=iyrange, IZRANGE=izrange
        
        ;Interchange indices
        ;   x -> -z
        ;   y -> -y
        ;   z -> +x
        'MAGNETOPAUSE': begin
            ;Interchange x and z
            self -> GetProperty, IXRANGE=izrange, IYRANGE=iyrange, IZRANGE=ixrange
            xDim = n_elements(*self.ZSim)
            yDim = n_elements(*self.YSim)
            zDim = n_elements(*self.XSim)
            
            ;Negate by picking from the opposite end of the simulation
            iyrange = yDim - iyrange - 1
            izrange = zDim - izrange - 1
            
            ;Order [min, max]
            iyrange = iyrange[sort(iyrange)]
            izrange = izrange[sort(izrange)]
        endcase
        
        else: ;Nothing
    endcase

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------
    ;
    ;Now loop through the 3D array in the data file and read only the
    ;needed y-slice of the data file - assume array is in FORTRAN ordering
    ;
    self -> GetInfo, NX=nx, NY=ny, NZ=nz
    
    if izrange[0] lt izrange[1] then begin
        izMin = izrange[0]
        izMax = izrange[1]
    endif else begin
        izMin = izrange[1]
        izMax = izrange[0]
    endelse
    
    if ixrange[0] lt ixrange[1] then begin
        ixMin = ixrange[0]
        ixMax = ixrange[1]
    endif else begin
        ixMin = ixrange[1]
        ixMax = ixrange[0]
    endelse
    
;    if iyrange[0] lt iyrange[1] then begin
;        iyMin = iyrange[0]
;        iyMax = iyrange[1]
;    endif else begin
;        iyMin = iyrange[1]
;        iyMax = iyrange[0]
;    endelse

    ;Open File and read data
    openr, lun, file, /GET_LUN

    ;Read x-z slice at constant y
    case self.orientation of
        'XZ': begin
            ;Declare output array for data
            slice = fltarr(ixMax-ixMin+1, izMax-izMin+1)

            ;Read this many elements from a row.
            temp = fltarr(ixMax-ixMin+1)

            ;Read the data
            ix = ulong64(ixMin)
            iy = ulong64(self.yslice)
            for iz = ulong64(izMin), izMax do begin
                ;Offset to first record.
                ;   ix       = column offset
                ;   iy*ny    = row offset
                ;   iz*nx*ny = 3rd dimension offset
                ;   4        = Bytes per element
                offset = 4ULL * (ix + iy*nx + iz*nx*ny)
                
                ;Read and save the records.
                point_lun, lun, offset
                readu, lun, temp
                slice[*,iz-izMin] = temp
            endfor
        endcase
        
        'XY': begin
            ;Declare array for data
            slice = fltarr(ixMax-ixMin+1, iyMax-iyMin+1)
            temp  = fltarr(ixMax-ixMin+1)
            
            ;Here, we must pick out single values of X within different rows.
            ;   - Jump to the proper z-slice
            ix      = ulong64(ixMin)
            iz      = ulong64(izMin)
            zOffset = iz*nx*ny*4ULL
            
            ;Read the data
            ;   - Jump to the proper row (y)
            ;   - Read all of the requested columns (x)
            for iy = ulong64(iyMin), iyMax do begin
                point_lun, lun, zOffset + 4ULL*(ix + iy*nx)
                readu, lun, temp
                slice[*,iy] = temp
            endfor
        endcase
        
        'YZ': begin
            ;Declare array for data
            slice = fltarr(iyMax-iyMin+1, izMax-izMin+1)
            temp  = fltarr(1)
            
            ;Here, we must pick out single values of X within different rows.
            ix     = ulong64(ixMin)
            for iz = ulong64(izMin), izMin+izMax-1 do begin
                for iy = ulong64(iyMin), iyMin+iyMax-1 do begin
                    ;Offset to first record.
                    ;   ix       = column offset
                    ;   iy*ny    = row offset
                    ;   iz*nx*ny = 3rd dimension offset
                    ;   4        = Bytes per element
                    offset = 4ULL * (ix + iy*nx + iz*nx*ny)
                
                    ;Read and save the records.
                    point_lun, lun, offset
                    readu, lun, temp
                    slice[iy-iyMin, iz-izMin] = temp
                endfor
            endfor
        endcase
    
        else: ;Nothing
    endcase  

    ;Close file
    free_lun, lun

;-------------------------------------------------------
; Change from SIMULATION Coordinates? //////////////////
;-------------------------------------------------------
    case self.coord_system of
        'SIMULATION':   ;Do nothing
        
        'MAGNETOPAUSE': begin
            ;Re-orient the axes
            case self.orientation of
                'XY': ;Do nothing
                'XZ': slice = transpose(slice)
                else: message, 'Orienatation "' + self.orientation + '" not recognized.'
            endcase
            
            ;Now reverse the  z-axes
            ;   single negate Pe-xz, Pe-yz, etc.
            ;   double negate Pe-zz, etc.
            case 1 of
                stregex(name, '[A-Z][a-z]*z$',   /BOOLEAN): slice = -slice
                stregex(name, '-(z[yx]|[xy]z)$', /BOOLEAN): slice = -slice
                else: ;Do nothing
            endcase
        end
        
        else: ;Nothing
    endcase
    
    ;Smooth data 
    if self.nSmooth gt 0 then slice = smooth(slice, self.nSmooth, /EDGE_TRUNCATE)  

    ;Store the data so that it does not have to be read again.
    return, slice
end


;+
;   The purpose of this method is to read data. If the data file does not exist, an
;   informational message will be printed and no data will be read.
;
;   Original method of reading field and moment data.
;
; :Params:
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;-
function MrSim3D::ReadMoment_Original, name, tIndex, xrange, yslice, zrange
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMsg()
        return, !Null
    endif

    ;Defaults
    if n_elements(nSmooth) eq 0 then nSmooth = self.nSmooth
    if n_elements(tindex)  eq 0 then tindex  = self.time
    if n_elements(xrange)  eq 0 then xrange  = self.xrange
    if n_elements(yslice)  eq 0 then yslice  = self.yslice
    if n_elements(zrange)  eq 0 then zrange  = self.zrange

    ;Create file name and test for the file
    file = filepath(name + '_' + string(tindex, FORMAT='(i0)') + '.gda', ROOT_DIR=self.directory)
    if file_test(file) eq 0 then begin
        file = cgPickFile(TITLE='Choose .gda Field or Moment File.', /READ)
        if file eq '' then return, !Null
    endif

;-------------------------------------------------------
; Coordinate System -> Simulation Coordinates //////////
;-------------------------------------------------------
    ;Subset of data to read
    ;   - Data is still in SIMULATION coordinates.
    ;   - Ranges are in COORD_SYSTEM coordinates.
    ;   - Must traslate from COORD_SYSTEM to SIMULATION.
    case self.coord_system of
        'SIMULATION':   self -> GetProperty, IXRANGE=ixrange, IYRANGE=iyrange, IZRANGE=izrange
        
        ;Interchange indices
        ;   x -> -z
        ;   y -> -y
        ;   z -> +x
        'MAGNETOPAUSE': begin
            ;Interchange x and z
            self -> GetProperty, IXRANGE=izrange, IYRANGE=iyrange, IZRANGE=ixrange
            xDim = n_elements(*self.ZSim)
            yDim = n_elements(*self.YSim)
            zDim = n_elements(*self.XSim)
            
            ;Negate by picking from the opposite end of the simulation
            iyrange = yDim - iyrange - 1
            izrange = zDim - izrange - 1
            
            ;Order [min, max]
            iyrange = iyrange[sort(iyrange)]
            izrange = izrange[sort(izrange)]
        endcase
        
        else: ;Nothing
    endcase

;-------------------------------------------------------
; Read Data ////////////////////////////////////////////
;-------------------------------------------------------
    ;
    ;Now loop through the 3D array in the data file and read only the
    ;needed y-slice of the data file - assume array is in FORTRAN ordering
    ;
    
    ;Dimension sizes.
    self -> GetInfo, NX=nx, NY=ny, NZ=nz

    ;Open File and read data
    openr, lun, file, /GET_LUN

    ;Read x-z slice at constant y
    case self.orientation of
        'XZ': begin
            if izrange[0] lt izrange[1] then begin
                izMin = izrange[0]
                izMax = izrange[1]
            endif else begin
                izMin = izrange[1]
                izMax = izrange[0]
            endelse
        
            ;Declare array for data
            slice = fltarr(nx,izMax-izMin+1)

            ;Get data from the file
            template = fltarr(nx)
            dfile = assoc(lun, template)

            ;Get the desired subset of data. DFILE access cannot be vectorized
            ;because it is a file... X values seem to be stored in the rows
            ;while y and z combinations make up the row numbers.
            j = self.yslice
            for k = izMin, izMax do begin
                n = (j + ny*k)
                slice[*,k-izMin] = dfile[n]
            endfor
            
            ;Trim the X component
            slice = slice[ixrange[0]:ixrange[1],*]
        endcase
        
        'XY': begin
            ;Declare array for data
            slice = fltarr(nx,ny)
        
            stemplate = fltarr(1)
            dfile = assoc(lun, stemplate)
            
            ;Here, we must pick out single values of X within different rows.
            i = self.yslice
            for k = 0, nz-1 do begin
                for j = 0, ny-1 do begin
                    n = (i + nx*(j + ny*k))
                    slice[j,k] = dfile[n]
                endfor
            endfor
        endcase
    
        else: ;Nothing
    endcase   

    ;Close file
    free_lun, lun

;-------------------------------------------------------
; Change from SIMULATION Coordinates? //////////////////
;-------------------------------------------------------
    case self.coord_system of
        'SIMULATION':   ;Do nothing
        
        'MAGNETOPAUSE': begin
            ;Re-orient the axes
            case self.orientation of
                'XY': ;Do nothing
                'XZ': slice = transpose(slice)
                else: message, 'Orienatation "' + self.orientation + '" not recognized.'
            endcase
            
            ;Now reverse the  z-axes
            ;   single negate Pe-xz, Pe-yz, etc.
            ;   double negate Pe-zz, etc.
            case 1 of
                stregex(name, '[A-Z][a-z]*z$',   /BOOLEAN): slice = -slice
                stregex(name, '-(z[yx]|[xy]z)$', /BOOLEAN): slice = -slice
                else: ;Do nothing
            endcase
        end
        
        else: ;Nothing
    endcase
    
    ;Smooth data 
    if self.nSmooth gt 0 then slice = smooth(slice, self.nSmooth, /EDGE_TRUNCATE)

    return, slice
end


;+
;   The purpose of this method is to read data. If the data file does not exist, an
;   informational message will be printed and no data will be read.
;
; :Params:
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;-
pro MrSim3D::ReadData, name
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMsg()
        return
    endif

    ;Read the data
    data = self -> ReadMoment(name)
    
    ;Store the data so the file does not have to be accessed again.
    self -> SetData, name, temporary(data)
end


;+
;   The purpose of this method is to read data from 3 adjacent slices in y so that
;   3D derivatives in the XZ plane can be taken.
;
; :Params:
;       NAME:           in, required, type=string
;                       Name of the data product to be read.
;-
function MrSim3D::Read3D, name
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if n_elements(lun) gt 0 then free_lun, lun
        if arg_present(the_error) eq 0 then void = cgErrorMsg()
        return, !Null
    endif

    
    self -> GetInfo, NX=nx, NY=ny
    self -> GetProperty, TIME=time, DIRECTORY=directory, $
                         ORIENTATION=orientation, YSLICE=yslice

    ;Check if NAME as an underscore. Replace it with a hyphen (i.e. "Pe_xx" -> "Pe-xx")
    if stregex(name, '_', /BOOLEAN) then name = strjoin(strsplit(name, '_', /EXTRACT), '-')

    ;Create file name and test for the file
    file = filepath(name + '_' + string(time, FORMAT='(i0)') + '.gda', ROOT_DIR=directory)
    if file_test(file) eq 0 then message, 'File does not exist: "' + file + '"'

    ;Open File and read data
    openr, lun, file, /GET_LUN
    
    ;Size of the data space to read
    self -> GetProperty, IXRANGE=ixrange, IZRANGE=izrange
    nz = izrange[1] - izrange[0] + 1

    ;Declare array for data
    slice3 = fltarr(nx, 3, izrange[1]-izrange[0]+1)
    buff   = fltarr(nx)

    ;Declare integer index
    n = 0L

    ;Associate buffer space with specified data file
    dfile = assoc(lun, buff)

    ;Make sure we don't request data outside of array
    j = self.yslice
    if (yslice eq 0)    then j = 1
    if (yslice eq ny-1) then j = ny-2

    ;Now loop through the 3D array in the data file and read only the
    ;needed y-slice of the data file - assume array is in FORTRAN ordering
    for k=0,nz-1 do begin
        n             = (j - 1 + ny*k)
        slice3[*,0,k] = dfile[n]
        n             = (j + ny*k)
        slice3[*,1,k] = dfile[n]
        n             = (j + 1 + ny*k)
        slice3[*,2,k] = dfile[n]
     endfor

    ;Close the file
    free_lun, lun
    
    ;Trim the data to the desired x-range
    slice3 = slice3[ixrange[0]:ixrange[1],*,*]

    ;Smooth data --- Data is already smoothed!
    if self.nSmooth gt 0 then slice3 = smooth(slice3, self.nSmooth, /EDGE_TRUNCATE)
    
    ;Return the data
    return, slice3
end


;+
;   The purpose of this program is to calculate x-component of the divergence of the
;   electron pressure tensor.
;       (\nabla \cdot P_{e})_{y} = \frac{1} {n_{e}}
;                                  \left( \frac{\partial (P_{e})_{xy}} {\partial y} +
;                                         \frac{\partial (P_{e})_{yy}} {\partial y} +
;                                         \frac{\partial (P_{e})_{zy}} {\partial y} \right)
;
; :Private:
;
; :Returns:
;       divPe_y:                The y-component of the divergence of the electron pressure
;                                   tensor.
;-
function MrSim3D::divPe_y
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
	;Pe-yx are unavailable, so I am assuming the matrix is symmetric
    Pe_yx = self -> Read3D('Pe-xy')
    Pe_yy = self -> Read3D('Pe-yy')
    Pe_yz = self -> Read3D('Pe-zy')
    n_e   = self -> getData('ne')
    
    dims = size(Pe_yy, /DIMENSIONS)
    
    ;Y-size of a simulation cell, in electron skin depths
    self -> GetInfo, DY_DE=dy_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    dPe_yx = fltarr(dims[[0,2]])
    dPe_yy = fltarr(dims[[0,2]])
    dPe_yz = fltarr(dims[[0,2]])
    divPe_y = fltarr(dims)
    
    ;Compute the derivative
    delta_Pe_xy = reform(Pe_xy[*,2,*] - Pe_xy[*,0,*]) / (2.0 * dy_de)
    delta_Pe_yy = reform(Pe_yy[*,2,*] - Pe_yy[*,0,*]) / (2.0 * dy_de)
    delta_Pe_zy = reform(Pe_zy[*,2,*] - Pe_zy[*,0,*]) / (2.0 * dy_de)
    
    ;Sum the terms, dividing by electron density
    divPe_y = (dPe_xy + dPe_yy + dPe_zy) / n_e
    
    return, divPe_y
end


;+
;   The purpose of this program is to calculate x-component of the divergence of the
;   electron pressure tensor.
;       (\nabla \cdot P_{i})_{y} = \frac{1} {n_{i}}
;                                  \left( \frac{\partial (P_{i})_{xy}} {\partial y} +
;                                         \frac{\partial (P_{i})_{yy}} {\partial y} +
;                                         \frac{\partial (P_{i})_{zy}} {\partial y} \right)
;
; :Private:
;
; :Returns:
;       divPi_y:                The y-component of the divergence of the ion pressure
;                                   tensor.
;-
function MrSim3D::divPi_y
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
	;Pi-yx are unavailable, so I am assuming the matrix is symmetric
    Pi_yx = self -> Read3D('Pi-xy')
    Pi_yy = self -> Read3D('Pi-yy')
    Pi_yz = self -> Read3D('Pi-zy')
    n_i   = self -> getData('ni')
    
    dims = size(Pi_yy, /DIMENSIONS)
    
    ;Y-size of a simulation cell, in electron skin depths
    self -> GetInfo, DY_DE=dy_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory
    dPi_yx = fltarr(dims[[0,2]])
    dPi_yy = fltarr(dims[[0,2]])
    dPi_yz = fltarr(dims[[0,2]])
    divPi_y = fltarr(dims)
    
    ;Compute the derivative
    delta_Pi_xy = reform(Pi_xy[*,2,*] - Pi_xy[*,0,*]) / (2.0 * dy_de)
    delta_Pi_yy = reform(Pi_yy[*,2,*] - Pi_yy[*,0,*]) / (2.0 * dy_de)
    delta_Pi_zy = reform(Pi_zy[*,2,*] - Pi_zy[*,0,*]) / (2.0 * dy_de)
    
    ;Sum the terms, dividing by electron density
    divPi_y = (dPi_xy + dPi_yy + dPi_zy) / n_i
    
    return, divPi_y
end


;+
;   The purpose of this program is to calculate the gradient of the scalar pressure.
;       \nabla p_{e} = \frac{1} {n_{e}}
;                                  \left( \frac{\partial p_{e}} {\partial x} +
;                                         \frac{\partial p_{e}} {\partial y} +
;                                         \frac{\partial p_{e}} {\partial z} \right)
;
; :Private:
;
; :Returns:
;       gradPe:         The x-component of the divergence of the electron pressure
;                           tensor.
;-
function MrSim3D::gradPe
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
	;Pe-yx, Pe-zx are unavailable, so I am assuming the matrix is symmetric
    Pe  = self -> GetData('Pe')
    n_e = self -> GetData('ne')
    
    dims = size(Pe, /DIMENSIONS)
    
    ;Get the size of each cell, in electron skin depths    
    self -> GetInfo, DX_DE=dx, DY_DE=dy, DZ_DE=dz
;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------
    ;Allocate memory (centered difference will have two fewer points)
    gradPe_x = fltarr(dims[0]-2, dims[2]-2)
    gradPe_y = fltarr(dims[0]-2, dims[2]-2)
    gradPe_z = fltarr(dims[0]-2, dims[2]-2)
    gradPe   = fltarr(dims[0]-2, dims[1]-2)
    
    ;Compute centered difference
    gradPe_x = reform((Pe[2:dims[0]-1,1,1:dims[2]-2] - Pe[0:dims[0]-3,1,1:dims[2]-2])) / (2.0*dx)
    gradPe_y = reform((Pe[1:dims[0]-2,2,1:dims[2]-2] - Pe[1:dims[0]-2,0,1:dims[2]-2])) / (2.0*dy)
    gradPe_z = reform((Pe[1:dims[0]-2,1,2:dims[2]-1] - Pe[1:dims[0]-2,1,0:dims[2]-3])) / (2.0*dz)
    Pe = 0B
    
    ;Sum the terms, dividing by electron density
    gradPe = sqrt((temporary(gradPe_x^2) + temporary(gradPe_y^2) + temporary(gradPe_z)^2) / (3.0*n_e^2))
    
    return, gradPe
end


;+
;   The purpose of this program is to calculate the x-component of the gradient of the
;   scalar pressure.
;       \left( \nabla p_{e} \right)_{x} = \frac{1} {n_{e}}
;                                         \left( \frac{\partial p_{e}} {\partial x} \right)
;
; :Private:
;
; :Returns:
;       gradPe_x:         The x-component of the divergence of the electron pressure
;                           tensor.
;-
function MrSim3D::gradPe_x
	compile_opt strictarr, hidden
	on_error, 2
	
	;Get the electron density and scalar pressure tensor
    Pe  = self -> GetData('Pe')
    n_e = self -> GetData('ne')
    dims = size(Pe, /DIMENSIONS)
    
    ;Get the size of each cell, in electron skin depths    
    self -> GetInfo, DX_DE=dx
    
    ;Compute centered difference
    gradPe_x = reform((Pe[2:dims[0]-1,1,1:dims[2]-2] - Pe[0:dims[0]-3,1,1:dims[2]-2])) / (2.0*dx*n_e)
    Pe = 0B
    
    return, gradPe_x
end


;+
;   The purpose of this program is to calculate the x-component of the gradient of the
;   scalar pressure.
;       \left( \nabla p_{e} \right)_{y} = \frac{1} {n_{e}}
;                                         \left( \frac{\partial p_{e}} {\partial y} \right)
;
; :Private:
;
; :Returns:
;       gradPe_y:         The y-component of the divergence of the electron pressure
;                           tensor.
;-
function MrSim3D::gradPe_y
	compile_opt strictarr, hidden
	on_error, 2
	
	;Get the electron density and scalar pressure tensor
    Pe  = self -> GetData('Pe')
    n_e = self -> GetData('ne')
    dims = size(Pe, /DIMENSIONS)
    
    ;Get the size of each cell, in electron skin depths    
    self -> GetInfo, DY_DE=dy
    
    ;Compute centered difference
    gradPe_y = reform((Pe[1:dims[0]-2,2,1:dims[2]-2] - Pe[1:dims[0]-2,0,1:dims[2]-2])) / (2.0*dy*n_e)
    Pe = 0B
    
    return, gradPe_y
end


;+
;   The purpose of this program is to calculate the z-component of the gradient of the
;   scalar pressure.
;       \left( \nabla p_{e} \right)_{z} = \frac{1} {n_{e}}
;                                         \left( \frac{\partial p_{e}} {\partial z} \right)
;
; :Private:
;
; :Returns:
;       gradPe_z:         The z-component of the divergence of the electron pressure
;                           tensor.
;-
function MrSim3D::gradPe_z
	compile_opt strictarr, hidden
	on_error, 2
	
	;Get the electron density and scalar pressure tensor
    Pe  = self -> GetData('Pe')
    n_e = self -> GetData('ne')
    dims = size(Pe, /DIMENSIONS)
    
    ;Get the size of each cell, in electron skin depths    
    self -> GetInfo, DZ_DE=dz
    
    ;Compute centered difference
    gradPe_z = reform((Pe[1:dims[0]-2,1,2:dims[2]-1] - Pe[1:dims[0]-2,1,0:dims[2]-3])) / (2.0*dz*n_e)
    Pe = 0B

    return, gradPe_z
end


;+
;   The purpose of this program is to calculate scalar pressure.
;       p_{e} = 1/3 * P_{e,xx} + P_{e,yy} + P_{e,zz}
;
; :Private:
;
; :Returns:
;       PE:         Electron scalar pressure.
;-
function MrSim3D::Pe
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
	;Pe-yx, Pe-zx are unavailable, so I am assuming the matrix is symmetric
    Pe_xx = self -> Read3D('Pe_xx')
    Pe_yy = self -> Read3D('Pe_yy')
    Pe_zz = self -> Read3D('Pe_zz')
    
    ;Calculate the scalar pressure
    Pe = (temporary(Pe_xx) + temporary(Pe_yy) + temporary(Pe_zz)) / 3.0
    
    return, Pe
end


;+
;   The purpose of this program is to calculate scalar pressure.
;       p_{i} = 1/3 * P_{i,xx} + P_{i,yy} + P_{i,zz}
;
; :Private:
;
; :Returns:
;       Pi:         Ion scalar pressure
;-
function MrSim3D::Pi
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
	;Pe-yx, Pe-zx are unavailable, so I am assuming the matrix is symmetric
    Pi_xx = self -> Read3D('Pi_xx')
    Pi_yy = self -> Read3D('Pi_yy')
    Pi_zz = self -> Read3D('Pi_zz')
    
    ;Calculate the scalar pressure
    Pi = (temporary(Pi_xx) + temporary(Pi_yy) + temporary(Pi_zz)) / 3.0
    
    return, Pi
end


;+
;   This method initializes the MrReadSim class.
;-
function MrSim3D::Init, theSim, time, yslice, $
BINARY=binary, $
ORIENTATION=orientation, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMsg()
        return, 0
    endif

    ;Superclass
    success = self -> MrSim::Init(theSim, time, yslice, BINARY=binary, _STRICT_EXTRA=extra)
    if success eq 0 then message, 'Unable to initialize the MrSim superclass.'

    ;
    ;From Bill Daughton:
    ;   The 3D data is so big, it is difficult to look at and/or move.
    ;   To make things easier, we usually apply a smoothing operator and then 
    ;   downsample by a factor of 2x in each direction - which gives a 
    ;   reduction of 8x in the file size.  This is why nx,ny,nz are 2x 
    ;   smaller in the data I sent. For example, in 1D the smoothing 
    ;   operator would be
    ;
    ;   B_j  =  0.25 *[ B(j-1) + 2*B(j) + B(j+1) ]
    ;
    ;   but instead of keeping every point, I would skip every other one for
    ;   the new "smoothed" array. 
    ;

    ;Report the actual simulation size
    self -> GetInfo, NX=nx, NY=ny, NZ=nz, LX_DE=lx_de, LY_DE=ly_de, LZ_DE=lz_de, MI_ME=mi_me
    message, string(FORMAT='(%"Simulation Size = %i x %i x %i")', nx, ny, nz), /INFORMATIONAL
    message, string(FORMAT='(%"Sim Size (de)   = %i x %i x %i")', lx_de, ly_de, lz_de), /INFORMATIONAL
    message, string(FORMAT='(%"m_i / m_e       = %f")', mi_me), /INFORMATIONAL
    
    ;But use the reduced size found in the binary file.
    if binary eq 0 then begin
        self -> ReadInfo_Binary, self.directory + '/info'
        self -> MakeSimDomain
    endif
    
    return, 1
end


;+
;   This is the class definition.
;-
pro MrSim3D__define, class
    class = { MrSim3D, $
              inherits MrSim $
            }
end