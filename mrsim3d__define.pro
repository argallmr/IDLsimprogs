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
;       2014/10/11  -   Oranization, comments, and flow to ::ReadMoments improved. - MRA
;       2014/10/14  -   ::ReadElectrons is now a function. File checking delegated
;                           elsewhere. - MRA
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
function MrSim3D::ReadGDA_FilePos, coords, $
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
        ix = getIndexRange(*self.xSim, self.xrange)
        iy = reform(getIndexRange(*self.ySim, self.yrange))
        iz = getIndexRange(*self.zSim, self.zrange)
    endif else begin
        ix = value_locate(*self.xSim, coords[0,*])
        iy = value_locate(*self.ySim, coords[1,*])
        iz = value_locate(*self.zSim, coords[2,*])
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
        ; The x- and z-axis need to be interchanged.
        ;   - The reverse occurs here.
        ;
        ;           MSP                         MSP
        ; +z  *--------------|       ++x *--------------|
        ;     |              |  ==>      |              |
        ; -z  |--------------+        +x |--------------+
        ;     +x    MSH    ++x           -z     MSH    +z
        ;
        ; To complete the transformation
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
        ; And finally
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
pro MrSim3D::ReadGDA_TransformData, name, data, $
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

    _name   = strupcase(name)
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
                stregex(_name, '[A-Z][a-z]*z$',   /BOOLEAN): data = -data
                stregex(_name, '-(z[yx]|[xy]z)$', /BOOLEAN): data = -data
                else: ;Do nothing
            endcase
        end
        
        'MAGNETOTAIL': begin
            case _orient of
                'XY': ;Do nothing
                'XZ': ;Do nothing
                else: message, 'Orienatation "' + _orient + '" not recognized.'
            endcase

            ;Now reverse the x-axes
            ;   Single negate Bx, Uix, Pe-xz, Pe-xy, etc.
            case 1 of
                stregex(_name, '[A-Z][a-z]*(x)$', /BOOLEAN): data = -data
                stregex(_name, '-(y|z)[x]$',      /BOOLEAN): data = -data
                stregex(_name, '-[x](y|z)$',      /BOOLEAN): data = -data
                else: ;Do nothing
            endcase
        end
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
function MrSim3D::ReadGDA_Cells, cells, name, tIndex, $
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
function MrSim3D::ReadGDA_XY, file, ixrange, iyrange, iz, $
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
function MrSim3D::ReadGDA_XZ, file, ixrange, izrange, iy, $
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

    ;How to read:
    ;   - Y is fixed
    ;   - Jump from z to z
    ;   - Read all columns (x) at each z
    ix = ulong64(ixMin)
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
function MrSim3D::ReadGDA_YZ, file, iyrange, izrange, ix, $
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
function MrSim3D::ReadGDA_XYZ, file, ixrange, iyrange, izrange, $
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
    ;   - When a full xy-plane is read, advance depth, repeat.
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
function MrSim3D::ReadGDA, name, tIndex, $
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
    if n_elements(directory)   eq 0 then directory   = self.directory
    if n_elements(nSmooth)     eq 0 then nSmooth     = self.nSmooth
    if n_elements(orientation) eq 0 then orientation = self.orientation
    if n_elements(tindex)      eq 0 then tindex      = self.time
    if n_elements(xrange)      eq 0 then xrange      = self.xrange
    if n_elements(yrange)      eq 0 then yrange      = self.yrange
    if n_elements(zrange)      eq 0 then zrange      = self.zrange
    coord_system = n_elements(coord_system) eq 0 ? self.coord_system : strupcase(coord_system)

    ;Create file name and test for the file
    file = filepath(name + '_' + string(tindex, FORMAT='(i0)') + '.gda', ROOT_DIR=directory)
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
        'XY': data = self -> ReadGDA_XY(file, cells[0,*], cells[1,*], cells[2,0])
        'XZ': data = self -> ReadGDA_XZ(file, cells[0,*], cells[2,*], cells[1,0])
        'YZ': data = self -> ReadGDA_YZ(file, cells[1,*], cells[2,*], cells[0,0])
        else: message, 'Orientation "' + orientation + '" not recognized.'
    endcase

    ;Manipulate data from SIMULATION to COORD_SYSTEM coordinates.
    self -> ReadGDA_TransformData, name, data

;-------------------------------------------------------
; Cleanup //////////////////////////////////////////////
;-------------------------------------------------------
    ;Smooth data 
    if nSmooth gt 0 then data = smooth(data, nSmooth, /EDGE_TRUNCATE)  

    ;Store the data so that it does not have to be read again.
    return, data
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
function MrSim3D::GridLine, r0, r1, $
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
    self -> GetInfo, DX_DE=dx_de, DY_DE=dy_de, DZ_DE=dz_de
    
    ;Number of grid cells between r0 and r1
    dr = sqrt( total( (r1 - r0)^2 ) )
    n  = dr * dx_de
    
    ;Coordinates of each point on the line.
    coords = MrLine3D(r0, r1, NPOINTS=n)
    
    ;Find the cells
    cells = coords
    cells[0,*] = self -> GetCell(coords[0,*], /X)
    cells[1,*] = self -> GetCell(coords[1,*], /Y)
    cells[2,*] = self -> GetCell(coords[2,*], /Z)

    return, cells
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
function MrSim3D::ReadElectrons, filename, $
ENERGY=energy, $
FMAP_DIR=fmap_dir, $
VELOCITY=velocity, $
VERBOSE=verbose, $
XRANGE=xrange, $
ZRANGE=zrange
    compile_opt strictarr
    
    ;catch errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        void = cgErrorMsg()
        return, !Null
    endif

    ;Get the y-grid-cell
    yslice = self -> GetCell(self.yrange[0], /Y)

    ;Create defaults    
    MrSim_Which, self.simname, EFILE=filename, FMAP_DIR=fMap_dir, TINDEX=self.time, YSLICE=yslice
    
    ;Set defaults
    if n_elements(filename) eq 0 then filename = eFile
    if n_elements(fmap_dir) eq 0 then fmap_dir = fmap
    if n_elements(xrange)   eq 0 then xrange   = self.xrange
    if n_elements(zrange)   eq 0 then zrange   = self.zrange
    
    ;Read the simulation info file
    data = MrSim_ReadParticles(filename, xrange, zrange, $
                               ENERGY           = energy, $
                               FMAP_DIR         = fMap_dir, $
                               VELOCITY         = velocity, $
                               VERBOSE          = verbose)
    
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
function MrSim3D::ReadMoment, name, tIndex, $
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
    if n_elements(xrange)      eq 0 then xrange      = self.xrange
    if n_elements(yrange)      eq 0 then yrange      = self.yrange
    if n_elements(zrange)      eq 0 then zrange      = self.zrange
    coord_system = n_elements(coord_system) eq 0 ? self.coord_system : strupcase(coord_system)

    ;Create file name and test for the file
    file = filepath(name + '_' + string(tindex, FORMAT='(i0)') + '.gda', ROOT_DIR=directory)
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
        'XY': data = self -> ReadGDA_XY(file, cells[0,*], cells[1,*], cells[2,0])
        'XZ': data = self -> ReadGDA_XZ(file, cells[0,*], cells[2,*], cells[1,0])
        'YZ': data = self -> ReadGDA_YZ(file, cells[1,*], cells[2,*], cells[0,0])
        else: message, 'Orientation "' + orientation + '" not recognized.'
    endcase

    ;Manipulate data from SIMULATION to COORD_SYSTEM coordinates.
    self -> ReadGDA_TransformData, name, data


;-------------------------------------------------------
; Cleanup //////////////////////////////////////////////
;-------------------------------------------------------
    ;Smooth data 
    if nSmooth gt 0 then data = smooth(data, nSmooth, /EDGE_TRUNCATE)  

    ;Store the data so that it does not have to be read again.
    return, data
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
    data = self -> ReadGDA(name)
    
    ;Store the data so the file does not have to be accessed again.
    self -> SetData, name, temporary(data)
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
pro MrSim3D::SetSim, sim_id, $
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
    if dimension ne '2D' then message, 'Cannot change dimensions from 3D to ' + dimension + '.'
    
    ;Call the superclass
    self -> MrSim::SetSim, sim_id, BINARY=binary, DIRECTORY=directory, INFO_FILE=info_file
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