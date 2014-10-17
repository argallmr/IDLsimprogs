; docformat = 'rst'
;
;+
;   The purpose of this program is to provide a base class for 2D and 3D simulation
;   data provided by Bill Daughton from Los Alamos National Laboratory.
;
;   MrSim2_Data is meant to be subclassed by a class that over-rides the ReadData method.
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
;       Copyright 2014 by the University of New Hampshire
;
; :History:
;   Modification History::
;
;       2014/09/06  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   This method initializes the MrSim2_Data class.
;-
function MrSim2_Data::INIT
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, 0
    endif

;-------------------------------------------------------
;Allocate Heap /////////////////////////////////////////
;-------------------------------------------------------

    ;Data Products
    self.electrons = ptr_new(/ALLOCATE_HEAP)
    self.Ay = ptr_new(/ALLOCATE_HEAP)
    self.Bx = ptr_new(/ALLOCATE_HEAP)
    self.By = ptr_new(/ALLOCATE_HEAP)
    self.Bz = ptr_new(/ALLOCATE_HEAP)
    self.Ex = ptr_new(/ALLOCATE_HEAP)
    self.Ey = ptr_new(/ALLOCATE_HEAP)
    self.Ez = ptr_new(/ALLOCATE_HEAP)
    self.n_e = ptr_new(/ALLOCATE_HEAP)
    self.n_i = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xx = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xy = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xz = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yx = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yy = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yz = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zx = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zy = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zz = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xx = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xy = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xz = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yx = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yy = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yz = ptr_new(/ALLOCATE_HEAP)
    self.Pi_zx = ptr_new(/ALLOCATE_HEAP)
    self.Pi_zy = ptr_new(/ALLOCATE_HEAP)
    self.Pi_zz = ptr_new(/ALLOCATE_HEAP)
    self.Uex = ptr_new(/ALLOCATE_HEAP)
    self.Uey = ptr_new(/ALLOCATE_HEAP)
    self.Uez = ptr_new(/ALLOCATE_HEAP)
    self.Uix = ptr_new(/ALLOCATE_HEAP)
    self.Uiy = ptr_new(/ALLOCATE_HEAP)
    self.Uiz = ptr_new(/ALLOCATE_HEAP)

    return, 1
end


;+
;   Clean up after the object is destroyed.
;-
pro MrSim2_Data::CLEANUP
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Data products
    ptr_free, self.Ay
    ptr_free, self.Bx
    ptr_free, self.By
    ptr_free, self.Bz
    ptr_free, self.electrons
    ptr_free, self.Ex
    ptr_free, self.Ey
    ptr_free, self.Ez
    ptr_free, self.n_e
    ptr_free, self.n_i
    ptr_free, self.Pe_xx
    ptr_free, self.Pe_xy
    ptr_free, self.Pe_xz
    ptr_free, self.Pe_yy
    ptr_free, self.Pe_yz
    ptr_free, self.Pe_zz
    ptr_free, self.Pi_xx
    ptr_free, self.Pi_xy
    ptr_free, self.Pi_xz
    ptr_free, self.Pi_yy
    ptr_free, self.Pi_yz
    ptr_free, self.Pi_zz
    ptr_free, self.Uex
    ptr_free, self.Uey
    ptr_free, self.Uez
    ptr_free, self.Uix
    ptr_free, self.Uiy
    ptr_free, self.Uiz
end



;+
;   Take the cross product of two quantities.
;
; :Params:
;       V1:             in, required, type=NxMx3 float or string
;                       A vectory quantity or the name of the vectory quantity to
;                           be crossed into `V2`.
;       V2:             in, required, type=NxMx3 float or string
;                       A vectory quantity or the name of the vectory quantity to
;                           be crossed with `V1`.
;
; :Keywords:
;       MAGNITUDE:      in, optional, type=boolean, default=0
;                       If set, the magnitude will be returned.
;       X:              in, optional, type=boolean, default=0
;                       If set, the x-component will be returned.
;       Y:              in, optional, type=boolean, default=0
;                       If set, the y-component will be returned.
;       Z:              in, optional, type=boolean, default=0
;                       If set, the z-component will be returned.
;
; :Returns:
;       V1XV2:          Cross product of `V1` and `V2`. If no keywords are set, an
;                           NxMx3 array will be retured, where the last demension
;                           corresponds to the x-, y-, and z-components.
;-
function MrSim2::CrossProduct, v1, v2, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_v1) then ptr_free, _v1
        if ptr_valid(_v2) then ptr_free, _v2
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Was a name or data given?
    if MrIsA(v1, 'STRING') $
        then _v1 = ptr_new(self -> GetData(v1)) $
        else _v1 = ptr_new(v1)
        
    if MrIsA(v2, 'STRING') $
        then _v2 = ptr_new(self -> GetData(v2)) $
        else _v2 = ptr_new(v2)
        
    ;Compute the cross product
    v1xv2 = [[[(*_v1)[*,*,1] * (*_v2)[*,*,2] - (*_v1)[*,*,2] * (*_v2)[*,*,1]]], $
             [[(*_v1)[*,*,2] * (*_v2)[*,*,0] - (*_v1)[*,*,0] * (*_v2)[*,*,2]]], $
             [[(*_v1)[*,*,0] * (*_v2)[*,*,1] - (*_v1)[*,*,1] * (*_v2)[*,*,0]]]]
    
    ;Free the pointers
    ptr_free, _v1
    ptr_free, _v2
    
    ;return
    case 1 of
        x:         return, v1xv2[*,*,0]
        y:         return, v1xv2[*,*,1]
        z:         return, v1xv2[*,*,2]
        magnitude: return, sqrt(total(v1xv2^2, 3))
        else:      return, v1xv2
    endcase
end


;+
;   Take the dot product of two quantities.
;
; :Params:
;       V1:             in, required, type=NxMx3 float or string
;                       A vectory quantity or the name of the vectory quantity to
;                           be dotted into `V2`.
;       V2:             in, required, type=NxMx3 float or string
;                       A vectory quantity or the name of the vectory quantity to
;                           be dotted with `V1`.
;
; :Keywords:
;       X:              in, optional, type=boolean, default=0
;                       If set, only the x-components will be used.
;       Y:              in, optional, type=boolean, default=0
;                       If set, only the y-components will be used.
;       Z:              in, optional, type=boolean, default=0
;                       If set, only the z-components will be used.
;
; :Returns:
;       V1DOTV2:        Dot product of `V1` and `V2`.
;-
function MrSim2::DotProduct, v1, v2, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_v1) then ptr_free, _v1
        if ptr_valid(_v2) then ptr_free, _v2
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Was a name or data given?
    if MrIsA(v1, 'STRING') $
        then _v1 = ptr_new(self -> GetData(v1)) $
        else _v1 = ptr_new(v1)
        
    if MrIsA(v2, 'STRING') $
        then _v2 = ptr_new(self -> GetData(v2)) $
        else _v2 = ptr_new(v2)
    
    ;return
    case 1 of
        x:    v1dotv2 =   (*_v1)[*,*,0] * (*_v2)[*,*,0]
        y:    v1dotv2 =   (*_v1)[*,*,1] * (*_v2)[*,*,1]
        z:    v1dotv2 =   (*_v1)[*,*,2] * (*_v2)[*,*,2]
        else: v1dotv2 = ( (*_v1)[*,*,0] * (*_v2)[*,*,0] + $
                          (*_v1)[*,*,1] * (*_v2)[*,*,1] + $
                          (*_v1)[*,*,2] * (*_v2)[*,*,2] )
    endcase
        
    ;Free the pointers
    ptr_free, _v1
    ptr_free, _v2
    
    return, v1dotv2
end


;+
;   Compute the x-derivative of a data array.
;
; :Params:
;       DATA:           in, required, type=string/NxM float
;                       A string naming the data product or the actual data of which
;                           the derivative is to be taken.
;
; :Returns:
;       D_DX:           The derivative of `DATA` with respect to X
;-
function MrSim2::D_DX, data
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_data) then ptr_free, _data
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get the data
    ;   - Avoid copying the data by creating a pointer
    if size(data, /TNAME) eq 'STRING' $
        then _data = ptr_new(self -> GetData(data)) $
        else _data = ptr_new(data)
    
    ;Get the x-size of the simulation in electron skin depths
    self -> GetInfo, DX_DE=dx_de
    dims = size(data, /DIMENSIONS)
    if n_elements(dims) ne 2 then $
        message, 'Data must be 2D in order to take the derivative.'
    
    ;Determine the orientation
    ;   'XY' -> data[X,Y]
    ;   'XZ' -> data[X,Z]
    ;   etc.
    xaxis = strmid(self.orientation, 0, 1)
    yaxis = strmid(self.orientation, 1, 2)        
    
    ;Allocate memory
    d_dx = fltarr(dims)

    ;Take the centered difference.
    case 'X' of
        xaxis: d_dx[1:dims[0]-2,*] = ((*_data)[2:dims[0]-1,*] - (*_data)[0:dims[0]-3,*]) / (2.0 * dx_de)
        yaxis: d_dx[*,1:dims[1]-2] = ((*_data)[*,2:dims[1]-1] - (*_data)[*,0:dims[1]-3]) / (2.0 * dx_de)
        else:  message, 'Orientation "' + self.orientation + '" does not allow the derivative with respect to x.'
    endcase
    
    ptr_free, _data
    return, d_dx
end


;+
;   Compute the y-derivative of a data array.
;
; :Params:
;       DATA:           in, required, type=string/NxM float
;                       A string naming the data product or the actual data of which
;                           the derivative is to be taken.
;
; :Returns:
;       D_DY:           The derivative of `DATA` with respect to Y.
;-
function MrSim2::D_DY, data
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_data) then ptr_free, _data
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get the data
    ;   - Avoid copying the data by creating a pointer
    if size(data, /TNAME) eq 'STRING' $
        then _data = ptr_new(self -> GetData(data)) $
        else _data = ptr_new(data)
    
    ;Get the x-size of the simulation in electron skin depths
    self -> GetInfo, DY_DE=dy_de
    dims = size(data, /DIMENSIONS)
    if n_elements(dims) ne 2 then $
        message, 'Data must be 2D in order to take the derivative.'
    
    ;Determine the orientation
    ;   'YZ' -> data[Y,Z]
    ;   'XY' -> data[X,Y]
    ;   etc.
    xaxis = strmid(self.orientation, 0, 1)
    yaxis = strmid(self.orientation, 1, 2)        
    
    ;Allocate memory
    d_dy = fltarr(dims)

    ;Take the centered difference.
    case 'Y' of
        xaxis: d_dy[1:dims[0]-2,*] = ((*_data)[2:dims[0]-1,*] - (*_data)[0:dims[0]-3,*]) / (2.0 * dy_de)
        yaxis: d_dy[*,1:dims[1]-2] = ((*_data)[*,2:dims[1]-1] - (*_data)[*,0:dims[1]-3]) / (2.0 * dy_de)
        else:  message, 'Orientation "' + self.orientation + '" does not allow the derivative with respect to x.'
    endcase
    
    ptr_free, _data
    return, d_dy
end


;+
;   Compute the z-derivative of a data array.
;
; :Params:
;       DATA:           in, required, type=string/NxM float
;                       A string naming the data product or the actual data of which
;                           the derivative is to be taken.
;
; :Returns:
;       D_DZ:           The derivative of `DATA` with respect to Z.
;-
function MrSim2::D_DZ, data
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_data) then ptr_free, _data
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Get the data
    ;   - Avoid copying the data by creating a pointer
    if size(data, /TNAME) eq 'STRING' $
        then _data = ptr_new(self -> GetData(data)) $
        else _data = ptr_new(data)
    
    ;Get the x-size of the simulation in electron skin depths
    self -> GetInfo, DZ_DE=dz_de
    dims = size(data, /DIMENSIONS)
    if n_elements(dims) ne 2 then $
        message, 'Data must be 2D in order to take the derivative.'
    
    ;Determine the orientation
    ;   'YZ' -> data[Y,Z]
    ;   'XZ' -> data[X,Z]
    ;   etc.
    xaxis = strmid(self.orientation, 0, 1)
    yaxis = strmid(self.orientation, 1, 2)        
    
    ;Allocate memory
    d_dz = fltarr(dims)

    ;Take the centered difference.
    case 'Z' of
        xaxis: d_dz[1:dims[0]-2,*] = ((*_data)[2:dims[0]-1,*] - (*_data)[0:dims[0]-3,*]) / (2.0 * dz_de)
        yaxis: d_dz[*,1:dims[1]-2] = ((*_data)[*,2:dims[1]-1] - (*_data)[*,0:dims[1]-3]) / (2.0 * dz_de)
        else:  message, 'Orientation "' + self.orientation + '" does not allow the derivative with respect to x.'
    endcase
    
    ptr_free, _data
    return, d_dz
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
function MrSim2::GetData, name, $
DX=dx, $
DY=dy, $
DZ=dz, $
CROSS=cross, $
DOT=dot, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
TENSOR=tensor, $
TXX=txx, $
TXY=txy, $
TXZ=Txz, $
TYY=Tyy, $
TYZ=Tyz, $
TZZ=Tzz, $
VECTOR=vec, $
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

    ;If no parameters were given, print a list of available data products.
    if n_params() eq 0 then begin
        MrSim2 -> listProducts
        return, !Null
    endif

    ;Defaults
    d_dx         = keyword_set(d_dx)
    d_dy         = keyword_set(d_dy)
    d_dz         = keyword_set(d_dz)
    magnitude    = keyword_set(magnitude)
    parallel     = keyword_set(parallel)
    perpendicuar = keyword_set(perpendicular)
    pointer      = keyword_set(pointer)
    Txx          = keyword_set(Txx)
    Txy          = keyword_set(Txy)
    Txz          = keyword_set(Txz)
    Tyy          = keyword_set(Tyy)
    Tyz          = keyword_set(Tyz)
    Tzz          = keyword_set(Tzz)
    x            = keyword_set(x)
    y            = keyword_set(y)
    z            = keyword_set(z)
    if (magnitude + x + y + z)             gt 1 then message, 'MAGNITUDE, X, Y, Z are mutually exclusive.'
    if (parallel + perpendicular)          gt 1 then message, 'PARALLEL and PERPENDICULAR are mutually exclusive.'
    if (d_dx + d_dz)                       gt 1 then message, 'D_DX and D_DZ are mutually exclusive.'
    if (Txx + Txy + Txz + Tyy + Tyz + Tzz) gt 1 then message, 'TXX, TXY, TXZ, TYY, TYZ, and TZZ are mutually exclusive.'

;-------------------------------------------------------
; Parse the NAME String ////////////////////////////////
;-------------------------------------------------------
    ;Operations are separated by spaces
    parts = strsplit(name, ' ', /EXTRACT, COUNT=nParts)
    _name = strupcase(parts[0])
    
    ;Perform an operation on the data?
    if nParts eq 3 then begin
        case strupcase(parts[1]) of
            'DOT':           data = self -> DotProduct(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'CROSS':         data = self -> CrossProduct(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'X':             data = self -> CrossProduct(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            '||':            data = self -> Parallel(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'PAR':           data = self -> Parallel(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'PARALLEL':      data = self -> Parallel(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'PERP':          data = self -> Perpendicular(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'PERPENDICULAR': data = self -> Perpendicular(parts[0], parts[2], X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            else: message, 'Operation not recognized: "' + parts[1] + '".'
        endcase
        
        return, data
    endif


;-------------------------------------------------------
;Had Data been Read? ///////////////////////////////////
;-------------------------------------------------------
    
    ;If the data has not been read, then try to read it from a file.
    if self -> HasData(name) eq 0 then self -> ReadData, name
    
    ;Check to see if the data has already been read first.
    case strupcase(data_product) of
        'AY':     data = self -> A(/Y)
        'BX':     data = self -> B(/X)
        'BY':     data = self -> B(/Y)
        'BZ':     data = self -> B(/Z)
        'EX':     data = self -> E(/X)
        'EY':     data = self -> E(/Y)
        'EZ':     data = self -> E(/Z)
        'NE':     data = self -> n_e()
        'NI':     data = self -> n_i()
        'PE-XX':  data = self -> Pe(/TXX)
        'PE-XY':  data = self -> Pe(/TXY)
        'PE-XZ':  data = self -> Pe(/TXZ)
        'PE-YX':  data = self -> Pe(/TXY)
        'PE-YY':  data = self -> Pe(/TYY)
        'PE-YZ':  data = self -> Pe(/TYZ)
        'PE-ZX':  data = self -> Pe(/TXZ)
        'PE-ZY':  data = self -> Pe(/TYZ)
        'PE-ZZ':  data = self -> Pe(/TZZ)
        'PI-XX':  data = self -> Pi(/TXX)
        'PI-XY':  data = self -> Pi(/TXY)
        'PI-XZ':  data = self -> Pi(/TXZ)
        'PI-YX':  data = self -> Pi(/TXY)
        'PI-YY':  data = self -> Pi(/TYY)
        'PI-YZ':  data = self -> Pi(/TYZ)
        'PI-ZX':  data = self -> Pi(/TXZ)
        'PI-ZY':  data = self -> Pi(/TYZ)
        'PI-ZZ':  data = self -> Pi(/TZZ)
        'UEX':    data = self -> Ue(/X)
        'UEY':    data = self -> Ue(/Y)
        'UEZ':    data = self -> Ue(/Z)
        'UIX':    data = self -> Ui(/X)
        'UIY':    data = self -> Ui(/Y)
        'UIZ':    data = self -> Ui(/Z)
        
        ;Custom Data Products
        'B':     data = self ->  B(X=x, Y=y, Z=z, MAG=magnitude, PAR=parallel, PERP=perpendicular)
        'E':     data = self ->  E(X=x, Y=y, Z=z, MAG=magnitude, PAR=parallel, PERP=perpendicular)
        'UE':    data = self -> Ue(X=x, Y=y, Z=z, MAG=magnitude, PAR=parallel, PERP=perpendicular)
        'UI':    data = self -> Ui(X=x, Y=y, Z=z, MAG=magnitude, PAR=parallel, PERP=perpendicular)
        'PE':    data = self -> Pe(TXX=Txx, TXY=Txy, TXZ=Txz, TYY=Tyy, TYZ=Tyz, TZZ=Tzz, MAG=magnitude, PAR=parallel, PERP=perpendicular)
        'PI':    data = self -> Pi(TXX=Txx, TXY=Txy, TXZ=Txz, TYY=Tyy, TYZ=Tyz, TZZ=Tzz, MAG=magnitude, PAR=parallel, PERP=perpendicular)
        'DNG_E': data = self -> Dng_e()
        else: message, 'Data product not available: "' + data_product + '".'
    endcase
    
    ;Take the derivative?
    if dx then data = self -> D_DX(data)
    if dy then data = self -> D_DY(data)
    if dz then data = self -> D_DZ(data)
    
    return, data
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
function MrSim2::HasData, data_product
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
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
        'PE-YX': if n_elements(*self.Pe_yx)     gt 0 then tf_has = 1B
        'PE-YY': if n_elements(*self.Pe_yy)     gt 0 then tf_has = 1B
        'PE-YZ': if n_elements(*self.Pe_yz)     gt 0 then tf_has = 1B
        'PE-ZX': if n_elements(*self.Pe_zx)     gt 0 then tf_has = 1B
        'PE-ZY': if n_elements(*self.Pe_zy)     gt 0 then tf_has = 1B
        'PE-ZZ': if n_elements(*self.Pe_zz)     gt 0 then tf_has = 1B
        'PI-XX': if n_elements(*self.Pi_xx)     gt 0 then tf_has = 1B
        'PI-XY': if n_elements(*self.Pi_xy)     gt 0 then tf_has = 1B
        'PI-XZ': if n_elements(*self.Pi_xz)     gt 0 then tf_has = 1B
        'PI-YX': if n_elements(*self.Pi_yx)     gt 0 then tf_has = 1B
        'PI-YY': if n_elements(*self.Pi_yy)     gt 0 then tf_has = 1B
        'PI-YZ': if n_elements(*self.Pi_yz)     gt 0 then tf_has = 1B
        'PI-ZX': if n_elements(*self.Pi_zx)     gt 0 then tf_has = 1B
        'PI-ZY': if n_elements(*self.Pi_zy)     gt 0 then tf_has = 1B
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
;   Compute the magnitude of a vector.
;
;   Calling Sequence:
;       data = oSim -> Magnitude('E')
;       data = oSim -> Magnitude(Ex, Ey, Ez)
;
; :Params:
;       VX:             in, required, type=string,fltarr(N\,M)
;                       Either the name or x-component of a vector quantity for which
;                           the magnitude is to be found.
;       VY:             in, optional, type=fltarr(N\,M)
;                       If `VX` is an array, then VY represents the y-component of the
;                           vector quantity for which the magnitude is to be computed.
;       VZ:             in, optional, type=fltarr(N\,M)
;                       If `VX` is an array, then VZ represents the z-component of the
;                           vector quantity for which the magnitude is to be computed.
;
; :Returns:
;       VMAG:           Magnitude of the vector quantity.
;-
function MrSim2::Magnitude, vx, vy, vz
    compile_opt strictarr
    on_error, 2
    
    ;Was a name or data given?
    if MrIsA(vx, 'STRING') then begin
        names = vx
    
        ;Get the data
        _vx = self -> GetData(name, /X)
        vy  = self -> GetData(name, /Y)
        vz  = self -> GetData(name, /Z)
    endif else begin
        _vx = vx
    endelse
    
    ;Compute the magnitude
    vmag = sqrt(_vx^2 + vy^2 + vz^2)
    
    return, vmag
end


;+
;   Compute the component of a vector parallel to another vector.
;
; :Params:
;       V1:             in, required, type=NxMx3 float or string
;                       A vectory quantity or the name of the vectory quantity to
;                           for which the component parallel to `V2` is to be determined.
;       V2:             in, required, type=NxMx3 float or string
;                       A vectory quantity or the name of the vectory quantity indicating
;                           the parallel direction.
;
; :Keywords:
;       X:              in, optional, type=boolean, default=0
;                       If set, only the x-components will be used.
;       Y:              in, optional, type=boolean, default=0
;                       If set, only the y-components will be used.
;       Z:              in, optional, type=boolean, default=0
;                       If set, only the z-components will be used.
;
; :Returns:
;       V1_PAR:         The component of `V1` parallel to `V2`.
;-
function MrSim2::Parallel, v1, v2, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_v1) then ptr_free, _v1
        if ptr_valid(_v2) then ptr_free, _v2
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Was a name or data given?
    if MrIsA(v1, 'STRING') $
        then _v1 = ptr_new(self -> GetData(v1)) $
        else _v1 = ptr_new(v1)
        
    if MrIsA(v2, 'STRING') $
        then _v2 = ptr_new(self -> GetData(v2)) $
        else _v2 = ptr_new(v2)
    
    ;Get Data
    v2_mag = sqrt(total(*_v2^2, 3))
    
    ;Components
    ;   par_x = vx * fx_hat
    ;   par_y = vy * fy_hat
    ;   par_z = vz * fz_hat
    case 1 of
        x:    v_par =   (*_v1)[*,*,0] * (*_v2)[*,*,0] / temporary(v2_mag))
        y:    v_par =   (*_v1)[*,*,1] * (*_v2)[*,*,1] / temporary(v2_mag))
        z:    v_par =   (*_v1)[*,*,2] * (*_v2)[*,*,2] / temporary(v2_mag))
        else: v_par = ( (*_v1)[*,*,0] * (*_v2)[*,*,0] + $
                        (*_v1)[*,*,1] * (*_v2)[*,*,1] + $
                        (*_v1)[*,*,2] * (*_v2)[*,*,2] ) / temporary(v2_mag)
    endcase
    
    ;Free the pointers
    ptr_free, _v1
    ptr_free, _v2
    
    return, v_par
end


;+
;   Compute the component of a vector perpendicular to another vector.
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
function MrSim2::Perpendicular, v1, v2, $
FIELD=field
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_v1) then ptr_free, _v1
        if ptr_valid(_v2) then ptr_free, _v2
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Compute the parallel component
    v1_par  = self -> Parallel(v1, v2)
    
    ;Get the data
    ;   - Was a name or data given?
    if MrIsA(v1, 'STRING') $
        then _v1 = ptr_new(self -> GetData(v1)) $
        else _v1 = ptr_new(v1)
        
    if MrIsA(v2, 'STRING') $
        then _v2 = ptr_new(self -> GetData(v2)) $
        else _v2 = ptr_new(v2)
        
    ;Find the perpendicular component
    v1_mag2 = total(*_v1^2, 3)
    v1_perp = sqrt(v_mag2 - v_par^2)
    
    ;Free pointers
    ptr_free, _v1
    ptr_free, _v2
    
    ;Perpendicular component
    return, v_perp 
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
pro MrSim2::ReadData, name
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
pro MrSim2::SetData, name, data
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
        'PE-YX': *self.Pe_yx = data
        'PE-YY': *self.Pe_yy = data
        'PE-YZ': *self.Pe_yz = data
        'PE-ZX': *self.Pe_zx = data
        'PE-ZY': *self.Pe_zy = data
        'PE-ZZ': *self.Pe_zz = data
        'PI-XX': *self.Pi_xx = data
        'PI-XY': *self.Pi_xy = data
        'PI-XZ': *self.Pi_xz = data
        'PI-YX': *self.Pi_yx = data
        'PI-YY': *self.Pi_yy = data
        'PI-YZ': *self.Pi_yz = data
        'PI-ZX': *self.Pi_zx = data
        'PI-ZY': *self.Pi_zy = data
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
;   The purpose of this program is to calculate the electron plasma beta::
;       \beta_{e} = \frac{Tr(P_{e})} {B^{2} / {2 \mu_{0}}
;
; :Private:
;
; :Returns:
;       BETA_E:                 The electron beta.
;-
function MrSim2_Data::A, $
CROSS=cross, $
DOT=dot, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    magnitude     = keyword_set(magnitude)
    parallel      = keyword_set(parallel)
    perpendicular = keyword_set(perpendicular)
    vec           = keyword_set(vec)
    x             = keyword_set(x)
    y             = keyword_set(y)
    z             = keyword_set(z)
    
    ;Ay is the only component
    if perpendicular + parallel + x + z + magnitude + vec gt 1 $
        then message, 'Only the y-component of the vector potential is available.'
    
    ;Ay
    if n_elements(*self.Ay) eq 0 then self -> ReadData, 'Ay'
    return, *self.Ay
end


;+
;   Return various data products associated with the magnetic field.
;
; :Private:
;
; :Keywords:
;       CROSS:              in, optional, type=boolean/string, default=0/''
;                           If set, the cross product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the cross product with the quantity specified will
;                               be returned.
;       DOT:                in, optional, type=boolean/string, default=0/''
;                           If set, the dot product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the dot product with the quantity specified will
;                               be returned.
;       PARALLEL:           in, optional, type=boolean/string, default=0/''
;                           If set, the component parallel to B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component parallel to the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the magnitude is returned.
;       VECTOR:             in, optional, type=boolean, default=0
;                           If set, all three components of the electric field will be
;                               returned in an NxMx3-dimensional array. If no other
;                               keywords are set, this is the default.
;       X:                  in, optional, type=boolean, default=0
;                           If set, the X-component is returned.
;       Y:                  in, optional, type=boolean, default=0
;                           If set, the Y-component is returned.
;       Z:                  in, optional, type=boolean, default=0
;                           If set, the Z-component is returned.
;
; :Returns:
;       DATA:               Magnetic field data.
;-
function MrSim2_Data::B, $
CROSS=cross, $
DOT=dot, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
VECTOR=vec
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    tf_cross  = keyword_set(cross)
    tf_dot    = keyword_set(dot)
    magnitude = keyword_set(magnitude)
    tf_par    = keyword_set(parallel)
    tf_perp   = keyword_set(perpendicular)
    vec       = keyword_set(vec)
    x         = keyword_set(x)
    y         = keyword_set(y)
    z         = keyword_set(z)
    
    ;If nothing was chosen, return the vector.
    vec = vec || (x + y + z + tf_cross + tf_dot + tf_perp + tf_par + magnitude) eq 0
    if vec then begin
        x = 1
        y = 1
        z = 1
    endif

    ;Parallel & perpendicular to what?
    dotName   = MrIsA(dot,           'STRING') ? dot           : 'B'
    crossName = MrIsA(cross,         'STRING') ? cross         : 'B'
    parName   = MrIsA(parallel,      'STRING') ? parellel      : 'B'
    perpName  = MrIsA(perpendicular, 'STRING') ? perpendicular : 'B'
    
    ;[XYZ]-Components
    if x then if n_elements(*self.Bx) eq 0 then self -> ReadData, 'Bx'
    if y then if n_elements(*self.By) eq 0 then self -> ReadData, 'By'
    if z then if n_elements(*self.Bz) eq 0 then self -> ReadData, 'Bz'
    
    ;Return the proper quantity
    case 1 of
        magnitude: return, self -> Magnitude('B')
        tf_cross:  return, self -> CrossProduct('B', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> DotProduct('B', dotName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Parallel( 'B', parName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Perpendicular('B', perpName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, [[[*self.Bx]], [[*self.By]], [[*self.Bz]]]
        x:         return, *self.Bx
        y:         return, *self.By
        z:         return, *self.Bz
    endcase
end


;+
;   The purpose of this program is to return various data products associated with
;   the electric field.
;
; :Private:
;
; :Keywords;
;       CROSS:              in, optional, type=boolean/string, default=0/''
;                           If set, the cross product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the cross product with the quantity specified will
;                               be returned.
;       DOT:                in, optional, type=boolean/string, default=0/''
;                           If set, the dot product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the dot product with the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the magnitude of E will be returned.
;       PARALLEL:           in, optional, type=boolean/string, default=0/''
;                           If set, the component parallel to B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component parallel to the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the component perpendicular to B (magnetic field) will
;                               be returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component perpendicular to the quantity specified
;                               will be returned.
;       VECTOR:             in, optional, type=boolean, default=0
;                           If set, all three components of the electric field will be
;                               returned in an NxMx3-dimensional array. If no other
;                               keywords are set, this is the default.
;       X:                  in, optional, type=boolean, default=0
;                           If set, the X-component of the electric field is returned.
;       Y:                  in, optional, type=boolean, default=0
;                           If set, the Y-component of the electric field is returned.
;       Z:                  in, optional, type=boolean, default=0
;                           If set, the Z-component of the electric field is returned.
;
; :Returns:
;       DATA:               Electric field data.
;-
function MrSim2_Data::E, $
CROSS=cross, $
DOT=dot, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    tf_cross  = keyword_set(cross)
    tf_dot    = keyword_set(dot)
    magnitude = keyword_set(magnitude)
    tf_par    = keyword_set(parallel)
    tf_perp   = keyword_set(perpendicular)
    vec       = keyword_set(vec)
    x         = keyword_set(x)
    y         = keyword_set(y)
    z         = keyword_set(z)
    
    ;If nothing was chosen, return the vector.
    vec = vec || (x + y + z + tf_cross + tf_dot + tf_par + tf_perp + magnitude) eq 0
    if vec then begin
        x = 1
        y = 1
        z = 1
    endif

    ;Parallel & perpendicular to what?
    dotName   = MrIsA(dot,           'STRING') ? dot           : 'B'
    crossName = MrIsA(cross,         'STRING') ? cross         : 'B'
    parName   = MrIsA(parallel,      'STRING') ? parellel      : 'B'
    perpName  = MrIsA(perpendicular, 'STRING') ? perpendicular : 'B'
    
    ;Components
    if x then if n_elements(*self.Ex) eq 0 then self -> ReadData, 'Ex'
    if y then if n_elements(*self.Ey) eq 0 then self -> ReadData, 'Ey'
    if z then if n_elements(*self.Ez) eq 0 then self -> ReadData, 'Ez'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        magnitude: return, self -> Magnitude('E')
        tf_cross:  return, self -> CrossProduct('E', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> DotProduct('E', dotName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Parallel('E', parName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Perpendicular('E', perpName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, [[[*self.Ex], [*self.Ey], [*self.Ez]]]
        x:         return, *self.Ex
        y:         return, *self.Ey
        z:         return, *self.Ez
        else:      ;Not possible
    endcase
end


;+
;   Return the ion density.
;
; :Private:
;
; :Returns:
;       n_i:                Ion density data.
;-
function MrSim2_Data::n_i
    compile_opt strictarr, hidden
    on_error, 2
    
    ;ni
    if n_elements(*self.n_i) eq 0 then self -> ReadData, 'ni'
    return, *self.n_i
end


;+
;   Return the ion density.
;
; :Private:
;
; :Returns:
;       n_e:                Electron density data.
;-
function MrSim2_Data::n_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;ne
    if n_elements(*self.n_e) eq 0 then self -> ReadData, 'ne'
    return, *self.n_e
end


;+
;   Return various data products associated with the electron pressure tensor.
;
; :Private:
;
; :Keywords;
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the magnitude is returned.
;       PARALLEL:           in, optional, type=boolean/string, default=0/''
;                           If set, the component parallel to B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component parallel to the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the component perpendicular to B (magnetic field) will
;                               be returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component perpendicular to the quantity specified
;                               will be returned.
;       TENSOR:             in, optional, type=boolean, default=0
;                           If set, an upper-diagonal matrix of tensor components will be
;                               returned in an NxMx3x3-dimensional array. If no other
;                               keywords are set, this is the default.
;       TXX:                in, optional, type=boolean, default=0
;                           If set, the Txx-component of the electron pressure is returned.
;       TXY:                in, optional, type=boolean, default=0
;                           If set, the Txy-component of the electron pressure is returned.
;       TXZ:                in, optional, type=boolean, default=0
;                           If set, the Txz-component of the electron pressure is returned.
;       TYY:                in, optional, type=boolean, default=0
;                           If set, the Tyy-component of the electron pressure is returned.
;       TYZ:                in, optional, type=boolean, default=0
;                           If set, the Tyz-component of the electron pressure is returned.
;       TZZ:                in, optional, type=boolean, default=0
;                           If set, the Tzz-component of the electron pressure is returned.
;
; :Returns:
;       DATA:               Electron pressure tensor data.
;-
function MrSim2_Data::Pe, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
TENSOR=tensor, $
TXX=Txx, $
TXY=Txy, $
TXZ=Txz, $
TYY=Tyy, $
TYZ=Tyz, $
TZZ=Tzz

    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    magnitude = keyword_set(magnitude)
    tf_par    = keyword_set(parallel)
    tf_perp   = keyword_set(perpendicular)
    tensor    = keyword_set(tensor)
    Txx       = keyword_set(Txx)
    Txy       = keyword_set(Txy)
    Txz       = keyword_set(Txz)
    Tyy       = keyword_set(Tyy)
    Tyz       = keyword_set(Tyz)
    Tzz       = keyword_set(Tzz)

    ;Parallel & perpendicular to what?
    parName  = MrIsA(parallel,      'STRING') ? parellel      : 'B'
    perpName = MrIsA(perpendicular, 'STRING') ? perpendicular : 'B'
    
    ;If nothing was chosen, return the tensor.
    tensor = tensor || (Txx + Txy + Txz + Tyy + Tyz + Tzz + tf_par + tf_perp + magnitude) eq 0
    if vec then begin
        Txx = 1
        Txy = 1
        Txz = 1
        Tyy = 1
        Tyz = 1
        Tzz = 1
    endif
    
    ;Tensor components
    if Txx then if n_elements(*self.Pe_xx) eq 0 then self -> ReadData, 'Pe-xx'
    if Txy then if n_elements(*self.Pe_xy) eq 0 then self -> ReadData, 'Pe-xy'
    if Txz then if n_elements(*self.Pe_xz) eq 0 then self -> ReadData, 'Pe-xz'
    if Tyy then if n_elements(*self.Pe_yy) eq 0 then self -> ReadData, 'Pe-yy'
    if Tyz then if n_elements(*self.Pe_yz) eq 0 then self -> ReadData, 'Pe-yz'
    if Tzz then if n_elements(*self.Pe_zz) eq 0 then self -> ReadData, 'Pe-zz'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Tensor_Par('Pe', X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Tensor_Perp('Pe', X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, self -> Tensor_Mag('Pe')
        Txx:           return, *self.Pe_xx
        Txy:           return, *self.Pe_xy
        Txz:           return, *self.Pe_xz
        Tyy:           return, *self.Pe_yy
        Tyz:           return, *self.Pe_yz
        Tzz:           return, *self.Pe_zz
        tensor: begin
            ;Allocate memory
            dims = size(*self.Pe_xx, /DIMENSIONS)
            data = fltarr(dims[0], dims[1], 3, 3)
            
            ;Create the data product
            data[*,0,0] = *self.Pe_xx
            data[*,0,1] = *self.Pe_xy
            data[*,0,2] = *self.Pe_xz
            data[*,1,1] = *self.Pe_yy
            data[*,1,2] = *self.Pe_yz
            data[*,2,2] = *self.Pe_zz
            return, data
        endcase
    endcase
end


;+
;   Return various data products associated with the electron pressure tensor.
;
; :Private:
;
; :Keywords;
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the magnitude is returned.
;       PARALLEL:           in, optional, type=boolean/string, default=0/''
;                           If set, the component parallel to B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component parallel to the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the component perpendicular to B (magnetic field) will
;                               be returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component perpendicular to the quantity specified
;                               will be returned.
;       TENSOR:             in, optional, type=boolean, default=0
;                           If set, an upper-diagonal matrix of tensor components will be
;                               returned in an NxMx3x3-dimensional array. If no other
;                               keywords are set, this is the default.
;       TXX:                in, optional, type=boolean, default=0
;                           If set, the Txx-component of the ion pressure is returned.
;       TXY:                in, optional, type=boolean, default=0
;                           If set, the Txy-component of the ion pressure is returned.
;       TXZ:                in, optional, type=boolean, default=0
;                           If set, the Txz-component of the ion pressure is returned.
;       TYY:                in, optional, type=boolean, default=0
;                           If set, the Tyy-component of the ion pressure is returned.
;       TYZ:                in, optional, type=boolean, default=0
;                           If set, the Tyz-component of the ion pressure is returned.
;       TZZ:                in, optional, type=boolean, default=0
;                           If set, the Tzz-component of the ion pressure is returned.
;
; :Returns:
;       DATA:               Ion pressure tensor data.
;-
function MrSim2_Data::Pi, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
TENSOR=tensor, $
TXX=Txx, $
TXY=Txy, $
TXZ=Txz, $
TYY=Tyy, $
TYZ=Tyz, $
TZZ=Tzz

    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    magnitude = keyword_set(magnitude)
    tf_par    = keyword_set(parallel)
    tf_perp   = keyword_set(perpendicular)
    tensor    = keyword_set(tensor)
    Txx       = keyword_set(Txx)
    Txy       = keyword_set(Txy)
    Txz       = keyword_set(Txz)
    Tyy       = keyword_set(Tyy)
    Tyz       = keyword_set(Tyz)
    Tzz       = keyword_set(Tzz)

    ;Parallel & perpendicular to what?
    parName  = MrIsA(parallel,      'STRING') ? parellel      : 'B'
    perpName = MrIsA(perpendicular, 'STRING') ? perpendicular : 'B'
    
    ;If nothing was chosen, return the tensor.
    tensor = tensor || (Txx + Txy + Txz + Tyy + Tyz + Tzz + tf_par + tf_perp + magnitude) eq 0
    if vec then begin
        Txx = 1
        Txy = 1
        Txz = 1
        Tyy = 1
        Tyz = 1
        Tzz = 1
    endif
    
    ;Tensor components
    if Txx then if n_elements(*self.Pi_xx) eq 0 then self -> ReadData, 'Pi-xx'
    if Txy then if n_elements(*self.Pi_xy) eq 0 then self -> ReadData, 'Pi-xy'
    if Txz then if n_elements(*self.Pi_xz) eq 0 then self -> ReadData, 'Pi-xz'
    if Tyy then if n_elements(*self.Pi_yy) eq 0 then self -> ReadData, 'Pi-yy'
    if Tyz then if n_elements(*self.Pi_yz) eq 0 then self -> ReadData, 'Pi-yz'
    if Tzz then if n_elements(*self.Pi_zz) eq 0 then self -> ReadData, 'Pi-zz'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Tensor_Par('Pi', X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Tensor_Pirp('Pi', X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, self -> Tensor_Mag('Pi')
        Txx:           return, *self.Pi_xx
        Txy:           return, *self.Pi_xy
        Txz:           return, *self.Pi_xz
        Tyy:           return, *self.Pi_yy
        Tyz:           return, *self.Pi_yz
        Tzz:           return, *self.Pi_zz
        tensor: begin
            ;Allocate memory
            dims = size(*self.Pi_xx, /DIMENSIONS)
            data = fltarr(dims[0], dims[1], 3, 3)
            
            ;Create the data product
            data[*,0,0] = *self.Pi_xx
            data[*,0,1] = *self.Pi_xy
            data[*,0,2] = *self.Pi_xz
            data[*,1,1] = *self.Pi_yy
            data[*,1,2] = *self.Pi_yz
            data[*,2,2] = *self.Pi_zz
            return, data
        endcase
    endcase
end


;+
;   Return various data products associated with the ion bulk velocity.
;
; :Private:
;
; :Keywords:
;       CROSS:              in, optional, type=boolean/string, default=0/''
;                           If set, the cross product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the cross product with the quantity specified will
;                               be returned.
;       DOT:                in, optional, type=boolean/string, default=0/''
;                           If set, the dot product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the dot product with the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the magnitude is returned.
;       PARALLEL:           in, optional, type=boolean/string, default=0/''
;                           If set, the component parallel to B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component parallel to the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the component perpendicular to B (magnetic field) will
;                               be returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component perpendicular to the quantity specified
;                               will be returned.
;       VECTOR:             in, optional, type=boolean, default=0
;                           If set, all three components of the electric field will be
;                               returned in an NxMx3-dimensional array. If no other
;                               keywords are set, this is the default.
;       X:                  in, optional, type=boolean, default=0
;                           If set, the X-component is returned.
;       Y:                  in, optional, type=boolean, default=0
;                           If set, the Y-component is returned.
;       Z:                  in, optional, type=boolean, default=0
;                           If set, the Z-component is returned.
;
; :Returns:
;       DATA:               Ion bulk velocity data.
;-
function MrSim2_Data::Ui, $
CROSS=cross, $
DOT=dot, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    tf_cross  = keyword_set(cross)
    tf_dot    = keyword_set(dot)
    magnitude = keyword_set(magnitude)
    tf_par    = keyword_set(parallel)
    tf_perp   = keyword_set(perpendicular)
    vec       = keyword_set(vec)
    x         = keyword_set(x)
    y         = keyword_set(y)
    z         = keyword_set(z)
    
    ;If nothing was chosen, return the vector.
    vec = vec || (x + y + z + tf_cross + tf_dot + tf_perp + tf_par + magnitude) eq 0
    if vec then begin
        x = 1
        y = 1
        z = 1
    endif

    ;Parallel & perpendicular to what?
    dotName   = MrIsA(dot,           'STRING') ? dot           : 'B'
    crossName = MrIsA(cross,         'STRING') ? cross         : 'B'
    parName   = MrIsA(parallel,      'STRING') ? parellel      : 'B'
    perpName  = MrIsA(perpendicular, 'STRING') ? perpendicular : 'B'
    
    ;[XYZ]-Components
    if x then if n_elements(*self.Uix) eq 0 then self -> ReadData, 'Uix'
    if y then if n_elements(*self.Uiy) eq 0 then self -> ReadData, 'Uiy'
    if z then if n_elements(*self.Uiz) eq 0 then self -> ReadData, 'Uiz'
    
    ;Which to return?
    ;   - The order is crucial.
    case 1 of
        magnitude: return, self -> Magnitude('Ui')
        tf_cross:  return, self -> CrossProduct('Ui', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> DotProduct('Ui', dotName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Parallel( 'Ui', parName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Perpendicular('Ui', perpName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, [[[*self.Uix]], [[*self.Uiy]], [[*self.Uiz]]]
        x:         return, *self.Uix
        y:         return, *self.Uiy
        z:         return, *self.Uiz
    endcase
end


;+
;   Return various data products associated with the electron bulk velocity.
;
; :Private:
;
; :Keywords:
;       CROSS:              in, optional, type=boolean/string, default=0/''
;                           If set, the cross product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the cross product with the quantity specified will
;                               be returned.
;       DOT:                in, optional, type=boolean/string, default=0/''
;                           If set, the dot product with B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the dot product with the quantity specified will
;                               be returned.
;       MAGNITUDE:          in, optional, type=boolean, default=0
;                           If set, the magnitude is returned.
;       PARALLEL:           in, optional, type=boolean/string, default=0/''
;                           If set, the component parallel to B (magnetic field) will be
;                               returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component parallel to the quantity specified will
;                               be returned.
;       PERPENDICULAR:      in, optional, type=boolean, default=0
;                           If set, the component perpendicular to B (magnetic field) will
;                               be returned. If a string is given, it must be the name of
;                               a vector quantity ('B' for the magnetic field), in which
;                               case the component perpendicular to the quantity specified
;                               will be returned.
;       VECTOR:             in, optional, type=boolean, default=0
;                           If set, all three components of the electric field will be
;                               returned in an NxMx3-dimensional array. If no other
;                               keywords are set, this is the default.
;       X:                  in, optional, type=boolean, default=0
;                           If set, the X-component is returned.
;       Y:                  in, optional, type=boolean, default=0
;                           If set, the Y-component is returned.
;       Z:                  in, optional, type=boolean, default=0
;                           If set, the Z-component is returned.
;
; :Returns:
;       DATA:               Electron bulk velocity data.
;-
function MrSim2_Data::Ue, $
CROSS=cross, $
DOT=dot, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Defaults
    magnitude = keyword_set(magnitude)
    tf_par    = keyword_set(parallel)
    tf_perp   = keyword_set(perpendicular)
    vec       = keyword_set(vec)
    x         = keyword_set(x)
    y         = keyword_set(y)
    z         = keyword_set(z)
    
    ;If nothing was chosen, return the vector.
    vec = vec || (x + y + z + tf_cross + tf_dot + tf_perp + tf_par + magnitude) eq 0
    if vec then begin
        x = 1
        y = 1
        z = 1
    endif

    ;Parallel & perpendicular to what?
    dotName   = MrIsA(dot,           'STRING') ? dot           : 'B'
    crossName = MrIsA(cross,         'STRING') ? cross         : 'B'
    parName   = MrIsA(parallel,      'STRING') ? parellel      : 'B'
    perpName  = MrIsA(perpendicular, 'STRING') ? perpendicular : 'B'
    
    if x then if n_elements(*self.Uex) eq 0 then self -> ReadData, 'Uex'
    if y then if n_elements(*self.Uey) eq 0 then self -> ReadData, 'Uey'
    if z then if n_elements(*self.Uez) eq 0 then self -> ReadData, 'Uez'
    
    ;Which to return?
    ;   - The order is crucial.
    case 1 of
        magnitude: return, self -> Magnitude('Ue')
        tf_cross:  return, self -> CrossProduct('B', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> DotProduct('B', dotName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Parallel( 'Ue', parName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Perpendicular('Ue', perpName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, [[[*self.Uex]], [[*self.Uey]], [[*self.Uez]]]
        x:         return, *self.Uex
        y:         return, *self.Uey
        z:         return, *self.Uez
    endcase
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
;   The purpose of this program is to calculate the x-component of the current density::
;       J_{x} = ( -n_{e} * U_{ex} ) + ( n_{i} * U_{ix} )
;
; :Private:
;
; :Returns:
;       Jx:                     The electric field dotted with the current density.
;-
function MrSim2_Data::J, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
MAGNITUDE=magnitude, $
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
	
	;Get the data
	;   - Number density
	;   - Bulk velocity
    n_e = self -> GetData('ne')
    n_e = self -> GetData('ni')
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, X=x, Y=y, Z=z, VECTOR=vec, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, X=x, Y=y, Z=z, VECTOR=vec, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)

    ;Total Current
    if vec then begin
        J = Ui
        J[*,*,0] = (-n_e * Ue[*,*,0]) + (n_i * Ui[*,*,0])
        J[*,*,1] = (-n_e * Ue[*,*,1]) + (n_i * Ui[*,*,1])
        J[*,*,2] = (-n_e * Ue[*,*,2]) + (n_i * Ui[*,*,2])
    endif else begin
        J = (-n_e * Ue) + (n_i * Ui)
    endelse
    
    return, J
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
function MrSim2_Data::Je, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
MAGNITUDE=magnitude, $
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	on_error, 2
	
	;Number density
	;   - ne is always stored as an object property.
	;   - Only the X, Y, and Z components of Ue are object properties. Since we
	;       do not know which was chosen without lots of checks, retrieve the data.
    n_e = self -> GetData('ne')
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, X=x, Y=y, Z=z, VECTOR=vec, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    
    ;Compute the current
    if vec then begin
        Je = Ue
        Je[*,*,0] = -n_e * Ue[*,*,0]
        Je[*,*,1] = -n_e * Ue[*,*,1]
        Je[*,*,2] = -n_e * Ue[*,*,2]
    endif else begin
        Je = -1.0 * temporary(n_e) * temporary(Ue)
    endelse
    
    return, Je
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
function MrSim2_Data::Ji, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	on_error, 2
	
	;Number density
    n_i = self -> GetData('ni')
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    
    ;Compute the current
    Ji = temporary(n_i) * temporary(Ui)
    
    return, Ji
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
function MrSim2_Data::V, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	
	;Get the Electric and Magnetic Field data
	self -> GetInfo, MI_ME=mi_me
	mi = mi_me
	me = 1.0
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    
    ;Calculate the magnitude of ExB
    V = (mi*Ui + me*Ue) / (mi + me)
    
    return, Vx
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
function MrSim2_Data::Fe, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	
	;Get the data
	E    = self -> GetData('E', X=x, Y=y, Z=z)
	UexB = self -> CrossProduct('Ue', 'B', X=x, Y=y, Z=z)
	
	;Compute the force
	Fe = -E + UexB
    
    ;Return the correct quantity
    case 1 of
        magnitude: return, sqrt(total(Fe^2, 3))
        x:         return, Fe[*,*,0]
        y:         return, Fe[*,*,1]
        z:         return, Fe[*,*,2]
        else:      return, Fe
    endcase
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
;       Fi:                     The lorentz force.
;-
function MrSim2_Data::Fi, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	
	;Get the data
	E    = self -> GetData('E', X=x, Y=y, Z=z)
	UixB = self -> CrossProduct('Ui', 'B', X=x, Y=y, Z=z)
	
	;Compute the force
	Fi = E + UixB
    
    ;Return the correct quantity
    case 1 of
        magnitude: return, sqrt(total(Fi^2, 3))
        x:         return, Fi[*,*,0]
        y:         return, Fi[*,*,1]
        z:         return, Fi[*,*,2]
        else:      return, Fi
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
function MrSim2_Data::v_ExB, $
MAGNITUDE=magnitude, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	on_error, 2
	
	;Get the Electric and Magnetic Field data
	ExB  = self -> CrossProduct('E', 'B')
	Bmag = self -> GetData('B', /MAGNITUDE)
    
    ;Compute the drift velocity
    v_ExB = ExB / Bmag^2
    
    ;Return the correct quantity
    case 1 of
        magnitude: return, sqrt(total(v_ExB^2, 3))
        x:         return, v_ExB[*,*,0]
        y:         return, v_ExB[*,*,1]
        z:         return, v_ExB[*,*,2]
        else:      return, v_ExB
    endcase
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
function MrSim2_Data::Beta_e
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
function MrSim2_Data::Beta_i
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
function MrSim2_Data::Beta_p
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
function MrSim2_Data::divPe_mag
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
function MrSim2_Data::divPe_x
	compile_opt strictarr, hidden
	on_error, 2

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
	
	;Get the electron density and the required components of the electron pressure tensor
    Pe_xx = self -> getData('Pe-xx')
    Pe_xy = self -> getData('Pe-xy')
    Pe_xz = self -> getData('Pe-xz')
    n_e = self -> getData('ne')
    
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
function MrSim2_Data::divPe_z
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
function MrSim2_Data::divPi_x
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
function MrSim2_Data::divPi_z
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
function MrSim2_Data::gradP
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
;   The purpose of this program is to calculate the electron inertial length::
;       \lambda_{e} = c / \omega_{p,e}
;
; :Private:
;
; :Returns:
;       LAMBDA_e:               Electron inertial length (skin depth).
;-
function MrSim2_Data::lambda_e
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
function MrSim2_Data::lambda_i
    compile_opt strictarr, hidden

    ;Get the electric field and current densitz data
    c = constants('c')
    w_pe = self -> getData('w_pi')

    ;Calculate the dot product
    lambda_i = c / w_pi
    
    return, lambda_i
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
function MrSim2_Data::v_Ae
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
function MrSim2_Data::v_A
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
function MrSim2_Data::w_pe
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
function MrSim2_Data::w_pi
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
function MrSim2_Data::w_ci
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
function MrSim2_Data::w_ci
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
;   This definition statement for the MrSim2_Data class.
;
; :Params:
;       CLASS           out, optional, type=named structure
;                       The class definition structure.
;
; :FIELDS:
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
pro MrSim2_Data__DEFINE, class
    compile_opt strictarr
    
    class = { MrSim2_Data, $

              ;Data
              Ay:        ptr_new(), $
              Bx:        ptr_new(), $
              By:        ptr_new(), $
              Bz:        ptr_new(), $
              electrons: ptr_new(), $
              Ex:        ptr_new(), $
              Ey:        ptr_new(), $
              Ez:        ptr_new(), $
              n_e:       ptr_new(), $
              n_i:       ptr_new(), $
              Pe_xx:     ptr_new(), $
              Pe_xy:     ptr_new(), $
              Pe_xz:     ptr_new(), $
              Pe_yy:     ptr_new(), $
              Pe_yz:     ptr_new(), $
              Pe_zz:     ptr_new(), $
              Pi_xx:     ptr_new(), $
              Pi_xy:     ptr_new(), $
              Pi_xz:     ptr_new(), $
              Pi_yy:     ptr_new(), $
              Pi_yz:     ptr_new(), $
              Pi_zz:     ptr_new(), $
              Uex:       ptr_new(), $
              Uey:       ptr_new(), $
              Uez:       ptr_new(), $
              Uix:       ptr_new(), $
              Uiy:       ptr_new(), $
              Uiz:       ptr_new() $
            }
end