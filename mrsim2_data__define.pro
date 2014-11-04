; docformat = 'rst'
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
; PURPOSE
;+
;   The purpose of this program is to provide a base class for 2D and 3D simulation
;   data provided by Bill Daughton from Los Alamos National Laboratory.
;
;   MrSim2_Data is meant to be subclassed by a class that over-rides the ReadGDA method.
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
    self.answer    = ptr_new(/ALLOCATE_HEAP)
    self.electrons = ptr_new(/ALLOCATE_HEAP)
    self.Ay        = ptr_new(/ALLOCATE_HEAP)
    self.Bx        = ptr_new(/ALLOCATE_HEAP)
    self.By        = ptr_new(/ALLOCATE_HEAP)
    self.Bz        = ptr_new(/ALLOCATE_HEAP)
    self.electrons = ptr_new(/ALLOCATE_HEAP)
    self.Ex        = ptr_new(/ALLOCATE_HEAP)
    self.Ey        = ptr_new(/ALLOCATE_HEAP)
    self.Ez        = ptr_new(/ALLOCATE_HEAP)
    self.n_e       = ptr_new(/ALLOCATE_HEAP)
    self.n_i       = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xx     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xy     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_xz     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yy     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_yz     = ptr_new(/ALLOCATE_HEAP)
    self.Pe_zz     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xx     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xy     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_xz     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yy     = ptr_new(/ALLOCATE_HEAP)
    self.Pi_yz     = ptr_new(/ALLOCATE_HEAP)
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
    ptr_free, self.answer
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
;   The purpose of this method is to release memory taken up by simulation data.
;
; :Private:
;
; :Params:
;       DATA_PRODUCT:           in, optional, type=string
;                               The name of the data product whose data is to be released.
;                                   If not given, the data for all data products is freed.
;-
pro MrSim2_Data::Clear_Data, data_product
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
        data_product = ['Ans', 'Ay', 'Bx', 'By', 'Bz', 'E-', 'Ex', 'Ey', 'Ez', 'ne', 'ni', $
                        'Pe_xx', 'Pe_xy', 'Pe_xz', 'Pe_yy', 'Pe_yz', 'Pe_zz', $
                        'Pi_xx', 'Pi_xy', 'Pi_xz', 'Pi_yy', 'Pi_yz', 'Pi_zz', $
                        'Uex', 'Uey', 'Uez', 'Uix', 'Uiy', 'Uiz']
    endif
    
    ;Step through all of the data products given and free the data
    foreach name, data_product do begin
        case name of
            'Ans':   *self.answer    = !Null
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
            'Pe_xx': *self.Pe_xx     = !Null
            'Pe_xy': *self.Pe_xy     = !Null
            'Pe_xz': *self.Pe_xz     = !Null
            'Pe_yy': *self.Pe_yy     = !Null
            'Pe_yz': *self.Pe_yz     = !Null
            'Pe_zz': *self.Pe_zz     = !Null
            'Pi_xx': *self.Pi_xx     = !Null
            'Pi_xy': *self.Pi_xy     = !Null
            'Pi_xz': *self.Pi_xz     = !Null
            'Pi_yy': *self.Pi_yy     = !Null
            'Pi_yz': *self.Pi_yz     = !Null
            'Pi_zz': *self.Pi_zz     = !Null
            'Uex':   *self.Uex       = !Null
            'Uey':   *self.Uey       = !Null
            'Uez':   *self.Uez       = !Null
            'Uix':   *self.Uix       = !Null
            'Uiy':   *self.Uiy       = !Null
            'Uiz':   *self.Uiz       = !Null
            else: message, 'Data product "' + name + '" does not exist.', /INFORMATIONAL
        endcase
    endforeach
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
function MrSim2_Data::CrossProduct, v1, v2, $
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
;   A helper method for the ::GetData method. Used to retrieve single data products when
;   no operation is specified.
;
; :Private:
;
; :Params:
;       NAME:               in, required, type=string, default=
;                           The name of the data product to be read. For a list of
;                               available data product, call mr_readSIM without any
;                               arguments.
;
; :Returns:
;       DATA:               The requested data. If the data product does not exist,
;                               then !Null will be returned.
;-
function MrSim2_Data::DataGet, name
    compile_opt strictarr
    on_error, 2

    ;If no parameters were given, print a list of available data products.
    if n_params() eq 0 then begin
        message, 'Use: data = mySim -> GetData(name)'
        MrSim -> ListProducts
        return, !Null
    endif

    ;Check to see if the data has already been read first.
    case strupcase(name) of
        'ANS':    data = *self.answer
        'AY':     data = self -> A(/Y)
        'BX':     data = self -> B(/X)
        'BY':     data = self -> B(/Y)
        'BZ':     data = self -> B(/Z)
        'E-':     data = *self.electrons
        'EX':     data = self -> E(/X)
        'EY':     data = self -> E(/Y)
        'EZ':     data = self -> E(/Z)
        'NE':     data = self -> n_i()
        'NI':     data = self -> n_e()
        'PE_XX':  data = self -> Pe(/TXX)
        'PE_XY':  data = self -> Pe(/TXY)
        'PE_XZ':  data = self -> Pe(/TXZ)
        'PE_YX':  data = self -> Pe(/TXY)
        'PE_YY':  data = self -> Pe(/TYY)
        'PE_YZ':  data = self -> Pe(/TYZ)
        'PE_ZX':  data = self -> Pe(/TXZ)
        'PE_ZY':  data = self -> Pe(/TYZ)
        'PE_ZZ':  data = self -> Pe(/TZZ)
        'PI_XX':  data = self -> Pi(/TXX)
        'PI_XY':  data = self -> Pi(/TXY)
        'PI_XZ':  data = self -> Pi(/TXZ)
        'PI_YX':  data = self -> Pi(/TXY)
        'PI_YY':  data = self -> Pi(/TYY)
        'PI_YZ':  data = self -> Pi(/TYZ)
        'PI_ZX':  data = self -> Pi(/TXZ)
        'PI_ZY':  data = self -> Pi(/TYZ)
        'PI_ZZ':  data = self -> Pi(/TZZ)
        'UEX':    data = self -> Ue(/X)
        'UEY':    data = self -> Ue(/Y)
        'UEZ':    data = self -> Ue(/Z)
        'UIX':    data = self -> Ui(/X)
        'UIY':    data = self -> Ui(/Y)
        'UIZ':    data = self -> Ui(/Z)
        
        ;Custom Data Products
        'B':      data = self ->  B(/VECTOR)
        'E':      data = self ->  E(/VECTOR)
        'PE':     data = self -> Pe(/TENSOR)
        'PI':     data = self -> Pi(/TENSOR)
        'Ue':     data = self -> Ue(/VECTOR)
        'Ui':     data = self -> Ui(/VECTOR)
        'A0_E':   data = self -> A0_e()
        'AN_E':   data = self -> An_e()
        'DNG_E':  data = self -> Dng_e()
        else: message, 'Data product not available: "' + name + '".'
    endcase
    
    return, data
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
function MrSim2_Data::DotProduct, v1, v2, $
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
function MrSim2_Data::d_dx, data
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
function MrSim2_Data::d_dy, data
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
function MrSim2_Data::d_dz, data
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
function MrSim2_Data::GetData, expression, $
SHOW=show, $
TEST=test, $
DIAGONAL=diagonal, $
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
        if n_elements(expr) gt 0 then message, 'Error evaluating: "' + expr + '".', /INFORMATIONAL
        void = cgErrorMSG()
        return, !Null
    endif

    ;If no parameters were given, print a list of available data products.
    if n_params() eq 0 then begin
        MrSim2 -> listProducts
        return, !Null
    endif

    ;Defaults
    show = keyword_set(show)
    test = keyword_set(test)

    ;Determine order of operations
    expr = self -> Op_Parser(expression, OPERATIONS=ops, COUNT=nOps)
    if expr eq '' $
        then message, 'Unable to parse expression: "' + expression + '".' $
        else if show then print, expr

    ;Test expression and return
    if test then begin
        print, expr
        return, !Null
    endif

    ;If no operations are present, then EXPR is a data product
    if nOps eq 0 then return, self -> DataGet(expr)
    
    ;Step through each operation
    for i = 0, nOps - 1 do begin
        lhs = self -> DataGet(ops[0,i])
        rhs = self -> DataGet(ops[2,i])
        
        ;Show the computation
        if show then begin
            ;The "what"
            print, FORMAT='(%"ans = %s %s %s")', ops[*,i]
            
            ;The "how"
            lhs_dims = strjoin(strtrim(size(lhs, /DIMENSIONS), 2), 'x')
            rhs_dims = strjoin(strtrim(size(rhs, /DIMENSIONS), 2), 'x')
            print, FORMAT='(%"ans = %s %s %s")', lhs_dims, ops[1,i], rhs_dims
        endif
        
        ;Execute the operation
        case ops of
            '+': *self.answer = lhs + rhs
            '-': *self.answer = lhs - rhs
            '*': *self.answer = lhs * rhs
            '/': *self.answer = lhs / rhs
            '^': *self.answer = lhs ^ rhs
            '.': *self.answer = self -> DotProduct(  lhs, rhs, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            'X': *self.answer = self -> CrossProduct(lhs, rhs, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
            else: message, 'Operation invalid: "' + ops[1,i] + '".'
        endcase
    endfor
    
    ;Return the result
    return, *self.answer
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
function MrSim2_Data::HasData, data_product
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
        'ANSWER': if n_elements(*self.answer)   gt 0 then tf_has = 1B
        'AY':     if n_elements(*self.Ay)        gt 0 then tf_has = 1B
        'BX':     if n_elements(*self.Bx)        gt 0 then tf_has = 1B
        'BY':     if n_elements(*self.By)        gt 0 then tf_has = 1B
        'BZ':     if n_elements(*self.Bz)        gt 0 then tf_has = 1B
        'E-':     if n_elements(*self.electrons) gt 0 then tf_has = 1B
        'EX':     if n_elements(*self.Ex)        gt 0 then tf_has = 1B
        'EY':     if n_elements(*self.Ey)        gt 0 then tf_has = 1B
        'EZ':     if n_elements(*self.Ez)        gt 0 then tf_has = 1B
        'NE':     if n_elements(*self.n_e)       gt 0 then tf_has = 1B
        'NI':     if n_elements(*self.n_i)       gt 0 then tf_has = 1B
        'PE_XX':  if n_elements(*self.Pe_xx)     gt 0 then tf_has = 1B
        'PE_XY':  if n_elements(*self.Pe_xy)     gt 0 then tf_has = 1B
        'PE_XZ':  if n_elements(*self.Pe_xz)     gt 0 then tf_has = 1B
        'PE_YX':  if n_elements(*self.Pe_yx)     gt 0 then tf_has = 1B
        'PE_YY':  if n_elements(*self.Pe_yy)     gt 0 then tf_has = 1B
        'PE_YZ':  if n_elements(*self.Pe_yz)     gt 0 then tf_has = 1B
        'PE_ZX':  if n_elements(*self.Pe_zx)     gt 0 then tf_has = 1B
        'PE_ZY':  if n_elements(*self.Pe_zy)     gt 0 then tf_has = 1B
        'PE_ZZ':  if n_elements(*self.Pe_zz)     gt 0 then tf_has = 1B
        'PI_XX':  if n_elements(*self.Pi_xx)     gt 0 then tf_has = 1B
        'PI_XY':  if n_elements(*self.Pi_xy)     gt 0 then tf_has = 1B
        'PI_XZ':  if n_elements(*self.Pi_xz)     gt 0 then tf_has = 1B
        'PI_YX':  if n_elements(*self.Pi_yx)     gt 0 then tf_has = 1B
        'PI_YY':  if n_elements(*self.Pi_yy)     gt 0 then tf_has = 1B
        'PI_YZ':  if n_elements(*self.Pi_yz)     gt 0 then tf_has = 1B
        'PI_ZX':  if n_elements(*self.Pi_zx)     gt 0 then tf_has = 1B
        'PI_ZY':  if n_elements(*self.Pi_zy)     gt 0 then tf_has = 1B
        'PI_ZZ':  if n_elements(*self.Pi_zz)     gt 0 then tf_has = 1B
        'UEX':    if n_elements(*self.Uex)       gt 0 then tf_has = 1B
        'UEY':    if n_elements(*self.Uey)       gt 0 then tf_has = 1B
        'UEZ':    if n_elements(*self.Uez)       gt 0 then tf_has = 1B
        'UIX':    if n_elements(*self.Uix)       gt 0 then tf_has = 1B
        'UIY':    if n_elements(*self.Uiy)       gt 0 then tf_has = 1B
        'UIZ':    if n_elements(*self.Uiz)       gt 0 then tf_has = 1B
        else: tf_has = -1
    endcase
    
    return, tf_has
end


;+
;   Print a list of data products to the command window.
;-
pro MrSim2_Data::ListProducts
    compile_opt strictarr
    on_error, 2
    
    ;Operations
    ops = [['+', 'Addition',       'Bx + By'], $
           ['-', 'Subtraction',    'ne - ni'], $
           ['*', 'Multiplication', 'ni * Ui'], $
           ['/', 'Division',       'ni / ne'], $
           ['^', 'Exponentiation', 'Uex^2'], $
           ['.', 'Dot Product',    'E . B'], $
           ['X', 'Cross Product',  'Ue X B']]
    
    ;Data products
    products = [['Ay',    'Vector potential'], $
                ['B',     'Magnetic field'], $
                ['Bx',    'Magnetic field, x-component'], $
                ['By',    'Magnetic field, y-component'], $
                ['Bz',    'Magnetic field, z-component'], $
                ['E',     'Electric field'], $
                ['Ex',    'Electric field, x-component'], $
                ['Ey',    'Electric field, y-component'], $
                ['Ez',    'Electric field, z-component'], $
                ['ne',    'Electron density'], $
                ['ni',    'Ion density'], $
                ['Pe',    'Electron pressure tensor'], $
                ['Pe_xx', 'Electron pressure tensor, xx-component'], $
                ['Pe_xy', 'Electron pressure tensor, xy-component'], $
                ['Pe_xz', 'Electron pressure tensor, xz-component'], $
                ['Pe_yy', 'Electron pressure tensor, yy-component'], $
                ['Pe_yz', 'Electron pressure tensor, yz-component'], $
                ['Pe_zz', 'Electron pressure tensor, zz-component'], $
                ['Pi',    'Ion pressure tensor'], $
                ['Pi_xx', 'Ion pressure tensor, xx-component'], $
                ['Pi_xy', 'Ion pressure tensor, xy-component'], $
                ['Pi_xz', 'Ion pressure tensor, xz-component'], $
                ['Pi_yy', 'Ion pressure tensor, yy-component'], $
                ['Pi_yz', 'Ion pressure tensor, yz-component'], $
                ['Pi_zz', 'Ion pressure tensor, zz-component'], $
                ['Ue',    'Electron bulk velocity'], $
                ['Uex',   'Electron bulk velocity, x-component'], $
                ['Uey',   'Electron bulk velocity, y-component'], $
                ['Uez',   'Electron bulk velocity, z-component'], $
                ['Ui',    'Ion bulk velocity'], $
                ['Uix',   'Ion bulk velocity, x-component'], $
                ['Uiy',   'Ion bulk velocity, y-component'], $
                ['Uiz',   'Ion bulk velocity, z-component'], $
                ['An_e',  'Electron anisotropy'], $
                ['A0_e',  'Electron agyrotropy'], $
                ['Dng_e', 'Electron non-gyrotropy']]
    
    ;Print Data Products
    print, FORMAT='(%"  %s      %s")', 'NAME', 'DESCRIPTION'
    print, FORMAT='(%"  %-6s     %s")', products
    print, ''
    
    ;Print operations
    print, FORMAT='(%"  %s      %s       %s")', 'OPERATION', 'DESCRIPTION', 'EXAMPLE'
    print, FORMAT='(%"      %s          %-14s    %s")', ops
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
function MrSim2_Data::Magnitude, vx, vy, vz
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
;   Parse the next value from the string of operations.
;
; :Private:
;
; :Params:
;       OPERATOR:       in, required, type=string
;                       The operator for which the associativity is to be determined.
;
; :Returns:
;       ASSOCIATIVITY:  Associativity of the given operator. Returns 1 for left-associative
;                           and 2 for right-associative.
;-
function MrSim2_Data::Op_GetAssoc, operator
    on_error, 2

    ;Define order of operations.
    case operator of
        'EQ': associativity = 'LEFT'
        '+':  associativity = 'LEFT'
        '-':  associativity = 'LEFT'
        '*':  associativity = 'LEFT'
        '/':  associativity = 'LEFT'
        '^':  associativity = 'LEFT'
        '.':  associativity = 'RIGHT'
        'X':  associativity = 'RIGHT'
        '':   associativity = ''
        else: message, 'Operator not recognized: "' + operator + '".'
    endcase

    return, associativity
end


;+
;   Parse the next value from the string of operations.
;
; :Private:
;
; :Params:
;       OPERATOR:       in, required, type=string
;                       The operator for which the precedence is to be determined.
;
; :Returns:
;       PRECEDENCE:     Precedence of the given operator.
;-
function MrSim2_Data::Op_GetPrecedence, operator
    on_error, 2

    ;Define order of operations.
    case operator of
        'EQ': precedence =  0   ;Equal to
        '+':  precedence =  1   ;Addition
        '-':  precedence =  1   ;Subtraction
        '*':  precedence =  2   ;Multiplication
        '/':  precedence =  2   ;Division
        '^':  precedence =  3   ;Exponentiation
        '.':  precedence =  4   ;Dot Product
        'X':  precedence =  5   ;Cross Product
        '':   precedence = -1
        else: message, 'Operator not recognized: "' + operator + '".'
    endcase

    return, precedence
end


;+
;   Parse the next value from the string of operations.
;
; :Private:
;
; :Params:
;       OPSTR:          in, required, type=string
;                       A string from which the value is extracted. The operator is
;                           expected to begin at the first character in the string.
;       REMAINDER:      out, optional, type=string
;                       The remainder of `OPSTR` after the operator has been extracted.
;
; :Returns:
;       ATOM:           The extracted value.
;-
function MrSim2_Data::Op_NextAtom, opStr, remainder, $
COUNT=count, $
OPERATIONS=operations
    on_error, 2

    _opStr = strtrim(opStr, 2)
;-------------------------------------------------------
; Subexpression ////////////////////////////////////////
;-------------------------------------------------------
    ;Is the first character a parentheses?
    char = strmid(_opStr, 0, 1)

    if char eq '(' then begin
        pos   = 1
        nOpen = 1
        nChars = strlen(_opStr)
    
        ;Find the matching close parens
        ;   - Continue until
        ;       o The parenthesis is closed
        ;       o We reach the end of the string
        while (char ne ')' && nOpen gt 0) || (pos lt nChars) do begin
            ;Get the next character
            char = strmid(_opStr, pos, pos+1)
            
            ;A parenthesis?
            case char of
                '(': nOpen += 1
                ')': nOpen -= 1
                else: ;Ignore
            endcase
            
            ;Next character
            pos += 1
        endwhile
        
        ;Parentheses not balanced?
        if nOpen ne 0 then message, 'Parentheses are not balanced.'

        ;Evaluate the subexpression within parentheses
        ;   - Extract the expression
        ;   - Evaluate the expression -- it is the next atom.
        expression = strmid(_opStr, 1, pos-2)
        atom       = self -> Op_Parser(expression, COUNT=count, OPERATIONS=operations)

;-------------------------------------------------------
; Any Non-Operator Sequence ////////////////////////////
;-------------------------------------------------------
    endif else begin
        ;Find the next operator
        opRegEx = '(EQ|\+|-|\*|/|\^|\.|X)'
        pos     = stregex(_opStr, opRegEx)
    
        ;Extract the atom
        ;   - Anything that is not an operator
        ;   - If no operator was found, take the whole string
        ;   - Otherwise, take up to the next operator.
        if pos eq -1 $
            then atom = _opStr $
            else atom = strtrim(strmid(_opStr, 0, pos), 2)
    endelse
    
    ;Extract the remainder
    remainder = strtrim(strmid(_opStr, pos), 2)
    return, atom
end


;+
;   Parse the next operation from the string.
;
; :Private:
;
; :Params:
;       OPSTR:          in, required, type=string
;                       A string from which the operator is extracted. The operator is
;                           expected to begin at the first character in the string.
;       REMAINDER:      out, optional, type=string
;                       The remainder of `OPSTR` after the operator has been extracted.
;
; :Returns:
;       NEXTOP:         The extracted operator.
;-
function MrSim2_Data::Op_NextOp, opStr, remainder, $
PRECEDENCE=precedence, $
ASSOCIATIVITY=associativity
    on_error, 2
    
    ;Extract the next operator
    opRegEx = '(EQ|\+|-|\*|/|\^|\.|X)'
    pos     = stregex(strtrim(opStr, 2), opRegEx, LEN=len)

    ;Extract the operator
    if pos eq -1 $
        then nextOp = '' $
        else nextOp = strmid(opStr, pos, len)

    ;Extract the remainder
    remainder = strmid(opStr, len)
    
    ;Get other features?
    if arg_present(precedence)    then precedence    = self -> Op_GetPrecedence(nextOp)
    if arg_present(associativity) then associativity = self -> Op_GetAssoc(nextOp)

    return, nextOp
end


;+
;   Evaluate the expression.
;
; :Params:
;       EXPRESSION:     in, required, type=string
;                       The expression to be evaluated.
;       MIN_PREC:       in, private, required, type=string
;                       Minimum precedence at which operations should be evaluated.
;       REMAINDER:      in, out, private, required, type=string
;                       When recursing, the part of `EXPRESSION` yet to be evaluated.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of operations evaluated.
;       OPERATIONS:     out, optional, type=strarr(3\,`COUNT`)
;                       The [lhs, op, rhs] of each operation, in the order they are to
;                           be evaluated. "ans" refers to the answer computed in the
;                           previous operation.
;
; :Returns:
;       LHS:            A string joining each operation, parenthesized to demonstrate
;                           the order in which operations should be executed.
;-
function MrSim2_Data::Op_Parser, expression, min_prec, remainder, $
COUNT=count, $
OPERATIONS=operations
    compile_opt strictarr
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        void = cgErrorMSG(/QUIET)
        return, ''
    endif

;-------------------------------------------------------
; Initial Conditions ///////////////////////////////////
;-------------------------------------------------------
    ;Initial conditions of recursion.
    ;   - Extract the left-hand side of the first expression.
    ;   - Precedence starts at 0.
    if n_elements(min_prec) eq 0 then begin
        lhs      = self -> Op_NextAtom(expression, remainder)
        min_prec = 0
        count    = 0
        if n_elements(operations) ne 0 then void = temporary(operations)
    endif else begin
        lhs = expression
    endelse
    
    ;Initial conditions of while loop
    ;   - First operator and the value on which it acts.
    this_op = self -> Op_NextOp(remainder, tempRemain, $
                                PRECEDENCE    = this_prec, $
                                ASSOCIATIVITY = this_assoc)
    rhs     = self -> Op_NextAtom(tempRemain, rhs_remain, $
                                  COUNT=count, OPERATIONS=operations)

;-------------------------------------------------------
; Evaluate Expressions /////////////////////////////////
;-------------------------------------------------------
    ;Example: 2 + 1 + 2
    ;   - Since the minimum precedence of all operators is 0, we will
    ;       loop through all operators until we find one that is < 0
    ;       (i.e. the end-of-expression empty string).
    while this_prec ge min_prec do begin
        ;Look ahead to the next operator.
        next_op = self -> Op_NextOp(rhs_remain, $
                                    PRECEDENCE    = next_prec, $
                                    ASSOCIATIVITY = next_assoc)

    ;-------------------------------------------------------
    ; Next Operation Takes Precedence to Current? //////////
    ;-------------------------------------------------------
        ;Example: 2 + 1 * 2
        ;   - Evaluate as 2 + (1 * 2)
        while next_prec gt this_prec do begin
            ;Example: 2 * 3 * 4
            ;   - 'LEFT':  => ((2 * 3) * 4) = 24
            ;   - 'RIGHT': => (2 * (3 * 2)) = 12
            ;   - Bump up the precedence for left-associative operators so that
            ;       a look-ahead addition does not supercede an earlier addition.
            if this_assoc eq 'LEFT' $
                then next_min_prec = this_prec + 1 $
                else next_min_prec = this_prec
        
            ;Recurse and evaluate the next operation
            rhs = self -> Op_Parser(rhs, next_min_prec, rhs_remain, $
                                    COUNT=count, OPERATIONS=operations)

            ;Next iteration
            ;   - "Look Ahead" to the next operator.
            la_op   = self -> Op_NextOp(rhs_remain, $
                                        PRECEDENCE    = la_prec, $
                                        ASSOCIATIVITY = la_assoc)

            ;Example: 2 + 3^2 * 6
            ;   - Both ^ and * must be evaluated before returning to +
            this_prec = next_prec
            next_prec = la_prec
            if la_prec gt next_prec then lhs = rhs
        endwhile

    ;-------------------------------------------------------
    ; Record Results ///////////////////////////////////////
    ;-------------------------------------------------------
        ;Count and store operations in order.
        if n_elements(operations) eq 0 then begin
            operations = [lhs, this_op, rhs]
            count      = 1
        endif else begin
            ;Substitute "ans" for the part of the expression that is already evaluated.
            if stregex(rhs, '[()]', /BOOLEAN) $
                then operations = [[operations], [ lhs,  this_op, 'ans']] $
                else operations = [[operations], ['ans', this_op,  rhs ]]
            count++
        endelse

        ;Form the operation
        lhs = '(' + lhs + ' ' + this_op + ' ' + rhs + ')'

    ;-------------------------------------------------------
    ; Next Iteration ///////////////////////////////////////
    ;-------------------------------------------------------
        ;Update the remainder to chop of parts that have been evaluated.
        remainder = rhs_remain

        ;Next iteration
        ;   - Get the next operator and value
        this_op   = self -> Op_NextOp(rhs_remain, tempRemain, $
                                      PRECEDENCE    = this_prec, $
                                      ASSOCIATIVITY = this_assoc)
        rhs       = self -> Op_NextAtom(tempRemain, rhs_remain, $
                                        COUNT=count, OPERATIONS=operations)
    endwhile
    
    return, lhs
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
function MrSim2_Data::Parallel, v1, v2, $
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
        x:    v_par =   (*_v1)[*,*,0] * (*_v2)[*,*,0] / temporary(v2_mag)
        y:    v_par =   (*_v1)[*,*,1] * (*_v2)[*,*,1] / temporary(v2_mag)
        z:    v_par =   (*_v1)[*,*,2] * (*_v2)[*,*,2] / temporary(v2_mag)
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
function MrSim2_Data::Perpendicular, v1, v2, $
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
pro MrSim2_Data::ReadGDA, name
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
;       DATA:           in, optional, type=any
;                       The data to be stored. If not provided, data will be read from
;                           the appropriate .gda file.
;-
pro MrSim2_Data::SetData, name, data
    compile_opt strictarr

    ;Catch any errors
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif

    ;Read data from file
    if n_elements(data) eq 0 then data = self -> ReadGDA(name)
    
    ;Store the data in the proper location
    case strupcase(name) of
        'ANSWER': *self.answer    = data
        'AY':     *self.Ay        = data
        'BX':     *self.Bx        = data
        'BY':     *self.By        = data
        'BZ':     *self.Bz        = data
        'E-':     *self.electrons = data
        'EX':     *self.Ex        = data
        'EY':     *self.Ey        = data
        'EZ':     *self.Ez        = data
        'NE':     *self.n_e       = data
        'NI':     *self.n_i       = data
        'PE_XX':  *self.Pe_xx     = data
        'PE_XY':  *self.Pe_xy     = data
        'PE_XZ':  *self.Pe_xz     = data
        'PE_YX':  *self.Pe_yx     = data
        'PE_YY':  *self.Pe_yy     = data
        'PE_YZ':  *self.Pe_yz     = data
        'PE_ZX':  *self.Pe_zx     = data
        'PE_ZY':  *self.Pe_zy     = data
        'PE_ZZ':  *self.Pe_zz     = data
        'PI_XX':  *self.Pi_xx     = data
        'PI_XY':  *self.Pi_xy     = data
        'PI_XZ':  *self.Pi_xz     = data
        'PI_YX':  *self.Pi_yx     = data
        'PI_YY':  *self.Pi_yy     = data
        'PI_YZ':  *self.Pi_yz     = data
        'PI_ZX':  *self.Pi_zx     = data
        'PI_ZY':  *self.Pi_zy     = data
        'PI_ZZ':  *self.Pi_zz     = data
        'UEX':    *self.Uex       = data
        'UEY':    *self.Uey       = data
        'UEZ':    *self.Uez       = data
        'UIX':    *self.Uix       = data
        'UIY':    *self.Uiy       = data
        'UIZ':    *self.Uiz       = data
        else: message, 'Data cannot be set: "' + name + '".'
    endcase
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
function MrSim2_Data::Tensor_Grad, tensor
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_T) then ptr_free, _T
        if ptr_valid(_V) then ptr_free, _V
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Was a name or data given?
    if MrIsA(tensor, 'STRING') $
        then _T = ptr_new(self -> GetData(T)) $
        else _T = ptr_new(T)
    
    return, div_T
end


;+
;   Compute the divergence of a tensor::
;
;       \nabla_{i} T_{ij} = \frac {\partial Tij} {\partial xi}
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
function MrSim2_Data::Tensor_Div, tensor, $
VECTOR=vec, $
X=x, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(T) then ptr_free, T
        void = cgErrorMSG()
        return, !Null
    endif
    
    vec = keyword_set(vec)
    x   = keyword_set(x)
    z   = keyword_set(z)
    
    vec = vec || x + z gt 0
    if vec then begin
        x = 1
        z = 1
    endif
    
    ;Was a name or data given?
    if MrIsA(tensor, 'STRING') $
        then T = ptr_new(self -> GetData(tensor)) $
        else T = ptr_new(T)

;-------------------------------------------------------
;Get Data //////////////////////////////////////////////
;-------------------------------------------------------
    
    dims = size(*T, /DIMENSIONS)
    
    ;Get the x-size of a grid cell in electron skin depth
    self -> GetInfo, DX_DE=dx_de, DY_DE=dy_de, DZ_DE=dz_de

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------

    ;Allocate memory
    divT = self.dimension eq '2D' ? fltarr(dims[0:1],2) : fltarr(dims[0:1],3)
    
    ;Component 1
    divT[1:dims[0]-2,*,0] = ((*T)[2:dims[0]-1,*,0] - (*T)[0:dims[0]-3,*,0]) + $
                            ((*T)[2:dims[0]-1,*,1] - (*T)[0:dims[0]-3,*,1]) + $
                            ((*T)[2:dims[0]-1,*,2] - (*T)[0:dims[0]-3,*,2])
    
    ;Component 2
    divT[1:dims[0]-2,*,1] = ((*T)[*,2:dims[1]-1,1] - (*T)[*,0:dims[1]-3,1]) + $
                            ((*T)[*,2:dims[1]-1,3] - (*T)[*,0:dims[1]-3,3]) + $
                            ((*T)[*,2:dims[1]-1,4] - (*T)[*,0:dims[1]-3,4])
    
    ;Divide by the proper width.idl
    case self.orientation of
        'XY': begin
            divT[*,*,0] /= (2.0 * dx_de)
            divT[*,*,1] /= (2.0 * dy_de)
        endcase
        
        'XZ': begin
            divT[*,*,0] /= (2.0 * dx_de)
            divT[*,*,1] /= (2.0 * dz_de)
        endcase
        
        'YZ': begin
            divT[*,*,0] /= (2.0 * dy_de)
            divT[*,*,1] /= (2.0 * dz_de)
        endcase
    endcase
    
    ;Return
    case 1 of
        vec: return, divT
        x:   return, divT[*,*,0]
        z:   return, divT[*,*,1]
    endcase
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
function MrSim2_Data::Tensor_Par, T, V
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_T) then ptr_free, _T
        if ptr_valid(_V) then ptr_free, _V
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Was a name or data given?
    if MrIsA(T, 'STRING') $
        then _T = ptr_new(self -> GetData(T)) $
        else _T = ptr_new(T)
        
    if MrIsA(V, 'STRING') $
        then _V = ptr_new(self -> GetData(V)) $
        else _V = ptr_new(V)
    
    ;Get Data
    v_hat = _V / sqrt(total(*_V^2, 3))
    ptr_free, _V
    
    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    T_par = _T[*,*,0] * v_hat[*,*,0]^2 + 2.0 * _T[*,*,1] * v_hat[*,*,0] * v_hat[*,*,1] + 2.0 * _T[*,*,2] * v_hat[*,*,0] * v_hat[*,*,2] + $
            _T[*,*,3] * v_hat[*,*,1]^2 + 2.0 * _T[*,*,4] * v_hat[*,*,1] * v_hat[*,*,2] + $
            _T[*,*,5] * v_hat[*,*,2]^2
    ptr_free, _T
    
    return, T_par
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
function MrSim2_Data::Tensor_Perp, T, V
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_T) then ptr_free, _T
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Was a name or data given?
    T_par = self -> Tensor_Par(T, V)    
    
    if MrIsA(T, 'STRING') $
        then _T = ptr_new(self -> GetData(T, /DIAGONAL)) $
        else _T = ptr_new(T)

    ;Perpendicular pressure
    ;   - Assuming the pressure tensor is mostly gyrotropic, P = P_par + 2 P_perp
    T_perp = (_T[*,*,0] + _T[*,*,0] + _T[*,*,0] - T_par) / 2.0
    
    ;Free the pointers
    ptr_free, _T
    
    return, T_perp
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
    if n_elements(*self.Ay) eq 0 then self -> SetData, 'Ay'
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
    if x then if n_elements(*self.Bx) eq 0 then self -> SetData, 'Bx'
    if y then if n_elements(*self.By) eq 0 then self -> SetData, 'By'
    if z then if n_elements(*self.Bz) eq 0 then self -> SetData, 'Bz'
    
    ;Return the proper quantity
    case 1 of
        magnitude: return, self -> Magnitude(    'B')
        tf_cross:  return, self -> CrossProduct( 'B', crossName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_dot:    return, self -> DotProduct(   'B', dotName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_par:    return, self -> Parallel(     'B', parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Perpendicular('B', perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
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
    if x then if n_elements(*self.Ex) eq 0 then self -> SetData, 'Ex'
    if y then if n_elements(*self.Ey) eq 0 then self -> SetData, 'Ey'
    if z then if n_elements(*self.Ez) eq 0 then self -> SetData, 'Ez'
    
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
    if n_elements(*self.n_i) eq 0 then self -> SetData, 'ni'
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
    if n_elements(*self.n_e) eq 0 then self -> SetData, 'ne'
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
DIAGONAL=diagonal, $
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
    diagonal  = keyword_set(diagonal)
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
    tensor = tensor || (diagonal + Txx + Txy + Txz + Tyy + Tyz + Tzz + tf_par + tf_perp + magnitude) eq 0
    if diagonal then begin
        Txx = 1
        Tyy = 1
        Tzz = 1
    endif
    if tensor then begin
        Txx = 1
        Txy = 1
        Txz = 1
        Tyy = 1
        Tyz = 1
        Tzz = 1
    endif
    
    ;Tensor components
    if Txx then if n_elements(*self.Pe_xx) eq 0 then self -> SetData, 'Pe_xx'
    if Txy then if n_elements(*self.Pe_xy) eq 0 then self -> SetData, 'Pe_xy'
    if Txz then if n_elements(*self.Pe_xz) eq 0 then self -> SetData, 'Pe_xz'
    if Tyy then if n_elements(*self.Pe_yy) eq 0 then self -> SetData, 'Pe_yy'
    if Tyz then if n_elements(*self.Pe_yz) eq 0 then self -> SetData, 'Pe_yz'
    if Tzz then if n_elements(*self.Pe_zz) eq 0 then self -> SetData, 'Pe_zz'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        tf_par:    return, self -> Tensor_Par('Pe', parName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Tensor_Perp('Pe', perpName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude: return, self -> Tensor_Mag('Pe')
        Txx:       return, *self.Pe_xx
        Txy:       return, *self.Pe_xy
        Txz:       return, *self.Pe_xz
        Tyy:       return, *self.Pe_yy
        Tyz:       return, *self.Pe_yz
        Tzz:       return, *self.Pe_zz
        diagonal: begin
            ;Allocate memory
            dims = size(*self.Pe_xx, /DIMENSIONS)
            data = fltarr(dims[0], dims[1], 3)
            
            ;Return the diagonal
            data[*,*,0] = *self.Pe_xx
            data[*,*,1] = *self.Pe_yy
            data[*,*,2] = *self.Pe_zz
            return, data
        endcase
        tensor: begin
            ;Allocate memory
            dims = size(*self.Pe_xx, /DIMENSIONS)
            data = fltarr(dims[0], dims[1], 6)
            
            ;Create the data product
            data[*,*,0] = *self.Pe_xx
            data[*,*,1] = *self.Pe_xy
            data[*,*,2] = *self.Pe_xz
            data[*,*,3] = *self.Pe_yy
            data[*,*,4] = *self.Pe_yz
            data[*,*,5] = *self.Pe_zz
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
DIAGONAL=diagonal, $
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
    diagonal  = keyword_set(diagonal)
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
    tensor = tensor || (diagonal + Txx + Txy + Txz + Tyy + Tyz + Tzz + tf_par + tf_perp + magnitude) eq 0
    if diagonal then begin
        Txx = 1
        Tyy = 1
        Tzz = 1
    endif
    if tensor then begin
        Txx = 1
        Txy = 1
        Txz = 1
        Tyy = 1
        Tyz = 1
        Tzz = 1
    endif
    
    ;Tensor components
    if Txx then if n_elements(*self.Pi_xx) eq 0 then self -> SetData, 'Pi_xx'
    if Txy then if n_elements(*self.Pi_xy) eq 0 then self -> SetData, 'Pi_xy'
    if Txz then if n_elements(*self.Pi_xz) eq 0 then self -> SetData, 'Pi_xz'
    if Tyy then if n_elements(*self.Pi_yy) eq 0 then self -> SetData, 'Pi_yy'
    if Tyz then if n_elements(*self.Pi_yz) eq 0 then self -> SetData, 'Pi_yz'
    if Tzz then if n_elements(*self.Pi_zz) eq 0 then self -> SetData, 'Pi_zz'
    
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
        diagonal: begin
            ;Allocate memory
            dims = size(*self.Pi_xx, /DIMENSIONS)
            data = fltarr(dims[0], dims[1], 3)
            
            ;Return the diagonal
            data[*,*,0] = *self.Pi_xx
            data[*,*,1] = *self.Pi_yy
            data[*,*,2] = *self.Pi_zz
            return, data
        endcase
        tensor: begin
            ;Allocate memory
            dims = size(*self.Pi_xx, /DIMENSIONS)
            data = fltarr(dims[0], dims[1], 3, 3)
            
            ;Create the data product
            data[*,*,0] = *self.Pi_xx
            data[*,*,1] = *self.Pi_xy
            data[*,*,2] = *self.Pi_xz
            data[*,*,3] = *self.Pi_yy
            data[*,*,4] = *self.Pi_yz
            data[*,*,5] = *self.Pi_zz
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
    if x then if n_elements(*self.Uix) eq 0 then self -> SetData, 'Uix'
    if y then if n_elements(*self.Uiy) eq 0 then self -> SetData, 'Uiy'
    if z then if n_elements(*self.Uiz) eq 0 then self -> SetData, 'Uiz'
    
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
    
    if x then if n_elements(*self.Uex) eq 0 then self -> SetData, 'Uex'
    if y then if n_elements(*self.Uey) eq 0 then self -> SetData, 'Uey'
    if z then if n_elements(*self.Uez) eq 0 then self -> SetData, 'Uez'
    
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
function MrSim2_Data::An_e
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
function MrSim2_Data::A0_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx = self -> GetData('Bx')
    By = self -> GetData('By')
    Bz = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe_xx')
    Pe_xy = self -> GetData('Pe_xy')
    Pe_xz = self -> GetData('Pe_xz')
    Pe_yy = self -> GetData('Pe_yy')
    Pe_yz = self -> GetData('Pe_yz')
    Pe_zz = self -> GetData('Pe_zz')
    
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
function MrSim2_Data::Dng_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Gather data
    Bx    = self -> GetData('Bx')
    By    = self -> GetData('By')
    Bz    = self -> GetData('Bz')
    Pe_xx = self -> GetData('Pe_xx')
    Pe_xy = self -> GetData('Pe_xy')
    Pe_xz = self -> GetData('Pe_xz')
    Pe_yy = self -> GetData('Pe_yy')
    Pe_yz = self -> GetData('Pe_yz')
    Pe_zz = self -> GetData('Pe_zz')
    
    ;Magnitude and direction of magnetic field
    Bmag   = sqrt(Bx^2 + By^2 + Bz^2)
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
    on_error, 2
	
	;Get the data
	;   - Number density
	;   - Bulk velocity
    n_e = self -> GetData('ne')
    n_i = self -> GetData('ni')
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, X=x, Y=y, Z=z, VECTOR=vec, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, X=x, Y=y, Z=z, VECTOR=vec, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)

    ;Total Current
    J = (-temporary(n_e) * temporary(Ue)) + $
        ( temporary(n_i) * temporary(Ui))

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
    Je = -1.0 * temporary(n_e) * temporary(Ue)
    
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
    
    return, V
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
              answer:    ptr_new(), $
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