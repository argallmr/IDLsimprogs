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
;       2014/11/20  -   Keywords properly checked in many routines. Divergence of a tensor
;                           is computed correctly. - MRA
;       2015/04/04  -   Added the ::OuterProduct method. - MRA
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
;   A helper method for ::GetDataByName method. Used to retrieve GDA data and ensure it
;   is returned in the proper coordinate system.
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
function MrSim2_Data::GetGDAByName, name
    compile_opt strictarr
    on_error, 2

    ;Translate between coordinate systems
    ;   - NAME is given in the current coordinate system
    ;   - GDA names reference the simulation coordinate system
    ;   - Vector and Tensor products (without [xyz] subscripts) are not affected
    ;       o Must be re-ordered below.
    _name = MrSim_GDA_Rename(name, self.coord_system)

    ;Read and/or get the data.
    case strupcase(_name) of
        ;Vectors
        'B':      data = self -> B(/VECTOR)
        'E':      data = self -> E(/VECTOR)
        'UE':     data = self -> Ue(/VECTOR)
        'UI':     data = self -> Ui(/VECTOR)
        
        ;Tensors
        'PE':     data = self -> Pe(/TENSOR)
        'PI':     data = self -> Pi(/TENSOR)
        
        ;Components
        'AY':     return, self -> A(/Y)
        'BX':     return, self -> B(/X)
        'BY':     return, self -> B(/Y)
        'BZ':     return, self -> B(/Z)
        'E-':     return, *self.electrons
        'EX':     return, self -> E(/X)
        'EY':     return, self -> E(/Y)
        'EZ':     return, self -> E(/Z)
        'NE':     return, self -> n_i()
        'NI':     return, self -> n_e()
        'PE_XX':  return, self -> Pe(/TXX)
        'PE_XY':  return, self -> Pe(/TXY)
        'PE_XZ':  return, self -> Pe(/TXZ)
        'PE_YX':  return, self -> Pe(/TXY)
        'PE_YY':  return, self -> Pe(/TYY)
        'PE_YZ':  return, self -> Pe(/TYZ)
        'PE_ZX':  return, self -> Pe(/TXZ)
        'PE_ZY':  return, self -> Pe(/TYZ)
        'PE_ZZ':  return, self -> Pe(/TZZ)
        'PI_XX':  return, self -> Pi(/TXX)
        'PI_XY':  return, self -> Pi(/TXY)
        'PI_XZ':  return, self -> Pi(/TXZ)
        'PI_YX':  return, self -> Pi(/TXY)
        'PI_YY':  return, self -> Pi(/TYY)
        'PI_YZ':  return, self -> Pi(/TYZ)
        'PI_ZX':  return, self -> Pi(/TXZ)
        'PI_ZY':  return, self -> Pi(/TYZ)
        'PI_ZZ':  return, self -> Pi(/TZZ)
        'UEX':    return, self -> Ue(/X)
        'UEY':    return, self -> Ue(/Y)
        'UEZ':    return, self -> Ue(/Z)
        'UIX':    return, self -> Ui(/X)
        'UIY':    return, self -> Ui(/Y)
        'UIZ':    return, self -> Ui(/Z)
        else:     message, 'Not a GDA data product: "' + name + '".'
    endcase
    
    ;
    ;If we get to here, a tensor or vector was requested.
    ;   - Rearrange the components so that the match the current coordinate system
    ;       not the simulation coordinate system.
    ;
    
    ;Get the dimensions
    dims = size(data, /DIMENSIONS)

    ;Vector
    if dims[2] eq 3 then begin
        case self.coord_system of
            'SIMULATION':   ;Do nothing
            
            ; Vx -> Vz
            ; Vy -> Vy
            ; Vz -> Vx
            'MAGNETOPAUSE': data = data[*,*,[2,1,0]]
            'MAGNETOTAIL':  ;Do nothing
        endcase
        
    ;Tensor
    endif else if dims[2] eq 6 then begin
        case self.coord_system of
            'SIMULATION':   ;Do nothing
            
            ; Txx -> Tzz
            ; Txy -> Tzy = Tyz
            ; Txz -> Tzx = Txz
            ; Tyy -> Tyy
            ; Tyz -> Tyx = Txy
            ; Tzz -> Txx
            'MAGNETOPAUSE': data = data[*,*,[5,4,2,3,1,0]]
            'MAGNETOTAIL':  ;Do nothing
        endcase
    endif
    
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
    x    = keyword_set(x)
    y    = keyword_set(y)
    z    = keyword_set(z)
    Txx  = keyword_set(Txx)
    Txy  = keyword_set(Txy)
    Txz  = keyword_set(Txz)
    Tyy  = keyword_set(Tyy)
    Tyz  = keyword_set(Tyz)
    Tzz  = keyword_set(Tzz)
    dx   = keyword_set(dx)
    dy   = keyword_set(dy)
    dz   = keyword_set(dz)

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
    if nOps eq 0 then begin
        data = self -> GetDataByName(expr)

        ;Get the correct component
        case 1 of
            X:    return, data[*,*,0]
            Y:    return, data[*,*,1]
            Z:    return, data[*,*,2]
            TXX:  return, data[*,*,0]
            TXY:  return, data[*,*,1]
            TXZ:  return, data[*,*,2]
            TYY:  return, data[*,*,3]
            TYZ:  return, data[*,*,4]
            TZZ:  return, data[*,*,5]
            DX:   return, self -> Scalar_Derivative(data, /D_DX)
            DY:   return, self -> Scalar_Derivative(data, /D_DY)
            DZ:   return, self -> Scalar_Derivative(data, /D_DZ)
            else: return, data
        endcase
    endif

    ;Step through each operation
    for i = 0, nOps - 1 do begin
        lhs = self -> GetDataByName(ops[0,i])
        rhs = self -> GetDataByName(ops[2,i])
        
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
        case ops[1,i] of
            '+': *self.answer = lhs + rhs
            '-': *self.answer = lhs - rhs
            '*': *self.answer = lhs * rhs
            '/': *self.answer = lhs / rhs
            '^': *self.answer = lhs ^ rhs
            '.': *self.answer = self -> Vector_DotProduct(  lhs, rhs )
            'X': *self.answer = self -> Vector_CrossProduct(lhs, rhs )
            else: message, 'Operation invalid: "' + ops[1,i] + '".'
        endcase
    endfor
    
    ;Take a component of the answer?
    case 1 of
        x:   *self.answer = (*self.answer)[*,*,0]
        y:   *self.answer = (*self.answer)[*,*,1]
        z:   *self.answer = (*self.answer)[*,*,2]
        Txx: *self.answer = (*self.answer)[*,*,0]
        Txy: *self.answer = (*self.answer)[*,*,1]
        Txz: *self.answer = (*self.answer)[*,*,2]
        Tyy: *self.answer = (*self.answer)[*,*,3]
        Tyz: *self.answer = (*self.answer)[*,*,4]
        Tzz: *self.answer = (*self.answer)[*,*,5]
        DX:  *self.answer = self -> Scalar_Derivative(*self.answer, /D_DX)
        DY:  *self.answer = self -> Scalar_Derivative(*self.answer, /D_DY)
        DZ:  *self.answer = self -> Scalar_Derivative(*self.answer, /D_DZ)
        else: ;Do nothing
    endcase
    
    ;Return the result
    return, *self.answer
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
function MrSim2_Data::GetDataByName, name
    compile_opt strictarr
    on_error, 2

    ;If no parameters were given, print a list of available data products.
    if n_params() eq 0 then begin
        message, 'Use: data = mySim -> GetData(name)'
        MrSim -> ListProducts
        return, !Null
    endif

    ;Need to remove LMN coordinates
    if self.mva_frame $
        then theName = MrSim_GDA_Rename(name, self.coord_system, self.coord_system, /MVA_FRAME) $
        else theName = name

    ;Check to see if the data has already been read first.
    case strupcase(theName) of
        'ANS':    return, *self.answer
        'AY':     return, self -> GetGDAByName('Ay')
        'B':      return, self -> GetGDAByName('B')
        'BX':     return, self -> GetGDAByName('Bx')
        'BY':     return, self -> GetGDAByName('By')
        'BZ':     return, self -> GetGDAByName('Bz')
        'E-':     return, *self.electrons
        'E':      return, self -> GetGDAByName('E')
        'EX':     return, self -> GetGDAByName('Ex')
        'EY':     return, self -> GetGDAByName('Ey')
        'EZ':     return, self -> GetGDAByName('Ez')
        'NE':     return, self -> GetGDAByName('Ne')
        'NI':     return, self -> GetGDAByName('Ni')
        'PE':     return, self -> GetGDAByName('Pe')
        'PE_XX':  return, self -> GetGDAByName('Pe_xx')
        'PE_XY':  return, self -> GetGDAByName('Pe_xy')
        'PE_XZ':  return, self -> GetGDAByName('Pe_xz')
        'PE_YX':  return, self -> GetGDAByName('Pe_yx')
        'PE_YY':  return, self -> GetGDAByName('Pe_yy')
        'PE_YZ':  return, self -> GetGDAByName('Pe_yz')
        'PE_ZX':  return, self -> GetGDAByName('Pe_zx')
        'PE_ZY':  return, self -> GetGDAByName('Pe_zy')
        'PE_ZZ':  return, self -> GetGDAByName('Pe_zz')
        'Pi':     return, self -> GetGDAByName('Pi')
        'PI_XX':  return, self -> GetGDAByName('Pi_xx')
        'PI_XY':  return, self -> GetGDAByName('Pi_xy')
        'PI_XZ':  return, self -> GetGDAByName('Pi_xz')
        'PI_YX':  return, self -> GetGDAByName('Pi_yx')
        'PI_YY':  return, self -> GetGDAByName('Pi_yy')
        'PI_YZ':  return, self -> GetGDAByName('Pi_yz')
        'PI_ZX':  return, self -> GetGDAByName('Pi_zx')
        'PI_ZY':  return, self -> GetGDAByName('Pi_zy')
        'PI_ZZ':  return, self -> GetGDAByName('Pi_zz')
        'UE':     return, self -> GetGDAByName('Ue')
        'UEX':    return, self -> GetGDAByName('Uex')
        'UEY':    return, self -> GetGDAByName('Uey')
        'UEZ':    return, self -> GetGDAByName('Uez')
        'UI':     return, self -> GetGDAByName('Ui')
        'UIX':    return, self -> GetGDAByName('Uix')
        'UIY':    return, self -> GetGDAByName('Uiy')
        'UIZ':    return, self -> GetGDAByName('Uiz')
        
        ;Custom Data Products
        'EIE':     data = self -> EIe()
        'EIE_X':   data = self -> EIe(/X)
        'EIE_Y':   data = self -> EIe(/Y)
        'EIE_Z':   data = self -> EIe(/Z)
        'EI':      data = self -> EI()
        'EI_X':    data = self -> EI(/X)
        'EI_Y':    data = self -> EI(/Y)
        'EI_Z':    data = self -> EI(/Z)
        'EI_JJ':   data = self -> EI_JJ()
        'EI_JJ_X': data = self -> EI_JJ(/X)
        'EI_JJ_Y': data = self -> EI_JJ(/Y)
        'EI_JJ_Z': data = self -> EI_JJ(/Z)
        'EI_VJ':   data = self -> EI_VJ()
        'EI_VJ_X': data = self -> EI_VJ(/X)
        'EI_VJ_Y': data = self -> EI_VJ(/Y)
        'EI_VJ_Z': data = self -> EI_VJ(/Z)
        'EI_JV':   data = self -> EI_JV()
        'EI_JV_X': data = self -> EI_JV(/X)
        'EI_JV_Y': data = self -> EI_JV(/Y)
        'EI_JV_Z': data = self -> EI_JV(/Z)
        'DIVPE':   data = self -> Tensor_Divergence('Pe')
        'DIVPE_X': data = self -> Tensor_Divergence('Pe', /X)
        'DIVPE_Y': data = self -> Tensor_Divergence('Pe', /Y)
        'DIVPE_Z': data = self -> Tensor_Divergence('Pe', /Z)
        'DIVPI':   data = self -> Tensor_Divergence('Pi')
        'DIVPI_X': data = self -> Tensor_Divergence('Pi', /X)
        'DIVPI_Y': data = self -> Tensor_Divergence('Pi', /Y)
        'DIVPI_Z': data = self -> Tensor_Divergence('Pi', /Z)
        'J':       data = self -> J(/VECTOR)
        'JX':      data = self -> J(/X)
        'JY':      data = self -> J(/Y)
        'JZ':      data = self -> J(/Z)
        'JE':      data = self -> Je(/VECTOR)
        'JEX':     data = self -> Je(/X)
        'JEY':     data = self -> Je(/Y)
        'JEZ':     data = self -> Je(/Z)
        'JI':      data = self -> Ji(/VECTOR)
        'JIX':     data = self -> Ji(/X)
        'JIY':     data = self -> Ji(/Y)
        'JIZ':     data = self -> Ji(/Z)
        'JXB':     data = self -> Vector_CrossProduct('J', 'B')
        'JXB_X':   data = self -> Vector_CrossProduct('J', 'B', /X)
        'JXB_Y':   data = self -> Vector_CrossProduct('J', 'B', /Y)
        'JXB_Z':   data = self -> Vector_CrossProduct('J', 'B', /Z)
        'PE_PAR':  data = self -> Pe(/PARALLEL)
        'PE_PERP': data = self -> Pe(/PERPENDICULAR)
        'PE_TR':   data = self -> Pe(/TRACE)
        'PI_PAR':  data = self -> Pi(/PARALLEL)
        'PI_PERP': data = self -> Pi(/PERPENDICULAR)
        'PI_TR':   data = self -> Pi(/TRACE)
        'TE':      data = self -> Te(/TENSOR)
        'TE_PAR':  data = self -> Te(/PARALLEL)
        'TE_PERP': data = self -> Te(/PERPENDICULAR)
        'TE_MAG':  data = self -> Te(/MAGNITUDE)
        'TE_TR':   data = self -> Te(/TRACE)
        'TE_XX':   data = self -> Te(/TXX)
        'TE_XY':   data = self -> Te(/TXY)
        'TE_XZ':   data = self -> Te(/TXZ)
        'TE_YX':   data = self -> Te(/TXY)
        'TE_YY':   data = self -> Te(/TYY)
        'TE_YZ':   data = self -> Te(/TYZ)
        'TE_ZX':   data = self -> Te(/TXZ)
        'TE_ZY':   data = self -> Te(/TYZ)
        'TE_ZZ':   data = self -> Te(/TZZ)
        'TI':      data = self -> Ti(/TENSOR)
        'TI_PAR':  data = self -> Ti(/PARALLEL)
        'TI_PERP': data = self -> Ti(/PERPENDICULAR)
        'TI_MAG':  data = self -> Ti(/MAGNITUDE)
        'TI_TR':   data = self -> Ti(/TRACE)
        'TI_XX':   data = self -> Ti(/TXX)
        'TI_XY':   data = self -> Ti(/TXY)
        'TI_XZ':   data = self -> Ti(/TXZ)
        'TI_YX':   data = self -> Ti(/TXY)
        'TI_YY':   data = self -> Ti(/TYY)
        'TI_YZ':   data = self -> Ti(/TYZ)
        'TI_ZX':   data = self -> Ti(/TXZ)
        'TI_ZY':   data = self -> Ti(/TYZ)
        'TI_ZZ':   data = self -> Ti(/TZZ)
        'V':       data = self -> V(/VECTOR)
        'VX':      data = self -> V(/X)
        'VY':      data = self -> V(/Y)
        'VZ':      data = self -> V(/Z)
        'VXB':     data = self -> Vector_CrossProduct('V', 'B')
        'VXB_X':   data = self -> Vector_CrossProduct('V', 'B', /X)
        'VXB_Y':   data = self -> Vector_CrossProduct('V', 'B', /Y)
        'VXB_Z':   data = self -> Vector_CrossProduct('V', 'B', /Z)
        'A0_E':    data = self -> A0_e()
        'AN_E':    data = self -> An_e()
        'DNG_E':   data = self -> Dng_e()
        else: message, 'Data product not available: "' + name + '".'
    endcase
    
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
        'ANSWER': if n_elements(*self.answer)    gt 0 then tf_has = 1B
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
        'PI_YX':  if n_elements(*self.Pi_xy)     gt 0 then tf_has = 1B
        'PI_YY':  if n_elements(*self.Pi_yy)     gt 0 then tf_has = 1B
        'PI_YZ':  if n_elements(*self.Pi_yz)     gt 0 then tf_has = 1B
        'PI_ZX':  if n_elements(*self.Pi_xz)     gt 0 then tf_has = 1B
        'PI_ZY':  if n_elements(*self.Pi_yz)     gt 0 then tf_has = 1B
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
        'PE_YX':  *self.Pe_xy     = data
        'PE_YY':  *self.Pe_yy     = data
        'PE_YZ':  *self.Pe_yz     = data
        'PE_ZX':  *self.Pe_xz     = data
        'PE_ZY':  *self.Pe_yz     = data
        'PE_ZZ':  *self.Pe_zz     = data
        'PI_XX':  *self.Pi_xx     = data
        'PI_XY':  *self.Pi_xy     = data
        'PI_XZ':  *self.Pi_xz     = data
        'PI_YX':  *self.Pi_xy     = data
        'PI_YY':  *self.Pi_yy     = data
        'PI_YZ':  *self.Pi_yz     = data
        'PI_ZX':  *self.Pi_xz     = data
        'PI_ZY':  *self.Pi_yz     = data
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
;   Compute the derivative of a data scalar quantity.
;
; :Params:
;       DATA:           in, required, type=string/NxM float
;                       A string naming the data product or the actual data of which
;                           the derivative is to be taken.
;
; :Returns:
;       DER:            The derivative of `DATA`.
;-
function MrSim2_Data::Scalar_Derivative, data, $
D_DX=d_dx, $
D_DY=d_dy, $
D_DZ=d_dz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_data) then ptr_free, _data
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    d_dx = keyword_set(d_dx)
    d_dy = keyword_set(d_dy)
    d_dz = keyword_set(d_dz)
    if d_dx + d_dy + d_dz gt 1 then d_dx = 1
    if d_dx + d_dy + d_dz gt 1 then message, 'D_DX, D_DY, and D_DZ are mutually exclusive.'
    
    ;Get the data
    ;   - Avoid copying the data by creating a pointer
    if size(data, /TNAME) eq 'STRING' $
        then _data = ptr_new(self -> GetData(data)) $
        else _data = ptr_new(data)
    
    ;Ensure the data is a scalar quantity (i.e. not a vector or tensor)
    dims = size(*_data, /DIMENSIONS)
    if n_elements(dims) ne 2 then $
        message, 'Data must be 2D in order to take the derivative.'
    
    ;Get the grid spacing of the simulation in electron skin depths
    self -> GetInfo, DX_DE=dx_de, DY_DE=dy_de, DZ_DE=dz_de
    case self.coord_system of
        'SIMULATION':   begin
            dx = dx_de
            dy = dy_de
            dz = dz_de
        endcase
        'MAGNETOPAUSE': begin
            dx = -dz_de
            dy =  dy_de
            dz =  dx_de
        endcase
        'MAGNETOTAIL':  begin
            dx = dx_de
            dy = dy_de
            dz = dz_de
        endcase
    endcase
    
    ;Determine the orientation
    ;   'XY' -> data[X,Y]
    ;   'XZ' -> data[X,Z]
    ;   etc.
    dim1 = strmid(self.orientation, 0, 1)
    dim2 = strmid(self.orientation, 1, 2)
    
    ;Over which dimension will the derivative be taken?
    if d_dx then begin
        delta = dx
        dimension = dim1 eq 'X' ? 1 : $
                    dim2 eq 'X' ? 2 : 0
    endif else if d_dy then begin
        delta = dy
        dimension = dim1 eq 'Y' ? 1 : $
                    dim2 eq 'Y' ? 2 : 0
    endif else begin
        delta = dz
        dimension = dim1 eq 'Z' ? 1 : $
                    dim2 eq 'Z' ? 2 : 0
    endelse
    
    ;Allocate memory and take centered difference
    ;  - D/DY in 2D simulations is zero.
    ;  - !!! THIS DEPENDS ON A PROPERTY OF THE SUBCLASS MRSIM2__DEFINE !!!
    der = fltarr(dims)
    if ~(d_dy && self.dimension eq '2D') then begin
        case dimension of
            0: message, 'Derivative not compatible with orientation "' + self.orientation + '".'
            1: der[1:dims[0]-2,*] = ((*_data)[2:dims[0]-1,*] - (*_data)[0:dims[0]-3,*]) / (2.0 * delta)
            2: der[*,1:dims[1]-2] = ((*_data)[*,2:dims[1]-1] - (*_data)[*,0:dims[1]-3]) / (2.0 * delta)
        endcase
    endif
    
    ptr_free, _data
    return, der
end


;+
;   Compute the derivative of a data scalar quantity.
;
; :Params:
;       DATA:           in, required, type=string/NxM float
;                       A string naming the data product or the actual data of which
;                           the derivative is to be taken.
;
; :Returns:
;       DER:            The derivative of `DATA`.
;-
function MrSim2_Data::Scalar_Gradient, data, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_data) then ptr_free, _data
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    x = keyword_set(x)
    y = keyword_set(y)
    z = keyword_set(z)
    
    ;Take the derivatives
    ds_dx = self -> Scalar_Derivative(data, /D_DX)
    ds_dy = self -> Scalar_Derivative(data, /D_DY)
    ds_dz = self -> Scalar_Derivative(data, /D_DZ)
    
    ;Return the desired part
    case 1 of
        x:    return, ds_dx
        y:    return, ds_dy
        z:    return, ds_sz
        else: return, [[[temporary(ds_dx)]], [[temporary(ds_dy)]], [[temporary(ds_dz)]]]
    endcase
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
function MrSim2_Data::Tensor_Divergence, tensor, $
VECTOR=vec, $
X=x, $
Y=y, $
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
    y   = keyword_set(y)
    z   = keyword_set(z)
    if x + z eq 0 then vec = 1

    ;Was a name or data given?
    if MrIsA(tensor, 'STRING') $
        then T = ptr_new(self -> GetData(tensor)) $
        else T = ptr_new(tensor)

;-------------------------------------------------------
;Take Derivative ///////////////////////////////////////
;-------------------------------------------------------

    ;Allocate memory
    dims = size(*T, /DIMENSIONS)
    divT = fltarr([dims[0:1],3])
    
    ;
    ; This routine expects the tensor to be symmetic and for only the upper-diagonal
    ; elements to be passed in.
    ;   - T[*,*,0] = Txx
    ;   - T[*,*,1] = Txy
    ;   - T[*,*,2] = Txz
    ;   - T[*,*,3] = Tyy
    ;   - T[*,*,4] = Tyz
    ;   - T[*,*,5] = Tzz
    ;
    
    ;Component 1
    ;   - ( div Pe)_x = d/dx Txx + d/dy Tyx + d/dz Tzx
    ;   - d/dy = 0
    divT[0,0,0] = self -> Scalar_Derivative( (*T)[*,*,0], /D_DX) + $      ;Txx = Txx
                  self -> Scalar_Derivative( (*T)[*,*,2], /D_DZ)          ;Txz = Tzx
    
    ;Component 2
    ;   - ( div Pe)_y = d/dx Txy + d/dy Tyy + d/dz Tzy
    ;   - d/dy = 0
    divT[0,0,1] = self -> Scalar_Derivative( (*T)[*,*,1], /D_DX) + $      ;Txy
                  self -> Scalar_Derivative( (*T)[*,*,4], /D_DZ)          ;Tyz = Tzy
    
    ;Component 3
    ;   - ( div Pe)_z = d/dx Txz + d/dy Tyz + d/dz Tzz
    ;   - d/dy = 0
    divT[0,0,2] = self -> Scalar_Derivative( (*T)[*,*,2], /D_DX) + $      ;Txz
                  self -> Scalar_Derivative( (*T)[*,*,5], /D_DZ)          ;Tzz
    
    ;Return
    case 1 of
        vec: return, divT
        x:   return, divT[*,*,0]
        y:   return, difT[*,*,1]
        z:   return, divT[*,*,2]
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
function MrSim2_Data::Tensor_Parallel, T, V
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
    v_hat         = *_V
    v_hat[*,*,0] /= sqrt(total(*_V^2, 3))
    v_hat[*,*,1] /= sqrt(total(*_V^2, 3))
    v_hat[*,*,2] /= sqrt(total(*_V^2, 3))
    ptr_free, _V
    
    ;
    ; This routine expects the tensor to be symmetic and for only the upper-diagonal
    ; elements to be passed in.
    ;   - T[*,*,0] = Txx
    ;   - T[*,*,1] = Txy
    ;   - T[*,*,2] = Txz
    ;   - T[*,*,3] = Tyy
    ;   - T[*,*,4] = Tyz
    ;   - T[*,*,5] = Tzz
    ;

    ;Parallel pressure.
    ;   - Pressure tensor is symmetric, hence, the multiple of 2.0
    ;   - P_par = integral f (v dot b) (v dot b) d3v
    T_par = (*_T)[*,*,0] * v_hat[*,*,0]^2 + 2.0 * (*_T)[*,*,1] * v_hat[*,*,0] * v_hat[*,*,1] + 2.0 * (*_T)[*,*,2] * v_hat[*,*,0] * v_hat[*,*,2] + $
            (*_T)[*,*,3] * v_hat[*,*,1]^2 + 2.0 * (*_T)[*,*,4] * v_hat[*,*,1] * v_hat[*,*,2] + $
            (*_T)[*,*,5] * v_hat[*,*,2]^2
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
function MrSim2_Data::Tensor_Perpendicular, T, V
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
    T_par = self -> Tensor_Parallel(T, V)    
    
    if MrIsA(T, 'STRING') $
        then _T = ptr_new(self -> GetData(T, /DIAGONAL)) $
        else _T = ptr_new(T)
    
    ;
    ; This routine expects the tensor to be symmetic and for only the upper-diagonal
    ; elements to be passed in.
    ;   - T[*,*,0] = Txx
    ;   - T[*,*,1] = Txy
    ;   - T[*,*,2] = Txz
    ;   - T[*,*,3] = Tyy
    ;   - T[*,*,4] = Tyz
    ;   - T[*,*,5] = Tzz
    ;

    ;Perpendicular pressure
    ;   - Assuming the pressure tensor is mostly gyrotropic, P = P_par + 2 P_perp
    T_perp = ((*_T)[*,*,0] + (*_T)[*,*,3] + (*_T)[*,*,5] - T_par) / 2.0
    
    ;Free the pointers
    ptr_free, _T

    return, T_perp
end


;+
;   The purpose of this method is to read the ASCII "info" file relating to Bill
;   Daughton's simulations.
;
; :Private:
;
; :Params:
;       NAME:               in, required, type=string
;                           Name of the vector quantity for which a unit vector
;                               is to be found.
;
; :Keywords:
;       TIME:               in, optional, type=integer
;                           Time index at which to create the unit vector.
;       XRANGE:             in, optional, type=fltarr(2)
;                           X-range (in de) over which particle data is to be kept.
;       ZRANGE:             in, optional, type=fltarr(2)
;                           Z-range (in de) over which particle data is to be kept.
;
; :Returns:
;       DATA:               Electron particle data.
;-
function MrSim2_Data::Unit_Vector, name, $
TIME=time, $
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

    ;Set defaults
    if n_elements(time)   eq 0 then time   = self.time
    if n_elements(xrange) eq 0 then xrange = self.xrange
    if n_elements(zrange) eq 0 then zrange = self.zrange
    
    ;Read vector component data
    Vx = self -> ReadMoment(name + 'x', time, XRANGE=xrange, ZRANGE=zrange)
    Vy = self -> ReadMoment(name + 'y', time, XRANGE=xrange, ZRANGE=zrange)
    Vz = self -> ReadMoment(name + 'z', time, XRANGE=xrange, ZRANGE=zrange)
    
    ;Take the average over the given range
    Vx = mean(Bx)
    Vy = mean(By)
    Vz = mean(Bz)
    
    ;Magnitude
    Vmag = sqrt(Vx^2 + Vy^2 + Vz^2)
    
    ;Unit vector
    V_hat = [Vx, Vy, Vz] / Vmag
    
    return, V_hat
end


;+
;   The purpose of this method is to read the ASCII "info" file relating to Bill
;   Daughton's simulations.
;
; :Private:
;
; :Params:
;       NAME:               in, required, type=string
;                           Name of the vector quantity for which a unit vector
;                               is to be found.
;
; :Keywords:
;       TIME:               in, optional, type=integer
;                           Time index at which to create the unit vector.
;       XRANGE:             in, optional, type=fltarr(2)
;                           X-range (in de) over which particle data is to be kept.
;       ZRANGE:             in, optional, type=fltarr(2)
;                           Z-range (in de) over which particle data is to be kept.
;
; :Returns:
;       DATA:               Electron particle data.
;-
function MrSim2_Data::Unit_Perp1, $
B_HAT=B_hat, $
E_HAT=E_hat, $
TIME=time, $
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

    ;E and B unit vectors.
    B_hat = self -> Unit_Vector('B', TIME=time, XRANGE=xrange, ZRANGE=zrange)
    E_hat = self -> Unit_Vector('E', TIME=time, XRANGE=xrange, ZRANGE=zrange)
    
    ;ExB direction
    perp1  = [ E_hat[1]*B_hat[2] - E_hat[2]*B_hat[1], $
               E_hat[2]*B_hat[0] - E_hat[0]*B_hat[2], $
               E_hat[0]*B_hat[1] - E_hat[1]*B_hat[0] ]

    ;Unit vector
    p1_hat = perp1 / sqrt(total(perp1^2))
    
    return, p1_hat
end


;+
;   The purpose of this method is to read the ASCII "info" file relating to Bill
;   Daughton's simulations.
;
; :Private:
;
; :Params:
;       NAME:               in, required, type=string
;                           Name of the vector quantity for which a unit vector
;                               is to be found.
;
; :Keywords:
;       TIME:               in, optional, type=integer
;                           Time index at which to create the unit vector.
;       XRANGE:             in, optional, type=fltarr(2)
;                           X-range (in de) over which particle data is to be kept.
;       ZRANGE:             in, optional, type=fltarr(2)
;                           Z-range (in de) over which particle data is to be kept.
;
; :Returns:
;       DATA:               Electron particle data.
;-
function MrSim2_Data::Unit_Perp2, $
B_HAT=B_hat, $
E_HAT=E_hat, $
P1_HAT=p1_hat, $
TIME=time, $
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

    ;Perp1 and B unit vectors.
    p1_hat = self -> Unit_Perp1(B_HAT=B_hat, TIME=time, XRANGE=xrange, ZRANGE=zrange)
    
    ;Bx(ExB) direction
    perp2  = [ B_hat[1]*p1_hat[2] - B_hat[2]*p1_hat[1], $
               B_hat[2]*p1_hat[0] - B_hat[0]*p1_hat[2], $
               B_hat[0]*p1_hat[1] - B_hat[1]*p1_hat[0] ]

    ;Unit vector
    p2_hat = perp2 / sqrt(total(perp2^2))
    
    return, p2_hat
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
function MrSim2_Data::Vector_CrossProduct, v1, v2, $
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
    
    ;Defaults
    magnitude = keyword_set(magnitude)
    x = keyword_set(x)
    y = keyword_set(y)
    z = keyword_set(z)
    
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
function MrSim2_Data::Vector_Divergence, vec, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_V) then ptr_free, _V
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    x = keyword_set(x)
    y = keyword_set(y)
    z = keyword_set(z)
    
    ;Was a name or data given?
    if MrIsA(v1, 'STRING') $
        then _V = ptr_new(self -> GetData(vec)) $
        else _V = ptr_new(vec)
    
    ;Take the derivative of each component.
    dv_dx = self -> Scalar_Derivative((*_V)[*,*,0], /D_DX)
    dv_dy = self -> Scalar_Derivative((*_V)[*,*,1], /D_DY)
    dv_dz = self -> Scalar_Derivative((*_V)[*,*,2], /D_DZ)
    
    ;Divergence terms
    case 1 of
        x:    div = temporary(dv_dx)
        y:    div = temporary(dv_dy)
        z:    div = temporary(dv_dz)
        else: div = temporary(dv_dx) + temporary(dv_dy) + temporary(dv_dz)
    endcase
        
    ;Free the pointers
    ptr_free, _V
    
    return, div
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
function MrSim2_Data::Vector_DotProduct, v1, v2, $
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
function MrSim2_Data::Vector_Gradient, vec, $
TXX=Txx, $
TXY=Txy, $
TXZ=Txz, $
TYX=Tyx, $
TYY=Tyy, $
TYZ=Tyz, $
TZX=Tzx, $
TZY=Tzy, $
TZZ=Tzz
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if ptr_valid(_V) then ptr_free, _V
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Defaults
    Txx = keyword_set(Txx)
    Txy = keyword_set(Txy)
    Txz = keyword_set(Txz)
    Tyx = keyword_set(Tyx)
    Tyy = keyword_set(Tyy)
    Tyz = keyword_set(Tyz)
    Tzx = keyword_set(Tzx)
    Tzy = keyword_set(Tzy)
    Tzz = keyword_set(Tzz)
    
    ;Was a name or data given?
    if MrIsA(v1, 'STRING') $
        then _V = ptr_new(self -> GetData(vec)) $
        else _V = ptr_new(vec)
    
    ;
    ;Compute the outer product
    ;           | d/dx vx    d/dx vy    d/dx vz |
    ;  grad V = | d/dy vx    d/dy vy    d/dy vz | = T
    ;           | d/dz vx    d/dz vy    d/dz vz |
    ;
    ;This is an asymmetric matrix
    ;   - T[*,*,0] = Txx
    ;   - T[*,*,1] = Txy
    ;   - T[*,*,2] = Txz
    ;   - T[*,*,3] = Tyx
    ;   - T[*,*,4] = Tyy
    ;   - T[*,*,5] = Tyz
    ;   - T[*,*,6] = Tzx
    ;   - T[*,*,7] = Tzy
    ;   - T[*,*,8] = Tzz
    ;
    T = fltarr( [size(*_V), 9] )
    T[*,*,0] = self -> Scalar_Derivative((*_V)[*,*,0], /D_DX)
    T[*,*,1] = self -> Scalar_Derivative((*_V)[*,*,1], /D_DX)
    T[*,*,2] = self -> Scalar_Derivative((*_V)[*,*,2], /D_DX)
    T[*,*,3] = self -> Scalar_Derivative((*_V)[*,*,0], /D_DY)
    T[*,*,4] = self -> Scalar_Derivative((*_V)[*,*,1], /D_DY)
    T[*,*,5] = self -> Scalar_Derivative((*_V)[*,*,2], /D_DY)
    T[*,*,6] = self -> Scalar_Derivative((*_V)[*,*,0], /D_DZ)
    T[*,*,7] = self -> Scalar_Derivative((*_V)[*,*,1], /D_DZ)
    T[*,*,8] = self -> Scalar_Derivative((*_V)[*,*,2], /D_DZ)
        
    ;Free the pointers
    ptr_free, _V
    
    ;Return the proper term
    case 1 of
        Txx:  return, T[*,*,0]
        Txy:  return, T[*,*,1]
        Txz:  return, T[*,*,2]
        Tyx:  return, T[*,*,3]
        Tyy:  return, T[*,*,4]
        Tyz:  return, T[*,*,5]
        Tzx:  return, T[*,*,6]
        Tzy:  return, T[*,*,7]
        Tzz:  return, T[*,*,8]
        else: return, T
    endcase
end


;+
;   Compute the magnitude of a vector.
;
;   Calling Sequence:
;       data = oSim -> Vector_Magnitude('E')
;       data = oSim -> Vector_Magnitude(Ex, Ey, Ez)
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
function MrSim2_Data::Vector_Magnitude, vx, vy, vz
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
;   Take the outer product of two quantities.
;
; :Params:
;       V1:             in, required, type=string or NxMx3 float
;                       A vectory quantity or the name of the vectory quantity.
;       V2:             in, required, type=string or NxMx3 float
;                       A vectory quantity or the name of the vectory quantity.
;
; :Returns:
;       OUTERPRODUCT:   The outer product of `V1` and `V2`.
;-
function MrSim2_Data::Vector_OuterProduct, v1, v2
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
    
    ;Allocate memory
    dims = size(*_v1, /DIMENSIONS)
    outerProduct = fltarr( [dims[0:1], 6] )
    
    ;
    ;Compute the outer product
    ;              | v1x * v2x    v1y * v2x    v1z * v2x |
    ;  v1 (X) v2 = | v1x * v2y    v1y * v2y    v1z * v2y | = T
    ;              | v1x * v2z    v1y * v2z    v1z * v2z |
    ;
    ;This is a symmetric matrix, like the pressure tensor
    ;   - Include only the upper-diagonal terms
    ;   - T[*,*,0] = Txx
    ;   - T[*,*,1] = Txy
    ;   - T[*,*,2] = Txz
    ;   - T[*,*,3] = Tyy
    ;   - T[*,*,4] = Tyz
    ;   - T[*,*,5] = Tzz
    ;
    outerProduct[0,0,0] = (*_v1)[*,*,0] * (*_v2)[*,*,0]
    outerProduct[0,0,1] = (*_v1)[*,*,1] * (*_v2)[*,*,0]
    outerProduct[0,0,2] = (*_v1)[*,*,2] * (*_v2)[*,*,0]
    outerProduct[0,0,3] = (*_v1)[*,*,1] * (*_v2)[*,*,1]
    outerProduct[0,0,4] = (*_v1)[*,*,2] * (*_v2)[*,*,1]
    outerProduct[0,0,5] = (*_v1)[*,*,2] * (*_v2)[*,*,2]
        
    ;Free the pointers
    ptr_free, _v1
    ptr_free, _v2
    
    return, outerProduct
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
function MrSim2_Data::Vector_Parallel, v1, v2, $
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
function MrSim2_Data::Vector_Perpendicular, v1, v2, $
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
    v1_par  = self -> Vector_Parallel(v1, v2)
    
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
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    'B')
        tf_cross:  return, self -> Vector_CrossProduct( 'B', crossName, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_dot:    return, self -> Vector_DotProduct(   'B', dotName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_par:    return, self -> Vector_Parallel(     'B', parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular('B', perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
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
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    'E')
        tf_cross:  return, self -> Vector_CrossProduct( 'E', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   'E', dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     'E', parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular('E', perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
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
;       TRACE:              in, optional, type=boolean, default=0
;                           If set, the trace of the pressure tensor is returned.
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
TRACE=Ttrace, $
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
    Ttrace    = keyword_set(Ttrace)
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
    tensor = tensor || (diagonal + Ttrace + Txx + Txy + Txz + Tyy + Tyz + Tzz + tf_par + tf_perp + magnitude) eq 0
    if diagonal || Ttrace then begin
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
    ;   - TENSOR and DIAGONAL must come before T[XYZ][XYZ].
    ;   - TTRACE must come before DIAGONAL
    case 1 of
        tf_par:    return, self -> Tensor_Parallel('Pe', parName)
        tf_perp:   return, self -> Tensor_Perpendicular('Pe', perpName)
        magnitude: return, self -> Tensor_Magnitude('Pe')
        Ttrace:    return, ( *self.Pe_xx + *self.Pe_yy + *self.Pe_zz ) / 3.0
        diagonal:  return, [[[*self.Pe_xx]], [[*self.Pe_yy]], [[*self.Pe_zz]]]
        tensor:    return, [[[*self.Pe_xx]], [[*self.Pe_xy]], [[*self.Pe_xz]], $
                            [[*self.Pe_yy]], [[*self.Pe_yz]], [[*self.Pe_zz]]]
        Txx:       return, *self.Pe_xx
        Txy:       return, *self.Pe_xy
        Txz:       return, *self.Pe_xz
        Tyy:       return, *self.Pe_yy
        Tyz:       return, *self.Pe_yz
        Tzz:       return, *self.Pe_zz
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
;       TRACE:              in, optional, type=boolean, default=0
;                           If set, the trace of the pressure tensor is returned.
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
TRACE=Ttrace, $
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
    Ttrace    = keyword_set(Ttrace)
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
    tensor = tensor || (diagonal + Ttrace + Txx + Txy + Txz + Tyy + Tyz + Tzz + tf_par + tf_perp + magnitude) eq 0
    if diagonal || Ttrace then begin
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
    ;   - TENSOR and DIAGONAL must come before T[XYZ][XYZ].
    ;   - TTRACE must come before DIAGONAL
    case 1 of
        tf_par:    return, self -> Tensor_Pararallel('Pi')
        tf_perp:   return, self -> Tensor_Perpendicular('Pi')
        magnitude: return, self -> Tensor_Magnitude('Pi')
        Ttrace:    return, ( *self.Pe_xx + *self.Pe_yy + *self.Pe_zz ) / 3.0
        diagonal:  return, [[[*self.Pi_xx]], [[*self.Pi_yy]], [[*self.Pi_zz]]]
        tensor:    return, [[[*self.Pi_xx]], [[*self.Pi_xy]], [[*self.Pi_xz]], $
                            [[*self.Pi_yy]], [[*self.Pi_yz]], [[*self.Pi_zz]]]
        Txx:       return, *self.Pi_xx
        Txy:       return, *self.Pi_xy
        Txz:       return, *self.Pi_xz
        Tyy:       return, *self.Pi_yy
        Tyz:       return, *self.Pi_yz
        Tzz:       return, *self.Pi_zz
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
        magnitude: return, self -> Vector_Magnitude(    'Ui')
        tf_cross:  return, self -> Vector_CrossProduct( 'Ui', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   'Ui', dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     'Ui', parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular('Ui', perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
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
    tf_cross  = keyword_set(crossproduct)
    tf_dot    = keyword_set(dotproduct)
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
        magnitude: return, self -> Vector_Magnitude(    'Ue')
        tf_cross:  return, self -> Vector_CrossProduct( 'Ue', crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   'Ue', dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     'Ue', parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular('Ue', perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
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
;   Compute the spacial gradient of the inertial term in Generalized Ohm's law assuming
;   that electons make up the current density. In this case, cJV = cVJ = -cJJ, and two
;   of the three contributing terms cancel, leaving
;
;      E_{inert,e} = \frac{m_{e}} {e} \nabla \cdot \vec{v_{e}} \vec{v_{e}}
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
function MrSim2_Data::EIe, $
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
    
    ;Get electron velocity and total current information
    ve_ve = self -> Vector_OuterProduct('Ue', 'Ue')
    EIe   = self -> Tensor_Divergence(ve_ve)
    
    ;Return
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    EIe)
        tf_cross:  return, self -> Vector_CrossProduct( EIe, crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   EIe, dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     EIe, parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular(EIe, perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, EIe
        x:         return, EIe[*,*,0]
        y:         return, EIe[*,*,1]
        z:         return, EIe[*,*,2]
        else:      ;Not possible
    endcase
end


;+
;   Compute the spacial gradient of the inertial term in Generalized Ohm's law assuming
;   that electons make up the current density. In this case, cJV = cVJ = -cJJ, and two
;   of the three contributing terms cancel, leaving
;
;      E_{inert,e} = \frac{m_{e}} {e} \nabla \cdot \vec{v_{e}} \vec{v_{e}}
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
function MrSim2_Data::EIe_UeDotGradUe, $
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
    
    ;Get electron velocity and total current information
    Ue     = self -> getData('Ue')
    gradUe = self -> Scalar_Derivative('Ue')
    ve_ve = self -> Vector_OuterProduct('Ue', 'Ue')
    EIe   = self -> Tensor_Divergence(ve_ve)
    
    ;Return
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    EIe)
        tf_cross:  return, self -> Vector_CrossProduct( EIe, crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   EIe, dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     EIe, parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular(EIe, perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, EIe
        x:         return, EIe[*,*,0]
        y:         return, EIe[*,*,1]
        z:         return, EIe[*,*,2]
        else:      ;Not possible
    endcase
end


;+
;   The purpose of this program is to compute the inertial term of the
;   generalized Ohm's Law.
;
;     E_{i} = \frac{m_{e}} { e^{2} n_{e} } \frac{ \partial V_{i} J_{j} } { \partial x_{i} }
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
function MrSim2_Data::EI, $
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
    
    ;Get electron velocity and total current information
    EI_JV = self -> GetData('EI_JV')
    EI_VJ = self -> GetData('EI_VJ')
    EI_JJ = self -> GetData('EI_JJ')
    EI    = EI_JV + EI_VJ + EI_JJ
    
    ;Return
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    EI)
        tf_cross:  return, self -> Vector_CrossProduct( EI, crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   EI, dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     EI, parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular(EI, perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, EI
        x:         return, EI[*,*,0]
        y:         return, EI[*,*,1]
        z:         return, EI[*,*,2]
        else:      ;Not possible
    endcase
end


;+
;   The purpose of this program is to compute the VJ term in the inertial term of the
;   generalized Ohm's Law.
;
;     E_{i} = \frac{m_{e}} { e^{2} n_{e} } \frac{ \partial V_{i} J_{j} } { \partial x_{i} }
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
function MrSim2_Data::EI_VJ, $
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
    
    ;Get electron velocity and total current information
    ; - VJ_ij = Vi * Jj
    ; - div VJ_ij = d/dx_i VJ_ij
    VJ     = self -> Vector_OuterProduct('V', 'J')
    div_VJ = self -> Tensor_Divergence(VJ)
    
    ;Electron Density
    n_e = self -> GetData('ne')
    
    ;Compute
    ;  - me / (e^2 * ne) * { d/dxi * (Vi * Jj) }
    ;  - me = 1
    ;  - e  = 1
    Ei_vj = div_VJ / n_e
    
    ;Return
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    Ei_vj)
        tf_cross:  return, self -> Vector_CrossProduct( Ei_vj, crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   Ei_vj, dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     Ei_vj, parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular(Ei_vj, perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, Ei_vj
        x:         return, Ei_vj[*,*,0]
        y:         return, Ei_vj[*,*,1]
        z:         return, Ei_vj[*,*,2]
        else:      ;Not possible
    endcase
end


;+
;   The purpose of this program is to compute the JV term in the inertial term of the
;   generalized Ohm's Law.
;
;     E_{i} = \frac{m_{e}} { e^{2} n_{e} } \frac{ \partial J_{i} V_{j} } { \partial x_{i} }
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
function MrSim2_Data::EI_JV, $
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
    
    ;Get electron velocity and total current information
    ; - VJ_ij = Vi * Jj
    ; - div VJ_ij = d/dx_i VJ_ij
    JV     = self -> Vector_OuterProduct('J', 'V')
    div_JV = self -> Tensor_Divergence(JV)
    
    ;Electron Density
    n_e = self -> GetData('ne')
    
    ;Compute
    ;  - me / (e^2 * ne) * { d/dxi * (Vi * Jj) }
    ;  - me = 1
    ;  - e  = 1
    Ei_jv = div_JV / n_e
    
    ;Return
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    Ei_jv)
        tf_cross:  return, self -> Vector_CrossProduct( Ei_jv, crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   Ei_jv, dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     Ei_jv, parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular(Ei_jv, perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, Ei_jv
        x:         return, Ei_jv[*,*,0]
        y:         return, Ei_jv[*,*,1]
        z:         return, Ei_jv[*,*,2]
        else:      ;Not possible
    endcase
end


;+
;   The purpose of this program is to compute the VJ term in the inertial term of the
;   generalized Ohm's Law.
;
;     E_{i} = \frac{m_{e}} { e^{3} n_{e}^{2} } \frac{ - \partial V_{i} J_{j} } { \partial x_{i} }
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
function MrSim2_Data::EI_JJ, $
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
    
    ;Get electron velocity and total current information
    ; - VJ_ij = Vi * Jj
    ; - div VJ_ij = d/dx_i VJ_ij
    JJ     = self -> Vector_OuterProduct('J', 'J')
    div_JJ = self -> Tensor_Divergence(JJ)
    
    ;Electron Density
    n_e = self -> GetData('ne')
    
    ;Compute
    ;  - me / (e^2 * ne) * { d/dxi * -1 / (e * ne) * (Ji * Jj) }
    ;  - me = 1
    ;  - e  = 1
    Ei_jj = -div_JJ / n_e^2
    
    ;Return
    ;   - VECTOR must come before [XYZ].
    case 1 of
        magnitude: return, self -> Vector_Magnitude(    Ei_jj)
        tf_cross:  return, self -> Vector_CrossProduct( Ei_jj, crossName, X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_dot:    return, self -> Vector_DotProduct(   Ei_jj, dotName,   X=x, Y=y, Z=z, MAGNITUDE=magitude)
        tf_par:    return, self -> Vector_Parallel(     Ei_jj, parName,   X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        tf_perp:   return, self -> Vector_Perpendicular(Ei_jj, perpName,  X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        vec:       return, Ei_jj
        x:         return, Ei_jj[*,*,0]
        y:         return, Ei_jj[*,*,1]
        z:         return, Ei_jj[*,*,2]
        else:      ;Not possible
    endcase
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
    if size(Ue, /N_DIMENSIONS) eq 3 then begin
        dims = size(Ue, /DIMENSIONS)
        J    = fltarr(dims)
        J[0,0,0] = -n_e * Ue[*,*,0] + n_i * Ui[*,*,0]
        J[0,0,1] = -n_e * Ue[*,*,1] + n_i * Ui[*,*,1]
        J[0,0,2] = -n_e * Ue[*,*,2] + n_i * Ui[*,*,2]
    endif else begin
        J = (-temporary(n_e) * temporary(Ue)) + $
            ( temporary(n_i) * temporary(Ui))
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

    ;Total Current
    if size(Ue, /N_DIMENSIONS) eq 3 then begin
        dims                  = size(Ue, /DIMENSIONS)
        Je                    = fltarr(dims)
        Je[0,0,0] = -n_e * Ue[*,*,0]
        Je[0,0,1] = -n_e * Ue[*,*,1]
        Je[0,0,2] = -n_e * Ue[*,*,2]
    endif else begin
        Je = (-temporary(n_e) * temporary(Ue))
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
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	on_error, 2
	
	;Number density
    n_i = self -> GetData('ni')
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, VECTOR=vec, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    
    ;Ion Current
    if size(Ue, /N_DIMENSIONS) eq 3 then begin
        dims                  = size(Ui, /DIMENSIONS)
        Ji                    = fltarr(dims)
        Ji[0,0,0] = -n_i * Ui[*,*,0]
        Ji[0,0,1] = -n_i * Ui[*,*,1]
        Ji[0,0,2] = -n_i * Ui[*,*,2]
    endif else begin
        Ji = (-temporary(n_e) * temporary(Ue))
    endelse
    
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
VECTOR=vec, $
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	
	;Get the ion and electron bulk velocities and masses
	self -> GetInfo, MI_ME=mi_me
	me = 1.0
	mi = me * mi_me
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, VECTOR=vec, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, VECTOR=vec, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    
    ;Calculate the fluid velocity
    V = (mi*Ui + me*Ue) / (mi + me)
    
    return, V
end


;+
;   Return various data products associated with the electron temperature tensor.
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
function MrSim2_Data::Te, $
DIAGONAL=diagonal, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
TENSOR=tensor, $
TRACE=Ttrace, $
TXX=Txx, $
TXY=Txy, $
TXZ=Txz, $
TYY=Tyy, $
TYZ=Tyz, $
TZZ=Tzz

    compile_opt strictarr, hidden
    on_error, 2

    ;Get the electron density and pressure
    n_e = self -> n_e()
    Pe  = self -> Pe(DIAGONAL=diagonal, $
                     MAGNITUDE=magnitude, $
                     PARALLEL=parallel, $
                     PERPENDICULAR=perpendicular, $
                     TENSOR=tensor, $
                     TRACE=Ttrace, $
                     TXX=Txx, $
                     TXY=Txy, $
                     TXZ=Txz, $
                     TYY=Tyy, $
                     TYZ=Tyz, $
                     TZZ=Tzz)
    
    ;Calculate the temperature
    Te = temporary(Pe) / temporary(n_e)
    
    return, Te
end


;+
;   Return various data products associated with the ion temperature tensor.
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
function MrSim2_Data::Ti, $
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
    
    ;Get the electron density and pressure
    n_i = self -> n_i()
    Pi  = self -> Ti(DIAGONAL=diagonal, $
                     MAGNITUDE=magnitude, $
                     PARALLEL=parallel, $
                     PERPENDICULAR=perpendicular, $
                     TENSOR=tensor, $
                     TXX=Txx, $
                     TXY=Txy, $
                     TXZ=Txz, $
                     TYY=Tyy, $
                     TYZ=Tyz, $
                     TZZ=Tzz)
    
    ;Calculate the temperature
    Ti = temporary(Pi) / temporary(n_i)
    
    return, Ti
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
	UexB = self -> Vector_CrossProduct('Ue', 'B', X=x, Y=y, Z=z)
	
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
	UixB = self -> Vector_CrossProduct('Ui', 'B', X=x, Y=y, Z=z)
	
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
	
	;Defaults
	magnitude = keyword_set(magnitude)
	x = keyword_set(x)
	y = keyword_set(y)
	z = keyword_set(z)
	
	;Get the Electric and Magnetic Field data
	ExB  = self -> Vector_CrossProduct('E', 'B')
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