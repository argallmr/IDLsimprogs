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
;   The purpose of this program is to calculate the electron plasma beta::
;       \beta_{e} = \frac{Tr(P_{e})} {B^{2} / {2 \mu_{0}}
;
; :Private:
;
; :Returns:
;       BETA_E:                 The electron beta.
;-
function MrSim2_Data::A, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Ay is the only component
    if perpendicular + parallel + x + z + magnitude gt 1 $
        then message, 'Only the y-component of the vector potential is available.'
    
    ;Ay
    if n_elements(*self.Ay) eq 0 then self -> ReadData, 'Ay'
    return, self.Ay
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
function MrSim2_Data::B, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    vec = (x + y + z + perpendicular + parallel + magnitude) eq 0
    if vec then begin
        x = 1
        y = 1
        z = 1
    endif
    
    ;B cannot be perpendicular to B
    ;   - Keep in mind that the direction parallel to B is B
    if perpendicular then message, 'The direction of B perpendicular to B is 0 always.'
    
    ;Bx
    if magnitude || parallel || x $
        then if n_elements(*self.Bx) eq 0 then self -> ReadData, 'Bx'
    
    ;By
    if magnitude || parallel || y $
        then if n_elements(*self.By) eq 0 then self -> ReadData, 'By'
    
    ;Bz
    if magnitude || parallel || z $
        then if n_elements(*self.Bz) eq 0 then self -> ReadData, 'Bz'
    
    ;Magnitude
    if magnitude || parallel then Bmag = sqrt(*self.Bx^2 + *self.By^2 + *self.Bz^2)
    
    ;Return the proper quantity
    case 1 of
        vec:       return, [[[*self.Bx]], [[*self.By]], [[*self.Bz]]]
        magnitude: return, Bmag
        parallel:  return, Bmag
        x:         return, *self.Bx
        y:         return, *self.By
        z:         return, *self.Bz
        else:      ;Do nothing
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
function MrSim2_Data::E, $
DESTROY=destroy, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    
    vec = (x + y + z + perpendicular + parallel + magnitude) eq 0
    if vec then begin
        x = 1
        y = 1
        z = 1
    endif
    
    ;Ex
    if magnitude || parallel || perpendicular || x $
        then if n_elements(*self.Ex) eq 0 then self -> ReadData, 'Ex'
    
    ;Ey
    if magnitude || parallel || perpendicular || y $
        then if n_elements(*self.Ey) eq 0 then self -> ReadData, 'Ey'
    
    ;Ez
    if magnitude || parallel || perpendicular || z $
        then if n_elements(*self.Ez) eq 0 then self -> ReadData, 'Ez'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Parallel(Ex, Ey, Ez, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Perpendicular(Ex, Ey, Ez, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, sqrt(*self.Ex^2 + *self.Ey^2 + *self.Ez^2)
        x:             return, *self.Ex
        y:             return, *self.Ey
        z:             return, *self.Ez
        else:          ;Not possible
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
function MrSim2_Data::n_i
    compile_opt strictarr, hidden
    on_error, 2
    
    ;ni
    if n_elements(*self.n_i) eq 0 then self -> ReadData, 'ni'
    return, *self.n_i
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
function MrSim2_Data::n_e
    compile_opt strictarr, hidden
    on_error, 2
    
    ;ne
    if n_elements(*self.n_e) eq 0 then self -> ReadData, 'ne'
    return, *self.n_e
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
function MrSim2_Data::Pe, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
TXX=Txx, $
TXY=Txy, $
TXZ=Txz, $
TYY=Tyy, $
TYZ=Tyz, $
TZZ=Tzz

    compile_opt strictarr, hidden
    on_error, 2
    
    ;Pe-xx
    if magnitude || parallel || perpendicular || Txx $
        then if n_elements(*self.Pe_xx) eq 0 then self -> ReadData, 'Pe-xx'
    
    ;Pe-xy
    if magnitude || parallel || perpendicular || Txy $
        then if n_elements(*self.Pe_xy) eq 0 then self -> ReadData, 'Pe-xy'
    
    ;Pe-xz
    if magnitude || parallel || perpendicular || Txz $
        then if n_elements(*self.Pe_xz) eq 0 then self -> ReadData, 'Pe-xz'
    
    ;Pe-yy
    if magnitude || parallel || perpendicular || Tyy $
        then if n_elements(*self.Pe_yy) eq 0 then self -> ReadData, 'Pe-yy'
    
    ;Pe-yz
    if magnitude || parallel || perpendicular || Tyz $
        then if n_elements(*self.Pe_yz) eq 0 then self -> ReadData, 'Pe-yz'
    
    ;Pe-zz
    if magnitude || parallel || perpendicular || Tzz $
        then if n_elements(*self.Pe_zz) eq 0 then self -> ReadData, 'Pe-zz'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Tensor_Par(Pe_xx, Pe_xy, Pe_xz, Pe_yy, Pe_yz, Pe_zz, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Tensor_Perp(Pe_xx, Pe_xy, Pe_xz, Pe_yy, Pe_yz, Pe_zz, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, (*self.Pe_xx + *self.Pe_yy + *self.Pe_zz) / 3.0
        Txx:           return, *self.Pe_xx
        Txy:           return, *self.Pe_xy
        Txz:           return, *self.Pe_xz
        Tyy:           return, *self.Pe_yy
        Tyz:           return, *self.Pe_yz
        Tzz:           return, *self.Pe_zz
        else:          ;Not possible
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
function MrSim2_Data::Pi, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
TXX=Txx, $
TXY=Txy, $
TXZ=Txz, $
TYY=Tyy, $
TYZ=Tyz, $
TZZ=Tzz

    compile_opt strictarr, hidden
    on_error, 2
    
    ;Pi-xx
    if magnitude || parallel || perpendicular || Txx $
        then if n_elements(*self.Pi_xx) eq 0 then self -> ReadData, 'Pi-xx'
    
    ;Pi-xy
    if magnitude || parallel || perpendicular || Txy $
        then if n_elements(*self.Pi_xy) eq 0 then self -> ReadData, 'Pi-xy'
    
    ;Pi-xz
    if magnitude || parallel || perpendicular || Txz $
        then if n_elements(*self.Pi_xz) eq 0 then self -> ReadData, 'Pi-xz'
    
    ;Pi-yy
    if magnitude || parallel || perpendicular || Tyy $
        then if n_elements(*self.Pi_yy) eq 0 then self -> ReadData, 'Pi-yy'
    
    ;Pi-yz
    if magnitude || parallel || perpendicular || Tyz $
        then if n_elements(*self.Pi_yz) eq 0 then self -> ReadData, 'Pi-yz'
    
    ;Pi-zz
    if magnitude || parallel || perpendicular || Tzz $
        then if n_elements(*self.Pi_zz) eq 0 then self -> ReadData, 'Pi-zz'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Tensor_Par(Pi_xx, Pi_xy, Pi_xz, Pi_yy, Pi_yz, Pi_zz, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Tensor_Pirp(Pi_xx, Pi_xy, Pi_xz, Pi_yy, Pi_yz, Pi_zz, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, (*self.Pi_xx + *self.Pi_yy + *self.Pi_zz) / 3.0
        Txx:           return, *self.Pi_xx
        Txy:           return, *self.Pi_xy
        Txz:           return, *self.Pi_xz
        Tyy:           return, *self.Pi_yy
        Tyz:           return, *self.Pi_yz
        Tzz:           return, *self.Pi_zz
        else:          ;Not possible
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
function MrSim2_Data::Ui, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Uix
    if magnitude || parallel || perpendicular || x $
        then if n_elements(*self.Uix) eq 0 then self -> ReadData, 'Uix'
    
    ;Uiy
    if magnitude || parallel || perpendicular || y $
        then if n_elements(*self.Uiy) eq 0 then self -> ReadData, 'Uiy'
    
    ;Uiz
    if magnitude || parallel || perpendicular || z $
        then if n_elements(*self.Uiz) eq 0 then self -> ReadData, 'Uiz'
    
    ;Which to return?
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Parallel(*self.Uix, *self.Uiy, *self.Uiz, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Perpendicular(*self.Uix, *self.Uiy, *self.Uiz, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, sqrt(*self.Uix^2 + *self.Uiy^2 + *self.Uiz^2)
        x:             return, *self.Uix
        y:             return, *self.Uiy
        z:             return, *self.Uiz
        else:          ;Not possible
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
function MrSim2_Data::Ue, $
MAGNITUDE=magnitude, $
PARALLEL=parallel, $
PERPENDICULAR=perpendicular, $
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
    on_error, 2
    
    ;Uex
    if magnitude || parallel || perpendicular || x then begin
        then if n_elements(*self.Uex) eq 0 then self -> ReadData, 'Uex'
    
    ;Uey
    if magnitude || parallel || perpendicular || y then begin
        then if n_elements(*self.Uey) eq 0 then self -> ReadData, 'Uey'
    
    ;Uez
    if magnitude || parallel || perpendicular || z then begin
        then if n_elements(*self.Uez) eq 0 then self -> ReadData, 'Uez'
    
    ;Return
    ;   - The order is crucial.
    case 1 of
        parallel:      return, self -> Parallel(*self.Uex, *self.Uey, *self.Uez, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        perpendicular: return, self -> Perpendicular(*self.Uex, *self.Uey, *self.Uez, X=x, Y=y, Z=z, MAGNITUDE=magnitude)
        magnitude:     return, sqrt(*self.Uex^2 + *self.Uey^2 + *self.Uez^2)
        x:             return, *self.Uex
        y:             return, *self.Uey
        z:             return, *self.Uez
        else:          ;Not possible
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
X=x, $
Y=y, $
Z=z
    compile_opt strictarr, hidden
	
	;Get the data
	;   - Number density
	;   - Bulk velocity
    n_e = self -> GetData('ne')
    n_e = self -> GetData('ni')
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)
    Ui  = self -> GetData('Ui', MAGNITUDE=magnitude, X=x, Y=y, Z=z, $
                                PARALLEL=parallel, PERPENDICULAR=perpendicular)

    ;Total Current
    J = (-n_e * Ue) + (n_i * Ui)
    
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
X=x, $
Y=y, $
Z=z
	compile_opt strictarr, hidden
	on_error, 2
	
	;Number density
	;   - ne is always stored as an object property. Retrieve the pointer
	;   - Only the X, Y, and Z components of Ue are object properties. Since we
	;       do not know which was chosen without lots of checks, retrieve the data.
    n_e = self -> GetData('ne')
    Ue  = self -> GetData('Ue', MAGNITUDE=magnitude, X=x, Y=y, Z=z, $
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
    
    return, Vx
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
function MrSim2_Data::E_dot_J
    compile_opt strictarr, hidden
    on_error, 2

    ;Get the electric field and current density data
    E_dot_J = self -> DotProduct('E', 'J', X=x, Y=y, Z=y)
    return, E_dot_J
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
;   The purpose of this program is to calculate the x-component of the JxB portion
;   of the Generalized Ohm's Law::
;       (J_{e} \times B)_{x} = \frac{1} {n_{e}} (J_{e,y} B_{z} - J_{e,z} B_{y})
;
; :Private:
;
; :Returns:
;       JexB_x:             The x-component of the JxB force.
;-
function MrSim2_Data::JexB_x
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
function MrSim2_Data::JexB_y
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
function MrSim2_Data::JexB_z
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
;   The purpose of this program is to calculate the x-component of the JxB portion
;   of the Generalized Ohm's Law::
;       (J_{i} \times B)_{x} = \frac{1} {n_{i}} (J_{i,y} B_{z} - J_{i,z} B_{y})
;
; :Private:
;
; :Returns:
;       JixB_x:             The x-component of the JxB force.
;-
function MrSim2_Data::JixB_x
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
function MrSim2_Data::JixB_y
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
function MrSim2_Data::JixB_z
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
function MrSim2_Data::Jpar
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
function MrSim2_Data::JxB_x
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
function MrSim2_Data::JxB_y
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
function MrSim2_Data::JxB_z
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
function MrSim2_Data::UexB_mag
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
function MrSim2_Data::UexB_x
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
function MrSim2_Data::UexB_y
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
function MrSim2_Data::UexB_z
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
function MrSim2_Data::UixB_mag
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
function MrSim2_Data::UixB_x
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
function MrSim2_Data::UixB_y
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
function MrSim2_Data::UixB_z
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
function MrSim2_Data::VxB_x
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
function MrSim2_Data::VxB_y
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
function MrSim2_Data::VxB_z
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