; docformat = 'rst'
;
; NAME:
;    MrSim_GDA_Rename
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
; PURPOSE:
;+
;   Translate data product names to simulation coordinates.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Examples:
;   Try the main level program at the end of this document::
;       IDL> .r MrSim_GDAnames
;
; :Params:
;       NAME:               in, required, type=string
;                           Name of a GDA data product.
;       COORD_SYSTEM:       in, optional, type=string, default='SIMULATION'
;                           Coordinate system in which the data will reside after
;                               transformation. Options include::
;                                   'SIMULATION'
;                                   'MAGNETOPAUSE'
;                                   'MAGNETOTAIL'
;
; :Keywords:
;       MVA_FRAME:          in, optional, type=boolean, default=0
;                           If set and `COORD_SYSTEM` is "Magnetopause" or "Magnetotail",
;                               then the minimum variance coordinates N, M, and L will
;                               be assigned appropriately.
;
; :Author:
;    Matthew Argall::
;    University of New Hampshire
;    Morse Hall Room 113
;    8 College Road
;    Durham, NH 03824
;    matthew.argall@wildcats.unh.edu
;
; :History:
;    Modification History::
;       2014/11/12  -   Written by Matthew Argall
;-
function MrSim_GDA_Rename, name, coord_system, $
MVA_FRAME=mva_frame
    compile_opt strictarr
    on_error, 2
    
    ;Defaults
    _coord_system = n_elements(coord_system) eq 0 ? 'SIMULATION' : strupcase(coord_system)
    _mva_frame    = keyword_set(mva_frame)

    ;Coordinates in use now
    if _mva_frame then begin
        case _coord_system of
            'SIMULATION':   coords_in = ['L', 'M', 'N']
            'MAGNETOPAUSE': coords_in = ['N', 'M', 'L']
            'MAGNETOTAIL':  coords_in = ['L', 'M', 'N']
            else: message, 'Coordinate system "' + coord_system + '" not recognized.'
        endcase
    endif else coords_in = ['X', 'Y', 'Z']
    
    ;How current system relates to the simulation coordinate system:
    case _coord_system of
        'SIMULATION':   coords_out = ['x', 'y', 'z']
        'MAGNETOPAUSE': coords_out = ['z', 'y', 'x']
        'MAGNETOTAIL':  coords_out = ['x', 'y', 'z']
        else: message, 'Coordinate system "' + coord_system + '" not recognized.'
    endcase
            
    ;Regex: ^(any character)(Not coordinates)*(coordinate?)(coordinate?)$
    ;   - Careful of "ne" with "LMN" coordinates.
    jcoord = strjoin(coords_in)
    regex  = '(^.[^' + jcoord + ']*)([' + jcoord + ']?)([' + jcoord + ']?$)'
    
    ;Find the coordinates in the name
    namePos   = stregex(name, regex, /SUBEXP, /EXTRACT, /FOLD_CASE)
    nameParts = stregex(name, regex, /SUBEXP, /EXTRACT, /FOLD_CASE)
    new_name  = reform(nameParts[1,*])
    nameParts = strupcase(nameParts)
    
    ;Change coordinates.
    for i = 2, 3 do begin
        ix = where(nameParts[i,*] eq coords_in[0], nx)
        iy = where(nameParts[i,*] eq coords_in[1], ny)
        iz = where(nameParts[i,*] eq coords_in[2], nz)
        if nx gt 0 then new_name[ix] += coords_out[0]
        if ny gt 0 then new_name[iy] += coords_out[1]
        if nz gt 0 then new_name[iz] += coords_out[2]
    endfor
    
    ;Return a scalar
    if n_elements(new_name) eq 1 then new_name = new_name[0]
    return, new_name
end


;-----------------------------------------------------
; Main Level Program: IDL> .r MrSim_Rename \\\\\\\\\\\
;-----------------------------------------------------
;EXAMPLE 1 -- Magnetopause system
names = ['Ay', 'Bx', 'By', 'Bz', 'E-', 'Ex', 'Ey', 'Ez', 'ne', 'ni', $
         'Pe-xx', 'Pe-xy', 'Pe-xz', 'Pe-yx', 'Pe-yy', 'Pe-yz', $
         'Pe-zx', 'Pe-zy', 'Pe-zz', 'Pi-xx', 'Pi-xy', 'Pi-xz', $
         'Pi-yx', 'Pi-yy', 'Pi-yz', 'Pi-zx', 'Pi-zy', 'Pi-zz', $
         'Uex', 'Uey', 'Uez', 'Uix', 'Uiy', 'Uiz']

new_names = MrSim_GDA_Rename(names, 'Magnetopause')
print, '-------------------------------------'
print, 'Magnetopause Names:'
print, names, FORMAT='(5(3x, a5))'
print, 'Simulation Names:'
print, new_names, FORMAT='(5(3x, a5))'


;EXAMPLE 1 -- Magnetopause MVA system
names = ['Am', 'Bn', 'Bm', 'Bl', 'E-', 'En', 'Em', 'El', 'ne', 'ni', $
         'Pe-nn', 'Pe-nm', 'Pe-nl', 'Pe-mn', 'Pe-mm', 'Pe-ml', $
         'Pe-ln', 'Pe-lm', 'Pe-ll', 'Pi-nn', 'Pi-nm', 'Pi-nl', $
         'Pi-mn', 'Pi-mm', 'Pi-ml', 'Pi-ln', 'Pi-lm', 'Pi-ll', $
         'Uen', 'Uem', 'Uel', 'Uin', 'Uim', 'Uil']

new_names = MrSim_GDA_Rename(names, 'Magnetopause', /MVA_FRAME)
print, '-------------------------------------'
print, 'Magnetopause Names:'
print, names, FORMAT='(5(3x, a5))'
print, 'Simulation Names:'
print, new_names, FORMAT='(5(3x, a5))'
end