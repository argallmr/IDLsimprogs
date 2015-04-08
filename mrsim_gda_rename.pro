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
;   Translate data product names from one coordinate system to another.
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
;       CS_OUT:             in, optional, type=string, default='SIMULATION'
;                           Destination coordinate system.
;
; :Keywords:
;       MVA_FRAME:          in, optional, type=boolean, default=0
;                           If set and `COORD_SYSTEM` is "Magnetopause" or "Magnetotail",
;                               then the minimum variance coordinates N, M, and L will
;                               be assigned appropriately.
;       MVA_OUT:            in, optional, type=boolean, default=0
;                           If set and `NEW_NAME` will be labeled with the appropriate
;                               minimum variance coordinate.
;
; :Returns:
;       NEW_NAME:           The same as `NAME` but with coordinate labels converted to
;                                the equivalent labels in `CS_OUT`, with `MVA_OUT` taken
;                                into consideration.
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
;       2015/04/03  -   Added the CS_OUT parameter and MVA_OUT keywords. - MRA
;       2015/04/04  -   Reworked coordinate finding algorithm
;-
function MrSim_GDA_Rename, name, coord_system, cs_out, $
MVA_FRAME=mva_frame, $
MVA_OUT=mva_out
    compile_opt strictarr
    on_error, 2
    
    ;Defaults
    _cs_out       = n_elements(cs_out)       eq 0 ? 'SIMULATION' : strupcase(cs_out)
    _coord_system = n_elements(coord_system) eq 0 ? 'SIMULATION' : strupcase(coord_system)
    _mva_frame    = keyword_set(mva_frame)
    _mva_out      = keyword_set(mva_out)

;-----------------------------------------------------
; Input CS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
    if _mva_frame then begin
        case _coord_system of
            'SIMULATION':   coords_in = ['L', 'M', 'N']
            'MAGNETOPAUSE': coords_in = ['N', 'M', 'L']
            'MAGNETOTAIL':  coords_in = ['L', 'M', 'N']
            else: message, 'Coordinate system "' + coord_system + '" not recognized.'
        endcase
    endif else coords_in = ['X', 'Y', 'Z']
    
;-----------------------------------------------------
; Output in MVA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
    if _mva_out then begin
        case _cs_out of
            'SIMULATION': begin
                case _coord_system of
                    'SIMULATION':   coords_out = ['L', 'M', 'N']
                    'MAGNETOPAUSE': coords_out = ['N', 'M', 'L']
                    'MAGNETOTAIL':  coords_out = ['L', 'M', 'N']
                endcase
            endcase
        
            'MAGNETOPAUSE': begin
                case _coord_system of
                    'SIMULATION':   coords_out = ['L', 'M', 'N']
                    'MAGNETOPAUSE': coords_out = ['N', 'M', 'L']
                    'MAGNETOTAIL':  coords_out = ['L', 'M', 'N']
                endcase
            endcase
        
            'MAGNETOTAIL': begin
                case _coord_system of
                    'SIMULATION':   coords_out = ['L', 'M', 'N']
                    'MAGNETOPAUSE': coords_out = ['N', 'M', 'L']
                    'MAGNETOTAIL':  coords_out = ['L', 'M', 'N']
                endcase
            endcase
        endcase
    
;-----------------------------------------------------
; Output in XYZ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
    endif else begin
        case _cs_out of
            'SIMULATION': begin
                case _coord_system of
                    'SIMULATION':   coords_out = ['x', 'y', 'z']
                    'MAGNETOPAUSE': coords_out = ['z', 'y', 'x']
                    'MAGNETOTAIL':  coords_out = ['x', 'y', 'z']
                endcase
            endcase
        
            'MAGNETOPAUSE': begin
                case _coord_system of
                    'SIMULATION':   coords_out = ['z', 'y', 'x']
                    'MAGNETOPAUSE': coords_out = ['x', 'y', 'z']
                    'MAGNETOTAIL':  coords_out = ['z', 'y', 'x']
                endcase
            endcase
        
            'MAGNETOTAIL': begin
                case _coord_system of
                    'SIMULATION':   coords_out = ['x', 'y', 'z']
                    'MAGNETOPAUSE': coords_out = ['z', 'y', 'x']
                    'MAGNETOTAIL':  coords_out = ['x', 'y', 'z']
                endcase
            endcase
        endcase
    endelse

;-----------------------------------------------------
; Translate CS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
    ;Regex: Find any coordinate at the end of the name
    regex  = '[' + strjoin(coords_in) + ']$'

    ;Split into the part before the coordinate and the part after the coordinate.
    pre  = name
    post = replicate('', n_elements(name))
    
    ;There will be at most 2 coordinates
    for i = 0, 1 do begin
        ;Find the coordinates in the name
        pos    = stregex(pre, regex, /FOLD_CASE)
        iCoord = where(pos ne -1, nCoord)

        ;Convert the coordinate
        if nCoord gt 0 then begin
            ;Extract the coordinate
            coords      = strmid(pre[iCoord], transpose(pos[iCoord]))
            pre[iCoord] = strmid(pre[iCoord], 0, transpose(pos[iCoord]))

            ;Convert the coordinate
            iSortCoords = sort(coords_in)
            iNewCoord   = value_locate(coords_in[iSortCoords], strupcase(coords))
            newCoord    = coords_out[iSortCoords[iNewCoord]]

            ;Append coordinate to post
            post[iCoord] = newCoord + post[iCoord]
        endif
    endfor
    
    ;Combine the results
    new_name = pre + post
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