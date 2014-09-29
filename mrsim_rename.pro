; docformat = 'rst'
;
; NAME:
;    MrSim_RenameData
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
;   The purpose of this program is to rename a data product to reflect the coordinate
;   system in which it is displayed. Names are given in simulation coordinates (x,y,z)
;   and are renamed appropriately for the transformation they have to undergo to reach
;   the destination system.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       NAME:               in, required, type=string
;                           The name of the vector quantity to be plotted in color.
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
;       2014/06/03  -   Written by Matthew Argall
;-
function MrSim_Rename, name, coord_system, $
MVA_FRAME=mva_frame, $
SUBSCRIPT=subscript
    compile_opt strictarr
    on_error, 2
    
    ;Defaults
    new_name      = name
    _coord_system = n_elements(coord_system) eq 0 ? 'SIMULATION' : strupcase(coord_system)
    _mva_frame    = keyword_set(mva_frame)
    _subscript    = keyword_set(subscript)

    case _coord_system of
        'SIMULATION': ;Do nothing
        
        'MAGNETOPAUSE': begin
            ;Which system?
            coords = mva_frame ? ['N', 'M', 'L'] : ['x', 'y', 'z']
            if _subscript then coords = '$\down' + coords + '$'
        
            ;Find the (data product) (first coordinate) (second coordinate)
            nameParts = stregex(name, '([^x-z]*)([x-z]?)([x-z]?$)', /SUBEXP, /EXTRACT)

            ;Begin the name as the data product
            new_name  = reform(nameParts[1,*])
            
            ;Was the first coordinate found? Replace it.
            iX = where(nameParts[2,*] eq 'x', nX)
            iY = where(nameParts[2,*] eq 'y', nY)
            iZ = where(nameParts[2,*] eq 'z', nZ)
            if nX gt 0 then new_name[iX] += coords[2]
            if nY gt 0 then new_name[iY] += coords[1]
            if nZ gt 0 then new_name[iZ] += coords[0]
        
            ;Was the second coordinate found? Replace it.
            iX = where(nameParts[3,*] eq 'x', nX)
            iY = where(nameParts[3,*] eq 'y', nY)
            iZ = where(nameParts[3,*] eq 'z', nZ)
            if nX gt 0 then new_name[iX] += coords[2]
            if nY gt 0 then new_name[iY] += coords[1]
            if nZ gt 0 then new_name[iZ] += coords[0]
        endcase
        
        'MAGNETOTAIL': begin
            ;Which system?
            coords = mva_frame ? ['L', 'M', 'N'] : ['x', 'y', 'z']
            if _subscript then coords = '$\down' + coords + '$'
        
            ;Find the (data product) (first coordinate) (second coordinate)
            nameParts = stregex(name, '(.*)([x-z]?)([x-z]?$)', /SUBEXP, /EXTRACT)
            
            ;Begin the name as the data product
            new_name  = reform(nameParts[1,*])
            
            ;Was the first coordinate found? Replace it.
            iX = where(nameParts[2,*] eq 'x', nX)
            iY = where(nameParts[2,*] eq 'y', nY)
            iZ = where(nameParts[2,*] eq 'z', nZ)
            if nX gt 0 then new_name[iX] += coords[0]
            if nY gt 0 then new_name[iY] += coords[1]
            if nZ gt 0 then new_name[iZ] += coords[2]
        
            ;Was the second coordinate found? Replace it.
            iX = where(nameParts[3,*] eq 'x', nX)
            iY = where(nameParts[3,*] eq 'y', nY)
            iZ = where(nameParts[3,*] eq 'z', nZ)
            if nX gt 0 then new_name[iX] += coords[0]
            if nY gt 0 then new_name[iY] += coords[1]
            if nZ gt 0 then new_name[iZ] += coords[2]
        endcase
        
        else: message, 'Coordinate system "' + coord_system + '" not recognized.'
    endcase
    
    if _subscript then begin
        ;Convert "e" to a subscript
        nameParts = stregex(new_name, '([^e]*)(e?)([^e]*)', /EXTRACT, /SUBEXP)
        iChange = where(nameParts[2,*] ne '', nChange, COMPLEMENT=iEmpty, NCOMPLEMENT=nEmpty)
        if nChange gt 0 then new_name[iChange] = nameParts[1,iChange] + '$\downe$' + nameParts[3,iChange]
        if nEmpty  gt 0 then new_name[iEmpty]  = nameParts[1,iEmpty]  + nameParts[3,iEmpty]
        
        ;Convert "i" to a subscript
        nameParts = stregex(new_name, '([^i]*)(i?)([^i]*)', /EXTRACT, /SUBEXP)
        iChange = where(nameParts[2,*] ne '', nChange, COMPLEMENT=iEmpty, NCOMPLEMENT=nEmpty)
        if nChange gt 0 then new_name[iChange] = nameParts[1,iChange] + '$\downi$' + nameParts[3,iChange]
        if nEmpty  gt 0 then new_name[iEmpty]  = nameParts[1,iEmpty]  + nameParts[3,iEmpty]
    endif
    
    if n_elements(new_name) eq 1 then new_name = new_name[0]
    return, new_name
end


;-----------------------------------------------------
; Main Level Program: IDL> .r MrSim_Rename \\\\\\\\\\\
;-----------------------------------------------------
