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
;   The purpose of this program is to convert coordinates and species symbols within
;   names of GDA data products to subscripts for plotting purposes. Names are expected to
;   be in the following format::
;       [Quantity][species][optional '-' or '_'][coordinate][optional coordinate]
;
;   where::
;       Species:     'e' or 'i'
;       Coordinates: 'x', 'y', 'z', 'l', 'm', 'n', 'tr', 'para', 'perp'
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Examples:
;   See the main-level program at the end of this file.
;       IDL> .r MrSim_GDA_Format
;
; :Params:
;       NAME:               in, required, type=string
;                           Name of the GDA dataset.
;
; :Keywords:
;       UPCASE:             in, optional, type=boolean, default=0
;                           If set, coordinates will be converted to uppercase.
;       LOWCASE:            in, optional, type=boolean, default=0
;                           If set, coordinates will be converted to lowercase.
;       SYMBOLS:            in, optional, type=boolean, default=0
;                           If set, coordinate names will be converted to symbols::
;                               'PARA' -> '||'
;                               'TR'   -> ''
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
;       2014/01/16  -   Written by Matthew Argall
;       2015/04/03  -   Forgot to ignoring case when searching for coordinates. Fixed. - MRA
;-
function MrSim_GDA_Format, name, $
UPCASE = upcase, $
LOWCASE = lowcase, $
SYMBOLS = symbols
	compile_opt strictarr
	on_error, 2

	;Defaults
	upcase  = keyword_set(upcase)
	lowcase = keyword_set(lowcase)
	symbols = keyword_set(symbols)
	if upcase && lowcase then message, 'UPCASE and LOWCASE are mutually exclusive.'

	;Copy the input. Coordinates and species will be moved from the head to the tail.
	nNames = n_elements(name)
	head   = name
	tail   = strarr(nNames)

	;Possible combinations    
	species = ['e', 'i']
	coords  = ['x', 'y', 'z', 'l', 'm', 'n', 'tr', 'para', 'perp']

	;Step through each name
	for i = 0, nNames - 1 do begin
		thisHead = head[i]
		thisTail = tail[i]
	
	;-----------------------------------------------------
	; Subscript Coordinates \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Check for a coordinate
		regex    = '(' + strjoin(coords, '|') + ')$'
		tf_coord = stregex(thisHead, regex, /FOLD_CASE, /BOOLEAN)
		while tf_coord do begin
			;Location of the coordinate.
			len_head  = strlen(thisHead)
			pos_coord = stregex(thisHead, regex, LENGTH=len_coord, /FOLD_CASE)
			theCoord  = strmid(thisHead, pos_coord)
			
			;Upper or lower case?
			if upcase  then theCoord = strupcase(theCoord)
			if lowcase then theCoord = strlowcase(theCoord)
			
			;Convert to a symbol?
			case strupcase(theCoord) of
				'PARA': theCoord = '||'
				'TR':   theCoord = ''
				else: ;Do nothing
			endcase
		
			;Trim the head and add to the tail, making it a subscript.
			if theCoord ne '' then thisTail = '$\down' + theCoord + '$' + thisTail
			thisHead = strmid(thisHead, 0, pos_coord)

			;Remaining coordinates
			tf_coord = stregex(thisHead, regex, /FOLD_CASE, /BOOLEAN)
		endwhile
		
		;Remove any '[-_]' that appears between coordinates and species
		if stregex(thisHead, '[-_]$', /BOOLEAN) $
			then thisHead = strmid(thisHead, 0, strlen(thisHead)-1)
	
	;-----------------------------------------------------
	; Subscript Species \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Check for species
		regex = '(' + strjoin(species, '|') + ')$'
		tf_species = stregex(thisHead, regex, /FOLD_CASE, /BOOLEAN)
		
		;Make species a subscript
		;   - Electric field Ex must not convert E to a subscript.
		if tf_species && strlen(thisHead) gt 1 then begin
			;Location of the coordinate.
			len_head   = strlen(thisHead)
			pos_coord  = stregex(thisHead, regex, LENGTH=len_coord)
			theSpecies = strmid(thisHead, pos_coord)
		
			;Trim the head and add to the tail, making it a subscript.
			if thisTail eq '' $
				then thisTail = '$\down' + theSpecies + '$'  + thisTail $
				else thisTail = '$\down' + theSpecies + ',$' + thisTail
			thisHead = strmid(thisHead, 0, pos_coord)

			;Remaining coordinates
			tf_coord = stregex(thisHead, regex, /FOLD_CASE, /BOOLEAN)
		endif
		
		head[i] = thisHead
		tail[i] = thisTail
	endfor

	return, head + tail

end


;-----------------------------------------------------
; Main Level Program: IDL> .r MrSim_Rename \\\\\\\\\\\
;-----------------------------------------------------
;EXAMPLE 1 -- XYZ coordinates
names = ['Ay', 'Bx', 'By', 'Bz', 'E-', 'Ex', 'Ey', 'Ez', 'ne', 'ni', $
         'Pe-xx', 'Pe-xy', 'Pe-xz', 'Pe-yx', 'Pe-yy', 'Pe-yz', $
         'Pe-zx', 'Pe-zy', 'Pe-zz', 'Pi-xx', 'Pi-xy', 'Pi-xz', $
         'Pi-yx', 'Pi-yy', 'Pi-yz', 'Pi-zx', 'Pi-zy', 'Pi-zz', $
         'Uex', 'Uey', 'Uez', 'Uix', 'Uiy', 'Uiz']

new_names = MrSim_GDA_Format(names)
print, '-------------------------------------'
print, 'Names:'
print, names, FORMAT='(5(3x, a-5))'
print, 'Names with subscripts:'
print, new_names, FORMAT='(5(3x, a-26))'
print, ''

;EXAMPLE 2 -- LMN coordinates
names = ['Am', 'Bn', 'Bm', 'Bl', 'E-', 'En', 'Em', 'El', 'ne', 'ni', $
         'Pe-nn', 'Pe-nm', 'Pe-nl', 'Pe-mn', 'Pe-mm', 'Pe-ml', $
         'Pe-ln', 'Pe-lm', 'Pe-ll', 'Pi-nn', 'Pi-nm', 'Pi-nl', $
         'Pi-mn', 'Pi-mm', 'Pi-ml', 'Pi-ln', 'Pi-lm', 'Pi-ll', $
         'Uen', 'Uem', 'Uel', 'Uin', 'Uim', 'Uil']

new_names = MrSim_GDA_Format(names)
print, '-------------------------------------'
print, 'Names:'
print, names, FORMAT='(5(3x, a-5))'
print, 'Names with subscripts:'
print, new_names, FORMAT='(5(3x, a-26))'
print, ''


end