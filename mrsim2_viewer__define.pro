; docformat = 'rst'
;
;+
;   The purpose of this program is to provide a base class for 2D and 3D simulation
;   data provided by Bill Daughton from Los Alamos National Laboratory.
;
;   MrSim2 is meant to be subclassed by a class that over-rides the ReadData method.
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
;       Copyright 2013 by the University of New Hampshire
;
; :History:
;   Modification History::
;       2014/09/12  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   This method initializes the MrSim2_Viewer class.
;
; :Params:
;       THESIM:         in, optional, type=string/integer
;                       Name or index of the simulation to be created. If not present,
;                           a list of all available simulations will be printed to the
;                           display window.
;
; :Keywords:
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrSim_Create is also accepted via
;                           keyword inheritance.
;-
function MrSim2_Viewer::BuildGUI, theSim, $
GROUP_LEADER=group_leader
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, 0
    endif
    
;---------------------------------------------------------------------
;Make the Top Level Base /////////////////////////////////////////////
;---------------------------------------------------------------------

	;Make a top-level base with or without a groupleader. cdf_read_gui2 is called by other
	;blocking widgets, so if a group_leader is given, then make cdf_read_gui2 modal.
	if n_elements(group_leader) ne 0 then begin
	    no_block = 0
	    tlb = widget_base(GROUP_LEADER=group_leader, TITLE='MrSim Viewer', /COLUMN, $
                          XOFFSET=100, YOFFSET=0, UNAME='tlb', /BASE_ALIGN_CENTER, $
                          /MODAL)
	endif else begin
	    no_block = 1
	    tlb = widget_base(TITLE='MrSim Viewer', /COLUMN, XOFFSET=100, YOFFSET=0, $
                          UNAME='tlb', /BASE_ALIGN_CENTER)
	endelse
    
    self.tlb = tlb

    return, 1
end


;+
;   Clean up after the object is destroyed.
;-
pro MrSim2_Viewer::CLEANUP
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Free all of the pointers
    if obj_valid(self.oSim) then obj_destroy, oSim
end


;+
;   This method initializes the MrSim2_Viewer class.
;
; :Params:
;       THESIM:         in, optional, type=string/integer
;                       Name or index of the simulation to be created. If not present,
;                           a list of all available simulations will be printed to the
;                           display window.
;
; :Keywords:
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrSim_Create is also accepted via
;                           keyword inheritance.
;-
function MrSim2_Viewer::INIT, theSim, $
_REF_EXTRA=extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, 0
    endif

;-------------------------------------------------------
; Create the Simulation Object /////////////////////////
;-------------------------------------------------------
    tname = size(theSim, /TNAME)
    
    ;Simulation object?
    if tname eq 'OBJREF' then begin
        ;Make sure it is a MrSim2 object.
        if obj_class(oSim, 'MrSim2') eq 0 then $
            message, 'SIM_OBJECT must be a subclass of MrSim2.'
        
        ;Store the object as a property
        self.oSim = 
    
    ;Simulation name/number?
    endif else begin
        sim_object = MrSim_Create(theSim, _STRICT_EXTRA=extra)
        if obj_valid(sim_object) eq 0 $
            then return, 0 $
            else self.oSim = sim_object
    endelse

;-------------------------------------------------------
; Create the GUI ///////////////////////////////////////
;-------------------------------------------------------
    self -> BuildGUI

    return, 1
end


;+
;   This definition statement for the MrSim2 class.
;
; :Params:
;       CLASS           out, optional, type=named structure
;                       The class definition structure.
;
; :FIELDS:
;-
pro MrSim2_Viewer__DEFINE, class
    compile_opt strictarr
    
    class = { MrSim2_Viewer, $
              inherits IDL_Object, $
             
              ;Object Properties
              oSim: obj_new() $
            }
end