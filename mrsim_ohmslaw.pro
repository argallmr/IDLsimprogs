; docformat = 'rst'
;
; NAME:
;    MrSim_ColorSlab
;
;*****************************************************************************************
;   Copyright (c) 2013, Matthew Argall                                                   ;
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
;   The purpose of this program is to create a color image of a single data product with
;   the option of overlaying contours and vertical or horizontal lines.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Uses:
;   Uses the following external programs::
;       cgErrorMSG.pro
;       cgConLevels.pro
;       GetMrWindows.pro
;       MrImage.pro
;       MrPlot.pro
;       MrColorbar.pro
;       MrContour.pro
;       MrSim2D__Define.pro
;       MrSim3D__Define.pro
;       MrWindow.pro
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
;       2014/03/28  -   Written by Matthew Argall
;       2014/04/08  -   Include separate plots for each term in Ohm's law. Added the IONS
;                           and ELECTRONS keywords. - MRA
;       2014/04/13  -   Forgot to divide the Hall terms by the sum of the densities. Fixed. - MRA
;       2014/04/16  -   Added the TOTAL, VXB, JXB, GRADPE keywords. Created individual
;                           function for each term in Ohm's Law. - MRA
;-
;*****************************************************************************************
;+
;   Plot the components of the pressure gradient term in the Generalized Ohm's Law
;
; :Params:
;       COMPONENT:          in, required, type=string
;                           The component of the Hall electric field to be plotted.
;       OSIM:               in, required, type=objref
;                           A "MRSIM2D" or "MRSIM3D" object containing information about
;                               the simulation domain.
; :Keywords:
;       ELECTRONS:          in, optional, type=boolean, default=0
;                           If set, the electron current will be plotted instead of the
;                               total current.
;       IONS:               in, optional, type=boolean, default=0
;                           If set, the ion current will be plotted instead of the
;                               total current.
;-
pro MrSim_OhmsLaw_Total, oSim, component, cut, $
ELECTRONS = electrons, $
HORIZONTAL = horizontal, $
IONS = ions
    compile_opt strictarr
    on_error, 2

;-------------------------------------------------------
; Defaults ////////////////////////////////////////////
;-------------------------------------------------------
    _comp     = strupcase(component)
    electrons = keyword_set(electrons)
    ions      = keyword_set(ions)
    Sim3D     = keyword_set(Sim3D)
    if electrons + ions eq 0 then fluid = 1 else fluid = 0
    if electrons + ions gt 1 then message, 'Keywords ELECTRONS and IONS are mutually exclusive.'
    
;-------------------------------------------------------
; Prepare to Plot //////////////////////////////////////
;-------------------------------------------------------
    ;Get the simulation size and time
    oSim -> GetProperty, MVA_FRAME=mva_frame, AXIS_LABELS=axlabels, TIME=time, XSIM=XSim, ZSIM=ZSim
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units

    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't=' + string(time*dtxwci, FORMAT='(f0.1)') + ' $\Omega$$\downc,i$$\up-1$' $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')
        
    ;Title
    if obj_class(oSim) eq 'MRSIM3D' then begin
        oSim -> GetProperty, YRANGE=yrange
        title += string(FORMAT='(%"  [%s,%s]=(%0.1f,%0.1f)%s")', axlabels[0:1], cut, yrange[0], units)
    endif else begin
        title += string(FORMAT='(%"  %s=%0.1f%s")', axlabels[0], cut, units)
    endelse
    
    ;Species to be plotted
    case 1 of
        electrons: begin
            n_name = 'ne'
            V_name = 'Ue'
            P_name = 'Pe'
        endcase
        
        ions: begin
            n_name = 'ni'
            V_name = 'Ui'
            P_name = 'Pi'
        endcase
        
        fluid: begin
            n_name = 'ne'
            V_name = 'V'
            P_name = 'Pe'
        endcase
    endcase

;-------------------------------------------------------
; Plot Ohm's Law ///////////////////////////////////////
;-------------------------------------------------------
    case _comp of
        'X': begin
            E     = MrSim_LineCut(oSim, 'Ex',      cut, /CURRENT, HORIZONTAL=horizontal)
            VxB   = MrSim_LineCut(oSim, 'VxB_x',   cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
            JxB   = MrSim_LineCut(oSim, 'JxB_x',   cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Forest Green')
            divPe = MrSim_LineCut(oSim, 'divPe_x', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Red')
        endcase
        'Y': begin
            E     = MrSim_LineCut(oSim, 'Ey',      cut, /CURRENT, HORIZONTAL=horizontal)
            VxB   = MrSim_LineCut(oSim, 'VxB_y',   cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
            JxB   = MrSim_LineCut(oSim, 'JxB_y',   cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Forest Green')
            divPe = MrSim_LineCut(oSim, 'divPe_y', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Red')
        endcase
        'Z': begin
            E     = MrSim_LineCut(oSim, 'Ez',      cut, /CURRENT, HORIZONTAL=horizontal)
            VxB   = MrSim_LineCut(oSim, 'VxB_z',   cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
            JxB   = MrSim_LineCut(oSim, 'JxB_z',   cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Forest Green')
            divPe = MrSim_LineCut(oSim, 'divPe_z', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Red')
        endcase
        else: message, 'Component "' + component + '" not rectognized. Choose from {"X" | "Y" | "Z"}'
    endcase

;-------------------------------------------------------
; Adjust the Plot //////////////////////////////////////
;-------------------------------------------------------
    
    ;MVA labels?
    if mva_frame then begin
        case _comp of
            'X': _comp = axlabels[0]
            'Y': _comp = axlabels[1]
            'Z': _comp = axlabels[2]
        endcase
    endif
    
    ;Rename
    E.name     = 'Total E'       + component
    VxB.name   = 'Total (VxB)'   + component
    JxB.name   = 'Total (JxB)'   + component
    divPe.name = 'Total div(Pe)' + component
    
    ;Set the range
    E     -> GetData, Edata
    VxB   -> GetData, VxBdata
    JxB   -> GetData, JxBdata
    divPe -> GetData, divPedata
    
    ;Find the min and max
    maxRange = max([max(temporary(Edata),   MIN=Emin),   max(-temporary(VxBdata),  MIN=VxBmin), $
                    max(temporary(JxBdata), MIN=JxBmin), max(temporary(divPedata), MIN=divPemin)])
    range = [min([Emin, VxBmin, JxBmin, divPemin]), maxRange]
    E.YRANGE = range

    ;Plot -VxB
    VxB -> GetData, x, y
    VxB -> SetData, temporary(x), -temporary(y)
    
    legend_titles = ['E', '-(VxB)', '(JxB)', '(divPe)'] + '$\down' + _comp +'$'
    legend_titles[2] += '/en'
    
    ;Create a legend
    ohmLegend = MrLegend(TARGET=E, LOCATION=4, LENGTH=0, NAME="Ohm's Law", $
                         TCOLORS=['Black', 'Blue', 'Forest Green', 'Red'], $
                         TITLE=legend_titles)

end


;+
;   Plot the components of the pressure gradient term in the Generalized Ohm's Law
;
; :Params:
;       COMPONENT:          in, required, type=string
;                           The component of the Hall electric field to be plotted.
;       OSIM:               in, required, type=objref
;                           A "MRSIM2D" or "MRSIM3D" object containing information about
;                               the simulation domain.
; :Keywords:
;       ELECTRONS:          in, optional, type=boolean, default=0
;                           If set, the electron current will be plotted instead of the
;                               total current.
;       IONS:               in, optional, type=boolean, default=0
;                           If set, the ion current will be plotted instead of the
;                               total current.
;-
pro MrSim_OhmsLaw_GradPe, oSim, component, cut, $
ELECTRONS = electrons, $
HORIZONTAL = horizontal, $
IONS = ions
    compile_opt strictarr
    on_error, 2

;-------------------------------------------------------
; Defaults ////////////////////////////////////////////
;-------------------------------------------------------
    electrons = keyword_set(electrons)
    ions      = keyword_set(ions)
    if electrons + ions eq 0 then fluid = 1 else fluid = 0
    if electrons + ions gt 1 then message, 'Keywords ELECTRONS and IONS are mutually exclusive.'
    
;-------------------------------------------------------
; Prepare to Plot //////////////////////////////////////
;-------------------------------------------------------
    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, XSIM=XSim, ZSIM=ZSim
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    
    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't=' + string(time*dtxwci, FORMAT='(f0.1)') + ' $\Omega$$\downc,i$$\up-1$' $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')
        
    ;Title
    if obj_class(oSim) eq 'MRSIM3D' then begin
        oSim -> GetProperty, YSLICE=yslice
        title += '  [x,y]=(' + string(cut, FORMAT='(f0.1)')    + ', ' $
                             + string(yslice, FORMAT='(f0.1)') + ')' + units
    endif else begin
        title += '  x=' + string(cut, FORMAT='(f0.1)') + units
    endelse
    
    ;Species to be plotted
    case 1 of
        electrons: begin
            n_name = 'ne'
            V_name = 'Ue'
            P_name = 'Pe'
        endcase
        
        ions: begin
            n_name = 'ni'
            V_name = 'Ui'
            P_name = 'Pi'
        endcase
        
        fluid: begin
            n_name = 'ne'
            V_name = 'V'
            P_name = 'Pe'
        endcase
    endcase

;-------------------------------------------------------
; Plot the Inertial Term ///////////////////////////////
;-------------------------------------------------------
    subX = '$\downX$'
    subY = '$\downY$'
    subZ = '$\downZ$'
    case strupcase(component) of
        'X': begin
            ;Ex and Divergence of Pressure (x-component)
            E     = MrSim_LineCut(oSim, 'Ex', cut, /CURRENT, HORIZONTAL=horizontal)
            divPe = MrSim_LineCut(oSim, 'div' + P_name + '_x', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
                                
            ;Get Data
            n_e     = oSim -> LineCuts(n_name,         cut, pos,      HORIZONTAL=horizontal)
            divPe_1 = oSim -> LineCuts(P_name + '-xx', cut, pos, /DX, HORIZONTAL=horizontal)
            divPe_2 = oSim -> LineCuts(P_name + '-xy', cut,      /DX, HORIZONTAL=horizontal)
            divPe_3 = oSim -> LineCuts(P_name + '-xz', cut,      /DX, HORIZONTAL=horizontal)
                        
            ;Titles for the legend
            titles = '(div ' + P_name + ')' + [subX, subX+subX, subX+subY, subX+subZ]
        endcase
        'Y': begin
            ;Ey and Divergence of Pressure (y-component)
            E       = MrSim_LineCut(oSim, 'Ey', cut, /CURRENT, HORIZONTAL=horizontal)
            divPe   = MrSim_LineCut(oSim, 'div' + P_name + '_y', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
            
            ;Get data
            n_e     = oSim -> LineCuts(n_name,         cut, pos,      HORIZONTAL=horizontal)
            divPe_1 = oSim -> LineCuts(P_name + '-xy', cut, pos, /DY, HORIZONTAL=horizontal)
            divPe_2 = oSim -> LineCuts(P_name + '-yy', cut,      /DY, HORIZONTAL=horizontal)
            divPe_3 = oSim -> LineCuts(P_name + '-yz', cut,      /DY, HORIZONTAL=horizontal)
            
            ;Title for legend
            titles = '(div' + Pe_name + ')' + [subY, subY+subX, subY+subY, subY+subZ]
        endcase
        'Z': begin
            ;Ez and Divergence of Pressure (z-component)
            E       = MrSim_LineCut(oSim, 'Ez', cut, /CURRENT, HORIZONTAL=horizontal)
            divPe   = MrSim_LineCut(oSim, 'div' + P_name + '_z', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
            
            ;Get data
            n_e     = oSim -> LineCuts(n_name,         cut, pos,      HORIZONTAL=horizontal)
            divPe_1 = oSim -> LineCuts(P_name + '-xz', cut, pos, /DZ, HORIZONTAL=horizontal)
            divPe_2 = oSim -> LineCuts(P_name + '-yz', cut,      /DZ, HORIZONTAL=horizontal)
            divPe_3 = oSim -> LineCuts(P_name + '-zz', cut,      /DZ, HORIZONTAL=horizontal)
            
            ;Title for legend
            titles = '(div ' + P_name + ')' + [subZ, subX+subZ, subY+subZ, subZ+subZ]
        endcase
        else: message, 'Component "' + component + '" not rectognized. Choose from {"X" | "Y" | "Z"}'
    endcase

;-------------------------------------------------------
; Adjust the Plot //////////////////////////////////////
;-------------------------------------------------------
    
    ;Rename
    E.name       = 'E' + component + ' vs. E inert'
    divPe.name   = 'Inertial div(Pe)' + component
    
    ;Find the min and max
    E     -> GetData, Edata
    divPe -> GetData, Pdata

    maxRange = max([max(temporary(Edata), MIN=Emin), max(temporary(Pdata), MIN=Pmin), $
                    max(divPe_1, MIN=P1min), max(divPe_2, MIN=P2min), $
                    max(divPe_3, MIN=P3min)])
    range = [min([Emin, Pmin, P1min, P2min, P3min]), maxRange]
    E.YRANGE=range

    ;Overplot the individual terms
    !Null = MrPlot(pos, divPe_1/n_e, OVERPLOT=E, NAME='Inertial div(Pe)' + component + 'X', COLOR='Forest Green')
    !NULL = MrPlot(pos, divPe_2/n_e, OVERPLOT=E, NAME='Inertial div(Pe)' + component + 'Y', COLOR='Red')
    !NULL = MrPlot(pos, divPe_3/n_e, OVERPLOT=E, NAME='Inertial div(Pe)' + component + 'Z', COLOR='Purple')
    
    ;Create a legend
    ohmLegend = MrLegend(TARGET=E, TITLE=['E'+subZ, titles], $
                         TCOLORS=['Black', 'Blue', 'Forest Green', 'Red', 'Purple'], $
                         LOCATION=4, LENGTH=0, NAME="Ohm's Law: div(Pe) term")
end


;+
;   Plot the components of the convective term in the Generalized Ohm's Law
;
; :Params:
;       COMPONENT:          in, required, type=string
;                           The component of the Hall electric field to be plotted.
;       OSIM:               in, required, type=objref
;                           A "MRSIM2D" or "MRSIM3D" object containing information about
;                               the simulation domain.
; :Keywords:
;       ELECTRONS:          in, optional, type=boolean, default=0
;                           If set, the electron current will be plotted instead of the
;                               total current.
;       IONS:               in, optional, type=boolean, default=0
;                           If set, the ion current will be plotted instead of the
;                               total current.
;-
pro MrSim_OhmsLaw_VxB, oSim, component, cut, $
ELECTRONS = electrons, $
HORIZONTAL = horizontal, $
IONS = ions
    compile_opt strictarr
    on_error, 2

;-------------------------------------------------------
; Defaults ////////////////////////////////////////////
;-------------------------------------------------------
    electrons = keyword_set(electrons)
    ions      = keyword_set(ions)
    if electrons + ions eq 0 then fluid = 1 else fluid = 0
    if electrons + ions gt 1 then message, 'Keywords ELECTRONS and IONS are mutually exclusive.'
    
;-------------------------------------------------------
; Prepare to Plot //////////////////////////////////////
;-------------------------------------------------------
    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, XSIM=XSim, ZSIM=ZSim
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    
    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't=' + string(time*dtxwci, FORMAT='(f0.1)') + ' $\Omega$$\downc,i$$\up-1$' $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')
        
    ;Title
    if obj_class(oSim) eq 'MRSIM3D' then begin
        oSim -> GetProperty, YSLICE=yslice
        title += '  (x,y)=(' + string(cut, FORMAT='(f0.1)')    + ', ' $
                             + string(yslice, FORMAT='(f0.1)') + ')' + units
    endif else begin
        title += '  x=' + string(cut, FORMAT='(f0.1)') + units
    endelse
    
    ;Species to be plotted
    case 1 of
        electrons: begin
            n_name = 'ne'
            V_name = 'Ue'
            P_name = 'Pe'
        endcase
        
        ions: begin
            n_name = 'ni'
            V_name = 'Ui'
            P_name = 'Pi'
        endcase
        
        fluid: begin
            n_name = 'ne'
            V_name = 'V'
            P_name = 'Pe'
        endcase
    endcase

;-------------------------------------------------------
; Plot the Convective Term /////////////////////////////
;-------------------------------------------------------
    subX = '$\downX$'
    subY = '$\downY$'
    subZ = '$\downZ$'
    case strupcase(component) of
        'X': begin
            ;Ex and Convective Electric Field (x-component)
            E   = MrSim_LineCut(oSim, 'Ex', cut, /CURRENT, HORIZONTAL=horizontal)
            VxB = MrSim_LineCut(oSim, V_name + 'xB_x', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
                                
            ;Get Data
            Vy  = oSim -> LineCuts(V_name + 'y', cut, pos, HORIZONTAL=horizontal)
            Vz  = oSim -> LineCuts(V_name + 'z', cut,      HORIZONTAL=horizontal)
            By  = oSim -> LineCuts('By',         cut,      HORIZONTAL=horizontal)
            Bz  = oSim -> LineCuts('Bz',         cut,      HORIZONTAL=horizontal)
            
            ;-(VxB)_x = -(Vy*Bz) + (Vz*By)
            vB1 = -Vy * Bz
            vB2 =  Vz * By
            
            ;Titles for the legend
            titles = V_name + ['xB', subY+'B'+subZ, subZ+'B'+subY]
            titles[2] = '-' + titles[2]
        endcase
        'Y': begin
            ;Ex and Convective Electric Field (x-component)
            E   = MrSim_LineCut(oSim, 'Ey', cut, /CURRENT, HORIZONTAL=horizontal)
            VxB = MrSim_LineCut(oSim, V_name + 'xB_y', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
                                
            ;Get Data
            Vx  = oSim -> LineCuts(V_name + 'x', cut, pos, HORIZONTAL=horizontal)
            Vz  = oSim -> LineCuts(V_name + 'z', cut,      HORIZONTAL=horizontal)
            Bx  = oSim -> LineCuts('Bx',         cut,      HORIZONTAL=horizontal)
            Bz  = oSim -> LineCuts('Bz',         cut,      HORIZONTAL=horizontal)
            
            ;-(VxB)_y = -(Vx*Bz) + (Vz*Bx)
            vB1 = -Vz * Bx
            vB2 =  Vx * Bz
            
            ;Titles for the legend
            titles = V_name + ['xB', subX+'B'+subZ, subZ+'B'+subX]
            titles[2] = '-' + titles[2]
        endcase
        'Z': begin
            ;Ex and Convective Electric Field (x-component)
            E   = MrSim_LineCut(oSim, 'Ez', cut, /CURRENT, HORIZONTAL=horizontal)
            VxB = MrSim_LineCut(oSim, V_name + 'xB_z', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
                                
            ;Get Data
            Vx  = oSim -> LineCuts(V_name + 'x', cut, pos, HORIZONTAL=horizontal)
            Vy  = oSim -> LineCuts(V_name + 'y', cut,      HORIZONTAL=horizontal)
            Bx  = oSim -> LineCuts('Bx',         cut,      HORIZONTAL=horizontal)
            By  = oSim -> LineCuts('By',         cut,      HORIZONTAL=horizontal)
            
            ;-(VxB)_z = -(Vx*By) + (Vy*Bx)
            vB1 = -Vx * By
            vB2 =  Vy * Bx
            
            ;Titles for the legend
            titles = V_name + ['xB)'+subZ, subX+'B'+subY, subY+'B'+subX]
            titles[0]   = '(' + titles[0]
            titles[0:1] = '-' + titles[0:1]
        endcase
        else: message, 'Component "' + component + '" not rectognized. Choose from {"X" | "Y" | "Z"}'
    endcase

;-------------------------------------------------------
; Adjust the Plot //////////////////////////////////////
;-------------------------------------------------------
    ;Rename
    E.name     = 'E' + component + ' vs. Ec'
    VxB.name   = 'Convective (VxB)'   + component
    
    ;Find the min and max
    E   -> GetData, Edata
    VxB -> GetData, vData
    maxRange = max([max(temporary(Edata), MIN=Emin), max(-vDATA, MIN=Vmin), $
                    max(vB1, MIN=vB1min), max(vB2, MIN=vB2min)])
    range = [min([Emin, Vmin, vB1min, vB2min]), maxRange]
    E.YRANGE=range
    
    ;Plot -VxB
    VxB -> SetData, -vData
    
    ;Overplot the individual terms
    !Null = MrPlot(pos, vB1, OVERPLOT=E, NAME='Convective -Vx*By', COLOR='Forest Green')
    !NULL = MrPlot(pos, vB2, OVERPLOT=E, NAME='Convective Vy*Bx', COLOR='Red')
    
    ;Create a legend
    ohmLegend = MrLegend(TARGET=E, TITLE=['E'+subZ, titles], $
                         TCOLORS=['Black', 'Blue', 'Forest Green', 'Red'], $
                         LOCATION=4, LENGTH=0, NAME="Ohm's Law: VxB term")
end


;+
;   Plot the components of the Hall term in the Generalized Ohm's Law
;
; :Params:
;       COMPONENT:          in, required, type=string
;                           The component of the Hall electric field to be plotted.
;       OSIM:               in, required, type=objref
;                           A "MRSIM2D" or "MRSIM3D" object containing information about
;                               the simulation domain.
; :Keywords:
;       ELECTRONS:          in, optional, type=boolean, default=0
;                           If set, the electron current will be plotted instead of the
;                               total current.
;       IONS:               in, optional, type=boolean, default=0
;                           If set, the ion current will be plotted instead of the
;                               total current.
;-
pro MrSim_OhmsLaw_JxB, oSim, component, cut, $
ELECTRONS = electrons, $
HORIZONTAL = horizontal, $
IONS = ions
    compile_opt strictarr
    on_error, 2

;-------------------------------------------------------
; Defaults ////////////////////////////////////////////
;-------------------------------------------------------
    electrons = keyword_set(electrons)
    ions      = keyword_set(ions)
    Sim3D     = keyword_set(Sim3D)
    if electrons + ions eq 0 then fluid = 1 else fluid = 0
    if electrons + ions gt 1 then message, 'Keywords ELECTRONS and IONS are mutually exclusive.'
    
;-------------------------------------------------------
; Prepare to Plot //////////////////////////////////////
;-------------------------------------------------------

    ;Get the simulation size and time
    oSim -> GetProperty, TIME=time, XSIM=XSim, ZSIM=ZSim
    oSim -> GetInfo, DTXWCI=dtxwci, UNITS=units
    
    ;Time is inverse gyro-time?
    if n_elements(dtxwci) gt 0 $
        then title = 't=' + string(time*dtxwci, FORMAT='(f0.1)') + ' $\Omega$$\downc,i$$\up-1$' $
        else title = 't$\downindex$=' + string(time, FORMAT='(i0)')
        
    ;Title
    if obj_class(oSim) eq 'MRSIM3D' then begin
        oSim -> GetProperty, YSLICE=yslice
        title += '  [x,y]=(' + string(cut, FORMAT='(f0.1)')    + ', ' $
                             + string(yslice, FORMAT='(f0.1)') + ')' + units
    endif else begin
        title += '  x=' + string(cut, FORMAT='(f0.1)') + units
    endelse
    
    ;Species to be plotted
    case 1 of
        electrons: begin
            n_name = 'ne'
            V_name = 'Ue'
            P_name = 'Pe'
        endcase
        
        ions: begin
            n_name = 'ni'
            V_name = 'Ui'
            P_name = 'Pi'
        endcase
        
        fluid: begin
            n_name = 'ne'
            V_name = 'V'
            P_name = 'Pe'
        endcase
    endcase

;-------------------------------------------------------
; Plot the Hall Term ///////////////////////////////////
;-------------------------------------------------------
    subX = '$\downX$'
    subY = '$\downY$'
    subZ = '$\downZ$'
    case strupcase(component) of
        'X': begin
            ;Ex and Hall Electric Field (x-component)
            E   = MrSim_LineCut(oSim, 'Ex', cut, /CURRENT, HORIZONTAL=horizontal)
            JxB = MrSim_LineCut(oSim, 'JxB_x', cut, OVERPLOT=E, COLOR='Blue', HORIZONTAL=horizontal)
                                
            ;Get Data
            n_i  = oSim -> LineCuts('ni', cut,      HORIZONTAL=horizontal)
            n_e  = oSim -> LineCuts('ne', cut,      HORIZONTAL=horizontal)
            Jy   = oSim -> LineCuts('Jy', cut, pos, HORIZONTAL=horizontal)
            Jz   = oSim -> LineCuts('Jz', cut,      HORIZONTAL=horizontal)
            By   = oSim -> LineCuts('By', cut,      HORIZONTAL=horizontal)
            Bz   = oSim -> LineCuts('Bz', cut,      HORIZONTAL=horizontal)
            
            ;E_N = (Jy*Bz) - (Jz*By) = JB1 + JB2
            JB1 = ( Jy * Bz) / (n_i + n_e)
            JB2 = (-Jz * By) / (n_i + n_e)
            
            ;Titles for the legend
            titles = ['(JxB)'+subX, 'J'+subY+'B'+subZ, 'J'+subZ+'B'+subY] + '/en'
        endcase
        'Y': begin
            ;Ey and Hall Electric Field (y-component)
            E   = MrSim_LineCut(oSim, 'Ey', cut, /CURRENT, HORIZONTAL=horizontal)
            JxB = MrSim_LineCut(oSim, 'JxB_y', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
                                
            ;Get Data
            n_i = oSim -> LineCuts('ni', cut,      HORIZONTAL=horizontal)
            n_e = oSim -> LineCuts('ne', cut,      HORIZONTAL=horizontal)
            Jx  = oSim -> LineCuts('Jx', cut, pos, HORIZONTAL=horizontal)
            Jz  = oSim -> LineCuts('Jz', cut,      HORIZONTAL=horizontal)
            Bx  = oSim -> LineCuts('Bx', cut,      HORIZONTAL=horizontal)
            Bz  = oSim -> LineCuts('Bz', cut,      HORIZONTAL=horizontal)
            
            ;(JxB)_y = (Jx*Bz) - (Jz*Bx) = JB1 + JB2
            JB1 = ( Jz * Bx) / (n_i + n_e)
            JB2 = (-Jx * Bz) / (n_i + n_e)
            
            ;Titles for the legend
            titles = ['(JxB)'+subY, 'J'+subZ+'B'+subX, '-J'+subX+'B'+subZ] + '/en'
        endcase
        'Z': begin
            ;Ez and Hall Electric Field (z-component)
            E   = MrSim_LineCut(oSim, 'Ez', cut, /CURRENT, HORIZONTAL=horizontal)
            JxB = MrSim_LineCut(oSim, 'JxB_z', cut, OVERPLOT=E, HORIZONTAL=horizontal, COLOR='Blue')
                                
            ;Get Data
            Jx  = oSim -> LineCuts('Jx', cut, pos, HORIZONTAL=horizontal)
            n_i = oSim -> LineCuts('ni', cut,      HORIZONTAL=horizontal)
            n_e = oSim -> LineCuts('ne', cut,      HORIZONTAL=horizontal)
            Jy  = oSim -> LineCuts('Jy', cut,      HORIZONTAL=horizontal)
            Bx  = oSim -> LineCuts('Bx', cut,      HORIZONTAL=horizontal)
            By  = oSim -> LineCuts('By', cut,      HORIZONTAL=horizontal)
            
            ;(JxB)_z = (Jx*By) - (Jy*Bx) = JB1 + JB2
            JB1 = ( Jx * By) / (n_i + n_e)
            JB2 = (-Jy * Bx) / (n_i + n_e)
            
            ;Titles for the legend
            titles = ['(JxB)'+subZ, 'J'+subX+'B'+subY, '-J'+subY+'B'+subX] + '/en'
        endcase
        else: message, 'Component "' + component + '" not rectognized. Choose from {"X" | "Y" | "Z"}'
    endcase
    
    ;Rename
    E.name     = 'E' + component + ' vs. Hall E'
    JxB.name   = 'Hall (JxB)'   + component
    
    ;Find the min and max
    E   -> GetData, Edata
    JxB -> GetData, JData
    maxRange = max([max(temporary(Edata), MIN=Emin), max(Jdata, MIN=Jmin), $
                    max(JB1, MIN=JB1min), max(JB2, MIN=JB2min)])
    range = [min([Emin, Jmin, JB1min, JB2min]), maxRange]
    E.YRANGE=range
    
    ;Overplot the individual terms
    !Null = MrPlot(pos, JB1, OVERPLOT=E, NAME='Hall Jx*By', COLOR='Forest Green')
    !NULL = MrPlot(pos, JB2, OVERPLOT=E, NAME='Hall -Jy*Bx', COLOR='Red')
    
    ;Create a legend
    ohmLegend = MrLegend(TARGET=E, TITLE=['E'+subZ, titles], $
                         TCOLORS=['Black', 'Blue', 'Forest Green', 'Red'], $
                         LOCATION=4, LENGTH=0, NAME="Ohm's Law: JxB term")
end


;+
;   Plot the components of Ohm's Law
;
; :Params:
;       THESIM:             in, required, type=string/integer/object
;                           The name or number of the simulation to be used, or a
;                               valid simulation object. See MrSim_Create.pro.
;       COMPONENT:          in, required, type=string, default='X'
;                           The component of the Generalized Ohm's Law to be plotted.
;                               Options are "X", "Y", and "Z".
;       CUT:                in, required, type=float
;                           Location along the horizontal axis at which to take a vertical
;                               cut through the simulation domain.
;       TIME:               in, optional, type=int
;                           The simulation time for which a velocity vector field
;                               is to be plotted. Only used if `SIM_OBJECT` is not given.
; :Keywords:
;       CURRENT:            in, optional, type=boolean, default=0
;                           Added the plots to the current MrWindow graphics window.
;       ELECTRONS:          in, optional, type=boolean, default=0
;                           If set, the electron current will be plotted instead of the
;                               total current.
;       GRADPE:             in, optional, type=boolean, default=0
;                           If set, the pressure tensor term will be plotted.
;       HORIZONTAL:         in, optional, type=boolean, default=0
;                           If set, a horizontal cut will be taken. The default is
;                               to take a vertical cut. For an "XY" orientation, "X" is
;                               horizontal, "Y" is vertical. Similar for "XZ", etc.
;       IONS:               in, optional, type=boolean, default=0
;                           If set, the ion current will be plotted instead of the
;                               total current.
;       JXB:                in, optional, type=boolean, default=0
;                           If set, the Hall term will be plotted.
;       OFILENAME:          in, optional, type=string, default=''
;                           Name of a file to which the graphics will be saved.
;       SIM_OBJECT:         in, optional, type=objref
;                           A "MrSim2D" or a "MrSim3D" object.
;       TOTAL:              in, optional, type=boolean, default=0
;                           If set, the total E, VxB, JxB, and gradP will be plotted.
;       VXB:                in, optional, type=boolean, default=0
;                           If set, the convective term will be plotted.
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by MrSim2D__Define or MrSim3D__Define,
;                               depending on the nature of the simulation.
;
; :Returns:
;       OHMWIN:             MrWindow graphic window containing the requested graphics.
;                               If the image data does not exist, an invalid object
;                               will be returned.
;-
function MrSim_OhmsLaw, theSim, component, cut, time, $
CURRENT = current, $
ELECTRONS = electrons, $
GRADPE = gradPe, $
HORIZONTAL = horizontal, $
IONS = ions, $
JXB = JxB, $
OFILENAME = ofilename, $
TOTAL = ETotal, $
VXB = VxB, $
_REF_EXTRA = extra
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        if osim_created        && arg_present(oSim)    eq 0 then obj_destroy, oSim
        if obj_valid(colorWin) && keyword_set(current) eq 0 then obj_destroy, colorWin
        void = cgErrorMSG()
        return, obj_new
    endif

;-------------------------------------------------------
; Check Simulation /////////////////////////////////////
;-------------------------------------------------------
    osim_created = 0B
    
    ;Simulation name or number?
    if MrIsA(theSim, 'STRING') || MrIsA(theSim, 'INTEGER') then begin
        oSim = MrSim_Create(theSim, time, yslice, _STRICT_EXTRA=extra)
        if obj_valid(oSim) eq 0 then return, obj_new()
        osim_created = 1B
        
    ;Object?
    endif else if MrIsA(theSim, 'OBJREF') then begin
        if obj_isa(theSim, 'MRSIM') eq 0 $
            then message, 'THESIM must be a subclass of the MrSim class.' $
            else oSim = theSim
            
    ;Somthing else
    endif else begin
        MrSim_Which
        message, 'THESIM must be a simulation name, number, or object.'
    endelse
    sim_class = obj_class(oSim)

;-------------------------------------------------------
; Defaults ////////////////////////////////////////////
;-------------------------------------------------------
    current   = keyword_set(current)
    electrons = keyword_set(electrons)
    ETotal    = keyword_set(ETotal)
    gradPe    = keyword_set(gradPe)
    ions      = keyword_set(ions)
    JxB       = keyword_set(JxB)
    Sim3D     = keyword_set(Sim3D)
    VxB       = keyword_set(VxB)
    if n_elements(ofilename) eq 0 then ofilename = ''
    
    ;Dependencies
    if electrons + ions eq 0 then fluid = 1 else fluid = 0
    if electrons + ions gt 1 then message, 'Keywords ELECTRONS and IONS are mutually exclusive.'
    
    if (ETotal + gradPe + JxB + VxB) eq 0 then begin
        ETotal = 1
        gradPe = 1
        JxB    = 1
        VxB    = 1
    endif
    
    ;Buffer the output?
    if current eq 0 then $
        if ofilename eq '' then buffer = 0 else buffer = 1
    
;-------------------------------------------------------
; Prepare to Plot //////////////////////////////////////
;-------------------------------------------------------
    ;Create a window
    if current $
        then ohmWin = GetMrWindows(/CURRENT) $
        else ohmWin = MrWindow(REFRESH=0, YGAP=0.5, YSIZE=690, BUFFER=buffer)

;-------------------------------------------------------
; Plot Ohm's Law ///////////////////////////////////////
;-------------------------------------------------------
    if ETotal then MrSim_OhmsLaw_Total,  oSim, component, cut, ELECTRONS=electrons, HORIZONTAL=horizontal, IONS=ions
    if VxB    then MrSim_OhmsLaw_VxB,    oSim, component, cut, ELECTRONS=electrons, HORIZONTAL=horizontal, IONS=ions
    if JxB    then MrSim_OhmsLaw_JxB,    oSim, component, cut, ELECTRONS=electrons, HORIZONTAL=horizontal, IONS=ions
    if gradPe then MrSim_OhmsLaw_gradPe, oSim, component, cut, ELECTRONS=electrons, HORIZONTAL=horizontal, IONS=ions
    
;-------------------------------------------------------
;Output ////////////////////////////////////////////////
;-------------------------------------------------------
    nPlots = ETotal + VxB + JxB + gradPe
    switch nPlots of
        4: begin
            middle = ohmWin -> FindByPIndex(3)
            middle -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
        endswitch
        
        3: begin
            middle = ohmWin -> FindByPIndex(2)
            middle -> SetProperty, TITLE='', XTITLE='', XTICKFORMAT='(a1)'
        endswitch 
        
        2: begin
            top = ohmWin -> FindByPIndex(1)
            top -> SetProperty, XTITLE='', XTICKFORMAT='(a1)'
        
            bottom = ohmWin -> FindByPIndex(nPlots)
            bottom -> SetProperty, TITLE=''
        endswitch
        
        1: ;Do nothing
    endswitch

    ;Destroy the object
    if osim_created && arg_present(oSim) eq 0 then obj_destroy, oSim
    
    ;Refresh and output, if requested.
    if current eq 0 then begin
        ohmWin -> Refresh
        if ofilename ne '' then ohmWin -> Save, ofilename
    endif

    return, ohmWin
end
