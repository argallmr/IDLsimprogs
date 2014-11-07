; docformat = 'rst'
;
; NAME:
;    MrSim_Create_fMap
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
;   Create a map into the binary particles files.
;
;   NOTE:
;       Bill Daughton's Description of the file (u is momentum):
;
;           This is just raw binary (using 4 byte reals)
;
;               x, z, ux, uy, uz
;
;           for each particle. This just repeats for all the particles in the file.
;           Whatever tool you are using, just declare single precision (4 byte) real
;           numbers, open the file, and then begin reading the particles in within a loop.
;
; :Examples:
;   See the main level program at the end of this document::
;       IDL> .r MrSim_Create_fMap
;
;   For an example of how to use the fMap, see MrSim_ReadParticles.pro
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       FILENAME:           in, required, type=string
;                           Name of the file containing particle data.
;
; :Keywords:
;       N_RECS_PER_CHUNK:   in, optional, type=ulong64, default=1000000ULL
;                           Number of records per data chunk.
;       RANGE1:             out, optional, type=fltarr(2)
;                           Range, in data coordinates, of the first spacial dimension.
;       RANGE2:             out, optional, type=fltarr(2)
;                           Range, in data coordinates, of the second spacial dimension.
;       RANGE3:             out, optional, type=fltarr(2)
;                           Range, in data coordinates, of the third spacial dimension.
;       REC_SAMPLE:         in, optional, type=any, default=fltarr(5)
;                           An example of how records are to be read. The default assumes
;                               one record consists of five 4-byte (16-bit) floating point
;                               data values of the form [x, z, ux, uy, uz], where x and
;                               z are two spatial locations and u[xyz] are the three
;                               components of the momentum.
;       SAVE_DIR:           in, optional, type=string, default=pwd
;                           Directory in which the fMap will be saved.
;       SIM3D:              in, optional, type=boolean, default=0
;                           If set, position coordinates are assumed to be provided in
;                               3 dimensions: x, y, and z.
;       VERBOSE:            in, optional, type=boolean, default=0
;                           If set, information about particle will be printed data to
;                               the command window.
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
;       2014/09/05  -   Written by Matthew Argall
;       2014/10/30  -   Added the SIM3D keyword. - MRA
;       2014/11/06  -   Added the OVERWRITE keyword and changed the output file name to
;                           create fewer conflicts. - MRA
;-
pro MrSim_Create_fMap, filename, $
N_RECS_PER_CHUNK=n_rec_per_chunk, $
OVERWRITE=overwrite, $
RANGE1=range1, $
RANGE2=range2, $
RANGE3=range3, $
REC_SAMPLE=rec_sample, $
SAVE_DIR=save_dir, $
SIM3D=sim3d, $
VERBOSE=verbose
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return
    endif

    ;Defaults
    ;   - There are five 4-byte numbers associated with each particle.
    overwrite = keyword_set(overwrite)
    sim3d     = keyword_set(sim3d)
    verbose   = keyword_set(verbose)
    if n_elements(n_recs_per_chunk) eq 0 then n_recs_per_chunk = 1000000ULL
    if n_elements(rec_sample)       eq 0 then rec_sample       = sim3d ? fltarr(6) : fltarr(5)
    if n_elements(save_dir)         eq 0 then void             = cgRootName(filename, DIRECTORY=save_dir)
    
    ;Make sure the file exists
    if file_test(filename) eq 0 then $
        message, 'File does not exist: "' + filename + '".'
    
    ;Form the output file name
    fbase    = cgRootName(filename, DIRECTORY=directory)
    sname    = idl_validname(file_basename(directory), /CONVERT_ALL)
    file_out = filepath('MrFMap' + (sim3d ? '3d' : '2d') + '_' + sname + '_' + fbase + '.sav', ROOT_DIR=save_dir)
    
    ;Check if the output file exists already
    if file_test(file_out) then begin
        if ~overwrite then $
            message, 'Output file name already exists. Use /OVERWRITE. "' + file_out + '".'
    endif

;---------------------------------------------------------------------
;Record Info /////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    n_per_rec = n_elements(rec_sample)
    tname     = size(rec_sample, /TNAME)
    type      = size(rec_sample, /TYPE)
    
    ;Bytes for IDL data types
    case tname of
        'BYTE':    nBytes = 1
        'INT':     nBytes = 1
        'LONG':    nBytes = 8
        'FLOAT':   nBytes = 4
        'DOUBLE':  nBytes = 8
        'UINT':    nBytes = 2
        'ULONG':   nBytes = 8
        'LONG64':  nBytes = 16
        'ULONG64': nBytes = 16
        else: message, 'Data type "' + tname + '" not recognized.'
    endcase
        
    ;Size of the record
    rec_size = n_per_rec * nBytes
    
;---------------------------------------------------------------------
; Open File //////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    
    ;Open the file
    openr, lun, filename, /GET_LUN
    finfo       = fstat(lun)
    n_recs      = ulong64(finfo.size / rec_size)       ; number of records in the file
    n_chunks    = n_recs / n_recs_per_chunk            ; integer divide
    n_recs_last = n_recs mod n_recs_per_chunk
    if n_recs_last ne 0 then ++n_chunks

    ;Print information about the file
    if verbose then begin
        print, 'Information about ' + filename + ':'
        print, FORMAT='(%"   # per rec          : %i")', n_per_rec
        print, FORMAT='(%"   record type        : %s")', tname
        print, FORMAT='(%"   record size        : %i")', rec_size
        print, FORMAT='(%"   filesize           : %i")', finfo.size
        print, FORMAT='(%"   #records in file   : %i")', n_recs
        print, FORMAT='(%"   #records per chunk : %i")', n_recs_per_chunk
        print, FORMAT='(%"   #chunks to read    : %i")', n_chunks
        if n_recs_last gt 0 then $
            print, FORMAT='(%"   last chunk size    : %i")', n_recs_last
    endif

;---------------------------------------------------------------------
; Allocate Memory ////////////////////////////////////////////////////
;---------------------------------------------------------------------
    map_entry = { pos:      0ULL, $
                  nRecs:    0ULL, $
                  d1_range: fltarr(2), $
                  d2_range: fltarr(2), $
                  d3_range: fltarr(2) $
                }
    map_entry = replicate(map_entry, n_chunks)
    data      = make_array(n_per_rec, n_recs_per_chunk, TYPE=type)
    
;---------------------------------------------------------------------
; Read Chunks ////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Save the maximum ranges
    range1 = [!values.f_infinity, -!values.f_infinity]
    range2 = [!values.f_infinity, -!values.f_infinity]
    range3 = [!values.f_infinity, -!values.f_infinity]
    
    ;Loop through each data chunk
    ;   - Read data in chunks so that all xxGB are not in memory at once.
    n_recs_read = 0ULL
    for i = 0ULL, n_chunks - 1ULL do begin
        ;The last chunk might not be a complete chunk.
        ;   - Allocate the correct amount of memory.
        if i eq n_chunks - 1 && n_recs_last gt 0 then begin
            n_recs_per_chunk = n_recs_last
            data = make_array(n_per_rec, n_recs_last, TYPE=type)
        endif
        map_entry[i].nRecs = n_recs_per_chunk
        
        ;Store the current file position
        point_lun, -lun, filepos
        map_entry[i].pos = filepos

        ;Print info
        if verbose then if ((i+1) mod 10 eq 0) || (i eq n_chunks-1) then begin
            print, FORMAT='(%"Chunk %i of %i. Bytes %i-%i of %i")', $
                   i+1, n_chunks, n_recs_read*rec_size, (n_recs_read+n_recs_per_chunk)*rec_size, finfo.size
        endif
        n_recs_read += n_recs_per_chunk

        ;Read the data
        readu, lun, data
        
        ;3D electron data?
        if sim3d then begin
            ;Get the minimum and maximum values in the two spacial dimensions
            d1_min = min(data[0,*], MAX=d1_max)
            d2_min = min(data[1,*], MAX=d2_max)
            d3_min = min(data[2,*], MAX=d3_max)
            map_entry[i].d1_range = [d1_min, d1_max]
            map_entry[i].d2_range = [d2_min, d2_max]
            map_entry[i].d3_range = [d3_min, d3_max]
        
            ;Record the domain being read
            range1[0] <= d1_min
            range1[1] >= d1_max
            range2[0] <= d2_min
            range2[1] >= d2_max
            range3[0] <= d3_min
            range3[1] >= d3_max
        
        ;2D electron data
        endif else begin
            ;Get the minimum and maximum values in the two spacial dimensions
            d1_min = min(data[0,*], MAX=d1_max)
            d3_min = min(data[1,*], MAX=d3_max)
            map_entry[i].d1_range = [d1_min, d1_max]
            map_entry[i].d3_range = [d3_min, d3_max]
        
            ;Record the domain being read
            range1[0] <= d1_min
            range1[1] >= d1_max
            range2[0]  = 0
            range2[1]  = 0
            range3[0] <= d3_min
            range3[1] >= d3_max
        endelse
    endfor

    ;Close the file
    free_lun, lun

    ;Display the range
    if verbose then begin
        print, FORMAT='(%"Range1 = [%f, %f]")', range1
        print, FORMAT='(%"Range2 = [%f, %f]")', range2
        print, FORMAT='(%"Range3 = [%f, %f]")', range3
    endif

    ;Save data
    save, map_entry, FILENAME=file_out, DESCRIPTION='Map into binary electron distribution files.'
    print, 'File saved to: "' + file_out + '".'
end



;---------------------------------------------------------------------
; Main-Level Example Program: IDL> .r MrSim_Create_fMap //////////////
;---------------------------------------------------------------------
;Files and directories
filename = '/data2/Asymm-2D-Large/electrons-t90.bin'
save_dir = '/home/argall/simulations/fmaps/'

;Read the particle data
MrSim_Create_fMap, filename, SAVE_DIR=save_dir, /VERBOSE

end