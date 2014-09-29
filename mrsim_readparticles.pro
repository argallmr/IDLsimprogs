; docformat = 'rst'
;
; NAME:
;    MrSim_ReadParticles
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
;   Read particle data from binary files.
;
;   NOTES:
;       Reading binary file: Bill Daughton's Description of the file (u is momentum):
;
;           This is just raw binary (using 4 byte reals)
;
;               x,z,ux,uy,uz
;
;           for each particle. This just repeats for all the particles in the file.
;           Whatever tool you are using, just declare single precision (4 byte) real
;           numbers, open the file, and then begin reading the particles in within a loop.
;
;       Counts:
;           Bill states that only every-other particle is saved. When dealing with
;           electron counts for distribution functions, etc., counts must be multiplied
;           by two.
;
;
; :Categories:
;    Bill Daughton, Simulation
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
;       2014/08/28  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Examples of using MrSim_ReadParticles
;-
function MrSim_ReadParticles_Examples, example
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if obj_valid(win)  then obj_destroy, win
        if obj_valid(oSim) then obj_destroy, oSim
        void = cgErrorMSG()
        return, !Null
    endif

    case example of
    
;---------------------------------------------------------------------
; Histogrom of Particle Counts vs. Denstiy: Asymm-3D /////////////////
;---------------------------------------------------------------------
        1: begin
            ;Define in puts
            time     = 30
            xrange   = [420, 450]
            zrange   = [-10, 10]
            binsize  = [0.5, 0.5]
            filename = '/data2/Asymm-Scan/By1/electrons-twci60-it70200.bin'

            ;Histogram particle data.
            data = MrSim_ReadParticles(filename, xrange, zrange, /VERBOSE, $
                                       BINSIZE=binSize, XLOC=x, ZLOC=z)

            ;Create a window
            win = MrWindow(OXMARGIN=[10,15], NAME='e- Counts vs. Density.')

            ;Plot an image of the histogram.
            img   = MrImage(data, x, z, /AXES, /SCALE, CTINDEX=13, /CURRENT, /LOG, $
                            XTITLE='x (de)', YTITLE='z (de)', TITLE='Particle Counts', $
                            MISSING_VALUE=0)
            !Null = MrColorbar(TARGET=img, TITLE='Counts', /CURRENT)

            ;Compare it with density.
            oSim  = MrSim_Create('Asymm-Scan/By1', time, XRANGE=xrange, ZRANGE=zrange)
            !Null = MrSim_ColorSlab('ne', /CURRENT, SIM_OBJECT=oSim)
            obj_destroy, oSim
            
            return, win
        endcase
        
;---------------------------------------------------------------------
; Histogrom of Particle Counts vs. Denstiy: Asymm-Large-2D-NEW ///////
;---------------------------------------------------------------------
        2: begin
            ;Define inputs
            time     = 18
            xrange   = 1473 + [-40,40]
            zrange   = [-25, 10]
            binsize  = [0.5, 0.5]
            filename = '/data2/Asymm-Large-2D-NEW/electrons-t' + strtrim(5*time, 2) + '.bin'

            ;Read the particle data
            data = MrSim_ReadParticles(filename, xrange, zrange, FMAP_DIR='', /VERBOSE)
            data = data[0:1,*]
            dims = size(data, /DIMENSIONS)

            ;Histogram the data
            hist = hist_nd(data, binsize, MIN=[xrange[0], zrange[0]], MAX=[xrange[1], zrange[1]])

            ;Create the data arrays
            x = linspace(xrange[0], xrange[1], binsize[0], /INTERVAL)
            z = linspace(zrange[0], zrange[1], binsize[1], /INTERVAL)

            ;Create a window
            win2 = MrWindow(XSIZE=500, YSIZE=550, OXMARGIN=[10,15])

            ;Plot the histogram
            img   = MrImage(hist, x, z, /AXES, /SCALE, CTINDEX=13, /CURRENT, $
                            XTITLE='x (de)', YTITLE='z (de)', TITLE='Particle Counts')
            !Null = MrColorbar(TARGET=img, TITLE='Counts', /CURRENT)

            ;Get rid of all the data
            data  = !Null
            x     = !Null
            z     = !Null

            ;
            ; Plot density.
            ;
            dir   = '/data2/Asymm-Large-2D-NEW/'
            oSim  = obj_new('MrSim2D', time, XRANGE=xrange, ZRANGE=zrange, $
                           DIRECTORY=dir+'data', INFO_FILE=dir+'info')
            !Null = MrSim_ColorSlab('ne', /CURRENT, SIM_OBJECT=oSim)
            obj_destroy, oSim

            win2 -> Refresh
        endcase
        
;---------------------------------------------------------------------
; Histogrom of Particle Counts vs. Denstiy: Asymm-3D /////////////////
;---------------------------------------------------------------------
        3: begin
            ;Total size
            ;   TIME   = 108090
            ;   YSLICE = 650, 905, 1400
            ;   XRANGE = [250, 700]
            ;   ZRANGE = [-80,  80]
            t        = 108090
            yslice   = 650
            xrange   = [400, 600]
            zrange   = [-10,  10]
            binsize  = [0.5, 0.5]
            filename = '/home/daughton/Asymm-3D/electrons-y' + strtrim(yslice, 2) + '.bin'

            ;Read the particle data
            data = MrSim_ReadParticles(filename, xrange, zrange, /VERBOSE)
            data = data[0:1,*]
            dims = size(data, /DIMENSIONS)

            ;Histogram the data
            hist = hist_nd(data, binsize, MIN=[xrange[0], zrange[0]], MAX=[xrange[1], zrange[1]])

            ;Create the data arrays
            x = linspace(xrange[0], xrange[1], binsize[0], /INTERVAL)
            z = linspace(zrange[0], zrange[1], binsize[1], /INTERVAL)

            ;Create a window
            win = MrWindow(XSIZE=500, YSIZE=550, OXMARGIN=[10,15])

            ;Plot the histogram
            img   = MrImage(hist, x, z, /AXES, /SCALE, CTINDEX=13, /CURRENT, $
                            XTITLE='x (de)', YTITLE='z (de)', TITLE='Particle Counts')
            !Null = MrColorbar(TARGET=img, TITLE='Counts', /CURRENT)

            ;Get rid of all the data
            data  = !Null
            x     = !Null
            z     = !Null

            ;
            ; Plot density.
            ;
            dir    = '/home/daughton/Asymm-3D/'
            oSim   = obj_new('MrSim3D', time, yslice, XRANGE=xrange, ZRANGE=zrange, $
                             DIRECTORY=dir+'data', INFO_FILE=dir+'info')
            !Null   = MrSim_ColorSlab('ne', /CURRENT, SIM_OBJECT=oSim)
            obj_destroy, oSim

            win -> Refresh
        endcase
        
;---------------------------------------------------------------------
; Use fMap file //////////////////////////////////////////////////////
;---------------------------------------------------------------------
        else: message, 'Invalid example choise: ' + strtrim(example, 2) + '.'
    endcase
end


;+
;   Read particle data.
;
; :Params:
;       FILENAME:           in, required, type=string
;                           Name of the file containing particle data.
;       XRANGE:             in, required, type=fltarr(2)
;                           X-range (in de) over which particle data is to be kept.
;       ZRANGE:             in, required, type=fltarr(2)
;                           Z-range (in de) over which particle data is to be kept.
;
; :Keywords:
;       BINSIZE:            in, optional, type=integer/intarr(2)
;                           If provided, electrons will be histogrammed into spacial bins
;                               of this size (units of de), and `DATA` will contain the
;                               particle counts in each bin. If a scalar, the x- and z-
;                               dimensions will have the same binsize.
;       FMAP_DIR:           in, optional, type=string, default=pwd
;                           Directory in which to find an fMap. See MrSim_Create_fMap.pro.
;       N_RECS_PER_CHUNK:   in, optional, type=ulong64, default=1000000ULL
;                           Number of records per data chunk.
;       NBINS:              in, optional, type=integer/intarr(2)
;                           If provided, electrons will be histogrammed into this number
;                               of spacial bins, and `DATA` will contain the particle
;                               counts in each bin. If a scalar, the same number of bins
;                               will be used for x and z.
;       REC_SAMPLE:         in, optional, type=any, default=fltarr(5)
;                           An example of how records are to be read. The default assumes
;                               one record consists of five 4-byte (16-bit) floating point
;                               data values of the form [x, z, ux, uy, uz], where x and
;                               z are two spatial locations and u[xyz] are the three
;                               components of the momentum.
;       VERBOSE:            in, optional, type=boolean, default=0
;                           If set, information about particle will be printed data to
;                               the command window.
;       XPOS:               out, optional, type=fltarr
;                           If `BINSIZE` or `NBINS` is given, then XPOS is the spacial
;                               location of the center of each column in `DATA`.
;       ZPOS:               out, optional, type=fltarr
;                           If `BINSIZE` or `NBINS` is given, then ZPOS is the spacial
;                               location of the center of each row in `DATA`.
;
; :Returns:
;       DATA:               MrWindow graphic window containing the requested graphics.
;                                   If the image data does not exist, an invalid object
;                                   will be returned.
;-
function MrSim_ReadParticles, filename, xrange, zrange, $
BINSIZE=binSize, $
FMAP_DIR=fMap_dir, $
N_RECS_PER_CHUNK=n_rec_per_chunk, $
NBINS=nBins, $
REC_SAMPLE=rec_sample, $
VERBOSE=verbose, $
XLOC=xloc, $
ZLOC=zloc
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif

    ;Make sure the file exists
    if file_test(filename) eq 0 then $
        message, 'File does not exist: "' + filename + '".'

    ;Defaults
    ;   - There are five 4-byte numbers associated with each particle.
    verbose = keyword_set(verbose)
    if n_elements(fMap_dir)         eq 0 then fMap_dir         = file_dirname(filename)
    if n_elements(n_recs_per_chunk) eq 0 then n_recs_per_chunk = 1000000ULL
    if n_elements(rec_sample)       eq 0 then rec_sample       = fltarr(5)
    
    ;Make a histogram?
    if n_elements(binsize) + n_elements(nBins) gt 0 $
        then make_hist = 1 $
        else make_hist = 0
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
    
    ;Open the file and get its info
    openr, lun, filename, /GET_LUN
    finfo = fstat(lun)  

;---------------------------------------------------------------------
; How to Read File ///////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Check for a save file
    fbase = cgRootName(filename)
    ftest = filepath(fbase + '_mrfMap.sav', ROOT_DIR=fMap_dir)
    if file_test(ftest) then begin
        restore, ftest
        iChunks = where(map_entry.d1_range[0] lt xrange[1] and $
                        map_entry.d1_range[1] gt xrange[0] and $
                        map_entry.d2_range[0] lt zrange[1] and $
                        map_entry.d2_range[1] gt zrange[0],    $
                        n_chunks)
        
        n_per_chunk  = map_entry[iChunks].nRecs
        byte_offsets = map_entry[iChunks].pos
        map_entry    = !Null
        
    ;Determine how many chunks to use.
    endif else begin
        message, 'fMap not found: "' + ftest + '".', /INFORMATIONAL
    
        ;Number and size of chunks
        n_recs       = ulong64(finfo.size / rec_size)       ; number of records in the file
        n_chunks     = n_recs / n_recs_per_chunk            ; integer divide
        n_last_chunk = n_recs mod n_recs_per_chunk
        if n_last_chunk ne 0 then ++n_chunks
        
        iChunks      = lindgen(n_chunks)
        byte_offsets = iChunks * n_recs_per_chunk * rec_size
        n_per_chunk  = replicate(n_recs_per_chunk, n_chunks)
        if n_last_chunk gt 0 then n_per_chunk[-1] = n_last_chunk
    endelse

    ;Allocate memory
    region  = fltarr(n_per_rec,   n_per_chunk[0])
    data    = fltarr(n_per_rec, 2*n_per_chunk[0])
    iStart  = 0ULL
    nEmpty  = 2*n_per_chunk[0]

;---------------------------------------------------------------------
; Read Each Chunk ////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Loop through each data chunk
    ;   - Read data in chunks so that all 9GB are not in memory at once.
    for i = 0ULL, n_chunks - 1ULL do begin
        if verbose then if ((i+1) mod 10 eq 0) || (i eq n_chunks - 1) $
            then print, FORMAT='(%"Chunk %i of %i")', i+1, n_chunks
        
        ;Last chunk may be a different size
        if i gt 0 then if n_per_chunk[i] ne n_per_chunk[i-1] $
            then region = fltarr(n_per_rec, n_per_chunk[i])
    
        ;Read the data
        point_lun, lun, byte_offsets[i]
        readu, lun, region

        ;Are we within the range?
        iGood = where(region[0,*] ge xrange[0] and $
                      region[0,*] le xrange[1] and $
                      region[1,*] ge zrange[0] and $
                      region[1,*] le zrange[1], nGood)
        
        ;Were there good points?
        if nGood gt 0 then begin
            ;Expand by two chunks
            if nGood gt nEmpty then begin
                nExpand = nGood > 2*n_per_chunk[i]
                data    = [[data], [fltarr(n_per_rec, nExpand)]]
                nEmpty  = nEmpty + nExpand
            endif
        
            ;Fill in the data
            iStop = iStart + nGood - 1
            data[*, iStart:iStop] = region[*, iGood]
            
            ;Proceed to the next chunk
            nEmpty -= nGood
            iStart  = iStop + 1
        endif
    endfor

    ;Close the file
    free_lun, lun
    data = data[*,0:iStart-1]

;---------------------------------------------------------------------
; Histogram into Particle Counts? ////////////////////////////////////
;---------------------------------------------------------------------
    if make_hist then begin
        ;We only need the spacial locations
        ;   - Throw away momentum
        data = data[0:1,*]
    
        ;Make a histrogram.
        data = hist_nd(data, binSize, NBINS=nBins, $
                       MIN=[xrange[0], zrange[0]], $
                       MAX=[xrange[1], zrange[1]])

        ;Position of the histogram bin centers.
        if arg_present(xloc) then xloc = linspace(xrange[0], xrange[1], binsize[0], /INTERVAL)
        if arg_present(zloc) then zloc = linspace(zrange[0], zrange[1], binsize[1], /INTERVAL)
    endif
    
    return, data
end



;---------------------------------------------------------------------
; Main-Level Example Program: IDL> .r MrSim_ReadParticles ////////////
;---------------------------------------------------------------------

;Choose example 1-3
example = 1
win = MrSim_ReadParticles_Examples(example)

end