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
;       2014/10/14  -   Added the VELOCITY and ENERGY keywords. - MRA
;       2014/10/30  -   Added DIST3D keyword. - MRA
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
; Histogrom of Particle Counts vs. Denstiy: Asymm-Scan/By1 ///////////
;---------------------------------------------------------------------
        1: begin
            ;Define in puts
            theSim   = 'Asymm-Scan/By1'
            time     = 30
            xrange   = [420, 450]
            zrange   = [-10, 10]
            binsize  = [0.5, 0.5]
            MrSim_Which, theSim, EFILE=eFile, FMAP_DIR=fMap_dir, TINDEX=tIndex

            ;Histogram particle data.
            data = MrSim_ReadParticles(eFile, xrange, yrange, zrange, FMAP_DIR=fMap_dir, /VERBOSE)

            ;Create a window
            win = MrWindow(OXMARGIN=[10,15], NAME='e- Counts vs. Density.')

            ;Plot an image of the histogram.
            img   = MrImage(data, x, z, /AXES, /SCALE, CTINDEX=13, /CURRENT, /LOG, $
                            XTITLE='x (de)', YTITLE='z (de)', TITLE='Particle Counts', $
                            MISSING_VALUE=0, NAME='Counts')
            !Null = MrColorbar(TARGET=img, TITLE='Counts', NAME='CB: Counts', /CURRENT)

            ;Compare it with density.
            oSim  = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange)
            !Null = MrSim_ColorSlab(oSim, 'ne', /CURRENT)
            obj_destroy, oSim
            
            return, win
        endcase
        
;---------------------------------------------------------------------
; Histogrom of Particle Counts vs. Denstiy: Asymm-Large-2D-NEW ///////
;---------------------------------------------------------------------
        2: begin
            ;Define inputs
            theSim   = 'Asymm-Large-2D-NEW'
            tIndex   = 18
            xrange   = 1473 + [-40,40]
            zrange   = [-25, 10]
            binsize  = [0.5, 0.5]
            MrSim_Which, theSim, EFILE=eFile, FMAP_DIR=fMap_dir, TINDEX=tIndex

            ;Read the particle data
            data = MrSim_ReadParticles(eFile, xrange, yrange, zrange, FMAP_DIR=fMap_dir, /VERBOSE)
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
            oSim  = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange)
            !Null = MrSim_ColorSlab(oSim, 'ne', /CURRENT)
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
            theSim   = 'Asymm-3D'
            t        = 108090
            yslice   = 650
            xrange   = [400, 600]
            zrange   = [-10,  10]
            binsize  = [0.5, 0.5]
            MrSim_Which, theSim, EFILE=eFile, FMAP_DIR=fMap_dir, TINDEX=tIndex

            ;Read the particle data
            data = MrSim_ReadParticles(eFile, xrange, zrange, FMAP_DIR=fMap_dir, /VERBOSE)
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
                            XTITLE='x (de)', YTITLE='z (de)', TITLE='Particle Counts', $
                            NAME='Counts')
            !Null = MrColorbar(TARGET=img, TITLE='Counts', NAME='CB: Counts', /CURRENT)

            ;Get rid of all the data
            data  = !Null
            x     = !Null
            z     = !Null

            ;
            ; Plot density.
            ;
            oSim  = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange)
            !Null = MrSim_ColorSlab(oSim, 'ne', /CURRENT)
            obj_destroy, oSim

            win -> Refresh
        endcase
        
;---------------------------------------------------------------------
; Histogrom of Particle Counts vs. Denstiy: Asymm-3D /////////////////
;---------------------------------------------------------------------
        4: begin
            theSim   = 'data-by0.03-NEW'
            tIndex   = 22
            yrange   = [0, 0]
            
            ;Particle domains for the various times
            case tIndex of
                17: begin
                    xrange   = 807 + [-30, 30]
                    zrange   = [-15,  15]
                endcase
                18: begin
                    xrange   = 806 + [-30, 30]
                    zrange   = [-15,  15]
                endcase
                19: begin
                    xrange   = 807 + [-40, 40]
                    zrange   = [-15,  15]
                endcase
                22: begin
                    xrange   = 799 + [-60, 60]
                    zrange   = [-15,  15]
                endcase
                29: begin
                    xrange   = 771 + [-30, 30]
                    zrange   = [-15,  15]
                endcase
                else: message, 'No particle data for tIndex ' + strtrim(tIndex, 2) + '.'
            endcase
            
            binsize  = [0.5, 0.5]
            MrSim_Which, theSim, EFILE=eFile, FMAP_DIR=fMap_dir, TINDEX=tIndex

            ;Read the particle data
            data = MrSim_ReadParticles(eFile, xrange, yrange, zrange, FMAP_DIR=fMap_dir, /VERBOSE)
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
                            XTITLE='x (de)', YTITLE='z (de)', TITLE='Particle Counts', $
                            NAME='Counts')
            !Null = MrColorbar(TARGET=img, TITLE='Counts', NAME='CB: Counts', /CURRENT)

            ;Get rid of all the data
            data  = !Null
            x     = !Null
            z     = !Null

            ;
            ; Plot density.
            ;
            oSim  = MrSim_Create(theSim, tIndex, XRANGE=xrange, ZRANGE=zrange, /BINARY)
            !Null = MrSim_ColorSlab(oSim, 'ne', /CURRENT)
            obj_destroy, oSim

            win -> Refresh
        endcase
        
;---------------------------------------------------------------------
; Use fMap file //////////////////////////////////////////////////////
;---------------------------------------------------------------------
        else: message, 'Invalid example choice: ' + strtrim(example, 2) + '.'
    endcase
    
    return, win
end


;+
;   Read particle data.
;
;   CALLING SEQUENCE:
;       data = MrSim_ReadParticles(filename, xrange, zrange)
;       data = MrSim_ReadParticles(filename, xrange, yrange, zrange)
;
; :Params:
;       FILENAME:           in, required, type=string
;                           Name of the file containing particle data.
;       XRANGE:             in, required, type=fltarr(2)
;                           X-range (in de) over which particle data is to be kept.
;       YRANGE:             in, optional, type=fltarr(2), default=[0,0]
;                           Y-range (in de) over which particle data is to be kept.
;       ZRANGE:             in, required, type=fltarr(2)
;                           Z-range (in de) over which particle data is to be kept.
;
; :Keywords:
;       BINSIZE:            in, optional, type=integer/intarr(2)
;                           If provided, electrons will be histogrammed into spacial bins
;                               of this size (units of de), and `DATA` will contain the
;                               particle counts in each bin. If a scalar, the x- and z-
;                               dimensions will have the same binsize.
;       DIST3D:             in, optioanl, type=boolean, default=0
;                           If set, the data file is ordered as [x, y, z, ux, uy, uz], i.e.
;                               it has 3D coordinates in both space and momentum.
;       ENERGY:             in, optional, type=boolean, default=0
;                           If set, momentum will be converted to energy.
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
;       VELOCITY:           in, optional, type=boolean, default=0
;                           If set, momentum will be converted to velocity.
;       VERBOSE:            in, optional, type=boolean, default=0
;                           If set, information about particle will be printed data to
;                               the command window.
;       XLOC:               out, optional, type=fltarr
;                           If `BINSIZE` or `NBINS` is given, then XLOC is the spacial
;                               location of the center of each column in `DATA`.
;       ZLOC:               out, optional, type=fltarr
;                           If `BINSIZE` or `NBINS` is given, then ZLOC is the spacial
;                               location of the center of each row in `DATA`.
;
; :Returns:
;       DATA:               Electron distribution data.
;-
function MrSim_ReadParticles, filename, xrange, yrange, zrange, $
BINSIZE=binSize, $
DIST3D=dist3D, $
ENERGY=energy, $
FMAP_DIR=fMap_dir, $
N_RECS_PER_CHUNK=n_rec_per_chunk, $
NBINS=nBins, $
REC_SAMPLE=rec_sample, $
VELOCITY=velocity, $
VERBOSE=verbose, $
XLOC=xloc, $
ZLOC=zloc
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(lun) gt 0 then free_lun, lun
        if nParams eq 3 then yrange = zrange
        void = cgErrorMSG()
        return, !Null
    endif

;---------------------------------------------------------------------
; Inputs /////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    nParams = n_params()
    if nParams eq 3 then begin
        zrange = yrange
        yrange = [0,0]
    endif

    ;Make sure the file exists
    ;   - If it does not, open a dialog.
    ;   - If cancelled, return
    if n_elements(filename) eq 0 then filename = ''
    if file_test(filename) eq 0 then begin
        if filename ne '' then $
            message, 'Choose a file. Given file does not exist. "' + filename + '".', /INFORMATIONAL
        if filename eq '' then filename = dialog_pickfile(/READ, TITLE='Pick electron file.')
        if filename eq '' then return, !Null
    endif

    ;Defaults
    ;   - There are five 4-byte numbers associated with each particle.
    energy   = keyword_set(energy)
    dist3D   = keyword_set(dist3D)
    velocity = keyword_set(velocity)
    verbose  = keyword_set(verbose)
    if n_elements(fMap_dir)         eq 0 then fMap_dir         = file_dirname(filename)
    if n_elements(n_recs_per_chunk) eq 0 then n_recs_per_chunk = 1000000ULL
    if n_elements(rec_sample)       eq 0 then rec_sample       = dist3D ? fltarr(6) : fltarr(5)
    if n_elements(yrange)           eq 0 then yrange           = [0, 0]
    
    ;Make a histogram?
    if n_elements(binsize) + n_elements(nBins) gt 0 $
        then make_hist = 1 $
        else make_hist = 0
    
    ;Dependencies
    if energy + velocity gt 1 then message, 'ENERGY and VELOCITY are mutually exclusive.'
;---------------------------------------------------------------------
; Record Info ////////////////////////////////////////////////////////
;---------------------------------------------------------------------
    ;Number and type of records.
    n_per_rec = n_elements(rec_sample)
    tname     = size(rec_sample, /TNAME)
    type      = size(rec_sample, /TYPE)
    
    ;Bytes of storage used by IDL data types
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
    ;fMap?
    fbase = cgRootName(filename, DIRECTORY=directory)
    sname = idl_validname(file_basename(directory), /CONVERT_ALL)
    ftest = filepath('MrFMap' + (dist3D ? '3d' : '2d') + '_' + sname + '_' + fbase + '.sav', ROOT_DIR=fMap_dir)
    if file_test(ftest) then begin
        restore, ftest
        iChunks = where((map_entry.d1_range[0] le xrange[1]  and $
                         map_entry.d1_range[1] ge xrange[0]) and $
                        (map_entry.d2_range[0] le yrange[1]  and $
                         map_entry.d2_range[1] ge yrange[0]) and $
                        (map_entry.d3_range[0] le zrange[1]  and $
                         map_entry.d3_range[1] ge zrange[0]),    $
                        n_chunks)
        
        ;No chunks found
        if n_chunks eq 0 then begin
            msg = string(FORMAT='(%"No particles found in range x=[%i, %i], y=[%i, %i], z=[%i, %i]")', $
                         xrange, yrange, zrange)
            message, msg
        endif
        
        n_per_chunk  = map_entry[iChunks].nRecs
        byte_offsets = map_entry[iChunks].pos
        map_entry    = !Null
        
    ;Read through the entire file.
    endif else begin
        message, 'fMap not found: "' + ftest + '".', /INFORMATIONAL
    
        ;Number and size of data chunks to read per loop
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
        if dist3D then begin
            iGood = where(region[0,*] ge xrange[0] and $
                          region[0,*] le xrange[1] and $
                          region[1,*] ge yrange[0] and $
                          region[1,*] le yrange[1] and $
                          region[2,*] ge zrange[0] and $
                          region[2,*] le zrange[1], nGood)
        endif else begin
            iGood = where(region[0,*] ge xrange[0] and $
                          region[0,*] le xrange[1] and $
                          region[1,*] ge zrange[0] and $
                          region[1,*] le zrange[1], nGood)
        endelse
        
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
    
    nParticles = iStart
    data = data[*,0:iStart-1]

;---------------------------------------------------------------------
; Convert Units //////////////////////////////////////////////////////
;---------------------------------------------------------------------

    ;Convert to velocity
    ;   u     = gamma * v
    ;   gamma = lorentz factor: sqrt( 1 - v^2 / c^2)^-1
    if velocity then begin
        if dist3D then begin
            e_gamma = sqrt(1.0 + total(data[3:5, *]^2, 1))
            data[2:4, *] /= rebin(reform(temporary(e_gamma), 1, nParticles), 3, nParticles)
        endif else begin
            e_gamma = sqrt(1.0 + total(data[2:4, *]^2, 1))
            data[2:4, *] /= rebin(reform(temporary(e_gamma), 1, nParticles), 3, nParticles)
        endelse
    endif
    
    ;Convert to energy
    if energy then begin
        message, 'Converting to energy not implemented.'
    endif

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
    
    if nParams eq 3 then yrange = zrange
    return, data
end



;---------------------------------------------------------------------
; Main-Level Example Program: IDL> .r MrSim_ReadParticles ////////////
;---------------------------------------------------------------------

;Choose example 1-4
example = 4
win = MrSim_ReadParticles_Examples(example)

end