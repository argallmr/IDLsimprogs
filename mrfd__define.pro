; docformat = 'rst'
;
; NAME:
;    MrLine3D
;
; PURPOSE:
;+
;   Create a series of point along a line in 3D space, given to points that define the line.
;
; :Categories:
;    Bill Daughton, Simulation
;
; :Params:
;       R0:             in, required, type=fltarr(3)
;                       A point in cartesian space through which the line should pass.
;       R1:             in, required, type=fltarr(3)
;                       A point in cartesian space through which the line should pass.
;
; :Keywords:
;       NPOINTS:        in, optional, type=integer, default=100
;                       Number of points to pass between `X0` and `X1`.
;
; :Returns:
;       LINE:           out, required, type=flaot(3\,`NPOINTS`)
;                       Points along the line connecting `X0` and `X1`.
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
;       2014/10/15  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Retrieve data from the info file.
;-
function MrFD::GetProbe, x, y, z
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;Pick out the grid
    coords = *self.coords
    xgrid  = coords[0, uniq(coords[0,*], sort(coords[0,*]))]
    ygrid  = coords[1, uniq(coords[1,*], sort(coords[1,*]))]
    zgrid  = coords[2, uniq(coords[2,*], sort(coords[2,*]))]
    
    ;Get the cell locations in the x, y, and z-directions
    grid   = value_locate(xgrid, x) > 0
    ypoint = value_locate(ygrid, y) > 0
    probe  = value_locate(zgrid, z) > 0

    ;Return the cell
    return, [grid, ypoint, probe]
end


;+
;   Retrieve object properties.
;-
pro MrFD::GetProperty, $
COORDS=coords, $
DIRECTORY=directory, $
INFO_FILE=info_file, $
NTIMES=nTimes
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Retrieve properties
    if arg_present(directory) then directory =  self.directory
    if arg_present(coords)    then coords    = *self.coords
    if arg_present(info_file) then info_file =  self.info_file
    if arg_present(nTimes)    then nTimes    =  self.nTimes
end


;+
;   Retrieve data from the info file.
;-
pro MrFD::GetInfo, $
EY_SAVE_INTERVAL=ey_save_interval, $
LX=Lx, $
LY=Ly, $
LZ=Lz, $
MI_ME=mi_me, $
NX=nx, $
NY=ny, $
NZ=nz, $
SAVE_DT=save_dt, $
TOPOLOGY_X=topology_x, $
TOPOLOGY_Y=topology_y, $
TOPOLOGY_Z=topology_z, $
VTHI=vthi, $
VTHE=vthe, $
WPE_WCE=wpe_wce
    compile_opt strictarr
    
    ;Error handling
    catch, the_error
    if the_error ne 0 then begin
        catch, /cancel
        void = cgErrorMSG()
        return
    endif
    
    ;Binary Info 
    if arg_present(topology_x)       then topology_x       = (*self.info).topology_x
    if arg_present(topology_y)       then topology_y       = (*self.info).topology_y
    if arg_present(topology_z)       then topology_z       = (*self.info).topology_z
    if arg_present(ey_save_interval) then ey_save_interval = (*self.info).ey_save_interval
    if arg_present(nx)               then nx               = (*self.info).nx
    if arg_present(ny)               then ny               = (*self.info).ny
    if arg_present(nz)               then nz               = (*self.info).nz
    if arg_present(Lx)               then Lx               = (*self.info).Lx
    if arg_present(Ly)               then Ly               = (*self.info).Ly
    if arg_present(Lz)               then Lz               = (*self.info).Lz
    if arg_present(wpe_wce)          then wpe_wce          = (*self.info).wpe_wce
    if arg_present(mi_me)            then mi_me            = (*self.info).mi_me
    if arg_present(vthi)             then vthi             = (*self.info).vthi
    if arg_present(vthe)             then vthe             = (*self.info).vthe
end


;+
;   Create a file name.
;-
function MrFD::Make_Filename, name, grid, probe
    compile_opt strictarr
    on_error, 2

    ;Name is lowercase
    _name = strlowcase(name)

    ;Create the file name
    strGrid  = string(FORMAT='(%"grid.%i")', grid)
    strProbe = string(FORMAT='(%"%sk.%i")', _name, probe)
    filename = filepath(strProbe, SUBDIRECTORY=[_name, strGrid], ROOT_DIR=self.directory)
    
    return, filename
end


;+
;   Read the binary info file.
;-
function MrFD::ReadInfo_Binary, filename
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif
    
    ;If the file does not exist or none was given, open a dialog
    if n_elements(filename) eq 0 then filename = ''
    if file_test(filename) eq 0 then begin
        filename = dialog_pickfile(/READ, TITLE='Choose binary info file.')
        if filename eq '' then return, !Null
    endif
    
    ;Allocate Memory
    topology_x = 1.0D
    topology_y = 1.0D
    topology_z = 1.0D
    Lx         = 1.0D
    Ly         = 1.0D
    Lz         = 1.0D
    nx         = 1.0D
    ny         = 1.0D
    nz         = 1.0D
    sim_dt     = 1.0D
    mime       = 1.0D
    wpewce     = 1.0D
    vthe       = 1.0D
    vthi       = 1.0D
    ey_save_interval = 1L

    ;Open the file
    openr, lun, filename, /GET_LUN

    ;Read quantities
    readu, lun, topology_x, topology_y, topology_z, Lx, Ly, Lz, nx, ny, nz
    readu, lun, sim_dt, mime, wpewce
    readu, lun, ey_save_interval

    ;Close the file and free the LUN
    free_lun, lun

    ;Store the information in a structure.
    info = { topology_x       : topology_x, $
             topology_y       : topology_y, $
             topology_z       : topology_z, $
             Lx               : Lx, $
             Ly               : Ly, $
             Lz               : Lz, $
             nx               : nx, $
             ny               : ny, $
             nz               : nz, $
             sim_dt           : sim_dt, $
             mi_me            : mime, $
             wpe_wce          : wpewce, $
             vthe             : vthe, $
             vthi             : vthi, $
             ey_save_interval : ey_save_interval $
           }
    
    return, info
end


;+
;   Read the coordinates of each probe.
;-
function MrFD::ReadCoords, filename, $
XGRID=xgrid, $
YGRID=ygrid, $
ZGRID=zgrid
    compile_opt idl2
    on_error, 2
    
    ;If the file does not exist or none was given, open a dialog
    if n_elements(filename) eq 0 then filename = ''
    if file_test(filename) eq 0 then begin
        filename = dialog_pickfile(/READ, TITLE='Choose binary info file.')
        if filename eq '' then return, !Null
    endif

    ;Read the file
    data = MrRead_Ascii(filename, $
                        COLUMN_NAMES = ['COORDS', 'COORDS', 'COORDS'], $
                        COLUMN_TYPES = ['INT', 'INT', 'INT'], $
                        GROUPS       = [0, 0, 0])

    ;Store coordinates as [x,y,z] instead of [y,x,z]
    coords = data.coords[[1,0,2],*]
    data   = !Null
    
    ;Pick out the grid
    xgrid = uniq(coords[sort(coords[*,0]),0])
    ygrid = uniq(coords[sort(coords[*,1]),1])
    zgrid = uniq(coords[sort(coords[*,2]),2])

    ;Return the 
    return, coords
end


;+
;   Read data from a particular probe.
;-
function MrFD::ReadProbe, name, grid, probe, index, $
TIME=time
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(lun) gt 0 then free_lun, lun
        void = cgErrorMSG()
        return, !Null
    endif

    ;Defaults
    time = keyword_set(time)

    ;Create & test the file name
    ;   - Open a dialog if file does not exist.
    filename = self -> Make_Filename(name, grid, probe)
    if file_test(filename) eq 0 then begin
        message, 'Choose data file. File not found: "' + filename + '".', /INFORMATIONAL
        filename = dialog_pickfile(/READ, TITLE='Choose FD file.')
        if filename eq '' then return, !Null
    endif
    
    ;Info about the file
    self -> GetInfo, NY=ny
    nTimes = self.nTimes
    
    ;Open file
    openr, lun, filename, /GET_LUN
    
;---------------------------------------------------------------------
; Time at Fixed Y ////////////////////////////////////////////////////
;---------------------------------------------------------------------
    if time then begin
        data = dblarr(nTimes)
        temp = 0D
    
        ;Step through each time
        for it = 0ULL, nTimes - 1 do begin
            ;Jump ahead to the correct point along y
            offset = 8ULL * (index + it*ny)
            point_lun, lun, offset
            
            ;Read the time point
            readu, lun, temp
            
            ;Sore the data
            data[it] = temp
        endfor
;---------------------------------------------------------------------
; Y at Fixed Time ////////////////////////////////////////////////////
;---------------------------------------------------------------------
    endif else begin
        data = dblarr(ny)
    
        ;Jump ahead in time to the correct y
        offset = 8ULL * index * ny
        point_lun, lun, offset
        
        ;Read the data along y
        readu, lun, data
    endelse
    
    ;Close file and free LUN
    free_lun, lun

    ;Return the data    
    return, data
end


;+
;   Initialize the FD object.
;-
function MrFD::Init, $
COORDS_FILE=coords_file, $
DIRECTORY=directory, $
INFO_FILE=info_file
    compile_opt idl2
    
    catch, the_error
    if the_error ne 0 then begin
        catch, /CANCEL
        if n_elements(pwd) gt 0 then cd, pwd
        void = cgErrorMSG()
        return, 0
    endif
    
    ;Defaults
    if n_elements(directory) eq 0 then cd, CURRENT=directory
    
    ;Get directory one up from the data directory
    cd, CURRENT=pwd
    cd, directory
    cd, '..'
    cd, CURRENT=dir_upOne
    cd, pwd
    if n_elements(info_file)   eq 0 then info_file   = filepath('info.bin',        ROOT_DIR=dir_upOne)
    if n_elements(coords_file) eq 0 then coords_file = filepath('coordinates.dat', ROOT_DIR=dir_upOne)
    
    ;Allocate Heap
    self.coords = ptr_new(/ALLOCATE_HEAP)
    self.info   = ptr_new(/ALLOCATE_HEAP)
    
    ;Set Properties
    self.directory = directory
    *self.info   = self -> ReadInfo_Binary(info_file)
    *self.coords = self -> ReadCoords(coords_file)
    
    ;Print information about the simulation
    self -> GetInfo, NX=nx, NY=ny, NZ=nz, MI_ME=mi_me
    nTimes = ny*8
    
    ;Number of time slices
    ;   - Get a sample file name.
    ;   - Numbers are doubles (8 bytes)
    ;   - T0 at all y, T1 at all y, etc.
    ey_file = self -> Make_Filename('ey', 0, 0)
    ey_info = file_info(ey_file)
    nTimes  = ey_info.size / (8D*ny)
    self.nTimes = nTimes
    
    message, '', /INFORMATIONAL
    print, FORMAT='(%"   Sim Size: %ix%ix%i")', nx, ny, nz
    print, FORMAT='(%"   mi/me:    %0.1f")', mi_me
    print, FORMAT='(%"   nTimes:   %i")', nTimes
    
    return, 1
end


;+
;
;-
pro MrFD__Define, class
    class = { MrFD, $
              coords:    ptr_new(), $
              directory: '', $
              info_file: '', $
              info:      ptr_new(), $
              nTimes:    0L $
            }
end