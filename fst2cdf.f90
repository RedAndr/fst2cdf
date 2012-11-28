  program fst2cdf

  use Gridder

!!$/*
!!$ * The fst2cdf package converts RPN Standard Files (FST) to NetCDF files.  
!!$ * Copyright (C) 2003 Ron McTaggart-Cowan and John Gyakum,
!!$ * Copyright (C) 2006 Andrew Ryzhkov and Michel Bourqui,
!!$ * at McGill University.
!!$ *
!!$ * This library is free software; you can redistribute it and/or
!!$ * modify it under the terms of the GNU Lesser General Public
!!$ * License as published by the Free Software Foundation; either
!!$ * version 2.1 of the License, or (at your option) any later version.
!!$ *
!!$ * This library is distributed in the hope that it will be useful,
!!$ * but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$ * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!$ * Lesser General Public License for more details.
!!$ *
!!$ * You should have received a copy of the GNU Lesser General Public
!!$ * License along with this library; if not, write to the Free Software
!!$ * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!!$ *
!!$ */

!  This program converts RPN Standard Files (FST) to NetCDF files
!  The primary jobs of this program are to sort through the input file for 
!  projection, variable, level, and time information which is then packaged 
!  to NetCDF. There are a few limitations (esp for projections):

!  Projections:
!    Unrotated Polar Stereographic:  The central meridian of the grid
                            !  points due north (i.e. towards the 
                            !  pole).  This type of grid is passed
                            !  directly to the Vis5D output file.
!    Rotated Polar Stereographic:  The central meridian of the grid is
                            !  NOT oriented north/south and therefore
                            !  does not point directly at the pole.
                            !  In this case, the data is interpolated
                            !  onto the minimally-sized lat/long grid
                            !  that will fit the entire dataset.  The
                            !  resolution of the lat/long grid is
                            !  approx "resBoost" (set in v5df.h) * the
                            !  resolution of the original grid to take
                            !  the subdomain rotation into account.
!    Unrotated Mercator:  The bounds of this projection are lat/long lines,
                            !  and so a simple interpolation to the
                            !  perfect-fit lat/long grid is performed.
                            !  Note that there is no resolution boost 
                            !  in this operation.  This is the cleanest
                            !  type of grid to view.
!    Rotated Mercator:  The minimally-size lat/long box that will fit the
                            !  full dataset is determined, and, as for
                            !  the rotated stereographic grid, the
                            !  resolution of the output grid is increased
                            !  by the "resBoost" factor.  The bounds of
                            !  the map displayed by Vis5D are of course
                            !  extended beyond those of the data.

!  A note about the domain extensions:  all values in the Vis5D domain which
!  extend beyond those of the input domain are set to the minimum value
!  on the horizontal domain of the variable, less 5% of the variability
!  on the level.  For more information of the interpolation of data
!  in this procedure, see the (sparse) package documentation for "ezscint",
!  a very nice set of interpolation procedures written by Yves Chartier.

!  Program Name:  fst2cdf
!  Version: 1.0
!  Programmer: Vis5D+ template file
!  Date: 11 July 2003
!  Modified by: Ron McTaggart-Cowan
!  Date: 20 Nov 2006
!  Modified by: Andrew Ryzhkov
!  Date: 15 Nov 2007
!  Modified by: AR
!  Date: 10 Jan 2010
!  Modified by: AR

  implicit none

  include  "netcdf.inc"

  ! Missing values
  real    MISSING
  integer IMISSING
  parameter (MISSING =-999.99)
  parameter (IMISSING=-999   )

  ! grid limits
  integer MAXVARS, MAXTIMES, MAXLEVELS

  parameter (MAXVARS=256)
  parameter (MAXTIMES=1024)
  parameter (MAXLEVELS=64)

  ! Local variables
  integer :: m!,n
  character (len=100) inname, outname
  integer :: fileUnit=20
  integer :: maxnl,nFields,Pig1,Pig2,Pig3,maxLevelCount,ip1,LLGrid
  integer :: Pig4,Rig1,Rig2,Rig3,Rig4,Lig1,Lig2,Lig3,Lig4
  integer :: i,j,k,err,junk,dMark,varFields,grid
  integer :: key,fieldCount,LevelCount,ip2,ip3,date
  integer :: myNr,myNc,deet,npas,yyyy,mm,dd
  integer :: hh,mn,ss,hits,nMax
  real    :: tMark,cMark,fMark
  character(len=8) ::  cjunk
  integer, dimension(2) :: fullDate
  integer, dimension(:), allocatable :: fieldList,newLevel,allLevels,keyList
  character(len=1)   :: gridType,varType
  character(len=4)   :: name
  character(len=4), dimension(:), allocatable :: newField
  character(len=256) :: TimeStr
  character(len=16)  :: etiket
  logical :: foundField,foundLevel,fieldOK,interpToLL,expandToLL,expandToPS=.false.,gotIt
  real :: xg1,xg2,xg3,xg4
  real, dimension(:), pointer :: LLlat,LLlong
  real, dimension(:,:), allocatable :: myField,readField
  real (kind=4), dimension(:,:,:), allocatable :: G
  real :: ptop,pref,rcoef
  integer :: ierr, cdfid, zdim(2), latid, lonid, levid, timid, ig1,ig2,ig3,ig4
  integer :: start(4), count(4), latvar, lonvar, kind 
  integer :: levvar, idim, timvar, lv0id, lv0var, aklvar, bklvar
  integer, dimension(:), allocatable :: var_id
  real    :: Eta(MAXLEVELS), Ak(MAXLEVELS), Bk(MAXLEVELS), plev(MAXLEVELS)

  logical write_lat, write_lon, append, info
  data write_lat /.false./, write_lon /.false./, append /.false./, info /.false./

!     External functions
  integer, external :: fnom,fstouv,fstfrm,fstinl,fstnbr,fstprm,     &
       fstinf,fstluk,fstlir,gdllfxy,gdxyfll,ezqkdef,gdll,           & 
       ezgxprm,ezsint,ezdefset,ezsetopt,ezgdef,fstecr,ezsetval,     &
       newdate, read_decode_hyb

  integer nr, nc, nl(MAXVARS)
  integer NumTimes, NumVars, NumTimes0
  data NumTimes0 / 0 /
  character*16 varname(MAXVARS)
  real Time, Times(MAXTIMES)
  integer projection
  real proj_args(128)
  real vert_args(MAXLEVELS)
  integer ios

  ! initialize the variables to missing values
  data nr,nc / IMISSING, IMISSING /
  data (nl(i), i=1,MAXVARS) / MAXVARS*IMISSING /
  data NumTimes,numvars / IMISSING, IMISSING /
  data (varname(i), i=1,MAXVARS) / MAXVARS*"          " /
  data (Times(i), i=1,MAXTIMES) / MAXTIMES*IMISSING /
  data projection / IMISSING /
  data (proj_args(i), i=1,128) / 128*MISSING /
  data (vert_args(i), i=1,MAXLEVELS) / MAXLEVELS*MISSING /

  ! level name and units in the FST format
!  character*64 LevelNames (0:6)
  character*16  LevelUnits (0:6)
!  data LevelNames / "height", "sigma", "pressure", "arbitrary code", &
!                    "height with respect to ground level",           &
!                    "hybrid coordinates", "theta" /
  data LevelUnits / "m", "sigma", "mb", "??", "M", "hy", "th" /

  character*256 CDF_Title
  data CDF_Title  / 'NetCDF file created by the FST2CDF program, Andrew Ryzhkov <http://RedAndr.ca>' /

  character*256 Usage_msg
  data Usage_msg / "Usage: fst2cdf -i input_file.fst -o output_file.cdf [options]\n"&
  "options:\n"&
  "-la and -lo - output latitude and longitude variables" /		!C

  ! arguments variables
  integer iargc
  character (len=256) arg


  include  "FST_var.h" 


  ! get command line arguments
  if (iargc() == 0) then
     print '(a)',Usage_msg
     call exit(1)
  endif

  i=1
  do while (iargc().ge.i)
       call getarg(i,arg)
       i=i+1
       select case (arg)
       case ('-i')
           call getarg(i,inname)
           i=i+1
       case ('-o')
           call getarg(i,outname)
           i=i+1
       case ('-la')
           write_lat = .true.
       case ('-lo')
           write_lon = .true.
       case ('-h')
           print *,Usage_msg
           stop
       case ('-f')
           info = .true.
       end select
  enddo

  if ( len_trim(inname) == 0 ) then
    print *, 'No input file specified with "-i" argument'
    stop
  end if  

  if ( len_trim(outname) == 0 ) then
    print *, 'No output file specified with "-o" argument'
    stop
  end if  

  print '(2a)',"Input file : ", inname
  print '(2a)',"Output file: ", outname
!  print *, 'write_lat, write_lon: ', write_lat, write_lon

  ! Open data file
  err = fnom(fileUnit, trim(inname), 'STD+RND',0)
  if (err < 0) then
     print *, 'ERROR: cannot open file <',trim(inname),'>'
     stop
  endif

  err = fstouv(fileUnit,'RND')
  if (err < 0) then
     print *, 'ERROR: cannot open <',trim(inname),'> in random mode'
     stop
  endif

  ! Read hybryd level parameters
  err = read_decode_hyb(fileUnit,'HY',-1,-1,' ',-1,ptop,pref,rcoef) !date
!  print *,err
  print '(a,f8.4,f8.1,f8.4)', 'Ptop, Pref, Rcoef=', Ptop, Pref, Rcoef

  ! Get grid dimensions and specifications
  key = fstinf(fileUnit,nc,junk,junk,-1,' ',-1,-1,-1,' ','>>')
  err = fstprm(key,date,deet,npas,junk,junk,junk,                   &
       junk,junk,Pig1,Pig2,Pig3,cjunk,name,etiket,gridType,Rig1,Rig2, &
       Rig3,Rig4,junk,junk,junk,junk,junk,junk,junk)
  key = fstinf(fileUnit,junk,nr,junk,-1,' ',Pig1,Pig2,Pig3,' ','^^')
  print '(2a)','Grid type:',gridType
  print '(a,i3,a,i3)','Grid size:',nc,'x',nr

  ! Define input grid
  err = ezsetopt('VERBOSE','NON')
  grid = ezqkdef(nc,nr,'Z',Pig1,Pig2,Pig3,Pig4,fileUnit)

  ! Initialize grid requests
  expandToPS = .false.; expandToLL = .false.; interpToLL = .false.

  select case (gridType)
    case ('N')                ! polar stereographic input
       call gridPolar(fileUnit,grid,LLlat,LLlong,projection,          &
            proj_args,myNc,myNr,expandToLL)
    
    case ('E')                ! rotated mercator input
       call gridMercator(fileUnit,grid,LLlat,LLlong,projection,       &
            proj_args,myNc,myNr,expandToLL,interpToLL)

    case DEFAULT
       print*, "WARNING: unsupported grid type ",gridType
       print*, "attempting to continue ..."
       projection = 0
       proj_args(1) =   90.
       proj_args(2) = -180.
       proj_args(3) = 1.
       proj_args(4) = 1.
  end select

  if (info) stop							! print info and exit

  ! Get date and field information
  nFields = fstnbr(fileUnit)
  allocate( newField(nFields) )
  allocate( newLevel(nFields) )
  allocate( fieldList(nFields) )
  err = fstinl(fileUnit,junk,junk,junk,-1,' ',-1,-1,-1,' ',' ',fieldList,nFields,nFields+1)

  FieldCount=0; NumTimes=0

  do i=1,nFields
     err = fstprm(fieldList(i), date,deet,npas,junk,junk,junk,              &
          junk,junk,junk,ip2,ip3,varType,name,etiket,cjunk,junk,junk,        &
          junk,junk,junk,junk,junk,junk,junk,junk,junk)

     if (i == 1) then 
       err  = newdate(date, fullDate(1), fullDate(2), -3)		! origin time
       yyyy = fullDate(1)/10000
       mm   = mod(fullDate(1)/100,100)
       dd   = mod(fullDate(1), 100)
       hh   = fullDate(2)/1000000
       mn   = mod(fullDate(2)/10000,100)
       ss   = mod(fullDate(2),100)
     endif

     Time = deet*npas/3600.						! (deet*npas+1800.)/3600.

     if ( (i==3) .or. (i>3 .and. Time>Times(NumTimes)) ) then
       NumTimes = NumTimes+1
       Times(NumTimes) = Time
     endif

     ! Exit on special names
     fieldOK = .true.
     if (name == '>>' .or. name == '^^' .or. name == 'HT')          &
          fieldOK = .false.

     ! Field names
     foundField=.false.
     j=1
     do while (.not.foundField .and. j <= fieldCount .and. fieldOK)
        if (name == newField(j)) foundField=.true.
        j=j+1
     enddo
     if (.not.foundField .and. fieldOK) then
        newField(j) = name
        fieldCount = j
     endif

  enddo

!     Enter values for field names and time list
  numvars = fieldCount
  varname = newField(1:fieldCount)

print *,'Variables      : ',NumVars
print *,'Variables names: ',VarName(1:NumVars)
print *,'Number of times: ',NumTimes
print *,'Times array    : ',Times(1:NumTimes)

  ! Grid levels for each variable
  allocate( allLevels(nFields) ); maxLevelCount=0

  do i=1,fieldCount
     
     err = fstinl(fileUnit,junk,junk,junk,dMark,' ',-1,tMark,-1,' ',&
          newField(i), fieldList,varFields,nFields+1)
     if (varFields == 0)                                            &
          err = fstinl(fileUnit,junk,junk,junk,dMark,' ',-1,-1,-1,  &
          ' ',newField(i), fieldList,varFields,nFields+1)    
     if (varFields == 0)                                            &
          err = fstinl(fileUnit,junk,junk,junk,-1,' ',-1,-1,-1,     &
          ' ',newField(i), fieldList,varFields,nFields+1)
     if (varFields == 0) print *, "WARNING: cannot find ",newField(i)
     
print *,'field=',i,'varFields=',varFields

     levelCount = 0

     do j=1,varFields

        err = fstprm(fieldList(j), date,  junk,junk,junk,junk,junk, &
             junk,junk,ip1,junk,junk,cjunk,cjunk,cjunk,cjunk,       &
             ig1,ig2,ig3,ig4, junk,junk,junk,junk, junk,junk,junk)

        foundLevel=.false.
        k=1

        do while (.not.foundLevel .and. k <= levelCount)
           if (ip1 == newLevel(k)) foundLevel=.true.
           k=k+1
        enddo

        if (.not. foundLevel) then
           newLevel(k) = ip1
           levelCount = k
        endif

!        if ( i==1 .and. j==1) then
!          call convip( ip1, Ptop, kind, -1, "", .false.)              ! to get level name and dimension from HY record
!          pref = ig1
!          rcoef = ig2/1000.0
!print '(a,i6,f8.4,f8.1,f8.4)', 'ip1, Ptop, Pref, Rcoef=', ip1, Ptop, Pref, Rcoef
!        endif
     enddo

     if (levelCount > maxLevelCount) then
        allLevels(1:levelCount) = newLevel(1:levelCount)
        maxLevelCount = levelCount
     endif
     
     ! Enter number of levels
     nl(i) = levelCount
     
  enddo

  ! Enter values for levels
  vert_args = allLevels(1:maxLevelCount)

print *,'vert_args=',vert_args(1:maxLevelCount)
print *,'nl=',nl(1:numvars)

  deallocate(fieldList)

  ! Calculate number of levels      
  maxnl = nl(1)
  do i=2,numvars
     if (nl(i) > maxnl) then
        maxnl = nl(i)
     endif
  enddo

  ! =================================== !
  !     Create/open the NetCDF file.    !
  ! =================================== !

  open ( unit=10, file=outname, status='old', iostat=ios)
  if ( ios == 0 ) then
    close ( unit=10 )
    ierr = nf_open ( outname, NF_WRITE, cdfid )
    call check_err(ierr,'Open NetCDF file ')
    print '(3a)','NetCDF file <',trim(outname), '> was opened'          ! the NetCDF file is exist
    append = .true.							! just append data
  else
     ierr = nf_create ( outname, NF_NOCLOBBER, cdfid )
    call check_err(ierr,'Create NetCDF file ')
    print '(3a)','NetCDF file <',trim(outname), '> was created'         ! create new NetCDF file
  endif

  if ( .not. append ) then

    ierr = nf_put_att_text ( cdfid, NF_GLOBAL, 'title', LEN_TRIM(CDF_Title), CDF_Title )
    call check_err(ierr,'Define title ')

    ! Define dimensions of the data in the NetCDF file: 
    ! latitude and longitude, level, level0 and time
    ierr = nf_def_dim ( cdfid, 'lon'   , myNc , lonid )
    call check_err(ierr,'Define longitude dimension ')

    ierr = nf_def_dim ( cdfid, 'lat'   , myNr , latid )
    call check_err(ierr,'Define latitude dimension ')

    ierr = nf_def_dim ( cdfid, 'level' , maxnl, levid )
    call check_err(ierr,'Define level dimension ')

    ierr = nf_def_dim ( cdfid, 'level0', 1    , lv0id )
    call check_err(ierr,'Define level0 dimension ')

    ierr = nf_def_dim ( cdfid, 'time'  , NF_UNLIMITED, timid )
    call check_err(ierr,'Define time dimension ')

    ! create variables in the NetCDF

    ! define longitude variable
    ierr = nf_def_var ( cdfid, 'lon', NF_FLOAT, 1, lonid, lonvar )
    call check_err(ierr,'Define longitude variable ')
    ierr = nf_put_att_text ( cdfid, lonvar, 'long_name', 9, 'Longitude' )
    call check_err(ierr,'Define long_name ')
    ierr = nf_put_att_text ( cdfid, lonvar, 'units'    , 7, 'degrees'   )
    call check_err(ierr,'Define lon units ')

    ! define latitude variable
    ierr = nf_def_var ( cdfid, 'lat', NF_FLOAT, 1, latid, latvar )
    call check_err(ierr,'Define latitude variable ')
    ierr = nf_put_att_text ( cdfid, latvar, 'long_name', 8, 'Latitude'  )
    call check_err(ierr,'Define long_name ')
    ierr = nf_put_att_text ( cdfid, latvar, 'units'    , 7, 'degrees'   )
    call check_err(ierr,'Define lat units ')

    ! define level variable
    ierr = nf_def_var ( cdfid, "level", NF_FLOAT, 1, levid, levvar )
    call check_err(ierr,'Define level variable ')
    ierr = nf_put_att_text ( cdfid, levvar, 'long_name',14, 'Pressure level' )
    call check_err(ierr,'Define long_name ')
    ierr = nf_put_att_text ( cdfid, levvar, 'units'    , 4, 'mbar'    )		! hPa
    call check_err(ierr,'Define level units ')

    ! define level0 variable
    ierr = nf_def_var ( cdfid, "level0", NF_FLOAT, 1, lv0id, lv0var )
    call check_err(ierr,'Define level0 variable ')

    ! define A(k) and B(k) variables
    ierr = nf_def_var ( cdfid, "aklay", NF_FLOAT, 1, levid, aklvar )
    call check_err(ierr,'Define A(k) variable ')
    ierr = nf_put_att_text ( cdfid, aklvar, 'long_name', 4, 'A(k)' )
    call check_err(ierr,'Define long_name ')
    ierr = nf_def_var ( cdfid, "bklay", NF_FLOAT, 1, levid, bklvar )
    call check_err(ierr,'Define B(k) variable ')
    ierr = nf_put_att_text ( cdfid, bklvar, 'long_name', 4, 'B(k)' )
    call check_err(ierr,'Define long_name ')

    ! define time variable
    ierr = nf_def_var ( cdfid, "time", NF_FLOAT, 1, timid, timvar )
    call check_err(ierr,'Define time variable ')
    TimeStr = "time"
    ierr = nf_put_att_text ( cdfid, timvar, 'long_name', LEN_TRIM(TimeStr), TimeStr )
    call check_err(ierr,'Define long_name ')
    write (TimeStr,'("hours since ",i4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)')  &
           yyyy, mm, dd, hh, mn, ss
    ierr = nf_put_att_text ( cdfid, timvar, 'units'    , LEN_TRIM(TimeStr), TimeStr )
    call check_err(ierr,'Define time units ')

  endif

  ! define other variables
  allocate( var_id(numvars) )

  if ( .not. append ) then

    idim = 4						! 4 dimensions: lat, lon, lev, tim
    zdim(1) = lonid;      zdim(2) = latid;      zdim(3) = levid;      zdim(4) = timid

    do j=1, NumVars

      if (trim(varname(j))=='LA' .and. .not. write_lat) cycle
      if (trim(varname(j))=='LO' .and. .not. write_lon) cycle
      if (trim(varname(j))=='HY') cycle			! hybrid levels defenition

      if ( nl(j) == 1 ) then                          	! only one level => no level dependant variable
        zdim(3) = lv0id;
      else
        zdim(3) = levid;
      endif
    
!      ierr = nf_def_var ( cdfid, trim(varname(j)), NF_FLOAT, idim, zdim, var_id(j) )
!      call check_err(ierr,'Define variables ')

print *,"Test variable: ",trim(varname(j))
      do i=1, FST_VARS_COUNT
         if( trim(vars(i)%shortname) .eq. trim(varname(j)) ) then
           ierr = nf_def_var ( cdfid, trim(vars(i)%cdfname), NF_FLOAT, idim, zdim, var_id(j) )
           call check_err(ierr,'Define variables ')
           ierr = nf_put_att_text ( cdfid, var_id(j), 'long_name'    , LEN_TRIM(vars(i)%longname), vars(i)%longname )
           call check_err(ierr,'Define long_name ')
           ierr = nf_put_att_text ( cdfid, var_id(j), 'units'        , LEN_TRIM(vars(i)%units)   , vars(i)%units    )
           call check_err(ierr,'Define units ')
           ierr = nf_put_att_real ( cdfid, var_id(j), 'missing_data' , NF_FLOAT, 1, MISSING )
           call check_err(ierr,'Define missing data ')
           ierr = nf_put_att_real ( cdfid, var_id(j), 'missing_value', NF_FLOAT, 1, MISSING )
           call check_err(ierr,'Define missing value ')
print *,"NetCDF variable: ",TRIM(vars(i)%cdfname)," was created."
         endif
      enddo

    enddo

  else

    ! get number of times
    ierr = nf_inq_dimid ( cdfid, 'time', timid ) 
    call check_err(ierr,'Get time dim ID ')
print *, "timid=", timid
    ierr = nf_inq_dim ( cdfid, timid, name, NumTimes0 ) 
    call check_err(ierr,'Get time len ')
print *, "NumTimes0=", NumTimes0
    ierr = nf_inq_varid ( cdfid, 'time', timvar ) 
print *, "timvar=", timvar

    ! get variables IDs
    do j=1, NumVars

      if (trim(varname(j))=='LA' .and. .not. write_lat) cycle
      if (trim(varname(j))=='LO' .and. .not. write_lon) cycle
      if (trim(varname(j))=='HY') cycle				! hybrid levels defenition

      do i=1, FST_VARS_COUNT
        if( TRIM(vars(i)%shortname) .eq. TRIM(varname(j)) ) then
          ierr = nf_inq_varid ( cdfid, TRIM(vars(i)%cdfname), var_id(j) )
          call check_err(ierr,'Get variable ID ')
print *, trim(varname(j)), " => ", TRIM(vars(i)%cdfname), " #", var_id(j)
        endif
      enddo
    enddo

  endif

  if ( .not. append ) then
! takes an open netCDF file out of define mode
! to write variables
    ierr = nf_enddef ( cdfid  )
    call check_err(ierr,'End define mode ')
  endif

! read and convert data
  allocate( readField(nc,nr)   )
  allocate( myField(myNc,myNr) )

  nMax = fstnbr(fileUnit)
  allocate( keyList(nMax)      )
  if (interpToLL .or. expandToLL .or. expandToPS) then

     if (interpToLL) err = ezsetopt('EXTRAP_DEGREE','LINEAR')
     if (expandToLL) err = ezsetopt('EXTRAP_DEGREE','MINIMUM')
     if (expandToPS) err = ezsetopt('EXTRAP_DEGREE','MINIMUM')
     if (interpToLL .or. expandToLL) then
        xg1=0.; xg2=0.; xg3=1.; xg4=1.; gridType='L'
     elseif (expandToPS) then
        xg1=0.; xg2=0.; xg3=1000.; xg4=0.; gridType='N'
     endif
     call cxgaig(gridType,Lig1,Lig2,Lig3,Lig4,xg1,xg2,xg3,xg4)
     LLgrid = ezgdef( myNc,myNr, 'Z',gridType,Lig1,Lig2,Lig3,Lig4, LLlong,LLlat )


! write latitide and longitute variables
     print *,'Latitides:' ;  print '(8f8.2)',LLlat
     print *,'Longitutes:';  print '(8f8.2)',LLlong

     if ( .not. append ) then
       ierr = nf_put_vara_real (cdfid, latvar, 1, myNr, LLlat )
       call check_err(ierr,'Write latitide variable ')
       ierr = nf_put_vara_real (cdfid, lonvar, 1, myNc, LLlong )
       call check_err(ierr,'Write longitude variable ')
     endif

     deallocate(LLlat)
     deallocate(LLlong)
     err = ezsetopt('INTERP_DEGREE','CUBIC')
  endif

  if ( .not. append ) then

    ! calculate A & B
    print '(a,i4,a)','Levels [',maxLevelCount,'] (#,Eta,Ak,Bk,P):'
    do i=1,maxLevelCount
      call CONVIP( allLevels(i), Eta(i), kind, -1, "", .false.)
      Bk(i) = Eta(i)**Rcoef
      Ak(i) = ( Pref*(Eta(i)-Bk(i)) + Ptop*(1-Eta(i)) )
      Plev(i) = Ak(i) + Bk(i) * 1015;
      print '(i6,f10.6,3a,f8.2,a,f8.4,a,f8.2)',i,Eta(i),' ',LevelUnits(kind),' ',Ak(i),' ',Bk(i),' ',Plev(i)
    enddo

    Ak  (1:maxLevelCount) = Ak  (maxLevelCount:1:-1)					! reverse levels
    Bk  (1:maxLevelCount) = Bk  (maxLevelCount:1:-1)
    Plev(1:maxLevelCount) = Plev(maxLevelCount:1:-1)

    ! write level variable
!    ierr = nf_put_vara_real (cdfid, levvar, 1, maxLevelCount, Eta )
    ierr = nf_put_vara_real (cdfid, levvar, 1, maxLevelCount, Plev )
    call check_err(ierr,'Write level variable ')

    ! write level0 variable = 1.0
    ierr = nf_put_vara_real (cdfid, lv0var, 1, 1, 1.0 )
    call check_err(ierr,'Write level0 variable ')

    ! write a(k) and b(k) variables
    ierr = nf_put_vara_real (cdfid, aklvar, 1, maxLevelCount, Ak )
    call check_err(ierr,'Write A(k) variable ')
    ierr = nf_put_vara_real (cdfid, bklvar, 1, maxLevelCount, Bk )
    call check_err(ierr,'Write B(k) variable ')

    ! write time variable
    print '(a,i4,a)','Times [',NumTimes,'] :'
    do i=1,NumTimes
      print '(i6,f12.2,2a)',i,Times(i),' ',TRIM(TimeStr)
    enddo

    ierr = nf_put_vara_real (cdfid, timvar,           1, NumTimes, Times )
    call check_err(ierr,'Write time variable ')
  
  else
    
    print *,'Appends time variable (NumTimes, Times): ',NumTimes, Times(1:NumTimes)
    ierr = nf_put_vara_real (cdfid, timvar, NumTimes0+1, NumTimes, Times )
    call check_err(ierr,'Append time variable ')

  endif

  ! write other variables

  allocate( G(myNc,myNr,maxnl) )
 
  do i=1, NumTimes

     do j=1, NumVars

        dMark = -1; tMark = Times(i)

        if (trim(varname(j))=='LA' .and. .not. write_lat) cycle
        if (trim(varname(j))=='LO' .and. .not. write_lon) cycle
        if (trim(varname(j))=='HY') cycle						! hybrid levels defenition

        do k=1, nl(j)

           err = fstinl(fileUnit,junk,junk,junk,dMark,' ',int(vert_args(k)), -1,-1,' ',trim(varname(j)),keyList,hits,nMax)

           if (hits == 0)                                           &
              err = fstinl(fileUnit,junk,junk,junk,dMark,' ',-1,-1,-1,' ',trim(varname(j)), keyList,hits,nMax)

           m=0; gotIt=.false.
           do while (.not.gotIt .and. m < hits)
              m=m+1
              err = fstprm(keyList(m), date,junk,npas,junk,junk,junk,&
                   junk,junk,ip1,junk,junk,cjunk,cjunk,cjunk,       &
                   gridType,junk,junk,junk,junk,junk,junk,junk,junk,   &		!cjunk
                   junk,junk,junk)
              cMark = real(npas)*real(deet)/3600.
              fMark = abs(cMark-tMark)
              if ( fMark < 0.001 ) then
                gotIt = .true.
!                print '(a,3f16.8)',"Found: " ,cMark,tMark,fMark
!              else
!                print '(a,3f16.8)',"Missed: ",cMark,tMark,fMark
              endif
           enddo

           if (gridType == 'G') cycle

           if (gotIt) then
!              print *,"gridType: ",gridType
              err = fstluk(readField,keyList(m), junk,junk,junk)
           else
              print *, 'WARNING: Cannot find ',trim(varname(j)), ' at step ',tMark,' ... padding'
              err = fstlir(readField,fileUnit,junk,junk,junk,  dMark,' ',-1,-1,-1,' ',trim(varname(j)))
           endif

           if (interpToLL .or. expandToLL .or. expandToPS) then
              err = ezdefset(LLgrid,grid)
              err = ezsint(myField,readField)
           else
              myField = readField
           endif

           if ( trim(varname(j))=='UU' .or. trim(varname(j))=='UV' .or. trim(varname(j))=='VV') then
             G(1:myNc,1:myNr,k) = myField * 0.51444444			! convert from knots to m/s
           else
             G(1:myNc,1:myNr,k) = myField
           endif

!           do m=1,myNr              ! nj
!              do n=1,myNc           ! ni
!                 if ( trim(varname(j))=='UU' .or.      &
!                      trim(varname(j))=='UV' .or.      &
!                      trim(varname(j))=='VV')     then
!                   G(n,m,k) = myField(n,m) * 0.51444444			! convert from knots to m/s
!                 else
!                   G(n,m,k) = myField(n,m)
!                 endif
!              enddo
!           enddo

        enddo                       ! k - levels

        ! variables start and count: 
        ! lon,         lat,           level,          time
        start(1)=1   ; start(2)=1;    start(3)=1;     start(4)=i + NumTimes0;
        count(1)=myNc; count(2)=myNr; count(3)=nl(j); count(4)=1

        G(1:myNc,1:myNr,1:nl(j)) = G(1:myNc,1:myNr,nl(j):1:-1)		! reverse levels

        print '(a,a4,a,i4,a,4i4,a,4i4)', 'Writing variable <',trim(varname(j)),   &
                            '> at time:', i, ' start:', start, ' count:',count
        ierr = nf_put_vara_real ( cdfid, var_id(j), start, count, G )
        call check_err(ierr,'Write variable: "'//trim(varname(j))//'"')

     enddo			    ! j - vars

  enddo                             ! i - times

! free memory
  deallocate(myField,readField,G,keyList,var_id)
        
! close the input file
  err = fstfrm(fileUnit)

! close the output NetCDF file
  call NCCLOS (cdfid, ierr )
  call check_err(ierr,'Close file ')

  end

  subroutine check_err(iret,msg)
    implicit none
    integer iret
    character*(*) msg
    include 'netcdf.inc'
    if ( iret /= NF_NOERR ) then
      print '("NetCDF error <",a,"> [",a,"]")', trim(msg), trim(nf_strerror(iret))
      stop
    endif
  end

