subroutine gridMercator(unit,grid,LLlat,LLlong,projection,proj_args,    &
     myNc,myNr,expandToLL,interpToLL)

!!$/*
!!$ * The fst2v5d package converts RPN Standard Files (FST) to Vis5D
!!$ * input files.  Copyright (C) 2003 Ron McTaggart-Cowan and John Gyakum,
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
!!$ * You should have received a copy of the GNU General Public License
!!$ * along with this program; if not, write to the Free Software
!!$ * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!!$ *
!!$ */

  implicit none

  ! This subroutine handles the mapping of rotated Mercator input grids.

  ! Program Name: gridMercator
  ! Version: 1.0
  ! Programmer: Ron McTaggart-Cowan
  ! Date Modified: 14 July 2003
  ! Modified by:

  real :: resBoost
  parameter (resBoost=1.)

  ! Universal constants
  real :: pi
  parameter (pi=3.141592654)

  integer MAXCOLUMNS, MAXROWS
  parameter (MAXROWS=400)
  parameter (MAXCOLUMNS=400)

  ! Input variables
  integer, intent(in) :: grid,unit

  ! Input/output variables
  logical, intent(inout) :: expandToLL,interpToLL

  ! Output variables
  integer, intent(out) :: projection,myNc,myNr
  real, dimension(100), intent(out) :: proj_args
  real, dimension(:), pointer :: LLlat,LLlong

  ! Internal variables
  integer :: err,key,ig1,ig2,ig3,junk,rig1,rig2,rig3,rig4,i,nc,nr
  real :: long1,lat1,long2,lat2,latInc,longInc
  real :: minRot=0.5
  real, dimension(:), allocatable :: lat,long
  real, dimension(:,:), allocatable :: fullLat,fullLong
  character(len=8) :: cjunk

  ! External functions
  integer, external :: fstinf,fstprm,gdll,gdllfxy,ezsetval,fstluk

  ! Get grid parameters
  key = fstinf(unit,nc,junk,junk,-1,' ',-1,-1,-1,' ','>>')
  err = fstprm(key,junk,junk,junk,junk,junk,junk,junk,junk,             &
       ig1,ig2,ig3,cjunk,cjunk,cjunk,cjunk,rig1,rig2,                   &
       rig3,rig4,junk,junk,junk,junk,junk,junk,junk)
  allocate(long(nc))
  err = fstluk(long,key,junk,junk,junk)
  key = fstinf(unit,junk,nr,junk,-1,' ',ig1,ig2,ig3,' ','^^')
  allocate(lat(nr))
  err = fstluk(lat,key,junk,junk,junk)

  ! Compute map lat/longs
  allocate( fullLat(nc,nr)  )
  allocate( fullLong(nc,nr) )
  err = gdll(grid,fullLat,fullLong)

  ! Compute rotation of the grid
  long1 = atan((fullLong(nc/2,nr)-fullLong(nc/2,1)) /                   &
       (fullLat(nc/2,nr)-fullLat(nc/2,1))) * (180./pi)

  ! Unrotated Mercator - remap to lat/long
  if (abs(long1) < minRot) then  !rotation is less than minRot degrees
     if (interpToLL) then
        print*, 'INTERPOLATING MERCATOR GRID - request'
     else
        print*, 'DETECTED UNROTATED MERCATOR: projecting to ',          &
             'proportional unrotated Lat/Long...'
        interpToLL = .true.
     endif

     ! Domain corner parameters
     err = gdllfxy(grid,lat1,long1,1.,float(nr),1)

     ! Lat/Long increments
     latInc = (lat(nr)-lat(1))/(float(nr)-1.)
     longInc = (long(nc)-long(1))/(float(nc)-1.)

     ! Domain dimensions
     myNr = nr; myNc = nc

     ! Define new domain descriptors
     allocate( LLlat(myNr)  )
     allocate( LLlong(myNc) )
     do i=1,nr
        LLlat(i) = fullLat(1,1)+latInc*float(i-1)
     enddo
     do i=1,nc
        LLlong(i) = fullLong(1,1)+longInc*float(i-1)
     enddo
     
     ! Fill V5D projection vector
     projection = 1  !unrotated lat/long
     proj_args(1) = lat1
     proj_args(2) = long1
     if (proj_args(2) < 180.) proj_args(2)=-proj_args(2)
     if (proj_args(2) > 180.) proj_args(2)=360.-proj_args(2)
     proj_args(3) = latInc
     proj_args(4) = longInc

  else  !rotated mercator - convert to lat/long grid
     if (expandToLL) then
        print *, 'EXPANDING MERCATOR GRID - request'
     else
        print '(a,f6.2,a)', 'DETECTED ROTATED MERCATOR (',long1,' degrees): projecting to resized Lat/Long...'
        expandToLL = .true.
     endif
          
     ! Find domain extent
     lat1 = maxval(fullLat)
     lat2 = minval(fullLat)
     long1 = min(fullLong(1,1),fullLong(1,nr))
     long2 = max(fullLong(nc,1),fullLong(nc,nr))
     if (long2 < long1) long2=long2+360. !greenwich protection
     
     ! Domain dimensions (scale to resolution)
     myNr = abs(int(float(nr)*(lat1-lat2)*resBoost/             &
          (lat(nr)-lat(1))))
     myNr = min(myNr,MAXROWS)
     myNc = abs(int(float(nc)*(long1-long2)*resBoost/           &
          (long(1)-long(nc))))
     myNc = min(myNc,MAXCOLUMNS)

     ! Lat/Long increments
     latInc = abs(lat1-lat2)/(float(myNr)-1)
     longInc = abs(long1-long2)/(float(myNc)-1)

     ! Define new domain descriptors
     allocate( LLlat(myNr)  )
     allocate( LLlong(myNc) )
     do i=1,myNr
        LLlat(i) = lat2+latInc*(float(i-1))
     enddo
     do i=1,myNc
        LLlong(i) = long1+longInc*(float(i-1))
     enddo

     ! Fill V5D projection vector
     projection = 1 ! unrotated lat/long
     proj_args(1) = lat1
     proj_args(2) = long1
     if (proj_args(2) < 180.) proj_args(2)=-proj_args(2)
     if (proj_args(2) > 180.) proj_args(2)=360.-proj_args(2)
     proj_args(3) = latInc
     proj_args(4) = longInc

     ! Set extrapolation value for ezscint package
!print *,'ezsetval - 1'
!     err = ezsetval('EXTRAP_VALUE', 0.0)
!print *,'ezsetval - 2'

  endif

  ! End of subroutine
  return
end subroutine gridMercator
