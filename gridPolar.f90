subroutine gridPolar(unit,grid,LLlat,LLlong,projection,proj_args,       &
     myNc,myNr,expandToLL)

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

  ! This subroutine handles the mapping of Polar Steregraphic input grids.

  ! Program Name: gridPolar
  ! Version: 1.0
  ! Programmer: Ron McTaggart-Cowan
  ! Date Modified: 14 July 2003
  ! Modified by:
  
  real :: resBoost
  parameter (resBoost=1.)

  ! Universal constants
  real :: pi,dtr,rtd,rEarth
  parameter (pi=3.141592654)
  parameter (dtr=pi/180.)
  parameter (rtd=180./pi)
  parameter (rEarth=6371.)

  integer MAXCOLUMNS, MAXROWS
  parameter (MAXROWS=400)
  parameter (MAXCOLUMNS=400)

  ! Input variables
  integer, intent(in) :: grid,unit

  ! Input/output variables
  logical, intent(inout) :: expandToLL

  ! Output variables
  integer, intent(out) :: projection,myNc,myNr
  real, dimension(100), intent(out) :: proj_args
  real, dimension(:), pointer :: LLlat,LLlong

  ! Internal variables
  integer :: err,key,ig1,ig2,ig3,junk,rig1,rig2,rig3,rig4,i,nc,nr
  real :: rjunk,long1,lat1,long2,lat2,centIndX,centIndY,d60,dgrw,       &
       latInc,longInc
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

  ! Get map parameters
  call cigaxg('N',rjunk,rjunk,d60,dgrw,rig1,rig2,rig3,rig4)

  ! Compute map lat/longs
  allocate( fullLat(nc,nr)  )
  allocate( fullLong(nc,nr) )
  err = gdll(grid,fullLat,fullLong)

  ! Input map coordinates and rotation
  centIndX = float(nc)/2.
  centIndY = float(nr)/2.
  d60 = abs(long(2)-long(1))
  long1 = fullLong(int(centIndX),nr)-fullLong(int(centIndX),1)
  if (mod(centIndX,1.) > .25) long1 =                            &
              (fullLong(int(centIndX),nr)+                              &
              fullLong(int(centIndX+1.),nr))/2. -                       &
              (fullLong(int(centIndX),1)+                               &
              fullLong(int(centIndX+1.),1))/2.
  
  ! Unrotated Polar Stereographic
  if (abs(long1) < minRot .and. .not. expandToLL) then  !rotation is less than minRot degrees
     print*, 'DETECTED UNROTATED STEREOGRAPHIC PROJECTION'
     
     ! Central lat/long
     err = gdllfxy(grid,lat1,long1,centIndX,centIndY,1)

     ! Initialize output variables
     myNr = nr; myNc = nc
     allocate( LLlat(myNr)  )
     allocate( LLlong(myNc) )
     LLlat = 0.
     LLlong = 0.
     
     ! Fill V5D projection vector
     projection = 3  !polar stereographic
     proj_args(1) = lat1
     proj_args(2) = long1
     if (proj_args(2) < 180.) proj_args(2)=-proj_args(2)
     if (proj_args(2) > 180.) proj_args(2)=360.-proj_args(2)
     proj_args(3) = centIndY
     proj_args(4) = centIndX
     proj_args(5) = d60 /                                               &
                 ((1.+sin(60.*dtr))/(1.+sin(lat1*dtr)))

  else  !rotated polar - convert to lat/long grid
     if (expandToLL) then
        print*, 'EXPANDING POLAR STEREOGRPHIC GRID - request'
     else
        print*, 'DETECTED ROTATED POLAR STEREOGRAPHIC (',long1,         &
             'degrees): projecting to resized Lat/Long'
        expandToLL = .true.
     endif
     
     ! Find domain extent
     lat1 = maxval(fullLat)
     lat2 = minval(fullLat)
     long1 = minval(fullLong(1,:))
     long2 = maxval(fullLong(nc,:))

     ! Domain dimensions (scale to resolution)
     myNr = nr *                                                        &
          abs(rEarth*(lat1-lat2)/(d60*nr*rtd))*resBoost
     myNr = min(myNr,MAXROWS)
     myNc = nc *                                                        &
          abs(rEarth*(lat1-lat2)/(d60*nr*rtd))*resBoost
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
     projection = 1  !unrotated lat/long
     proj_args(1) = lat1
     proj_args(2) = long1
     if (proj_args(2) < 180.) proj_args(2)=-proj_args(2)
     if (proj_args(2) > 180.) proj_args(2)=360.-proj_args(2)
     proj_args(3) = latInc
     proj_args(4) = longInc

     ! Set extrapolation value for ezscint package
     err = ezsetval('EXTRAP_VALUE',0.)

  endif

  ! End of subroutine
  return
end subroutine gridPolar
