module Gridder

  !  This module contains all of the grid operations for the RPN - Vis5D
  !  converter.

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
!!$ */

  !  Interface block
  interface
     subroutine gridPolar(unit,grid,LLlat,LLlong,projection,        &
          proj_args,myNc,myNr,expandToLL)
       integer, intent(in) :: grid,unit
       integer, intent(out) :: projection,myNc,myNr
       real, dimension(100), intent(out) :: proj_args
       real, dimension(:), pointer :: LLlat,LLlong
       logical, intent(inout) :: expandToLL
     end subroutine gridPolar
     subroutine gridMercator(unit,grid,LLlat,LLlong,projection,     &
          proj_args,myNc,myNr,expandToLL,interpToLL)
       integer, intent(in) :: grid,unit
       integer, intent(out) :: projection,myNc,myNr
       real, dimension(100), intent(out) :: proj_args
       real, dimension(:), pointer :: LLlat,LLlong
       logical, intent(inout) :: expandToLL,interpToLL
     end subroutine gridMercator
  end interface

end module Gridder
