! -----------------------------------------------------------------------------
! This file is part of Fortran-YAML: a lightweight YAML parser written in
! object-oriented Fortran.
!
! Official repository: https://github.com/BoldingBruggeman/fortran-yaml
!
! Copyright 2013-2016 Bolding & Bruggeman ApS.
!
! This is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation (https://www.gnu.org/licenses/gpl.html). It is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! A copy of the license is provided in the COPYING file.
! -----------------------------------------------------------------------------

program test_yaml

   use yaml_types
   use yaml
   use, intrinsic :: iso_fortran_env

   character(error_length) :: error
   character(256) :: path
   class (type_node),pointer :: root

   call get_command_argument(1, path)
   if (path=='') then
      write (*,*) 'ERROR: path to YAML file not provided.'
      stop 2
   end if
   root => parse(path,unit=100,error=error)
   if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
      stop 1
   end if
   call root%dump(unit=output_unit,indent=0)

end program
