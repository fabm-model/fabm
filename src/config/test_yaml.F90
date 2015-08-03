program test_yaml

   use fabm_config_types
   use fabm_yaml
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