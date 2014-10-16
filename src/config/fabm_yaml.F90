module fabm_yaml

! This module parses a subset of YAML: http://yaml.org/
! Only block mappings are supported; flow mappings are not, nor are
! block or flow sequences. All scalars are left as strings; interpreting
! them as native data types is left to the caller. Thus, keys in the
! (key : value) pairs of a mapping are strings by definition.
! Comments (starting with #) are allowed. As per the YAML specification,
! indentation must consist of spaces only (no tabs!)
!
! For instance:
!
! instances:
!   P1:
!     model: pml/ersem/vphyt
!     parameters:
!       mu_max: 2.2
!       K: 1
!     coupling:
!       R: pom/R6
!
! NB This is not an attempt to write a complete YAML parser in Fortran!
! It handles a narrowly defined subset of YAML only. However, it IS meant
! to only accept documents that are valid YAML - if you find that this parser
! permits constructs that the YAML specification disallows (but not vice versa),
! please contact the author.
!
! Original author(s): Jorn Bruggeman

   use fabm_config_types

   implicit none

   private

   public parse,error_length

   integer,parameter :: line_length  = 2048
   integer,parameter :: error_length = 2048

   type type_file
      integer                 :: unit   = -1
      character(line_length)  :: line   = ''
      integer                 :: indent = 0
      logical                 :: eof    = .false.
      integer                 :: iline  = 0
      character(error_length) :: error_message = ''
      logical                 :: has_error     = .false.
   contains
      procedure :: next_line
      procedure :: set_error
   end type

contains

   function parse(path,unit,error) result(root)
      integer,                intent(in)  :: unit
      character(len=*),       intent(in)  :: path
      character(error_length),intent(out) :: error
      class (type_node),pointer           :: root

      type (type_file) :: file
      logical          :: already_open

      nullify(root)
      error = ''

      inquire(unit=unit, opened=already_open)
      if (.not.already_open) open(unit=unit,file=path,status='old',action='read',err=90)
      file%unit = unit
      file%eof = .false.
      call file%next_line()
      if (.not.file%has_error) root => read_value(file)
      if (.not.already_open) close(file%unit)
      if (file%has_error) then
         write (error,'(a,a,i0,a,a)') trim(path),', line ',file%iline,': ',trim(file%error_message)
      elseif (.not.file%eof) then
         if (associated(root)) then
            select type (root)
               class is (type_dictionary)
                  write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': unexpected decrease in indentation.'
               class is (type_scalar)
                  write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': expected end of file after reading &
                                             &one scalar value.'
               class default
                  write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': expected end of file.'
            end select
         else
            write (error,'(a,a,i0,a)') trim(path),', line ',file%iline,': expected end of file.'
         end if
      end if

      if (associated(root)) call root%set_path('')

      return

90    error = 'Unable to open '//trim(path)//' for reading.'
   end function

   subroutine next_line(file)
      class (type_file),intent(inout) :: file
      integer                         :: i
      logical                         :: done

      done = .false.
      do while (.not.done)
         ! Read entire line
         read (file%unit,'(A)',end=91) file%line
         file%iline = file%iline + 1

         ! Determine indentation and strip this.
         file%indent = len(file%line)
         do i=1,len(file%line)
            if (file%line(i:i)==achar(9)) then
               ! Found tabs in indentation: not allowed.
               call file%set_error('tab in indentation is not allowed.')
               return
            elseif (file%line(i:i)/=' ') then
               ! Found non-space: indentation ends here.
               file%indent = i-1
               exit
            end if
         end do
         file%line = file%line(file%indent+1:)

         ! If the line starts with comment character; move to next.
         if (file%line(1:1)=='#') cycle

         ! Search for whitespace delimited comment within the string; remove if found.
         do i=1,len_trim(file%line)-1
            if (is_whitespace(file%line(i:i)).and.file%line(i+1:i+1)=='#') then
               file%line = file%line(:i-1)
               exit
            end if
         end do

         ! Strip trailing whitespace
         do i=len(file%line),1,-1
            if (.not.is_whitespace(file%line(i:i))) then
               ! We found a non-whitespace character. Strip trailing whitespace and report we have a valid line.
               file%line = file%line(:i)
               done = .true.
               exit
            end if
         end do
      end do

      ! Check for unsupported YAML features.
      do i=1,len_trim(file%line)
         if (file%line(i:i)=='['.or.file%line(i:i)==']'.or.file%line(i:i)=='{'.or.file%line(i:i)=='}') then
            call file%set_error('flow mappings and sequences using []{} are not supported.')
            return
         end if
         if (file%line(i:i)=='"'.or.file%line(i:i)=='''') then
            call file%set_error('single- and double-quoted strings are not supported.')
            return
         end if
      end do
      do i=1,len(file%line)-1
         if (file%line(i:i+1)=='- ') then
            call file%set_error('block sequences using "- " are not supported.')
            return
         end if
      end do

      return

91    file%indent = 0
      file%eof = .true.
   end subroutine

   recursive function read_value(file) result(node)
      class (type_file),intent(inout) :: file
      class (type_node),pointer       :: node

      integer                    :: icolon,icolon_stop,firstindent
      type (type_key_value_pair) :: pair

      nullify(node)
      if (file%eof) return

      ! Find the first colon (if any)
      call find_mapping_character(file%line,icolon,icolon_stop)

      if (icolon==-1) then
         ! No colon found: item is a value
         allocate(type_scalar::node)
         select type (node)
            class is (type_scalar)
               node%string = trim(file%line)
         end select
         call file%next_line()
      else
         ! Colon found: item starts a mapping
         allocate(type_dictionary::node)
         firstindent = file%indent
         do
            pair = read_key_value_pair(file,icolon,icolon_stop)
            if (file%has_error) return
            select type (node)
               class is (type_dictionary)
                  call node%set(pair%key,pair%value)
            end select

            ! Check indentation of next line.
            if (file%indent>firstindent) then
               call file%set_error('unexpected increase in indentation following key-value pair "'//trim(pair%key)//'".')
               return
            elseif (file%eof .or. file%indent<firstindent) then
               ! End-of-file or decrease in indentation signifies that the mapping has ended.
               exit
            end if

            ! We are expecting a new key-value pair, since indentation has not changed. Find position of colon.
            call find_mapping_character(file%line,icolon,icolon_stop)
            if (icolon==-1) then
               call file%set_error('expected a key indicated by inline ": " or trailing :')
               return
            end if
         end do
      end if
   end function

   recursive function read_key_value_pair(file,icolon,icolon_stop) result(pair)
      class (type_file),intent(inout) :: file
      integer,          intent(in)    :: icolon,icolon_stop
      type (type_key_value_pair)      :: pair

      integer :: istop,baseindent

      istop = len_trim(file%line)

      pair%key = file%line(:icolon-1)
      if (icolon_stop==istop) then
         ! Colon ends the line; we need to read the value from the next line.
         baseindent = file%indent
         call file%next_line()
         if (file%has_error) return
         if (file%eof .or. file%indent<=baseindent) then
            ! Indentation equal to, or below, that of label (or file ends after label).
            ! That implies the value of the key-value pair is null.
            ! See YAML specification, section 7.2. Empty Nodes.
            allocate(type_null::pair%value)
         else
            ! Value on next line with higher indentation - read it.
            pair%value => read_value(file)
         end if
      else
         ! Value follows colon-space. Read the value and proceed to next line.
         allocate(type_scalar::pair%value)
         select type (node=>pair%value)
            class is (type_scalar)
               node%string = file%line(icolon_stop+1:istop)
         end select
         call file%next_line()
      end if
   end function

   subroutine find_mapping_character(string,istart,istop)
      character(len=*),intent(in)  :: string
      integer,         intent(out) :: istart,istop
      integer                      :: i,length

      ! Default: mapping indicator not found.
      istart = -1
      istop = -1

      ! Search for mapping indicator
      length = len_trim(string)
      do i=1,length-1
         if (string(i:i+1)==': ') then
            ! Found "colon space" mapping indicator
            istart = i
            exit
         end if
      end do

      ! No mapping indicator found yet; check whether string ends with colon.
      if (istart==-1 .and. string(length:length)==':') istart = length

      ! If we have not found a mapping indicator by now, there isn't one: return.
      if (istart==-1) return

      ! Eliminate all trailing whitespace
      istop = istart
      do i=istart+1,length
         if (.not.is_whitespace(string(i:i))) then
            istop = i-1
            exit
         end if
      end do

      ! Eliminate all preceding whitespace
      do i=istart-1,1,-1
         if (.not.is_whitespace(string(i:i))) then
            istart = i+1
            exit
         end if
      end do
   end subroutine

   logical function is_whitespace(string)
      character(len=*),intent(in) :: string
      ! White space in YAML includes spaces and tabs only (NB tabs are not allowed in indentation!)
      is_whitespace = (string(1:1)==' '.or.string(1:1)==achar(9))
   end function

   subroutine set_error(file,error)
      class (type_file),intent(inout) :: file
      character(len=*), intent(in)    :: error
      file%error_message = error
      file%has_error = .true.
   end subroutine

end module fabm_yaml
