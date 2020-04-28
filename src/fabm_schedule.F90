module fabm_schedule

   use fabm_types, only: type_base_model, source_unknown, type_link, prefill_previous_value, operator_assign, presence_internal
   use fabm_parameters, only: rke

   implicit none

   public type_schedules, schedule_pattern_monthly

   private

   type type_logical_pointer
      logical, pointer                     :: p    => null()
      type (type_logical_pointer), pointer :: next => null()
   end type

   integer, parameter :: schedule_pattern_none    = -1
   integer, parameter :: schedule_pattern_monthly = 1

   type type_schedule
      class (type_base_model),     pointer :: model          => null()
      integer                              :: source         =  source_unknown
      integer                              :: pattern        =  schedule_pattern_none
      integer                              :: last_year      = -1
      integer                              :: last_month     = -1
      integer                              :: last_day       = -1
      type (type_logical_pointer), pointer :: first_affected => null()
      type (type_schedule),        pointer :: next           => null()
   end type

   type type_schedules
      type (type_schedule), pointer :: first => null()
   contains
      procedure :: add
      procedure :: attach
      procedure :: update
   end type

contains

   subroutine add(self, model, source, pattern)
      class (type_schedules), intent(inout) :: self
      class (type_base_model), target :: model
      integer, intent(in) :: source
      integer,                intent(in)    :: pattern

      type (type_schedule), pointer :: schedule
      type (type_link),     pointer :: link

      allocate(schedule)
      schedule%model => model
      schedule%source = source
      schedule%pattern = pattern
      schedule%next => self%first
      self%first => schedule

      link => model%links%first
      do while (associated(link))
         if (associated(link%original, link%target) .and. link%target%source == source &
             .and. link%target%presence == presence_internal .and. link%target%write_operator == operator_assign) &
               link%target%prefill = prefill_previous_value
         link => link%next
      end do
   end subroutine

   subroutine attach(self, model, source, active)
      class (type_schedules), intent(inout) :: self
      class (type_base_model), target :: model
      integer, intent(in) :: source
      logical,                 target  :: active

      type (type_schedule),        pointer :: schedule
      type (type_logical_pointer), pointer :: pactive

      schedule => self%first
      do while (associated(schedule))
         if (associated(schedule%model, model) .and. schedule%source == source) exit
         schedule => schedule%next
      end do

      if (associated(schedule)) then
         allocate(pactive)
         pactive%p => active
         pactive%next => schedule%first_affected
         schedule%first_affected => pactive
      end if
   end subroutine

   subroutine update(self, year, month, day, seconds)
      class (type_schedules), intent(inout) :: self
      integer,                intent(in)    :: year, month, day
      real(rke),              intent(in)    :: seconds

      type (type_schedule),        pointer :: schedule
      logical                              :: active
      type (type_logical_pointer), pointer :: pactive

      schedule => self%first
      do while (associated(schedule))
         select case (schedule%pattern)
         case (schedule_pattern_monthly)
            active = month > schedule%last_month .or. year > schedule%last_year
         end select
         if (active) then
            schedule%last_year  = year
            schedule%last_month = month
            schedule%last_day   = day
         end if
         pactive => schedule%first_affected
         do while (associated(pactive))
            pactive%p = active
            pactive => pactive%next
         end do
         schedule => schedule%next
      end do
   end subroutine

end module
