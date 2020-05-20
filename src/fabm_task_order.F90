#include "fabm_driver.h"

module fabm_task_order

   use fabm_graph
   use fabm_types, only: source_unknown, source_do_bottom, source_do_surface, source_do_horizontal, source2string
   use fabm_driver, only: driver

   implicit none

   private

   public find_best_order, type_step, type_graph_subset_node_pointer

   type type_graph_subset_node_pointer
      type (type_graph_subset_node),         pointer :: p    => null()
      type (type_graph_subset_node_pointer), pointer :: next => null()
   end type

   type type_graph_subset_node_set
      type (type_graph_subset_node_pointer), pointer :: first => null()
      integer                                        :: log_unit = -1
   contains
      procedure :: collect            => graph_subset_node_set_collect
      procedure :: collect_and_branch => graph_subset_node_set_collect_and_branch
      procedure :: branch             => graph_subset_node_set_branch
      procedure :: finalize           => graph_subset_node_set_finalize
   end type

   type type_graph_subset_node
      type (type_node), pointer         :: graph_node => null()
      integer                           :: outgoing_count = 0
      type (type_graph_subset_node_set) :: incoming
   end type

   type type_task_tree_node
      integer                             :: operation = source_unknown
      type (type_task_tree_node), pointer :: next_sibling => null()
      type (type_task_tree_node), pointer :: first_child  => null()
      type (type_task_tree_node), pointer :: parent       => null()
   contains
      procedure :: to_string         => task_tree_node_to_string
      procedure :: get_minimum_depth => task_tree_node_get_minimum_depth
      procedure :: get_leaf_at_depth => task_tree_node_get_leaf_at_depth
      procedure :: finalize          => task_tree_node_finalize
   end type

   type, extends(type_graph_subset_node_set) :: type_step
      integer :: operation = source_unknown
   end type

contains

   function find_best_order(graph, first_operation, log_unit) result(steps)
      type (type_graph), intent(in) :: graph
      integer,           intent(in) :: first_operation
      integer,           intent(in) :: log_unit
      type (type_step), allocatable :: steps(:)

      type (type_graph_subset_node_set)    :: subset
      class (type_task_tree_node), pointer :: leaf
      type (type_task_tree_node)           :: root
      integer                              :: ntasks
      integer                              :: itask

      ! Create tree that describes all possible task orders (root is at the right of the call list, i.e., at the end of the very last call)
      if (log_unit /= -1) write (log_unit,'(a)') '  possible task orders:'
      call create_graph_subset_node_set(graph, subset, log_unit)
      root%operation = first_operation
      if (root%operation /= source_unknown) then
         call subset%collect_and_branch(root)
      else
         call subset%branch(root)
      end if

      ! Determine optimum task order and use it to create the task objects.
      ntasks = root%get_minimum_depth()
      leaf => root%get_leaf_at_depth(ntasks)
      _ASSERT_(associated(leaf), 'find_best_order', 'BUG: get_leaf_at_depth did not return valid leaf.')

      if (root%operation == source_unknown) ntasks = ntasks - 1
      if (log_unit /= -1) write (log_unit,'(a,i0,a,a)') '  best task order (', ntasks, ' tasks): ', trim(leaf%to_string())

      allocate(steps(ntasks))
      do itask = 1, ntasks
         _ASSERT_(associated(leaf), 'find_best_order', 'BUG: leaf pointer not associated.')
         _ASSERT_(leaf%operation /= source_unknown, 'find_best_order', 'BUG: step with operation source_unknown.')
         steps(itask)%operation = leaf%operation
         leaf => leaf%parent
      end do
      call root%finalize()
      do itask = ntasks, 1, -1
         call subset%collect(steps(itask)%operation, steps(itask))
      end do
      _ASSERT_(.not. associated(subset%first), 'find_best_order', 'BUG: graph subset should be empty after create_tasks.')
   end function find_best_order

   subroutine create_graph_subset_node_set(graph, set, log_unit)
      type (type_graph),                 intent(in)  :: graph
      type (type_graph_subset_node_set), intent(out) :: set
      integer,                           intent(in)  :: log_unit

      type (type_node_list_member),         pointer :: graph_node
      type (type_node_set_member),          pointer :: dependency
      type (type_graph_subset_node_pointer),pointer :: set_node, set_node2, set_node3

      set%log_unit = log_unit

      ! Create representatives for each original graph node.
      graph_node => graph%first
      do while (associated(graph_node))
         allocate(set_node)
         allocate(set_node%p)
         set_node%p%graph_node => graph_node%p
         set_node%next => set%first
         set%first => set_node
         graph_node => graph_node%next
      end do

      set_node => set%first
      do while (associated(set_node))
         ! Enumerate dependencies of the original graph node
         dependency => set_node%p%graph_node%dependencies%first
         do while (associated(dependency))
            ! Find set node that corresponds to the dependency [set_node2]
            ! Note that this node does not exist if the dependency could be fulfilled by an earlier processed job/graph!
            set_node2 => set%first
            do while (associated(set_node2))
               if (associated(set_node2%p%graph_node, dependency%p)) exit
               set_node2 => set_node2%next
            end do

            if (associated(set_node2)) then
               ! Add dependency as "incoming" node.
               allocate(set_node3)
               set_node3%p => set_node2%p
               set_node3%next => set_node%p%incoming%first
               set_node%p%incoming%first => set_node3

               ! Increment the number of outgoing edges of the dependency itself.
               set_node2%p%outgoing_count = set_node2%p%outgoing_count + 1
            end if
            dependency => dependency%next
         end do
         set_node => set_node%next
      end do
   end subroutine create_graph_subset_node_set

   subroutine graph_subset_node_set_finalize(self)
      class (type_graph_subset_node_set), intent(inout) :: self

      type (type_graph_subset_node_pointer), pointer :: node, next_node, dependency, next_dependency

      node => self%first
      do while (associated(node))
         next_node => node%next
         dependency => node%p%incoming%first
         do while (associated(dependency))
            next_dependency => dependency%next
            deallocate(dependency)
            dependency => next_dependency
         end do
         deallocate(node%p)
         deallocate(node)
         node => next_node
      end do
      self%first => null()
   end subroutine graph_subset_node_set_finalize

   subroutine graph_subset_node_set_collect(self, operation, removed)
      class (type_graph_subset_node_set), intent(inout) :: self
      integer,                            intent(in)    :: operation
      class (type_graph_subset_node_set), intent(inout) :: removed

      type (type_graph_subset_node_pointer), pointer :: pnode, pnode2, pnode_previous, pnode_next
      logical                                        :: found, source_match
      integer                                        :: new_operation

      do
         found = .false.
         pnode_previous => null()
         pnode => self%first
         do while (associated(pnode))
            ! Store pointer to next node as we may remove the node from the set (thus destroying its original next pointer)
            pnode_next => pnode%next

            ! Determine whether the node's source matches
            new_operation = operation
            source_match = is_source_compatible(new_operation,pnode%p%graph_node%source)
            source_match = source_match .and. new_operation == operation ! Source is considered compatible only if it does not change the original task

            if (pnode%p%outgoing_count == 0 .and. source_match) then
               found = .true.

               ! Remove node from active set...
               if (associated(pnode_previous)) then
                  ! This is NOT the first node in the set.
                  pnode_previous%next => pnode%next
               else
                  ! This is the first node in the set.
                  self%first => pnode%next
               end if

               ! ...disconnect it from remaining nodes
               pnode2 => pnode%p%incoming%first
               do while (associated(pnode2))
                  pnode2%p%outgoing_count = pnode2%p%outgoing_count - 1
                  pnode2 => pnode2%next
               end do

               ! ...and add it to the set of removed nodes.
               pnode%next => removed%first
               removed%first => pnode
            else
               pnode_previous => pnode
            end if
            pnode => pnode_next
         end do
         if (.not. found) exit
      end do
   end subroutine graph_subset_node_set_collect

   recursive subroutine graph_subset_node_set_collect_and_branch(self, tree_node)
      class (type_graph_subset_node_set), intent(inout) :: self
      type (type_task_tree_node),target,  intent(inout) :: tree_node

      type (type_graph_subset_node_pointer), pointer :: pnode, pnode2, pnode_next
      type (type_graph_subset_node_set)              :: removed

      call self%collect(tree_node%operation, removed)

      if (self%log_unit /= -1 .and. .not. associated(self%first)) write (self%log_unit,'(a,a)') '  - ', trim(tree_node%to_string())

      ! We have processed all graph end points with compatible sources.
      ! Now process end points with other sources.
      call self%branch(tree_node)

      ! Restore removed nodes: add each removed node back to the set, and increment edge count for all nodes that connect to it.
      pnode => removed%first
      do while (associated(pnode))
         pnode_next => pnode%next
         pnode2 => pnode%p%incoming%first
         do while (associated(pnode2))
            pnode2%p%outgoing_count = pnode2%p%outgoing_count + 1
            pnode2 => pnode2%next
         end do
         pnode%next => self%first
         self%first => pnode
         pnode => pnode_next
      end do
   end subroutine graph_subset_node_set_collect_and_branch

   function task_tree_node_to_string(self) result(string)
      class (type_task_tree_node), intent(in) :: self
      character(len=512) :: string

      type (type_task_tree_node), pointer :: node

      string = source2string(self%operation)
      node => self%parent
      do while (associated(node))
         if (node%operation /= source_unknown) string = trim(string) // ', ' // trim(source2string(node%operation))
         node => node%parent
      end do
   end function task_tree_node_to_string

   recursive function task_tree_node_get_minimum_depth(self) result(depth)
      class (type_task_tree_node),intent(in) :: self
      integer :: depth

      type (type_task_tree_node), pointer :: child

      if (.not.associated(self%first_child)) then
         depth = 0
      else
         depth = huge(depth)
         child => self%first_child
         do while (associated(child))
            depth = min(depth, child%get_minimum_depth())
            child => child%next_sibling
         end do
      end if
      depth = depth + 1
   end function task_tree_node_get_minimum_depth

   recursive function task_tree_node_get_leaf_at_depth(self, length) result(node)
      class (type_task_tree_node), target, intent(in) :: self
      integer,                             intent(in) :: length
      class (type_task_tree_node), pointer :: node

      type (type_task_tree_node), pointer :: child

      node => null()
      if (length == 1) then
         if (.not. associated(self%first_child)) node => self
         return
      end if

      child => self%first_child
      do while (associated(child))
         node => child%get_leaf_at_depth(length - 1)
         if (associated(node)) return
         child => child%next_sibling
      end do
   end function

   recursive subroutine task_tree_node_finalize(self)
      class (type_task_tree_node), intent(inout) :: self

      type (type_task_tree_node), pointer :: child, child2

      child => self%first_child
      do while (associated(child))
         child2 => child%next_sibling
         call child%finalize()
         deallocate(child)
         child => child2
      end do
      self%first_child => null()
   end subroutine

   recursive subroutine graph_subset_node_set_branch(self, parent)
      class (type_graph_subset_node_set), intent(inout) :: self
      type (type_task_tree_node),target,  intent(inout) :: parent

      type (type_graph_subset_node_pointer), pointer :: pnode
      type (type_task_tree_node),            pointer :: pbranch
      integer                                        :: operation

      pnode => self%first
      do while (associated(pnode))
         if (pnode%p%outgoing_count==0) then
            operation = source2operation(pnode%p%graph_node%source)
            pbranch => parent%first_child
            do while (associated(pbranch))
               if (pbranch%operation == operation) exit
               pbranch => pbranch%next_sibling
            end do
            if (.not. associated(pbranch)) then
               allocate(pbranch)
               pbranch%operation = operation
               pbranch%parent => parent
               pbranch%next_sibling => parent%first_child
               parent%first_child => pbranch
            end if
         end if
         pnode => pnode%next
      end do

      pbranch => parent%first_child
      do while (associated(pbranch))
         call self%collect_and_branch(pbranch)
         pbranch => pbranch%next_sibling
      end do
   end subroutine graph_subset_node_set_branch

   logical function is_source_compatible(operation, new_source)
      integer, intent(inout) :: operation
      integer, intent(in)    :: new_source

      integer :: new_operation

      new_operation = source2operation(new_source)

      if (operation == source_unknown) then
         ! Current operation is still undefined. Adopt the current and return.
         is_source_compatible = .true.
         operation = new_operation
         return
      end if

      select case (new_operation)
      case (source_do_bottom)
         is_source_compatible = operation == source_do_horizontal .or. operation == source_do_bottom
         if (operation == source_do_horizontal) operation = source_do_bottom
      case (source_do_surface)
         is_source_compatible = operation == source_do_horizontal .or. operation == source_do_surface
         if (operation == source_do_horizontal) operation = source_do_surface
      case (source_do_horizontal)
         is_source_compatible = operation == source_do_horizontal .or. operation == source_do_surface &
            .or. operation == source_do_bottom
      case default
         is_source_compatible = operation == new_operation
      end select
   end function is_source_compatible

end module
