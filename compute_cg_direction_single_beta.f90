program xcompute_cg_direction
! ADIOS implementation: Ebru & Matthieu
! Princeton, September 2013

  use mpi
  use adios_read_mod
  use adios_write_mod
  use adios_helpers_mod

  implicit none

  include '../constants.h'
  include '../values_from_mesher.h'
  include '../precision.h'

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE
  integer,parameter:: NKERNEL=4
  integer:: myrank, sizeprocs,ier
  integer:: iker,ispec,i,j,k

  character(len=500):: direction_0_dir, direction_1_dir, gradient_0_dir, gradient_1_dir
  character(len=500):: direction_0_file, direction_1_file, gradient_0_file, gradient_1_file
  character(len=250):: kernel_name(NKERNEL)

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: direction_0, gradient_0,gradient_1,gradient_sq
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL):: direction_1
  real(kind=CUSTOM_REAL)::beta,beta_upper,beta_down,beta_upper_all_tmp,beta_down_all_tmp
  real(kind=CUSTOM_REAL),dimension(NKERNEL)::beta_upper_all,beta_down_all

  !--- Local parameters for ADIOS ---
  character(len=256) :: output_name,direction_bp_0_file,direction_bp_1_file
  character(len=256) :: grad_bp_0_file,grad_bp_1_file
  !character(len=64), parameter :: group_name  = "KERNELS_GROUP"
  character(len=150) :: group_name,kernels_bp(1)
  integer(kind=8) :: group, handle(NKERNEL), read_handle(NKERNEL)
  integer(kind=8) :: groupsize(NKERNEL), totalsize
  integer :: local_dim

  integer(kind=8) :: sel
  integer(kind=8), dimension(1) :: start, count_ad

  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  call adios_init_noxml (MPI_COMM_WORLD, ier);
  !sizeMB = 200 ! TODO 200MB is surely not the right size for the adios buffer
  call adios_allocate_buffer (ADIOS_BUFFER_SIZE_IN_MB, ier)

  call getarg(1,direction_0_dir)
  call getarg(2,direction_1_dir)
  call getarg(3,gradient_0_dir)
  call getarg(4,gradient_1_dir)

  if (trim(direction_0_dir) == '' .or. trim(direction_1_dir) == '' &
        .or. trim(gradient_0_dir) == '' .or. trim(gradient_1_dir) == '') then
        call exit_MPI(myrank,'USAGE: xcompute_cg_direction direction_0_dir direction_1_dir gradient_0_dir gradient_1_dir')
  end if

  !kernels_bp=(/"kernels_smooth.bp"/)
  kernels_bp=(/"kernels_precond.bp"/)

  kernel_name=(/character(len=150) :: "bulk_betah_kl_crust_mantle",&
                                      "bulk_betav_kl_crust_mantle",&
                                      "bulk_c_kl_crust_mantle",&
                                      "eta_kl_crust_mantle"/)

  !kernel_name=(/"reg1_bulk_betah_kernel_precond_smooth","reg1_bulk_betav_kernel_precond_smooth","reg1_bulk_c_kernel_precond_smooth","reg1_eta_kernel_precond_smooth"/)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD,"verbose=1", ier)

  write(grad_bp_0_file,'(a,i6.6,a)') trim(gradient_0_dir)//'/'//trim(kernels_bp(1))
  write(grad_bp_1_file,'(a,i6.6,a)') trim(gradient_1_dir)//'/'//trim(kernels_bp(1))
  !write(direction_bp_0_file,'(a,i6.6,a)') trim(direction_0_dir)//'/'//trim(kernels_bp(1))
  write(direction_bp_0_file,'(a,i6.6,a)') trim(direction_0_dir)//'/'//'cg_direction.bp'
  write(output_name,'(a,i6.6,a)') trim(direction_1_dir)//'/'//'cg_direction.bp'

  groupsize(1:NKERNEL) = 0
  group_name = "KERNELS_GROUP"

  call adios_declare_group(group, group_name, "", 1, ier)
  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)
  call define_adios_scalar(group, groupsize(1), "", "nspec", nspec)

  do iker=1,NKERNEL
    !-----------------------------------.
    ! Setup ADIOS for the current group |
    !-----------------------------------'
    !------------------------.
    ! Define ADIOS Variables |
    !------------------------'
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_array1D(group, groupsize(1), local_dim, &
                                   "", trim(kernel_name(iker)), direction_1)
  enddo

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  call adios_open(handle(1), group_name, output_name, "w", &
                  MPI_COMM_WORLD, ier)

  call adios_group_size (handle(1), groupsize(1), totalsize, ier)

  call adios_write(handle(1), "nspec", nspec, ier)


  do iker = 1,NKERNEL
  !        direction_1=0._CUSTOM_REAL

  !        write(direction_0_file,'(a,i6.6,a)') trim(direction_0_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
!        write(direction_1_file,'(a,i6.6,a)') trim(direction_1_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
!!        write(gradient_0_file,'(a,i6.6,a)') trim(gradient_0_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin' 
!!        write(gradient_1_file,'(a,i6.6,a)') trim(gradient_1_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin' 

    gradient_0=0.
    call adios_read_open_file (read_handle(iker), trim(grad_bp_0_file), 0, &
                               MPI_COMM_WORLD, ier)
    call adios_get_scalar(read_handle(iker), trim(kernel_name(iker))//"/local_dim",&
                          local_dim, ier)
    start(1) = local_dim * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec

    !print *, myrank, "box", start, count_ad
    call adios_selection_boundingbox(sel, 1, start, count_ad)

    call adios_schedule_read(read_handle(iker), sel, trim(kernel_name(iker))//"/array", &
                             0, 1, gradient_0, ier)
    call adios_perform_reads(read_handle(iker), ier)

    call adios_read_close(read_handle(iker),ier)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    gradient_1=0.
    call adios_read_open_file (read_handle(iker), trim(grad_bp_1_file), 0, &
                               MPI_COMM_WORLD, ier)
    call adios_get_scalar(read_handle(iker), trim(kernel_name(iker))//"/local_dim",&
                          local_dim, ier)

    start(1) = local_dim * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec

    !print *, myrank, "box", start, count_ad
    call adios_selection_boundingbox(sel, 1, start, count_ad)

    call adios_schedule_read(read_handle(iker), sel, trim(kernel_name(iker))//"/array", &
                             0, 1, gradient_1, ier)
    call adios_perform_reads(read_handle(iker), ier)

    call adios_read_close(read_handle(iker),ier)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    direction_0=0.
    call adios_read_open_file (read_handle(iker), trim(direction_bp_0_file), 0, &
                               MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle(iker), trim(kernel_name(iker))//"/local_dim",&
                          local_dim, ier)

    start(1) = local_dim * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec

    !print *, myrank, "box", start, count_ad
    call adios_selection_boundingbox(sel, 1, start, count_ad)

    call adios_schedule_read(read_handle(iker), sel, trim(kernel_name(iker))//"/array", &
                             0, 1, direction_0, ier)
    call adios_perform_reads(read_handle(iker), ier)

    call adios_read_close(read_handle(iker),ier)


!        if ( myrank == 0) print*,'reading direction0:',trim(direction_0_file)
!        if (ier /= 0 ) then 
!                print*, 'error opening:',trim(direction_0_file) 
!                call exit_mpi(myrank,'file not found') 
!        end if 
!        read(1001) direction_0(:,:,:,1:NSPEC) 
!        close(1001) 

!        open(1001,file=trim(gradient_0_file),status='old',form='unformatted',iostat=ier)
!        if (myrank == 0) print*, 'reading gradient0:',trim(gradient_0_file)
!        if ( ier /=0) then 
!                print*, 'error opening:',trim(gradient_0_file)
!                call exit_mpi(myrank,'file not found') 
!        end if 
!        read(1001) gradient_0(:,:,:,1:NSPEC)
!        close(1001)

!        open(1001,file=trim(gradient_1_file),status='old',form='unformatted',iostat=ier) 
!        if (myrank == 0) print*, 'reading gradient1:',trim(gradient_1_file)
!        if (ier/=0) then 
!                print*, 'error opening:',trim(gradient_1_file)
!                call exit_mpi(myrank,'file not found')
!        end if 
!        read(1001) gradient_1(:,:,:,1:NSPEC)
!        close(1001) 


!        beta_upper =sum(gradient_1(:,:,:,:)*(gradient_1(:,:,:,:)-gradient_0(:,:,:,:)))
!        beta_down=sum(gradient_0(:,:,:,:)*gradient_0(:,:,:,:)) 
    beta_upper=sum(gradient_1*(gradient_1-gradient_0))
    beta_down=sum(gradient_0*gradient_0)

    call mpi_barrier(MPI_COMM_WORLD,ier)
    call mpi_allreduce(beta_upper,beta_upper_all_tmp,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
    call mpi_allreduce(beta_down,beta_down_all_tmp,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)

    beta_upper_all(iker)=beta_upper_all_tmp
    beta_down_all(iker)=beta_down_all_tmp
  end do

  beta=sum(beta_upper_all)/sum(beta_down_all)

  if (myrank < 10 ) then
    print*,'before zero',myrank,beta
  end if

  if ( beta < 0.0 ) then
    beta=0.0
  end if

  beta=0.0 ! Ebru: command out for steepest descent

  if (myrank < 10 ) then
    print*,myrank,beta
  end if


  do iker = 1,NKERNEL

    direction_1=0._CUSTOM_REAL

    !=== read current gradient
    gradient_1=0.

    call adios_read_open_file (read_handle(iker), trim(grad_bp_1_file), 0, &
                               MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle(iker), trim(kernel_name(iker))//"/local_dim",&
                          local_dim, ier)

    start(1) = local_dim * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
    !print *, myrank, "box", start, count_ad

    call adios_selection_boundingbox(sel, 1, start, count_ad)

    call adios_schedule_read(read_handle(iker), sel, trim(kernel_name(iker))//"/array", &
                           0, 1, gradient_1, ier)

    call adios_perform_reads(read_handle(iker), ier)

    call adios_read_close(read_handle(iker),ier)


    !=== read previous direction
    direction_0=0.

    call adios_read_open_file (read_handle(iker), trim(direction_bp_0_file), 0, &
                               MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle(iker), trim(kernel_name(iker))//"/local_dim",&
                          local_dim, ier)

    start(1) = local_dim * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
    !print *, myrank, "box", start, count_ad

    call adios_selection_boundingbox(sel, 1, start, count_ad)

    call adios_schedule_read(read_handle(iker), sel, trim(kernel_name(iker))//"/array", &
                           0, 1, direction_0, ier)

    call adios_perform_reads(read_handle(iker), ier)

    call adios_read_close(read_handle(iker),ier)

    !=== calculate new direction
    do ispec=1,NSPEC
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            direction_1(i,j,k,ispec,iker) = -gradient_1(i,j,k,ispec) + beta * direction_0(i,j,k,ispec)
          end do
        end do
      end do
    end do ! ispec

    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(handle(1), myrank, sizeprocs, local_dim, &
                                     trim(kernel_name(iker)), direction_1(:,:,:,:,iker))

  end do ! kernel type

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle(1), "", ier)
  call adios_close(handle(1), ier)

  if (myrank==0) write(*,*) 'done calculting the search direction for all the kernels'

  call adios_finalize (myrank, ier)
  call MPI_FINALIZE(ier)

end program xcompute_cg_direction
