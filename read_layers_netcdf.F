        PROGRAM read_layers

        use netcdf

        implicit none

        INTEGER,PARAMETER :: nx = 240
        INTEGER,PARAMETER :: ny = 240
        INTEGER,PARAMETER :: lvl = 16
        INTEGER,PARAMETER :: lh=nx/2+1,ld=2*lh
        INTEGER,PARAMETER :: tini = 111000 
        INTEGER,PARAMETER :: tend = 160000
        INTEGER,PARAMETER :: length = (tend-tini)
        INTEGER,PARAMETER :: output_width = 1 
        INTEGER,PARAMETER :: avg_width = 100
        INTEGER,PARAMETER :: new_length = length/avg_width
        INTEGER,PARAMETER :: skip_length = length/100

        INTEGER :: ncid_w,varid_w,ierr_w,varid_t
        INTEGER :: ncid_r,varid_r,ierr_r 
        INTEGER :: dimids(3)
        INTEGER,ALLOCATABLE :: time(:)

        REAL,PARAMETER :: k0 = 4.0

        REAL(KIND=8) :: mean_w
        REAL(KIND=8) :: cov_w,var_w

        REAL(KIND=8),DIMENSION(nx,ny) :: w

        REAL(KIND=8),DIMENSION(nx,ny) :: w_slice,w_pert
        REAL(KIND=8),DIMENSION(nx,ny) :: w_pert_shift,w_pert_copy
        REAL(KIND=8),DIMENSION(nx,ny) :: correlation_w

        REAL(KIND=8),DIMENSION(nx,ny,new_length) :: w_full

        REAL(KIND=8),DIMENSION(nx) :: kx,ky

        REAL(KIND=8),DIMENSION(4*nx+15) :: wsave

        COMPLEX(KIND=8),DIMENSION(nx) :: cdat_w

        INTEGER i,j,k,l,m,n,r,s

        CHARACTER(len=28) :: fname,fname_w,fname_r

        ! Initialize the size of the time array based on the condition
        if (avg_width < 100) then
          allocate(time(skip_length))  ! Allocate the array with skip_length
        else
          allocate(time(new_length))  ! Allocate the array with new_length
        endif

        ! Create a NetCDF file
        WRITE (fname_w,'(a,I2.2,a,I4.4,a,I3.3,a)') &
          'w_z',(lvl),'_t',(avg_width),'_k',(int(k0)),'.nc'
        ierr_w = nf90_create(fname_w, NF90_CLOBBER, ncid_w)
        if (ierr_w /= NF90_NOERR) stop 'Error creating NetCDF file'
        ! Define the dimensions of the data (nx, ny, nt)
        ierr_w = nf90_def_dim(ncid_w, 'nx', nx, dimids(1))
        ierr_w = nf90_def_dim(ncid_w, 'ny', ny, dimids(2))
        if (avg_width<100) then
          ierr_w = nf90_def_dim(ncid_w, 'nt', skip_length, dimids(3))
        else
          ierr_w = nf90_def_dim(ncid_w, 'nt', new_length, dimids(3))
        endif

        ! Create the variable to store the data (3D array)
        ierr_w = nf90_def_var(ncid_w, 'w', NF90_DOUBLE, dimids, varid_w)
        ierr_w = nf90_def_var(ncid_w, 'time',NF90_INT,dimids(3),varid_t)
        ! End the definition mode (important to allow data to be written)
        ierr_w = nf90_enddef(ncid_w)

        WRITE (fname_r,'(a,I2.2,a,I4.4,a,I3.3,a)') &
          'r_z',(lvl),'_t',(avg_width),'_k',(int(k0)),'.nc'
        ierr_r = nf90_create(fname_r, NF90_CLOBBER, ncid_r)
        if (ierr_r /= NF90_NOERR) stop 'Error creating NetCDF file'
        ! Define the dimensions of the data (nx, ny, nt)
        ierr_r = nf90_def_dim(ncid_r, 'nx', nx, dimids(1))
        ierr_r = nf90_def_dim(ncid_r, 'ny', ny, dimids(2))
        if (avg_width<100) then
          ierr_r = nf90_def_dim(ncid_r, 'nt', skip_length, dimids(3))
        else
          ierr_r = nf90_def_dim(ncid_r, 'nt', new_length, dimids(3))
        endif
        ! Create the variable to store the data (3D array)
        ierr_r = nf90_def_var(ncid_r, 'r', NF90_DOUBLE, dimids, varid_r)
        ierr_r = nf90_def_var(ncid_r, 'time',NF90_INT,dimids(3),varid_t)
        ! End the definition mode (important to allow data to be written)
        ierr_r = nf90_enddef(ncid_r)

        s=0.0
        do n=1,new_length
          if (mod(tini + avg_width * (n - 1), 100) /= 0) then
            cycle
          else
            s=s+1      
            w_slice(:,:)=0.0
            do r=1,avg_width
              w(:,:) = 0.0           
              WRITE (fname,'(a,I3.3,a,I2.2,a,I7.7,a)') &
              'w_layer_z',(lvl),'_t',(output_width),'_', &
              (tini+(r-1)+(avg_width*(n-1))),'.out'
              open(10,file=fname,form='unformatted')
              READ (10) w(:,:)
              close(10)
              do i=1,nx
                do j=1,ny
                  w_slice(i,j)=w_slice(i,j)+w(i,j)
                enddo
              enddo
            enddo
            w_slice(:,:)=w_slice(:,:)/avg_width
            ! Compute 2D Fourier transform and shift zero frequency component to center
            do i=1,nx
              do j=1,ny
                if (j <= (nx/2)+1) then
                  ky(j) = j-1
                else
                  ky(j) = nx-j+1
                endif
              enddo
              cdat_w(:)=0.0
              do j=1,ny
                cdat_w(j)=CMPLX(w_slice(i,j),0.0)
              enddo
              call cffti(ny,wsave)
              call cfftf(ny,cdat_w,wsave)
              do j=1,ny
                if (ABS(ky(j)) > k0) then
                  cdat_w(j) = 0.0
                endif
              enddo
              call cfftb(ny,cdat_w,wsave)
              do j=1,ny
                w_slice(i,j)=dble(cdat_w(j))/dble(nx)
              enddo
            enddo
            do j=1,ny
              do i=1,nx
                if (i <= (nx/2)+1) then
                  kx(i) = i-1
                else
                  kx(i) = nx-i+1
                endif
              enddo
              cdat_w(:) = 0.0
              do i=1,nx
                cdat_w(i)=CMPLX(w_slice(i,j),0.0)
              enddo
              call cffti(nx,wsave)
              call cfftf(nx,cdat_w,wsave)
              do i=1,nx
                if (ABS(kx(i)) > k0) then
                  cdat_w(i) = 0.0
                endif
              enddo
              call cfftb(nx,cdat_w,wsave)
              do i=1,nx
                w_slice(i,j)=dble(cdat_w(i))/dble(ny)
              enddo
            enddo
            time(s)=(tini+(avg_width*(n-1)))
            ierr_w=nf90_put_var(ncid_w,varid_w,w_slice,start = [1, 1,s],count=[nx, ny, 1])
            mean_w=0.0
            do i=1,nx
              do j=1,ny
                mean_w=mean_w+w_slice(i,j)
              enddo
            enddo
            mean_w=mean_w/(nx*ny)
            w_pert=w_slice-mean_w
            do i=1,nx
              do l=1,nx
                do m=1,ny
                  w_pert_copy(l,m)=w_pert(l,m)
                enddo
              enddo
              w_pert_shift=w_pert_copy
              call circular_shift_x(w_pert_shift,nx,ny,i-(nx/2))
              do j=1,ny
                if (j==1) then
                  call circular_shift_y(w_pert_shift,nx,ny,j-(ny/2))
                else
                  call circular_shift_y(w_pert_shift,nx,ny,1)
                endif
                cov_w = 0.0
                var_w = 0.0
                do l=1,nx
                  do m=1,ny
                    cov_w=cov_w+(w_pert(l,m)*w_pert_shift(l,m))
                    var_w=var_w+(w_pert(l,m)**2)
                  enddo
                enddo
                cov_w=cov_w/(nx*ny)
                var_w=var_w/(nx*ny)
                correlation_w(i,j) = cov_w/var_w
              enddo
            enddo
            ierr_r=nf90_put_var(ncid_r,varid_r,correlation_w,start = [1,1,s],count=[nx, ny, 1])
          endif
        enddo
        print*,time
        ierr_w=nf90_put_var(ncid_w,varid_t,time)
        ierr_r=nf90_put_var(ncid_r,varid_t,time)
        ierr_w=nf90_close(ncid_w)
        ierr_r=nf90_close(ncid_r)
        END
        subroutine circular_shift_x(vector,nx,ny,shift)
        implicit none
        integer, intent(in) :: nx,ny           ! Size of array
        integer, intent(in) :: shift          ! Shift amount (positive or negative)
        integer :: i,new_shift
        real(kind=8), intent(inout) :: vector(nx,ny)
        real(kind=8) :: temp(nx,ny)

        ! Adjust shift to be within [1, size]
        new_shift = mod(shift,nx)
        ! Allocate temporary array
        ! Perform circular shift
        if (new_shift > 0) then
          temp(1:(nx-new_shift),:) = vector(new_shift+1:nx,:)
          temp((nx-new_shift)+1:nx,:) = vector(1:new_shift,:)
        elseif (new_shift < 0) then
          new_shift = nx - abs(new_shift)
          temp(1:(nx-new_shift),:) = vector(new_shift+1:nx,:)
          temp((nx-new_shift)+1:nx,:) = vector(1:new_shift,:)
        else
          temp = vector  ! No shift needed
        endif

        ! Copy shifted elements back to original vector
         vector = temp

        ! Deallocate temporary array

        end subroutine circular_shift_x
        subroutine circular_shift_y(vector,nx,ny,shift)
        implicit none
        integer, intent(in) :: nx,ny           ! Size of array
        integer, intent(in) :: shift          ! Shift amount (positive or negative)
        integer :: i,new_shift
        real(kind=8), intent(inout) :: vector(nx,ny)
        real(kind=8) :: temp(nx,ny)

        ! Adjust shift to be within [1, size]
        new_shift = mod(shift,ny)
        ! Allocate temporary array
        ! Perform circular shift
        if (new_shift > 0) then
          temp(:,1:(ny-new_shift)) = vector(:,new_shift+1:ny)
          temp(:,(ny-new_shift)+1:ny) = vector(:,1:new_shift)
        elseif (new_shift < 0) then
          new_shift = ny - abs(new_shift)
          temp(:,1:(ny-new_shift)) = vector(:,new_shift+1:ny)
          temp(:,(ny-new_shift)+1:ny) = vector(:,1:new_shift)
        else
          temp = vector  ! No shift needed
        endif

        ! Copy shifted elements back to original vector
         vector = temp

        ! Deallocate temporary array

        end subroutine circular_shift_y
