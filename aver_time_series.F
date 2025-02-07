        PROGRAM aver_time_series

        IMPLICIT NONE

        ! Parameters
        INTEGER,PARAMETER :: nx = 240 ! x-direction number of grid points
        INTEGER,PARAMETER :: nz = 161 ! z-direction number of grid points
        INTEGER,PARAMETER :: tend = 160000 ! total number of timesteps
        INTEGER,PARAMETER :: vf = 1000 ! number of timesteps averaged  in aver.out
        INTEGER,PARAMETER :: nt = tend/vf ! number of time data points in aver.out  
        REAL(KIND=8),PARAMETER :: Lz = 4000 ! vertical domain
        REAL(KIND=8),PARAMETER :: T_scale = 300 ! temperature scale (param.f90)
        REAL(KIND=8),PARAMETER :: wt_s = 0.05 ! domain average surface heat flux
        REAL(KIND=8),PARAMETER :: ug = 15 ! geopstrophic wind magnitude
        REAL(KIND=8),PARAMETER :: g = 9.81 ! accelatiion due to gravity

        ! loop integers
        INTEGER :: i,j,k,f

        ! Defining variable names and tpes
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_u,aver_v,aver_w
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_uw,aver_vw
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_u2,aver_v2,aver_w2
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_dudz,aver_dvdz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_txz,aver_tyz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_theta,aver_dTdz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_wt,aver_sgs_t3
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_txx,aver_tyy,aver_tzz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_dissip

        REAL(KIND=8),DIMENSION(nt,nz) :: u_yx,v_yx,w_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: uw_yx,vw_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: u2_yx,v2_yx,w2_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: dudz_yx,dvdz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: txz_yx,tyz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: theta_yx,dTdz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: wt_yx,sgs_t3_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: txx_yx,tyy_yx,tzz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: dissip_yx

        REAL(KIND=8),DIMENSION(nt,nz) :: uw_prime,vw_prime
        REAL(KIND=8),DIMENSION(nt,nz) :: bouy_TKE_prod,heat_flux
        REAL(KIND=8),DIMENSION(nt,nz) :: ushear_TKE_prod,vshear_TKE_prod
        REAL(KIND=8),DIMENSION(nt,nz) :: theta_temp,u_temp,v_temp
        REAL(KIND=8),DIMENSION(nt,nz) :: u2_prime,v2_prime,w2_prime,TKE
        REAL(KIND=8),DIMENSION(nt,nz) :: duwdz,dvwdz
        REAL(KIND=8),DIMENSION(nt) :: height,height_index
        REAL(KIND=8),DIMENSION(nt) :: w_star,u_star
        REAL(KIND=8),DIMENSION(nt) :: time

        CHARACTER(len=30) :: fname

        REAL(KIND=8) :: ABL_T

        ! Calculated parameters
        REAL(KIND=8) :: dz = Lz/DBLE(nz-1)

        ! Reading in data from aver.out files
        fname = 'aver_u.out'
        open(unit=1,file=fname,status='old',action='read')
        fname = 'aver_v.out'
        open(unit=2,file=fname,status='old',action='read')
        fname = 'aver_w.out'
        open(unit=3,file=fname,status='old',action='read')
        fname = 'aver_uw.out'
        open(unit=4,file=fname,status='old',action='read')
        fname = 'aver_vw.out'
        open(unit=5,file=fname,status='old',action='read')
        fname = 'aver_u2.out'
        open(unit=6,file=fname,status='old',action='read')
        fname = 'aver_v2.out'
        open(unit=7,file=fname,status='old',action='read')
        fname = 'aver_w2.out'
        open(unit=8,file=fname,status='old',action='read')
        fname = 'aver_dudz.out'
        open(unit=9,file=fname,status='old',action='read')
        fname = 'aver_dvdz.out'
        open(unit=10,file=fname,status='old',action='read')
        fname = 'aver_txz.out'
        open(unit=11,file=fname,status='old',action='read')
        fname = 'aver_tyz.out'
        open(unit=12,file=fname,status='old',action='read')
        fname = 'aver_theta.out'
        open(unit=13,file=fname,status='old',action='read')
        fname = 'aver_dTdz.out'
        open(unit=14,file=fname,status='old',action='read')
        fname = 'aver_wt.out'
        open(unit=15,file=fname,status='old',action='read')
        fname = 'aver_sgs_t3.out'
        open(unit=16,file=fname,status='old',action='read')
        fname = 'aver_txx.out'
        open(unit=17,file=fname,status='old',action='read')
        fname = 'aver_tyy.out'
        open(unit=18,file=fname,status='old',action='read')
        fname = 'aver_tzz.out'
        open(unit=19,file=fname,status='old',action='read')
        fname = 'aver_dissip.out'
        open(unit=20,file=fname,status='old',action='read')

        do i=1,nt
          do j=1,nz-1
            read(1,*)time(i),(aver_u(i,j,k),k=1,nx)
            read(2,*)time(i),(aver_v(i,j,k),k=1,nx)
            read(3,*)time(i),(aver_w(i,j,k),k=1,nx)
            read(4,*)time(i),(aver_uw(i,j,k),k=1,nx)
            read(5,*)time(i),(aver_vw(i,j,k),k=1,nx)
            read(6,*)time(i),(aver_u2(i,j,k),k=1,nx)
            read(7,*)time(i),(aver_v2(i,j,k),k=1,nx)
            read(8,*)time(i),(aver_w2(i,j,k),k=1,nx)
            read(9,*)time(i),(aver_dudz(i,j,k),k=1,nx)
            read(10,*)time(i),(aver_dvdz(i,j,k),k=1,nx)
            read(11,*)time(i),(aver_txz(i,j,k),k=1,nx)
            read(12,*)time(i),(aver_tyz(i,j,k),k=1,nx)
            read(13,*)time(i),(aver_theta(i,j,k),k=1,nx)
            read(14,*)time(i),(aver_dTdz(i,j,k),k=1,nx)
            read(15,*)time(i),(aver_wt(i,j,k),k=1,nx)
            read(16,*)time(i),(aver_sgs_t3(i,j,k),k=1,nx)
            read(17,*)time(i),(aver_txx(i,j,k),k=1,nx)
            read(18,*)time(i),(aver_tyy(i,j,k),k=1,nx)
            read(19,*)time(i),(aver_tzz(i,j,k),k=1,nx)
            read(20,*)time(i),(aver_dissip(i,j,k),k=1,nx)
          end do
        end do
        
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
        close(7)
        close(8)
        close(9)
        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
        close(17)
        close(18)
        close(19)
        close(20)

        ! averaging data over x
        u_yx(:,:) = sum(aver_u(:,:,:),3)/nx
        v_yx(:,:) = sum(aver_v(:,:,:),3)/nx
        w_yx(:,:) = sum(aver_w(:,:,:),3)/nx
        uw_yx(:,:) = sum(aver_uw(:,:,:),3)/nx
        vw_yx(:,:) = sum(aver_vw(:,:,:),3)/nx
        u2_yx(:,:) = sum(aver_u2(:,:,:),3)/nx
        v2_yx(:,:) =sum(aver_v2(:,:,:),3)/nx
        w2_yx(:,:) = sum(aver_w2(:,:,:),3)/nx
        dudz_yx(:,:) = sum(aver_dudz(:,:,:),3)/nx
        dvdz_yx(:,:) = sum(aver_dvdz(:,:,:),3)/nx
        txz_yx(:,:) = sum(aver_txz(:,:,:),3)/nx
        tyz_yx(:,:) = sum(aver_tyz(:,:,:),3)/nx
        theta_yx(:,:) = sum(aver_theta(:,:,:),3)/nx
        dTdz_yx(:,:) = sum(aver_dTdz(:,:,:),3)/nx
        wt_yx(:,:) = sum(aver_wt(:,:,:),3)/nx
        sgs_t3_yx(:,:) = sum(aver_sgs_t3(:,:,:),3)/nx
        txx_yx(:,:) = sum(aver_txx(:,:,:),3)/nx
        tyy_yx(:,:) = sum(aver_tyy(:,:,:),3)/nx
        tzz_yx(:,:) = sum(aver_tzz(:,:,:),3)/nx
        dissip_yx(:,:) = sum(aver_dissip(:,:,:),3)/nx
        
        ! u,v,theta on w-grid
        theta_temp(:,2:nz) = (theta_yx(:,2:nz)+theta_yx(:,1:(nz-1)))/2
        theta_temp(:,1) = ((theta_yx(:,1)*T_scale)- &
        ((dTdz_yx(:,1)*(T_scale/Lz))*dz/2))/T_scale
        u_temp(:,2:nz) = (u_yx(:,2:nz)+u_yx(:,1:(nz-1)))/2
        u_temp(:,1) = 0.0
        v_temp(:,2:nz) = (v_yx(:,2:nz)+v_yx(:,1:(nz-1)))/2
        v_temp(:,1) = 0.0

        ! ABL depth
        height_index = maxloc(dTdz_yx,2)
        height = (height_index-1)*dz

        ! convective velocity calculation
        do i=1,nt
          ABL_T = 0.0
          ABL_T = T_scale*(sum(theta_temp(i,1:height_index(i))) & ! Average ABL Potential Temp
          /height_index(i)) 
          w_star(i)=(((9.8/ABL_T)*height(i)*wt_s)**(1./3.))/ug !Nondimensional
        enddo
        height = height/Lz !Nondimensional

        ! horizontal turbulent momentum fluxes
        uw_prime = uw_yx(:,:)-(u_temp(:,:)*w_yx(:,:))+txz_yx(:,:)
        vw_prime = vw_yx(:,:)-(v_temp(:,:)*w_yx(:,:))+tyz_yx(:,:)

        ! Shear production of TKE
        ushear_TKE_prod = uw_prime(:,:)*dudz_yx(:,:)
        vshear_TKE_prod = vw_prime(:,:)*dvdz_yx(:,:)

        ! TKE calculation
        u2_prime = u2_yx(:,:)-(u_yx(:,:)*u_yx(:,:))+txx_yx(:,:)
        v2_prime = v2_yx(:,:)-(v_yx(:,:)*v_yx(:,:))+tyy_yx(:,:)
        w2_prime = w2_yx(:,:)-(w_yx(:,:)*w_yx(:,:))+tzz_yx(:,:)
        TKE = 0.5*(u2_prime+v2_prime+w2_prime)

        ! Bouyant production of TKE
        heat_flux = wt_yx(:,:)-(theta_temp(:,:)*w_yx(:,:)) &
        +sgs_t3_yx(:,:)
        bouy_TKE_prod = (g*heat_flux(:,:))/theta_yx(:,:)

        ! horizontal turbulent momentum flux vertical divergence
        duwdz(:,1:nz-1) = (uw_prime(:,2:nz)-uw_prime(:,1:nz-1))/dz
        duwdz(:,nz) = 0.0
        duwdz = duwdz*Lz
        dvwdz(:,1:nz-1) = (vw_prime(:,2:nz)-vw_prime(:,1:nz-1))/dz
        dvwdz(:,nz) = 0.0
        dvwdz = dvwdz*Lz

        ! friction velocity
        u_star(:)=((uw_prime(:,1)**2)+(vw_prime(:,1)**2))**(1./4.)

        ! output data
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dissip_prod_TKE_timeseries.txt'
        open(1,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(1,'(800F20.14)')((dissip_yx(f,k)),f=1,nt)
        end do
        close(1)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'ushear_prod_TKE_timeseries.txt'
        open(2,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(2,'(800F20.14)')((ushear_TKE_prod(f,k)),f=1,nt)
        end do
        close(2)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'vshear_prod_TKE_timeseries.txt'
        open(3,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(3,'(800F20.14)')((vshear_TKE_prod(f,k)),f=1,nt)
        end do
        close(3)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'bouy_prod_TKE_timeseries.txt'
        open(4,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(4,'(800F20.14)')((bouy_TKE_prod(f,k)),f=1,nt)
        end do
        close(4)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'TKE_profile_timeseries.txt'
        open(5,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(5,'(800F20.14)')((TKE(f,k)),f=1,nt)
        end do
        close(5)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dvwdz_profile_timeseries.txt'
        open(6,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(6,'(800F20.14)')((dvwdz(f,k)),f=1,nt)
        end do
        close(6)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'duwdz_profile_timeseries.txt'
        open(7,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(7,'(800F20.14)')((duwdz(f,k)),f=1,nt)
        end do
        close(7)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'u_profile_timeseries.txt'
        open(8,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(8,'(800F20.14)')((u_yx(f,k)),f=1,nt)
        end do
        close(8)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'v_profile_timeseries.txt'
        open(9,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(9,'(800F20.14)')((v_yx(f,k)),f=1,nt)
        end do
        close(9)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'zi_timeseries.txt'
        open(10,file=fname,status='unknown',action='write')
        do f = 1,nt
          write(10,'(F20.14)')(height(f))
        end do
        close(10)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'w_star_timeseries.txt'
        open(11,file=fname,status='unknown',action='write')
        do f = 1,nt
          write(11,'(F20.14)')(w_star(f))
        end do
        close(11)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'u_star_timeseries.txt'
        open(12,file=fname,status='unknown',action='write')
        do f = 1,nt
          write(12,'(F20.14)')(u_star(f))
        end do
        close(12)
        END PROGRAM aver_time_series
