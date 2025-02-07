        PROGRAM ensemble_average_profile

        IMPLICIT NONE

        INTEGER,PARAMETER :: nx = 240
        INTEGER,PARAMETER :: ny = 240
        INTEGER,PARAMETER :: nz = 161
        INTEGER,PARAMETER :: tini = 108000
        INTEGER,PARAMETER :: tend = 160000
        INTEGER,PARAMETER :: vf = 1000
        INTEGER,PARAMETER :: nt = tend/vf
        INTEGER,PARAMETER :: num = (tend-tini)/vf+1

        INTEGER :: i,j,k

        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_u,aver_v,aver_w
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_u2,aver_v2,aver_w2
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_u3,aver_v3,aver_w3
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_uw,aver_vw,aver_uv
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_uww,aver_uuw
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_txz,aver_tyz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_txx,aver_tyy,aver_tzz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_Cs,aver_beta_sgs
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_p
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_theta,aver_DTdz,aver_t2
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_wt,aver_sgs_t3
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_dudz,aver_dvdz
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_dudt,aver_dwdt
        REAL(KIND=8),DIMENSION(nt,nz,nx) :: aver_dissip

        REAL(KIND=8),DIMENSION(nz) :: u_yxt,v_yxt,w_yxt
        REAL(KIND=8),DIMENSION(nz) :: u2_yxt,v2_yxt,w2_yxt
        REAL(KIND=8),DIMENSION(nz) :: u3_yxt,v3_yxt,w3_yxt
        REAL(KIND=8),DIMENSION(nz) :: uw_yxt,vw_yxt,uv_yxt
        REAL(KIND=8),DIMENSION(nz) :: uww_yxt,uuw_yxt
        REAL(KIND=8),DIMENSION(nz) :: txz_yxt,tyz_yxt
        REAL(KIND=8),DIMENSION(nz) :: txx_yxt,tyy_yxt,tzz_yxt
        REAL(KIND=8),DIMENSION(nz) :: Cs_yxt,beta_sgs_yxt
        REAL(KIND=8),DIMENSION(nz) :: p_yxt
        REAL(KIND=8),DIMENSION(nz) :: theta_yxt,DTdz_yxt,t2_yxt
        REAL(KIND=8),DIMENSION(nz) :: wt_yxt,sgs_t3_yxt
        REAL(KIND=8),DIMENSION(nz) :: dudz_yxt,dvdz_yxt
        REAL(KIND=8),DIMENSION(nz) :: dudt_yxt,dwdt_yxt
        REAL(KIND=8),DIMENSION(nz) :: dissip_yxt

        REAL(KIND=8),DIMENSION(nt,nz) :: u_yx,v_yx,w_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: u2_yx,v2_yx,w2_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: u3_yx,v3_yx,w3_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: uw_yx,vw_yx,uv_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: uww_yx,uuw_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: txz_yx,tyz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: txx_yx,tyy_yx,tzz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: Cs_yx,beta_sgs_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: p_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: theta_yx,DTdz_yx,t2_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: wt_yx,sgs_t3_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: dudz_yx,dvdz_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: dudt_yx,dwdt_yx
        REAL(KIND=8),DIMENSION(nt,nz) :: dissip_yx

        REAL(KIND=8),DIMENSION(nt) :: time

        CHARACTER(len=30) :: fname

        u_yxt(:) = 0.0
        fname='aver_u.out'
        open(unit=1,file=fname,status='old',action='read')
        v_yxt(:) = 0.0
        fname='aver_v.out'
        open(unit=2,file=fname,status='old',action='read')
        w_yxt(:) = 0.0
        fname='aver_w.out'
        open(unit=3,file=fname,status='old',action='read')
        u2_yxt(:) = 0.0
        fname='aver_u2.out'
        open(unit=4,file=fname,status='old',action='read')
        v2_yxt(:) = 0.0
        fname='aver_v2.out'
        open(unit=5,file=fname,status='old',action='read')
        w2_yxt(:) = 0.0
        fname='aver_w2.out'
        open(unit=6,file=fname,status='old',action='read')
        u3_yxt(:) = 0.0
        fname='aver_u3.out'
        open(unit=7,file=fname,status='old',action='read')
        v3_yxt(:) = 0.0
        fname='aver_v3.out'
        open(unit=8,file=fname,status='old',action='read')
        w3_yxt(:) = 0.0
        fname='aver_w3.out'
        open(unit=9,file=fname,status='old',action='read')
        uw_yxt(:) = 0.0
        fname='aver_uw.out'
        open(unit=10,file=fname,status='old',action='read')
        vw_yxt(:) = 0.0
        fname='aver_vw.out'
        open(unit=11,file=fname,status='old',action='read')
        uv_yxt(:) = 0.0
        fname='aver_uv.out'
        open(unit=12,file=fname,status='old',action='read')
        uww_yxt(:) = 0.0
        fname='aver_uww.out'
        open(unit=13,file=fname,status='old',action='read')
        uuw_yxt(:) = 0.0
        fname='aver_uuw.out'
        open(unit=14,file=fname,status='old',action='read')
        txz_yxt(:) = 0.0
        fname='aver_txz.out'
        open(unit=15,file=fname,status='old',action='read')
        tyz_yxt(:) = 0.0
        fname='aver_tyz.out'
        open(unit=16,file=fname,status='old',action='read')
        txx_yxt(:) = 0.0
        fname='aver_txx.out'
        open(unit=17,file=fname,status='old',action='read')
        tyy_yxt(:) = 0.0
        fname='aver_tyy.out'
        open(unit=18,file=fname,status='old',action='read')
        tzz_yxt(:) = 0.0
        fname='aver_tzz.out'
        open(unit=19,file=fname,status='old',action='read')
        Cs_yxt(:) = 0.0
        fname='aver_Cs.out'
        open(unit=20,file=fname,status='old',action='read')
        beta_sgs_yxt(:) = 0.0
        fname='aver_beta_sgs.out'
        open(unit=21,file=fname,status='old',action='read')
        p_yxt(:) = 0.0
        fname='aver_p.out'
        open(unit=22,file=fname,status='old',action='read')
        theta_yxt(:) = 0.0
        fname='aver_theta.out'
        open(unit=23,file=fname,status='old',action='read')
        dTdz_yxt(:) = 0.0
        fname='aver_dTdz.out'
        open(unit=24,file=fname,status='old',action='read')
        t2_yxt(:)=0.0
        fname='aver_t2.out'
        open(unit=25,file=fname,status='old',action='read')
        wt_yxt(:)=0.0
        fname='aver_wt.out'
        open(unit=26,file=fname,status='old',action='read')
        sgs_t3_yxt(:) = 0.0
        fname='aver_sgs_t3.out'
        open(unit=27,file=fname,status='old',action='read')
        dudz_yxt(:)= 0.0
        fname='aver_dudz.out'
        open(unit=30,file=fname,status='old',action='read')
        dvdz_yxt(:) = 0.0
        fname='aver_dvdz.out'
        open(unit=31,file=fname,status='old',action='read')
        dudt_yxt(:) = 0.0
        fname='aver_dudt.out'
        open(unit=32,file=fname,status='old',action='read')
        dwdt_yxt(:) = 0.0
        fname='aver_dwdt.out'
        open(unit=33,file=fname,status='old',action='read')
        dissip_yxt(:) = 0.0
        fname='aver_dissip.out'
        open(unit=34,file=fname,status='old',action='read')

        do i=1,nt
          do j=1,nz-1
            read(1,*)time(i),(aver_u(i,j,k),k=1,nx)
            read(2,*)time(i),(aver_v(i,j,k),k=1,nx)
            read(3,*)time(i),(aver_w(i,j,k),k=1,nx)
            read(4,*)time(i),(aver_u2(i,j,k),k=1,nx)
            read(5,*)time(i),(aver_v2(i,j,k),k=1,nx)
            read(6,*)time(i),(aver_w2(i,j,k),k=1,nx)
            read(7,*)time(i),(aver_u3(i,j,k),k=1,nx)
            read(8,*)time(i),(aver_v3(i,j,k),k=1,nx)
            read(9,*)time(i),(aver_w3(i,j,k),k=1,nx)
            read(10,*)time(i),(aver_uw(i,j,k),k=1,nx)
            read(11,*)time(i),(aver_vw(i,j,k),k=1,nx)
            read(12,*)time(i),(aver_uv(i,j,k),k=1,nx)
            read(13,*)time(i),(aver_uww(i,j,k),k=1,nx)
            read(14,*)time(i),(aver_uuw(i,j,k),k=1,nx)
            read(15,*)time(i),(aver_txz(i,j,k),k=1,nx)
            read(16,*)time(i),(aver_tyz(i,j,k),k=1,nx)
            read(17,*)time(i),(aver_txx(i,j,k),k=1,nx)
            read(18,*)time(i),(aver_tyy(i,j,k),k=1,nx)
            read(19,*)time(i),(aver_tzz(i,j,k),k=1,nx)
            read(20,*)time(i),(aver_Cs(i,j,k),k=1,nx)
            read(21,*)time(i),(aver_beta_sgs(i,j,k),k=1,nx)
            read(22,*)time(i),(aver_p(i,j,k),k=1,nx)
            read(23,*)time(i),(aver_theta(i,j,k),k=1,nx)
            read(24,*)time(i),(aver_dTdz(i,j,k),k=1,nx)
            read(25,*)time(i),(aver_t2(i,j,k),k=1,nx)
            read(26,*)time(i),(aver_wt(i,j,k),k=1,nx)
            read(27,*)time(i),(aver_sgs_t3(i,j,k),k=1,nx)
            read(30,*)time(i),(aver_dudz(i,j,k),k=1,nx)
            read(31,*)time(i),(aver_dvdz(i,j,k),k=1,nx)
            read(32,*)time(i),(aver_dudt(i,j,k),k=1,nx)
            read(33,*)time(i),(aver_dwdt(i,j,k),k=1,nx)
            read(34,*)time(i),(aver_dissip(i,j,k),k=1,nx)
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
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(27)
        close(30)
        close(31)
        close(32)
        close(33)
        close(34)
        u_yx(:,:)=sum(aver_u(:,:,:),3)/nx
        v_yx(:,:)=sum(aver_v(:,:,:),3)/nx
        w_yx(:,:)=sum(aver_w(:,:,:),3)/nx
        u2_yx(:,:)=sum(aver_u2(:,:,:),3)/nx
        v2_yx(:,:)=sum(aver_v2(:,:,:),3)/nx
        w2_yx(:,:)=sum(aver_w2(:,:,:),3)/nx
        u3_yx(:,:)=sum(aver_u3(:,:,:),3)/nx
        v3_yx(:,:)=sum(aver_v3(:,:,:),3)/nx
        w3_yx(:,:)=sum(aver_w3(:,:,:),3)/nx
        uw_yx(:,:)=sum(aver_uw(:,:,:),3)/nx
        vw_yx(:,:)=sum(aver_vw(:,:,:),3)/nx
        uv_yx(:,:)=sum(aver_uv(:,:,:),3)/nx
        uww_yx(:,:)=sum(aver_uww(:,:,:),3)/nx
        uuw_yx(:,:)=sum(aver_uuw(:,:,:),3)/nx
        txz_yx(:,:)=sum(aver_txz(:,:,:),3)/nx
        tyz_yx(:,:)=sum(aver_tyz(:,:,:),3)/nx
        txx_yx(:,:)=sum(aver_txx(:,:,:),3)/nx
        tyy_yx(:,:)=sum(aver_tyy(:,:,:),3)/nx
        tzz_yx(:,:)=sum(aver_tzz(:,:,:),3)/nx
        Cs_yx(:,:)=sum(aver_Cs(:,:,:),3)/nx
        beta_sgs_yx(:,:)=sum(aver_beta_sgs(:,:,:),3)/nx
        p_yx(:,:)=sum(aver_p(:,:,:),3)/nx
        theta_yx(:,:)=sum(aver_theta(:,:,:),3)/nx
        dTdz_yx(:,:)=sum(aver_dTdz(:,:,:),3)/nx
        t2_yx(:,:)=sum(aver_t2(:,:,:),3)/nx
        wt_yx(:,:)=sum(aver_wt(:,:,:),3)/nx
        sgs_t3_yx(:,:)=sum(aver_sgs_t3(:,:,:),3)/nx
        dudz_yx(:,:)=sum(aver_dudz(:,:,:),3)/nx
        dvdz_yx(:,:)=sum(aver_dvdz(:,:,:),3)/nx
        dudt_yx(:,:)=sum(aver_dudt(:,:,:),3)/nx
        dwdt_yx(:,:)=sum(aver_dwdt(:,:,:),3)/nx
        dissip_yx(:,:)=sum(aver_dissip(:,:,:),3)/nx
        do j=1,nz-1
          do i=tini/vf,tend/vf
            u_yxt(j) = u_yxt(j) + u_yx(i,j)
            v_yxt(j) = v_yxt(j) + v_yx(i,j)
            w_yxt(j) = w_yxt(j) + w_yx(i,j)
            u2_yxt(j) = u2_yxt(j) + u2_yx(i,j)
            v2_yxt(j) = v2_yxt(j) + v2_yx(i,j)
            w2_yxt(j) = w2_yxt(j) + w2_yx(i,j)
            u3_yxt(j) = u3_yxt(j) + u3_yx(i,j)
            v3_yxt(j) = v3_yxt(j) + v3_yx(i,j)
            w3_yxt(j) = w3_yxt(j) + w3_yx(i,j)
            uw_yxt(j) = uw_yxt(j) + uw_yx(i,j)
            vw_yxt(j) = vw_yxt(j) + vw_yx(i,j)
            uv_yxt(j) = uv_yxt(j) + uv_yx(i,j)
            uww_yxt(j) = uww_yxt(j) + uww_yx(i,j)
            uuw_yxt(j) = uuw_yxt(j) + uuw_yx(i,j)
            txz_yxt(j) = txz_yxt(j) + txz_yx(i,j)
            tyz_yxt(j) = tyz_yxt(j) + tyz_yx(i,j)
            txx_yxt(j) = txx_yxt(j) + txx_yx(i,j)
            tyy_yxt(j) = tyy_yxt(j) + tyy_yx(i,j)
            tzz_yxt(j) = tzz_yxt(j) + tzz_yx(i,j)
            Cs_yxt(j) = Cs_yxt(j) + Cs_yx(i,j)
            beta_sgs_yxt(j) = beta_sgs_yxt(j) + beta_sgs_yx(i,j)
            p_yxt(j) = p_yxt(j) + p_yx(i,j)
            theta_yxt(j) = theta_yxt(j) + theta_yx(i,j)
            dTdz_yxt(j) = dTdz_yxt(j) + dTdz_yx(i,j)
            t2_yxt(j) = t2_yxt(j) + t2_yx(i,j)
            wt_yxt(j) = wt_yxt(j) + wt_yx(i,j)
            sgs_t3_yxt(j) = sgs_t3_yxt(j) + sgs_t3_yx(i,j)
            dudz_yxt(j) = dudz_yxt(j) + dudz_yx(i,j)
            dvdz_yxt(j) = dvdz_yxt(j) + dvdz_yx(i,j)
            dudt_yxt(j)= dudt_yxt(j) + dudt_yx(i,j)
            dwdt_yxt(j)= dwdt_yxt(j) + dwdt_yx(i,j)
            dissip_yxt(j) = dissip_yxt(j) + dissip_yx(i,j)
          end do
        end do
        u_yxt(:) = u_yxt(:)/num
        v_yxt(:) = v_yxt(:)/num
        w_yxt(:) = w_yxt(:)/num
        u2_yxt(:) = u2_yxt(:)/num
        v2_yxt(:) = v2_yxt(:)/num
        w2_yxt(:) = w2_yxt(:)/num
        u3_yxt(:) = u3_yxt(:)/num
        v3_yxt(:) = v3_yxt(:)/num
        w3_yxt(:) = w3_yxt(:)/num
        uw_yxt(:) = uw_yxt(:)/num
        vw_yxt(:) = vw_yxt(:)/num
        uv_yxt(:) = uv_yxt(:)/num
        uww_yxt(:) = uww_yxt(:)/num
        uuw_yxt(:) = uuw_yxt(:)/num
        txz_yxt(:) = txz_yxt(:)/num
        tyz_yxt(:) = tyz_yxt(:)/num
        txx_yxt(:) = txx_yxt(:)/num
        tyy_yxt(:) = tyy_yxt(:)/num
        tzz_yxt(:) = tzz_yxt(:)/num
        Cs_yxt(:) = Cs_yxt(:)/num
        beta_sgs_yxt(:) = beta_sgs_yxt(:)/num
        p_yxt(:) = p_yxt(:)/num
        theta_yxt(:) = theta_yxt(:)/num
        dTdz_yxt(:) = dTdz_yxt(:)/num
        t2_yxt(:)= t2_yxt(:)/num
        wt_yxt(:) = wt_yxt(:)/num
        sgs_t3_yxt(:) = sgs_t3_yxt(:)/num
        dudz_yxt(:) = dudz_yxt(:)/num
        dvdz_yxt(:) = dvdz_yxt(:)/num
        dudt_yxt(:) = dudt_yxt(:)/num
        dwdt_yxt(:) = dwdt_yxt(:)/num
        dissip_yxt(:) = dissip_yxt(:)/num
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'u_z_',tini,'_',tend,'.txt'
        open(1,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(1,'(F20.14)')(u_yxt(k))
        end do
        close(1)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'v_z_',tini,'_',tend,'.txt'
        open(2,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(2,'(F20.14)')(v_yxt(k))
        end do
        close(2)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'w_z_',tini,'_',tend,'.txt'
        open(3,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(3,'(F20.14)')(w_yxt(k))
        end do
        close(3)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'u2_z_',tini,'_',tend,'.txt'
        open(4,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(4,'(F20.14)')(u2_yxt(k))
        end do
        close(4)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'v2_z_',tini,'_',tend,'.txt'
        open(5,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(5,'(F20.14)')(v2_yxt(k))
        end do
        close(5)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'w2_z_',tini,'_',tend,'.txt'
        open(6,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(6,'(F20.14)')(w2_yxt(k))
        end do
        close(6)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'u3_z_',tini,'_',tend,'.txt'
        open(7,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(7,'(F20.14)')(u3_yxt(k))
        end do
        close(7)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'v3_z_',tini,'_',tend,'.txt'
        open(8,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(8,'(F20.14)')(v3_yxt(k))
        end do
        close(8)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'w3_z_',tini,'_',tend,'.txt'
        open(9,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(9,'(F20.14)')(w3_yxt(k))
        end do
        close(9)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'uw_z_',tini,'_',tend,'.txt'
        open(10,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(10,'(F20.14)')(uw_yxt(k))
        end do
        close(10)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'vw_z_',tini,'_',tend,'.txt'
        open(11,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(11,'(F20.14)')(vw_yxt(k))
        end do
        close(11)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'uv_z_',tini,'_',tend,'.txt'
        open(12,file=fname,status='unknown',action='write')
        do k=1,nz
          write(12,'(F20.14)')(uv_yxt(k))
        enddo
        close(12)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'uww_z_',tini,'_',tend,'.txt'
        open(13,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(13,'(F20.14)')(uww_yxt(k))
        end do
        close(13)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'uuw_z_',tini,'_',tend,'.txt'
        open(14,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(14,'(F20.14)')(uuw_yxt(k))
        end do
        close(14)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'txz_z_',tini,'_',tend,'.txt'
        open(15,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(15,'(F20.14)')(txz_yxt(k))
        end do
        close(15)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'tyz_z_',tini,'_',tend,'.txt'
        open(16,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(16,'(F20.14)')(tyz_yxt(k))
        end do
        close(16)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'txx_z_',tini,'_',tend,'.txt'
        open(17,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(17,'(F20.14)')(txx_yxt(k))
        end do
        close(17)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'tyy_z_',tini,'_',tend,'.txt'
        open(18,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(18,'(F20.14)')(tyy_yxt(k))
        end do
        close(18)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'tzz_z_',tini,'_',tend,'.txt'
        open(19,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(19,'(F20.14)')(tzz_yxt(k))
        end do
        close(19)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'Cs_z_',tini,'_',tend,'.txt'
        open(20,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(20,'(F20.14)')(Cs_yxt(k))
        end do
        close(20)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'beta_sgs_z_',tini,'_',tend,'.txt'
        open(21,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(21,'(F20.14)')(beta_sgs_yxt(k))
        end do
        close(21)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'p_z_',tini,'_',tend,'.txt'
        open(22,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(22,'(F20.14)')(p_yxt(k))
        end do
        close(22)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'theta_z_',tini,'_',tend,'.txt'
        open(23,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(23,'(F20.14)')(theta_yxt(k))
        end do
        close(23)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dTdz_z_',tini,'_',tend,'.txt'
        open(24,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(24,'(F20.14)')(dTdz_yxt(k))
        end do
        close(24)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             't2_z_',tini,'_',tend,'.txt'
        open(25,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(25,'(F20.14)')(t2_yxt(k))
        end do
        close(25)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'wt_z_',tini,'_',tend,'.txt'
        open(26,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(26,'(F20.14)')(wt_yxt(k))
        end do
        close(26)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'sgs_t3_z_',tini,'_',tend,'.txt'
        open(27,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(27,'(F20.14)')(sgs_t3_yxt(k))
        end do
        close(27)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dudz_z_',tini,'_',tend,'.txt'
        open(30,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(30,'(F20.14)')(dudz_yxt(k))
        end do
        close(30)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dvdz_z_',tini,'_',tend,'.txt'
        open(31,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(31,'(F20.14)')(dvdz_yxt(k))
        end do
        close(31)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dudt_z_',tini,'_',tend,'.txt'
        open(32,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(32,'(F20.14)')(dudt_yxt(k))
        end do
        close(32)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dwdt_z_',tini,'_',tend,'.txt'
        open(33,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(33,'(F20.14)')(dwdt_yxt(k))
        end do
        close(33)
        WRITE(fname,'(A,I7.7,A,I7.7,A)'), &
             'dissip_z_',tini,'_',tend,'.txt'
        open(34,file=fname,status='unknown',action='write')
        do k = 1,nz
          write(34,'(F20.14)')(dissip_yxt(k))
        end do
        close(34)

        END PROGRAM ensemble_average_profile
