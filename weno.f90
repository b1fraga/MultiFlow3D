!######################################################################
      subroutine HJ_WENO_dx(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          INTEGER :: I,J,K,L,k_star,npp,nq,nr,op,ib,pll
          REAL(dp), pointer, DIMENSION(:,:,:)::FI
          REAL(dp) :: e,Q1,Q2,Q3,c_star,v1,v2,v3,v4,v5,v11,v22,v33,v44,v55,S1,S2 &
    ,S3,a1,a2,a3,w1,w2,w3
          pll=pl_ex+1
          do ib=1,nbp

              select case (op)
                Case (1)
                  fi => dom(ib)%u
                  npp=dom(ib)%ieu+pll; nq=dom(ib)%jeu+pll; nr=dom(ib)%keu+pll
                Case (2)
                  fi => dom(ib)%v
                  npp=dom(ib)%iev+pll; nq=dom(ib)%jev+pll; nr=dom(ib)%kev+pll
                Case (3)
                  fi => dom(ib)%w
                  npp=dom(ib)%iew+pll; nq=dom(ib)%jew+pll; nr=dom(ib)%kew+pll
              end select

!
!.....Compute Dividierte Differenzen
!
              DO K=1,NR
                  DO I=1,npp-1
                      DO J=1,NQ
                          dom(ib)%D1(i,j,k) = (FI(i+1,j,k) - FI(i,j,k))/dom(ib)%dx
                      END DO
                  END DO
              END DO
!
!..... COMPUTE DERIVATIVE OF PHI (5-th Ord.)
!

              dom(ib)%dphi_dxplus = 0.0_dp
              dom(ib)%dphi_dxminus = 0.0_dp

              DO K=1,NR
                  DO I=1,npp-6
                      DO J=1,NQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          L=I

                          v1 = dom(ib)%D1(L+5,j,k)
                          v2 = dom(ib)%D1(L+4,j,k)
                          v3 = dom(ib)%D1(L+3,j,k)
                          v4 = dom(ib)%D1(L+2,j,k)
                          v5 = dom(ib)%D1(L+1,j,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0E-6_dp) * MAX(v11,v22,v33,v44,v55) + 1.0E-99_dp

                          S1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          S2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          S3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+S1)**2
                          a2 = 0.6_dp/(e+S2)**2
                          a3 = 0.3_dp/(e+S3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dxplus(i+3,j,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          L=I-1

                          v1 = dom(ib)%D1(L+1,j,k)
                          v2 = dom(ib)%D1(L+2,j,k)
                          v3 = dom(ib)%D1(L+3,j,k)
                          v4 = dom(ib)%D1(L+4,j,k)
                          v5 = dom(ib)%D1(L+5,j,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0E-6_dp) * MAX(v11,v22,v33,v44,v55) + 1.0E-99_dp

                          S1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          S2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          S3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+S1)**2
                          a2 = 0.6_dp/(e+S2)**2
                          a3 = 0.3_dp/(e+S3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dxminus(i+3,j,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      END DO
                  END DO
              END DO

          end do

          RETURN
      end subroutine HJ_WENO_dx
!######################################################################
      subroutine HJ_WENO_dy(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          INTEGER :: I,J,K,L,k_star,npp,nq,nr,op,ib,pll
          REAL(dp), pointer, DIMENSION(:,:,:)::FI
          REAL(dp) :: e,Q1,Q2,Q3,c_star,v1,v2,v3,v4,v5,v11,v22,v33,v44,v55,S1,S2 &
    ,S3,a1,a2,a3,w1,w2,w3
          pll=pl_ex+1
          do ib=1,nbp

              select case (op)
                Case (1)
                  fi => dom(ib)%u
                  npp=dom(ib)%ieu+pll; nq=dom(ib)%jeu+pll; nr=dom(ib)%keu+pll
                Case (2)
                  fi => dom(ib)%v
                  npp=dom(ib)%iev+pll; nq=dom(ib)%jev+pll; nr=dom(ib)%kev+pll
                Case (3)
                  fi => dom(ib)%w
                  npp=dom(ib)%iew+pll; nq=dom(ib)%jew+pll; nr=dom(ib)%kew+pll
              end select

!
!.....Compute Dividierte Differenzen
!
              DO K=1,NR
                  DO I=1,npp
                      DO J=1,NQ-1
                          dom(ib)%D1(i,j,k) = (FI(i,j+1,k) - FI(i,j,k))/dom(ib)%dy
                      END DO
                  END DO
              END DO
!
!..... COMPUTE DERIVATIVE OF PHI (5-th Ord.)
!

              dom(ib)%dphi_dyplus = 0.0_dp
              dom(ib)%dphi_dyminus = 0.0_dp

              DO K=1,NR
                  DO I=1,npp
                      DO J=1,NQ-6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          L=J

                          v1 = dom(ib)%D1(i,L+5,k)
                          v2 = dom(ib)%D1(i,L+4,k)
                          v3 = dom(ib)%D1(i,L+3,k)
                          v4 = dom(ib)%D1(i,L+2,k)
                          v5 = dom(ib)%D1(i,L+1,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0E-6_dp) * MAX(v11,v22,v33,v44,v55) + 1.0E-99_dp

                          S1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          S2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          S3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+S1)**2
                          a2 = 0.6_dp/(e+S2)**2
                          a3 = 0.3_dp/(e+S3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dyplus(i,j+3,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          L=J-1

                          v1 = dom(ib)%D1(i,L+1,k)
                          v2 = dom(ib)%D1(i,L+2,k)
                          v3 = dom(ib)%D1(i,L+3,k)
                          v4 = dom(ib)%D1(i,L+4,k)
                          v5 = dom(ib)%D1(i,L+5,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0E-6_dp) * MAX(v11,v22,v33,v44,v55) + 1.0E-99_dp

                          S1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          S2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          S3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+S1)**2
                          a2 = 0.6_dp/(e+S2)**2
                          a3 = 0.3_dp/(e+S3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dyminus(i,j+3,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      END DO
                  END DO
              END DO

          end do

          RETURN
      end subroutine HJ_WENO_dy
!######################################################################
      subroutine HJ_WENO_dz(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          INTEGER :: I,J,K,L,k_star,npp,nq,nr,op,ib,pll
          REAL(dp), pointer, DIMENSION(:,:,:)::FI
          REAL(dp) :: e,Q1,Q2,Q3,c_star,v1,v2,v3,v4,v5,v11,v22,v33,v44,v55,S1,S2 &
    ,S3,a1,a2,a3,w1,w2,w3
          pll=pl_ex+1
          do ib=1,nbp

              select case (op)
                Case (1)
                  fi => dom(ib)%u
                  npp=dom(ib)%ieu+pll; nq=dom(ib)%jeu+pll; nr=dom(ib)%keu+pll
                Case (2)
                  fi => dom(ib)%v
                  npp=dom(ib)%iev+pll; nq=dom(ib)%jev+pll; nr=dom(ib)%kev+pll
                Case (3)
                  fi => dom(ib)%w
                  npp=dom(ib)%iew+pll; nq=dom(ib)%jew+pll; nr=dom(ib)%kew+pll
              end select

!
!.....Compute Dividierte Differenzen
!
              DO K=1,NR-1
                  DO I=1,npp
                      DO J=1,NQ
                          dom(ib)%D1(i,j,k) = (FI(i,j,k+1) - FI(i,j,k))/dom(ib)%dz
                      END DO
                  END DO
              END DO
!
!..... COMPUTE DERIVATIVE OF PHI (5-th Ord.)
!

              dom(ib)%dphi_dzplus = 0.0_dp
              dom(ib)%dphi_dzminus = 0.0_dp

              DO K=1,NR-6
                  DO I=1,npp
                      DO J=1,NQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          L=K

                          v1 = dom(ib)%D1(i,j,L+5)
                          v2 = dom(ib)%D1(i,j,L+4)
                          v3 = dom(ib)%D1(i,j,L+3)
                          v4 = dom(ib)%D1(i,j,L+2)
                          v5 = dom(ib)%D1(i,j,L+1)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0E-6_dp) * MAX(v11,v22,v33,v44,v55) + 1.0E-99_dp

                          S1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          S2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          S3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+S1)**2
                          a2 = 0.6_dp/(e+S2)**2
                          a3 = 0.3_dp/(e+S3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dzplus(i,j,k+3) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          L=K-1

                          v1 = dom(ib)%D1(i,j,L+1)
                          v2 = dom(ib)%D1(i,j,L+2)
                          v3 = dom(ib)%D1(i,j,L+3)
                          v4 = dom(ib)%D1(i,j,L+4)
                          v5 = dom(ib)%D1(i,j,L+5)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0E-6_dp) * MAX(v11,v22,v33,v44,v55) + 1.0E-99_dp

                          S1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          S2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          S3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+S1)**2
                          a2 = 0.6_dp/(e+S2)**2
                          a3 = 0.3_dp/(e+S3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dzminus(i,j,k+3) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      END DO
                  END DO
              END DO

          end do

          RETURN
      end subroutine HJ_WENO_dz
!######################################################################
