      program main
      implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

       PetscInt is,js,iw,jw
       PetscInt one,three
       PetscErrorCode ierr
       KSP ksp
       DM dm
       external ComputeRHS,ComputeMatrix,ComputeInitialGuess

       one = 1
       three = 3

       call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
       call KSPCreate(MPI_COMM_WORLD,ksp,ierr)
       call DMDACreate2D(MPI_COMM_WORLD, DM_BOUNDARY_NONE,              &
     &    DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,-three,-three,            &
     &    PETSC_DECIDE,PETSC_DECIDE,one,one, PETSC_NULL_INTEGER,        &
     &    PETSC_NULL_INTEGER, dm, ierr)
       call KSPSetDM(ksp,dm,ierr)
       call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,             &
     &                                PETSC_NULL_OBJECT,ierr)
       call KSPSetComputeRHS(ksp,ComputeRHS,PETSC_NULL_OBJECT,ierr)
       call KSPSetComputeOperators(ksp,ComputeMatrix,                   &
     &      PETSC_NULL_OBJECT,ierr)
       call DMDAGetCorners(dm,is,js,PETSC_NULL_INTEGER,iw,jw,             &
     &                     PETSC_NULL_INTEGER,ierr)
       call KSPSetFromOptions(ksp,ierr)
       call KSPSetUp(ksp,ierr)
       call KSPSolve(ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
       call KSPDestroy(ksp,ierr)
       call DMDestroy(dm,ierr)
       call PetscFinalize(ierr)
       end


       subroutine ComputeInitialGuess(ksp,b,ctx,ierr)
       implicit none
       PetscErrorCode  ierr
       KSP ksp
       PetscInt ctx(*)
       Vec b
       PetscScalar  h

       h=0.0
       call VecSet(b,h,ierr)
       end subroutine

       subroutine ComputeRHS(ksp,b,dummy,ierr)
       implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
       PetscErrorCode  ierr
       KSP ksp
       Vec b
       integer dummy(*)
       PetscScalar  h,Hx,Hy
       PetscInt  mx,my
       DM dm

       call KSPGetDM(ksp,dm,ierr)
       call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,  &
     &                     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &                     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &                     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &                     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &                     PETSC_NULL_INTEGER,ierr)

       Hx = 1.0 / (mx-1)
       Hy = 1.0 / (my-1)
       h=Hx*Hy
       call VecSet(b,h,ierr)
       end subroutine

       subroutine ComputeMatrix(ksp,A,B,dummy,ierr)
       implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
       PetscErrorCode  ierr
       KSP ksp
       Mat A,B
       integer dummy(*)
       DM dm

      PetscInt    i,j,mx,my,xm
      PetscInt    ym,xs,ys,i1,i5
      PetscScalar  v(5),Hx,Hy
      PetscScalar  HxdHy,HydHx
      MatStencil   row(4),col(4,5)

      i1 = 1
      i5 = 5
      call KSPGetDM(ksp,dm,ierr)
      call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,  &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,             &
     &               PETSC_NULL_INTEGER,ierr)

      Hx = 1.0 / (mx-1)
      Hy = 1.0 / (my-1)
      HxdHy = Hx/Hy
      HydHx = Hy/Hx
      call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,            &
     &     PETSC_NULL_INTEGER,ierr)
      do 10,j=ys,ys+ym-1
        do 20,i=xs,xs+xm-1
          row(MatStencil_i) = i
          row(MatStencil_j) = j
          if (i.eq.0 .or. j.eq.0 .or. i.eq.mx-1 .or. j.eq.my-1 ) then
            v(1) = 2.0*(HxdHy + HydHx)
            call MatSetValuesStencil(B,i1,row,i1,row,v,                 &
     &           INSERT_VALUES,ierr)
          else
            v(1) = -HxdHy
            col(MatStencil_i,1) = i
            col(MatStencil_j,1) = j-1
            v(2) = -HydHx
            col(MatStencil_i,2) = i-1
            col(MatStencil_j,2) = j
            v(3) = 2.0*(HxdHy + HydHx)
            col(MatStencil_i,3) = i
            col(MatStencil_j,3) = j
            v(4) = -HydHx
            col(MatStencil_i,4) = i+1
            col(MatStencil_j,4) = j
            v(5) = -HxdHy
            col(MatStencil_i,5) = i
            col(MatStencil_j,5) = j+1
            call MatSetValuesStencil(B,i1,row,i5,col,v,                 &
     &           INSERT_VALUES,ierr)
            endif
 20      continue
 10   continue
       call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
       call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
       if ( A .ne. B) then
         call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
       endif
       end subroutine
