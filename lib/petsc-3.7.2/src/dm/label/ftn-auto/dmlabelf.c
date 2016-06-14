#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dmlabel.c */
/* Fortran interface file */

/*
* This file was generated automatically by bfort from the C source
* file.  
 */

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(PetscFortranAddr *)(a))
#define PetscFromPointer(a) (PetscFortranAddr)(a)
#define PetscRmPointer(a)
#endif

#include "petscdmlabel.h"
#include "petscis.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelhasvalue_ DMLABELHASVALUE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelhasvalue_ dmlabelhasvalue
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelhaspoint_ DMLABELHASPOINT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelhaspoint_ dmlabelhaspoint
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelstratumhaspoint_ DMLABELSTRATUMHASPOINT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelstratumhaspoint_ dmlabelstratumhaspoint
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelgetvalue_ DMLABELGETVALUE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelgetvalue_ dmlabelgetvalue
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelsetvalue_ DMLABELSETVALUE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelsetvalue_ dmlabelsetvalue
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelclearvalue_ DMLABELCLEARVALUE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelclearvalue_ dmlabelclearvalue
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelinsertis_ DMLABELINSERTIS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelinsertis_ dmlabelinsertis
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dmlabelgather_ DMLABELGATHER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dmlabelgather_ dmlabelgather
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void PETSC_STDCALL  dmlabelhasvalue_(DMLabel label,PetscInt *value,PetscBool *contains, int *__ierr ){
*__ierr = DMLabelHasValue(
	(DMLabel)PetscToPointer((label) ),*value,contains);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelhaspoint_(DMLabel label,PetscInt *point,PetscBool *contains, int *__ierr ){
*__ierr = DMLabelHasPoint(
	(DMLabel)PetscToPointer((label) ),*point,contains);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelstratumhaspoint_(DMLabel label,PetscInt *value,PetscInt *point,PetscBool *contains, int *__ierr ){
*__ierr = DMLabelStratumHasPoint(
	(DMLabel)PetscToPointer((label) ),*value,*point,contains);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelgetvalue_(DMLabel label,PetscInt *point,PetscInt *value, int *__ierr ){
*__ierr = DMLabelGetValue(
	(DMLabel)PetscToPointer((label) ),*point,value);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelsetvalue_(DMLabel label,PetscInt *point,PetscInt *value, int *__ierr ){
*__ierr = DMLabelSetValue(
	(DMLabel)PetscToPointer((label) ),*point,*value);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelclearvalue_(DMLabel label,PetscInt *point,PetscInt *value, int *__ierr ){
*__ierr = DMLabelClearValue(
	(DMLabel)PetscToPointer((label) ),*point,*value);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelinsertis_(DMLabel label,IS is,PetscInt *value, int *__ierr ){
*__ierr = DMLabelInsertIS(
	(DMLabel)PetscToPointer((label) ),
	(IS)PetscToPointer((is) ),*value);
}
PETSC_EXTERN void PETSC_STDCALL  dmlabelgather_(DMLabel label,PetscSF sf,DMLabel *labelNew, int *__ierr ){
*__ierr = DMLabelGather(
	(DMLabel)PetscToPointer((label) ),
	(PetscSF)PetscToPointer((sf) ),labelNew);
}
#if defined(__cplusplus)
}
#endif
