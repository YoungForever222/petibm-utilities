/*! Converts the numerical solution from one format to another.
 * \file convert.cpp
 */

#include <string>

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"

#ifndef DIMENSIONS
#define DIMENSIONS 2
#endif


struct AppCtx
{
	char source[PETSC_MAX_PATH_LEN];
	char destination[PETSC_MAX_PATH_LEN];
	PetscBool hdf52binary = PETSC_FALSE;
}; // AppCtx


PetscErrorCode AppGetOptions(const char prefix[], AppCtx *ctx)
{
	PetscErrorCode ierr;
	char path[PETSC_MAX_PATH_LEN];
	PetscBool found;

	PetscFunctionBeginUser;

	// get path of configuration file
	ierr = PetscOptionsGetString(nullptr, prefix, "-config_file",
	                             path, sizeof(path), &found); CHKERRQ(ierr);
	if (found)
	{
		ierr = PetscOptionsInsertFile(
			PETSC_COMM_WORLD, nullptr, path, PETSC_FALSE); CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetString(nullptr, prefix, "-source", ctx->source,
	                             sizeof(ctx->source), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetString(
		nullptr, prefix, "-destination", ctx->destination,
		sizeof(ctx->destination), &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, prefix, "-hdf52binary", &ctx->hdf52binary, &found); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // AppGetOptions


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	PetscInt dim = DIMENSIONS;
	DM da;
	PetibmGridCtx gridCtx;
	PetibmFieldCtx fieldCtx;
	AppCtx appCtx;
	PetibmField field;
	DMBoundaryType bType_x, bType_y, bType_z;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	// parse command-line options
	ierr = PetibmGridGetOptions(nullptr, &gridCtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions(nullptr, &fieldCtx); CHKERRQ(ierr);
	ierr = AppGetOptions(nullptr, &appCtx); CHKERRQ(ierr);

	// create DMDA object
	bType_x = (fieldCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (fieldCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (fieldCtx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	if (dim == 2)
	{
		ierr = DMDACreate2d(PETSC_COMM_WORLD,
		                    bType_x, bType_y,
		                    DMDA_STENCIL_STAR,
		                    gridCtx.nx, gridCtx.ny,
		                    PETSC_DECIDE, PETSC_DECIDE,
		                    1, 1, nullptr, nullptr,
		                    &da); CHKERRQ(ierr);
	}
	else if (dim == 3)
	{
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
		                    bType_x, bType_y, bType_z,
		                    DMDA_STENCIL_STAR,
		                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
		                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
		                    1, 1, nullptr, nullptr, nullptr,
		                    &da); CHKERRQ(ierr);
	}
	else
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
		        "Support only provided for 2D or 3D DMDA objects");
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);

	// initialize, read, and write
	ierr = PetibmFieldInitialize(da, field); CHKERRQ(ierr);
	if (appCtx.hdf52binary)
	{
		ierr = PetibmFieldHDF5Read(
			appCtx.source, fieldCtx.name, field); CHKERRQ(ierr);
		ierr = PetibmFieldBinaryWrite(
			appCtx.destination, fieldCtx.name, field); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetibmFieldBinaryRead(
			appCtx.source, fieldCtx.name, field); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Write(
			appCtx.destination, fieldCtx.name, field); CHKERRQ(ierr);
	}

	ierr = PetibmFieldDestroy(field); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
} // main
