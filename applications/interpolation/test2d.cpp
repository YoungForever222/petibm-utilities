/*! Tests 2D interpolation.
 * \file test2d.cpp
 */

#include <cmath>
#include <string>

#include <petscsys.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	DMBoundaryType bType_x,
	               bType_y;
	PetscReal xstart = 0.0, xend = 1.0,
	          ystart = 1.0, yend = 2.0,
	          h,
	          gammaA = 1.01, gammaB = 1.02;
	PetscInt i, j, ite, nt = 100;
	PetscBool found;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	ierr = PetscOptionsGetInt(
		nullptr, nullptr, "-nt", &nt, &found); CHKERRQ(ierr);

	PetibmField fieldA;
	PetibmGridCtx gridACtx;
	gridACtx.nx = 8;
	gridACtx.ny = 8;
	bType_x = (gridACtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (gridACtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    gridACtx.nx, gridACtx.ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &fieldA.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldA.da, nullptr, "-fieldA_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldA); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridACtx, fieldA.grid); CHKERRQ(ierr);
	ierr = VecSet(fieldA.global, 1.2345);
	ierr = DMGlobalToLocalBegin(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(
		fieldA.da, fieldA.global, INSERT_VALUES, fieldA.local); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(fieldA, 1.2345); CHKERRQ(ierr);

	PetscReal *xA, *yA;
	ierr = VecGetArray(fieldA.grid.x, &xA); CHKERRQ(ierr);
	ierr = VecGetArray(fieldA.grid.y, &yA); CHKERRQ(ierr);
	h = (xend - xstart) * (gammaA - 1.0) / (std::pow(gammaA, gridACtx.nx) - 1.0);
	xA[0] = xstart + 0.5 * h;
	for (i=0; i<gridACtx.nx-1; i++)
		xA[i+1] = xA[i] + h * std::pow(gammaA, i);
	h = (yend - ystart) * (gammaA - 1.0) / (std::pow(gammaA, gridACtx.ny) - 1.0);
	yA[0] = ystart + 0.5 * h;
	for (j=0; j<gridACtx.ny-1; j++)
		yA[j+1] = yA[j] + h * std::pow(gammaA, j);
	ierr = VecRestoreArray(fieldA.grid.x, &xA); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldA.grid.y, &yA); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write("gridA.h5", fieldA.grid); CHKERRQ(ierr);

	PetibmField fieldB;
	PetibmGridCtx gridBCtx;
	gridBCtx.nx = 4;
	gridBCtx.ny = 4;
	bType_x = (gridBCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (gridBCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    gridBCtx.nx, gridBCtx.ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &fieldB.da); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) fieldB.da, nullptr, "-fieldB_dmda_view"); CHKERRQ(ierr);
	ierr = PetibmFieldInitialize(fieldB); CHKERRQ(ierr);
	ierr = PetibmGridInitialize(gridBCtx, fieldB.grid); CHKERRQ(ierr);
	ierr = PetibmFieldExternalGhostPointsSet(fieldB, 1.2345); CHKERRQ(ierr);

	PetscReal *xB, *yB;
	ierr = VecGetArray(fieldB.grid.x, &xB); CHKERRQ(ierr);
	ierr = VecGetArray(fieldB.grid.y, &yB); CHKERRQ(ierr);
	h = (xend - xstart) * (gammaB - 1.0) / (std::pow(gammaB, gridBCtx.nx) - 1.0);
	xB[0] = xstart + 0.5 * h;
	for (i=0; i<gridBCtx.nx-1; i++)
		xB[i+1] = xB[i] + h * std::pow(gammaB, i);
	h = (yend - ystart) * (gammaB - 1.0) / (std::pow(gammaB, gridBCtx.ny) - 1.0);
	yB[0] = ystart + 0.5 * h;
	for (j=0; j<gridBCtx.ny-1; j++)
		yB[j+1] = yB[j] + h * std::pow(gammaB, j);
	ierr = VecRestoreArray(fieldB.grid.x, &xB); CHKERRQ(ierr);
	ierr = VecRestoreArray(fieldB.grid.y, &yB); CHKERRQ(ierr);
	ierr = PetibmGridHDF5Write("gridB.h5", fieldB.grid); CHKERRQ(ierr);

	for (ite=0; ite<nt; ite++)
	{
		if (ite % 4 == 0 or ite % 4 == 1)
		{
			ierr = PetibmFieldInterpolate2D(fieldA, fieldB, 1.2345); CHKERRQ(ierr);
		}
		else
		{
			ierr = PetibmFieldInterpolate2D(fieldB, fieldA, 1.2345); CHKERRQ(ierr);
		}
	}		

	ierr = PetibmFieldHDF5Write("fieldA.h5", "phi", fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldHDF5Write("fieldB.h5", "phi", fieldB); CHKERRQ(ierr);

	ierr = PetibmFieldDestroy(fieldA); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(fieldB); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main