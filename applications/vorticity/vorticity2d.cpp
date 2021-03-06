/*! Computes the vorticity field from the 2D velocity vector field.
 * \file vorticity2d.cpp
 */

#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <functional>

#include <petscsys.h>
#include <petscdmda.h>

#include "petibm-utilities/field.h"
#include "petibm-utilities/grid.h"
#include "petibm-utilities/misc.h"
#include "petibm-utilities/timestep.h"
#include "petibm-utilities/vorticity.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	std::string directory, outdir, gridpath;
	PetibmGrid grid, gridux, griduy, gridwz;
	PetibmGridCtx gridCtx, griduxCtx, griduyCtx;
	PetibmField ux, uy, wz;
	PetibmFieldCtx fieldCtx;
	PetibmTimeStepCtx stepCtx;
	DM da;
	const PetscInt *plx, *ply;
	PetscInt *lx, *ly;
	PetscInt M, N, m, n;
	DMBoundaryType bType_x, bType_y;
	PetscInt ite;
	PetscMPIInt rank;
	PetscBool found = PETSC_FALSE;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	// parse command-line options
	ierr = PetibmGetDirectory(&directory); CHKERRQ(ierr);
	ierr = PetibmTimeStepGetOptions(nullptr, &stepCtx); CHKERRQ(ierr);
	ierr = PetibmGridGetOptions(nullptr, &gridCtx); CHKERRQ(ierr);
	ierr = PetibmFieldGetOptions(nullptr, &fieldCtx); CHKERRQ(ierr);
	{
		char dir[PETSC_MAX_PATH_LEN];
		ierr = PetscOptionsGetString(nullptr, nullptr, "-output_directory",
		                             dir, sizeof(dir), &found); CHKERRQ(ierr);
		outdir = (!found) ? directory : dir;
		mkdir((outdir).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	{
		char path[PETSC_MAX_PATH_LEN];
		ierr = PetscOptionsGetString(nullptr, nullptr, "-grid_path",
		                             path, sizeof(path), &found); CHKERRQ(ierr);
		gridpath = (!found) ? directory+"/grid.h5" : path;
	}

	// read cell-centered gridline stations
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.nx, &grid.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.ny, &grid.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		gridpath, "p", "x", grid.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		gridpath, "p", "y", grid.y.coords); CHKERRQ(ierr);
	// read staggered gridline stations for x-velocity
	griduxCtx.nx = (fieldCtx.periodic_x) ? gridCtx.nx : gridCtx.nx-1;
	griduxCtx.ny = gridCtx.ny;
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduxCtx.nx, &gridux.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduxCtx.ny, &gridux.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		gridpath, "u", "x", gridux.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		gridpath, "u", "y", gridux.y.coords); CHKERRQ(ierr);
	// read staggered gridline stations for y-velocity
	griduyCtx.nx = gridCtx.nx;
	griduyCtx.ny = (fieldCtx.periodic_y) ? gridCtx.ny : gridCtx.ny-1;
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduyCtx.nx, &griduy.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduyCtx.ny, &griduy.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		gridpath, "v", "x", griduy.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(
		gridpath, "v", "y", griduy.y.coords); CHKERRQ(ierr);
	// create grid for z-vorticity
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.nx-1, &gridwz.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.ny-1, &gridwz.y.coords); CHKERRQ(ierr);
	ierr = PetibmVorticityZComputeGrid(gridux, griduy, gridwz); CHKERRQ(ierr);
	if (rank == 0)
	{
		gridpath = outdir + "/grid.h5";
		ierr = PetibmGridHDF5Write(gridpath, "wz", gridwz); CHKERRQ(ierr);
	}
	// create base DMDA object
	bType_x = (fieldCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (fieldCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny,
	                    PETSC_DECIDE, PETSC_DECIDE, 1, 1, nullptr, nullptr,
	                    &da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);
	// get info from base DMDA for velocity components and z-vorticity
	ierr = DMDAGetOwnershipRanges(da, &plx, &ply, nullptr); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,
	                   nullptr,
	                   nullptr, nullptr, nullptr,
	                   &m, &n, nullptr,
	                   nullptr, nullptr,
	                   &bType_x, &bType_y, nullptr,
	                   nullptr); CHKERRQ(ierr);
	// create DMDA and vector for velocity in x-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	M = gridCtx.nx;
	N = gridCtx.ny;
	if (!fieldCtx.periodic_x)
	{
	  lx[m-1]--;
	  M--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    M, N, m, n, 1, 1, lx, ly,
	                    &ux.da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(ux.da); CHKERRQ(ierr);
	ierr = DMSetUp(ux.da);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ux.da, &ux.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(ux.da, &ux.local); CHKERRQ(ierr);
	// create DMDA and vector for velocity in y-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	M = gridCtx.nx;
	N = gridCtx.ny;
	if (!fieldCtx.periodic_y)
	{
	  ly[n-1]--;
	  N--;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_BOX,
	                    M, N, m, n, 1, 1, lx, ly,
	                    &uy.da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(uy.da); CHKERRQ(ierr);
	ierr = DMSetUp(uy.da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(uy.da, &uy.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(uy.da, &uy.local); CHKERRQ(ierr);
	// create DMDA and vector for z-vorticity
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	lx[m-1]--;
	ly[n-1]--;
	ierr = DMDACreate2d(PETSC_COMM_WORLD,
	                    bType_x, bType_y,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx-1, gridCtx.ny-1, m, n, 1, 1, lx, ly,
	                    &wz.da); CHKERRQ(ierr);
	ierr = DMSetFromOptions(wz.da); CHKERRQ(ierr);
	ierr = DMSetUp(wz.da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(wz.da, &wz.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(wz.da, &wz.local); CHKERRQ(ierr);

	// loop over the time steps to compute the z-vorticity
	for (ite=stepCtx.start; ite<=stepCtx.end; ite+=stepCtx.step)
	{
		ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %d]\n", ite); CHKERRQ(ierr);
		// get name of time-step file
		std::stringstream ss;
		ss << std::setfill('0') << std::setw(7) << ite << ".h5";
		std::string filename(ss.str());
		// read velocity field
		ierr = PetibmFieldHDF5Read(directory+"/"+filename, "u", ux); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Read(directory+"/"+filename, "v", uy); CHKERRQ(ierr);
		// compute the z-vorticity field
		ierr = PetibmVorticityZComputeField(
			gridux, griduy, ux, uy, wz); CHKERRQ(ierr);
		ierr = PetibmFieldHDF5Write(outdir+"/"+filename, "wz", wz); CHKERRQ(ierr);
	}

	ierr = PetibmGridDestroy(gridux); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(griduy); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridwz); CHKERRQ(ierr);
	ierr = DMDestroy(&da); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(ux); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uy); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(wz); CHKERRQ(ierr);
	
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
