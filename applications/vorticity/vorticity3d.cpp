/*! Computes the vorticity field from the 3D velocity vector field.
 * \file vorticity3d.cpp
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
	std::string directory, outdir, griddir, gridpath;
	PetibmGrid grid, gridux, griduy, griduz, gridwx, gridwz;
	PetibmGridCtx gridCtx, griduxCtx, griduyCtx, griduzCtx;
	PetibmField ux, uy, uz, wx, wz;
	PetibmFieldCtx fieldCtx;
	PetibmTimeStepCtx stepCtx;
	DM da;
	const PetscInt *plx, *ply, *plz;
	PetscInt *lx, *ly, *lz;
	PetscInt M, N, P, m, n, p;
	DMBoundaryType bType_x, bType_y, bType_z;
	PetscInt ite;
	PetscMPIInt rank;
	PetscBool found,
	          compute_wx=PETSC_FALSE,
	          compute_wz=PETSC_FALSE,
	          binary_format=PETSC_FALSE;

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
	}
	{
		char dir[PETSC_MAX_PATH_LEN];
		ierr = PetscOptionsGetString(nullptr, nullptr, "-grid_directory",
		                             dir, sizeof(dir), &found); CHKERRQ(ierr);
		griddir = (!found) ? directory+"/grids" : dir;
	}
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-compute_wx", &compute_wx, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-compute_wz", &compute_wz, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(
		nullptr, nullptr, "-binary_format", &binary_format, &found); CHKERRQ(ierr);

	// read cell-centered gridline stations
	grid.dim = 3;
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.nx, &grid.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.ny, &grid.y.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, gridCtx.nz, &grid.z.coords); CHKERRQ(ierr);
	gridpath = griddir+"/cell-centered.h5";
	ierr = PetibmGridlineHDF5Read(gridpath, "x", grid.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "y", grid.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "z", grid.z.coords); CHKERRQ(ierr);
	// read staggered gridline stations for x-velocity
	griduxCtx.nx = (fieldCtx.periodic_x) ? gridCtx.nx : gridCtx.nx-1;
	griduxCtx.ny = gridCtx.ny;
	griduxCtx.nz = gridCtx.nz;
	gridux.dim = 3;
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduxCtx.nx, &gridux.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduxCtx.ny, &gridux.y.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduxCtx.nz, &gridux.z.coords); CHKERRQ(ierr);
	gridpath = griddir+"/staggered-x.h5";
	ierr = PetibmGridlineHDF5Read(gridpath, "x", gridux.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "y", gridux.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "z", gridux.z.coords); CHKERRQ(ierr);
	// read staggered gridline stations for y-velocity
	griduyCtx.nx = gridCtx.nx;
	griduyCtx.ny = (fieldCtx.periodic_y) ? gridCtx.ny : gridCtx.ny-1;
	griduyCtx.nz = gridCtx.nz;
	griduy.dim = 3;
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduyCtx.nx, &griduy.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduyCtx.ny, &griduy.y.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduyCtx.nz, &griduy.z.coords); CHKERRQ(ierr);
	gridpath = griddir+"/staggered-y.h5";
	ierr = PetibmGridlineHDF5Read(gridpath, "x", griduy.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "y", griduy.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "z", griduy.z.coords); CHKERRQ(ierr);
	// read staggered gridline stations for z-velocity
	griduzCtx.nx = gridCtx.nx;
	griduzCtx.ny = gridCtx.ny;
	griduzCtx.nz = (fieldCtx.periodic_z) ? gridCtx.nz : gridCtx.nz-1;
	griduz.dim = 3;
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduzCtx.nx, &griduz.x.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduzCtx.ny, &griduz.y.coords); CHKERRQ(ierr);
	ierr = VecCreateSeq(
		PETSC_COMM_SELF, griduzCtx.nz, &griduz.z.coords); CHKERRQ(ierr);
	gridpath = griddir+"/staggered-z.h5";
	ierr = PetibmGridlineHDF5Read(gridpath, "x", griduz.x.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "y", griduz.y.coords); CHKERRQ(ierr);
	ierr = PetibmGridlineHDF5Read(gridpath, "z", griduz.z.coords); CHKERRQ(ierr);
	// create grid for x-vorticity
	if (compute_wx)
	{
		gridwx.dim = 3;
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, gridCtx.nx, &gridwx.x.coords); CHKERRQ(ierr);
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, gridCtx.ny-1, &gridwx.y.coords); CHKERRQ(ierr);
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, gridCtx.nz-1, &gridwx.z.coords); CHKERRQ(ierr);
		ierr = PetibmVorticityXComputeGrid(griduy, griduz, gridwx); CHKERRQ(ierr);
		if (rank == 0)
		{
			gridpath = outdir+"/grids/wx.h5";
			mkdir((outdir+"/grids").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			ierr = PetibmGridHDF5Write(gridpath, gridwx); CHKERRQ(ierr);
		}
	}
	// create grid for z-vorticity
	if (compute_wz)
	{
		gridwz.dim = 3;
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, gridCtx.nx-1, &gridwz.x.coords); CHKERRQ(ierr);
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, gridCtx.ny-1, &gridwz.y.coords); CHKERRQ(ierr);
		ierr = VecCreateSeq(
			PETSC_COMM_SELF, gridCtx.nz, &gridwz.z.coords); CHKERRQ(ierr);
		ierr = PetibmVorticityZComputeGrid(gridux, griduy, gridwz); CHKERRQ(ierr);
		if (rank == 0)
		{
			gridpath = outdir+"/grids/wz.h5";
			mkdir((outdir+"/grids").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			ierr = PetibmGridHDF5Write(gridpath, gridwz); CHKERRQ(ierr);
		}
	}
	// create base DMDA object
	bType_x = (fieldCtx.periodic_x) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_y = (fieldCtx.periodic_y) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	bType_z = (fieldCtx.periodic_z) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_STAR,
	                    gridCtx.nx, gridCtx.ny, gridCtx.nz,
	                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	                    1, 1, nullptr, nullptr, nullptr,
	                    &da); CHKERRQ(ierr);
	// get info from base DMDA for velocity components and z-vorticity
	ierr = DMDAGetOwnershipRanges(da, &plx, &ply, &plz); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,
	                   nullptr,
	                   nullptr, nullptr, nullptr,
	                   &m, &n, &p,
	                   nullptr, nullptr,
	                   &bType_x, &bType_y, &bType_z,
	                   nullptr); CHKERRQ(ierr);
	// create DMDA and vector for velocity in x-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	M = gridCtx.nx;
	N = gridCtx.ny;
	P = gridCtx.nz;
	if (!fieldCtx.periodic_x)
	{
		lx[m-1]--;
		M--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    M, N, P, m, n, p, 1, 1, lx, ly, lz,
	                    &ux.da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ux.da, &ux.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(ux.da, &ux.local); CHKERRQ(ierr);
	// create DMDA and vector for velocity in y-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	M = gridCtx.nx;
	N = gridCtx.ny;
	P = gridCtx.nz;
	if (!fieldCtx.periodic_y)
	{
		ly[n-1]--;
		N--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    M, N, P, m, n, p, 1, 1, lx, ly, lz,
	                    &uy.da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(uy.da, &uy.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(uy.da, &uy.local); CHKERRQ(ierr);
	// create DMDA and vector for velocity in z-direction
	ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
	ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
	ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
	M = gridCtx.nx;
	N = gridCtx.ny;
	P = gridCtx.nz;
	if (!fieldCtx.periodic_z)
	{
		lz[p-1]--;
		P--;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
	                    bType_x, bType_y, bType_z,
	                    DMDA_STENCIL_BOX,
	                    M, N, P, m, n, p, 1, 1, lx, ly, lz,
	                    &uz.da); CHKERRQ(ierr);
	ierr = PetscFree(lx); CHKERRQ(ierr);
	ierr = PetscFree(ly); CHKERRQ(ierr);
	ierr = PetscFree(lz); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(uz.da, &uz.global); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(uz.da, &uz.local); CHKERRQ(ierr);
	// create DMDA and vector for x-vorticity
	if (compute_wx)
	{
		ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
		ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
		ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
		ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
		ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
		ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
		ly[n-1]--;
		lz[p-1]--;
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
		                    bType_x, bType_y, bType_z,
		                    DMDA_STENCIL_STAR,
		                    gridCtx.nx, gridCtx.ny-1, gridCtx.nz-1, m, n, p,
		                    1, 1, lx, ly, lz,
		                    &wx.da); CHKERRQ(ierr);
		ierr = PetscFree(lx); CHKERRQ(ierr);
		ierr = PetscFree(ly); CHKERRQ(ierr);
		ierr = PetscFree(lz); CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(wx.da, &wx.global); CHKERRQ(ierr);
		ierr = DMCreateLocalVector(wx.da, &wx.local); CHKERRQ(ierr);
	}
	// create DMDA and vector for z-vorticity
	if (compute_wz)
	{
		ierr = PetscMalloc(m*sizeof(*lx), &lx); CHKERRQ(ierr);
		ierr = PetscMalloc(n*sizeof(*ly), &ly); CHKERRQ(ierr);
		ierr = PetscMalloc(p*sizeof(*lz), &lz); CHKERRQ(ierr);
		ierr = PetscMemcpy(lx, plx, m*sizeof(*lx)); CHKERRQ(ierr);
		ierr = PetscMemcpy(ly, ply, n*sizeof(*ly)); CHKERRQ(ierr);
		ierr = PetscMemcpy(lz, plz, p*sizeof(*lz)); CHKERRQ(ierr);
		lx[m-1]--;
		ly[n-1]--;
		ierr = DMDACreate3d(PETSC_COMM_WORLD,
		                    bType_x, bType_y, bType_z,
		                    DMDA_STENCIL_STAR,
		                    gridCtx.nx-1, gridCtx.ny-1, gridCtx.nz, m, n, p,
		                    1, 1, lx, ly, lz,
		                    &wz.da); CHKERRQ(ierr);
		ierr = PetscFree(lx); CHKERRQ(ierr);
		ierr = PetscFree(ly); CHKERRQ(ierr);
		ierr = PetscFree(lz); CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(wz.da, &wz.global); CHKERRQ(ierr);
		ierr = DMCreateLocalVector(wz.da, &wz.local); CHKERRQ(ierr);
	}

	std::string extension = ".h5";
	PetscViewerType viewerType = PETSCVIEWERHDF5;
	if (binary_format == PETSC_TRUE)
	{
		extension = ".dat";
		viewerType = PETSCVIEWERBINARY;
	}

	// loop over the time steps to compute the z-vorticity
	for (ite=stepCtx.start; ite<=stepCtx.end; ite+=stepCtx.step)
	{
		ierr = PetscPrintf(
			PETSC_COMM_WORLD, "[time-step %d]\n", ite); CHKERRQ(ierr);
		// get time-step directory
		std::stringstream ss;
		ss << directory << "/" << std::setfill('0') << std::setw(7) << ite;
		std::string folder(ss.str());
		// read velocity field
		ierr = PetibmFieldRead(
			folder+"/uy"+extension, "uy", viewerType, uy); CHKERRQ(ierr);
		if (compute_wz)
		{
			ierr = PetibmFieldRead(
				folder+"/ux"+extension, "ux", viewerType, ux); CHKERRQ(ierr);
		}
		if (compute_wx)
		{
			ierr = PetibmFieldRead(
				folder+"/uz"+extension, "uz", viewerType, uz); CHKERRQ(ierr);
		}
		// get time-step directory to save
		std::stringstream ssout;
		ssout << outdir << "/" << std::setfill('0') << std::setw(7) << ite;
		std::string outfolder(ssout.str());
		mkdir(outfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		// compute the x-vorticity field
		if (compute_wx)
		{
			ierr = PetibmVorticityXComputeField(
				griduy, griduz, uy, uz, wx); CHKERRQ(ierr);
			ierr = PetibmFieldHDF5Write(outfolder+"/wx.h5", "wx", wx); CHKERRQ(ierr);
		}
		// compute the z-vorticity field
		if (compute_wz)
		{
			ierr = PetibmVorticityZComputeField(
				gridux, griduy, ux, uy, wz); CHKERRQ(ierr);
			ierr = PetibmFieldHDF5Write(outfolder+"/wz.h5", "wz", wz); CHKERRQ(ierr);
		}
	}

	ierr = PetibmGridDestroy(gridux); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(griduy); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(griduz); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridwx); CHKERRQ(ierr);
	ierr = PetibmGridDestroy(gridwz); CHKERRQ(ierr);
	ierr = DMDestroy(&da); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(ux); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uy); CHKERRQ(ierr);
	ierr = PetibmFieldDestroy(uz); CHKERRQ(ierr);
	if (compute_wx)
	{
		ierr = PetibmFieldDestroy(wx); CHKERRQ(ierr);
	}
	if (compute_wz)
	{
		ierr = PetibmFieldDestroy(wz); CHKERRQ(ierr);
	}
	
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
