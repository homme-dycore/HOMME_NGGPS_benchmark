///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereGrid.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _CUBEDSPHEREGRID_H_
#define _CUBEDSPHEREGRID_H_

///////////////////////////////////////////////////////////////////////////////

#include "netcdfcpp.h"

#include "Preferences.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

#include "CubedSphereTrans.h"

///////////////////////////////////////////////////////////////////////////////

class CubedSphereGrid {

public:
	///	<summary>
	///		Calculate reconstruction coefficients from given data.
	///	</summary>
	virtual void ReconstructHighOrder(
		const DataVector<double> & dU
	) = 0;

	///	<summary>
	///		Calculate reconstruction coefficients from given data.
	///	</summary>
	virtual void ReconstructLowOrderMonotone(
		const DataVector<double> & dU
	) = 0;

	///	<summary>
	///		Evaluate the reconstruction at the given point.
	///	</summary>
	virtual void EvaluateReconstruction(
		int iP,
		int i,
		int j,
		double dX,
		double dY
	) {
	}

	///	<summary>
	///		Compute some data statistics from a pointwise data object.
	///	</summary>
	virtual void Checksum(
		const DataVector<double> & dU,
		double & dSum,
		double & dMin,
		double & dMax
	) const {
	}

public:
	///	<summary>
	///		Determine the element index from the given gnomonic coordinates.
	///	</summary>
	void IndexFromXY(
		double dX,
		double dY,
		int & iX,
		int & iY
	) const;

public:
	///	<summary>
	///		Accessors.
	///	</summary>
	inline int GetResolution() const {
		return m_nResolution;
	}

	inline int GetGhostElements() const {
		return m_nGhostElements;
	}

	inline int GetTotalElements() const {
		return m_nTotalElements;
	}

	inline double GetDeltaA() const {
		return m_dDeltaA;
	}

	inline const DataMatrix<double>& GetElementArea() const {
		return m_matElementArea;
	}

	inline const DataMatrix<double>& GetEdgeLength() const {
		return m_dEdgeLength;
	}

public:
	inline const DataMatrix3D<double>& GetGeoCentroidX() const {
		return m_matGeoCentroidX;
	}

public:
	inline const DataVector<double>& GetCentroidA() const {
		return m_dCentroidA;
	}

	inline const DataVector<double>& GetEdgeA() const {
		return m_dEdgeA;
	}

	inline const DataVector<double>& GetCentroidX() const {
		return m_dCentroidX;
	}

	inline const DataVector<double>& GetEdgeX() const {
		return m_dEdgeX;
	}

public:
	///	<summary>
	///		Combined accessors.
	///	</summary>
	inline int GetInteriorBegin() const {
		return 0;
	}

	inline int GetInteriorEnd() const {
		return m_nResolution;
	}

public:
	///	<summary>
	///		Get the reconstruction coefficients.
	///	</summary>
	inline double GetReconstruction(
		int iP,
		int iX,
		int iY,
		int iCoeff
	) const {
		return m_dReconstruction[iP][iX][iY][iCoeff];
	}

	///	<summary>
	///		Get an array of reconstruction coefficients.
	///	</summary>
	inline const double * GetReconstruction(
		int iP,
		int iX,
		int iY
	) const {
		return m_dReconstruction[iP][iX][iY];
	}

protected:
	///	<summary>
	///		Number of elements in each direction.
	///	</summary>
	int m_nResolution;

	///	<summary>
	///		Number of ghost elements on each side.
	///	</summary>
	int m_nGhostElements;

	///	<summary>
	///		Total number of elements in each direction.
	///	</summary>
	int m_nTotalElements;

	///	<summary>
	///		Horizontal and vertical edge lengths in terms of angle.
	///	</summary>
	double m_dDeltaA;

	///	<summary>
	///		Calculated per-element area.
	///	</summary>
	DataMatrix<double> m_matElementArea;

	///	<summary>
	///		Calculated line segment length.
	///	</summary>
	DataMatrix<double> m_dEdgeLength;

protected:
	///	<summary>
	///		Calculated element geometric centroids.
	///	</summary>
	DataMatrix3D<double> m_matGeoCentroidX;

protected:
	///	<summary>
	///		Calculated element centroids, stored as a vector of angles.
	///	</summary>
	DataVector<double> m_dCentroidA;

	///	<summary>
	///		Calculated edge coordinates, stored as a vector of angles.
	///	</summary>
	DataVector<double> m_dEdgeA;

	///	<summary>
	///		Calculated element centroids, stored as a vector of gnomonic
	///		X coordinates.
	///	</summary>
	DataVector<double> m_dCentroidX;

	///	<summary>
	///		Calculated edge coordinates, stored as a vector of gnomonic
	///		X coordinates.
	///	</summary>
	DataVector<double> m_dEdgeX;

protected:
	///	<summary>
	///		Reconstruction coefficients.
	///	</summary>
	DataMatrix4D<double> m_dReconstruction;

};

///////////////////////////////////////////////////////////////////////////////

#endif

