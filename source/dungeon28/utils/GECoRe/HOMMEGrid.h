///////////////////////////////////////////////////////////////////////////////
///
///	\file    HOMMEGrid.h
///	\author  Paul Ullrich
///	\version August 13, 2010
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

#ifndef _HOMMECONVERTER_H_
#define _HOMMECONVERTER_H_

///////////////////////////////////////////////////////////////////////////////

#include "CubedSphereGrid.h"

#include "netcdfcpp.h"

#include "Preferences.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "DataMatrix4D.h"

#include "CubedSphereTrans.h"

///////////////////////////////////////////////////////////////////////////////

class HOMMEGrid : public CubedSphereGrid {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		HOMMEGrid(
			NcFile & ncdf_file,
			int nNp
		);

	protected:
		///	<summary>
		///		Generate edge positions.
		///	</summary>
		void GenerateEdgePositions();

		///	<summary>
		///		Generate corresponding node list.
		///	</summary>
		void GenerateNodeList();

		///	<summary>
		///		Generate reconstruction.
		///	</summary>
		void GenerateReconstructionMatrix();

	public:
		///	<summary>
		///		Compute some data statistics from a pointwise data object.
		///	</summary>
		virtual void Checksum(
			const DataVector<double> & dU,
			double & dSum,
			double & dMin,
			double & dMax
		) const;

		///	<summary>
		///		Generate a reconstruction based on pointwise data.
		///	</summary>
		virtual void ReconstructHighOrder(
			const DataVector<double> & dU
		);

		///	<summary>
		///		Generate a monotone reconstruction based on pointwise data.
		///	</summary>
		virtual void ReconstructLowOrderMonotone(
			const DataVector<double> & dU
		);

		///	<summary>
		///		Evaluate the reconstruction at the given point.
		///	</summary>
		void EvaluateReconstruction(
			int iP,
			int i,
			int j,
			double dA,
			double dB
		) {
		}

	protected:
		///	<summary>
		///		Degree of polynomial.
		///	</summary>
		int m_nNp;

		///	<summary>
		///		Number of elements along an edge (resolution).
		///	</summary>
		int m_nNe;

		///	<summary>
		///		Number of nodes along an edge.
		///	</summary>
		int m_nNodesPerEdge;

		///	</summary>
		///	<summary>
		///		Longitudes from NetCDF file.
		///	</summary>
		DataVector<double> m_dLon;

		///	<summary>
		///		Latitudes from NetCDF file.
		///	</summary>
		DataVector<double> m_dLat;

		///	<summary>
		///		Area of elements on the dual grid.
		///	</summary>
		DataVector<double> m_dArea;

		///	<summary>
		///		Element centroid positions in equiangular coordinates.
		///	</summary>
		DataVector<double> m_dElementCentroidA;

		///	<summary>
		///		Node list.
		///	</summary>
		DataMatrix3D<int> m_matNodeList;

		///	<summary>
		///		Reconstruction matrix.
		///	</summary>
		DataMatrix<double> m_dReconsMatrix;
};

///////////////////////////////////////////////////////////////////////////////

#endif
