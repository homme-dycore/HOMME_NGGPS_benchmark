///////////////////////////////////////////////////////////////////////////////
///
///	\file    HOMMEGrid.cpp
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

#include "HOMMEGrid.h"

#include "Announce.h"
#include "CubedSphereTrans.h"

#include "LinearAlgebra.h"
#include "GaussLobattoQuadrature.h"
#include "MathHelper.h"

#include <map>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

class LonLatIndex {

public:
	static const double TINY;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	LonLatIndex(
		double dLon,
		double dLat
	) {
		// Determine the branch of this longitude value
		int nBranch = static_cast<int>((dLon + TINY) / (360.0));

		// Remap the longitude into the range [0, 360 - TINY]
		m_dLon = dLon - 360.0 * static_cast<double>(nBranch);

		// Store latitude and index
		m_dLat = dLat;
	}

	///	<summary>
	///		Accessors.
	///	</summary>
	double GetLongitude() const {
		return m_dLon;
	}

	double GetLatitude() const {
		return m_dLat;
	}

public:
	///	<summary>
	///		Comparator.
	///	</summary>
	int operator< (const LonLatIndex & lonlatix) const {
		if (fabs(m_dLat - lonlatix.m_dLat) < TINY) {
			if (90.0 - fabs(m_dLat) < TINY) {
				return 0;
			}

			if (m_dLon < lonlatix.m_dLon + TINY) {
				return 0;

			} else {
				return 1;
			}

		} else if (m_dLat < lonlatix.m_dLat) {
			return 1;

		} else {
			return 0;
		}
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	int operator== (const LonLatIndex & lonlatix) const {
		if (((fabs(m_dLon - lonlatix.m_dLon) < TINY)
				|| (fabs(m_dLat) - 90.0 < TINY))
			&& (fabs(m_dLat - lonlatix.m_dLat) < TINY)
		) {
			return 1;
		}

		return 0;
	}

	///	<summary>
	///		Negative comparator.
	///	</summary>
	int operator!= (const LonLatIndex & lonlatix) const {
		return !(*this == lonlatix);
	}

protected:
	///	<summary>
	///		Longitude of this point.
	///	</summary>
	double m_dLon;

	///	<summary>
	///		Latitude of this point.
	///	</summary>
	double m_dLat;
};

///////////////////////////////////////////////////////////////////////////////

const double LonLatIndex::TINY = 1.0e-10;

///////////////////////////////////////////////////////////////////////////////

typedef std::map<LonLatIndex, int> LonLatIndexMap;
typedef std::pair<LonLatIndex, int> LonLatIndexPair;
typedef LonLatIndexMap::const_iterator LonLatIndexIterator;

///////////////////////////////////////////////////////////////////////////////

HOMMEGrid::HOMMEGrid(
	NcFile & ncdf_file,
	int nNp
) {
	// Get point count
	NcDim *ncDimCol = ncdf_file.get_dim("ncol");
	if (ncDimCol == NULL) {
		_EXCEPTIONT("NetCDF file has no dimension \"ncol\".");
	}
	long nCol = ncDimCol->size();

	// Load in array of latitudes
	m_dLat.Initialize(nCol);

	NcVar *ncLat = ncdf_file.get_var("lat");
	if (ncLat == NULL) {
		_EXCEPTIONT("Input NetCDF file has no variable \"lat\".");
	}
	ncLat->get(m_dLat, nCol);

	// Load in array of longitudes
	m_dLon.Initialize(nCol);

	NcVar *ncLon = ncdf_file.get_var("lon");
	if (ncLat == NULL) {
		_EXCEPTIONT("Input NetCDF file has no variable \"lon\".");
	}
	ncLon->get(m_dLon, nCol);

	// Load in areas
	m_dArea.Initialize(nCol);

	NcVar *ncArea = ncdf_file.get_var("area");
	if (ncArea == NULL) {
		_EXCEPTIONT("Input NetCDF file has no variable \"area\".");
	}
	ncArea->get(m_dArea, nCol);

	// Get the polynomial degree
	m_nNp = nNp;

	// Calculate number of elements along each face
	//   - obtained from m_nNe = sqrt(6*ncol-12) / (6 * (m_nNp - 1))
	m_nNe = ISqrt(6 * static_cast<int>(nCol) - 12);
	if (m_nNe * m_nNe != 6* static_cast<int>(nCol) - 12) {
		_EXCEPTIONT("Number of columns does not match a known resolution!");
	}
	if ((m_nNe % (m_nNp - 1)) != 0) {
		_EXCEPTIONT("Number of columns does not match a known resolution!");
	}
	m_nNe /= (m_nNp - 1);

	if ((m_nNe % 6) != 0) {
		_EXCEPTIONT("Number of columns does not match a known resolution!");
	}
	m_nNe /= 6;

	// Number of nodes per edge
	m_nNodesPerEdge = m_nNe * (m_nNp - 1) + 1;

	// Number of elements per edge
	m_nResolution = m_nNodesPerEdge - 1;

	// Output
	Announce("Detected resolution: %i", m_nNe);
	Announce("Nodes per edge: %i", m_nNodesPerEdge);

	// Generate edge positions
	GenerateEdgePositions();

	// Generate node list
	GenerateNodeList();

	// Generate the reconstruction matrix
	GenerateReconstructionMatrix();
}

///////////////////////////////////////////////////////////////////////////////

void HOMMEGrid::GenerateEdgePositions() {

	// Grid spacing per element
	double dDeltaA = 0.5 * M_PI / static_cast<double>(m_nNe);

	// Gauss-Lobatto quadrature points
	DataVector<double> dG;
	DataVector<double> dW;

	GaussLobattoQuadrature::GetPoints(m_nNp, 0.0, dDeltaA, dG, dW);

	// Initialize the array of edge positions
	m_dEdgeA.Initialize(m_nNodesPerEdge);
	m_dEdgeX.Initialize(m_nNodesPerEdge);

	for (int a = 0; a < m_nNe; a++) {
	for (int i = 0; i < m_nNp-1; i++) {
		int ix = a * (m_nNp-1) + i;

		m_dEdgeA[ix] = - M_PI / 4.0 + dDeltaA * static_cast<double>(a) + dG[i];
		m_dEdgeX[ix] = tan(m_dEdgeA[ix]);
	}
	}

	m_dEdgeA[m_nNodesPerEdge-1] = M_PI / 4.0;
	m_dEdgeX[m_nNodesPerEdge-1] = 1.0;

	// Calculate element centroids
	m_dElementCentroidA.Initialize(m_nNe);

	for (int i = 0; i < m_nNe; i++) {
		m_dElementCentroidA[i] =
			(static_cast<double>(i) + 0.5) * dDeltaA - 0.25 * M_PI;
	}

	// Sub-element centroids are simply given by element centroids
	m_dCentroidA.Initialize(m_nNodesPerEdge-1);
	m_dCentroidX.Initialize(m_nNodesPerEdge-1);

	for (int i = 0; i < m_nNodesPerEdge-1; i++) {
		m_dCentroidA[i] = m_dElementCentroidA[i / (m_nNp-1)];

		m_dCentroidX[i] = tan(m_dCentroidA[i]);
	}
}

///////////////////////////////////////////////////////////////////////////////

void HOMMEGrid::GenerateNodeList() {

	// Generate node indices
	m_matNodeList.Initialize(6, m_nNodesPerEdge, m_nNodesPerEdge);

	// Put all nodes into a latlon map
	LonLatIndexMap mapLonLatIx;
	for (int k = 0; k < m_dLon.GetRows(); k++) {
		mapLonLatIx.insert(
			LonLatIndexPair(
				LonLatIndex(m_dLon[k], m_dLat[k]), k));
	}

	// Tiny quantity
	const double TINY = LonLatIndex::TINY;

	// Loop over all nodes
	for (int iP = 0; iP < 6; iP++) {
	for (int i = 0; i < m_nNodesPerEdge; i++) {
	for (int j = 0; j < m_nNodesPerEdge; j++) {

		// Convert node positions to RLL
		double dPtLon;
		double dPtLat;

		CubedSphereTrans::RLLFromXYP(
			m_dEdgeX[i],
			m_dEdgeX[j],
			iP,
			dPtLon, dPtLat);

		dPtLat = dPtLat * 180.0 / M_PI;
		dPtLon = dPtLon * 180.0 / M_PI;

		LonLatIndex llix(dPtLon, dPtLat);

		LonLatIndexIterator iter = mapLonLatIx.find(llix);

		if (iter != mapLonLatIx.end()) {
			m_matNodeList[iP][i][j] = iter->second;

			if (llix != iter->first) {
				std::cout << "I: " << iP << ", " << i << ", " << j << std::endl;
				std::cout << "Pt: " << dPtLon << ", " << dPtLat << std::endl;
				std::cout << "iter: " << iter->first.GetLongitude()
					<< ", " << iter->first.GetLatitude() << std::endl;
				_EXCEPTIONT("Mismatch error.");
			}

		} else {
			_EXCEPTION5("\nUnable to find node CS(%i, %f, %f) RLL(%f, %f)",
				iP, m_dEdgeA[i], m_dEdgeA[j],
				dPtLon, dPtLat);
		}
	}
	}
	}

}

///////////////////////////////////////////////////////////////////////////////

void HOMMEGrid::GenerateReconstructionMatrix() {

	// Vandermonde matrix
	DataMatrix<double> dInvVandermonde;
	dInvVandermonde.Initialize(m_nNp * m_nNp, m_nNp * m_nNp);

	// Integration matrix
	DataMatrix<double> dIntegration;
	dIntegration.Initialize(m_nNp * m_nNp, m_nNp * m_nNp);

	// Buffers used for LAPACK calls
	DataVector<int> iPIV;
	iPIV.Initialize(m_nNp * m_nNp);

	DataVector<double> dWork;
	dWork.Initialize(m_nNp * m_nNp);

	// Element centroid
	double dElementAlpha0 = m_dElementCentroidA[0];
	double dElementBeta0 = m_dElementCentroidA[0];

	// Construct the Vandermonde matrix and invert
	dInvVandermonde.Zero();
	for (int i = 0; i < m_nNp; i++) {
	for (int j = 0; j < m_nNp; j++) {
		int iCV = i * m_nNp + j;

		double dDA = m_dEdgeA[i] - dElementAlpha0;
		double dDB = m_dEdgeA[j] - dElementBeta0;

		int iMom = 0;
		for (int p = 0; p < m_nNp; p++) {
		for (int q = 0; q < m_nNp; q++) {
			dInvVandermonde[iMom][iCV] = IPow(dDA, p) * IPow(dDB, q);

			iMom++;
		}
		}
	}
	}

	// Invert the Vandermonde matrix
	LAPACK::DGETRF(dInvVandermonde, iPIV);
	LAPACK::DGETRI(dInvVandermonde, iPIV, dWork);

	// Copy inverse Vandermonde matrix to reconstruction matrix
	m_dReconsMatrix = dInvVandermonde;
}

///////////////////////////////////////////////////////////////////////////////

void HOMMEGrid::Checksum(
	const DataVector<double> & dU,
	double & dSum,
	double & dMin,
	double & dMax
) const {

	// Verify data is the right size
	if (dU.GetRows() != m_dLon.GetRows()) {
		_EXCEPTION2("Invalid data size: Found %i, Expected %i",
			dU.GetRows(), m_dLon.GetRows());
	}

	// Reset data statistics
	dSum = 0.0;
	dMax = dU[0];
	dMin = dU[0];

	// Loop through all elements
	for (int i = 0; i < dU.GetRows(); i++) {
		dSum += dU[i] * m_dArea[i];

		if (dU[i] > dMax) {
			dMax = dU[i];
		}
		if (dU[i] < dMin) {
			dMin = dU[i];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void HOMMEGrid::ReconstructHighOrder(
	const DataVector<double> & dU
) {
	// Verify data is the right size
	if (dU.GetRows() != m_dLon.GetRows()) {
		_EXCEPTION2("Invalid data size: Found %i, Expected %i",
			dU.GetRows(), m_dLon.GetRows());
	}

	// Number of sub-elements along each edge
	int nNumSubElements = m_nNe * (m_nNp - 1);

	// Initialize reconstruction coefficients
	m_dReconstruction.Initialize(
		6, nNumSubElements, nNumSubElements, m_nNp * m_nNp);

	// Temporary pointwise values for this element
	DataVector<double> dPointValues;
	dPointValues.Initialize(m_nNp * m_nNp);

	// Reconstruction coefficients for the entire element
	DataVector<double> dElementReconsCoeffs;
	dElementReconsCoeffs.Initialize(m_nNp * m_nNp);

	// Loop through all elements
	for (int iP = 0; iP < 6; iP++) {
	for (int a = 0; a < m_nNe; a++) {
	for (int b = 0; b < m_nNe; b++) {

		int iAElement = a * (m_nNp - 1);
		int iBElement = b * (m_nNp - 1);

		// Copy pointwise data into point values vector
		for (int i = 0; i < m_nNp; i++) {
		for (int j = 0; j < m_nNp; j++) {

			int iA = iAElement + i;
			int iB = iBElement + j;

			// Define the metric at this point
			double dX = m_dEdgeX[iA];
			double dY = m_dEdgeX[iB];
			double dDelta = sqrt(1.0 + dX * dX + dY * dY);
			double dSqrtG = (1.0 + dX * dX) * (1.0 + dY * dY)
				/ (dDelta * dDelta * dDelta);

			int iNode = m_matNodeList[iP][iA][iB];

			double dValue = dU[iNode];

			// DEBUG:  Override file information
			//dValue = 100.0;
			/*
			if (iP == 0) {
				double dR = sqrt(
					m_dElementCentroidA[i] * m_dElementCentroidA[i]
					+ m_dElementCentroidA[j] * m_dElementCentroidA[j]);

				dValue = 100.0 * dR;
			} else {
				dValue = 0.0;
			}
			*/

			// Multiply by metric weight
			dPointValues[i * m_nNp + j] = dValue * dSqrtG;
		}
		}

		// Multiply to obtain elementwise reconstruction coefficients
		dElementReconsCoeffs.Zero();

		for (int m = 0; m < m_nNp * m_nNp; m++) {
		for (int n = 0; n < m_nNp * m_nNp; n++) {
			dElementReconsCoeffs[m] +=
				m_dReconsMatrix[n][m] * dPointValues[n];
		}
		}
/*
		// DEBUG: FLATTEN RECONSTRUCTION
        for (int k = 1; k < 16; k++) {
            dElementReconsCoeffs[k] = 0.0;
        }

        double dW[4];
        dW[0] = 0.166666666666667;
        dW[1] = 0.833333333333333;
        dW[2] = 0.833333333333333;
        dW[3] = 0.166666666666667;

		dElementReconsCoeffs[0] = 0.0;
        for (int n = 0; n < m_nNp * m_nNp; n++) {
            dElementReconsCoeffs[0] +=
                dW[n/m_nNp] * dW[n%m_nNp] * dPointValues[n] / 4.0;
        }
*/
/*
		// Element centroid
		double dA0 = m_dElementCentroidA[i];
		double dB0 = m_dElementCentroidA[j];

		//printf("%1.5e %1.5e\n\n", dA0, dB0);

		// Shift elementwise reconstruction coefficients to subelements
		// Constant term
		m_dReconstruction[iP][iE][jE][0] =
			dElementReconsCoeffs[0]
			- dA0 * dElementReconsCoeffs[4]
			- dB0 * dElementReconsCoeffs[1]
			+ dA0 * dA0 * dElementReconsCoeffs[8]
			+ dA0 * dB0 * dElementReconsCoeffs[5]
			+ dB0 * dB0 * dElementReconsCoeffs[2]
			- dA0 * dA0 * dA0 * dElementReconsCoeffs[12]
			- dA0 * dA0 * dB0 * dElementReconsCoeffs[9]
			- dA0 * dB0 * dB0 * dElementReconsCoeffs[6]
			- dB0 * dB0 * dB0 * dElementReconsCoeffs[3];

		m_dReconstruction[iP][iE][jE][0] +=
			+ dA0 * dB0 * dB0 * dB0 * dElementReconsCoeffs[7]
			+ dA0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[10]
			+ dA0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[13]
			- dA0 * dA0 * dB0 * dB0 * dB0 * dElementReconsCoeffs[11]
			- dA0 * dA0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[14]
			+ dA0 * dA0 * dA0 * dB0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Beta derivative
		m_dReconstruction[iP][iE][jE][1] =
			dElementReconsCoeffs[1]
			- 2.0 * dB0 * dElementReconsCoeffs[2]
			- dA0 * dElementReconsCoeffs[5]
			+ 3.0 * dB0 * dB0 * dElementReconsCoeffs[3]
			+ 2.0 * dA0 * dB0 * dElementReconsCoeffs[6]
			+ dA0 * dA0 * dElementReconsCoeffs[9]
			- 3.0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[7]
			- 2.0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[10]
			- dA0 * dA0 * dA0 * dElementReconsCoeffs[13]
			+ 3.0 * dA0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[11]
			+ 2.0 * dA0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[14]
			- 3.0 * dA0 * dA0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Alpha derivative
		m_dReconstruction[iP][iE][jE][4] =
			dElementReconsCoeffs[4]
			- dB0 * dElementReconsCoeffs[5]
			- 2.0 * dA0 * dElementReconsCoeffs[8]
			+ dB0 * dB0 * dElementReconsCoeffs[6]
			+ 2.0 * dA0 * dB0 * dElementReconsCoeffs[9]
			+ 3.0 * dA0 * dA0 * dElementReconsCoeffs[12]
			- dB0 * dB0 * dB0 * dElementReconsCoeffs[7]
			- 2.0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[10]
			- 3.0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[13]
			+ 2.0 * dA0 * dB0 * dB0 * dB0 * dElementReconsCoeffs[11]
			+ 3.0 * dA0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[14]
			- 3.0 * dA0 * dA0 * dB0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Beta, Beta derivative
		m_dReconstruction[iP][iE][jE][2] =
			dElementReconsCoeffs[2]
			- 3.0 * dB0 * dElementReconsCoeffs[3]
			- dA0 * dElementReconsCoeffs[6]
			+ 3.0 * dA0 * dB0 * dElementReconsCoeffs[7]
			+ dA0 * dA0 * dElementReconsCoeffs[10]
			- 3.0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[11]
			- dA0 * dA0 * dA0 * dElementReconsCoeffs[14]
			+ 3.0 * dA0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[15];

		// Alpha, Beta derivative
		m_dReconstruction[iP][iE][jE][5] =
			dElementReconsCoeffs[5]
			- 2.0 * dB0 * dElementReconsCoeffs[6]
			- 2.0 * dA0 * dElementReconsCoeffs[9]
			+ 3.0 * dB0 * dB0 * dElementReconsCoeffs[7]
			+ 4.0 * dA0 * dB0 * dElementReconsCoeffs[10]
			+ 3.0 * dA0 * dA0 * dElementReconsCoeffs[13]
			- 6.0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[11]
			- 6.0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[14]
			+ 9.0 * dA0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Alpha, Alpha derivative
		m_dReconstruction[iP][iE][jE][8] =
			dElementReconsCoeffs[8]
			- dB0 * dElementReconsCoeffs[9]
			- 3.0 * dA0 * dElementReconsCoeffs[12]
			+ dB0 * dB0 * dElementReconsCoeffs[10]
			+ 3.0 * dA0 * dB0 * dElementReconsCoeffs[13]
			- dB0 * dB0 * dB0 * dElementReconsCoeffs[11]
			- 3.0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[14]
			+ 3.0 * dA0 * dB0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Beta, Beta, Beta derivative
		m_dReconstruction[iP][iE][jE][3] =
			dElementReconsCoeffs[3]
			- dA0 * dElementReconsCoeffs[7]
			+ dA0 * dA0 * dElementReconsCoeffs[11]
			- dA0 * dA0 * dA0 * dElementReconsCoeffs[15];

		// Alpha, Beta, Beta derivative
		m_dReconstruction[iP][iE][jE][6] =
			dElementReconsCoeffs[6]
			- 3.0 * dB0 * dElementReconsCoeffs[7]
			- 2.0 * dA0 * dElementReconsCoeffs[10]
			+ 6.0 * dA0 * dB0 * dElementReconsCoeffs[11]
			+ 3.0 * dA0 * dA0 * dElementReconsCoeffs[14]
			- 9.0 * dA0 * dA0 * dB0 * dElementReconsCoeffs[15];

		// Alpha, Alpha, Beta derivative
		m_dReconstruction[iP][iE][jE][9] =
			dElementReconsCoeffs[9]
			- 2.0 * dB0 * dElementReconsCoeffs[10]
			- 3.0 * dA0 * dElementReconsCoeffs[13]
			+ 3.0 * dB0 * dB0 * dElementReconsCoeffs[11]
			+ 6.0 * dA0 * dB0 * dElementReconsCoeffs[14]
			- 9.0 * dA0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Alpha, Alpha, Alpha derivative
		m_dReconstruction[iP][iE][jE][12] =
			dElementReconsCoeffs[12]
			- dB0 * dElementReconsCoeffs[13]
			+ dB0 * dB0 * dElementReconsCoeffs[14]
			- dB0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Beta, Beta, Beta, Alpha derivative
		m_dReconstruction[iP][iE][jE][7] =
			dElementReconsCoeffs[7]
			- 2.0 * dA0 * dElementReconsCoeffs[11]
			+ 3.0 * dA0 * dA0 * dElementReconsCoeffs[15];

		// Beta, Beta, Alpha, Alpha derivative
		m_dReconstruction[iP][iE][jE][10] =
			dElementReconsCoeffs[10]
			- 3.0 * dB0 * dElementReconsCoeffs[11]
			- 3.0 * dA0 * dElementReconsCoeffs[14]
			+ 9.0 * dA0 * dB0 * dElementReconsCoeffs[15];

		// Beta, Alpha, Alpha, Alpha derivative
		m_dReconstruction[iP][iE][jE][13] =
			dElementReconsCoeffs[13]
			- 2.0 * dB0 * dElementReconsCoeffs[14]
			+ 3.0 * dB0 * dB0 * dElementReconsCoeffs[15];

		// Beta, Beta, Beta, Alpha, Alpha derivative
		m_dReconstruction[iP][iE][jE][11] =
			dElementReconsCoeffs[11]
			- 3.0 * dA0 * dElementReconsCoeffs[15];

		// Beta, Beta, Alpha, Alpha, Alpha derivatives
		m_dReconstruction[iP][iE][jE][14] =
			dElementReconsCoeffs[14]
			- 3.0 * dB0 * dElementReconsCoeffs[15];

		// Beta, Beta, Beta, Alpha, Alpha, Alpha derivatives
		m_dReconstruction[iP][iE][jE][15] =
			dElementReconsCoeffs[15];
*/
		for (int p = 0; p < m_nNp * m_nNp; p++) {
			m_dReconstruction[iP][iAElement][iBElement][p] =
				dElementReconsCoeffs[p];
		}

		// Since sub-elements are considered to have independent
		// reconstructions we must store the corresponding reconstruction
		// for these elements as well.
		for (int i = 0; i < m_nNp-1; i++) {
		for (int j = 0; j < m_nNp-1; j++) {
			if ((i == 0) && (j == 0)) {
				continue;
			}

			int iA = iAElement + i;
			int iB = iBElement + j;

			for (int p = 0; p < m_nNp * m_nNp; p++) {
				m_dReconstruction[iP][iA][iB][p] =
					m_dReconstruction[iP][iAElement][iBElement][p];
			}
		}
		}
/*
		for (m = 0; m < m_nNp * m_nNp; m++) {
			printf("[%1.5e, %1.5e],\n", dPointValues[m], m_dReconstruction[iP][iE][jE][p]);
		}

		_EXCEPTION();
*/
	}
	}
	}
/*
	// Total mass
	double dMass = 0.0;
	for (int iP = 0; iP < 6; iP++) {
	for (int i = 0; i < nNumSubElements; i++) {
	for (int j = 0; j < nNumSubElements; j++) {

		double dArea =
			  (m_dEdgeA[i+1] - m_dEdgeA[i])
			* (m_dEdgeA[j+1] - m_dEdgeA[j]);

		dMass += m_dReconstruction[iP][i][j][0] * dArea;
	}
	}
	}
	printf("HIGHORDER MASS: %1.15e\n", dMass / (4.0 * M_PI));

	_EXCEPTION();
*/
}

///////////////////////////////////////////////////////////////////////////////

void HOMMEGrid::ReconstructLowOrderMonotone(
	const DataVector<double> & dU
) {
	// Verify data is the right size
	if (dU.GetRows() != m_dLon.GetRows()) {
		_EXCEPTION2("Invalid data size: Found %i, Expected %i",
			dU.GetRows(), m_dLon.GetRows());
	}

	// Number of sub-elements along each edge
	int nNumSubElements = m_nNe * (m_nNp - 1);

	// Initialize reconstruction coefficients
	m_dReconstruction.Initialize(
		6, nNumSubElements, nNumSubElements, m_nNp * m_nNp);

	// Temporary pointwise values for this element
	DataVector<double> dPointValues;
	dPointValues.Initialize(m_nNp * m_nNp);

	// Reconstruction coefficients for the entire element
	DataVector<double> dElementReconsCoeffs;
	dElementReconsCoeffs.Initialize(m_nNp * m_nNp);

	// Jacobian at nodes within a finite element
	DataMatrix<double> dJacobian;
	dJacobian.Initialize(m_nNp, m_nNp);

	// Reconstructed sub-element areas
	DataMatrix<double> dSubElementArea;
	dSubElementArea.Initialize(m_nNp, m_nNp);

	// Grid spacing
	double dElementDeltaA = 0.5 * M_PI / static_cast<double>(m_nNe);

	// Quadrature weights at each node
	DataVector<double> dG;
	DataVector<double> dW;

	GaussLobattoQuadrature::GetPoints(m_nNp, 0.0, 1.0, dG, dW);
/*
	DataMatrix<double> dWeight;
	dWeight.Initialize(m_nNp, m_nNp);

	for (int i = 0; i < m_nNp; i++) {
	for (int j = 0; j < m_nNp; j++) {
		dWeight[i][j] = dW[i] * dW[j] * dElementDeltaA * dElementDeltaA;
	}
	}
	for (int i = 1; i < m_nNp-1; i++) {
	for (int j = 0; j < m_nNp; j++) {
		dWeight[i][j] /= 2.0;
	}
	}
	for (int i = 0; i < m_nNp; i++) {
	for (int j = 1; j < m_nNp-1; j++) {
		dWeight[i][j] /= 2.0;
	}
	}
*/
	// Remapping matrix
	DataMatrix<double> dR;
	dR.Initialize(3, 3);
	dR[0][0] = 0.090903957942705;
	dR[0][1] = 0.251678634003417;
	dR[0][2] = 0.405738774050462;
	dR[1][1] = 0.125362728546156;
	dR[1][2] = 0.374637271453844;
	dR[2][2] = 0.25;

	// Loop through all elements
	for (int iP = 0; iP < 6; iP++) {
	for (int a = 0; a < m_nNe; a++) {
	for (int b = 0; b < m_nNe; b++) {

		int iAElement = a * (m_nNp - 1);
		int iBElement = b * (m_nNp - 1);

		// Compute metric at sub-grid points
		for (int i = 0; i < m_nNp; i++) {
		for (int j = 0; j < m_nNp; j++) {

			int iA = iAElement + i;
			int iB = iBElement + j;

			// Define the metric at this point
			double dX = m_dEdgeX[iA];
			double dY = m_dEdgeX[iB];
			double dDelta = sqrt(1.0 + dX * dX + dY * dY);

			dJacobian[i][j] =
				(1.0 + dX * dX) * (1.0 + dY * dY)
				/ (dDelta * dDelta * dDelta);
		}
		}

		if (m_nNp != 4) {
			_EXCEPTIONT("UNIMPLEMENTED");
		}

		// Copy pointwise data into point values vector
		for (int i = 0; i < m_nNp; i++) {
		for (int j = 0; j < m_nNp; j++) {

			int iA = iAElement + i;
			int iB = iBElement + j;

			// Define the metric at this point
			double dX = m_dEdgeX[iA];
			double dY = m_dEdgeX[iB];
			double dDelta = sqrt(1.0 + dX * dX + dY * dY);
			double dSqrtG = (1.0 + dX * dX) * (1.0 + dY * dY)
				/ (dDelta * dDelta * dDelta);

			// Multiply by metric weight
			dPointValues[i * m_nNp + j] = dSqrtG;
		}
		}

		// Multiply to obtain elementwise reconstruction coefficients
		dElementReconsCoeffs.Zero();

		for (int m = 0; m < m_nNp * m_nNp; m++) {
		for (int n = 0; n < m_nNp * m_nNp; n++) {
			dElementReconsCoeffs[m] +=
				m_dReconsMatrix[n][m] * dPointValues[n];
		}
		}

		// Since sub-elements are considered to have independent
		// reconstructions we must store the corresponding reconstruction
		// for these elements as well.
		for (int i = 0; i < m_nNp-1; i++) {
		for (int j = 0; j < m_nNp-1; j++) {
			int iA = iAElement + i;
			int iB = iBElement + j;

			for (int p = 0; p < m_nNp * m_nNp; p++) {
				m_dReconstruction[iP][iA][iB][p] =
					dElementReconsCoeffs[p];
			}
		}
		}
/*
		// Calculate sub-element masses using line integrals
		for (int i = 0; i < m_nNp-1; i++) {
		for (int j = 0; j < m_nNp-1; j++) {
			dSubElementArea[i][j] = 0.0;

			int iA = iAElement + i;
			int iB = iBElement + j;

			double dA0 = m_dCentroidA[iA];
			double dB0 = m_dCentroidA[iB];

			double dA1 = m_dEdgeA[iA];
			double dB1 = m_dEdgeA[iB];

			double dA2 = m_dEdgeA[iA+1];
			double dB2 = m_dEdgeA[iB+1];

			int w = 0;
			for (int p = 0; p < m_nNp; p++) {
			for (int q = 0; q < m_nNp; q++) {
				dSubElementArea[i][j] +=
					m_dReconstruction[iP][iA][iB][w]
						/ static_cast<double>(p + 1)
						/ static_cast<double>(q + 1)
						* (IPow(dA2-dA0, p + 1) * IPow(dB2-dB0, q + 1)
							- IPow(dA1-dA0, p + 1) * IPow(dB2-dB0, q + 1));

				dSubElementArea[i][j] -=
					m_dReconstruction[iP][iA][iB][w]
						/ static_cast<double>(p + 1)
						/ static_cast<double>(q + 1)
						* (IPow(dA2-dA0, p + 1) * IPow(dB1-dB0, q + 1)
							- IPow(dA1-dA0, p + 1) * IPow(dB1-dB0, q + 1));

				w++;
			}
			}
		}
		}
*/
		// Compute mass in each subcell
		for (int i = 0; i < m_nNp-1; i++) {
		for (int j = 0; j < m_nNp-1; j++) {

			int iA = iAElement + i;
			int iB = iBElement + j;

			double dValue[2][2];
			dValue[0][0] = dU[m_matNodeList[iP][iA  ][iB  ]];
			dValue[1][0] = dU[m_matNodeList[iP][iA+1][iB  ]];
			dValue[0][1] = dU[m_matNodeList[iP][iA  ][iB+1]];
			dValue[1][1] = dU[m_matNodeList[iP][iA+1][iB+1]];
/*
			// DEBUG
			dValue[0][0] = 1.0;
			dValue[0][1] = 1.0;
			dValue[1][0] = 1.0;
			dValue[1][1] = 1.0;
*/
			// Reconstruction
			double dElementMass;
			double dElementArea;
			if (i == 0) {
				if (j == 0) {
					dElementMass =
						+ dR[0][0] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[0][1] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[0][1] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[0][2] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[0][0] * dJacobian[i  ][j  ]
						+ dR[0][1] * dJacobian[i+1][j  ]
						+ dR[0][1] * dJacobian[i  ][j+1]
						+ dR[0][2] * dJacobian[i+1][j+1];

				} else if (j == 1) {
					dElementMass =
						+ dR[1][1] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[1][2] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[1][1] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[1][2] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[1][1] * dJacobian[i  ][j  ]
						+ dR[1][2] * dJacobian[i+1][j  ]
						+ dR[1][1] * dJacobian[i  ][j+1]
						+ dR[1][2] * dJacobian[i+1][j+1];

				} else if (j == 2) {
					dElementMass =
						+ dR[0][1] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[0][2] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[0][0] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[0][1] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[0][1] * dJacobian[i  ][j  ]
						+ dR[0][2] * dJacobian[i+1][j  ]
						+ dR[0][0] * dJacobian[i  ][j+1]
						+ dR[0][1] * dJacobian[i+1][j+1];
				}

			} else if (i == 1) {
				if (j == 0) {
					dElementMass =
						+ dR[1][1] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[1][1] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[1][2] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[1][2] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[1][1] * dJacobian[i  ][j  ]
						+ dR[1][1] * dJacobian[i+1][j  ]
						+ dR[1][2] * dJacobian[i  ][j+1]
						+ dR[1][2] * dJacobian[i+1][j+1];

				} else if (j == 1) {
					dElementMass =
						+ dR[2][2] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[2][2] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[2][2] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[2][2] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[2][2] * dJacobian[i  ][j  ]
						+ dR[2][2] * dJacobian[i+1][j  ]
						+ dR[2][2] * dJacobian[i  ][j+1] 
						+ dR[2][2] * dJacobian[i+1][j+1];

				} else if (j == 2) {
					dElementMass =
						+ dR[1][2] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[1][2] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[1][1] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[1][1] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[1][2] * dJacobian[i  ][j  ]
						+ dR[1][2] * dJacobian[i+1][j  ]
						+ dR[1][1] * dJacobian[i  ][j+1] 
						+ dR[1][1] * dJacobian[i+1][j+1];
				}

			} else if (i == 2) {
				if (j == 0) {
					dElementMass =
						+ dR[0][1] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[0][0] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[0][2] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[0][1] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[0][1] * dJacobian[i  ][j  ]
						+ dR[0][0] * dJacobian[i+1][j  ]
						+ dR[0][2] * dJacobian[i  ][j+1]
						+ dR[0][1] * dJacobian[i+1][j+1];

				} else if (j == 1) {
					dElementMass =
						+ dR[1][2] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[1][1] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[1][2] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[1][1] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[1][2] * dJacobian[i  ][j  ]
						+ dR[1][1] * dJacobian[i+1][j  ]
						+ dR[1][2] * dJacobian[i  ][j+1]
						+ dR[1][1] * dJacobian[i+1][j+1];

				} else if (j == 2) {
					dElementMass =
						+ dR[0][2] * dJacobian[i  ][j  ] * dValue[0][0]
						+ dR[0][1] * dJacobian[i+1][j  ] * dValue[1][0]
						+ dR[0][1] * dJacobian[i  ][j+1] * dValue[0][1]
						+ dR[0][0] * dJacobian[i+1][j+1] * dValue[1][1];

					dElementArea =
						+ dR[0][2] * dJacobian[i  ][j  ]
						+ dR[0][1] * dJacobian[i+1][j  ]
						+ dR[0][1] * dJacobian[i  ][j+1]
						+ dR[0][0] * dJacobian[i+1][j+1];
				}
			}
/*
			m_dReconstruction[iP][iA][iB][0] = dElementMass;
			for (int p = 1; p < m_nNp * m_nNp; p++) {
				m_dReconstruction[iP][iA][iB][p] = 0.0;
			}
*/
/*
			dElementArea = dSubElementArea[i][j]
				/ (m_dEdgeA[iA+1] - m_dEdgeA[iA])
				/ (m_dEdgeA[iB+1] - m_dEdgeA[iB]);
*/
			for (int p = 0; p < m_nNp * m_nNp; p++) {
				m_dReconstruction[iP][iA][iB][p] *= dElementMass / dElementArea;
			}

/*
			double dUMass =
				+ dJacobian[i  ][j  ] * dU[iNode[0][0]] * dWeight[i  ][j  ]
				+ dJacobian[i+1][j  ] * dU[iNode[1][0]] * dWeight[i+1][j  ]
				+ dJacobian[i  ][j+1] * dU[iNode[0][1]] * dWeight[i  ][j+1]
				+ dJacobian[i+1][j+1] * dU[iNode[1][1]] * dWeight[i+1][j+1];

			m_dReconstruction[iP][iA][iB][0] = dUMass;

			double dArea = 
				  (m_dEdgeA[iA+1] - m_dEdgeA[iA])
				* (m_dEdgeA[iB+1] - m_dEdgeA[iB]);

			m_dReconstruction[iP][iA][iB][0] /= dArea;
*/
/*
			// DEBUG
			double dA0 = m_dCentroidA[iA];
			double dB0 = m_dCentroidA[iB];

			double dA1 = m_dEdgeA[iA];
			double dB1 = m_dEdgeA[iB];

			double dA2 = m_dEdgeA[iA+1];
			double dB2 = m_dEdgeA[iB+1];

			int w = 0;
			double dNewElementMass = 0.0;
			for (int p = 0; p < m_nNp; p++) {
			for (int q = 0; q < m_nNp; q++) {
				dNewElementMass +=
					m_dReconstruction[iP][iA][iB][w]
						/ static_cast<double>(p + 1)
						/ static_cast<double>(q + 1)
						* (IPow(dA2-dA0, p + 1) * IPow(dB2-dB0, q + 1)
							- IPow(dA1-dA0, p + 1) * IPow(dB2-dB0, q + 1));

				dNewElementMass -=
					m_dReconstruction[iP][iA][iB][w]
						/ static_cast<double>(p + 1)
						/ static_cast<double>(q + 1)
						* (IPow(dA2-dA0, p + 1) * IPow(dB1-dB0, q + 1)
							- IPow(dA1-dA0, p + 1) * IPow(dB1-dB0, q + 1));

				w++;
			}
			}

			dNewElementMass /= (dA2-dA1) * (dB2-dB1);

			printf("%i %i : %1.14e %1.14e\n",
				i, j, dElementMass, dNewElementMass);
*/
		}
		}
		//_EXCEPTION();

	}
	}
	}
/*
	// Sum mass in cell
	for (int iP = 0; iP < 6; iP++) {
	for (int a = 0; a < m_nNe; a++) {
	for (int b = 0; b < m_nNe; b++) {

		int iAElement = a * (m_nNp - 1);
		int iBElement = b * (m_nNp - 1);

		double dMass = 0.0;
		for (int i = 0; i < m_nNp-1; i++) {
		for (int j = 0; j < m_nNp-1; j++) {

			int iA = iAElement + i;
			int iB = iBElement + j;

			dMass += m_dReconstruction[iP][iA][iB][0];
		}
		}
		for (int i = 0; i < m_nNp-1; i++) {
		for (int j = 0; j < m_nNp-1; j++) {

			int iA = iAElement + i;
			int iB = iBElement + j;

			m_dReconstruction[iP][iA][iB][0] = dMass;
		}
		}
	}
	}
	}

	std::cout << "LOWORDER: " << m_dReconstruction[0][0][0][0] << std::endl;
*/
/*
	// Total mass
	double dMass = 0.0;
	for (int iP = 0; iP < 6; iP++) {
	for (int i = 0; i < nNumSubElements; i++) {
	for (int j = 0; j < nNumSubElements; j++) {

		double dArea =
			  (m_dEdgeA[i+1] - m_dEdgeA[i])
			* (m_dEdgeA[j+1] - m_dEdgeA[j]);

		dMass += m_dReconstruction[iP][i][j][0] * dArea;
	}
	}
	}
	printf("\nLOWORDER MASS: %1.15e\n", dMass / (4.0 * M_PI));

	_EXCEPTION();
*/
}

///////////////////////////////////////////////////////////////////////////////

