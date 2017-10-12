///////////////////////////////////////////////////////////////////////////////
///
///	\file    Interpolator.cpp
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

#include "Interpolator.h"

#include "MathHelper.h"

////////////////////////////////////////////////////////////////////////////////

const double Interpolator::TINY = 1e-12;

////////////////////////////////////////////////////////////////////////////////

void Interpolator::AddLineSegment(
	const CubedSphereGrid & gridCS,
	const LatitudeLongitudeGrid & gridRLL,
	LineType nLineType,
	double dX1,
	double dY1,
	double dX2,
	double dY2,
	int nPanel,
	int ixLon,
	int ixLat,
	int ixX,
	int ixY
) {
	int i;
	int j;

	double dWeight[16];
	double dWeightAux[16];

	int ixPrev;

	// Element centroids
	const DataVector<double> & dCentroidA = gridCS.GetCentroidA();

	// Check for coincident begin and endpoints
	if ((fabs(dX1 - dX2) < TINY) && (fabs(dY1 - dY2) < TINY)) {
		return;
	}

	// Calculate the weights associated with this line segment
	CalculateLineSegmentWeight(
		nLineType, dCentroidA[ixX], dCentroidA[ixY],
		dX1, dY1, dX2, dY2, nPanel, dWeight);

	// Check for zero weights
	for (i = 0; i < m_nWeights; i++) {
		if (fabs(dWeight[i]) > TINY) {
			break;
		}
	}
	if (i == m_nWeights) {
		return;
	}

	// Lines of constant Y
	if (nLineType == LineType_Beta) {

		// Line passes through a RLL element
		if (ixY != gridCS.GetInteriorBegin()) {
			CalculateLineSegmentWeight(
				nLineType, dCentroidA[ixX], dCentroidA[ixY-1],
				dX1, dY1, dX2, dY2, nPanel, dWeightAux);

			weights[ixLat][ixLon].AddWeight(
				ixX, ixY-1, nPanel, dWeightAux, m_nWeights, 1.0);
		}
		if (ixY != gridCS.GetInteriorEnd()) {
			weights[ixLat][ixLon].AddWeight(
				ixX, ixY, nPanel, dWeight, m_nWeights, -1.0);
		}

	// Lines of constant latitude
	} else if (nLineType == LineType_Lat) {

		// Determine the previous latitude element
		if ((ixLat == 0) || (ixLat == gridRLL.GetLatitudes())) {
			_EXCEPTION();
		}
		ixPrev = ixLat - 1;

		// Adjust the weights of the adjacent RLL elements, giving special
		// consideration to the equator which may be coincident with a line
		// of constant Y.
		if ((nPanel < 4) && (ixY == m_jMiddle) &&
			(fabs(dY1) < TINY) && (fabs(dY2) < TINY)
		) {
			CalculateLineSegmentWeight(
				nLineType, dCentroidA[ixX], dCentroidA[ixY-1],
				dX1, dY1, dX2, dY2, nPanel, dWeightAux);

			weights[ixLat][ixLon].AddWeight(
				ixX, ixY, nPanel, dWeight, m_nWeights, -1.0);
			weights[ixPrev][ixLon].AddWeight(
				ixX, ixY-1, nPanel, dWeightAux, m_nWeights, 1.0);
		} else {
			weights[ixLat][ixLon].AddWeight(
				ixX, ixY, nPanel, dWeight, m_nWeights, -1.0);
			weights[ixPrev][ixLon].AddWeight(
				ixX, ixY, nPanel, dWeight, m_nWeights, 1.0);
		}

	// Lines of constant longitude
	} else if (nLineType == LineType_Lon) {

		// Determine the previous longitude element
		if (ixLon == 0) {
			ixPrev = gridRLL.GetLongitudes()-1;
		} else {
			ixPrev = ixLon - 1;
		}

		// If we are on a longitude line coincident with a line of constant Y
		// then we must apply special consideration.
		if ((nPanel > 3) && (m_jMiddle != static_cast<int>(-1)) &&
			(fabs(dY1) < TINY) && (fabs(dY2) < TINY)
		) {

			// Longitude = pi/2
			if (fabs(gridRLL.GetLonEdge()[ixLon] - 0.5 * M_PI) < TINY) {
				CalculateLineSegmentWeight(
					nLineType, dCentroidA[ixX], dCentroidA[ixY-1],
					dX1, dY1, dX2, dY2, nPanel, dWeightAux);

				if (nPanel == 4) {
					weights[ixLat][ixLon].AddWeight(
						ixX, ixY, nPanel, dWeight, m_nWeights, 1.0);
					weights[ixLat][ixPrev].AddWeight(
						ixX, ixY-1, nPanel, dWeightAux, m_nWeights, -1.0);

				} else {
					weights[ixLat][ixLon].AddWeight(
						ixX, ixY-1, nPanel, dWeightAux, m_nWeights, 1.0);
					weights[ixLat][ixPrev].AddWeight(
						ixX, ixY, nPanel, dWeight, m_nWeights, -1.0);
				}

			// Longitude = 3 pi / 2
			} else {
				CalculateLineSegmentWeight(
					nLineType, dCentroidA[ixX], dCentroidA[ixY+1],
					dX1, dY1, dX2, dY2, nPanel, dWeightAux);

				if (nPanel == 4) {
					weights[ixLat][ixLon].AddWeight(
						ixX, ixY, nPanel, dWeight, m_nWeights, 1.0);
					weights[ixLat][ixPrev].AddWeight(
						ixX, ixY+1, nPanel, dWeightAux, m_nWeights, -1.0);

				} else {
					weights[ixLat][ixLon].AddWeight(
						ixX, ixY+1, nPanel, dWeightAux, m_nWeights, 1.0);
					weights[ixLat][ixPrev].AddWeight(
						ixX, ixY, nPanel, dWeight, m_nWeights, -1.0);
				}
			}

		// Adjust the weights of the adjacent RLL elements
		} else {
			weights[ixLat][ixLon].AddWeight(
				ixX, ixY, nPanel, dWeight, m_nWeights, 1.0);
			weights[ixLat][ixPrev].AddWeight(
				ixX, ixY, nPanel, dWeight, m_nWeights, -1.0);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void Interpolator::InitializeWeights(
	const CubedSphereGrid &gridCS,
	const LatitudeLongitudeGrid &gridRLL
) {
	int i;
	int j;
	int k;
	int n;

	// Initialize the weight storage array
	weights.Initialize(gridRLL.GetLatitudes(), gridRLL.GetLongitudes());

	// Determine the "middle" index for lines of constant beta
	m_jMiddle = static_cast<int>(-1);
	for (j = gridCS.GetInteriorBegin(); j < gridCS.GetInteriorEnd(); j++) {
		if (fabs(gridCS.GetEdgeX()[j]) < TINY) {
			m_jMiddle = j;
		}
	}

	// There exists a "middle" index.  Check for coincident longitude lines
	// (on polar panels) and coincident latitude lines (on equatorial panels).
	bool fHasMiddleLon1 = false;
	bool fHasMiddleLon2 = false;
	bool fHasMiddleLat = false;

	if (m_jMiddle != static_cast<int>(-1)) {
		for (j = 0; j < gridRLL.GetLongitudes(); j++) {
			if (fabs(gridRLL.GetLonEdge()[j] - 0.5 * M_PI) < TINY) {
				fHasMiddleLon1 = true;
			}
			if (fabs(gridRLL.GetLonEdge()[j] - 1.5 * M_PI) < TINY) {
				fHasMiddleLon2 = true;
			}
		}

		for (j = 0; j < gridRLL.GetLatitudes(); j++) {
			if (fabs(gridRLL.GetLatEdge()[j]) < TINY) {
				fHasMiddleLat = true;
				break;
			}
		}
	}

	// No implementation for split middle longitudes yet
	if (((fHasMiddleLon1) && (!fHasMiddleLon2)) ||
		((!fHasMiddleLon1) && (fHasMiddleLon2))
	) {
		_EXCEPTIONT("Grid must either have longitude lines at both "
			"pi/2 and 3pi/2 or at neither (implementation).");
	}

	// Loop over all longitude segments on the grid, careful to ignore
	// longitude line segments at pi/2 and 3pi/2, which correspond to lines
	// of constant Y.
	for (i = 0; i < gridRLL.GetLongitudes(); i++) {
		InitializeLongitudeSegment(gridCS, gridRLL, i);
	}

	// Loop over all latitude segments on the grid, careful to ignore the
	// equatorial line segment, which corresponds to a line of constant Y.
	for (i = 0; i < gridRLL.GetLatitudes(); i++) {
		InitializeLatitudeSegment(gridCS, gridRLL, i);
	}

	// Loop over all beta segments (lines of constant Y) on the grid
	for (i = 0; i < 6; i++) {
		for (j = gridCS.GetInteriorBegin(); j <= gridCS.GetInteriorEnd(); j++) {
			if (j == m_jMiddle) {
				if ((i < 4) && (fHasMiddleLat)) {
					continue;
				}
				if ((i > 3) && (fHasMiddleLon1) && (fHasMiddleLon2)) {
					continue;
				}
			}
			InitializeBetaSegment(gridCS, gridRLL, j, i);
		}
	}

	// (DEBUG)
	//InitializeBetaSegment(gridCS, gridRLL, 6, 4);
	//InitializeLatitudeSegment(gridCS, gridRLL, 5);
	//InitializeLongitudeSegment(gridCS, gridRLL, 4);

/*
	// Output the results in a given grid cell
	double dSum = 0.0;
	//for (i = 0; i < gridRLL.GetLongitudes(); i++) {
	//for (j = 0; j < gridRLL.GetLatitudes(); j++) {
		ElementWeightVector &vec = weights[0][0];
		for (k = 0; k < vec.size(); k++) {
			std::cout << vec[k].ixX << ", " << vec[k].ixY << ", " << vec[k].ixP << " : " << vec[k].dWeight[4] << std::endl;
			//dSum += vec[k].dWeight[1];
		}
	//}
	//}
*/

/*
	double dSum = 0.0;
	for (i = 0; i < gridRLL.GetLatitudes(); i++) {
	for (j = 0; j < gridRLL.GetLongitudes(); j++) {
		ElementWeightVector &vec = weights[i][j];
		for (k = 0; k < vec.size(); k++) {
			//std::cout << vec[k].ixX << ", " << vec[k].ixY << ", " << vec[k].ixP << " : " << vec[k].dWeight[0] << std::endl;
			dSum += vec[k].dWeight[4];
		}
	}
	}
	std::cout << "Area: " << dSum << std::endl;
*/

/*
	// (DEBUG) Check weights
	for (i = 0; i < gridRLL.GetLatitudes(); i++) {
	for (j = 0; j < gridRLL.GetLongitudes(); j++) {
		for (k = 0; k < weights[i][j].size(); k++) {
			// Check for negative overlap areas
			if (weights[i][j][k].dWeight[0] < - TINY) {
				std::cout << "In: " << i << ", " << j << std::endl;
				_EXCEPTIONT("Negative overlap area detected.");
			}

			if (m_nOrder == 1) {
				continue;
			}

			if (weights[i][j][k].dWeight[1] *
				gridCS.GetCentroidX()[weights[i][j][k].ixX] < - TINY
			) {
				std::cout << "In: " << i << ", " << j << std::endl;
				std::cout << "CS: " << weights[i][j][k].ixX
					<< ", " << weights[i][j][k].ixY
					<< ", " << weights[i][j][k].ixP << std::endl;
				std::cout << "Weight: " << weights[i][j][k].dWeight[1] << std::endl;

				_EXCEPTION();
			}

			if (weights[i][j][k].dWeight[2] *
				gridCS.GetCentroidX()[weights[i][j][k].ixY] < - TINY
			) {
				std::cout << "In: " << i << ", " << j << std::endl;
				std::cout << "CS: " << weights[i][j][k].ixX
					<< ", " << weights[i][j][k].ixY
					<< ", " << weights[i][j][k].ixP << std::endl;
				std::cout << "Weight: " << weights[i][j][k].dWeight[2] << std::endl;

				_EXCEPTION();
			}

			if (m_nOrder == 2) {
				continue;
			}

			if (weights[i][j][k].dWeight[3] < - TINY) {
				std::cout << "In: " << i << ", " << j << std::endl;
				std::cout << "CS: " << weights[i][j][k].ixX
					<< ", " << weights[i][j][k].ixY
					<< ", " << weights[i][j][k].ixP << std::endl;

				_EXCEPTIONT("Negative weight[3] detected.");
			}

			if (weights[i][j][k].dWeight[5] < - TINY) {
				std::cout << "In: " << i << ", " << j << std::endl;
				std::cout << "CS: " << weights[i][j][k].ixX
					<< ", " << weights[i][j][k].ixY
					<< ", " << weights[i][j][k].ixP << std::endl;

				_EXCEPTIONT("Negative weight[5] detected.");
			}

			if (weights[i][j][k].dWeight[4] < - TINY) {
				std::cout << "In: " << i << ", " << j << std::endl;
				std::cout << "CS: " << weights[i][j][k].ixX
					<< ", " << weights[i][j][k].ixY
					<< ", " << weights[i][j][k].ixP << std::endl;
				std::cout << "Wt: " << weights[i][j][k].dWeight[4] << std::endl;

				_EXCEPTIONT("Negative weight[4] detected.");
			}

		}
	}
	}
*/

}

////////////////////////////////////////////////////////////////////////////////

void Interpolator::Initialize(
	const CubedSphereGrid &gridCS,
	const LatitudeLongitudeGrid &gridRLL,
	int nOrder
) {

/*
	if (nOrder == 1) {
		m_nWeights = 1;
	} else if (nOrder == 2) {
		m_nWeights = 3;
	} else if (nOrder == 3) {
		m_nWeights = 6;
	}
*/
	// Initialize the weights for this grid
	InitializeWeights(gridCS, gridRLL);
}

////////////////////////////////////////////////////////////////////////////////

void Interpolator::RemapCStoRLL(
	const CubedSphereGrid & gridCS,
	const LatitudeLongitudeGrid & gridRLL,
	int nInVar,
	DataMatrix<double> & dataRLL,
	bool fUseRemappedArea
) {
	int i;
	int j;

	int k;

	int iX;
	int iY;
	int iP;

	int iW;

	// Element area
	const DataMatrix<double> * pArea;
	if (fUseRemappedArea) {
		pArea = &(gridRLL.GetRemappedElementArea());
	} else {
		pArea = &(gridRLL.GetElementArea());
	}

	// Loop through all elements on the RLL grid
	for (i = 0; i < gridRLL.GetLatitudes(); i++) {
	for (j = 0; j < gridRLL.GetLongitudes(); j++) {

		// Obtain the element weight vector
		ElementWeightVector &vec = weights[i][j];

		// Remap
		dataRLL[i][j] = 0.0;
		for (k = 0; k < vec.size(); k++) {

			// Obtain panel and CS element index
			iX = vec[k].ixX;
			iY = vec[k].ixY;
			iP = vec[k].ixP;

			// Get the reconstruction coefficients
			const double * dRecons = gridCS.GetReconstruction(iP, iX, iY);

			for (iW = 0; iW < m_nWeights; iW++) {
				dataRLL[i][j] +=
					dRecons[iW] * vec[k].dWeight[iW];
			}
		}

		// Take the element average
		dataRLL[i][j] /= (*pArea)[i][j];
	}
	}
}

////////////////////////////////////////////////////////////////////////////////

void Interpolator::CropExtremeValues(
	const LatitudeLongitudeGrid &gridRLL,
	DataMatrix<double> & dataRLL,
	double dMin,
	double dMax
) {
	
	// Loop through all elements on the RLL grid
	for (int i = 0; i < gridRLL.GetLatitudes(); i++) {
	for (int j = 0; j < gridRLL.GetLongitudes(); j++) {
		if (dataRLL[i][j] < dMin) {
			dataRLL[i][j] = dMin;
		}
		if (dataRLL[i][j] > dMax) {
			dataRLL[i][j] = dMax;
		}
	}
	}
}

////////////////////////////////////////////////////////////////////////////////

