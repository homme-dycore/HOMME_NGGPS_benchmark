///////////////////////////////////////////////////////////////////////////////
///
///	\file    HOMMEInterpolator.cpp
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

#include "HOMMEInterpolator.h"

#include "MathHelper.h"

////////////////////////////////////////////////////////////////////////////////

void HOMMEInterpolator::CalculateLineSegmentWeight(
	LineType nLineType,
	double dA0,
	double dB0,
	double dX1,
	double dY1,
	double dX2,
	double dY2,
	int nPanel,
	double * dWeight
) {

	// Gaussian quadrature points
	static const int QuadraturePoints = 5;

	double dG[QuadraturePoints];
	dG[0] = - 0.90617984593866399282;
	dG[1] = - 0.53846931010568309105;
	dG[2] =   0.0;
	dG[3] = + 0.53846931010568309105;
	dG[4] = + 0.90617984593866399282;

	double dW[QuadraturePoints];
	dW[0] = 0.23692688505618908752;
	dW[1] = 0.47862867049936646804;
	dW[2] = 0.56888888888888888889;
	dW[3] = 0.47862867049936646804;
	dW[4] = 0.23692688505618908752;

	// Calculate deltas
	double dDelta1 = sqrt(1.0 + dX1 * dX1 + dY1 * dY1);
	double dDelta2 = sqrt(1.0 + dX2 * dX2 + dY2 * dY2);

	// Obtain latitude and longitude coordinates
	double dLat1;
	double dLon1;

	double dLat2;
	double dLon2;

	CubedSphereTrans::RLLFromXYP(dX1, dY1, nPanel, dLon1, dLat1);
	CubedSphereTrans::RLLFromXYP(dX2, dY2, nPanel, dLon2, dLat2);

	if (dLon1 > 2.0 * M_PI - TINY) {
		dLon1 = 0.0;
	}
	if (dLon2 > 2.0 * M_PI - TINY) {
		dLon2 = 0.0;
	}

	// Panel indicator (for polar panels)
	double dP = (nPanel == 4)?(1.0):(-1.0);

	// Buffer variables
	double dInt1;
	double dInt2;

	// Weight id
	//int w;

	//int p;
	//int q;

	// Buffer data for computing Gaussian quadrature
	//int m;

	double dA1 = atan(dX1);
	double dA2 = atan(dX2);

	double dB1 = atan(dY1);
	double dB2 = atan(dY2);

	double dA;
	double dB;

	// Loop through all weights
	for (int w = 0; w < m_nOrder * m_nOrder; w++) {
		dWeight[w] = 0.0;
	}

	// Lines of constant beta
	if (nLineType == LineType_Beta) {

		// Loop through all weights
		int w = 0;
		for (int p = 0; p < m_nOrder; p++) {
		for (int q = 0; q < m_nOrder; q++) {
			dWeight[w] =
				1.0 / static_cast<double>(p + 1)
					/ static_cast<double>(q + 1)
				* (IPow(dA2-dA0, p + 1) * IPow(dB2-dB0, q + 1)
					- IPow(dA1-dA0, p + 1) * IPow(dB1-dB0, q + 1));

			w++;
		}
		}

	// Lines of constant longitude on polar panels
	} else if ((nPanel > 3) && (nLineType == LineType_Lon)) {

		// Fix longitude line to be constant even at the pole
		if ((fabs(dX2) < TINY) && (fabs(dY2) < TINY)) {
			dLon2 = dLon1;
		}

		if ((fabs(dX1) < TINY) && (fabs(dY1) < TINY)) {
			dLon1 = dLon2;
		}

		// Verify this is a line of constant longitude
		if (fabs(dLon1 - dLon2) > TINY) {
			_EXCEPTION2(
				"Line of constant longitude does not have "
				"constant longitude: %1.5e %1.5e", dLon1, dLon2);
		}

		// Ignore lines of constant alpha
		if ((fabs(dX1) < TINY) && (fabs(dX2) < TINY)) {
			return;
		}

		// Loop through all weights
		for (int w = 0; w < m_nOrder * m_nOrder; w++) {
			dWeight[w] = 0.0;
		}

		for (int m = 0; m < QuadraturePoints; m++) {

			dA = 0.5 * (dA1 + dA2) + dG[m] * 0.5 * (dA2 - dA1);
			dB = -dP * atan(tan(dA) / tan(dLon1));

			int w = 0;
			for (int p = 0; p < 4; p++) {
			for (int q = 0; q < 4; q++) {
				dWeight[w] +=
					dW[m] / static_cast<double>(q + 1)
					* IPow(dA-dA0, p)
					* IPow(dB-dB0, q + 1);

				w++;
			}
			}
		}

		for (int w = 0; w < m_nWeights; w++) {
			dWeight[w] *= 0.5 * (dA2 - dA1);
		}

	// Lines of constant latitude on equatorial panels
	} else if ((nPanel < 4) && (nLineType == LineType_Lat)) {

		// Verify this is, in fact, a line of constant latitude
		if (fabs(dLat1 - dLat2) > TINY) {
			_EXCEPTIONT(
				"Line of constant latitude does not have "
				"constant latitude.");
		}

		// Loop through all weights
		for (int w = 0; w < m_nWeights; w++) {
			dWeight[w] = 0.0;
		}

		for (int m = 0; m < QuadraturePoints; m++) {

			dA = 0.5 * (dA1 + dA2) + dG[m] * 0.5 * (dA2 - dA1);
			dB = atan(tan(dLat2) / cos(dA));

			int w = 0;
			for (int p = 0; p < m_nOrder; p++) {
			for (int q = 0; q < m_nOrder; q++) {
				dWeight[w] +=
					dW[m] / static_cast<double>(q + 1)
					* IPow(dA-dA0, p)
					* IPow(dB-dB0, q + 1);

				w++;
			}
			}
		}

		for (int w = 0; w < m_nWeights; w++) {
			dWeight[w] *= 0.5 * (dA2 - dA1);
		}

	// Lines of constant latitude on polar panels
	} else if ((nPanel > 3) && (nLineType == LineType_Lat)) {

		// Verify this is, in fact, a line of constant latitude
		if (fabs(dLat1 - dLat2) > TINY) {
			_EXCEPTIONT(
				"Line of constant latitude does not have "
				"constant latitude.");
		}

		// Need to special case the situation where beta crosses 0
		if (dB1 * dB2 < -TINY) {
			_EXCEPTIONT("\nLine of constant longitude that crosses beta = 0 "
				"detected.\nThis is not yet implemented.");
		}

		// Loop through all weights
		for (int w = 0; w < m_nWeights; w++) {
			dWeight[w] = 0.0;
		}

		for (int m = 0; m < QuadraturePoints; m++) {

			dA = 0.5 * (dA1 + dA2) + dG[m] * 0.5 * (dA2 - dA1);

			dB = atan(sqrt(
				1.0 / (tan(dLat1) * tan(dLat1))
				- tan(dA) * tan(dA)));

			if (dB * (dB1 + dB2) < -TINY) {
				dB *= -1.0;
			}

			int w = 0;
			for (int p = 0; p < m_nOrder; p++) {
			for (int q = 0; q < m_nOrder; q++) {
				dWeight[w] +=
					dW[m] / static_cast<double>(q + 1)
					* IPow(dA-dA0, p)
					* IPow(dB-dB0, q + 1);

				w++;
			}
			}
		}

		for (int w = 0; w < m_nWeights; w++) {
			dWeight[w] *= 0.5 * (dA2 - dA1);
		}

	// Lines of constant alpha
	} else {

		// Loop through all weights
		for (int w = 0; w < m_nWeights; w++) {
			dWeight[w] = 0.0;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void HOMMEInterpolator::Initialize(
	const CubedSphereGrid &gridCS,
	const LatitudeLongitudeGrid &gridRLL,
	int nOrder
) {
	// Initialize the order of this method and the associated number of weights
	m_nOrder = nOrder;
	m_nWeights = nOrder * nOrder;

	// Initialize the weights for this grid
	Interpolator::Initialize(gridCS, gridRLL, nOrder);
}

////////////////////////////////////////////////////////////////////////////////

