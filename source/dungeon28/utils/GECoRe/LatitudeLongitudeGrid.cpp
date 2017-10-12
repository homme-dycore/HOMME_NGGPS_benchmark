///////////////////////////////////////////////////////////////////////////////
///
///	\file    LatitudeLongitudeGrid.cpp
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

#include "LatitudeLongitudeGrid.h"

#include <cmath>
#include <cfloat>

///////////////////////////////////////////////////////////////////////////////

void LatitudeLongitudeGrid::Initialize(
	int nLatitudes,
	int nLongitudes,
	double dLongitudeShift
) {
	int m;
	int n;

	double dN;

	double dFirstA;
	
	double dDeltaA;

	// Number of ghost elements around static grid
	m_nLatitudes = nLatitudes;

	// Number of elements in each direction
	m_nLongitudes = nLongitudes;

	// Shift in longitude
	m_dLongitudeShift = fmod(dLongitudeShift, 2.0 * M_PI);

	// Initialize edge coordinate vectors
	m_matLatEdge.Initialize(m_nLatitudes + 1);

	m_matLonEdge.Initialize(m_nLongitudes + 1);

	// Calculate coordinates of all latitude edges
	dFirstA = - 0.5 * M_PI;
	dDeltaA = M_PI / static_cast<double>(m_nLatitudes);
	for (n = 0; n < m_nLatitudes + 1; n++) {
		dN = static_cast<double>(n);

		m_matLatEdge[n] = dFirstA + dDeltaA * dN;
	}

	// Calculate coordinates of all longitude edges
	dFirstA = 0.0;
	dDeltaA = 2.0 * M_PI / static_cast<double>(m_nLongitudes);
	for (n = 0; n < m_nLongitudes + 1; n++) {
		dN = static_cast<double>(n);

		m_matLonEdge[n] = dFirstA + dDeltaA * dN;
	}

	// Initialize centroid coordinate vectors
	m_matLatCentroid.Initialize(m_nLatitudes);

	m_matLonCentroid.Initialize(m_nLongitudes);

	// Calculate coordinates of all latitude centroids
	for (n = 0; n < m_nLatitudes; n++) {
		m_matLatCentroid[n] = 0.5 * (m_matLatEdge[n] + m_matLatEdge[n+1]);
	}

	// Calculate coordinates of all longitude centroids
	for (n = 0; n < m_nLongitudes; n++) {
		m_matLonCentroid[n] = 0.5 * (m_matLonEdge[n] + m_matLonEdge[n+1]);
	}

	// Initialize element area matrix
	m_matElementArea.Initialize(m_nLatitudes, m_nLongitudes);

	// Calculate element area of all elements
	for (m = 0; m < m_nLatitudes; m++) {
	for (n = 0; n < m_nLongitudes; n++) {
		m_matElementArea[m][n] =
			(m_matLonEdge[n+1] - m_matLonEdge[n]) *
			(sin(m_matLatEdge[m+1]) - sin(m_matLatEdge[m]));
	}
	}

	// Scale longitude edges and centroids to be in the range [0, 2*pi)
	for (n = 0; n < m_nLongitudes + 1; n++) {
		m_matLonEdge[n] += dLongitudeShift;
		m_matLonEdge[n] = fmod(m_matLonEdge[n], 2.0 * M_PI);
	}
	for (n = 0; n < m_nLongitudes; n++) {
		m_matLonCentroid[n] += dLongitudeShift;
		m_matLonCentroid[n] = fmod(m_matLonCentroid[n], 2.0 * M_PI);
	}

}

///////////////////////////////////////////////////////////////////////////////

void LatitudeLongitudeGrid::SetRemappedGridArea(
	const DataMatrix<double> & matRemappedElementArea
) {
	m_matRemappedElementArea = matRemappedElementArea;
}

///////////////////////////////////////////////////////////////////////////////

void LatitudeLongitudeGrid::IndexFromRLL(
	double dLon,
	double dLat,
	int &iLon,
	int &iLat
) const {
	const double TINY = 1.0e-12;

	int i;

	iLon = (-1);
	iLat = (-1);

	for (i = 0; i < GetLongitudes(); i++) {
		double dLon1 = GetLonEdge()[i];
		double dLon2 = GetLonEdge()[i+1];

		if (dLon2 < dLon1) {
			if ((dLon > dLon1) || (dLon < dLon2 + TINY)) {
				iLon = i;
				break;
			}

		} else {
			if ((dLon > dLon1) && (dLon < dLon2 + TINY)) {
				iLon = i;
				break;
			}
		}
	}
	for (i = 0; i < GetLatitudes(); i++) {
		if (dLat > GetLatEdge()[i]) {
			iLat = i;
		} else {
			break;
		}
	}

	if ((iLon == (-1)) || (iLat == (-1))) {
		std::cout << dLon << ", " << dLat << std::endl;
		for (i = 0; i < GetLongitudes()+1; i++) {
			printf("%1.5e\n", GetLonEdge()[i]);
		}
		_EXCEPTIONT("Error in IndexFromRLL");
	}
}

///////////////////////////////////////////////////////////////////////////////

