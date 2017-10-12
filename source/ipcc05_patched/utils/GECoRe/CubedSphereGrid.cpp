///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereGrid.cpp
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

#include "CubedSphereGrid.h"
#include "CubedSphereTrans.h"

#include <cmath>
#include <cfloat>

///////////////////////////////////////////////////////////////////////////////

void CubedSphereGrid::IndexFromXY(
	double dX,
	double dY,
	int & iX,
	int & iY
) const {
	int i;

	iX = static_cast<int>(-1);
	iY = static_cast<int>(-1);

	if (dX < GetEdgeX()[GetInteriorBegin()]) {
		iX = GetInteriorBegin();
	}
	if (dY < GetEdgeX()[GetInteriorBegin()]) {
		iY = GetInteriorBegin();
	}

	for (i = GetInteriorBegin(); i < GetInteriorEnd(); i++) {
		if ((dX >= GetEdgeX()[i]) && (dX <= GetEdgeX()[i+1])) {
			iX = i;
		}
		if ((dY >= GetEdgeX()[i]) && (dY <= GetEdgeX()[i+1])) {
			iY = i;
		}
	}

	if (iX == static_cast<int>(-1)) {
		iX = GetInteriorEnd()-1;
	}
	if (iY == static_cast<int>(-1)) {
		iY = GetInteriorEnd()-1;
	}
}

///////////////////////////////////////////////////////////////////////////////

