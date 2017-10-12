///////////////////////////////////////////////////////////////////////////////
///
///	\file    HOMMEInterpolator.h
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

#ifndef _HOMMEINTERPOLATOR_H_
#define _HOMMEINTERPOLATOR_H_

#include "Interpolator.h"
#include "CubedSphereGrid.h"
#include "LatitudeLongitudeGrid.h"
#include "SystemState.h"

#include "Exception.h"
#include "DataMatrix3D.h"
#include "MathHelper.h"

#include <vector>
#include <cfloat>

////////////////////////////////////////////////////////////////////////////////

class HOMMEInterpolator : public Interpolator {

	public:
		///	<summary>
		///		Default constructor.
		///	</summary>
		HOMMEInterpolator() :
			Interpolator()
		{ }

		///	<summary>
		///		Constructor.
		///	</summary>
		HOMMEInterpolator(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL,
			int nOrder
		) {
			// Initialize the interpolator
			Initialize(gridCS, gridRLL, nOrder);
		}

	protected: 
		///	<summary>
		///		Calculate the weights associated with a given line segment.
		///	</summary>
		virtual void CalculateLineSegmentWeight(
			LineType nLineType,
			double dA0,
			double dB0,
			double dX1,
			double dY1,
			double dX2,
			double dY2,
			int nPanel,
			double * dWeight
		);

	public:
		///	<summary>
		///		Initialize the interpolator for the given grid.
		///	</summary>
		virtual void Initialize(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL,
			int nOrder
		);
};

////////////////////////////////////////////////////////////////////////////////

#endif

