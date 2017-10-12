///////////////////////////////////////////////////////////////////////////////
///
///	\file    Interpolator_search.cpp
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

//#define VERBOSE1

////////////////////////////////////////////////////////////////////////////////

void Interpolator::InitializeLongitudeSegment(
	const CubedSphereGrid &gridCS,
	const LatitudeLongitudeGrid &gridRLL,
	int nLon
) {
	// Current longitude
	const double dLon = gridRLL.GetLonEdge()[nLon];
	if ((dLon < 0.0) || (dLon > 2.0 * M_PI + TINY)) {
		_EXCEPTIONT("Longitude must be defined on the range [0, 2*pi].");
	}

	// Current latitude
	double dLat;
	double dNextLat1;
	double dNextLat2;

	// Current latitude element
	int iLat = 0;

	// Current XYP coordinates
	double dX1 = 0.0;
	double dY1 = 0.0;

	double dX2 = 0.0;
	double dY2 = 0.0;

	int nP1;
	int nP2;

	// Current cubed sphere element
	int iX;
	int iY;

	// Next intersection with a cubed sphere edge
	int iNextX;
	int iNextY;

	// Flag indicating if a panel switch has occurred
	bool fPanelSwitch = false;

	// Buffer variables
	double dTemp;
	int nTemp;

	// Quadrant of the polar panels that this longitude line passes through
	int iQuad;
	if (dLon < 0.5 * M_PI) {
		iQuad = 1;
	} else if (dLon < M_PI) {
		iQuad = 2;
	} else if (dLon < 1.5 * M_PI) {
		iQuad = 3;
	} else {
		iQuad = 4;
	}

	// Begin searching
	dLat = - 0.5 *  M_PI;

	dX1 = 0.0;
	dY1 = 0.0;
	nP1 = 5;

	// Initial element
	CubedSphereTrans::XYPFromRLL(dLon, dLat + TINY, dX2, dY2, nP2);

	// Note that we must apply special consideration to longitude lines
	// with lon = pi/2 or 3pi/2 since these might be coincident with lines
	// of constant Y.  Hence, we enforce that we are always in the XY
	// element above pi/2 and below 3pi/2.
	if (fabs(dLon - 0.5 * M_PI) < TINY) {
		gridCS.IndexFromXY(dX2, TINY, iX, iY);
	} else if (fabs(dLon - 1.5 * M_PI) < TINY) {
		gridCS.IndexFromXY(dX2, -TINY, iX, iY);
	} else {
		gridCS.IndexFromXY(dX2, dY2, iX, iY);
	}

	// Begin searching
	for (;;) {

		// Panel indicator (for polar panels)
		double dP = (nP1 == 4)?(1.0):(-1.0);

		// Determine direction of movement through cubed sphere elements
		if (nP1 == 4) {
			if (iQuad == 1) {
				iNextX = iX;
				iNextY = iY + 1;

			} else if (iQuad == 2) {
				iNextX = iX;
				iNextY = iY;

			} else if (iQuad == 3) {
				iNextX = iX + 1;
				iNextY = iY;

			} else {
				iNextX = iX + 1;
				iNextY = iY + 1;
			}

		} else {
			if (iQuad == 1) {
				iNextX = iX + 1;
				iNextY = iY + 1;

			} else if (iQuad == 2) {
				iNextX = iX + 1;
				iNextY = iY;

			} else if (iQuad == 3) {
				iNextX = iX;
				iNextY = iY;

			} else {
				iNextX = iX;
				iNextY = iY + 1;
			}
		}

		// Calculate theta at the next X and Y intersection,
		// ensuring no intersections that cross panels.
		if ((nP1 == 4) && (dX1 * gridCS.GetEdgeX()[iNextX] < TINY)) {
			dNextLat1 = M_PI;
		} else if (fabs(gridCS.GetEdgeX()[iNextX]) < TINY) {
			dNextLat1 = M_PI;
		} else {
			dNextLat1 = dP * atan(sin(dLon) / gridCS.GetEdgeX()[iNextX]);
		}
		if ((nP1 == 4) && (dY1 * gridCS.GetEdgeX()[iNextY] < TINY)) {
			dNextLat2 = M_PI;
		} else if (fabs(gridCS.GetEdgeX()[iNextY]) < TINY) {
			dNextLat2 = M_PI;
		} else {
			dNextLat2 = -atan(cos(dLon) / gridCS.GetEdgeX()[iNextY]);
		}

		// Check for the end of this RLL element
		if ((dNextLat1 > gridRLL.GetLatEdge()[iLat+1]) &&
			(dNextLat2 > gridRLL.GetLatEdge()[iLat+1])
		) {
			CubedSphereTrans::XYPFromRLL(
				dLon, gridRLL.GetLatEdge()[iLat+1], dX2, dY2, nP2);

			if (nP2 != nP1) {
				CubedSphereTrans::XYPFromXYP(
					dX2, dY2, nP2, nP1, dX2, dY2);

				nP2 = nP1;
			}

#ifdef VERBOSE1
			printf("LONLINE1: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
				dX1, dY1, nP1, dX2, dY2, nP2);
#endif

			// Add the line segment weights
			AddLineSegment(
				gridCS, gridRLL, LineType_Lon,
				dX1, dY1, dX2, dY2, nP1,
				nLon, iLat, iX, iY);

			// Update the latitude index
			iLat++;

		// Intersection with a line of constant X first
		} else if (dNextLat1 < dNextLat2) {

			CubedSphereTrans::XYPFromRLL(
				dLon, dNextLat1, dX2, dY2, nP2);

			if (nP2 != nP1) {
				CubedSphereTrans::XYPFromXYP(
					dX2, dY2, nP2, nP1, dX2, dY2);

				nP2 = nP1;
			}

#ifdef VERBOSE1
			printf("LONLINE2: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
				dX1, dY1, nP1, dX2, dY2, nP2);
#endif

			// Add the line segment weights
			AddLineSegment(
				gridCS, gridRLL, LineType_Lon,
				dX1, dY1, dX2, dY2, nP1,
				nLon, iLat, iX, iY);

			// Update index
			if (nP1 == 4) {
				if ((iQuad == 1) || (iQuad == 2)) {
					iX = iX - 1;
				} else {
					iX = iX + 1;
				}
			} else {
				if ((iQuad == 1) || (iQuad == 2)) {
					iX = iX + 1;
				} else {
					iX = iX - 1;
				}
			}

		// Intersection with a line of constant Y first
		} else {
			CubedSphereTrans::XYPFromRLL(
				dLon, dNextLat2, dX2, dY2, nP2);

			if (nP2 != nP1) {
				CubedSphereTrans::XYPFromXYP(
					dX2, dY2, nP2, nP1, dX2, dY2);

				nP2 = nP1;
			}

#ifdef VERBOSE1
			printf("LONLINE3: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
				dX1, dY1, nP1, dX2, dY2, nP2);
#endif

			// Add the line segment weights
			AddLineSegment(
				gridCS, gridRLL, LineType_Lon,
				dX1, dY1, dX2, dY2, nP1,
				nLon, iLat, iX, iY);

			// Update index
			if (nP1 == 4) {
				if ((iQuad == 1) || (iQuad == 4)) {
					iY = iY + 1;
				} else {
					iY = iY - 1;
				}
			} else {
				if ((iQuad == 2) || (iQuad == 3)) {
					iY = iY - 1;
				} else {
					iY = iY + 1;
				}
			}
		}

		// Update position
		dX1 = dX2;
		dY1 = dY2;

		if ((nP1 == 4) && (fabs(dX1) < TINY) && (fabs(dY1) < TINY)) {
			break;
		}

		// Change to northern panel (equatorial panels have no weight for
		// longitude line segments and so are ignored)
		if (nP1 == 5) {
			// Determine direction of exit
			if (fabs(dX1 - 1.0) < TINY) {
				fPanelSwitch = true;
				dY1 = -dY1;

			} else if (fabs(dX1 + 1.0) < TINY) {
				fPanelSwitch = true;
				dY1 = -dY1;

			} else if (fabs(dY1 - 1.0) < TINY) {
				dY1 = -1.0;
				fPanelSwitch = true;

			} else if (fabs(dY1 + 1.0) < TINY) {
				dY1 = 1.0;
				fPanelSwitch = true;
			}

			// Switch to northern panel
			if (fPanelSwitch) {
				nP1 = 4;
				fPanelSwitch = false;

				if (fabs(dLon - 0.5 * M_PI) < TINY) {
					gridCS.IndexFromXY(dX1, TINY, iX, iY);
				} else if (fabs(dLon - 1.5 * M_PI) < TINY) {
					gridCS.IndexFromXY(dX1, -TINY, iX, iY);
				} else {
					gridCS.IndexFromXY(dX1, dY1, iX, iY);
				}

				// Determine new RLL index
				CubedSphereTrans::RLLFromXYP(dX1, dY1, nP1, dTemp, dLat);
				gridRLL.IndexFromRLL(dLon, dLat, nTemp, iLat);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void Interpolator::InitializeLatitudeSegment(
	const CubedSphereGrid &gridCS,
	const LatitudeLongitudeGrid &gridRLL,
	int nLat
) {
	int i;

	double dTemp;

	// Current latitude
	const double dLat = gridRLL.GetLatEdge()[nLat];

	// No line segments at either pole
	if (fabs(fabs(dLat) - 0.5 * M_PI) < TINY) {
		return;
	}

	// Current longitude
	double dLon;
	double dInitialLon;

	double dNextLon1;
	double dNextLon2;

	double dNextLon3;
	double dNextLon4;

	// Current XYP coordinates
	double dX1 = 0.0;
	double dY1 = 0.0;

	double dX2 = 0.0;
	double dY2 = 0.0;

	double dNextX1;
	double dNextX2;

	double dNextX3;
	double dNextX4;

	int nP1;
	int nP2;

	// Current cubed sphere element
	int iX;
	int iY;

	// Next intersection with a cubed sphere edge
	int iNextX;
	int iNextY;

	// Begin searching
	for (i = 0; i < gridRLL.GetLongitudes(); i++) {

		// Get current and next latitude
		dLon = gridRLL.GetLonEdge()[i];

		dInitialLon = dLon;

		// Get XYP coordinates of current position
		CubedSphereTrans::XYPFromRLL(dLon, dLat, dX1, dY1, nP1);

		// Immediate move up one panel if we're right at the boundary
		if (fabs(dY1 - 1.0) < TINY) {
			if ((nP1 < 4) && (dX1 > -TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 4, dX1, dY1);
				nP1 = 4;

			} else if ((nP1 == 4) && (dX1 > TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 2, dX1, dY1);
				nP1 = 2;

			} else if ((nP1 == 5) && (dX1 < -TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 0, dX1, dY1);
				nP1 = 0;
			}
		}
		if (fabs(dY1 + 1.0) < TINY) {
			if ((nP1 < 4) && (dX1 > -TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 5, dX1, dY1);
				nP1 = 5;

			} else if ((nP1 == 4) && (dX1 < -TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 0, dX1, dY1);
				nP1 = 0;

			} else if ((nP1 == 5) && (dX1 > TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 2, dX1, dY1);
				nP1 = 2;
			}
		}
		if (fabs(dX1 - 1.0) < TINY) {
			if (nP1 < 4) {
				CubedSphereTrans::XYPFromXYP(
					dX1, dY1, nP1, (nP1 + 1) % 4, dX1, dY1);
				nP1 = (nP1 + 1) % 4;

			} else if ((nP1 == 4) && (dY1 < -TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 1, dX1, dY1);
				nP1 = 1;

			} else if ((nP1 == 5) && (dY1 > TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 1, dX1, dY1);
				nP1 = 1;
			}
		}
		if (fabs(dX1 + 1.0) < TINY) {
			if ((nP1 == 4) && (dY1 > TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 3, dX1, dY1);
				nP1 = 3;

			} else if ((nP1 == 5) && (dY1 < -TINY)) {
				CubedSphereTrans::XYPFromXYP(dX1, dY1, nP1, 3, dX1, dY1);
				nP1 = 3;
			}
		}

		// Obtain the element index of the initial position
		if (dLat >= 0.0) {
			CubedSphereTrans::XYPFromRLL(
				dLon + TINY, dLat + TINY, dX2, dY2, nP2);

			gridCS.IndexFromXY(dX2, dY2, iX, iY);

		} else if (dLat < 0.0) {
			CubedSphereTrans::XYPFromRLL(
				dLon + TINY, dLat - TINY, dX2, dY2, nP2);

			gridCS.IndexFromXY(dX2, dY2, iX, iY);
		}

#ifdef VERBOSE1
		printf("BEGIN: (%1.4f, %1.4f, %d) - [%1.4f, %1.4f] x [%1.4f, %1.4f]\n",
			dX1, dY1, nP1,
			gridCS.GetEdgeX()[iX], gridCS.GetEdgeX()[iX+1],
			gridCS.GetEdgeX()[iY], gridCS.GetEdgeX()[iY+1]
		);
#endif

		for (;;) {
			// Equatorial panels
			if (nP1 < 4) {
				// There are two intersections for each line of latitude -
				// one with positive X and one with negative X.  Find all four
				// intersections, then take the intersection that is smallest,
				// but farther along than the current position.

				// Determine X coordinates of intersection with bottom
				if (fabs(dLat) < TINY) {
					dNextX1 = M_PI;
					dNextX2 = M_PI;
				} else {
					dNextX1 = gridCS.GetEdgeX()[iY] / tan(dLat);
					if (dNextX1 < 1.0 + TINY) {
						dNextX1 = M_PI;
						dNextX2 = M_PI;
					} else {
						dNextX1 = sqrt(dNextX1 * dNextX1 - 1.0);
						dNextX2 = -dNextX1;
					}
				}

				// Determine X coordinates of intersection with top
				if (fabs(dLat) < TINY) {
					dNextX3 = M_PI;
					dNextX4 = M_PI;
				} else {
					dNextX3 = gridCS.GetEdgeX()[iY+1] / tan(dLat);
					if (dNextX3 < 1.0 + TINY) {
						dNextX3 = M_PI;
						dNextX4 = M_PI;
					} else {
						dNextX3 = sqrt(dNextX3 * dNextX3 - 1.0);
						dNextX4 = -dNextX3;
					}
				}

				// Find the closest intersection between this latitude and
				// the bottom/top of the element
				if ((dNextX2 > dX1) && (dNextX2 < dNextX1)) {
					dNextX1 = dNextX2;
				}

				if ((dNextX4 > dX1) && (dNextX4 < dNextX3)) {
					dNextX2 = dNextX4;
				} else {
					dNextX2 = dNextX3;
				}

				// X coordinate of the right-edge of this CS element.
				dNextX3 = gridCS.GetEdgeX()[iX+1];

				// Calculate the X coordinate of the right-edge of this RLL
				// element.  Note that this formula doesn't take into account
				// problems associated with large RLL elements.
				dNextX4 = gridRLL.GetLonEdge()[i+1];
				dNextX4 = fmod(dNextX4 + 0.25 * M_PI, 0.5 * M_PI) - 0.25 * M_PI;
				dNextX4 = tan(dNextX4);

				//std::cout << dX1 << ", " << dNextX1 << ", " << iX << ", " << iY << std::endl;
				//printf("%1.14e : %1.14e\n", dLat, atan(dY1));

				// Find the closest major point of intersection
				if (dX1 - dNextX1 > -TINY) {
					dNextX1 = M_PI;
				}
				if (dX1 - dNextX2 > -TINY) {
					dNextX2 = M_PI;
				}
				if (dX1 - dNextX3 > TINY) {
					_EXCEPTIONT("Logic error.");
				}
				if (dX1 - dNextX4 > TINY) {
					dNextX4 = M_PI;
				}

				// First intersection is with the bottom of this CS element
				if ((dNextX1 < dNextX2) &&
					(dNextX1 < dNextX3) &&
					(dNextX1 < dNextX4)
				) {
					dX2 = dNextX1;
					dY2 = gridCS.GetEdgeX()[iY];
					nP2 = nP1;

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE1a: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					// Update index
					iY = iY - 1;

					// Update panel (go to the south panel)
					if (fabs(dY2 + 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP1, 5, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = 5;
					}

				// First intersection is with the top of this CS element
				} else if ((dNextX2 < dNextX3) && (dNextX2 < dNextX4)) {
					dX2 = dNextX2;
					dY2 = gridCS.GetEdgeX()[iY+1];
					nP2 = nP1;

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE1b: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					// Update index
					iY = iY + 1;

					// Update panel (go to the north panel)
					if (fabs(dY2 - 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP1, 4, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = 4;
					}

				// First intersection is with a line of constant X
				} else if (dNextX3 < dNextX4) {
					dX2 = dNextX3;
					dY2 = tan(dLat) * sqrt(1.0 + dX2 * dX2);
					nP2 = nP1;

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE2: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					// Update index
					iX = iX + 1;

					// Update panel (go to the next equatorial panel)
					if (fabs(dX2 - 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP1, (nP1 + 1) % 4, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = (nP1 + 1) % 4;
					}

				// First intersection is with the right side of the RLL cell
				} else {
					CubedSphereTrans::XYPFromRLL(
						gridRLL.GetLonEdge()[i+1], dLat, dX2, dY2, nP2);

					if (nP2 != nP1) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, nP1, dX2, dY2);

						nP2 = nP1;
					}

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE3: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					break;
				}

				// Update position
				dX1 = dX2;
				dY1 = dY2;
				nP1 = nP2;

			// Polar panels
			} else {

				// Panel indicator (for polar panels)
				double dP = (nP1 == 4)?(1.0):(-1.0);

				// Obtain the current longitude
				CubedSphereTrans::RLLFromXYP(dX1, dY1, nP1, dLon, dTemp);

				if (fabs(dTemp - dLat) > TINY) {
					_EXCEPTIONT("Logic error.");
				}

				// Current quadrant on the polar panels
				int iQuad;
				if (dLon < 0.5 * M_PI) {
					iQuad = 1;
				} else if (dLon < M_PI) {
					iQuad = 2;
				} else if (dLon < 1.5 * M_PI) {
					iQuad = 3;
				} else {
					iQuad = 4;
				}

				// Determine direction of movement through cubed sphere elements
				if (nP1 == 4) {
					if (iQuad == 1) {
						iNextX = iX + 1;
						iNextY = iY + 1;

					} else if (iQuad == 2) {
						iNextX = iX;
						iNextY = iY + 1;

					} else if (iQuad == 3) {
						iNextX = iX;
						iNextY = iY;

					} else {
						iNextX = iX + 1;
						iNextY = iY;
					}

				} else {
					if (iQuad == 1) {
						iNextX = iX + 1;
						iNextY = iY;

					} else if (iQuad == 2) {
						iNextX = iX;
						iNextY = iY;

					} else if (iQuad == 3) {
						iNextX = iX;
						iNextY = iY + 1;

					} else {
						iNextX = iX + 1;
						iNextY = iY + 1;
					}
				}

				// Calculate longitude at the next X and Y intersections,
				// and map onto the interval [lon, lon + 2 pi].  Note that there
				// are 4 possible intersections, given by dNextLon1, dNextLon2
				// (intersections with a line of constant X) and dNextLon3,
				// dNextLon4 (intersections with a line of constant Y
				dNextLon1 = dP * tan(dLat) * gridCS.GetEdgeX()[iNextX];

				if (fabs(dNextLon1) - 1.0 > -TINY) {
					dNextLon1 = 4.0 * M_PI;
					dNextLon2 = 4.0 * M_PI;
				} else {
					dNextLon1 = fmod(asin(dNextLon1) + 2.0 * M_PI, 2.0 * M_PI);
					if (dNextLon1 - dInitialLon < -TINY) {
						dNextLon1 += 2.0 * M_PI;
					}
					if (dNextLon1 - dInitialLon < -TINY) {
						dNextLon1 += 2.0 * M_PI;
					}

					dNextLon2 = fmod(3.0 * M_PI - dNextLon1, 2.0 * M_PI);
					if (dNextLon2 - dInitialLon < -TINY) {
						dNextLon2 += 2.0 * M_PI;
					}
					if (dNextLon2 - dInitialLon < -TINY) {
						dNextLon2 += 2.0 * M_PI;
					}
				}

				dNextLon3 = - tan(dLat) * gridCS.GetEdgeX()[iNextY];

				if (fabs(dNextLon3) - 1.0 > -TINY) {
					dNextLon3 = 4.0 * M_PI;
					dNextLon4 = 4.0 * M_PI;
				} else {
					dNextLon3 = fmod(acos(dNextLon3) + 2.0 * M_PI, 2.0 * M_PI);
					if (dNextLon3 - dInitialLon < -TINY) {
						dNextLon3 += 2.0 * M_PI;
					}
					if (dNextLon3 - dInitialLon < -TINY) {
						dNextLon3 += 2.0 * M_PI;
					}
					dNextLon4 = fmod(2.0 * M_PI - dNextLon3, 2.0 * M_PI);
					if (dNextLon4 - dInitialLon < -TINY) {
						dNextLon4 += 2.0 * M_PI;
					}
					if (dNextLon4 - dInitialLon < -TINY) {
						dNextLon4 += 2.0 * M_PI;
					}
				}

				// Determine the first intersection of each type X, Y
				if (dNextLon2 < dNextLon1) {
					dNextLon1 = dNextLon2;
				}
				if (dNextLon4 < dNextLon3) {
					dNextLon3 = dNextLon4;
				}

				if (dNextLon1 - dInitialLon < -TINY) {
					_EXCEPTION();
				}
				if (dNextLon3 - dInitialLon < -TINY) {
					printf("%1.5e %1.5e %1.5e %1.5e %1.5e\n",
						dLat,
						gridCS.GetEdgeX()[iNextY],
						dInitialLon,
						dNextLon1,
						dNextLon3);

					_EXCEPTION();
				}

				//std::cout << dLat << ", " << gridCS.GetEdgeX()[iNextX] << std::endl;
				//std::cout << dLon << " : " << dNextLon1 << ", " << dNextLon3 << ", " << gridRLL.GetLonEdge()[i+1] << std::endl;

				// Intersection with the right side of this RLL element
				double dNextLonEdge = gridRLL.GetLonEdge()[i+1];
				if (gridRLL.GetLonEdge()[i+1] < gridRLL.GetLonEdge()[i]) {
					dNextLonEdge += 2.0 * M_PI;
				}

				if ((dNextLonEdge - dNextLon1 < TINY) &&
					(dNextLonEdge - dNextLon3 < TINY)
				) {
					CubedSphereTrans::XYPFromRLL(
						gridRLL.GetLonEdge()[i+1], dLat, dX2, dY2, nP2);

					if (nP2 != nP1) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, nP1, dX2, dY2);

						nP2 = nP1;
					}

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE4: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					break;

				// Intersection with a line of constant X
				} else if (dNextLon1 < dNextLon3) {
					CubedSphereTrans::XYPFromRLL(
						dNextLon1, dLat, dX2, dY2, nP2);

					if (nP2 != nP1) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, nP1, dX2, dY2);

						nP2 = nP1;
					}

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE2: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					// Update index
					if (iNextX == iX) {
						iX = iX - 1;
					} else {
						iX = iX + 1;
					}

					// Switch panel ?
					if (fabs(dX2 - 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, 1, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = 1;

					} else if (fabs(dX2 + 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, 3, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = 3;
					}

				// Intersection with a line of constant Y
				} else {
					CubedSphereTrans::XYPFromRLL(
						dNextLon3, dLat, dX2, dY2, nP2);

					if (nP2 != nP1) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, nP1, dX2, dY2);

						nP2 = nP1;
					}

					// Add the line segment weights
					AddLineSegment(
						gridCS, gridRLL, LineType_Lat,
						dX1, dY1, dX2, dY2, nP1,
						i, nLat, iX, iY);

#ifdef VERBOSE1
					printf("LATLINE1: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
						dX1, dY1, nP1, dX2, dY2, nP2);
#endif

					// Update index
					if (iNextY == iY) {
						iY = iY - 1;
					} else {
						iY = iY + 1;
					}

					// Switch panel ?
					if (fabs(dY2 - dP * 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, 2, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = 2;

					} else if (fabs(dY2 + dP * 1.0) < TINY) {
						CubedSphereTrans::XYPFromXYP(
							dX2, dY2, nP2, 0, dX2, dY2);
						gridCS.IndexFromXY(dX2, dY2, iX, iY);
						nP2 = 0;
					}
				}

				// Update the position
				dX1 = dX2;
				dY1 = dY2;
				nP1 = nP2;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void Interpolator::InitializeBetaSegment(
	const CubedSphereGrid &gridCS,
	const LatitudeLongitudeGrid &gridRLL,
	const int nBeta,
	const int nP
) {
	// Current Y line
	const double dY = gridCS.GetEdgeX()[nBeta];
	if ((dY < -1.0 - TINY) || (dY > 1.0 + TINY)) {
		_EXCEPTIONT("Y must be defined on the range [-1, 1].");
	}

	// Current X coordinate
	double dX = -1.0;

	// Next X coordinates
	double dNextX1;
	double dNextX2;
	double dNextX3;
	double dNextX4;

	// Current cubed sphere element
	int iX;
	int iY;

	// Latlon coordinates
	double dLon;
	double dLat;
	double dRelLon;

	// Current latlon element
	int iLon;
	int iLat;

	// Temporary panel storage
	int nPtemp;

	// Get RLL coordinates of current position
	CubedSphereTrans::RLLFromXYP(dX, dY + TINY, nP, dLon, dLat);

	// Get the index of the current RLL element
	gridRLL.IndexFromRLL(dLon, dLat, iLon, iLat);

	if (iLon >= gridRLL.GetLongitudes()) {
		iLon = 0;
	}

	// Initial cubed sphere index
	iX = gridCS.GetInteriorBegin();
	iY = nBeta;

	// Equatorial panels
	if (nP < 4) {

		// Begin searching
		for (;;) {

			// Determine the intersection with the bottom and top of
			// this RLL element.  Note that there are two possible
			// intersections for each line of latitude.
			if (dY * gridRLL.GetLatEdge()[iLat] < TINY) {
				dNextX1 = 2.0;
				dNextX2 = 2.0;
			} else {
				dNextX1 = dY / tan(gridRLL.GetLatEdge()[iLat]);
				dNextX1 = dNextX1 * dNextX1 - 1.0;
				if (dNextX1 < TINY) {
					dNextX1 = 2.0;
					dNextX2 = 2.0;
				} else {
					dNextX1 = - sqrt(dNextX1);
					dNextX2 = - dNextX1;
				}
			}

			if (dY * gridRLL.GetLatEdge()[iLat+1] < TINY) {
				dNextX3 = 2.0;
				dNextX4 = 2.0;
			} else {
				dNextX3 = dY / tan(gridRLL.GetLatEdge()[iLat+1]);
				dNextX3 = dNextX3 * dNextX3 - 1.0;
				if (dNextX3 < TINY) {
					dNextX3 = 2.0;
					dNextX4 = 2.0;
				} else {
					dNextX3 = - sqrt(dNextX3);
					dNextX4 = - dNextX3;
				}
			}

			// Dispose of these intersections if they occur before the
			// current positions.  Store the resulting intersections in
			// dNextX1 (bottom) and dNextX2 (top).
			if (dNextX1 - dX < TINY) {
				if (dNextX2 - dX < TINY) {
					dNextX1 = 2.0;
				} else {
					dNextX1 = dNextX2;
				}
			}

			if (dNextX3 - dX < TINY) {
				if (dNextX4 - dX < TINY) {
					dNextX2 = 2.0;
				} else {
					dNextX2 = dNextX4;
				}
			} else {
				dNextX2 = dNextX3;
			}

			// Determine the intersection with the next longitude line,
			// and store in dNextX3.
			dRelLon = gridRLL.GetLonEdge()[iLon+1]
				- static_cast<double>(nP) * 0.5 * M_PI;

			dNextX3 = tan(dRelLon);

			// Store the end of this CS element in dNextX4
			dNextX4 = gridCS.GetEdgeX()[iX+1];

			//std::cout << dNextX1 << ", " << dNextX2 << ", " << dNextX3 << ", " << dNextX4 << std::endl;

			// We have reached the end of the CS element
			if ((dNextX4 < dNextX1) &&
				(dNextX4 < dNextX2) &&
				(dNextX4 < dNextX3)
			) {
				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX4, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE4: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX4, dY, nP);
#endif
				// Update the coordinates
				dX = dNextX4;

				// Update the index
				iX = iX + 1;

				if (iX == gridCS.GetInteriorEnd()) { 
					break;
				}

			// Intersection with a line of constant longitude
			} else if ((dNextX3 < dNextX1) && (dNextX3 < dNextX2)) {

				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX3, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE3: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX3, dY, nP);
#endif

				// Update the index
				iLon = iLon + 1;

				if (iLon == gridRLL.GetLongitudes()) {
					iLon = 0;
				}

				// Update the coordinates
				dX = dNextX3;

			// Intersection with a line of constant latitude (top)
			} else if (dNextX2 < dNextX1) {
			
				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX2, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE2: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX2, dY, nP);
#endif

				// Update the index
				iLat = iLat + 1;

				// Update the coordinates
				dX = dNextX2;

			// Intersection with a line of constant latitude (bottom)
			} else {

				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX1, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE1: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX1, dY, nP);
#endif

				// Update the index
				iLat = iLat - 1;

				// Update the coordinates
				dX = dNextX1;
			}
		}

	// Polar panels
	} else {
		// Panel indicator (for polar panels)
		double dP = (nP == 4)?(1.0):(-1.0);

		// Begin searching
		for (;;) {

			// Determine the intersection of this line with latitude lines
			if (fabs(fabs(gridRLL.GetLatEdge()[iLat]) - 0.5 * M_PI) < TINY) {
				if (fabs(dY) < TINY) {
					dNextX1 = 0.0;
					dNextX2 = 0.0;
				} else {
					dNextX1 = 2.0;
					dNextX2 = 2.0;
				}
			} else {
				dNextX1 = 1.0 / tan(gridRLL.GetLatEdge()[iLat]);
				dNextX1 = dNextX1 * dNextX1 - dY * dY;
				if (dNextX1 < TINY) {
					dNextX1 = 2.0;
					dNextX2 = 2.0;
				} else {
					dNextX1 = - sqrt(dNextX1);
					dNextX2 = - dNextX1;
				}
			}

			if (fabs(fabs(gridRLL.GetLatEdge()[iLat+1]) - 0.5 * M_PI) < TINY) {
				if (fabs(dY) < TINY) {
					dNextX3 = 0.0;
					dNextX4 = 0.0;
				} else {
					dNextX3 = 2.0;
					dNextX4 = 2.0;
				}
			} else {
				dNextX3 = 1.0 / tan(gridRLL.GetLatEdge()[iLat+1]);
				dNextX3 = dNextX3 * dNextX3 - dY * dY;
				if (dNextX3 < TINY) {
					dNextX3 = 2.0;
					dNextX4 = 2.0;
				} else {
					dNextX3 = - sqrt(dNextX3);
					dNextX4 = - dNextX3;
				}
			}


			// Dispose of these intersections if they occur before the
			// current positions.  Store the resulting intersections in
			// dNextX1 (bottom) and dNextX2 (top).
			if (dNextX1 - dX < TINY) {
				if (dNextX2 - dX < TINY) {
					dNextX1 = 2.0;
				} else {
					dNextX1 = dNextX2;
				}
			}

			if (dNextX3 - dX < TINY) {
				if (dNextX4 - dX < TINY) {
					dNextX2 = 2.0;
				} else {
					dNextX2 = dNextX4;
				}
			} else {
				dNextX2 = dNextX3;
			}

			// Determine the intersection with the next longitude line,
			// and store in dNextX3.
			if (fabs(dY) > TINY) {
				if (((nP == 4) && (dY < 0.0)) || ((nP == 5) && (dY > 0.0))) {
					dRelLon = gridRLL.GetLonEdge()[iLon+1];

					if ((dRelLon > 0.5 * M_PI - TINY) &&
						(dRelLon < 1.5 * M_PI + TINY)
					) {
						dNextX3 = 2.0;
					} else {
						dNextX3 = -dP * dY * tan(dRelLon);
					}

				} else {
					dRelLon = gridRLL.GetLonEdge()[iLon];

					if ((dRelLon < 0.5 * M_PI + TINY) ||
						(dRelLon > 1.5 * M_PI - TINY)
					) {
						dNextX3 = 2.0;
					} else {
						dNextX3 = -dP * dY * tan(dRelLon);
					}
				}
			} else {
				dNextX3 = 2.0;
			}

			if (dNextX3 < dX - TINY) {
				std::cout << dNextX3 << ", " << dX << std::endl;
				_EXCEPTION();
			}

			// Store the end of this CS element in dNextX4
			dNextX4 = gridCS.GetEdgeX()[iX+1];

			// We have reached the end of the CS element
			if ((dNextX4 < dNextX1) &&
				(dNextX4 < dNextX2) &&
				(dNextX4 < dNextX3)
			) {
				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX4, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE4: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX4, dY, nP);
#endif
				// Update the coordinates
				dX = dNextX4;

				// Update the index
				iX = iX + 1;

				if (iX == gridCS.GetInteriorEnd()) { 
					break;
				}

			// Intersection with a line of constant longitude
			} else if ((dNextX3 < dNextX1) && (dNextX3 < dNextX2)) {

				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX3, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE3: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX3, dY, nP);
#endif

				// Update the index
				if (((nP == 4) && (dY < 0.0)) || ((nP == 5) && (dY > 0.0))) {
					iLon = iLon + 1;
					if (iLon == gridRLL.GetLongitudes()) {
						iLon = 0;
					}
				} else {
					iLon = iLon - 1;
					if (iLon == 0) {
						iLon = gridRLL.GetLongitudes();
					}
				}

				// Update the coordinates
				dX = dNextX3;

			// Intersection with a line of constant latitude (top)
			} else if (dNextX2 < dNextX1) {
			
				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX2, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE2: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX2, dY, nP);
#endif

				// Update the index
				if (iLat != gridRLL.GetLatitudes()-1) {
					iLat = iLat + 1;
				} else {
					CubedSphereTrans::RLLFromXYP(1.0, TINY, nP, dLon, dLat);
					gridRLL.IndexFromRLL(dLon, dLat, iLon, iLat);

					iLat = gridRLL.GetLatitudes()-1;
					dLat = 0.5 * M_PI;
				}

				// Update the coordinates
				dX = dNextX2;

			// Intersection with a line of constant latitude (bottom)
			} else {

				// Add the line segment weights
				AddLineSegment(
					gridCS, gridRLL, LineType_Beta,
					dX, dY, dNextX1, dY, nP,
					iLon, iLat, iX, nBeta);

#ifdef VERBOSE1
				printf("YLINE1: (%1.4f, %1.4f, %d) - (%1.4f, %1.4f, %d)\n",
					dX, dY, nP, dNextX1, dY, nP);
#endif

				// Update the index
				if (iLat != 0) {
					iLat = iLat - 1;
				} else {
					CubedSphereTrans::RLLFromXYP(1.0, TINY, nP, dLon, dLat);
					gridRLL.IndexFromRLL(dLon, dLat, iLon, iLat);

					iLat = 0;
					dLat = -0.5 * M_PI;
				}

				// Update the coordinates
				dX = dNextX1;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

