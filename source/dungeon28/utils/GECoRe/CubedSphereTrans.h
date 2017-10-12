///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereTrans.h
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

#ifndef _CUBEDSPHERETRANS_H_
#define _CUBEDSPHERETRANS_H_

#include "Exception.h"

#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

class CubedSphereTrans {

private:
	///	<summary>
	///		No initialization for objects of this class.
	///	</summary>
	CubedSphereTrans()
	{ }

public:
	///	<summary>
	///		Determine the Cartesian (XYZ) coordinates of a point on a sphere
	///		given its Gnomonic (XYP) coordinate.
	///	</summary>
	///	<parameters>
	///		XX       - Gnomonic X coordinate
	///		YY       - Gnomonic Y coordinate
	///		np       - Panel coordinate (0-5, 4 = north, 5 = south)
	///		xx (OUT) - Calculated X coordinate
	///		yy (OUT) - Calculated Y coordinate
	///		zz (OUT) - Calculated Z coordinate
	///	</parameters>
	static void XYZFromXYP(
		double dXX,
		double dYY,
		int np,
		double &xx,
		double &yy,
		double &zz
	);

	///	<summary>
	///		Determine the (X, Y, idest) coordinates of a source point on 
	///		panel isource.
	///	</summary>
	///	<notes>
	///		It is safe to use this function with the same inputs (dX_in, dY_in)
	///		as outputs (dX_out, dY_out).
	///	</notes>
	///	<parameters>
	///		dX_in           - Gnomonic X coordinate on source panel
	///		dY_in           - Gnomonic Y coordinate on source panel
	///		isource         - Source panel (0-5, 4 = north, 5 = south)
	///		idest           - Destination panel (0-5, 4 = north, 5 = south)
	///		dX_out (OUT)    - Gnomonic X coordinate on destination panel
	///		dY_out (OUT)    - Gnomonic Y coordinate on destination panel
	///	</parameters>
	static void XYPFromXYP(
		double dX_in,
		double dY_in,
		int isource,
		int idest,
		double &dX_out,
		double &dY_out
	);

	///	<summary>
	///		Determine the (lat, lon) coordinates of a point in gnomonic XYP
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		dX     - Gnomonic X coordinate (IN)
	///		dY     - Gnomonic Y coordinate (IN)
	///		nP     - Gnomonic panel coordinate (0-5, 4 = north, 5 = south) (IN)
	///		lat    - Latitude (on the interval [-pi/2, pi/2]) (OUT)
	///		lon    - Longitude (on the interval [0, 2 pi]) (OUT)
	///	</parameters>
	static void RLLFromXYP(
		double dX,
		double dY,
		int nP,
		double &lon,
		double &lat
	);

	///	<summary>
	///		Determine the gnomonic XYP coordinates of a point in (lat, lon)
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		lat    - Latitude (IN)
	///		lon    - Longitude (IN)
	///		dX     - Gnomonic X coordinate (OUT)
	///		dY     - Gnomonic Y coordinate (OUT)
	///		nP     - Gnomonic panel coordinate (0-5, 4 = north, 5 = south) (OUT)
	///	</parameters>
	static void XYPFromRLL(
		double lon,
		double lat,
		double &dX,
		double &dY,
		int &nP
	);

	///	<summary>
	///		Translate a vector in spherical coordinates to a vector in
	///		equiangular coordinates.  The components of the vector field in
	///		spherical coordinates are expected in the unit basis, whereas the
	///		components in equiangular coordinates are given in the geometric
	///		basis.
	///	</summary>
	///	<parameters>
	///		dX      - Gnomonic X coordinate of translation point (IN)
	///		dY      - Gnomonic Y coordinate of translation point (IN)
	///		nP      - Gnomonic panel coordinate (0-5, 4 = north, 5 = south) of
	///		          translation point (IN)
	///		dUlon   - Longitudinal component of vector field (IN)
	///		dUlat   - Latitudinal component of vector field (IN)
	///		dUalpha - Alpha component of vector field (OUT)
	///		dUbeta  - Beta component of vector field (OUT)
	///	</parameters>
	static void VecTransABPFromRLL(
		double dX,
		double dY,
		int nP,
		double dUlon,
		double dUlat,
		double &dUalpha,
		double &dUbeta
	);

	///	<summary>
	///		Translate a vector in equiangular coordinates to a vector in
	///		spherical coordinates.  The components of the vector field in
	///		spherical coordinates are expected in the unit basis, whereas the
	///		components in equiangular coordinates are given in the geometric
	///		basis.
	///	</summary>
	///	<parameters>
	///		dX      - Gnomonic X coordinate of translation point (IN)
	///		dY      - Gnomonic Y coordinate of translation point (IN)
	///		nP      - Gnomonic panel coordinate (0-5, 4 = north, 5 = south) of
	///		          translation point (IN)
	///		dUalpha - Alpha component of vector field (IN)
	///		dUbeta  - Beta component of vector field (IN)
	///		dUlon   - Longitudinal component of vector field (OUT)
	///		dUlat   - Latitudinal component of vector field (OUT)
	///	</parameters>
	static void VecTransRLLFromABP(
		double dX,
		double dY,
		int nP,
		double dUalpha,
		double dUbeta,
		double &dUlon,
		double &dUlat
	);

	///	<summary>
	///		Convert equiangular derivatives to gnomonic derivatives.
	///	</summary>
	///	<parameters>
	///		X - Gnomonic X coordinate (IN)
	///		Y - Gnomonic Y coordinate (IN)
	///		DxU - First alpha / X derivative (IN/OUT)
	///		DyU - First beta / Y derivative (IN/OUT)
	///		DxxU - Second alpha / X derivative (IN/OUT)
	///		DxyU - First cross derivative derivative (IN/OUT)
	///		DyyU - Second beta / Y derivative (IN/OUT)
	///	</parameters>
	static void DerivTransXYPFromABP(
		double dX,
		double dY,
		double &dDxU,
		double &dDyU,
		double &dDxxU,
		double &dDxyU,
		double &dDyyU
	);

	///	<summary>
	///		Compute the angle between grid lines on the gnomonic projection.
	///	</summary>
	///	<parameters>
	///		X - Gnomonic X coordinate
	///		Y - Gnomonic Y coordinate
	///	</parameters>
	inline static double GnomonicGridAngle(
		double dX,
		double dY
	) {
		return acos(- dX * dY / (sqrt(1.0 + dX * dX) * sqrt(1.0 + dY * dY)));
	}

	///	<summary>
	///		Compute the area of a square region on the Gnomonic cubed sphere
	///		grid.
	///	</summary>
	///	<parameters>
	///		dX      - Gnomonic X coordinate of bottom-left corner
	///		dDeltaX - Delta X
	///		dY      - Gnomonic Y coordinate of the bottom-left corner
	///		dDeltaY - Delta Y
	///		dRadius - Radius of the sphere
	///	</parameters>
	inline static double GnomonicElementArea(
		double dX,
		double dDeltaX,
		double dY,
		double dDeltaY,
		double dRadius = 1.0
	) {
		return dRadius * dRadius *
		       (+ GnomonicGridAngle(dX          , dY   )
	    	    - GnomonicGridAngle(dX + dDeltaX, dY   )
				- GnomonicGridAngle(dX          , dY + dDeltaY)
				+ GnomonicGridAngle(dX + dDeltaX, dY + dDeltaY));

	}

	///	<summary>
	///		Calculate the length of a spherical arc along a line of
	///		constant X, in Gnomonic coordinates.
	///	</summary>
	///	<parameters>
	///		dX      - X coordinate of line
	///		dY1     - Y coordinate of initial point
	///		dY2     - Y coordinate of final point
	///		dRadius - Radius of the sphere
	///	</parameters>
	inline static double GnomonicXLineLength(
		double dX,
		double dY1,
		double dY2,
		double dRadius = 1.0
	) {
		return dRadius * (
	    	     + atan(dY2 / sqrt(1.0 + dX * dX))
				 - atan(dY1 / sqrt(1.0 + dX * dX)));
	}

	///	<summary>
	///		Calculate the length of a spherical arc along a line of
	///		constant Y, in Gnomonic coordinates.
	///	</summary>
	///	<parameters>
	///		dY      - Y coordinate of line
	///		dX1     - X coordinate of initial point
	///		dX2     - X coordinate of final point
	///		dRadius - Radius of the sphere
	///	</parameters>
	inline static double GnomonicYLineLength(
		double dY,
		double dX1,
		double dX2,
		double dRadius = 1.0
	) {
		return dRadius * (
		         + atan(dX2 / sqrt(1.0 + dY * dY))
	    	     - atan(dX1 / sqrt(1.0 + dY * dY)));
	}

	///	<summary>
	///		Determine the panel id in the given direction relative to
	///		a given panel id.
	///	</summary>
	///	<parameters>
	///		Nc      - Resolution of the cubed sphere grid (IN)
	///		p_src   - Source panel (IN)
	///		ix_src  - Index in Y direction of source element (IN)
	///		jx_src  - Index in X direction of source element (IN)
	///		p_dest  - Destination panel (OUT)
	///		ix_dest - Index in the Y direction of destination element (OUT)
	///		jx_dest - Index in the X direction of destination element (OUT)
	///	</parameters>
	static void RelativeCoord(
		int Nc,
		int p_src,
		int ix_src,
		int jx_src,
		int &p_dest,
		int &ix_dest,
		int &jx_dest,
		bool &switchAB,
		bool &switchDir
	) {
		// Internal points
		if ((ix_src >= 0) && (jx_src >= 0) && (ix_src < Nc) && (jx_src < Nc)) {
			ix_dest = ix_src;
			jx_dest = jx_src;
			p_dest = p_src;
			switchAB = false;
			switchDir = false;

		// Equatorial panel 0
		} else if (p_src == 0) {
			if (jx_src >= Nc) {
				p_dest = 1;
				ix_dest = ix_src;
				jx_dest = jx_src - Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src >= Nc) {
				p_dest = 4;
				ix_dest = ix_src - Nc;
				jx_dest = jx_src;
				switchAB = false;
				switchDir = false;

			} else if (jx_src < 0) {
				p_dest = 3;
				ix_dest = ix_src;
				jx_dest = jx_src + Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src < 0) {
				p_dest = 5;
				ix_dest = ix_src + Nc;
				jx_dest = jx_src;
				switchAB = false;
				switchDir = false;
			}

		// Equatorial panel 1
		} else if (p_src == 1) {
			if (jx_src >= Nc) {
				p_dest = 2;
				ix_dest = ix_src;
				jx_dest = jx_src - Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src >= Nc) {
				p_dest = 4;
				ix_dest = jx_src;
				jx_dest = 2 * Nc - 1 - ix_src;
				switchAB = true;
				switchDir = true;

			} else if (jx_src < 0) {
				p_dest = 0;
				ix_dest = ix_src;
				jx_dest = jx_src + Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src < 0) {
				p_dest = 5;
				ix_dest = Nc - 1 - jx_src;
				jx_dest = Nc + ix_src;
				switchAB = true;
				switchDir = false;
			}

		// Equatorial panel 2
		} else if (p_src == 2) {
			if (jx_src >= Nc) {
				p_dest = 3;
				ix_dest = ix_src;
				jx_dest = jx_src - Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src >= Nc) {
				p_dest = 4;
				ix_dest = 2 * Nc - 1 - ix_src;
				jx_dest = Nc - jx_src - 1;
				switchAB = false;
				switchDir = true;

			} else if (jx_src < 0) {
				p_dest = 1;
				ix_dest = ix_src;
				jx_dest = jx_src + Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src < 0) {
				p_dest = 5;
				ix_dest = - ix_src - 1;
				jx_dest = Nc - jx_src - 1;
				switchAB = false;
				switchDir = true;
			}

		// Equatorial panel 3
		} else if (p_src == 3) {
			if (jx_src >= Nc) {
				p_dest = 0;
				ix_dest = ix_src;
				jx_dest = jx_src - Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src >= Nc) {
				p_dest = 4;
				ix_dest = Nc - jx_src - 1;
				jx_dest = ix_src - Nc;
				switchAB = true;
				switchDir = false;

			} else if (jx_src < 0) {
				p_dest = 2;
				ix_dest = ix_src;
				jx_dest = jx_src + Nc;
				switchAB = false;
				switchDir = false;

			} else if (ix_src < 0) {
				p_dest = 5;
				ix_dest = jx_src;
				jx_dest = - ix_src - 1;
				switchAB = true;
				switchDir = true;
			}

		// North polar panel
		} else if (p_src == 4) {
			if (jx_src >= Nc) {
				p_dest = 1;
				ix_dest = 2 * Nc - 1 - jx_src;
				jx_dest = ix_src;
				switchAB = true;
				switchDir = true;

			} else if (ix_src >= Nc) {
				p_dest = 2;
				ix_dest = 2 * Nc - 1 - ix_src;
				jx_dest = Nc - 1 - jx_src;
				switchAB = false;
				switchDir = true;

			} else if (jx_src < 0) {
				p_dest = 3;
				ix_dest = jx_src + Nc;
				jx_dest = Nc - 1 - ix_src;
				switchAB = true;
				switchDir = false;

			} else if (ix_src < 0) {
				p_dest = 0;
				ix_dest = ix_src + Nc;
				jx_dest = jx_src;
				switchAB = false;
				switchDir = false;
			}

		// South polar panel
		} else if (p_src == 5) {
			if (jx_src >= Nc) {
				p_dest = 1;
				ix_dest = jx_src - Nc;
				jx_dest = Nc - 1 - ix_src;
				switchAB = true;
				switchDir = false;

			} else if (ix_src >= Nc) {
				p_dest = 0;
				ix_dest = ix_src - Nc;
				jx_dest = jx_src;
				switchAB = false;
				switchDir = false;

			} else if (jx_src < 0) {
				p_dest = 3;
				ix_dest = - jx_src - 1;
				jx_dest = ix_src;
				switchAB = true;
				switchDir = true;

			} else if (ix_src < 0) {
				p_dest = 2;
				ix_dest = - ix_src - 1;
				jx_dest = Nc - 1 - jx_src;
				switchAB = false;
				switchDir = true;
			}
		}
	}

	///	<summary>
	///		Remap vectors from the given source panel to the given destination
	///		panel.  All calculations are performed in-place.
	///	</summary>
	static void VecPanelTrans(
		int p_src,
		int p_dest,
		double &dAlpha,
		double &dBeta,
		double dXdest,
		double dYdest
	) {
		if ((p_dest < 4) && (p_src < 4)) {
			if ((p_dest + 7 - p_src) % 4 == 0) {
				VecToP2FromP1(dAlpha, dBeta, dXdest, dYdest);
			} else if ((p_dest + 5 - p_src) % 4 == 0) {
				VecToP1FromP2(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_src == 4) {
			if (p_dest == 0) {
				VecToP1FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 1) {
				VecToP2FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 2) {
				VecToP3FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 3) {
				VecToP4FromP5(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_dest == 4) {
			if (p_src == 0) {
				VecToP5FromP1(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 1) {
				VecToP5FromP2(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 2) {
				VecToP5FromP3(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 3) {
				VecToP5FromP4(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_src == 5) {
			if (p_dest == 0) {
				VecToP1FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 1) {
				VecToP2FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 2) {
				VecToP3FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_dest == 3) {
				VecToP4FromP6(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else if (p_dest == 5) {
			if (p_src == 0) {
				VecToP6FromP1(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 1) {
				VecToP6FromP2(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 2) {
				VecToP6FromP3(dAlpha, dBeta, dXdest, dYdest);
			} else if (p_src == 3) {
				VecToP6FromP4(dAlpha, dBeta, dXdest, dYdest);
			} else {
				_EXCEPTION();
			}

		} else {
			_EXCEPTION();
		}
	}

	///	<summary>
	///		Remap vectors to equatorial panel 2 from equatorial panel 1.  All
	///		calculations are performed in-place.
	///	</summary>
	///	<param name="dAlpha">
	///		Alpha component of the vector system.
	///	</param>
	///	<param name="dBeta">
	///		Beta component of the vector system.
	///	</param>
	///	<param name="dX2">
	///		Gnomonic X coordinate (panel 2) of point of rotation.
	///	</param>
	///	<param name="dY2">
	///		Gnomonic Y coordinate (panel 2) of point of rotation.
	///	</param>
	inline static void VecToP2FromP1(
		double &dAlpha,
		double &dBeta,
		double dX2,
		double dY2
	) {
		dBeta =
			dY2 * (1.0 + dX2 * dX2) / (dX2 * (1.0 + dY2 * dY2)) * dAlpha
			- (dX2 * dX2 + dY2 * dY2) / (dX2 * (1.0 + dY2 * dY2)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 1 from
	///		equatorial panel 2.
	///	</summary>
	inline static void VecToP1FromP2(
		double &dAlpha,
		double &dBeta,
		double dX1,
		double dY1
	) {
		dBeta =
			dY1 * (1.0 + dX1 * dX1) / (dX1 * (1.0 + dY1 * dY1)) * dAlpha
			+ (dX1 * dX1 + dY1 * dY1) / (dX1 * (1.0 + dY1 * dY1)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 1.
	///	</summary>
	inline static void VecToP5FromP1(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		dAlpha =
			- (dX5 * dX5 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dAlpha
			+ dX5 * (1.0 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 1
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP1FromP5(
		double &dAlpha,
		double &dBeta,
		double dX1,
		double dY1
	) {
		dAlpha =
			(dX1 * dX1 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dAlpha
			+ dX1 * (1.0 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 2.
	///	</summary>
	inline static void VecToP5FromP2(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		double dTemp = dBeta;

		dBeta =
			(dX5 * dX5 + dY5 * dY5) / (dX5 * (1.0 + dY5 * dY5)) * dAlpha
			- (dY5 * (dX5 * dX5 + 1.0)) / (dX5 * (1.0 + dY5 * dY5)) * dBeta;

		dAlpha = -dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 2
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP2FromP5(
		double &dAlpha,
		double &dBeta,
		double dX2,
		double dY2
	) {
		double dTemp = dAlpha;

		dAlpha =
			- (dX2 * (1.0 + dY2 * dY2)) / (dY2 * (1.0 + dX2 * dX2)) * dAlpha
			+ (dX2 * dX2 + dY2 * dY2) / (dY2 * (1.0 + dX2 * dX2)) * dBeta;

		dBeta = -dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 3.
	///	</summary>
	inline static void VecToP5FromP3(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		dAlpha =
			- (dX5 * dX5 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dAlpha
			- dX5 * (1.0 + dY5 * dY5) / (dY5 * (dX5 * dX5 + 1.0)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 3
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP3FromP5(
		double &dAlpha,
		double &dBeta,
		double dX3,
		double dY3
	) {
		dAlpha =
			- (dX3 * dX3 + dY3 * dY3) / (dY3 * (dX3 * dX3 + 1.0)) * dAlpha
			- dX3 * (dY3 * dY3 + 1.0) / (dY3 * (dX3 * dX3 + 1.0)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 5 (north pole)
	///		from equatorial panel 4.
	///	</summary>
	inline static void VecToP5FromP4(
		double &dAlpha,
		double &dBeta,
		double dX5,
		double dY5
	) {
		double dTemp = dBeta;

		dBeta =
			(dX5 * dX5 + dY5 * dY5) / (dX5 * (1.0 + dY5 * dY5)) * dAlpha
			+ dY5 * (dX5 * dX5 + 1.0) / (dX5 * (1.0 + dY5 * dY5)) * dBeta;

		dAlpha = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 4
	///		from panel 5 (north pole).
	///	</summary>
	inline static void VecToP4FromP5(
		double &dAlpha,
		double &dBeta,
		double dX4,
		double dY4
	) {
		double dTemp = dAlpha;

		dAlpha =
			(dX4 * (1.0 + dY4 * dY4)) / (dY4 * (1.0 + dX4 * dX4)) * dAlpha
			- (dX4 * dX4 + dY4 * dY4) / (dY4 * (1.0 + dX4 * dX4)) * dBeta;

		dBeta = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 1.
	///	</summary>
	inline static void VecToP6FromP1(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		dAlpha =
			(dX6 * dX6 + dY6 * dY6) / (dY6 * (1.0 + dX6 * dX6)) * dAlpha
			+ dX6 * (1.0 + dY6 * dY6) / (dY6 * (1.0 + dX6 * dX6)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 1
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP1FromP6(
		double &dAlpha,
		double &dBeta,
		double dX1,
		double dY1
	) {
		dAlpha =
			- (dX1 * dX1 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dAlpha
			+ dX1 * (1.0 + dY1 * dY1) / (dY1 * (1.0 + dX1 * dX1)) * dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 2.
	///	</summary>
	inline static void VecToP6FromP2(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		double dTemp = dBeta;

		dBeta =
			- (dX6 * dX6 + dY6 * dY6) / (dX6 * (1.0 + dY6 * dY6)) * dAlpha
			+ (dY6 * (1.0 + dX6 * dX6)) / (dX6 * (1.0 + dY6 * dY6)) * dBeta;

		dAlpha = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 2
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP2FromP6(
		double &dAlpha,
		double &dBeta,
		double dX2,
		double dY2
	) {
		double dTemp = dAlpha;
		
		dAlpha =
			(dX2 * (1.0 + dY2 * dY2)) / (dY2 * (1.0 + dX2 * dX2)) * dAlpha
			+ (dX2 * dX2 + dY2 * dY2) / (dY2 * (1.0 + dX2 * dX2)) * dBeta;

		dBeta = dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 3.
	///	</summary>
	inline static void VecToP6FromP3(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		dAlpha =
			(dX6 * dX6 + dY6 * dY6) / (dY6 * (1.0 + dX6 * dX6)) * dAlpha
			- (dX6 * (1.0 + dY6 * dY6)) / (dY6 * (1.0 + dX6 * dX6)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 3
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP3FromP6(
		double &dAlpha,
		double &dBeta,
		double dX3,
		double dY3
	) {
		dAlpha =
			(dX3 * dX3 + dY3 * dY3) / (dY3 * (dX3 * dX3 + 1.0)) * dAlpha
			- dX3 * (dY3 * dY3 + 1.0) / (dY3 * (dX3 * dX3 + 1.0)) * dBeta;

		dBeta = -dBeta;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to panel 6 (south pole)
	///		from equatorial panel 4.
	///	</summary>
	inline static void VecToP6FromP4(
		double &dAlpha,
		double &dBeta,
		double dX6,
		double dY6
	) {
		double dTemp = dBeta;

		dBeta =
			- (dX6 * dX6 + dY6 * dY6) / (dX6 * (1.0 + dY6 * dY6)) * dAlpha
			- (dY6 * (1.0 + dX6 * dX6)) / (dX6 * (1.0 + dY6 * dY6)) * dBeta;

		dAlpha = -dTemp;
	}

	///	<summary>
	///		As VecToP2FromP1 except remapping vectors to equatorial panel 4
	///		from panel 6 (south pole).
	///	</summary>
	inline static void VecToP4FromP6(
		double &dAlpha,
		double &dBeta,
		double dX4,
		double dY4
	) {
		double dTemp = dAlpha;

		dAlpha =
			- (dX4 * (1.0 + dY4 * dY4)) / (dY4 * (1.0 + dX4 * dX4)) * dAlpha
			- (dY4 * dY4 + dX4 * dX4) / (dY4 * (1.0 + dX4 * dX4)) * dBeta;

		dBeta = -dTemp;
	}
};

////////////////////////////////////////////////////////////////////////////////

#endif

