///////////////////////////////////////////////////////////////////////////////
///
///	\file    LatitudeLongitudeGrid.h
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

#ifndef _LATITUDELONGITUDEGRID_H_
#define _LATITUDELONGITUDEGRID_H_

///////////////////////////////////////////////////////////////////////////////

#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"

#include "CubedSphereTrans.h"

///////////////////////////////////////////////////////////////////////////////

class LatitudeLongitudeGrid {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	///	<parameters>
	///		nLongitudes - Number of longitude elements
	///		nLatitudes - Number of latitude elements
	///	</parameters>
	LatitudeLongitudeGrid(
		int nLatitudes,
		int nLongitudes,
		double dLongitudeShift = 0.0
	) {
		Initialize(nLatitudes, nLongitudes, dLongitudeShift);
	}

public:
	///	<summary>
	///		Initialize this latitude-longitude grid.
	///	</summary>
	///	<parameters>
	///		nLongitudes - Number of longitude elements
	///		nLatitudes - Number of latitude elements
	///	</parameters>
	void Initialize(
		int nLatitudes,
		int nLongitudes,
		double dLongitudeShift = 0.0
	);

	///	<summary>
	///		Set the remapped grid areas obtained from interpolation.
	///	</summary>
	void SetRemappedGridArea(
		const DataMatrix<double> & matRemappedElementArea
	);

public:
	///	<summary>
	///		Determine the element index from the given RLL coordinates.
	///	</summary>
	///	<parameters>
	///		dLon - Longitude coordinate (IN)
	///		dLat - Latitude coordinate (IN)
	///		iLon - Longitude index (OUT)
	///		iLat - Latitude index (OUT)
	///	</parameters>
	void IndexFromRLL(
		double dLon,
		double dLat,
		int &iLon,
		int &iLat
	) const;

public:
	///	<summary>
	///		Accessors.
	///	</summary>
	inline int GetLongitudes() const {
		return m_nLongitudes;
	}

	inline int GetLatitudes() const {
		return m_nLatitudes;
	}

	inline const DataMatrix<double>& GetElementArea() const {
		return m_matElementArea;
	}

	inline const DataMatrix<double>& GetRemappedElementArea() const {
		return m_matRemappedElementArea;
	}

	inline const DataMatrix<double>& GetEdgeLength() const {
		_EXCEPTION();
		//return m_matEdgeLength;
	}

	inline const DataMatrix3D<double>& GetGeoCentroidX() const {
		_EXCEPTION();
		//return m_matGeoCentroidX;
	}

	inline const DataVector<double>& GetLonCentroid() const {
		return m_matLonCentroid;
	}

	inline const DataVector<double>& GetLatCentroid() const {
		return m_matLatCentroid;
	}

	inline const DataVector<double>& GetLonEdge() const {
		return m_matLonEdge;
	}

	inline const DataVector<double>& GetLatEdge() const {
		return m_matLatEdge;
	}

private:
	///	<summary>
	///		Number of elements in each direction.
	///	</summary>
	int m_nLongitudes;

	///	<summary>
	///		Number of ghost elements on each side.
	///	</summary>
	int m_nLatitudes;

	///	<summary>
	///		Longitude shift.
	///	</summary>
	double m_dLongitudeShift;

	///	<summary>
	///		Calculated per-element area.
	///	</summary>
	DataMatrix<double> m_matElementArea;

	///	<summary>
	///		Calculated remapped per-element area, obtained by remapping the
	///		unit field to this grid.
	///	</summary>
	DataMatrix<double> m_matRemappedElementArea;

/*
	///	<summary>
	///		Calculated line segment length.
	///	</summary>
	DataMatrix<double> m_matEdgeLength;

	///	<summary>
	///		Calculated element geometric centroids.
	///	</summary>
	DataMatrix3D<double> m_matGeoCentroidX;
*/

	///	<summary>
	///		Calculated edge coordinates for lines of constant longitude.
	///	</summary>
	DataVector<double> m_matLonEdge;

	///	<summary>
	///		Calculated edge coordinates for lines of constant latitude.
	///	</summary>
	DataVector<double> m_matLatEdge;

	///	<summary>
	///		Calculated centroid coordinates for lines of constant longitude.
	///	</summary>
	DataVector<double> m_matLonCentroid;

	///	<summary>
	///		Calculated centroid coordinates for lines of constant latitude.
	///	</summary>
	DataVector<double> m_matLatCentroid;

};

///////////////////////////////////////////////////////////////////////////////

#endif

