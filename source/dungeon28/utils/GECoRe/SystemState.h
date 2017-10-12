///////////////////////////////////////////////////////////////////////////////
///
///	\file    SystemState.h
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

#ifndef _SYSTEMSTATE_H_
#define _SYSTEMSTATE_H_

#include "CubedSphereGrid.h"
#include "DataMatrix3D.h"

///////////////////////////////////////////////////////////////////////////////

class SystemState {

public:
	///	<summary>
	///		The data storage array associated with one rectangular grid.
	///	</summary>
	typedef DataMatrix3D<double> StaticStateData;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	SystemState() {
	}

	///	<summary>
	///		Constructor.
	///	</summary>
	SystemState(
		const CubedSphereGrid &grid,
		int nComponents
	) {
		int i;

		// Allocate states
		for (i = 0; i < 6; i++) {
			data[i].Initialize(
				grid.GetTotalElements(),
				grid.GetTotalElements(),
				nComponents);
		}
	}

public:
	///	<summary>
	///		Initializer.
	///	</summary>
	void Initialize(
		const CubedSphereGrid &grid,
		int nComponents
	) {
		int i;

		// Allocate states
		for (i = 0; i < 6; i++) {
			data[i].Initialize(
				grid.GetTotalElements(),
				grid.GetTotalElements(),
				nComponents);
		}
	}

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	SystemState & operator=(const SystemState &state) {
		int i;

		for (i = 0; i < 6; i++) {
			data[i] = state.data[i];
		}

		return (*this);
	}

	///	<summary>
	///		Zero operator.
	///	</summary>
	void Zero() {
		int i;

		for (i = 0; i < 6; i++) {
			data[i].Zero();
		}
	}

public:
	///	<summary>
	///		Accessors.
	///	</summary>
	inline StaticStateData & operator[](int n) {
		return data[n];
	}

	inline const StaticStateData & operator[](int n) const {
		return data[n];
	}

private:
	///	<summary>
	///		An array of StaticGridStates, one corresponding to each active
	///		region when performing simulations.
	///	</summary>
	StaticStateData data[6];
};

///////////////////////////////////////////////////////////////////////////////

#endif

