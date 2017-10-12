///////////////////////////////////////////////////////////////////////////////
///
///	\file    Interpolator.h
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

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "CubedSphereGrid.h"
#include "LatitudeLongitudeGrid.h"
#include "SystemState.h"

#include "Exception.h"
#include "DataMatrix3D.h"
#include "MathHelper.h"

#include <vector>
#include <cfloat>

////////////////////////////////////////////////////////////////////////////////

class Interpolator {
	protected:
		///	<summary>
		///		A small number, used in the search algorithm for determining
		///		when two floating point values are equal.
		///	</summary>
		static const double TINY;

		///	<summary>
		///		An enumeration of line types.
		///	</summary>
		enum LineType {
			LineType_Alpha = 0,
			LineType_Beta,
			LineType_Lat,
			LineType_Lon
		};

		///	<summary>
		///		A class for storing weights associated with a RLL element
		///		and CS element pair.
		///	</summary>
		class ElementWeight {
			public:
				///	<summary>
				///		Constructor.
				///	</summary>
				ElementWeight() :
					ixX(0),
					ixY(0),
					ixP(0)
				{ }

				///	<summary>
				///		Initializer.
				///	</summary>
				void Initialize(int nWeights) {
					int i;

					dWeight.resize(nWeights);

					for (i = 0; i < nWeights; i++) {
						dWeight[i] = 0.0;
					}
				}

			public:
				///	<summary>
				///		Coordinates of the associated CS element.
				///	</summary>
				int ixX;
				int ixY;
				int ixP;

				///	<summary>
				///		An array of weights associated with the reconstruction.
				///	</summary>
				std::vector<double> dWeight;
		};

		///	<summary>
		///		A class defining a vector of ElementWeights.
		///	</summary>
		class ElementWeightVector :
			public std::vector<ElementWeight>
		{
			public:
			///	<summary>
			///		Add a bond between this latitude-longitude element and the
			///		given cubed sphere element, with given weights.
			///	</summary>
			void AddWeight(
				int ixX,
				int ixY,
				int ixP,
				double *dWeight,
				int nWeights,
				double dDirection
			) {
				int i;
				int j;

				// Update an existing element weight
				for (i = 0; i < size(); i++) {
					if (((*this)[i].ixX == ixX) &&
						((*this)[i].ixY == ixY) &&
						((*this)[i].ixP == ixP)
					) {
						for (j = 0; j < nWeights; j++) {
							(*this)[i].dWeight[j] += dDirection * dWeight[j];
						}
						break;
					}
				}

				// Add a new element weight
				if (i == size()) {
					resize(size() + 1);

					(*this)[i].Initialize(nWeights);

					(*this)[i].ixX = ixX;
					(*this)[i].ixY = ixY;
					(*this)[i].ixP = ixP;
					for (j = 0; j < nWeights; j++) {
						(*this)[i].dWeight[j] += dDirection * dWeight[j];
					}
				}
			}
		};

		///	<summary>
		///		A class for storing weights associated with a RLL element.
		///	</summary>
		class RLLElementWeights {
			public:
				///	<summary>
				///		Constructor.
				///	</summary>
				RLLElementWeights() :
					dWeights(NULL)
				{ }

				///	<summary>
				///		Initializer.
				///	</summary>
				void Initialize(
					int nLon,
					int nLat
				) {
					int i;
					int j;

					// Allocate the weights array
					dWeights = new ElementWeightVector*[nLon];
					for (i = 0; i < nLon; i++) {
						dWeights[i] = new ElementWeightVector[nLat];
					}

					// Store size of resulting array
					m_nLon = nLon;
					m_nLat = nLat;
				}

				///	<summary>
				///		Destructor.
				///	</summary>
				~RLLElementWeights() {
					int i;

					if (dWeights != NULL) {
						for (i = 0; i < m_nLon; i++) {
							delete[] dWeights[i];
						}
						delete[] dWeights;
					}
				}

				///	<summary>
				///		Accessor operator.
				///	</summary>
				inline operator ElementWeightVector **() {
					return dWeights;
				}

			protected:
				///	<summary>
				///		The size of the weights array.
				///	</summary>
				int m_nLon;
				int m_nLat;

				///	<summary>
				///		An array of weights associated with each RLL element.
				///	</summary>
				ElementWeightVector ** dWeights;
		};

	public:
		///	<summary>
		///		Default constructor.
		///	</summary>
		Interpolator() :
			m_nOrder(0),
			m_nWeights(0)
		{ }

		///	<summary>
		///		Constructor.
		///	</summary>
		Interpolator(
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
		) = 0;

	protected:
		///	<summary>
		///		Calculate the weight assigned to the given line segment and
		///		update the element weights array.
		///	</summary>
		void AddLineSegment(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL,
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
		);

		///	<summary>
		///		Obtain line segments along a line of constant longitude.
		///	</summary>
		void InitializeLongitudeSegment(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL,
			int nLon
		);

		///	<summary>
		///		Obtain line segments along a line of constant latitude.
		///	</summary>
		void InitializeLatitudeSegment(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL,
			int nLon
		);

		///	<summary>
		///		Obtain line segments along a line of constant beta.
		///	</summary>
		void InitializeBetaSegment(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL,
			int nBeta,
			int nP
		);

		///	<summary>
		///		Initialize the interpolation weights in each element.
		///	</summary>
		void InitializeWeights(
			const CubedSphereGrid &gridCS,
			const LatitudeLongitudeGrid &gridRLL
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

		///	<summary>
		///		Remap from the given cubed sphere grid to a RLL grid.
		///	</summary>
		void RemapCStoRLL(
			const CubedSphereGrid & gridCS,
			const LatitudeLongitudeGrid & gridRLL,
			int nInVar,
			DataMatrix<double> & dataRLL,
			bool fUseRemappedArea = false
		);

	public:
		///	<summary>
		///		Crop extreme values, with the consequence of lost conservation.
		///	</summary>
		void CropExtremeValues(
			const LatitudeLongitudeGrid & gridRLL,
			DataMatrix<double> & dataRLL,
			double dMin,
			double dMax
		);

	private:
		///	<summary>
		///		An object for storing weight data associated with each
		///		RLL grid element.
		///	</summary>
		RLLElementWeights weights;

		///	<summary>
		///		"Middle" index for lines of constant beta.
		///	</summary>
		int m_jMiddle;

	protected:
		///	<summary>
		///		The order of this remapping scheme.
		///	</summary>
		int m_nOrder;

		///	<summary>
		///		The number of weights associated with the order of this scheme.
		///	</summary>
		int m_nWeights;
};

////////////////////////////////////////////////////////////////////////////////

#endif

