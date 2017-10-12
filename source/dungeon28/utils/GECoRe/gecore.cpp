///////////////////////////////////////////////////////////////////////////////
///
///	\file    gecore.cpp
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

#include "CommandLine.h"
#include "DataMatrix.h"
#include "Preferences.h"
#include "MathHelper.h"
#include "StringHelper.h"
#include "STLStringHelper.h"

#include "CubedSphereGrid.h"
#include "LatitudeLongitudeGrid.h"
#include "HOMMEInterpolator.h"

#include "HOMMEGrid.h"

#include "netcdfcpp.h"

#include <cfloat>
#include <iostream>
#include <memory>

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	try {

		// Order of accuracy (currently hardcoded to 4)
		const int nDefaultNp = 4;

		// Default output file
		const char * szDefaultOutputFile = "outRemapped.nc";

		// Default output file
		const char * szDefaultVariables = "*";

		// Input file
		std::string strInputFile;

		// Output file
		std::string strOutputFile = szDefaultOutputFile;

		// List of variables to remap
		std::string strVariables = szDefaultVariables;

		// Preferences file
		std::string strPreferences;

		// Weights file
		std::string strWeights;

		// Number of latitudes
		int nLatitudes = 360;

		// Number of longitudes
		int nLongitudes = 180;

		// Number of GLL nodes per cubed sphere element
		int nNp = nDefaultNp;

		// Parse the command line
		int argfirst = 1;
		if ((argc >= 2) && (argv[1][0] != '-')) {
			argfirst = 2;
			strPreferences = argv[1];
		}

		// Open the preferences file, if specified
		Preferences aPrefs;

		if (strPreferences != "") {
			aPrefs.ParsePreferences(strPreferences.c_str());

			// Number of latitudes from preferences file
			int nPrefLat = aPrefs.GetPreferenceAsInt_NoThrow("Latitudes");
			if (nPrefLat != 0) {
				nLatitudes = nPrefLat;
			}

			// Number of longitudes from preferences file
			int nPrefLon = aPrefs.GetPreferenceAsInt_NoThrow("Longitudes");
			if (nPrefLon != 0) {
				nLongitudes = nPrefLon;
			}

			// Input filename
			std::string strPrefInputFile =
				aPrefs.GetPreferenceAsString_NoThrow("InputFile");

			if (strPrefInputFile != "") {
				if (strInputFile != "") {
					_EXCEPTIONT("\n  Input file specified on command line"
						" and in Preferences file.");
				}
				strInputFile = strPrefInputFile;
			}

			// Output filename
			std::string strPrefOutputFile =
				aPrefs.GetPreferenceAsString_NoThrow("OutputFile");

			if (strPrefOutputFile != "") {
				if (strOutputFile != szDefaultOutputFile) {
					_EXCEPTIONT("\n  Output file specified on command line"
						" and in Preferences file.");
				}
				strOutputFile = strPrefOutputFile;
			}

			// Variable list
			std::string strPrefVariables =
				aPrefs.GetPreferenceAsString_NoThrow("Variables");

			if (strPrefVariables != "") {
				if (strVariables != "*") {
					_EXCEPTIONT("\n  Variables list specified on command line"
						" and in Preferences file.");
				}
				strVariables = strPrefVariables;
			}
		}

		if (argc == 1) {
			Announce("Usage: %s <Preferences File> <Parameters>", argv[0]);
			Announce("       %s <Parameters>\n", argv[0]);
		}

		BeginCommandLine()
			CommandLineString(strInputFile, "infile", strInputFile);
			CommandLineString(strOutputFile, "outfile", strOutputFile);
			CommandLineString(strVariables, "vars", strVariables);
			CommandLineString(strWeights, "weightsout", strWeights);
			CommandLineInt(nLatitudes, "lat", nLatitudes);
			CommandLineInt(nLongitudes, "lon", nLongitudes);
			CommandLineInt(nNp, "np", nNp);

			ParseCommandLine(argc, argv, argfirst);
		EndCommandLine(argv)

		if (argc == 1) {
			return (0);
		}

		// Check input file
		if (strInputFile == "") {
			return (0);
		}

		// System state
		int nResolution;
		int nGhostElements;

		// Load NetCDF file
		AnnounceBanner();
		AnnounceStartBlock("Loading input file");
		Announce("File: \"%s\"", strInputFile.c_str());

		NcFile ncdf_file(strInputFile.c_str());

		// Set error mode
		NcError err(NcError::verbose_nonfatal);

		// Determine dimensions of the given variable
		NcVar * var;
		NcDim * dim;

		long nTime = 1;
		bool fHasTime = false;

		dim = ncdf_file.get_dim("time");
		if (dim != NULL) {
			nTime = dim->size();
			fHasTime = true;
		} else {
			Announce("NOTE: Time dimension not detected");
		}

		dim = ncdf_file.get_dim("lev");
		if (dim == NULL) {
			_EXCEPTIONT("NetCDF file has no dimension \"lev\".");
		}
		long nLevel = dim->size();

		dim = ncdf_file.get_dim("ncol");
		if (dim == NULL) {
			_EXCEPTIONT("NetCDF file has no dimension \"ncol\".");
		}
		long nCol = dim->size();

		Announce("#Time %i / #Lev %i / #Col %i", nTime, nLevel, nCol);

		AnnounceEndBlock();

		// Initialize the cubed sphere grid
		AnnounceStartBlock("Initializing cubed sphere grid");

		std::auto_ptr<CubedSphereGrid> pGridCS;

		pGridCS.reset(new HOMMEGrid(ncdf_file, nNp));

		AnnounceEndBlock();

		// Initialize the latitude longitude grid
		AnnounceStartBlock("Initializing RLL grid");

		// Longitudinal shift on the output
		double dLongitudeShift =
			aPrefs.GetPreferenceAsDouble_NoThrow("OutputLongitudeShift");

		if ((dLongitudeShift < 0.0) || (dLongitudeShift > 360.0)) {
			_EXCEPTION1("OutputLongitudeShift must be in [0,360]. "
				"Given: %f", dLongitudeShift);
		}

		if (dLongitudeShift != 0.0) {
			Announce("Longitude shift: %1.5e", dLongitudeShift);
		}

		// Instantiate a RLL grid
		LatitudeLongitudeGrid gridRLL(
			nLatitudes, nLongitudes, dLongitudeShift * M_PI / 180.0);

		AnnounceEndBlock();

		// Initialize the interpolator
		AnnounceStartBlock("Initializing GECoRe");

		std::auto_ptr<Interpolator> pInterp;

		pInterp.reset(new HOMMEInterpolator(*pGridCS, gridRLL, nNp));

		AnnounceEndBlock("Done");

		// Output from remapping process
		DataMatrix<double> dataRLL;
		dataRLL.Initialize(nLatitudes, nLongitudes);

		// Temporary storage from netCDF file
		DataVector<float> dataTempFloat;
		dataTempFloat.Initialize(nCol);

		DataVector<double> dataTempDouble;
		dataTempDouble.Initialize(nCol);
/*
		// Build the reference grid and reset RLL areas
		AnnounceStartBlock("Computing reference grid areas");
		for (int i = 0; i < dataTempDouble.GetRows(); i++) {
			dataTempDouble[i] = 1.0;
		}
		pGridCS->ReconstructLowOrderMonotone(dataTempDouble);
		pInterp->RemapCStoRLL(*pGridCS, gridRLL, 0, dataRLL);
		double dOriginalArea = 0.0;
		double dRemappedArea = 0.0;
		for (int i = 0; i < dataRLL.GetRows(); i++) {
		for (int j = 0; j < dataRLL.GetColumns(); j++) {
			dOriginalArea += gridRLL.GetElementArea()[i][j];
			dataRLL[i][j] *= gridRLL.GetElementArea()[i][j];
			dRemappedArea += dataRLL[i][j];
		}
		}
		Announce("Original area: %1.15e", dOriginalArea);
		Announce("Remapped area: %1.15e", dRemappedArea);
		gridRLL.SetRemappedGridArea(dataRLL);
		AnnounceEndBlock("Done");
*/
		// Initialize NetCDF file
		AnnounceStartBlock("Initializing netCDF file");

		// Output the results to a netCDF file
		NcFile ncdf_out(
			strOutputFile.c_str(),
			NcFile::Replace,
			NULL, 0,
			NcFile::Offset64Bits);

		// Dimensions
		Announce("Dimensions");
		NcDim * dimTime;
		NcDim * dimLon = ncdf_out.add_dim("lon", nLongitudes);
		NcDim * dimLat = ncdf_out.add_dim("lat", nLatitudes);
		NcDim * dimLev = ncdf_out.add_dim("lev", nLevel);
		NcDim * dimILev = ncdf_out.add_dim("ilev", nLevel + 1);

		if (fHasTime) {
			dimTime = ncdf_out.add_dim("time");
		}

		// Buffer storage
		DataVector<double> dTemp;
		dTemp.Initialize(
			Max(static_cast<long int>(Max(nTime, nLevel + 1)),
				static_cast<long int>(Max(nLongitudes, nLatitudes))));

		// Latitudes
		Announce("Latitudes");
		for (int i = 0; i < nLatitudes; i++) {
			dTemp[i] = 180.0 / static_cast<double>(nLatitudes)
				* (static_cast<double>(i) + 0.5) - 90.0;
		}
		NcVar * latVar = ncdf_out.add_var("lat", ncDouble, dimLat);
		latVar->put(dTemp, nLatitudes);
		latVar->add_att("long_name", "latitude");
		latVar->add_att("units", "degrees_north");

		// Longitudes
		Announce("Longitudes");

		double dDeltaLongitude = 360.0 / static_cast<double>(nLongitudes);
		double dTruncatedLongitudeShift =
			fmod(dLongitudeShift, dDeltaLongitude);

		for (int i = 0; i < nLongitudes; i++) {
			dTemp[i] = dDeltaLongitude * (static_cast<double>(i) + 0.5)
				+ dTruncatedLongitudeShift;
		}

		NcVar * lonVar = ncdf_out.add_var("lon", ncDouble, dimLon);
		lonVar->put(dTemp, nLongitudes);
		lonVar->add_att("long_name", "longitude");
		lonVar->add_att("units", "degrees_east");

		// Levels
		Announce("Levels");
		NcVar * levVar = ncdf_file.get_var("lev");
		
		if (levVar != NULL) {
			levVar->get(dTemp, nLevel);

			levVar = ncdf_out.add_var("lev", ncDouble, dimLev);
			levVar->put(dTemp, nLevel);
			levVar->add_att("long_name", "hybrid level at midpoints (1000*(A+B))");
			levVar->add_att("units", "level");
			levVar->add_att("positive", "down");
			levVar->add_att("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate");
			levVar->add_att("lev:formula_terms", "a: hyam b: hybm p0: P0 ps: PS");
		} else {
			Announce("WARNING: Variable \"lev\" not found");
		}

		// Interface levels
		Announce("Interface Levels");
		NcVar * ilevVar = ncdf_file.get_var("ilev");

		if (ilevVar != NULL) {
			ilevVar->get(dTemp, nLevel+1);

			ilevVar = ncdf_out.add_var("ilev", ncDouble, dimILev);
			ilevVar->put(dTemp, nLevel + 1);
			ilevVar->add_att("long_name", "hybrid level at interfaces (1000*(A+B))");
			ilevVar->add_att("units", "level");
			ilevVar->add_att("positive", "down");
			ilevVar->add_att("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate");
			ilevVar->add_att("lev:formula_terms", "a: hyai b: hybi p0: P0 ps: PS");
		} else {
			Announce("WARNING: Variable \"ilev\" not found");
		}

		// Interface level hybrid coefficients
		Announce("Interface level hybrid coefficients");
		NcVar * hyaiVar = ncdf_file.get_var("hyai");
		
		if (hyaiVar != NULL) {
			hyaiVar->get(dTemp, nLevel+1);

			hyaiVar = ncdf_out.add_var("hyai", ncDouble, dimILev);
			hyaiVar->put(dTemp, nLevel + 1);
			hyaiVar->add_att("long_name", "hybrid A coefficient at layer interfaces");
		} else {
			Announce("WARNING: Variable \"hyai\" not found");
		}

		NcVar *hybiVar = ncdf_file.get_var("hybi");

		if (hybiVar != NULL) {
			hybiVar->get(dTemp, nLevel+1);

			hybiVar = ncdf_out.add_var("hybi", ncDouble, dimILev);
			hybiVar->put(dTemp, nLevel + 1);
			hybiVar->add_att("long_name", "hybrid B coefficient at layer interfaces");
		} else {
			Announce("WARNING: Variable \"hybi\" not found");
		}

		// Model level hybrid coefficients
		Announce("Model level hybrid coefficients");
		NcVar * hyamVar = ncdf_file.get_var("hyam");

		if (hyamVar != NULL) {
			hyamVar->get(dTemp, nLevel);

			hyamVar = ncdf_out.add_var("hyam", ncDouble, dimLev);
			hyamVar->put(dTemp, nLevel);
			hyamVar->add_att("long_name", "hybrid A coefficient at layer midpoints");
		} else {
			Announce("WARNING: Variable \"hyam\" not found");
		}
		
		NcVar * hybmVar = ncdf_file.get_var("hybm");
		
		if (hybmVar != NULL) {
			hybmVar->get(dTemp, nLevel);

			hybmVar = ncdf_out.add_var("hybm", ncDouble, dimLev);
			hybmVar->put(dTemp, nLevel);
			hybmVar->add_att("long_name", "hybrid B coefficient at layer midpoints");
		} else {
			Announce("WARNING: Variable \"hybm\" not found");
		}

		// Reference pressure
		Announce("Reference pressure");

		NcVar * p0Var = ncdf_file.get_var("P0");
		
		if (p0Var != NULL) {
			p0Var->get(dTemp, 1);

			p0Var = ncdf_out.add_var("P0", ncDouble);
			p0Var->put(dTemp, 1);
			p0Var->add_att("long_name", "reference pressure");
			p0Var->add_att("units", "Pa");
		} else {
			Announce("WARNING: Variable \"P0\" not found");
		}

		// Times
		Announce("Times");

		NcVar * timeVar = ncdf_file.get_var("time");
		
		if ((fHasTime) && (timeVar != NULL)) {
			timeVar->get(dTemp, nTime);

			timeVar = ncdf_out.add_var("time", ncDouble, dimTime);
			timeVar->put(dTemp, nTime);
			timeVar->add_att("long_name", "time");
			timeVar->add_att("unit", "days since 0005-09-01 00:00:00");
			timeVar->add_att("calender", "noleap");
			timeVar->add_att("bounds", "time_bnds");
		} else {
			Announce("WARNING: Variable \"time\" not found");
		}

		AnnounceEndBlock();

		// Build variable list
		if (strVariables == "*") {
			AnnounceStartBlock("Building variable list");
			_EXCEPTIONT("UNIMPLEMENTED");
			AnnounceEndBlock("Done");
		}

		// Begin remapping
		bool fFirstRemap = true;
		AnnounceBanner("BEGIN REMAPPING");

		const char * szPos1 = strVariables.c_str();
		const char * szPos2;

		char szVariable[50];

		for (int v = 0;; v++) {

			bool fHasLevels = false;
			bool fIsILev = false;

			// Parse
			szPos2 = StringHelper::FindWhitespace(szPos1);

			strncpy(szVariable, szPos1, szPos2 - szPos1);
			szVariable[szPos2 - szPos1] = '\0';

			if ((*szPos1 == '\0') || (*szPos1 == '\n') || (*szPos1 == '\r')) {
				break;
			}

			szPos1 = StringHelper::IgnoreWhitespace(szPos2);

			// Data
			var = ncdf_file.get_var(szVariable);
			if (var == NULL) {
				_EXCEPTION1("netCDF file does not contain variable \"%s\"",
					szVariable);
			}

			// Find special instructions for this variable
			std::string strSpecial =
				aPrefs.GetPreferenceAsString_NoThrow(szVariable);

			if (strSpecial != "") {
				STLStringHelper::ToLower(strSpecial);
			}

			// Rename longitude
			if (strcmp(szVariable, "lon") == 0) {
				strcat(szVariable, "x");
			}
			if (strcmp(szVariable, "lat") == 0) {
				strcat(szVariable, "x");
			}

			// Check dimensionality
			NcDim *firstDim = var->get_dim(0);
			if (firstDim == NULL) {
				_EXCEPTIONT("Cannot remap variable with no spatial dimension.");

			} else if (strcmp(firstDim->name(), "lev") == 0) {
				fHasTime = false;
				fHasLevels = true;
				nLevel = firstDim->size();

			} else if (strcmp(firstDim->name(), "ilev") == 0) {
				fHasTime = false;
				fHasLevels = true;
				fIsILev = true;
				nLevel = firstDim->size();

			} else if (strcmp(firstDim->name(), "ncol") == 0) {
				fHasTime = false;
				fHasLevels = false;
				nLevel = 1;

			} else if (strcmp(firstDim->name(), "time") == 0) {
				fHasTime = true;

				NcDim *secondDim = var->get_dim(1);
				if (secondDim == NULL) {
					_EXCEPTIONT("Cannot remap variable with no spatial dimension.");

				} else if (strcmp(secondDim->name(), "lev") == 0) {
					fHasLevels = true;
					nLevel = secondDim->size();

				} else if (strcmp(secondDim->name(), "ilev") == 0) {
					fHasLevels = true;
					fIsILev = true;
					nLevel = secondDim->size();

				} else if (strcmp(secondDim->name(), "ncol") == 0) {
					fHasLevels = false;
					nLevel = 1;

				} else {
					_EXCEPTIONT("Unknown second dimension:  "
						"Expected \"lev\" or \"ncol\".");
				}

			} else {
				_EXCEPTIONT("Unknown first dimension:  "
					"Expected \"time\", \"lev\" or \"ncol\".");
			}

			// Define the output variable
			NcVar *varOut;
			if (fHasTime) {
				if (fHasLevels) {
					if (fIsILev) {
						varOut = ncdf_out.add_var(szVariable, ncDouble,
							dimTime, dimILev, dimLat, dimLon);
					} else {
						varOut = ncdf_out.add_var(szVariable, ncDouble,
							dimTime, dimLev, dimLat, dimLon);
					}
				} else {
					varOut = ncdf_out.add_var(szVariable, ncDouble,
						dimTime, dimLat, dimLon);
				}
			} else {
				if (fHasLevels) {
					if (fIsILev) {
						varOut = ncdf_out.add_var(szVariable, ncDouble,
							dimILev, dimLat, dimLon);
					} else {
						varOut = ncdf_out.add_var(szVariable, ncDouble,
							dimLev, dimLat, dimLon);
					}
				} else {
					varOut = ncdf_out.add_var(szVariable, ncDouble,
						dimLat, dimLon);
				}
			}

			for (int i = 0; i < var->num_atts(); i++) {
				NcAtt *attA = var->get_att(i);
				varOut->add_att(attA->name(), attA->values()->as_string(0));
			}

			// Check for special instructions
			bool fMalformedSpecial = false;
			bool fLowOrderMonotone = false;
			bool fCropData = false;
			double dCropMin = 0.0;
			double dCropMax = 1.0;

			// No special instruction
			if (strSpecial == "") {

			// Monotone remapping
			} else if ((strSpecial.length() >= 4) &&
				(strncmp(strSpecial.c_str(), "mono", 4) == 0)
			) {
				fLowOrderMonotone = true;

			// Crop extrema
			} else if ((strSpecial.length() >= 4) &&
				(strncmp(strSpecial.c_str(), "crop", 4) == 0)
			) {
				if (strSpecial.length() > 4) {

					strSpecial[4] = '\0';
					const char * szSpecialPos =
						StringHelper::FindWhitespace(strSpecial.c_str() + 5);
					strSpecial[szSpecialPos - strSpecial.c_str()] = '\0';
					dCropMin = atof(strSpecial.c_str() + 5);
					dCropMax = atof(szSpecialPos+1);
				} else {
					fMalformedSpecial = true;
				}

			// Malformed special
			} else {
				fMalformedSpecial = true;
			}

		// Get the values array
		for (int t = 0; t < nTime; t++) {
		for (int l = 0; l < nLevel; l++) {

			if (fFirstRemap) {
				fFirstRemap = false;
			} else {
				AnnounceBanner();
			}

			AnnounceStartBlock("Remapping variable %s (t = %i, l = %i)",
				var->name(), t, l);

			// Unknown special instruction
			if (fMalformedSpecial) {
				Announce("MALFORMED SPECIAL: %s (IGNORED)",
					strSpecial.c_str());
			}

			// Load in data
			if (fHasTime) {
				if (fHasLevels) {
					var->set_cur(t, l, 0);

					if (var->type() == ncFloat) {
						var->get(dataTempFloat, 1, 1, nCol);
					} else if (var->type() == ncDouble) {
						var->get(dataTempDouble, 1, 1, nCol);
					}
				} else {
					var->set_cur(t, 0);

					if (var->type() == ncFloat) {
						var->get(dataTempFloat, 1, nCol);
					} else if (var->type() == ncDouble) {
						var->get(dataTempDouble, 1, nCol);
					}
				}
			} else {
				if (fHasLevels) {
					var->set_cur(l, 0);

					if (var->type() == ncFloat) {
						var->get(dataTempFloat, 1, nCol);
					} else if (var->type() == ncDouble) {
						var->get(dataTempDouble, 1, nCol);
					}
				} else {
					if (var->type() == ncFloat) {
						var->get(dataTempFloat, nCol);
					} else if (var->type() == ncDouble) {
						var->get(dataTempDouble, nCol);
					}
				}
			}

			// Copy data into double precision array
			if (var->type() == ncFloat) {
				AnnounceStartBlock("Converting floats to doubles");
				for (int i = 0; i < dataTempFloat.GetRows(); i++) {
					dataTempDouble[i] = static_cast<double>(dataTempFloat[i]);
				}
				AnnounceEndBlock("Done");
			}

			// Perform reconstruction
			AnnounceStartBlock("Performing reconstruction");
			if (fLowOrderMonotone) {
				Announce("SPECIAL: Use low-order monotone reconstruction");
				pGridCS->ReconstructLowOrderMonotone(dataTempDouble);
			} else {
				pGridCS->ReconstructHighOrder(dataTempDouble);
			}
			AnnounceEndBlock("Done");

			// Perform remapping
			AnnounceStartBlock("Performing remapping");
			pInterp->RemapCStoRLL(
				*pGridCS, gridRLL, 0, dataRLL);
			AnnounceEndBlock("Done");

			// No special instruction
			if (fCropData) {
				Announce("SPECIAL: Crop to %1.10e / %1.10e",
					dCropMin, dCropMax);

				pInterp->CropExtremeValues(
					gridRLL, dataRLL, dCropMin, dCropMax);
			}

			// Output the mass checksum on the cubed sphere grid
			double dSumCS;
			double dMinCS;
			double dMaxCS;

			pGridCS->Checksum(dataTempDouble, dSumCS, dMinCS, dMaxCS);

			Announce(" CS Min: %1.8e / Max: %1.8e / Sum: %1.14e",
				dMinCS, dMaxCS, dSumCS);

			// Output the mass checksum on the RLL grid
			double dSumRLL = 0.0;
			double dMaxRLL = -DBL_MAX;
			double dMinRLL =  DBL_MAX;

			int iMinLon;
			int iMinLat;

			for (int i = 0; i < gridRLL.GetLatitudes(); i++) {
			for (int j = 0; j < gridRLL.GetLongitudes(); j++) {
				dSumRLL += dataRLL[i][j] * gridRLL.GetElementArea()[i][j];

				if (dMaxRLL < dataRLL[i][j]) {
					dMaxRLL = dataRLL[i][j];
				}
				if (dMinRLL > dataRLL[i][j]) {
					dMinRLL = dataRLL[i][j];
					iMinLon = j;
					iMinLat = i;
				}
			}
			}

			Announce("RLL Min: %1.8e / Max: %1.8e / Sum: %1.14e",
				dMinRLL, dMaxRLL, dSumRLL);

			Announce("RLL Min LonLat: %i %i", iMinLon, iMinLat);

			// Output the mass difference
			Announce("Mass difference: %1.14e", (dSumCS - dSumRLL));

			// Output the results to a netCDF file
			AnnounceStartBlock("Writing results");

			if (fHasTime) {
				if (fHasLevels) {
					varOut->set_cur(t, l, 0, 0);
					varOut->put(&(dataRLL[0][0]), 1, 1, nLatitudes, nLongitudes);
				} else {
					varOut->set_cur(t, 0, 0);
					varOut->put(&(dataRLL[0][0]), 1, nLatitudes, nLongitudes);
				}
			} else {
				if (fHasLevels) {
					varOut->set_cur(l, 0, 0);
					varOut->put(&(dataRLL[0][0]), 1, nLatitudes, nLongitudes);
				} else {
					varOut->set_cur(0, 0);
					varOut->put(&(dataRLL[0][0]), nLatitudes, nLongitudes);
				}
			}
			AnnounceEndBlock("Done");

			// Done remapping this time / level
			AnnounceEndBlock("Done");
		}
		}
		}

		AnnounceBanner("END REMAPPING");

		// Return successfully
		return (0);

	// Exception
	} catch (Exception &e) {
		std::cout << std::endl << e.ToString() << std::endl;
		return (-1);
	}
}

////////////////////////////////////////////////////////////////////////////////

