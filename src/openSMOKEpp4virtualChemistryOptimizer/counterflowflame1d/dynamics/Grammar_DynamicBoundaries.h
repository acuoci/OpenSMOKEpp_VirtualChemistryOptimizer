/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Grammar_DynamicBoundaries_H
#define OpenSMOKE_Grammar_DynamicBoundaries_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_DynamicBoundaries : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Type of dynamic boundaries: TEMPERATURE_VELOCITY_SLOPES | FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE | FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE | FIXED_BALANCE_FUEL_VELOCITY_SLOPE | FIXED_BALANCE_OX_VELOCITY_SLOPE | SIN_VELOCITY_INPHASE",
				true));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SlopeFuelTemperature",
				OpenSMOKE::SINGLE_MEASURE,
				"Slope of fuel temperature (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SlopeOxidizerTemperature",
				OpenSMOKE::SINGLE_MEASURE,
				"Slope of oxidizer temperature (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SlopeOxidizerVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Slope of oxidizer velocity (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SlopeFuelVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Slope of oxidizer velocity (default: 0)",
				false));


			// ------------------------------------------------------------------------------------------------
			// Output options section
			// ------------------------------------------------------------------------------------------------

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EndTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Total time of simulation (default: 1e4 s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CompactOutput",
				OpenSMOKE::SINGLE_BOOL,
				"Compact output on file (only most relevant data)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StepsOutput",
				OpenSMOKE::SINGLE_DOUBLE,
				"Frequency of output (default: 5)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SnapshotTimes",
				OpenSMOKE::VECTOR_STRING,
				"Times for snapshots (example: @SnapshotTimes 10 20 30 s;)",
				false));


			// ------------------------------------------------------------------------------------------------
			// Sinusoidal signals section
			// ------------------------------------------------------------------------------------------------

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FrequencyFuel",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency of sinusoidal signal (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SemiAmplitudeFuel",
				OpenSMOKE::SINGLE_DOUBLE,
				"Semiamplitude (relative) of sinusoidal signal (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DelayTimeFuel",
				OpenSMOKE::SINGLE_MEASURE,
				"Delay time of sinusoidal signal (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FrequencyOxidizer",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency of sinusoidal signal (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SemiAmplitudeOxidizer",
				OpenSMOKE::SINGLE_DOUBLE,
				"Semiamplitude (relative) of sinusoidal signal (default: 0)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DelayTimeOxidizer",
				OpenSMOKE::SINGLE_MEASURE,
				"Delay time of sinusoidal signal (default: 0)",
				false));

		}
	};
}

#endif /* OpenSMOKE_Grammar_DynamicBoundaries_H */