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

#ifndef OpenSMOKE_Grammar_CounterFlowFlame1D_H
#define OpenSMOKE_Grammar_CounterFlowFlame1D_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_CounterFlowFlame1D : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder containing the kinetic scheme (XML Version)",
				true,
				"@KineticsPreProcessor",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Backup",
				OpenSMOKE::SINGLE_PATH,
				"Name of backup file (XML Version)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the list of kinetic files to be interpreted",
				true,
				"@KineticsFolder",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Type of simulation: CounterFlowDiffusion | BurnerStabilizedStagnation",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Grid",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the mesh",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelStream",
				OpenSMOKE::VECTOR_STRING,
				"Name of the dictionary/dictionaries defining the fuel stream composition, temperature and pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerStream",
				OpenSMOKE::VECTOR_STRING,
				"Name of the dictionary/dictionaries defining the oxidizer stream composition, temperature and pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PeakMixture",
				OpenSMOKE::VECTOR_STRING,
				"Name of the dictionary/dictionaries defining the peak mixture composition, temperature and pressure",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Fuel stream velocity",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Oxidizer stream velocity",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing additional options for solving the sensitivity analysis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DaeParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the stiff DAE system",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NlsParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the NL system",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FalseTransientParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the pseudo-transient phase",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Output",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder where to write the output files",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseDaeSolver",
				OpenSMOKE::SINGLE_BOOL,
				"Use the Dae solver (instead of NLS) to solve the steady state problems (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseNlsSolver",
				OpenSMOKE::SINGLE_BOOL,
				"Use the NLS solver to solve the steady state problems, after the DAE solver (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PolimiSoot",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the rules for analyzing soot calculated using the Polimi mechanism",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@HMOM",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the rules for applying the Hybrid Method of Moments (HMOM)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VirtualChemistry",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the rules for applying the Virtual Chemistry (VC)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FixedTemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary describing the fixed temperature profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialTemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary describing the initial temperature profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Soret",
				OpenSMOKE::SINGLE_BOOL,
				"Add Soret effect (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SimplifiedTransportProperties",
				OpenSMOKE::SINGLE_BOOL,
				"Simplified transport properties (default: false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesBundling",
				OpenSMOKE::SINGLE_DOUBLE,
				"Estimation of mass diffusion coefficients through the species bundling according the the specified maximum relative error (example: @SpeciesBundling 0.1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialProfiles",
				OpenSMOKE::SINGLE_STRING,
				"Type of initial profiles for species and temperature: triangular | linear (default) | plateau",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PeakPosition",
				OpenSMOKE::SINGLE_MEASURE,
				"Position of peak of temperature (if not provided it will be assumed equal to the position of the stagnation plane)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MixingZoneWidth",
				OpenSMOKE::SINGLE_MEASURE,
				"Width of the mixing zone",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EigenValueStartingGuess",
				OpenSMOKE::SINGLE_MEASURE,
				"Starting guess value for the eigen value: default -100 kg/m3/s2",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RadialGradientFuelSide",
				OpenSMOKE::SINGLE_MEASURE,
				"Radial gradient on the fuel side: default 0 1/s",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RadialGradientOxidizerSide",
				OpenSMOKE::SINGLE_MEASURE,
				"Radial gradient on the oxidizer side: default 0 1/s",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PlanarSymmetry",
				OpenSMOKE::SINGLE_BOOL,
				"Planar symmetry: default false (i.e. cylyndrical symmetry)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GasReactionRateMultiplier",
				OpenSMOKE::SINGLE_DOUBLE,
				"Optional multiplying factor applied to all the reaction rates (default is 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LewisNumbers",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary containing the list of Lewis numbers of species",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyPostProcessing",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary specifying the details for carrying out the post-processing analyses (on the fly)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Radiation",
				OpenSMOKE::SINGLE_BOOL,
				"Radiative heat transfer (default: none)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EnvironmentTemperature",
				OpenSMOKE::SINGLE_MEASURE,
				"Environment temperature (default: 298.15 K)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DepositionWall",
				OpenSMOKE::SINGLE_BOOL,
				"Deposition wall for burner stabilized stagnation flames",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DynamicBoundaries",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the definition of dynamic boundaries",
				false,
				"none",
				"@Backup",
				"none"));
		}
	};
}

#endif /* OpenSMOKE_Grammar_CounterFlowFlame1D_H */