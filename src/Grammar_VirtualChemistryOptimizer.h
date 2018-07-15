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
|   Copyright(C) 2014, 2013  Alberto Cuoci                                |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_VirtualChemistryOptimizer : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
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

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the list of kinetic files to be interpreted",
				true,
				"@KineticsFolder",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VirtualChemistry",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the rules for applying the Virtual Chemistry (VC)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfExperiments_PremixedFlames",
				OpenSMOKE::VECTOR_STRING,
				"List of input files for experiments: premixed 1D flames",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfExperiments_CounterFlowFlames",
				OpenSMOKE::VECTOR_STRING,
				"List of input files for experiments: counterflow 1D flames",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfExperiments_PlugFlowReactors",
				OpenSMOKE::VECTOR_STRING,
				"List of input files for experiments: plugflow reactors",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfParameters",
				OpenSMOKE::VECTOR_INT,
				"List of optimization parameters (from 0 to 22)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfFirstGuess",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of first guess values",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfRelativeMinima",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of relative minima",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfRelativeMaxima",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of relative maxima",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfAbsoluteMinima",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of absolute minima",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfAbsoluteMaxima",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of absolute maxima",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RelativeErrors",
				OpenSMOKE::SINGLE_BOOL,
				"Objective function of relative errors (default: true)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxIterations",
				OpenSMOKE::SINGLE_INT,
				"Max number of evaluations (default: 100000)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Algorithm",
				OpenSMOKE::SINGLE_STRING,
				"Algorithm OpenSMOKEpp-Simplex | DIRECT | CRS | MLSL | STOGO | ISRES | ESCH | COBYLA | BOBYQA |NEWUOA | PRAXIS | NELDERMEAD | SBPLX (default: OpenSMOKEpp-Simplex)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Variant",
				OpenSMOKE::SINGLE_STRING,
				"Algorithm variant (default: none)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationTarget",
				OpenSMOKE::SINGLE_STRING,
				"Optimization target: main | CO | NO (default: main)",
				true));

		}
	};
}

