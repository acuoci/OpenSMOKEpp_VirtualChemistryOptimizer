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

// Include OpenMP Header file
#if defined(_OPENMP)
#include <omp.h>
#endif

// Math library
#include <math.h>

// NLopt library
#if NLOPT_LIB == 1
#include <nlopt.h>
#endif

// Optim library
#if OPTIM_LIB == 1
#include "optim/optim.hpp"
#endif

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"

// PolimiSoot Analyzer
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Virtual Chemistry
#include "utilities/virtualchemistry/VirtualChemistry.h"

// Optimizer
#include "Grammar_VirtualChemistryOptimizer.h"
#include "math/native-minimization-solvers/MinimizationSimplex.h"

// Reactors/Flames
#include "pfr/PlugFlowReactorExperiment.h"
#include "premixed-1dflame/Premixed1DFlameExperiment.h"
#include "counterflow-1dflame/CounterFlow1DFlameExperiment.h"

const double UNFEASIBLE_BIG_NUMBER = 1.e16;

void FromMinimizationParametersToRealParameters(const Eigen::VectorXd& b, Eigen::VectorXd& parameters, const std::vector<bool>& flag);
void FromRealParametersToMinimizationParameters(const Eigen::VectorXd& parameters, Eigen::VectorXd& b, Eigen::VectorXd& bMin, Eigen::VectorXd& bMax, const std::vector<bool>& flag);
void WriteTables(const Eigen::VectorXd& parameters, const std::vector<bool>& flag);

double OptFunction(const Eigen::VectorXd& b);

#if NLOPT_LIB == 1
double NLOptFunction(unsigned n, const double *x, double *grad, void *my_func_data);
#endif

#if OPTIM_LIB == 1
double OPTIMFunction(const arma::vec& x, arma::vec* grad_out, void* opt_data);
#endif

OpenSMOKE::VirtualChemistry* virtual_chemistry;
OpenSMOKE::PlugFlowReactorExperiment* plugs;
OpenSMOKE::Premixed1DFlameExperiment* premixed_flames;
OpenSMOKE::CounterFlow1DFlameExperiment* counterflow_flames;

const double Rgas = 1987.;
std::vector<bool> flag;

std::vector<std::string> list_of_experiments_premixed_flames;
std::vector<std::string> list_of_experiments_plugflow_reactors;
std::vector<std::string> list_of_experiments_counterflow_flames;

bool obj_function_relative_errors;

int number_parameters;
bool central_difference;
std::string optimization_target;

double fobj_rel_best;
double fobj_abs_best;
std::ofstream fMonitoring;
std::ofstream fMonitoringBest;

std::vector<double> list_of_relative_minima;
std::vector<double> list_of_absolute_minima;
std::vector<double> list_of_relative_maxima;
std::vector<double> list_of_absolute_maxima;

unsigned int numberOfGradientEvaluations;
unsigned int numberOfFunctionEvaluations;

int main(int argc, char** argv)
{

	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_VirtualChemistryOptimizer", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "Optimizer";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_VirtualChemistryOptimizer");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"PlugFlowReactor\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				main_dictionary_name_ = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Defines the grammar rules
	OpenSMOKE::Grammar_VirtualChemistryOptimizer grammar_vco;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_vco);


	// Kinetic scheme
	boost::filesystem::path path_kinetics_output;
	if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", path_kinetics_output);
		OpenSMOKE::CheckKineticsFolder(path_kinetics_output);
	}
	else
	{
		std::string name_of_rapid_kinetics_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@KineticsPreProcessor") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@KineticsPreProcessor", name_of_rapid_kinetics_subdictionary);

		OpenSMOKE::Grammar_RapidKineticMechanism grammar_rapid_kinetics;
		dictionaries(name_of_rapid_kinetics_subdictionary).SetGrammar(grammar_rapid_kinetics);


		boost::filesystem::path path_input_thermodynamics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Thermodynamics") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Thermodynamics", path_input_thermodynamics);

		boost::filesystem::path path_input_kinetics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Kinetics") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Kinetics", path_input_kinetics);

		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_kinetics_output);

		OpenSMOKE::RapidKineticMechanismWithoutTransport(path_kinetics_output, path_input_thermodynamics.c_str(), path_input_kinetics.c_str());
	}


	// Read thermodynamics and kinetics maps
	OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML;
	OpenSMOKE::KineticsMap_CHEMKIN*				kineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN*	transportMapXML;
	{
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_output / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc);
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);

		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
	}

	// Optimization target
	optimization_target = "main";
	if (dictionaries(main_dictionary_name_).CheckOption("@OptimizationTarget") == true)
		dictionaries(main_dictionary_name_).ReadString("@OptimizationTarget", optimization_target);

	// Virtual Chemistry
	if (dictionaries(main_dictionary_name_).CheckOption("@VirtualChemistry") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@VirtualChemistry", name_of_subdictionary);
		virtual_chemistry = new OpenSMOKE::VirtualChemistry(*thermodynamicsMapXML, dictionaries(name_of_subdictionary));
		virtual_chemistry->SetOnTheFlyOptimization(optimization_target);
	}

	// List of experiments
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfExperiments_PremixedFlames") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfExperiments_PremixedFlames", list_of_experiments_premixed_flames);
	
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfExperiments_CounterFlowFlames") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfExperiments_CounterFlowFlames", list_of_experiments_counterflow_flames);
	
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfExperiments_PlugFlowReactors") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfExperiments_PlugFlowReactors", list_of_experiments_plugflow_reactors);

	// List of optimization parameters
	std::vector<int> list_of_parameters;
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfParameters") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfParameters", list_of_parameters);

	// List of optimization parameters
	std::vector<double> list_of_first_guess;
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfFirstGuess") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfFirstGuess", list_of_first_guess);

	// List of relative minima
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfRelativeMinima") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfRelativeMinima", list_of_relative_minima);

	// List of absolute minima
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfAbsoluteMinima") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfAbsoluteMinima", list_of_absolute_minima);

	// List of relative maxima
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfRelativeMaxima") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfRelativeMaxima", list_of_relative_maxima);

	// List of absolute maxima
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfAbsoluteMaxima") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfAbsoluteMaxima", list_of_absolute_maxima);

	// Objective function on relative errors
	obj_function_relative_errors = true;
	if (dictionaries(main_dictionary_name_).CheckOption("@RelativeErrors") == true)
		dictionaries(main_dictionary_name_).ReadBool("@RelativeErrors", obj_function_relative_errors);

	// Algorithm
	std::string algorithm = "OpenSMOKEpp-Simplex";
	if (dictionaries(main_dictionary_name_).CheckOption("@Algorithm") == true)
		dictionaries(main_dictionary_name_).ReadString("@Algorithm", algorithm);

	// Variant
	std::string variant = "";
	if (dictionaries(main_dictionary_name_).CheckOption("@Variant") == true)
		dictionaries(main_dictionary_name_).ReadString("@Variant", variant);

	// Central/Forward gradient
	central_difference = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@CentralGradient") == true)
		dictionaries(main_dictionary_name_).ReadBool("@CentralGradient", central_difference);

	// Max number iterations
	int max_eval = 1000000;
	if (dictionaries(main_dictionary_name_).CheckOption("@MaxIterations") == true)
		dictionaries(main_dictionary_name_).ReadInt("@MaxIterations", max_eval);

	// Experiments

	// Create premixed 1Dflames
	premixed_flames = new OpenSMOKE::Premixed1DFlameExperiment[list_of_experiments_premixed_flames.size()];
	for (unsigned int k = 0; k < list_of_experiments_premixed_flames.size(); k++)
	{
		premixed_flames[k].Setup(list_of_experiments_premixed_flames[k], thermodynamicsMapXML, kineticsMapXML, transportMapXML, virtual_chemistry);
		premixed_flames[k].Solve(true);
	}

	// Create counterflow 1Dflames
	counterflow_flames = new OpenSMOKE::CounterFlow1DFlameExperiment[list_of_experiments_counterflow_flames.size()];
	for (unsigned int k = 0; k < list_of_experiments_counterflow_flames.size(); k++)
	{
		counterflow_flames[k].Setup(list_of_experiments_counterflow_flames[k], thermodynamicsMapXML, kineticsMapXML, transportMapXML, virtual_chemistry);
		counterflow_flames[k].Solve(true);
	}

	// Creat plugflow reactors
	plugs = new OpenSMOKE::PlugFlowReactorExperiment[list_of_experiments_plugflow_reactors.size()];
	for (unsigned int k = 0; k < list_of_experiments_plugflow_reactors.size(); k++)
	{
		plugs[k].Setup(list_of_experiments_plugflow_reactors[k], thermodynamicsMapXML, kineticsMapXML, virtual_chemistry);
		plugs[k].Solve(true);
	}

	{
		// Create list of parameters
		if (optimization_target == "main") flag.resize(4);
		else if (optimization_target == "CO") flag.resize(0);
		else if (optimization_target == "NO") flag.resize(23);

		std::fill(flag.begin(), flag.end(), false);
		for (unsigned int k = 0; k < list_of_parameters.size(); k++)
			flag[list_of_parameters[k]] = true;
		number_parameters = std::count(flag.begin(), flag.end(), true);

		// Check
		if (list_of_first_guess.size() != 0 && list_of_first_guess.size() != number_parameters)
			OpenSMOKE::FatalErrorMessage("Error in specifying list of first guess values");
		if (list_of_absolute_maxima.size() != 0 && list_of_absolute_maxima.size() != number_parameters)
			OpenSMOKE::FatalErrorMessage("Error in specifying list of absolute maxima values");
		if (list_of_relative_maxima.size() != 0 && list_of_relative_maxima.size() != number_parameters)
			OpenSMOKE::FatalErrorMessage("Error in specifying list of relative maxima values");
		if (list_of_absolute_minima.size() != 0 && list_of_absolute_minima.size() != number_parameters)
			OpenSMOKE::FatalErrorMessage("Error in specifying list of absolute minima values");
		if (list_of_relative_minima.size() != 0 && list_of_relative_minima.size() != number_parameters)
			OpenSMOKE::FatalErrorMessage("Error in specifying list of relative minima values");

		// Get first guess values
		Eigen::VectorXd parameters0(number_parameters);
		if (list_of_first_guess.size() == 0)
		{
			virtual_chemistry->GetParameters(parameters0, flag);
		}
		else
		{
			for (unsigned int j = 0; j < list_of_first_guess.size(); j++)
				parameters0(j) = list_of_first_guess[j];
		}

		WriteTables(parameters0, flag);

		// Get minimization parameters
		Eigen::VectorXd b0(number_parameters), bMin(number_parameters), bMax(number_parameters);
		FromRealParametersToMinimizationParameters(parameters0, b0, bMin, bMax, flag);

		// Print data on screen
		std::cout << "Parameters: " << std::endl;
		for (unsigned int i = 0; i < number_parameters; i++)
			std::cout << i << " " << parameters0(i) << " " << b0(i) << " " << bMin(i) << " " << bMax(i) << std::endl;

		// Check
		for (unsigned int i = 0; i < number_parameters; i++)
		{
			if (bMin(i) >= bMax(i))
				OpenSMOKE::FatalErrorMessage("Error in min/max constraints: min > max");
			if (b0(i) <= bMin(i))
				OpenSMOKE::FatalErrorMessage("Error in min/max constraints: first guess value <= min");
			if (b0(i) <= bMin(i))
				OpenSMOKE::FatalErrorMessage("Error in min/max constraints: first guess value >= max");
		}

		// Calculates objective function in the starting point
		const double f0 = OptFunction(b0);

		// Start optimization procedure
		fobj_rel_best = 1.e20;
		fobj_abs_best = 1.e20;
		fMonitoring.open("log", std::ios::out);
		fMonitoring.setf(std::ios::scientific);
		fMonitoringBest.open("log.best", std::ios::out);
		fMonitoringBest.setf(std::ios::scientific);

		if(algorithm == "OpenSMOKEpp-Simplex")
		{
			// Initialize the optimization
			OpenSMOKE::MinimizationSimplex mr;
			mr(b0, OptFunction, bMin, bMax);
			
			// Solve the optimization
			mr();

			// Convert optimization parameters in real parameters
			Eigen::VectorXd parametersOpt(number_parameters);
			Eigen::VectorXd bOpt(number_parameters);
			mr.GetSolution(bOpt);
			FromMinimizationParametersToRealParameters(bOpt, parametersOpt, flag);

			// Write on the screen
			const double fOpt = OptFunction(bOpt);
			std::cout << "Objective function: " << f0 << " -> " << fOpt << std::endl;
			std::cout << "Optimal Parameters: " << std::endl;
			for (unsigned int i = 0; i < number_parameters; i++)
				std::cout << i << " " << parameters0(i) << " -> " << parametersOpt(i) << std::endl;

			// Write on file
			WriteTables(parametersOpt, flag);
		}

		// 1) DIRECT (DIviding RECTangles)

		// NLOPT_GN_DIRECT: DIviding RECTangles algorithm for global optimization:
		//                  D.R.Jones, C.D.Perttunen, and B.E.Stuckmann
		//                  "Lipschitzian optimization without the lipschitz constant"
		//                  J.Optimization Theory and Applications, vol. 79, p. 157 (1993)
		
		// NLOPT_GN_DIRECT_L: is the "locally biased" variant
		//                    J.M.Gablonsky and C.T.Kelley
		//                    "A locally-biased form of the DIRECT algorithm"
		//                    J.Global Optimization, vol. 21(1), p. 27-37 (2001)

		// NLOPT_GN_DIRECT_L_RAND: randomized variant of DIRECT-L

		// NLOPT_GNL_DIRECT_NOSCAL: unscaled variant
		// NLOPT_GNL_DIRECT_L_NOSCAL: unscaled variant
		// NLOPT_GNL_DIRECT_L_RAND_NOSCAL: unscaled variant

		// NLOPT_GN_ORIG_DIRECT: original implementation by J.M.Gablonsky

		// NLOPT_GN_ORIG_DIRECT_L: original implementation by J.M.Gablonsky


		// 2) CONTROLLED RANDOM SEARCH ALGORITHMS 

		// NLOPT_GN_CRS2_LM: Controlled Random Search (CRS) algorithm (CRS2 variant) with the "local mutation" modification
		//                   P.Kaelo and M.M.Ali
		//                   "Some variants of the controlled random search algorithm for global optimization"
		//                   J.Optim.Theory Appl. 130 (2), 253 - 264 (2006).
		//                   W.L.Price
		//                   "A controlled random search procedure for global optimization," in Towards Global Optimization 2, p. 71-84
		//                   Edited by L.C.W.Dixon and G.P.Szego(North - Holland Press, Amsterdam, 1978).
		//                   W.L.Price
		//                   "Global optimization by controlled random search"
		//                   J.Optim.Theory Appl. 40 (3), p. 333 - 348 (1983).


		// 3) MLSL (Multi-Level Single-Linkage)
		
		// NLOPT_G_MLSL: original MLSL algorithm (using pseudo-random numbers via the Mersenne twister algorithm)
		//                A.H.G. Rinnooy Kan and G. T. Timmer
		//                "Stochastic global optimization methods"
		//                Mathematical Programming, vol. 39, p. 27-78 (1987)

		// NLOPT_G_MLSL_LDS: LDS-based MLSL algorithm (Sobol's low-discrepancy sequence (LDS) instead of pseudorandom numbers)
		//                   S. Kucherenko and Y. Sytsko
		//                   "Application of deterministic low-discrepancy sequences in global optimization"
		//                   Computational Optimization and Applications, vol. 30, p. 297-318 (2005).
		
		// In both cases the local optimizer must be set: nlopt_opt_set_local_optimizer 

		#if NLOPT_LIB == 1

		else if (	algorithm == "DIRECT" || algorithm == "CRS" ||
					algorithm == "MLSL" || algorithm == "STOGO" ||
					algorithm == "ISRES" || algorithm == "ESCH" ||
					algorithm == "COBYLA" || algorithm == "BOBYQA" ||
					algorithm == "NEWUOA" || algorithm == "PRAXIS" ||
					algorithm == "NELDERMEAD" || algorithm == "SBPLX" ||
					algorithm == "SLSQP" || algorithm == "LBFGS" || 
					algorithm == "TNEWTON_PRECOND" || algorithm == "SLM_VAR" )
		{
			nlopt_opt opt;
			
			if (algorithm == "DIRECT")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_DIRECT, number_parameters);
				else if (variant == "L")
					opt = nlopt_create(NLOPT_GN_DIRECT_L, number_parameters);
				else if (variant == "L_RAND")
					opt = nlopt_create(NLOPT_GN_DIRECT_L_RAND, number_parameters);
				else if (variant == "NOSCAL")
					opt = nlopt_create(NLOPT_GN_DIRECT_NOSCAL, number_parameters);
				else if (variant == "L_NOSCAL")
					opt = nlopt_create(NLOPT_GN_DIRECT_L_NOSCAL, number_parameters);
				else if (variant == "L_RAND_NOSCAL")
					opt = nlopt_create(NLOPT_GN_DIRECT_L_RAND_NOSCAL, number_parameters);
				else if (variant == "ORIG")
					opt = nlopt_create(NLOPT_GN_ORIG_DIRECT, number_parameters);
				else if (variant == "L_ORIG")
					opt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for DIRECT algorithms: \
						NONE | L | L_RAND | NOSCAL | L_NOSCAL | L_RAND_NOSCAL | ORIG | L_ORIG");
			}
			else if (algorithm == "CRS")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_CRS2_LM, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for CRS algorithms: NONE ");
			}
			else if (algorithm == "MLSL")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_G_MLSL, number_parameters);
				else if (variant == "LDS")
					opt = nlopt_create(NLOPT_G_MLSL_LDS, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for MLSL algorithms: NONE | LDS");

				nlopt_opt local_opt;
				nlopt_create(NLOPT_LN_COBYLA, number_parameters);
				nlopt_set_min_objective(local_opt, NLOptFunction, NULL);
				nlopt_set_local_optimizer(opt, local_opt);
			}
			else if (algorithm == "STOGO")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GD_STOGO, number_parameters);
				else if (variant == "RAND")
					opt = nlopt_create(NLOPT_GD_STOGO_RAND, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for STOGO algorithms: NONE | RAND");
			}
			else if (algorithm == "ISRES")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_ISRES, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for ISRES algorithms: NONE");
			}
			else if (algorithm == "ESCH")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_ESCH, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for ESCH algorithms: NONE");
			}
			else if (algorithm == "COBYLA")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_COBYLA, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for COBYLA algorithms: NONE");
			}
			else if (algorithm == "BOBYQA")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_BOBYQA, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for BOBYQA algorithms: NONE");
			}
			else if (algorithm == "NEWUOA")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_NEWUOA_BOUND, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for NEWUOA algorithms: NONE");
			}
			else if (algorithm == "PRAXIS")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_PRAXIS, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for PRAXIS algorithms: NONE");
			}
			else if (algorithm == "NELDERMEAD")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_NELDERMEAD, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for NELDERMEAD algorithms: NONE");
			}
			else if (algorithm == "SBPLX")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_SBPLX, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for SBPLX algorithms: NONE");
			}
			else if (algorithm == "SLSQP")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LD_SLSQP, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for SLSQP algorithms: NONE");
			}
			else if (algorithm == "LBFGS")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LD_LBFGS, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for LBFGS algorithms: NONE");
			}
			else if (algorithm == "TNEWTON_PRECOND")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, number_parameters);
				else if (variant == "NO_RESTART")
					opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND, number_parameters);
				else if (variant == "NO_PRECOND")
					opt = nlopt_create(NLOPT_LD_TNEWTON_RESTART, number_parameters);
				else if (variant == "NO_RESTART_NO_PRECOND")
					opt = nlopt_create(NLOPT_LD_TNEWTON, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for TNEWTON_PRECOND algorithms: NONE | NO_RESTART | NO_PRECOND | NO_RESTART_NO_PRECOND");
			}
			else if (algorithm == "SLM_VAR")
			{
				if (variant == "VAR2")
					opt = nlopt_create(NLOPT_LD_VAR2, number_parameters);
				else if (variant == "VAR1")
					opt = nlopt_create(NLOPT_LD_VAR1, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for TNEWTON_PRECOND algorithms: VAR2 | VAR1");
			}
			

			double* lb = new double[number_parameters];// lower bounds
			double* ub = new double[number_parameters];// upper bounds
			double* x = new double[number_parameters];// first guess
			for (unsigned int i = 0; i < number_parameters; i++)
			{
				lb[i] = bMin(i);
				ub[i] = bMax(i);
				x[i] = b0(i);
			}

			nlopt_set_lower_bounds(opt, lb);
			nlopt_set_upper_bounds(opt, ub);
			nlopt_set_min_objective(opt, NLOptFunction, NULL);
			nlopt_set_maxeval(opt, max_eval);

			double fOpt;
			if (nlopt_optimize(opt, x, &fOpt) < 0)
			{
				std::cout << "NLopt failed" << std::endl;
				getchar();
				exit(-1);
			}
			else
			{
				Eigen::VectorXd bOpt(number_parameters);
				for (unsigned int i = 0; i < number_parameters; i++)
					bOpt(i) = x[i];

				const double fOpt = OptFunction(bOpt);

				Eigen::VectorXd parametersOpt(number_parameters);
				FromMinimizationParametersToRealParameters(bOpt, parametersOpt, flag);

				std::cout << "Objective function: " << f0 << " -> " << fOpt << std::endl;
				std::cout << "Optimal Parameters: " << std::endl;
				for (unsigned int i = 0; i < number_parameters; i++)
					std::cout << i << " " << parameters0(i) << " " << parametersOpt(i) << std::endl;

				WriteTables(parametersOpt, flag);
			}
		}
		#endif

		#if OPTIM_LIB == 1
		else if ( algorithm == "OPTIM-DE" || algorithm == "OPTIM-PSO" || algorithm == "OPTIM-NM")
		{
			arma::vec lb(number_parameters);
			arma::vec ub(number_parameters);
			arma::vec x(number_parameters);
			for (unsigned int i = 0; i < number_parameters; i++)
			{
				lb[i] = bMin(i);
				ub[i] = bMax(i);
				x[i] = b0(i);
			}

			optim::algo_settings_t settings;
			settings.vals_bound = true;
			settings.lower_bounds = lb;
			settings.upper_bounds = ub;

			std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

			bool success = false;
			if (algorithm == "OPTIM-DE")
				success = optim::de(x,OPTIMFunction,nullptr, settings);
			else if (algorithm == "OPTIM-PSO")
				success = optim::pso(x,OPTIMFunction,nullptr, settings);
			else if (algorithm == "OPTIM-NM")
				success = optim::nm(x,OPTIMFunction,nullptr, settings);			

			std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;

			if (success == true) 
			{
				Eigen::VectorXd bOpt(number_parameters);
				for (unsigned int i = 0; i < number_parameters; i++)
					bOpt(i) = x[i];

				const double fOpt = OptFunction(bOpt);

				Eigen::VectorXd parametersOpt(number_parameters);
				FromMinimizationParametersToRealParameters(bOpt, parametersOpt, flag);

				std::cout << "Objective function: " << f0 << " -> " << fOpt << std::endl;
				std::cout << "Optimal Parameters: " << std::endl;
				for (unsigned int i = 0; i < number_parameters; i++)
					std::cout << i << " " << parameters0(i) << " " << parametersOpt(i) << std::endl;

				WriteTables(parametersOpt, flag);
			} 
			else 
			{
				std::cout << "OPTIM failed" << std::endl;
				getchar();
				exit(-1);
			}
		}
		#endif
		else
		{
			OpenSMOKE::FatalErrorMessage("Error @Algorithm option:	OpenSMOKEpp-Simplex | DIRECT | CRS | MLSL | STOGO | ISRES | ESCH | \
																		COBYLA | BOBYQA |NEWUOA | PRAXIS | NELDERMEAD | SBPLX | \
																		SLSQP | LBFGS | TNEWTON_PRECOND | SLM_VAR");
		}

		fMonitoring.close();
		fMonitoringBest.close();
	}

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_VirtualChemistryOptimizer", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	std::cout << "Optimization completed! Press enter to exit..." << std::endl;
	getchar();

	return OPENSMOKE_SUCCESSFULL_EXIT;
}


double ReturnObjFunction(const Eigen::VectorXd parameters)
{
	double fobj_abs = 0.;
	double fobj_rel = 0.;

	// Plug flow reactors
	for (unsigned int k = 0; k < list_of_experiments_plugflow_reactors.size(); k++)
	{
		plugs[k].Solve();
		fobj_abs += plugs[k].norm2_abs_error();
		fobj_rel += plugs[k].norm2_rel_error();
	}

	// Premixed 1D flame
	for (unsigned int k = 0; k < list_of_experiments_premixed_flames.size(); k++)
	{
		premixed_flames[k].Solve();
		fobj_abs += premixed_flames[k].norm2_abs_error();
		fobj_rel += premixed_flames[k].norm2_rel_error();
	}

	// Counterflow flames
	for (unsigned int k = 0; k < list_of_experiments_counterflow_flames.size(); k++)
	{
		counterflow_flames[k].Solve();
		fobj_abs += counterflow_flames[k].norm2_abs_error();
		fobj_rel += counterflow_flames[k].norm2_rel_error();
	}

	std::cout << "f(abs) = " << fobj_abs << "   -   f(rel) = " << fobj_rel << std::endl;

	if (obj_function_relative_errors == true)
	{
		if (fobj_rel < fobj_rel_best)
		{
			fobj_rel_best = fobj_rel;

			fMonitoringBest << std::left << std::setw(16) << fobj_rel;
			for (unsigned int i = 0; i<number_parameters; i++)
				fMonitoringBest << std::left << std::setw(16) << parameters(i);
			fMonitoringBest << std::endl;
		}

		fMonitoring << std::left << std::setw(16) << fobj_rel;
		for (unsigned int i = 0; i<number_parameters; i++)
			fMonitoring << std::left << std::setw(16) << parameters(i);
		fMonitoring << std::endl;

		return fobj_rel;
	}
	else
	{
		if (fobj_abs < fobj_abs_best)
		{
			fobj_abs_best = fobj_abs;

			fMonitoringBest << std::left << std::setw(16) << fobj_abs;
			for (unsigned int i = 0; i<number_parameters; i++)
				fMonitoringBest << std::left << std::setw(16) << parameters(i);
			fMonitoringBest << std::endl;
		}

		fMonitoring<< std::left << std::setw(16) << fobj_abs;
		for (unsigned int i = 0; i<number_parameters; i++)
			fMonitoring << std::left << std::setw(16) << parameters(i);
		fMonitoring << std::endl;

		return fobj_abs;
	}
}

double OptFunction(const Eigen::VectorXd &b)
{
	Eigen::VectorXd parameters(b.size());
	FromMinimizationParametersToRealParameters(b, parameters, flag);
	virtual_chemistry->SetParameters(parameters, flag);
	
	return ReturnObjFunction(parameters);
}

#if NLOPT_LIB == 1
double NLOptFunction(unsigned n, const double *x, double *grad, void *my_func_data)
{
	Eigen::VectorXd b(number_parameters);
	for (unsigned int i = 0; i < number_parameters; i++)
		b(i) = x[i];

	const double f = OptFunction(b);

	if (grad)
	{
		std::cout << "Gradient evaluation..." << std::endl;

		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double ETA3 = std::pow(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE, 1. / 3.);
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_DOUBLE);

		// Dimensions of parameters vector
		Eigen::VectorXd b_dimensions(number_parameters);
		b_dimensions.setConstant(1.);

		// Choosing between forward and centrale difference
		double eta = ETA2;
		if (central_difference == true)
			eta = ETA3;

		// Estimate the increment for forward approximation
		Eigen::VectorXd deltab(number_parameters);
		for (unsigned int i = 0; i < number_parameters; i++)
		{
			if (b(i) < 0.)
				deltab(i) = -eta * std::max(std::fabs(b(i)), std::fabs(b_dimensions(i)));
			else
				deltab(i) = eta * std::max(std::fabs(b(i)), std::fabs(b_dimensions(i)));

			if (deltab(i) == 0.)
				deltab(i) = ZERO_DER;
		}

		// Forward gradient
		if (central_difference == false)
		{
			Eigen::VectorXd b_plus = b;
			for (unsigned int j = 0; j < number_parameters; j++)
			{
				b_plus(j) = b(j) + deltab(j);
				const double f_plus = OptFunction(b_plus);

				grad[j] = (f_plus - f) / deltab(j);

				b_plus(j) = b(j);
			}
		}

		// Central gradient
		if (central_difference == true)
		{
			Eigen::VectorXd b_star = b;
			for (unsigned int j = 0; j < number_parameters; j++)
			{
				b_star(j) = b(j) + deltab(j);
				const double f_plus = OptFunction(b_star);
				b_star(j) = b(j) - deltab(j);
				const double f_minus = OptFunction(b_star);

				grad[j] = (f_plus - f_minus) / (2.*deltab(j));

				b_star(j) = b(j);
			}
		}

		numberOfGradientEvaluations++;
	}

	numberOfFunctionEvaluations++;

	return f;
}
#endif

#if OPTIM_LIB == 1
double OPTIMFunction(const arma::vec& x, arma::vec* grad_out, void* opt_data)
{
	Eigen::VectorXd b(number_parameters);
	for (unsigned int i = 0; i < number_parameters; i++)
		b(i) = x[i];

	const double f = OptFunction(b);

	return f;
}
#endif

void FromMinimizationParametersToRealParameters(const Eigen::VectorXd& b, Eigen::VectorXd& parameters, const std::vector<bool>& flag)
{
	unsigned int k = 0;

	if (optimization_target == "main")
	{
		for (unsigned int j = 0; j < 4; j++)	// alpha1-4 
			if (flag[j] == true)
			{
				parameters(k) = b(k);
				k++;
			}
	}
	else if (optimization_target == "CO")
	{
		// TODO
	}
	if (optimization_target == "NO")
	{
		for (unsigned int j = 0; j < 5; j++)	// A3_NO_, A4_NO_, A5_NO_, A6_NO_, A7_NO_ 
			if (flag[j] == true)
			{
				parameters(k) = std::exp(b(k));
				k++;
			}

		for (unsigned int j = 5; j < 10; j++)	// E3_NO_, E4_NO_, E5_NO_, E6_NO_, E7_NO_ 
			if (flag[j] == true)
			{
				parameters(k) = b(k) * Rgas;
				k++;
			}

		if (flag[10] == true)				// Kc5_NO_
		{
			parameters(k) = std::exp(b(k));
			k++;
		}

		for (unsigned int j = 11; j < 14; j++)	// alpha_NO, beta_NO_, gamma_NO_ 
			if (flag[j] == true)
			{
				parameters(k) = b(k);
				k++;
			}

		for (unsigned int j = 14; j < 23; j++)	// reaction orders
			if (flag[j] == true)
			{
				parameters(k) = b(k);
				k++;
			}
	}
}

void FromRealParametersToMinimizationParameters(const Eigen::VectorXd& parameters, Eigen::VectorXd& b, Eigen::VectorXd& bMin, Eigen::VectorXd& bMax, const std::vector<bool>& flag)
{
	unsigned int k = 0;

	if (optimization_target == "main")
	{
		for (unsigned int j = 0; j < 4; j++)	// alpha1-4
			if (flag[j] == true)
			{
				b(k) = parameters(k);

				if (list_of_relative_minima.size() != 0)		bMin(k) = parameters(k) * list_of_relative_minima[k];
				else if (list_of_absolute_minima.size() != 0)	bMin(k) = list_of_absolute_minima[k];
				else                                            bMin(k) = 0.;

				if (list_of_relative_maxima.size() != 0)		bMax(k) = parameters(k)* list_of_relative_maxima[k];
				else  if (list_of_absolute_maxima.size() != 0)	bMax(k) = list_of_absolute_maxima[k];
				else                                            bMax(k) = 1.;

				k++;
			}
	}
	else if (optimization_target == "CO")
	{
		// TODO
	}
	if (optimization_target == "NO")
	{
		for (unsigned int j = 0; j < 5; j++)	// A3_NO_, A4_NO_, A5_NO_, A6_NO_, A7_NO_ 
			if (flag[j] == true)
			{
				b(k) = std::log(parameters(k));
				
				if (list_of_relative_minima.size() != 0)		bMin(k) = std::log(parameters(k) * list_of_relative_minima[k]);
				else if (list_of_absolute_minima.size() != 0)	bMin(k) = std::log(list_of_absolute_minima[k]);
				else                                            bMin(k) = std::log(parameters(k) / 1000.);

				if (list_of_relative_maxima.size() != 0)		bMax(k) = std::log(parameters(k) * list_of_relative_maxima[k]);
				else  if (list_of_absolute_maxima.size() != 0)	bMax(k) = std::log(list_of_absolute_maxima[k]);
				else                                            bMax(k) = std::log(parameters(k) * 1000.);

				k++;
			}

		for (unsigned int j = 5; j < 10; j++)	// E3_NO_, E4_NO_, E5_NO_, E6_NO_, E7_NO_ 
			if (flag[j] == true)
			{
				b(k) = parameters(k) / Rgas;

				if (list_of_relative_minima.size() != 0)		bMin(k) = (parameters(k) / Rgas) * list_of_relative_minima[k];
				else if (list_of_absolute_minima.size() != 0)	bMin(k) = list_of_absolute_minima[k] / Rgas;
				else                                            bMin(k) = (parameters(k) / Rgas) / 1.3;

				if (list_of_relative_maxima.size() != 0)		bMax(k) = (parameters(k) / Rgas) * list_of_relative_maxima[k];
				else  if (list_of_absolute_maxima.size() != 0)	bMax(k) = list_of_absolute_maxima[k] / Rgas;
				else                                            bMax(k) = (parameters(k) / Rgas) * 1.3;
				
				k++;
			}


		if (flag[10] == true)				// Kc5_NO_
		{
			b(k) = std::log(parameters(k));

			if (list_of_relative_minima.size() != 0)		bMin(k) = std::log(parameters(k) * list_of_relative_minima[k]);
			else if (list_of_absolute_minima.size() != 0)	bMin(k) = std::log(list_of_absolute_minima[k]);
			else                                            bMin(k) = std::log(parameters(k) / 1000.);

			if (list_of_relative_maxima.size() != 0)		bMax(k) = std::log(parameters(k) * list_of_relative_maxima[k]);
			else  if (list_of_absolute_maxima.size() != 0)	bMax(k) = std::log(list_of_absolute_maxima[k]);
			else                                            bMax(k) = std::log(parameters(k) * 1000.);

			k++;
		}

		for (unsigned int j = 11; j < 14; j++)	// alpha_NO, beta_NO_, gamma_NO_ 
			if (flag[j] == true)
			{
				b(k) = parameters(k);

				if (list_of_relative_minima.size() != 0)		bMin(k) = parameters(k) * list_of_relative_minima[k];
				else if (list_of_absolute_minima.size() != 0)	bMin(k) = list_of_absolute_minima[k];
				else                                            bMin(k) = 0.;

				if (list_of_relative_maxima.size() != 0)		bMax(k) = parameters(k)* list_of_relative_maxima[k];
				else  if (list_of_absolute_maxima.size() != 0)	bMax(k) = list_of_absolute_maxima[k];
				else                                            bMax(k) = 1.;

				k++;
			}

		for (unsigned int j = 14; j < 23; j++)	// reaction orders
			if (flag[j] == true)
			{
				b(k) = parameters(k);

				if (list_of_relative_minima.size() != 0)		bMin(k) = parameters(k) * list_of_relative_minima[k];
				else if (list_of_absolute_minima.size() != 0)	bMin(k) = list_of_absolute_minima[k];
				else                                            bMin(k) = -1.;

				if (list_of_relative_maxima.size() != 0)		bMax(k) = parameters(k)* list_of_relative_maxima[k];
				else  if (list_of_absolute_maxima.size() != 0)	bMax(k) = list_of_absolute_maxima[k];
				else                                            bMax(k) = 3.;

				k++;
			}
	}
}

void WriteTables(const Eigen::VectorXd& parameters, const std::vector<bool>& flag)
{
	Eigen::VectorXd p;
	virtual_chemistry->GetOriginalParameters(p);

	if (optimization_target == "main")
	{
		double alpha1 = p(0);		// #0
		double alpha2 = p(1);		// #1
		double alpha3 = p(2);		// #2
		double alpha4 = p(3);		// #3

		unsigned int k = 0;

		if (flag[0] == true)	alpha1 = parameters(k++);		// #0
		if (flag[1] == true)	alpha2 = parameters(k++);		// #1
		if (flag[2] == true)	alpha3 = parameters(k++);		// #2
		if (flag[3] == true)	alpha4 = parameters(k++);		// #3

		{
			std::ofstream fTable("Table.main", std::ios::out);
			fTable.setf(std::ios::scientific);
			fTable << "Y_N2[-]        alpha1[-]         alpha2[-]        alpha3[-]       alpha4[-]" << std::endl;
			fTable << "4 3" << std::endl;
			fTable << std::left << std::setw(16) << "6.00000E-01";
			fTable << std::left << std::setw(16) << alpha1;
			fTable << std::left << std::setw(16) << alpha2;
			fTable << std::left << std::setw(16) << alpha3;
			fTable << std::left << std::setw(16) << alpha4;
			fTable << std::endl;
			fTable << std::left << std::setw(16) << "7.00000E-01";
			fTable << std::left << std::setw(16) << alpha1;
			fTable << std::left << std::setw(16) << alpha2;
			fTable << std::left << std::setw(16) << alpha3;
			fTable << std::left << std::setw(16) << alpha4;
			fTable << std::endl;
			fTable << std::left << std::setw(16) << "8.00000E-01";
			fTable << std::left << std::setw(16) << alpha1;
			fTable << std::left << std::setw(16) << alpha2;
			fTable << std::left << std::setw(16) << alpha3;
			fTable << std::left << std::setw(16) << alpha4;
			fTable << std::endl;
			fTable << "END";
			fTable << std::endl;
			fTable.close();
		}
	}
	else if (optimization_target == "CO")
	{
		// TODO
	}
	if (optimization_target == "NO")
	{
		const double A3_NO_orig = p(0);		// #0
		const double A4_NO_orig = p(1);		// #1
		const double A5_NO_orig = p(2);		// #2
		const double A6_NO_orig = p(3);		// #3
		const double A7_NO_orig = p(4);		// #4

		double A3_NO_ = A3_NO_orig;		// #0
		double A4_NO_ = A4_NO_orig;		// #1
		double A5_NO_ = A5_NO_orig;		// #2
		double A6_NO_ = A6_NO_orig;		// #3
		double A7_NO_ = A7_NO_orig;		// #4

		double E3_NO_ = p(5);			// #5
		double E4_NO_ = p(6);			// #6
		double E5_NO_ = p(7);			// #7
		double E6_NO_ = p(8);			// #8
		double E7_NO_ = p(9);			// #9
		double Kc5_NO_ = p(10);			// #10

		double alpha_NO_ = p(11);		// #11
		double beta_NO_ = p(12);		// #12
		double gamma_NO_ = p(13);		// #13

		double nuF_3_NO_ = p(14);		// #14
		double nuOX_3_NO_ = p(15);		// #15
		double nuW1_4_NO_ = p(16);		// #16
		double nuNO_5f_NO_ = p(17);		// #17
		double nuW2_5f_NO_ = p(18);		// #18
		double nuNO_5b_NO_ = p(19);		// #19
		double nuW2_5b_NO_ = p(20);		// #20
		double nuW3_6_NO_ = p(21);		// #21
		double nuW3_7_NO_ = p(22);		// #22

		/*
		const double A3_NO_orig = 1.538479670000e+18;		// #0
		const double A4_NO_orig = 1.505543580620e+14;		// #1
		const double A5_NO_orig = 9.59355596366e+25;		// #2
		const double A6_NO_orig = 4.289769943040e+15;		// #3
		const double A7_NO_orig = 3.682428621600e+22;		// #4

		double A3_NO_ = A3_NO_orig;				// #0
		double A4_NO_ = A4_NO_orig;				// #1
		double A5_NO_ = A5_NO_orig;				// #2
		double A6_NO_ = A6_NO_orig;				// #3
		double A7_NO_ = A7_NO_orig;				// #4

		double E3_NO_ = 35075.1745;				// #5
		double E4_NO_ = 57094.3505553;			// #6
		double E5_NO_ = 188550.271919;			// #7
		double E6_NO_ = 175095.668571;			// #8
		double E7_NO_ = 46937.7896765;			// #9
		double Kc5_NO_ = 7.992241E-03;			// #10

		double alpha_NO_ = 0.1;					// #11
		double beta_NO_ = 1.407162E-03;			// #12
		double gamma_NO_ = 0.;					// #13

		double nuF_3_NO_ = 1.70998297;			// #14
		double nuOX_3_NO_ = 0.86862947;			// #15
		double nuW1_4_NO_ = 1.54432793325;		// #16
		double nuNO_5f_NO_ = 0.00887992386821;	// #17
		double nuW2_5f_NO_ = 2.95870228167;		// #18
		double nuNO_5b_NO_ = 1.00887992387;		// #19
		double nuW2_5b_NO_ = 1.95870228167;		// #20
		double nuW3_6_NO_ = 0.400450422151;		// #21
		double nuW3_7_NO_ = 2.76521242592;		// #22
		*/


		unsigned int k = 0;

		if (flag[0] == true)	A3_NO_ = parameters(k++);		// #0
		if (flag[1] == true)	A4_NO_ = parameters(k++);		// #1
		if (flag[2] == true)	A5_NO_ = parameters(k++);		// #2
		if (flag[3] == true)	A6_NO_ = parameters(k++);		// #3
		if (flag[4] == true)	A7_NO_ = parameters(k++);		// #4

		if (flag[5] == true)	E3_NO_ = parameters(k++);		// #5
		if (flag[6] == true)	E4_NO_ = parameters(k++);		// #6
		if (flag[7] == true)	E5_NO_ = parameters(k++);		// #7
		if (flag[8] == true)	E6_NO_ = parameters(k++);		// #8
		if (flag[9] == true)	E7_NO_ = parameters(k++);		// #9
		if (flag[10] == true)	Kc5_NO_ = parameters(k++);		// #10

		if (flag[11] == true)	alpha_NO_ = parameters(k++);	// #11
		if (flag[12] == true)	beta_NO_ = parameters(k++);		// #12
		if (flag[13] == true)	gamma_NO_ = parameters(k++);	// #13

		if (flag[14] == true)	nuF_3_NO_ = parameters(k++);	// #14
		if (flag[15] == true)	nuOX_3_NO_ = parameters(k++);	// #15
		if (flag[16] == true)	nuW1_4_NO_ = parameters(k++);	// #16
		if (flag[17] == true)	nuNO_5f_NO_ = parameters(k++);	// #17
		if (flag[18] == true)	nuW2_5f_NO_ = parameters(k++);	// #18
		if (flag[19] == true)	nuNO_5b_NO_ = parameters(k++);	// #19
		if (flag[20] == true)	nuW2_5b_NO_ = parameters(k++);	// #20
		if (flag[21] == true)	nuW3_6_NO_ = parameters(k++);	// #21
		if (flag[22] == true)	nuW3_7_NO_ = parameters(k++);	// #22

		{
			std::ofstream fTable1("Table_NO.1", std::ios::out);
			fTable1.setf(std::ios::scientific);
			fTable1 << "Y_N2[-]        f3[-] 		   f5f[-] 		   f5b[-] 		   f4[-] 		   f6[-] 		   f7[-]" << std::endl;
			fTable1 << "6 3" << std::endl;
			fTable1 << std::left << std::setw(16) << "6.00000E-01";
			fTable1 << std::left << std::setw(16) << A3_NO_ / A3_NO_orig;
			fTable1 << std::left << std::setw(16) << A5_NO_ / A5_NO_orig;
			fTable1 << std::left << std::setw(16) << A5_NO_ / A5_NO_orig;
			fTable1 << std::left << std::setw(16) << A4_NO_ / A4_NO_orig;
			fTable1 << std::left << std::setw(16) << A6_NO_ / A6_NO_orig;
			fTable1 << std::left << std::setw(16) << A7_NO_ / A7_NO_orig;
			fTable1 << std::endl;
			fTable1 << std::left << std::setw(16) << "7.00000E-01";
			fTable1 << std::left << std::setw(16) << A3_NO_ / A3_NO_orig;
			fTable1 << std::left << std::setw(16) << A5_NO_ / A5_NO_orig;
			fTable1 << std::left << std::setw(16) << A5_NO_ / A5_NO_orig;
			fTable1 << std::left << std::setw(16) << A4_NO_ / A4_NO_orig;
			fTable1 << std::left << std::setw(16) << A6_NO_ / A6_NO_orig;
			fTable1 << std::left << std::setw(16) << A7_NO_ / A7_NO_orig;
			fTable1 << std::endl;
			fTable1 << std::left << std::setw(16) << "8.00000E-01";
			fTable1 << std::left << std::setw(16) << A3_NO_ / A3_NO_orig;
			fTable1 << std::left << std::setw(16) << A5_NO_ / A5_NO_orig;
			fTable1 << std::left << std::setw(16) << A5_NO_ / A5_NO_orig;
			fTable1 << std::left << std::setw(16) << A4_NO_ / A4_NO_orig;
			fTable1 << std::left << std::setw(16) << A6_NO_ / A6_NO_orig;
			fTable1 << std::left << std::setw(16) << A7_NO_ / A7_NO_orig;
			fTable1 << std::endl;
			fTable1 << "END";
			fTable1 << std::endl;
			fTable1.close();
		}

		{
			std::ofstream fTable2("Table_NO.2", std::ios::out);
			fTable2.setf(std::ios::scientific);
			fTable2 << "Y_N2[-]        K_c[-]          beta[-]         alpha[-]        gamma[-]" << std::endl;
			fTable2 << "4 3" << std::endl;
			fTable2 << std::left << std::setw(16) << "6.00000E-01";
			fTable2 << std::left << std::setw(16) << Kc5_NO_;
			fTable2 << std::left << std::setw(16) << beta_NO_;
			fTable2 << std::left << std::setw(16) << alpha_NO_;
			fTable2 << std::left << std::setw(16) << gamma_NO_;
			fTable2 << std::endl;
			fTable2 << std::left << std::setw(16) << "7.00000E-01";
			fTable2 << std::left << std::setw(16) << Kc5_NO_;
			fTable2 << std::left << std::setw(16) << beta_NO_;
			fTable2 << std::left << std::setw(16) << alpha_NO_;
			fTable2 << std::left << std::setw(16) << gamma_NO_;
			fTable2 << std::endl;
			fTable2 << std::left << std::setw(16) << "8.00000E-01";
			fTable2 << std::left << std::setw(16) << Kc5_NO_;
			fTable2 << std::left << std::setw(16) << beta_NO_;
			fTable2 << std::left << std::setw(16) << alpha_NO_;
			fTable2 << std::left << std::setw(16) << gamma_NO_;
			fTable2 << std::endl;
			fTable2 << "END";
			fTable2 << std::endl;
			fTable2.close();
		}

		{
			std::ofstream fTable3("Table_NO.3", std::ios::out);
			fTable3.setf(std::ios::scientific);
			fTable3 << "Y_N2 [-] 		E3[-] 		E4[-] 		    E5[-] 		    E6[-] 		    E7[-] 	" << std::endl;
			fTable3 << "5 3" << std::endl;
			fTable3 << std::left << std::setw(16) << "6.00000E-01";
			fTable3 << std::left << std::setw(16) << E3_NO_;
			fTable3 << std::left << std::setw(16) << E4_NO_;
			fTable3 << std::left << std::setw(16) << E5_NO_;
			fTable3 << std::left << std::setw(16) << E6_NO_;
			fTable3 << std::left << std::setw(16) << E7_NO_;
			fTable3 << std::endl;
			fTable3 << std::left << std::setw(16) << "7.00000E-01";
			fTable3 << std::left << std::setw(16) << E3_NO_;
			fTable3 << std::left << std::setw(16) << E4_NO_;
			fTable3 << std::left << std::setw(16) << E5_NO_;
			fTable3 << std::left << std::setw(16) << E6_NO_;
			fTable3 << std::left << std::setw(16) << E7_NO_;
			fTable3 << std::endl;
			fTable3 << std::left << std::setw(16) << "8.00000E-01";
			fTable3 << std::left << std::setw(16) << E3_NO_;
			fTable3 << std::left << std::setw(16) << E4_NO_;
			fTable3 << std::left << std::setw(16) << E5_NO_;
			fTable3 << std::left << std::setw(16) << E6_NO_;
			fTable3 << std::left << std::setw(16) << E7_NO_;
			fTable3 << std::endl;
			fTable3 << "END";
			fTable3 << std::endl;
			fTable3.close();
		}

		{
			std::ofstream fTable4("Table_NO.4", std::ios::out);
			fTable4.setf(std::ios::scientific);
			fTable4 << "Y_N2[-] 		nuF_3 		nuOX_3 		nuW1_4 		    nuNO_5f			nuW2_5f			nuNO_5b			nuW2_5b			nuW3_6			nuW3_7" << std::endl;
			fTable4 << "9 3" << std::endl;
			fTable4 << std::left << std::setw(16) << "6.00000E-01";
			fTable4 << std::left << std::setw(16) << nuF_3_NO_;
			fTable4 << std::left << std::setw(16) << nuOX_3_NO_;
			fTable4 << std::left << std::setw(16) << nuW1_4_NO_;
			fTable4 << std::left << std::setw(16) << nuNO_5f_NO_;
			fTable4 << std::left << std::setw(16) << nuW2_5f_NO_;
			fTable4 << std::left << std::setw(16) << nuNO_5b_NO_;
			fTable4 << std::left << std::setw(16) << nuW2_5b_NO_;
			fTable4 << std::left << std::setw(16) << nuW3_6_NO_;
			fTable4 << std::left << std::setw(16) << nuW3_7_NO_;
			fTable4 << std::endl;
			fTable4 << std::left << std::setw(16) << "7.00000E-01";
			fTable4 << std::left << std::setw(16) << nuF_3_NO_;
			fTable4 << std::left << std::setw(16) << nuOX_3_NO_;
			fTable4 << std::left << std::setw(16) << nuW1_4_NO_;
			fTable4 << std::left << std::setw(16) << nuNO_5f_NO_;
			fTable4 << std::left << std::setw(16) << nuW2_5f_NO_;
			fTable4 << std::left << std::setw(16) << nuNO_5b_NO_;
			fTable4 << std::left << std::setw(16) << nuW2_5b_NO_;
			fTable4 << std::left << std::setw(16) << nuW3_6_NO_;
			fTable4 << std::left << std::setw(16) << nuW3_7_NO_;
			fTable4 << std::endl;
			fTable4 << std::left << std::setw(16) << "8.00000E-01";
			fTable4 << std::left << std::setw(16) << nuF_3_NO_;
			fTable4 << std::left << std::setw(16) << nuOX_3_NO_;
			fTable4 << std::left << std::setw(16) << nuW1_4_NO_;
			fTable4 << std::left << std::setw(16) << nuNO_5f_NO_;
			fTable4 << std::left << std::setw(16) << nuW2_5f_NO_;
			fTable4 << std::left << std::setw(16) << nuNO_5b_NO_;
			fTable4 << std::left << std::setw(16) << nuW2_5b_NO_;
			fTable4 << std::left << std::setw(16) << nuW3_6_NO_;
			fTable4 << std::left << std::setw(16) << nuW3_7_NO_;
			fTable4 << std::endl;
			fTable4 << "END";
			fTable4 << std::endl;
			fTable4.close();
		}

		{
			std::ofstream fTable5("Table_NO.5", std::ios::out);
			fTable5.setf(std::ios::scientific);
			fTable5 << "Y_N2 [-] 		A3[-] 				A4[-] 		    	A5[-] 		    	A6[-] 		   	 	A7[-]" << std::endl;
			fTable5 << "5 3" << std::endl;
			fTable5 << std::left << std::setw(16) << "6.00000E-01";
			fTable5 << std::left << std::setw(16) << A3_NO_orig;
			fTable5 << std::left << std::setw(16) << A4_NO_orig;
			fTable5 << std::left << std::setw(16) << A5_NO_orig;
			fTable5 << std::left << std::setw(16) << A6_NO_orig;
			fTable5 << std::left << std::setw(16) << A7_NO_orig;
			fTable5 << std::endl;
			fTable5 << std::left << std::setw(16) << "7.00000E-01";
			fTable5 << std::left << std::setw(16) << A3_NO_orig;
			fTable5 << std::left << std::setw(16) << A4_NO_orig;
			fTable5 << std::left << std::setw(16) << A5_NO_orig;
			fTable5 << std::left << std::setw(16) << A6_NO_orig;
			fTable5 << std::left << std::setw(16) << A7_NO_orig;
			fTable5 << std::endl;
			fTable5 << std::left << std::setw(16) << "8.00000E-01";
			fTable5 << std::left << std::setw(16) << A3_NO_orig;
			fTable5 << std::left << std::setw(16) << A4_NO_orig;
			fTable5 << std::left << std::setw(16) << A5_NO_orig;
			fTable5 << std::left << std::setw(16) << A6_NO_orig;
			fTable5 << std::left << std::setw(16) << A7_NO_orig;
			fTable5 << std::endl;
			fTable5 << "END";
			fTable5 << std::endl;
			fTable5.close();
		}
	}
}
