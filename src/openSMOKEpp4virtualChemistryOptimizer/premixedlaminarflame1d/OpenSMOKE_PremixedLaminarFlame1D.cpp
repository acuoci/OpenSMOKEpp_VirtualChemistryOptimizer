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

// Include OpenMP Header file
#if defined(_OPENMP)
#include <omp.h>
#endif

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/Utilities.h"
#include "idealreactors/utilities/Grammar_LewisNumbers.h"

// PolimiSoot Analyzer
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Hybrid Method of Moments
#include "utilities/soot/hmom/HMOM.h"

// Virtual Chemistry
#include "utilities/virtualchemistry/VirtualChemistry.h"

// 1D grid
#include "utilities/grids/adaptive/Grid1D.h"
#include "utilities/grids/adaptive/Grammar_Grid1D.h"
#include "utilities/grids/adaptive/Adapter_Grid1D.h"

// Laminar Flame 1D
#include "premixedlaminarflame1d/Grammar_PremixedLaminarFlame1D.h"
#include "premixedlaminarflame1d/OpenSMOKE_PremixedLaminarFlame1D.h"

// License
#include "licensegenerator/license/OpenSMOKELicenseUtilities.hpp"

// Sensitivity Analysis
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap_BlockTridiagonal_SteadyState.h"


int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKE_PremixedLaminarFlame1D", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	unsigned int max_number_allowed_species = 100000;
	OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "PremixedLaminarFlame1D", max_number_allowed_species);

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "PremixedLaminarFlame1D";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKE_PremixedLaminarFlame1D");
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
	OpenSMOKE::Grammar_LaminarFlame grammar_laminarflame;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_laminarflame);

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

		boost::filesystem::path path_input_transport;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Transport", path_input_transport);

		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_kinetics_output);

		OpenSMOKE::RapidKineticMechanismWithTransport(path_kinetics_output, path_input_transport.c_str(), path_input_thermodynamics.c_str(), path_input_kinetics.c_str());
	}

	// Species bundling 
	double species_bundling = 0.;
	if (dictionaries(main_dictionary_name_).CheckOption("@SpeciesBundling") == true)
		dictionaries(main_dictionary_name_).ReadDouble("@SpeciesBundling", species_bundling);

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

		// Transport properties
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		if (species_bundling > 0.)
			transportMapXML->ImportSpeciesBundlingFromXMLFile(doc, species_bundling);

		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

		if (thermodynamicsMapXML->NumberOfSpecies() > max_number_allowed_species)
		{
			std::stringstream tag; tag << max_number_allowed_species;
			OpenSMOKE::FatalErrorMessage("The OpenSMOKE++ license you are using is limited to " + tag.str() + " species");
		}
	}

	std::vector<double> inlet_T;
	std::vector<double> P_Pa;
	double inlet_velocity;
	std::vector<OpenSMOKE::OpenSMOKEVectorDouble> inlet_omega;
	std::vector<double> equivalence_ratios;

	double outlet_T;
	OpenSMOKE::OpenSMOKEVectorDouble outlet_omega;

	// Read inlet conditions
	{
		std::vector<std::string> list_of_strings;
		if (dictionaries(main_dictionary_name_).CheckOption("@InletStream") == true)
			dictionaries(main_dictionary_name_).ReadOption("@InletStream", list_of_strings);

		// If multiple inlet streams are specified
		if (list_of_strings.size() != 1)
		{
			inlet_T.resize(list_of_strings.size());
			P_Pa.resize(list_of_strings.size());
			inlet_omega.resize(list_of_strings.size());
			equivalence_ratios.resize(list_of_strings.size());
			for (unsigned int i = 0; i < list_of_strings.size(); i++)
				GetGasStatusFromDictionary(dictionaries(list_of_strings[i]), *thermodynamicsMapXML, inlet_T[i], P_Pa[i], inlet_omega[i]);
		}
		// If a single inlet stream is defined
		else
		{
			GetGasStatusFromDictionary(dictionaries(list_of_strings[0]), *thermodynamicsMapXML, inlet_T, P_Pa, inlet_omega, equivalence_ratios);
		}
	}
	
	// Read outlet conditions
	{
		double P_Pa_outlet;
		std::string name_of_gas_status_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@OutletStream") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@OutletStream", name_of_gas_status_subdictionary);

		GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML, outlet_T, P_Pa_outlet, outlet_omega);

		if (P_Pa_outlet != P_Pa[0])
			OpenSMOKE::FatalErrorMessage("The pressure of outlet stream does not match with the inlet stream");
	}

	// Read inlet velocity
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@InletVelocity") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@InletVelocity", value, units);
			if (units == "m/s")			inlet_velocity = value;
			else if (units == "cm/s")   inlet_velocity = value/100.;
			else if (units == "mm/s")	inlet_velocity = value/1000.;
			else if (units == "km/h")   inlet_velocity = value*10./36.;
			else OpenSMOKE::FatalErrorMessage("Unknown velocity units");
		}
	}

	Eigen::VectorXd w;
	
	// Adaptive grid
	OpenSMOKE::Grid1D* grid;
	{
		std::string name_of_adaptive_grid_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@Grid") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@Grid", name_of_adaptive_grid_subdictionary);

		grid = new OpenSMOKE::Grid1D(dictionaries(name_of_adaptive_grid_subdictionary), w);
	}

	flame = new OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *grid);
	
	// Output folder
	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
	{
		boost::filesystem::path output_folder;
		dictionaries(main_dictionary_name_).ReadPath("@Output", output_folder);
		flame->SetOutputFolder(output_folder);
	}

	// Solver type
	{
		std::string solver_type;
		if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@Type", solver_type);
			flame->SetSolverType(solver_type);
		}
	}

	// Soret effect
	{
		bool soret = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@Soret") == true)
			dictionaries(main_dictionary_name_).ReadBool("@Soret",soret);
		flame->SetSoret(soret);
	}

	// Simplified transport properties
	{
		bool is_simplfified_transport_properties = false;
		if (dictionaries(main_dictionary_name_).CheckOption("@SimplifiedTransportProperties") == true)
			dictionaries(main_dictionary_name_).ReadBool("@SimplifiedTransportProperties", is_simplfified_transport_properties);
		flame->SetSimplifiedTransportProperties(is_simplfified_transport_properties);
	}

	// Radiative heat transfer
	{
		bool radiative_heat_transfer = false;
		if (dictionaries(main_dictionary_name_).CheckOption("@Radiation") == true)
			dictionaries(main_dictionary_name_).ReadBool("@Radiation", radiative_heat_transfer);
		flame->SetRadiativeHeatTransfer(radiative_heat_transfer);
	}

	// Read environment temperature
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@EnvironmentTemperature", value, units);
			if (units == "K")			value *= 1.;
			else if (units == "C")		value += 273.15;
			else OpenSMOKE::FatalErrorMessage("Unknown temperature units");

			flame->SetEnvironmentTemperature(value);
		}
	}

	// Polimi soot
	OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@PolimiSoot") == true)
		{
			std::string name_of_polimisoot_analyzer_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@PolimiSoot", name_of_polimisoot_analyzer_subdictionary);
			polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML, dictionaries(name_of_polimisoot_analyzer_subdictionary));
			
			if (polimi_soot->number_sections() != 0)
				flame->SetPolimiSoot(polimi_soot);
		}
	}

	// Virtual Chemistry
	OpenSMOKE::VirtualChemistry* virtual_chemistry;
	if (dictionaries(main_dictionary_name_).CheckOption("@VirtualChemistry") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@VirtualChemistry", name_of_subdictionary);
		virtual_chemistry = new OpenSMOKE::VirtualChemistry(*thermodynamicsMapXML, dictionaries(name_of_subdictionary));
		flame->SetVirtualChemistry(virtual_chemistry);
	}

	// On the fly PostProcessing
	OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
	{
		on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, flame->output_folder());

		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
			on_the_fly_post_processing->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			flame->SetOnTheFlyPostProcessing(on_the_fly_post_processing);
		}
	}

	// Fixed temperature profile
	{
		std::string name_of_fixed_temperature_profile_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@FixedTemperatureProfile") == true)
		{
			dictionaries(main_dictionary_name_).ReadDictionary("@FixedTemperatureProfile", name_of_fixed_temperature_profile_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_fixed_temperature_profile_subdictionary), x, y, x_variable, y_variable);

			if (x_variable != "length")
				OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be defined versus spacee");
			if (y_variable != "temperature")
				OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be define the temperature profile");

			flame->SetFixedTemperatureProfile(x, y);
		}
	}

	// Fixed specific (i.e. per unit area) mass flow rate profile
	{
		std::string name_of_fixed_specific_mass_flow_rate_profile_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@FixedSpecificMassFlowRateProfile") == true)
		{
			dictionaries(main_dictionary_name_).ReadDictionary("@FixedSpecificMassFlowRateProfile", name_of_fixed_specific_mass_flow_rate_profile_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_fixed_specific_mass_flow_rate_profile_subdictionary), x, y, x_variable, y_variable);

			if (x_variable != "length")
				OpenSMOKE::FatalErrorMessage("The @FixedSpecificMassFlowRateProfile must be defined versus spacee");
			if (y_variable != "specific-mass-flow-rate")
				OpenSMOKE::FatalErrorMessage("The @FixedSpecificMassFlowRateProfile must be define the specific (i.e. per unit area) mass flow rate profile");

			flame->SetFixedSpecificMassFlowRateProfile(x, y);
		}
	}

	// Lewis numbers
	if (dictionaries(main_dictionary_name_).CheckOption("@LewisNumbers") == true)
	{
		std::string name_of_lewis_numbers_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@LewisNumbers", name_of_lewis_numbers_subdictionary);

		std::vector<double> lewis_numbers;
		OpenSMOKE::GetLewisNumbersFromDictionary(dictionaries(name_of_lewis_numbers_subdictionary), thermodynamicsMapXML[0], lewis_numbers);
		flame->SetLewisNumbers(lewis_numbers);
	}

	// Initialize from backup
	bool use_userdefined_grid_for_backup = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@Backup") == true)
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@DontUseBackupGrid") == true)
			dictionaries(main_dictionary_name_).ReadBool("@DontUseBackupGrid", use_userdefined_grid_for_backup);

		if (flame->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED)
		{
			boost::filesystem::path path_backup;
			dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

			// Set inlet and outlet values and first guess velocity according to user values
			flame->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame->SetOutlet(outlet_T, outlet_omega);
			flame->SetInletVelocity(inlet_velocity);

			// Setup the solution, accordingly to backup file
			flame->InitializeFromBackupFile(path_backup, use_userdefined_grid_for_backup);
		}
		else if (flame->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
		{
			boost::filesystem::path path_backup;
			dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

			// Set inlet and outlet values and first guess velocity according to user values
			flame->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame->SetOutlet(outlet_T, outlet_omega);
			flame->SetInletVelocity(inlet_velocity);

			// Setup the solution, accordingly to backup file
			flame->InitializeFromBackupFile(path_backup, use_userdefined_grid_for_backup);
		}
	}
	else
	{
		if (flame->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED)
		{
			flame->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame->SetOutlet(outlet_T, outlet_omega);
			flame->SetInletVelocity(inlet_velocity);
			flame->SetupForFlameSpeed(w);
		}
		else if (flame->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
		{
			flame->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame->SetOutlet(outlet_T, outlet_omega);
			flame->SetInletVelocity(inlet_velocity);
			flame->SetupForBurnerStabilized(w);
		}
	}

	// Use NLS Solver
	{
		bool use_nls_solver = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@UseNlsSolver") == true)
			dictionaries(main_dictionary_name_).ReadBool("@UseNlsSolver", use_nls_solver);
		flame->SetUseNlsSolver(use_nls_solver);
	}

	// Use DAE Solver
	{
		bool use_dae_solver = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@UseDaeSolver") == true)
			dictionaries(main_dictionary_name_).ReadBool("@UseDaeSolver", use_dae_solver);
		flame->SetUseDaeSolver(use_dae_solver);
	}

	// Sensitivity Options
	OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
	if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
	{
		sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
		std::string name_of_sensitivity_options_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);
		sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
		flame->EnableSensitivityAnalysis(*sensitivity_options);
	}

	// Dae Options
	DaeSMOKE::DaeSolver_Parameters* dae_parameters;
	dae_parameters = new DaeSMOKE::DaeSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@DaeParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@DaeParameters", name_of_subdictionary);
		dae_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}

	// Nls Options
	NlsSMOKE::NonLinearSolver_Parameters* nls_parameters;
	nls_parameters = new NlsSMOKE::NonLinearSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@NlsParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@NlsParameters", name_of_subdictionary);
		nls_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}
	
	// Pseudo Transient Options
	NlsSMOKE::FalseTransientSolver_Parameters* false_transient_parameters;
	false_transient_parameters = new NlsSMOKE::FalseTransientSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@FalseTransientParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@FalseTransientParameters", name_of_subdictionary);
		false_transient_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}

	if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeGasTemperature") == true)
	{
		std::string value;
		dictionaries(main_dictionary_name_).ReadString("@DerivativeGasTemperature", value);
		if (value == "upwind")					flame->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_UPWIND);
		else if (value == "backward")			flame->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_BACKWARD);
		else if (value == "forward")			flame->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_FORWARD);
		else if (value == "centered")			flame->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_CENTERED);
		else OpenSMOKE::FatalErrorMessage("Unknown derivative type for gas temperature");
	}

	if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeGasMassFractions") == true)
	{
		std::string value;
		dictionaries(main_dictionary_name_).ReadString("@DerivativeGasMassFractions", value);
		if (value == "upwind")					flame->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_UPWIND);
		else if (value == "backward")			flame->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_BACKWARD);
		else if (value == "forward")			flame->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_FORWARD);
		else if (value == "centered")			flame->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_CENTERED);
		else OpenSMOKE::FatalErrorMessage("Unknown derivative type for gas mass fractions");
	}

	// Hybrid Method of Moments
	OpenSMOKE::HMOM* hmom;
	if (dictionaries(main_dictionary_name_).CheckOption("@HMOM") == true)
	{
		hmom = new OpenSMOKE::HMOM();
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@HMOM", name_of_subdictionary);
		hmom->SetupFromDictionary(dictionaries(name_of_subdictionary));

		flame->SolveHMOMFromExistingSolution(*hmom, *dae_parameters, *nls_parameters, *false_transient_parameters);

		return(0);
	}

	if (flame->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED)
	{
		// Solve only for a single flame
		if (inlet_omega.size() == 1)
		{
			time_t timerStart;
			time_t timerEnd;
			
			time(&timerStart);
			flame->SolveFlameSpeedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
		// Solve for several flames
		else
		{
			boost::filesystem::path output_folder_root = flame->output_folder();
			std::ofstream fOutput((output_folder_root / "FlameSpeeds.out").string().c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);

			fOutput << std::left;
			fOutput << std::setw(8) << "Case";
			fOutput << std::setw(16) << "Speed[cm/s]";
			fOutput << std::setw(16) << "Eq.Ratio[-]";
			fOutput << std::setw(16) << "Temperature[K]";
			fOutput << std::setw(16) << "Pressure[atm]";
			fOutput << std::endl;

			for (unsigned int i = 0; i < inlet_omega.size(); i++)
			{
				std::stringstream index; index << i;
				std::string case_name = "Case" + index.str();

				flame->SetOutputFolder(output_folder_root / case_name);
				flame->ChangeInletConditions(inlet_T[i], P_Pa[i], inlet_omega[i]);
				flame->SolveFlameSpeedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);

				fOutput << std::setw(8) << std::fixed << std::setprecision(0) << i;
				fOutput << std::setw(16) << std::fixed << std::setprecision(4) << flame->flame_speed()*100.;
				fOutput << std::setw(16) << std::fixed << std::setprecision(4) << equivalence_ratios[i];
				fOutput << std::setw(16) << std::fixed << std::setprecision(2) << inlet_T[i];
				fOutput << std::setw(16) << std::fixed << std::setprecision(4) << P_Pa[i] / 101325.;
				fOutput << std::endl;
			}

			fOutput.close();
		}
	}
	else if (flame->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
	{
		flame->SolveBurnerStabilizedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
	}
	
	return(0);
}
