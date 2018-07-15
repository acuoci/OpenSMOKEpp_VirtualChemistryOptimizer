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

// Dynamic boundaries
#include "dynamics/DynamicBoundaries.h"

// Laminar Flame 1D
#include "counterflowflame1d/Grammar_CounterFlowFlame1D.h"
#include "counterflowflame1d/OpenSMOKE_CounterFlowFlame1D.h"

// License
#include "licensegenerator/license/OpenSMOKELicenseUtilities.hpp"

// Sensitivity Analysis
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap_BlockTridiagonal_SteadyState.h"


int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_CounterFlowFlame1D", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	unsigned int max_number_allowed_species = 100000;
	OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "CounterFlowFlame1D", max_number_allowed_species);

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "CounterFlowFlame1D";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_CounterFlowFlame1D");
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
	OpenSMOKE::Grammar_CounterFlowFlame1D grammar_laminarflame;

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
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);

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

	// Read fuel conditions
	std::vector<double> fuel_T;
	std::vector<double> fuel_P_Pa;
	std::vector<OpenSMOKE::OpenSMOKEVectorDouble> fuel_omega;
	{
		std::vector<std::string> list_of_strings;
		if (dictionaries(main_dictionary_name_).CheckOption("@FuelStream") == true)
			dictionaries(main_dictionary_name_).ReadOption("@FuelStream", list_of_strings);

		fuel_T.resize(list_of_strings.size());
		fuel_P_Pa.resize(list_of_strings.size());
		fuel_omega.resize(list_of_strings.size());
		for (unsigned int i = 0; i < list_of_strings.size(); i++)
			GetGasStatusFromDictionary(dictionaries(list_of_strings[i]), *thermodynamicsMapXML, fuel_T[i], fuel_P_Pa[i], fuel_omega[i]);
	}
	
	// Read oxidizer conditions
	std::vector<double> oxidizer_T;
	std::vector<double> oxidizer_P_Pa;
	std::vector<OpenSMOKE::OpenSMOKEVectorDouble> oxidizer_omega;
	{
		std::vector<std::string> list_of_strings;
		if (dictionaries(main_dictionary_name_).CheckOption("@OxidizerStream") == true)
			dictionaries(main_dictionary_name_).ReadOption("@OxidizerStream", list_of_strings);

		{
			oxidizer_T.resize(list_of_strings.size());
			oxidizer_P_Pa.resize(list_of_strings.size());
			oxidizer_omega.resize(list_of_strings.size());
			for (unsigned int i = 0; i < list_of_strings.size(); i++)
				GetGasStatusFromDictionary(dictionaries(list_of_strings[i]), *thermodynamicsMapXML, oxidizer_T[i], oxidizer_P_Pa[i], oxidizer_omega[i]);
		}
	}

	// Check over pressures

	// Read fuel stream velocity
	double fuel_velocity;
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@FuelVelocity") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@FuelVelocity", value, units);
			if (units == "m/s")			fuel_velocity = value;
			else if (units == "cm/s")   fuel_velocity = value/100.;
			else if (units == "mm/s")	fuel_velocity = value/1000.;
			else if (units == "km/h")   fuel_velocity = value*10./36.;
			else OpenSMOKE::FatalErrorMessage("Unknown fuel velocity units");
		}
	}

	// Read oxidizer stream velocity
	double oxidizer_velocity;
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@OxidizerVelocity") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@OxidizerVelocity", value, units);
			if (units == "m/s")			oxidizer_velocity = value;
			else if (units == "cm/s")   oxidizer_velocity = value / 100.;
			else if (units == "mm/s")	oxidizer_velocity = value / 1000.;
			else if (units == "km/h")   oxidizer_velocity = value*10. / 36.;
			else OpenSMOKE::FatalErrorMessage("Unknown oxidizer velocity units");
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

	flame_cfdf = new OpenSMOKE::OpenSMOKE_CounterFlowFlame1D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *grid);
	
	// Output folder
	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
	{
		boost::filesystem::path output_folder;
		dictionaries(main_dictionary_name_).ReadPath("@Output", output_folder);
		flame_cfdf->SetOutputFolder(output_folder);
	}

	// Solver type
	{
		std::string solver_type;
		if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@Type", solver_type);
			flame_cfdf->SetSolverType(solver_type);
		}
	}

	// Enable Soret effect
	{
		bool soret = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@Soret") == true)
			dictionaries(main_dictionary_name_).ReadBool("@Soret",soret);
		flame_cfdf->SetSoret(soret);
	}

	// Simplified transport properties
	{
		bool is_simplfified_transport_properties = false;
		if (dictionaries(main_dictionary_name_).CheckOption("@SimplifiedTransportProperties") == true)
			dictionaries(main_dictionary_name_).ReadBool("@SimplifiedTransportProperties", is_simplfified_transport_properties);
		flame_cfdf->SetSimplifiedTransportProperties(is_simplfified_transport_properties);
	}

	// Radiative heat transfer
	{
		bool radiative_heat_transfer = false;
		if (dictionaries(main_dictionary_name_).CheckOption("@Radiation") == true)
			dictionaries(main_dictionary_name_).ReadBool("@Radiation", radiative_heat_transfer);
		flame_cfdf->SetRadiativeHeatTransfer(radiative_heat_transfer);
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

			flame_cfdf->SetEnvironmentTemperature(value);
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
				flame_cfdf->SetPolimiSoot(polimi_soot);
		}
	}

	// Virtual Chemistry
	OpenSMOKE::VirtualChemistry* virtual_chemistry;
	if (dictionaries(main_dictionary_name_).CheckOption("@VirtualChemistry") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@VirtualChemistry", name_of_subdictionary);
		virtual_chemistry = new OpenSMOKE::VirtualChemistry(*thermodynamicsMapXML, dictionaries(name_of_subdictionary));
		flame_cfdf->SetVirtualChemistry(virtual_chemistry);
	}

	// On the fly PostProcessing
	OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
	{
		on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, flame_cfdf->output_folder());

		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
			on_the_fly_post_processing->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			flame_cfdf->SetOnTheFlyPostProcessing(on_the_fly_post_processing);
		}
	}

	// Check on the temperature profiles (if any)
	{
		if ( dictionaries(main_dictionary_name_).CheckOption("@FixedTemperatureProfile") == true &&
			 dictionaries(main_dictionary_name_).CheckOption("@InitialTemperatureProfile") == true)
			OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile and @InitialTemperatureProfile keywords are mutually exclusive");
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
				OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be defined versus space");
			if (y_variable != "temperature")
				OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be define the temperature profile");

			flame_cfdf->SetFixedTemperatureProfile(x, y);
		}
	}

	// Initial temperature profile
	{
		std::string name_of_fixed_temperature_profile_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@InitialTemperatureProfile") == true)
		{
			dictionaries(main_dictionary_name_).ReadDictionary("@InitialTemperatureProfile", name_of_fixed_temperature_profile_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_fixed_temperature_profile_subdictionary), x, y, x_variable, y_variable);

			if (x_variable != "length")
				OpenSMOKE::FatalErrorMessage("The @InitialTemperatureProfile must be defined versus space");
			if (y_variable != "temperature")
				OpenSMOKE::FatalErrorMessage("The @InitialTemperatureProfile must be define the temperature profile");

			flame_cfdf->SetInitialTemperatureProfile(x, y);
		}
	}

	// Read peak conditions (if any)
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@PeakMixture") == true)
		{
			std::vector<double> peak_T;
			std::vector<double> peak_P_Pa;
			std::vector<OpenSMOKE::OpenSMOKEVectorDouble> peak_omega;
			std::vector<std::string> list_of_strings;

			dictionaries(main_dictionary_name_).ReadOption("@PeakMixture", list_of_strings);

			peak_T.resize(list_of_strings.size());
			peak_P_Pa.resize(list_of_strings.size());
			peak_omega.resize(list_of_strings.size());
			for (unsigned int i = 0; i < list_of_strings.size(); i++)
				GetGasStatusFromDictionary(dictionaries(list_of_strings[i]), *thermodynamicsMapXML, peak_T[i], peak_P_Pa[i], peak_omega[i]);

			flame_cfdf->SetPeakMixture(peak_T[0], peak_omega[0]);
		}
	}

	// Peak position
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@PeakPosition") == true)
		{
			double value;
			std::string units;

			dictionaries(main_dictionary_name_).ReadMeasure("@PeakPosition", value, units);
			if (units == "m")			value = value;
			else if (units == "cm")		value = value / 100.;
			else if (units == "mm")		value = value / 1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown @PeakPosition units");

			flame_cfdf->SetPeakPosition(value);
		}
	}

	// Mixing zone width
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@MixingZoneWidth") == true)
		{
			double value;
			std::string units;

			dictionaries(main_dictionary_name_).ReadMeasure("@MixingZoneWidth", value, units);
			if (units == "m")			value = value;
			else if (units == "cm")		value = value / 100.;
			else if (units == "mm")		value = value / 1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown @MixingZoneWidth units");

			flame_cfdf->SetMixingZoneWidth(value);
		}
	}

	// Initial profile type
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@InitialProfiles") == true)
		{
			std::string value;

			dictionaries(main_dictionary_name_).ReadString("@InitialProfiles", value);

			flame_cfdf->SetInitialProfileType(value);
		}
	}

	// Starting guess of eigenvalue
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@EigenValueStartingGuess") == true)
		{
			double value;
			std::string units;

			dictionaries(main_dictionary_name_).ReadMeasure("@EigenValueStartingGuess", value, units);
			if (units == "kg/m3/s2")		value = value;
			else if (units == "g/cm3/s2")	value = value * 1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown @EigenValueStartingGuess units");

			flame_cfdf->SetEigenValueStartingGuess(value);
		}
	}

	// Radial gradient on the fuel side (G = rho*radial_gradient)
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@RadialGradientFuelSide") == true)
		{
			double value;
			std::string units;

			dictionaries(main_dictionary_name_).ReadMeasure("@RadialGradientFuelSide", value, units);
			if (units == "1/s")			value = value;
			else if (units == "Hz")		value = value;
			else OpenSMOKE::FatalErrorMessage("Unknown @RadialGradientFuelSide units");

			flame_cfdf->SetRadialGradientFuelSide(value);
		}
	}

	// Radial gradient on the oxidizer side (G = rho*radial_gradient)
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@RadialGradientOxidizerSide") == true)
		{
			double value;
			std::string units;

			dictionaries(main_dictionary_name_).ReadMeasure("@RadialGradientOxidizerSide", value, units);
			if (units == "1/s")			value = value;
			else if (units == "Hz")		value = value;
			else OpenSMOKE::FatalErrorMessage("Unknown @RadialGradientOxidizerSide units");

			flame_cfdf->SetRadialGradientOxidizerSide(value);
		}
	}

	// Set planar symmetry
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@PlanarSymmetry") == true)
		{
			bool flag;
			dictionaries(main_dictionary_name_).ReadBool("@PlanarSymmetry", flag);
			flame_cfdf->SetPlanarSymmetry(flag);
		}
	}

	// Set the reaction rate multiplier for all the reactions
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@GasReactionRateMultiplier") == true)
		{
			double value;
			dictionaries(main_dictionary_name_).ReadDouble("@GasReactionRateMultiplier", value);
			flame_cfdf->GasReactionRateMultiplier(value);
		}
	}

	if (dictionaries(main_dictionary_name_).CheckOption("@LewisNumbers") == true)
	{
		std::string name_of_lewis_numbers_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@LewisNumbers", name_of_lewis_numbers_subdictionary);

		std::vector<double> lewis_numbers;
		OpenSMOKE::GetLewisNumbersFromDictionary(dictionaries(name_of_lewis_numbers_subdictionary), thermodynamicsMapXML[0], lewis_numbers);
		flame_cfdf->SetLewisNumbers(lewis_numbers);
	}

	// Dynamic boundaries
	OpenSMOKE::DynamicBoundaries* dynamic_boundaries;
	{
		dynamic_boundaries = new OpenSMOKE::DynamicBoundaries(*thermodynamicsMapXML, flame_cfdf->output_folder());

		if (dictionaries(main_dictionary_name_).CheckOption("@DynamicBoundaries") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@DynamicBoundaries", name_of_options_subdictionary);
			dynamic_boundaries->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			flame_cfdf->SetDynamicBoundaries(dynamic_boundaries);
		}
	}

	// Initialize from backup
	if (dictionaries(main_dictionary_name_).CheckOption("@Backup") == true)
	{
		//else if (flame_cfdf->solver_type() == OpenSMOKE::OpenSMOKE_CounterFlowFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
		{
			boost::filesystem::path path_backup;
			dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

			// Set inlet and outlet values and first guess velocity according to user values
			flame_cfdf->SetFuelSide(fuel_T[0], fuel_P_Pa[0], fuel_omega[0], fuel_velocity);
			flame_cfdf->SetOxidizerSide(oxidizer_T[0], oxidizer_P_Pa[0], oxidizer_omega[0], oxidizer_velocity);

			// Setup the solution, accordingly to backup file
			flame_cfdf->InitializeFromBackupFile(path_backup);
		}
	}
	else
	{
	//	else if (flame_cfdf->solver_type() == OpenSMOKE::OpenSMOKE_CounterFlowFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
		{
			flame_cfdf->SetFuelSide(fuel_T[0], fuel_P_Pa[0], fuel_omega[0], fuel_velocity);
			flame_cfdf->SetOxidizerSide(oxidizer_T[0], oxidizer_P_Pa[0], oxidizer_omega[0], oxidizer_velocity);
			flame_cfdf->SetupForBurnerStabilized();
		}
	}

	// Use DAE Solver
	{
		bool use_dae_solver = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@UseDaeSolver") == true)
			dictionaries(main_dictionary_name_).ReadBool("@UseDaeSolver", use_dae_solver);
		flame_cfdf->SetUseDaeSolver(use_dae_solver);
	}

	// Use NLS Solver
	{
		bool use_nls_solver = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@UseNlsSolver") == true)
			dictionaries(main_dictionary_name_).ReadBool("@UseNlsSolver", use_nls_solver);
		flame_cfdf->SetUseNlsSolver(use_nls_solver);
	}

	// Sensitivity Options
	OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
	if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
	{
		sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
		std::string name_of_sensitivity_options_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);
		sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
		flame_cfdf->EnableSensitivityAnalysis(*sensitivity_options);
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

	// Hybrid Method of Moments
	OpenSMOKE::HMOM* hmom;
	if (dictionaries(main_dictionary_name_).CheckOption("@HMOM") == true)
	{
		hmom = new OpenSMOKE::HMOM();
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@HMOM", name_of_subdictionary);
		hmom->SetupFromDictionary(dictionaries(name_of_subdictionary));

		flame_cfdf->SolveHMOMFromExistingSolution(*hmom, *dae_parameters, *nls_parameters, *false_transient_parameters);

		return(0);
	}

	// Deposition wall in burner stabilized stagnation flames (BSS)
	if (dictionaries(main_dictionary_name_).CheckOption("@DepositionWall") == true)
	{
		bool flag;
		dictionaries(main_dictionary_name_).ReadBool("@DepositionWall", flag);

		if (polimi_soot->thermophoretic_effect() == false)
			OpenSMOKE::FatalErrorMessage("The @DepositionWall option can be used only if the thermophoretic effect is turned on");

		flame_cfdf->SetDepositionWall(flag);
	}

	if (flame_cfdf->solver_type() == OpenSMOKE::OpenSMOKE_CounterFlowFlame1D::SOLVER_TYPE_COUNTERFLOW_DIFFUSION)
	{
		// Solve only for a single flame
		if (fuel_omega.size() == 1)
		{
			time_t timerStart;
			time_t timerEnd;
			
			time(&timerStart);
			flame_cfdf->SolveFlameFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
		// TODO
		/*
		// Solve for several flames
		else
		{
			boost::filesystem::path output_folder_root = flame_cfdf->output_folder();
			std::ofstream fOutput((output_folder_root / "FlameSpeeds.out").string().c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);

			fOutput << std::left;
			fOutput << std::setw(8) << "Case";
			fOutput << std::setw(16) << "Speed[cm/s]";
			fOutput << std::setw(16) << "Eq.Ratio[-]";
			fOutput << std::setw(16) << "Temperature[K]";
			fOutput << std::setw(16) << "Pressure[atm]";
			fOutput << std::endl;

			for (unsigned int i = 0; i < fuel_omega.size(); i++)
			{
				std::stringstream index; index << i;
				std::string case_name = "Case" + index.str();

				flame_cfdf->SetOutputFolder(output_folder_root / case_name);
				flame_cfdf->ChangeFuelSide(fuel_T[i], fuel_P_Pa[i], fuel_omega[i], fuel_velocity);
				flame_cfdf->ChangeOxidizerSide(oxidizer_T[i], oxidizer_P_Pa[i], oxidizer_omega[i], oxidizer_velocity);
				flame_cfdf->SolveFlameFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
				
				fOutput << std::setw(8) << std::fixed << std::setprecision(0) << i;
				fOutput << std::setw(16) << std::fixed << std::setprecision(4) << flame_cfdf->flame_speed()*100.;
				fOutput << std::setw(16) << std::fixed << std::setprecision(4) << equivalence_ratios[i];
				fOutput << std::setw(16) << std::fixed << std::setprecision(2) << inlet_T[i];
				fOutput << std::setw(16) << std::fixed << std::setprecision(4) << P_Pa[i] / 101325.;
				fOutput << std::endl;
				
			}

			fOutput.close();
		}
		*/
	}
	else if (flame_cfdf->solver_type() == OpenSMOKE::OpenSMOKE_CounterFlowFlame1D::SOLVER_TYPE_BURNER_STABILIZED_STAGNATION)
	{
		flame_cfdf->SolveFlameFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
	}
	
	return(0);
}
