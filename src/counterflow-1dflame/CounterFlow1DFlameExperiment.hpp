/*-----------------------------------------------------------------------*\
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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2018  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

namespace OpenSMOKE
{
	void CounterFlow1DFlameExperiment::Setup(const std::string input_file_name,
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML,
		OpenSMOKE::KineticsMap_CHEMKIN*				kineticsMapXML,
		OpenSMOKE::TransportPropertiesMap_CHEMKIN*	transportMapXML,
		OpenSMOKE::VirtualChemistry*				virtual_chemistry)
	{
		// Pointers
		thermodynamicsMapXML_ = thermodynamicsMapXML;
		kineticsMapXML_ = kineticsMapXML;
		transportMapXML_ = transportMapXML;
		virtual_chemistry_ = virtual_chemistry;

		// Defines the grammar rules
		OpenSMOKE::Grammar_CounterFlow1DFlameExperiment grammar_1DFlame;

		// Define the dictionaries
		std::string main_dictionary_name_ = "CounterFlowFlame1D";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_1DFlame);



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
				else if (units == "cm/s")   fuel_velocity = value / 100.;
				else if (units == "mm/s")	fuel_velocity = value / 1000.;
				else if (units == "km/h")   fuel_velocity = value * 10. / 36.;
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
				else if (units == "km/h")   oxidizer_velocity = value * 10. / 36.;
				else OpenSMOKE::FatalErrorMessage("Unknown oxidizer velocity units");
			}
		}

		Eigen::VectorXd w;

		// Adaptive grid
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
				dictionaries(main_dictionary_name_).ReadBool("@Soret", soret);
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
		flame_cfdf->SetVirtualChemistry(virtual_chemistry);

		// On the fly PostProcessing
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
			if (dictionaries(main_dictionary_name_).CheckOption("@FixedTemperatureProfile") == true &&
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
		if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
		{
			sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
			std::string name_of_sensitivity_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);
			sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
			flame_cfdf->EnableSensitivityAnalysis(*sensitivity_options);
		}

		// Dae Options
		dae_parameters = new DaeSMOKE::DaeSolver_Parameters();
		if (dictionaries(main_dictionary_name_).CheckOption("@DaeParameters") == true)
		{
			std::string name_of_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@DaeParameters", name_of_subdictionary);
			dae_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
		}

		// Nls Options
		nls_parameters = new NlsSMOKE::NonLinearSolver_Parameters();
		if (dictionaries(main_dictionary_name_).CheckOption("@NlsParameters") == true)
		{
			std::string name_of_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@NlsParameters", name_of_subdictionary);
			nls_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
		}

		// Pseudo Transient Options
		false_transient_parameters = new NlsSMOKE::FalseTransientSolver_Parameters();
		if (dictionaries(main_dictionary_name_).CheckOption("@FalseTransientParameters") == true)
		{
			std::string name_of_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@FalseTransientParameters", name_of_subdictionary);
			false_transient_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
		}

		// Hybrid Method of Moments
		if (dictionaries(main_dictionary_name_).CheckOption("@HMOM") == true)
		{
			hmom = new OpenSMOKE::HMOM();
			std::string name_of_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@HMOM", name_of_subdictionary);
			hmom->SetupFromDictionary(dictionaries(name_of_subdictionary));

			flame_cfdf->SolveHMOMFromExistingSolution(*hmom, *dae_parameters, *nls_parameters, *false_transient_parameters);

			//return(0);
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

		// Optimization
		{
			optimization_ = new OpenSMOKE::OptimizationRules_CounterFlow1DFlame;

			std::string name_of_optimization_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@Optimization") == true)
			{
				dictionaries(main_dictionary_name_).ReadDictionary("@Optimization", name_of_optimization_subdictionary);
				optimization_->SetupFromDictionary(dictionaries(name_of_optimization_subdictionary), dictionaries);
			}
		}
	}

	void CounterFlow1DFlameExperiment::Solve(const bool verbose)
	{
		// Solve the flame
		flame_cfdf->SolveFlameFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
		
		// Extracting the profile
		const unsigned int npoints = flame_cfdf->T().size();

		// Errors (total)
		norm2_abs_error_ = 0.;
		norm2_rel_error_ = 0.;

		// Contribution from species
		if (optimization_->is_species() == true)
		{
			const unsigned int index = thermodynamicsMapXML_->IndexOfSpecies(optimization_->species_name()) - 1;

			Eigen::VectorXd x(npoints);
			Eigen::VectorXd y(npoints);
			for (unsigned int i = 0; i < npoints; i++)
			{
				x(i) = flame_cfdf->grid().x()[i];
				y(i) = flame_cfdf->Y()[i](index);
			}
			FixedProfile profile(npoints, x.data(), y.data());

			// Interpolating the profile
			Eigen::VectorXd y_interpolated(optimization_->species_np());
			profile.InterpolateFree(optimization_->species_x(), y_interpolated);

			// Calculating the errors
			double local_norm2_abs_error_ = 0.;
			double local_norm2_rel_error_ = 0.;
			for (unsigned int i = 0; i < optimization_->species_np(); i++)
			{
				const double abs_error = (optimization_->species_y(i) - y_interpolated(i));
				const double rel_error = (optimization_->species_y(i) - y_interpolated(i)) / optimization_->species_y(i);
				local_norm2_abs_error_ += abs_error * abs_error;
				local_norm2_rel_error_ += rel_error * rel_error;
			}

			local_norm2_abs_error_ /= static_cast<double>(optimization_->species_np());
			local_norm2_rel_error_ /= static_cast<double>(optimization_->species_np());

			norm2_abs_error_ += local_norm2_abs_error_;
			norm2_rel_error_ += local_norm2_rel_error_;
		}

		// Contribution from species
		if (optimization_->is_temperature() == true)
		{
			Eigen::VectorXd x(npoints);
			Eigen::VectorXd y(npoints);
			for (unsigned int i = 0; i < npoints; i++)
			{
				x(i) = flame_cfdf->grid().x()[i];
				y(i) = flame_cfdf->T()[i];
			}
			FixedProfile profile(npoints, x.data(), y.data());

			// Interpolating the profile
			Eigen::VectorXd y_interpolated(optimization_->t_np());
			profile.InterpolateFree(optimization_->t_x(), y_interpolated);

			// Calculating the errors
			double local_norm2_abs_error_ = 0.;
			double local_norm2_rel_error_ = 0.;
			for (unsigned int i = 0; i < optimization_->t_np(); i++)
			{
				const double abs_error = (optimization_->t_y(i) - y_interpolated(i));
				const double rel_error = (optimization_->t_y(i) - y_interpolated(i)) / optimization_->t_y(i);
				local_norm2_abs_error_ += abs_error * abs_error;
				local_norm2_rel_error_ += rel_error * rel_error;
			}

			local_norm2_abs_error_ /= static_cast<double>(optimization_->t_np());
			local_norm2_rel_error_ /= static_cast<double>(optimization_->t_np());

			norm2_abs_error_ += local_norm2_abs_error_;
			norm2_rel_error_ += local_norm2_rel_error_;
		}

		if (verbose == true)
		{
			std::cout << "Norm2(abs): " << norm2_abs_error_ << std::endl;
			std::cout << "Norm2(rel): " << norm2_rel_error_ << std::endl;
		}
	}
}