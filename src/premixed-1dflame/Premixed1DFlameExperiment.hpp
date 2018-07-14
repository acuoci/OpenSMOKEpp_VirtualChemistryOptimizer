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
	void Premixed1DFlameExperiment::Setup(const std::string input_file_name,
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
		OpenSMOKE::Grammar_Premixed1DFlameExperiment grammar_1DFlame;

		// Define the dictionaries
		std::string main_dictionary_name_ = "PremixedLaminarFlame1D";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_1DFlame);

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
				else if (units == "cm/s")   inlet_velocity = value / 100.;
				else if (units == "mm/s")	inlet_velocity = value / 1000.;
				else if (units == "km/h")   inlet_velocity = value * 10. / 36.;
				else OpenSMOKE::FatalErrorMessage("Unknown velocity units");
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
				dictionaries(main_dictionary_name_).ReadBool("@Soret", soret);
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
		flame->SetVirtualChemistry(virtual_chemistry);

		// On the fly PostProcessing
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
		if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
		{
			sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
			std::string name_of_sensitivity_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);
			sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
			flame->EnableSensitivityAnalysis(*sensitivity_options);
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
		if (dictionaries(main_dictionary_name_).CheckOption("@HMOM") == true)
		{
			hmom = new OpenSMOKE::HMOM();
			std::string name_of_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@HMOM", name_of_subdictionary);
			hmom->SetupFromDictionary(dictionaries(name_of_subdictionary));

			flame->SolveHMOMFromExistingSolution(*hmom, *dae_parameters, *nls_parameters, *false_transient_parameters);

//			return(0);
		}

		// Optimization
		{
			optimization_ = new OpenSMOKE::OptimizationRules_Premixed1DFlame;

			std::string name_of_optimization_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@Optimization") == true)
			{
				dictionaries(main_dictionary_name_).ReadDictionary("@Optimization", name_of_optimization_subdictionary);
				optimization_->SetupFromDictionary(dictionaries(name_of_optimization_subdictionary), dictionaries);
			}
		}

		// Solve
		//flame->SolveFlameSpeedFromScratchForOptimization(*dae_parameters, *nls_parameters, *false_transient_parameters);
	}

	void Premixed1DFlameExperiment::Solve(const bool verbose)
	{
		// Solve the flame
		flame->SolveFlameSpeedFromScratchForOptimization(*dae_parameters, *nls_parameters, *false_transient_parameters);
		
		// Extracting the profile
		const unsigned int npoints = flame->T().size();

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
				x(i) = flame->grid().x()[i];
				y(i) = flame->Y()[i](index);
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
				x(i) = flame->grid().x()[i];
				y(i) = flame->T()[i];
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

		// Contribution from velocity
		if (optimization_->is_velocity() == true)
		{
			Eigen::VectorXd x(npoints);
			Eigen::VectorXd y(npoints);
			for (unsigned int i = 0; i < npoints; i++)
			{
				x(i) = flame->grid().x()[i];
				y(i) = flame->v(i);
			}
			FixedProfile profile(npoints, x.data(), y.data());

			// Interpolating the profile
			Eigen::VectorXd y_interpolated(optimization_->u_np());
			profile.InterpolateFree(optimization_->u_x(), y_interpolated);

			// Calculating the errors
			double local_norm2_abs_error_ = 0.;
			double local_norm2_rel_error_ = 0.;
			for (unsigned int i = 0; i < optimization_->u_np(); i++)
			{
				const double abs_error = (optimization_->u_y(i) - y_interpolated(i));
				const double rel_error = (optimization_->u_y(i) - y_interpolated(i)) / optimization_->u_y(i);
				local_norm2_abs_error_ += abs_error * abs_error;
				local_norm2_rel_error_ += rel_error * rel_error;
			}

			local_norm2_abs_error_ /= static_cast<double>(optimization_->u_np());
			local_norm2_rel_error_ /= static_cast<double>(optimization_->u_np());

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