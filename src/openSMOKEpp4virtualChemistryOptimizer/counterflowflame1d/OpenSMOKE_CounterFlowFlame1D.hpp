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
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

#include "math/OpenSMOKEUtilities.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap_BlockTridiagonal_SteadyState.h"
#include "utilities/Utilities.h"

OpenSMOKE::OpenSMOKE_CounterFlowFlame1D* flame_cfdf;

// OpenSMOKE++ Solvers
#include "interfaces/Interface_OpenSMOKEppDae.h"
#include "interfaces/Interface_OpenSMOKEppNls.h"
#include "interfaces/Interface_FalseTransient_OpenSMOKEpp.h"

// Interfaces (Sundials)
#if OPENSMOKE_USE_SUNDIALS == 1
#include "interfaces/sundials_header.h"
#include "interfaces/Interface_Ida.h"
#include "interfaces/Interface_KinSol.h"
#include "interfaces/Interface_FalseTransient_KinSol.h"
#endif

// Interfaces (BzzMath)
#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#include "interfaces/Interface_BzzDae.h"
#include "interfaces/Interface_BzzNls.h"
#include "interfaces/Interface_FalseTransient_BzzNls.h"
#endif

// Interfaces (DASPK)
#if OPENSMOKE_USE_DASPK == 1
#include "interfaces/Interface_Daspk.h"
#endif

namespace OpenSMOKE
{
	OpenSMOKE_CounterFlowFlame1D::OpenSMOKE_CounterFlowFlame1D(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
																OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
																OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
																OpenSMOKE::Grid1D& grid) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap),
		transportMap_(transportMap),
		grid_(grid)
	{
		sensitivity_analysis_ = false;
		type_ = SIMULATION_TYPE_Y;
		solver_type_ = SOLVER_TYPE_COUNTERFLOW_DIFFUSION;
		soret_effect_ = false;
		n_steps_video_ = 10;
		n_geometry_ = 3.;

		// Default initial profiles: linear
		T_peak_ = 0.;
		x_peak_ = 0.;
		width_mixing_ = 0.;
		initial_profile_type_ = INITIAL_PROFILE_TYPE_LINEAR;

		// Default starting guess value for H [kg/m3/s2]
		H_starting_guess_ = -100.;

		// Radial gradients (by default they are assumed equal to zero)
		radial_gradient_fuel_ = 0.;
		radial_gradient_oxidizer_ = 0.;

		// Gas reaction rate multiplier
		gas_reaction_rate_multiplier_ = 1.;

		// Diffusion coefficients (by default they are calculated using the molecular theory of gases)
		mass_diffusion_coefficients_type_ = MASS_DIFFUSION_COEFFICIENTS_TYPE_MOLECULAR_THEORY_GASES;

		// Default derivatives for convective terms: upwind
		derivative_type_temperature_ = DERIVATIVE_1ST_UPWIND;
		derivative_type_species_ = DERIVATIVE_1ST_UPWIND;

		// Radiative heat transfer
		radiative_heat_transfer_ = false;
		environment_temperature_ = 298.15;

		count_video_ = n_steps_video_;
		iterations_ = 0;

		output_folder_ = "Output";
		use_dae_solver_ = true;
		use_nls_solver_ = true;

		timeFixedTemperature_ = 1.;
		timeFixedComposition_ = 1.;
		timeComplete_ = 1.;

		is_hmom_soot_ = false;
		is_polimi_soot_ = false;
		is_deposition_wall_ = false;
		soot_deposition_ = 0.;
		is_virtual_chemistry_ = false;
		is_simplified_transport_properties_ = false;

		is_on_the_fly_post_processing_ = false;
		is_userdefined_fixed_temperature_profile_ = false;
		is_userdefined_initial_temperature_profile_ = false;
		is_dynamic_boundaries_ = false;

		MemoryAllocation();
		SetAlgebraicDifferentialEquations();
	}

	void OpenSMOKE_CounterFlowFlame1D::SetSolverType(const std::string solver_type)
	{
		if (solver_type == "CounterFlowDiffusion")				solver_type_ = SOLVER_TYPE_COUNTERFLOW_DIFFUSION;
		else if (solver_type == "BurnerStabilizedStagnation")	solver_type_ = SOLVER_TYPE_BURNER_STABILIZED_STAGNATION;
		else FatalErrorMessage("Unknown solver type. The allowed solver types are: CounterFlowDiffusion | BurnerStabilizedStagnation");
	}

	void OpenSMOKE_CounterFlowFlame1D::SetSoret(const bool soret)
	{
		soret_effect_ = soret;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetSimplifiedTransportProperties(const bool flag)
	{
		is_simplified_transport_properties_ = flag;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetRadiativeHeatTransfer(const bool radiative_heat_transfer)
	{
		radiative_heat_transfer_ = radiative_heat_transfer;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetEnvironmentTemperature(const double environment_temperature)
	{
		environment_temperature_ = environment_temperature;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetPolimiSoot(OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_analyzer)
	{
		polimi_soot_analyzer_ = polimi_soot_analyzer;
		is_polimi_soot_ = true;

		if (polimi_soot_analyzer_->thermophoretic_effect() == true)
		{
			j_thermophoretic_star_.resize(grid_.Np() - 1);
			for (int i = 0; i < grid_.Np() - 1; i++)
				j_thermophoretic_star_[i].resize(polimi_soot_analyzer_->bin_indices().size());
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::SetVirtualChemistry(OpenSMOKE::VirtualChemistry* virtual_chemistry)
	{
		virtual_chemistry_ = virtual_chemistry;
		is_virtual_chemistry_ = true;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetDepositionWall(const bool flag)
	{
		is_deposition_wall_ = flag;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetOnTheFlyPostProcessing(OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing)
	{
		on_the_fly_post_processing_ = on_the_fly_post_processing;
		is_on_the_fly_post_processing_ = true;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetFixedTemperatureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& T)
	{
		is_userdefined_fixed_temperature_profile_ = true;
		is_userdefined_initial_temperature_profile_ = true;
		userdefined_temperature_profile_ = new FixedProfile(x.Size(), x.GetHandle(), T.GetHandle());
	}

	void OpenSMOKE_CounterFlowFlame1D::SetInitialTemperatureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& T)
	{
		is_userdefined_initial_temperature_profile_ = true;
		userdefined_temperature_profile_ = new FixedProfile(x.Size(), x.GetHandle(), T.GetHandle());
	}

	void OpenSMOKE_CounterFlowFlame1D::SetDynamicBoundaries(OpenSMOKE::DynamicBoundaries* dynamic_boundaries)
	{
		dynamic_boundaries_ = dynamic_boundaries;
		is_dynamic_boundaries_ = true;

		if (is_virtual_chemistry_ == true)
			dynamic_boundaries_->SetVirtualChemistry(virtual_chemistry_);
	}

	void OpenSMOKE_CounterFlowFlame1D::SetupForBurnerStabilized()
	{
		// In case the initial/fixed profile for temperature is provided by the user
		if (is_userdefined_initial_temperature_profile_ == true)
			userdefined_temperature_profile_->Interpolate(grid_.x(), T_);

		// Position of stagnation plane
		const double x_stagnation_plane = grid_.L() / (1. + (rho_oxidizer_ / rho_fuel_)*boost::math::pow<2>(v_oxidizer_ / v_fuel_));

		// Position of peaks
		if (x_peak_ == 0.)
		{
			if (is_userdefined_initial_temperature_profile_ == true)
			{
				int iMax;
				T_.maxCoeff(&iMax);
				x_peak_ = grid_.x()(iMax);
			}
			else
			{
				x_peak_ = x_stagnation_plane;
			}
		}

		// Width of mixing region
		if (width_mixing_ == 0.)
			width_mixing_ = 2.*std::min(grid_.L() - x_peak_, x_peak_);

		// First guess composition (triangular profiles)
		if (initial_profile_type_ == INITIAL_PROFILE_TYPE_TRIANGULAR)
		{
			if (T_peak_ == 0.)
				FatalErrorMessage("In order to use the triangular initial profiles, the peak mixture must be defined");

			for (int i = 0; i < grid_.Np(); i++)
				T_(i) = triangular_profile(T_fuel_, T_oxidizer_, T_peak_, x_peak_, width_mixing_, grid_.x()(i));

			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = triangular_profile(Y_fuel_(j), Y_oxidizer_(j), Y_peak_(j), x_peak_, width_mixing_, grid_.x()(i));
		}

		// First guess composition (linear profiles)
		if (initial_profile_type_ == INITIAL_PROFILE_TYPE_LINEAR)
		{
			if (is_userdefined_initial_temperature_profile_ == false)
				FatalErrorMessage("In order to use the linear initial profiles, a user-defined initial temperature profile must be provided");

			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = linear_profile(Y_fuel_(j), Y_oxidizer_(j), x_peak_, width_mixing_, grid_.x()(i));
		}

		// First guess composition (plateau profiles)
		if (initial_profile_type_ == INITIAL_PROFILE_TYPE_PLATEAU)
		{
			if (T_peak_ == 0.)
				FatalErrorMessage("In order to use the triangular initial profiles, the peak mixture must be defined");

			for (int i = 0; i < grid_.Np(); i++)
				T_(i) = plateau_profile(T_fuel_, T_oxidizer_, T_peak_, x_peak_, width_mixing_, grid_.L(), grid_.x()(i));

			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = plateau_profile(Y_fuel_(j), Y_oxidizer_(j), Y_peak_(j), x_peak_, width_mixing_, grid_.L(), grid_.x()(i));
		}

		// Momentum (linear profile)
		for (int i = 0; i < grid_.Np(); i++)
			U_(i) = U_fuel_ + (U_oxidizer_ - U_fuel_)*grid_.x()(i) / grid_.L();

		// Radial gradients G
		const double dU_over_dx = (U_oxidizer_ - U_fuel_) / grid_.L();
		G_(0) = rho_fuel_*radial_gradient_fuel_;
		for (int i = 1; i < grid_.Ni(); i++)
			G_(i) = dU_over_dx;
		G_(grid_.Ni()) = rho_oxidizer_*radial_gradient_oxidizer_;

		// Eigen value H (the starting guess can be chosen by the user)
		H_.setConstant(H_starting_guess_);

		// Calculates properties and diffusion fluxes
		Properties();
		DiffusionFluxes();
	}

	void OpenSMOKE_CounterFlowFlame1D::SetFuelSide(const double Tfuel, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaFuel, const double vFuel)
	{
		Y_fuel_.resize(thermodynamicsMap_.NumberOfSpecies());
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			Y_fuel_(j) = omegaFuel[j+1];

		P_ = P_Pa;

		T_fuel_ = Tfuel;
		v_fuel_ = vFuel;

		// Momentum and boundary conditions
		{
			// Molecular weight
			double mw;
			aux_Y.CopyFrom(Y_fuel_.data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw, aux_Y.GetHandle());

			if (is_virtual_chemistry_ == true)
			{
				mw = virtual_chemistry_->MWMix(aux_Y.GetHandle());
				virtual_chemistry_->MoleFractions(mw, aux_Y.GetHandle(), aux_X.GetHandle());
			}

			// Density
			rho_fuel_ = P_*mw / PhysicalConstants::R_J_kmol / T_fuel_;

			// Momentum variables
			U_fuel_ = rho_fuel_*v_fuel_ / (n_geometry_ - 1.);

			// Boundary conditions for species mass fractions
			boundary_fuel_mass_fractions_.resize(thermodynamicsMap_.NumberOfSpecies());
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				boundary_fuel_mass_fractions_(j) = rho_fuel_*v_fuel_*Y_fuel_(j);
		}
	}
	
	void OpenSMOKE_CounterFlowFlame1D::SetOxidizerSide(const double Toxidizer, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaOxidizer, const double vOxidizer)
	{
		Y_oxidizer_.resize(thermodynamicsMap_.NumberOfSpecies());
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			Y_oxidizer_(j) = omegaOxidizer[j + 1];

		P_ = P_Pa;

		T_oxidizer_ = Toxidizer;
		v_oxidizer_ = vOxidizer;

		// Momentum and boundary conditions
		{
			// Molecular weight
			double mw;
			aux_Y.CopyFrom(Y_oxidizer_.data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw, aux_Y.GetHandle());

			if (is_virtual_chemistry_ == true)
			{
				mw = virtual_chemistry_->MWMix(aux_Y.GetHandle());
				virtual_chemistry_->MoleFractions(mw, aux_Y.GetHandle(), aux_X.GetHandle());
			}

			// Density
			rho_oxidizer_ = P_*mw / PhysicalConstants::R_J_kmol / T_oxidizer_;

			// Momentum variables
			U_oxidizer_ = -rho_oxidizer_*v_oxidizer_ / (n_geometry_ - 1.);

			// Boundary conditions for species mass fractions
			boundary_oxidizer_mass_fractions_.resize(thermodynamicsMap_.NumberOfSpecies());
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				boundary_oxidizer_mass_fractions_(j) = -rho_oxidizer_*v_oxidizer_*Y_oxidizer_(j);
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::SetPeakMixture(const double Tpeak, const OpenSMOKE::OpenSMOKEVectorDouble& omegaPeak)
	{
		Y_peak_.resize(thermodynamicsMap_.NumberOfSpecies());
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			Y_peak_(j) = omegaPeak[j + 1];

		T_peak_ = Tpeak;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetPeakPosition(const double x_peak)
	{
		if (x_peak <= 0. || x_peak >= grid_.L())
			FatalErrorMessage("The peak position is outside the grid boundaries");

		x_peak_ = x_peak;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetMixingZoneWidth(const double width_mixing)
	{
		if (width_mixing <= 0. || width_mixing >= grid_.L())
			FatalErrorMessage("The width of the mixing zone is larger tha the grid");

		width_mixing_ = width_mixing;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetEigenValueStartingGuess(const double H)
	{
		if (H >= 0.)
			FatalErrorMessage("The provided starting guess eigenvalue H must be strictly negative");

		H_starting_guess_ = H;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetRadialGradientFuelSide(const double G)
	{
		radial_gradient_fuel_ = G;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetRadialGradientOxidizerSide(const double G)
	{
		radial_gradient_oxidizer_ = G;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetPlanarSymmetry(const bool flag)
	{
		if (flag == true)	n_geometry_ = 2.;	// planar
		else                n_geometry_ = 3.;	// cylyndrical
	}

	void OpenSMOKE_CounterFlowFlame1D::GasReactionRateMultiplier(const double value)
	{
		if (value < 0.)
			FatalErrorMessage("The provided gas reaction rate multiplier coefficient must be positive or equal to 0");
		gas_reaction_rate_multiplier_ = value;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetLewisNumbers(const std::vector<double> lewis_numbers)
	{
		if (lewis_numbers.size() != thermodynamicsMap_.NumberOfSpecies())
			FatalErrorMessage("The provided Lewis numbers are not correct");

		mass_diffusion_coefficients_type_ = MASS_DIFFUSION_COEFFICIENTS_TYPE_LEWIS_NUMBERS;
		lewis_numbers_ = lewis_numbers;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetInitialProfileType(const std::string type)
	{
		if (type == "linear")
			initial_profile_type_ = INITIAL_PROFILE_TYPE_LINEAR;
		else if (type == "triangular")
			initial_profile_type_ = INITIAL_PROFILE_TYPE_TRIANGULAR;
		else if (type == "plateau")
			initial_profile_type_ = INITIAL_PROFILE_TYPE_PLATEAU;
		else
			FatalErrorMessage("Wrong initial profile type. Allowed types are: triangular | linear | plateau");
	}

	void OpenSMOKE_CounterFlowFlame1D::SetUseNlsSolver(const bool use_nls_solver)
	{
		use_nls_solver_ = use_nls_solver;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetUseDaeSolver(const bool use_dae_solver)
	{
		use_dae_solver_ = use_dae_solver;
	}

	void OpenSMOKE_CounterFlowFlame1D::SetOutputFolder(const boost::filesystem::path output_folder)
	{
		output_folder_ = output_folder;
		if (!boost::filesystem::exists(output_folder_))
			OpenSMOKE::CreateDirectory(output_folder_);
	}

	int OpenSMOKE_CounterFlowFlame1D::BlockDimensions() const
	{
		if (type_ == SIMULATION_TYPE_YTM)
			return thermodynamicsMap_.NumberOfSpecies() + 4;
		else if (type_ == SIMULATION_TYPE_YM)
			return thermodynamicsMap_.NumberOfSpecies() + 3;
		else if (type_ == SIMULATION_TYPE_YT)
			return thermodynamicsMap_.NumberOfSpecies() + 1;
		else if (type_ == SIMULATION_TYPE_Y)
			return thermodynamicsMap_.NumberOfSpecies();
		else if (type_ == SIMULATION_TYPE_TM)
			return 4;
		else if (type_ == SIMULATION_TYPE_M)
			return 3;
		else if (type_ == SIMULATION_TYPE_HMOM)
			return hmom_->n_moments();
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
			return thermodynamicsMap_.NumberOfSpecies() + hmom_->n_moments();
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
			return thermodynamicsMap_.NumberOfSpecies() + 1 + hmom_->n_moments();
		else if (type_ == SIMULATION_TYPE_CSI)
			return 1;
		else
		{
			ErrorMessage("int OpenSMOKE_CounterFlowFlame1D::BlockDimensions() const", "Type not provided");
			return 0;
		}
	}

	int OpenSMOKE_CounterFlowFlame1D::NumberOfEquations() const
	{
		return BlockDimensions()*grid_.Np();
	}

	void OpenSMOKE_CounterFlowFlame1D::MemoryAllocation()
	{
		// Must be allocated only the first time
		if (T_.size() == 0)
		{
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_X, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_Y, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_C, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_R, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_prov, true);
		}

		T_.resize(grid_.Np());
		U_.resize(grid_.Np());
		G_.resize(grid_.Np());
		H_.resize(grid_.Np());

		M_.resize(grid_.Np());
		G_over_rho_.resize(grid_.Np());

		rho_.resize(grid_.Np());
		rho_.resize(grid_.Np());
		mw_.resize(grid_.Np());
		cp_.resize(grid_.Np());
		lambda_.resize(grid_.Np());
		mu_.resize(grid_.Np());
		Q_.resize(grid_.Np());

		cp_species_.resize(grid_.Np());
		for (int i = 0; i <grid_.Np(); i++)
			cp_species_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		gamma_fick_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			gamma_fick_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		gamma_fick_star_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			gamma_fick_star_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		gamma_soret_star_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			gamma_soret_star_[i].resize(transportMap_.iThermalDiffusionRatios().size());

		X_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			X_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		Y_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			Y_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		omega_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			omega_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		j_star_.resize(grid_.Np()-1);
		for (int i = 0; i < grid_.Np()-1; i++)
			j_star_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		j_fick_star_.resize(grid_.Np()-1);
		for (int i = 0; i < grid_.Np()-1; i++)
			j_fick_star_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		j_soret_star_.resize(grid_.Np()-1);
		for (int i = 0; i < grid_.Np()-1; i++)
			j_soret_star_[i].resize(transportMap_.iThermalDiffusionRatios().size());

		if (is_polimi_soot_ == true)
		{
			if (polimi_soot_analyzer_->thermophoretic_effect() == true)
			{
				j_thermophoretic_star_.resize(grid_.Np()-1);
				for (int i = 0; i < grid_.Np()-1; i++)
					j_thermophoretic_star_[i].resize(polimi_soot_analyzer_->bin_indices().size());
			}
		}

		jc_star_.resize(grid_.Np()-1);


		dX_over_dx_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dX_over_dx_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		dY_over_dx_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dY_over_dx_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		dT_over_dx_.resize(grid_.Np());
		dT_over_dx_centered_.resize(grid_.Np());
		lambda_d2T_over_dx2_.resize(grid_.Np());

		dM_over_dx_.resize(grid_.Np());
		dU_over_dx_.resize(grid_.Np());
		mu_d2G_over_rho_over_dx2_.resize(grid_.Np());

		// Time derivatives
		dY_over_dt_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dY_over_dt_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		dT_over_dt_.resize(grid_.Np());
		dU_over_dt_.resize(grid_.Np());
		dG_over_dt_.resize(grid_.Np());
		dH_over_dt_.resize(grid_.Np());

		// Radiative heat transfer (memory allocation)
		Q_radiation_.resize(grid_.Np());
		planck_mean_absorption_gas_.resize(grid_.Np());
		planck_mean_absorption_soot_.resize(grid_.Np());

		// Radiative heat transfer (initialization)
		Q_radiation_.setZero();
		planck_mean_absorption_gas_.setZero();
		planck_mean_absorption_soot_.setZero();

		// Mixture fraction
		csi_.resize(grid_.Np());
		dcsi_over_dt_.resize(grid_.Np());
		dcsi_over_dx_.resize(grid_.Np());
		rho_times_alpha_d2csi_over_dx2_.resize(grid_.Np());
		csi_.setConstant(0.);
		for (int i = 0; i < grid_.Np(); i++)
			csi_(i) = 1. - 1. * grid_.x()[i] / grid_.x()[grid_.Ni()];
	}

	void OpenSMOKE_CounterFlowFlame1D::EnableSensitivityAnalysis(OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		sensitivity_analysis_ = true;
		indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
		for (int i = 0; i<indices_of_sensitivity_species_.size(); i++)
			indices_of_sensitivity_species_[i] = thermodynamicsMap_.IndexOfSpecies(sensitivity_options.list_of_species()[i]) - 1;
	}

	void OpenSMOKE_CounterFlowFlame1D::Properties()
	{
		if (is_virtual_chemistry_ == false)
		{
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Thermodynamics
				{
					thermodynamicsMap_.SetPressure(P_);
					thermodynamicsMap_.SetTemperature(T_(i));

					aux_Y.CopyFrom(Y_[i].data());
					thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
					aux_X.CopyTo(X_[i].data());

					// Concentrations [kmol/m3]
					const double cTot = P_ / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
					Product(cTot, aux_X, &aux_C);

					// Mixture density
					rho_(i) = cTot * mw_(i);	// [kg/m3]

					// Species constant pressure specific heats
					thermodynamicsMap_.cpMolar_Species(aux_prov.GetHandle());

					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						cp_species_[i](j) = aux_prov[j + 1] / thermodynamicsMap_.MW(j); // [J/kg/K]

					// Mixture constant pressure specific heat
					cp_(i) = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(aux_X.GetHandle());
					cp_(i) /= mw_(i); // [J/kg/K]
				}

				// Transport properties
				{
					transportMap_.SetTemperature(T_(i));
					transportMap_.SetPressure(P_);

					// Mixture viscosity
					mu_(i) = transportMap_.DynamicViscosity(aux_X.GetHandle());

					// Mixture thermal conductivity
					lambda_(i) = transportMap_.ThermalConductivity(aux_X.GetHandle());

					if (is_simplified_transport_properties_ == true)
					{
						const double mu0 = 1.8405e-5;	// [kg/m/s]
						const double T0 = 300.;			// [K]
						const double Beta0 = 0.6759;
						const double Pr0 = 0.7;

						mu_(i) = mu0 * std::pow(T_(i) / T0, Beta0);
						lambda_(i) = mu_(i) * cp_(i) / Pr0;
					}

					// Mixture diffusion coefficients
					if (mass_diffusion_coefficients_type_ == MASS_DIFFUSION_COEFFICIENTS_TYPE_MOLECULAR_THEORY_GASES)
					{
						transportMap_.MassDiffusionCoefficients(aux_prov.GetHandle(), aux_X.GetHandle(), transportMap_.is_species_bundling());
						aux_prov.CopyTo(gamma_fick_[i].data());

						// Correct diffusion coefficients for soot particles
						if (is_polimi_soot_ == true)
						{
							if (polimi_soot_analyzer_->physical_diffusion() == true)
							{
								const double Dref = gamma_fick_[i](polimi_soot_analyzer_->bin_diffusivity_reference_species());
								for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
								{
									const unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
									gamma_fick_[i](index) = Dref * polimi_soot_analyzer_->bin_diffusivity_correction_factors()[jj];
								}
							}

							const double reduction_coefficient = polimi_soot_analyzer_->physical_diffusion_reduction_coefficient();
							if (reduction_coefficient != 1.)
							{
								for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
								{
									const unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
									gamma_fick_[i](index) *= reduction_coefficient;
								}
							}
						}

						for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
							gamma_fick_star_[i](j) = gamma_fick_[i](j)*thermodynamicsMap_.MW(j) / mw_(i);
					}
					else if (mass_diffusion_coefficients_type_ == MASS_DIFFUSION_COEFFICIENTS_TYPE_LEWIS_NUMBERS)
					{
						const double alpha = lambda_(i) / rho_(i) / cp_(i);
						for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
							gamma_fick_[i](j) = alpha / lewis_numbers_[j];

						for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
							gamma_fick_star_[i](j) = gamma_fick_[i](j)*thermodynamicsMap_.MW(j) / mw_(i);
					}

					// Thermal diffusion ratios
					if (soret_effect_ == true)
					{
						transportMap_.ThermalDiffusionRatios(aux_prov.GetHandle(), aux_X.GetHandle());
						for (unsigned int ii = 0; ii < transportMap_.iThermalDiffusionRatios().size(); ii++)
						{
							unsigned int index = transportMap_.iThermalDiffusionRatios()[ii];
							gamma_soret_star_[i](ii) = gamma_fick_[i](index - 1) * aux_prov[index] * thermodynamicsMap_.MW(index - 1) / mw_(i) / T_(i);
						}
					}
				}

				// Kinetics
				{
					// Sets the kinetic map
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_);

					// Calculates the reaction and formation rates (molar units)
					kineticsMap_.ReactionRates(aux_C.GetHandle());
					kineticsMap_.FormationRates(aux_R.GetHandle());

					// Corrects the formation rates
					if (gas_reaction_rate_multiplier_ != 1.)
						Product(gas_reaction_rate_multiplier_, &aux_R);

					// Reaction heat
					Q_(i) = kineticsMap_.HeatRelease(aux_R.GetHandle()); // [W/m3]

					// Formation rates in mass units
					ElementByElementProduct(aux_R.Size(), aux_R.GetHandle(), thermodynamicsMap_.MWs().data(), aux_R.GetHandle()); // [kg/m3/s]
					aux_R.CopyTo(omega_[i].data());
				}

				// Radiative heat transfer
				if (radiative_heat_transfer_ == true)
				{
					// Gaseous phase
					{
						planck_mean_absorption_gas_(i) = transportMap_.kPlanckMix(aux_X.GetHandle());
					}

					// Soot particles
					{
						if (is_hmom_soot_ == true)
						{
							if (hmom_->radiative_heat_transfer() == true)
							{
								hmom_->SetNormalizedMoments(hmom_M_[i](0), hmom_M_[i](1), hmom_M_[i](2), hmom_M_[i](3));
								planck_mean_absorption_soot_(i) = hmom_->planck_coefficient(T_(i), hmom_->SootVolumeFraction());
							}
						}
						else if (is_polimi_soot_ == true)
						{
							if (polimi_soot_analyzer_->radiative_heat_transfer() == true)
							{
								Eigen::VectorXd mass_fractions(thermodynamicsMap_.NumberOfSpecies());
								for (unsigned int j = 1; j <= thermodynamicsMap_.NumberOfSpecies(); j++)
									mass_fractions(j - 1) = aux_Y[j];

								planck_mean_absorption_soot_(i) = polimi_soot_analyzer_->planck_coefficient(rho_(i), T_(i), mass_fractions);
							}
						}
					}

					// Total Planck mean absorption coefficient [1/m]
					const double kPlanck = planck_mean_absorption_gas_(i) + planck_mean_absorption_soot_(i);

					// Total radiative heat transfer [W/m3]
					Q_radiation_(i) = 4.*PhysicalConstants::sigmaStefanBoltzmann*kPlanck * (boost::math::pow<4>(T_(i)) - boost::math::pow<4>(environment_temperature_));
				}
			}
		}

		if (is_virtual_chemistry_ == true)
		{
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Thermodynamics
				{
					// Molecular weight [kg/kmol] and mole fractions
					aux_Y.CopyFrom(Y_[i].data());
					mw_(i) = virtual_chemistry_->MWMix(aux_Y.GetHandle());
					virtual_chemistry_->MoleFractions(mw_(i), aux_Y.GetHandle(), aux_X.GetHandle());
					aux_X.CopyTo(X_[i].data());

					// Mixture density [kg/m3]
					rho_(i) = P_ * mw_(i) / PhysicalConstants::R_J_kmol / T_(i);

					// Constant pressure specific heat [J/kg/K]
					cp_(i) = virtual_chemistry_->CpMix(T_(i), P_, aux_Y.GetHandle());
					virtual_chemistry_->CpSpecies(T_(i), P_, cp_species_[i].data());
				}

				// Transport properties
				{
					// Dynamic viscosity [kg/m/s]
					mu_(i) = virtual_chemistry_->DynamicViscosity(T_(i));

					// Thermal conductivity [W/m/K]
					lambda_(i) = virtual_chemistry_->ThermalConductivity(mu_(i), cp_(i));

					// Mass diffusion coefficients [m2/s]
					virtual_chemistry_->MassDiffusionCoefficients(lambda_(i), rho_(i), cp_(i), gamma_fick_[i].data());
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						gamma_fick_star_[i](j) = gamma_fick_[i](j)*virtual_chemistry_->MW(j, aux_Y.GetHandle()) / mw_(i);
				}

				// Kinetics
				{
					// Concentrations [kmol/m3]
					const double cTot = P_ / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]

					// Calculates the formation rates of species [kg/m3/s]
					virtual_chemistry_->FormationRates(cTot, mw_(i), T_(i), aux_Y.GetHandle(), aux_R.GetHandle());
					aux_R.CopyTo(omega_[i].data());

					// Calculates the reaction heat [W/m3]
					Q_(i) = virtual_chemistry_->Qdot(T_(i), P_, aux_R.GetHandle());
				}

				if (radiative_heat_transfer_ == true)
					ErrorMessage("Properties", "The radiative heat transfer cannot be used in conjuction with virtual chemistry");

				if (soret_effect_ == true)
					ErrorMessage("Properties", "The soret effect cannot be used in conjuction with virtual chemistry");

				if (is_polimi_soot_ == true)
					ErrorMessage("Properties", "The POLIMI soot mechanism cannot be used in conjuction with virtual chemistry");
			}
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::DiffusionFluxes()
	{
		// Fick diffusion velocity (TODO)
		{
			grid_.Derivative(DERIVATIVE_1ST_BACKWARD, U_, X_, &dX_over_dx_);
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					j_fick_star_[i](j) = -0.50 * (rho_(i)*gamma_fick_star_[i](j) + rho_(i+1)*gamma_fick_star_[i + 1](j)) *dX_over_dx_[i + 1](j);
			}
		}

		// Soret diffusion velocity
		if (soret_effect_ == true)
		{
			grid_.Derivative(DERIVATIVE_1ST_BACKWARD, U_, T_, &dT_over_dx_);
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int jj = 0; jj < transportMap_.iThermalDiffusionRatios().size(); jj++)
					j_soret_star_[i](jj) = -0.50 * (rho_(i)*gamma_soret_star_[i](jj) + rho_(i+1)*gamma_soret_star_[i + 1](jj)) *dT_over_dx_[i + 1];
			}
		}

		// Thermophoretic diffusion velocity
		if (is_polimi_soot_ == true)
		{
			if (polimi_soot_analyzer_->thermophoretic_effect() == true)
			{
				grid_.Derivative(DERIVATIVE_1ST_BACKWARD, U_, T_, &dT_over_dx_);
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
					{
						unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
						j_thermophoretic_star_[i](jj) = -0.50 * (0.538*mu_(i) / T_(i)*Y_[i](index) + 0.538*mu_(i + 1) / T_(i + 1)*Y_[i + 1](index)) * dT_over_dx_[i + 1];
					}
				}
			}
		}

		// Correction diffusion velocity
		{
			jc_star_.setConstant(0.);

			// Fick contribution
			for (int i = 0; i < grid_.Ni(); i++)
			{
				if (is_virtual_chemistry_ == false)
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						jc_star_(i) -= j_fick_star_[i](j);
				}
				else
				{
					for (unsigned int j = 0; j < virtual_chemistry_->ns_main(); j++)
						jc_star_(i) -= j_fick_star_[i](j);
				}
			}

			// Soret contribution
			if (soret_effect_ == true)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int jj = 0; jj < transportMap_.iThermalDiffusionRatios().size(); jj++)
						jc_star_(i) -= j_soret_star_[i](jj);
				}
			}

			// Thermophoretic contribution
			if (is_polimi_soot_ == true)
			{
				if (polimi_soot_analyzer_->thermophoretic_effect() == true)
				{
					for (int i = 0; i < grid_.Ni(); i++)
					{
						for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
							jc_star_(i) -= j_thermophoretic_star_[i](jj);
					}
				}
			}

			if (is_polimi_soot_ == true)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					double omega_soot = 0.;
					for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
					{
						unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
						omega_soot += Y_[i](index);
					}
					const double omega_gas = 1. - omega_soot;
					jc_star_(i) /= omega_gas;
				}
			}
		}

		// Total diffusion velocity (TOSIMPLIFY)
		{
			// Fick + Correction velocity
			if (is_polimi_soot_ == false)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						j_star_[i](j) = j_fick_star_[i](j) + 0.50*(Y_[i](j) + Y_[i + 1](j))*jc_star_(i);
				}
			}
			else
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					{
						// The correction is applied only to gaseous species
						if (!OpenSMOKE::IsValuePresent(j, polimi_soot_analyzer_->bin_indices()))
							j_star_[i](j) = j_fick_star_[i](j) + 0.50*(Y_[i](j) + Y_[i + 1](j))*jc_star_(i);
						else if (OpenSMOKE::IsValuePresent(j, polimi_soot_analyzer_->bin_indices()))
							j_star_[i](j) = j_fick_star_[i](j);
					}
				}
			}

			// Soret contribution
			if (soret_effect_ == true)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int jj = 0; jj < transportMap_.iThermalDiffusionRatios().size(); jj++)
					{
						unsigned int index = transportMap_.iThermalDiffusionRatios()[jj];
						j_star_[i](index - 1) += j_soret_star_[i](jj);
					}
				}
			}

			// Thermophoretic contribution
			if (is_polimi_soot_ == true)
			{
				if (polimi_soot_analyzer_->thermophoretic_effect() == true)
				{
					for (int i = 0; i < grid_.Ni(); i++)
					{
						for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
						{
							unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
							j_star_[i](index) += j_thermophoretic_star_[i](jj);
						}
					}
				}
			}
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::ResidenceTime(Eigen::VectorXd& tau)
	{
		// Reconstruct the velocity field [m/s]
		Eigen::VectorXd v(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			v(i) = (n_geometry_ - 1.)*U_(i) / rho_(i);

		// Search for the stagnation plane
		unsigned int index_stagnation_plane = 0;
		for (int i = 0; i < grid_.Np(); i++)
			if (v(i) <= 0.)
			{
				index_stagnation_plane = i;
				break;
			}

		// Residence time
		tau.resize(grid_.Np());

		// Fuel side
		tau(0) = 0.;
		for (unsigned int i = 1; i < index_stagnation_plane; i++)
		{
			const double vmean = 0.50*(v(i) + v(i - 1));
			const double deltax = (grid_.x()(i) - grid_.x()(i - 1));
			tau(i) = tau(i - 1) + deltax/vmean;
		}

		// Oxidizer side
		tau(grid_.Np()-1) = 0.;
		for (unsigned int i = grid_.Np()-2; i >= index_stagnation_plane; i--)
		{
			const double vmean = 0.50*(v(i) + v(i+1));
			const double deltax = (grid_.x()(i+1) - grid_.x()(i));
			tau(i) = tau(i+1) + deltax / vmean;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations(const double t, const double* y, double* dy)
	{
		if (type_ == SIMULATION_TYPE_YTM)
			Equations_MassFractions_Temperature_Momentum(t, y, dy);
		if (type_ == SIMULATION_TYPE_YM)
			Equations_MassFractions_Momentum(t, y, dy);
		else if (type_ == SIMULATION_TYPE_YT)
			Equations_MassFractions_Temperature(t, y, dy);
		else if (type_ == SIMULATION_TYPE_Y)
			Equations_MassFractions(t, y, dy);
		else if (type_ == SIMULATION_TYPE_TM)
			Equations_Temperature_MassFlowRate(t, y, dy);
		else if (type_ == SIMULATION_TYPE_M)
			Equations_Momentum(t, y, dy);
		else if (type_ == SIMULATION_TYPE_HMOM)
			Equations_HMOM(t, y, dy);
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
			Equations_MassFractions_HMOM(t, y, dy);
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
			Equations_MassFractions_Temperature_HMOM(t, y, dy);
		else if (type_ == SIMULATION_TYPE_CSI)
			Equations_MixtureFraction(t, y, dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::SubEquations_MassFractions()
	{
		// Inlet boundary
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			dY_over_dt_[0](j) = ((n_geometry_ - 1.)*U_(0)*Y_[0](j) + j_star_[0](j)) - boundary_fuel_mass_fractions_(j);

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			{
				dY_over_dt_[i](j) = -(n_geometry_ - 1.)*U_(i)*dY_over_dx_[i](j)
					- (j_star_[i](j) - j_star_[i - 1](j)) / grid_.dxc_over_2()(i)
					+ omega_[i](j);
				dY_over_dt_[i](j) /= rho_(i);
			}

		// Outlet boundary
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			dY_over_dt_[grid_.Ni()](j) = ((n_geometry_ - 1.)*U_(grid_.Ni())*Y_[grid_.Ni()](j) + j_star_[grid_.Ni() - 1](j)) - boundary_oxidizer_mass_fractions_(j);

		// Outlet boundary (in case of deposition wall)
		if (is_deposition_wall_ == true)
		{
			if (is_polimi_soot_ == true)
			{
				if (polimi_soot_analyzer_->thermophoretic_effect() == true)
				{
					for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
					{
						unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
						dY_over_dt_[grid_.Ni()](index) = ((n_geometry_ - 1.)*U_(grid_.Ni())*Y_[grid_.Ni()](index) + j_star_[grid_.Ni() - 1](index)) - j_thermophoretic_star_[grid_.Ni() - 1](jj);
					}
				}
			}
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::SubEquations_Temperature()
	{
		// Inlet boundary
		dT_over_dt_(0) = T_(0) - T_fuel_;

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			double sumCp = 0.;
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				sumCp += cp_species_[i](j)* (j_star_[i - 1](j) + j_star_[i](j)) / 2.;
			sumCp *= dT_over_dx_centered_[i];

			dT_over_dt_(i) = -(n_geometry_ - 1.)*U_(i)*cp_(i)*dT_over_dx_(i)
				+ lambda_d2T_over_dx2_(i)
				- sumCp
				+ Q_(i);

			if (radiative_heat_transfer_ == true)
				dT_over_dt_(i) += -Q_radiation_(i);

			dT_over_dt_(i) /= (rho_(i) * cp_(i));
		}

		// Outlet boundary
		dT_over_dt_(grid_.Ni()) = T_(grid_.Ni()) - T_oxidizer_;
	}

	void OpenSMOKE_CounterFlowFlame1D::SubEquations_Momentum()
	{
		// Auxiliary fields
		for (int i = 0; i < grid_.Np(); i++)
		{
			G_over_rho_(i) = G_(i) / rho_(i);
			M_(i) = U_(i)*G_over_rho_(i);
		}
			
		// Derivatives
		grid_.Derivative(DERIVATIVE_1ST_UPWIND, U_, M_, &dM_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_BACKWARD, U_, U_, &dU_over_dx_);
		grid_.SecondDerivative(mu_, G_over_rho_, &mu_d2G_over_rho_over_dx2_);

		// Fuel boundary
		dU_over_dt_(0) = U_(0) - U_fuel_;
		dG_over_dt_(0) = rho_fuel_*radial_gradient_fuel_ - G_(0);
		dH_over_dt_(0) = H_(1) - H_(0);

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			dU_over_dt_(i) =	0.50*(G_(i-1)+G_(i)) - dU_over_dx_(i);
			dG_over_dt_(i) =	mu_d2G_over_rho_over_dx2_(i) + 
								n_geometry_*G_(i)*G_over_rho_(i) + H_(i) - (n_geometry_-1.)*dM_over_dx_(i);
			dH_over_dt_(i) =	H_(i+1) - H_(i);
		}

		// Oxidizer boundary
		dU_over_dt_(grid_.Ni()) = 0.50*(G_(grid_.Ni() - 1) + G_(grid_.Ni())) - dU_over_dx_(grid_.Ni());
		dG_over_dt_(grid_.Ni()) = G_(grid_.Ni()) - rho_oxidizer_*radial_gradient_oxidizer_;
		dH_over_dt_(grid_.Ni()) = U_oxidizer_ - U_(grid_.Ni());

		// Soot deposition on right wall (if any)
		SootDepositionOnTheWall();
	}

	void OpenSMOKE_CounterFlowFlame1D::SubEquations_HMOM()
	{
		const int jN2 = thermodynamicsMap_.IndexOfSpecies("N2") - 1;
		const int jOH = thermodynamicsMap_.IndexOfSpecies("OH") - 1;
		const int jH = thermodynamicsMap_.IndexOfSpecies("H") - 1;
		const int jH2O = thermodynamicsMap_.IndexOfSpecies("H2O") - 1;
		const int jH2 = thermodynamicsMap_.IndexOfSpecies("H2") - 1;
		const int jC2H2 = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
		const int jO2 = thermodynamicsMap_.IndexOfSpecies("O2") - 1;
		
		std::vector<int> jPAH(hmom_->pah_species().size());
		for (unsigned int i=0;i<hmom_->pah_species().size();i++)
			jPAH[i] = thermodynamicsMap_.IndexOfSpecies(hmom_->pah_species()[i]) - 1;
		
		// Inlet boundary
		for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			dhmom_M_over_dt_[0](j) = hmom_M_[0](j) - 1e-12;
		
		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			const double gamma = mu_(i)/rho_(i)/hmom_->schmidt_number();	// diffusion coefficient [m2/s]

			const double ctot = P_ / PhysicalConstants::R_J_kmol / T_[i];	// [kmol/m3]
			const double mass_fraction_OH = Y_[i](jOH);
			const double mass_fraction_H = Y_[i](jH);
			const double conc_OH = ctot*X_[i](jOH);
			const double conc_H = ctot*X_[i](jH);
			const double conc_H2O = ctot*X_[i](jH2O);
			const double conc_H2 = ctot*X_[i](jH2);
			const double conc_C2H2 = ctot*X_[i](jC2H2);
			const double conc_O2 = ctot*X_[i](jO2);
			double conc_PAH = 0.;
			for (unsigned int j = 0; j<hmom_->pah_species().size(); j++)
				conc_PAH += ctot*X_[i](jPAH[j]);

			hmom_->SetNormalizedMoments(hmom_M_[i](0), hmom_M_[i](1), hmom_M_[i](2), hmom_M_[i](3));
			hmom_->SetTemperatureAndPressure(T_[i], P_);
			hmom_->SetMassFractions(mass_fraction_OH, mass_fraction_H);
			hmom_->SetConcentrations("kmol/m3", conc_OH, conc_H, conc_H2O, conc_H2, conc_C2H2, conc_O2, conc_PAH);
			hmom_->SetViscosity(SutherlandViscosity(T_[i]));
			hmom_->CalculateSourceMoments();

			const double v = (n_geometry_ - 1.)*U_(i)/rho_(i);		// gas velocity [m/s]
			
			// Equations
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				dhmom_M_over_dt_[i](j) = -v*dhmom_M_over_dx_[i](j) + gamma*d2hmom_M_over_dx2_[i](j) + hmom_->sources()(j);
			
			// Thermophoretic effect
			if (hmom_->thermophoretic_model() != 0)
			{
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dhmom_M_over_dt_[i](j) += hmom_thermophoresis_[i](j) / rho_[i];
			}

			// Source terms gas phase
			if (hmom_->PAHConsumption() == true)
			{
				if (conc_PAH > 1.e-64)
				{
					const double R_PAH = hmom_->PAHConsumptionRate() / 1000.;	// [kmol/m3/s]
					
					double Omega_PAH = 0.;	// [kg/m3]
					for (unsigned int j = 0; j < hmom_->pah_species().size(); j++)
					{
						const double fraction = ctot*X_[i](jPAH[j]) / conc_PAH;
						const double omega = thermodynamicsMap_.MW(jPAH[j]) * R_PAH * fraction;
						dY_over_dt_[i](jPAH[j]) -= omega / rho_[i];
						Omega_PAH += omega;
					}

					dY_over_dt_[i](jN2) += Omega_PAH / rho_[i];
				}
			}
		}
		
		// Outlet boundary
		for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			dhmom_M_over_dt_[grid_.Ni()](j) = hmom_M_[grid_.Ni()](j) - 1e-12;
	}

	void OpenSMOKE_CounterFlowFlame1D::SubEquations_MixtureFraction()
	{
		// Inlet boundary
		dcsi_over_dt_(0) = csi_(0) - 1.;

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			dcsi_over_dt_(i) = -(n_geometry_ - 1.)*U_(i)*dcsi_over_dx_(i) + rho_times_alpha_d2csi_over_dx2_(i);
			dcsi_over_dt_(i) /= rho_(i);
		}

		// Outlet boundary
		dcsi_over_dt_(grid_.Ni()) = csi_(grid_.Ni()) - 0.;
	}

	void OpenSMOKE_CounterFlowFlame1D::Recover_Unknowns(const double* y)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];

				// Temperature
				T_(i) = y[count++];

				// Momentum
				U_(i) = y[count++];
				G_(i) = y[count++];
				H_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];

				// Momentum
				U_(i) = y[count++];
				G_(i) = y[count++];
				H_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];

				// Temperature
				T_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Temperature
				T_(i) = y[count++];

				// Momentum
				U_(i) = y[count++];
				G_(i) = y[count++];
				H_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				U_(i) = y[count++];
				G_(i) = y[count++];
				H_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];
				T_(i) = y[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			for (int i = 0; i < grid_.Np(); i++)
				csi_(i) = y[i];
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::Recover_Residuals(double* dy)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				dy[count++] = dT_over_dt_(i);
				dy[count++] = dU_over_dt_(i);
				dy[count++] = dG_over_dt_(i);
				dy[count++] = dH_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				dy[count++] = dU_over_dt_(i);
				dy[count++] = dG_over_dt_(i);
				dy[count++] = dH_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				dy[count++] = dT_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				dy[count++] = dT_over_dt_(i);
				dy[count++] = dU_over_dt_(i);
				dy[count++] = dG_over_dt_(i);
				dy[count++] = dH_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				dy[count++] = dU_over_dt_(i);
				dy[count++] = dG_over_dt_(i);
				dy[count++] = dH_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dy[count++] = dhmom_M_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dy[count++] = dhmom_M_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				dy[count++] = dT_over_dt_(i);
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dy[count++] = dhmom_M_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			for (int i = 0; i < grid_.Np(); i++)
				dy[i] = dcsi_over_dt_(i);
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MassFractions_Temperature_Momentum(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		{
			// Species
			grid_.Derivative(derivative_type_species_, U_, Y_, &dY_over_dx_);

			// Temperature
			grid_.Derivative(derivative_type_temperature_, U_, T_, &dT_over_dx_);
			grid_.Derivative(DERIVATIVE_1ST_CENTERED, U_, T_, &dT_over_dx_centered_);
			grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);
		}
		
		// Equations
		SubEquations_MassFractions();
		SubEquations_Temperature();
		SubEquations_Momentum();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MassFractions_Momentum(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives of species
		grid_.Derivative(derivative_type_species_, U_, Y_, &dY_over_dx_);

		// Equations
		SubEquations_MassFractions();
		SubEquations_Momentum();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_Momentum(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Equations
		SubEquations_Momentum();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MassFractions_Temperature(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		{
			// Species
			grid_.Derivative(derivative_type_species_, U_, Y_, &dY_over_dx_);

			// Temperature
			grid_.Derivative(derivative_type_temperature_, U_, T_, &dT_over_dx_);
			grid_.Derivative(DERIVATIVE_1ST_CENTERED, U_, T_, &dT_over_dx_centered_);
			grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);
		}

		// Equations
		SubEquations_MassFractions();
		SubEquations_Temperature();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_Temperature_MassFlowRate(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		{
			grid_.Derivative(derivative_type_temperature_, U_, T_, &dT_over_dx_);
			grid_.Derivative(DERIVATIVE_1ST_CENTERED, U_, T_, &dT_over_dx_centered_);
			grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);
		}

		// Equations
		SubEquations_Temperature();
		SubEquations_Momentum();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MassFractions(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives of Species
		grid_.Derivative(derivative_type_species_, U_, Y_, &dY_over_dx_);

		// Equations
		SubEquations_MassFractions();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MassFractions_Temperature_HMOM(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		grid_.Derivative(derivative_type_species_, U_, Y_, &dY_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_UPWIND, U_, hmom_M_, &dhmom_M_over_dx_);
		grid_.SecondDerivative(hmom_M_, &d2hmom_M_over_dx2_);
		grid_.Derivative(derivative_type_temperature_, U_, T_, &dT_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_CENTERED, U_, T_, &dT_over_dx_centered_);
		grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);

		if (hmom_->thermophoretic_model() != 0)
		{
			Eigen::VectorXd d2T_over_dx2(grid_.Np());
			grid_.SecondDerivative(T_, &d2T_over_dx2);
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				for (int i = 0; i < grid_.Np(); i++)
				{
					const double mu = SutherlandViscosity(T_(i));
					hmom_thermophoresis_[i](j) = hmom_M_[i](j)*0.55*mu / T_[i] * d2T_over_dx2[i];
				}
			}
		}

		// Equations
		SubEquations_MassFractions();
		SubEquations_Temperature();
		SubEquations_HMOM();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MassFractions_HMOM(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);
		
		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();
		
		// Derivatives
		grid_.Derivative(derivative_type_species_, U_, Y_, &dY_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_UPWIND, U_, hmom_M_, &dhmom_M_over_dx_);
		grid_.SecondDerivative(hmom_M_, &d2hmom_M_over_dx2_);

		if (hmom_->thermophoretic_model() != 0)
		{
			Eigen::VectorXd d2T_over_dx2(grid_.Np());
			grid_.SecondDerivative(T_, &d2T_over_dx2);
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				for (int i = 0; i < grid_.Np(); i++)
				{
					const double mu = SutherlandViscosity(T_(i));
					hmom_thermophoresis_[i](j) = hmom_M_[i](j)*0.55*mu / T_[i] * d2T_over_dx2[i];
				}
			}
		}

		// Equations
		SubEquations_MassFractions();
		SubEquations_HMOM();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_HMOM(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Derivatives of Species
		grid_.Derivative(DERIVATIVE_1ST_UPWIND, U_, hmom_M_, &dhmom_M_over_dx_);
		grid_.SecondDerivative(hmom_M_, &d2hmom_M_over_dx2_);

		if (hmom_->thermophoretic_model() != 0)
		{
			Eigen::VectorXd d2T_over_dx2(grid_.Np());
			grid_.SecondDerivative(T_, &d2T_over_dx2);
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				for (int i = 0; i < grid_.Np(); i++)
				{
					const double mu = SutherlandViscosity(T_(i));
					hmom_thermophoresis_[i](j) = hmom_M_[i](j)*0.55*mu / T_[i] * d2T_over_dx2[i];
				}
			}
		}

		// Equations
		SubEquations_HMOM();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Equations_MixtureFraction(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		//Properties();

		// Derivatives
		Eigen::VectorXd rho_times_alpha(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			rho_times_alpha(i) = lambda_(i) / cp_(i);

		grid_.Derivative(DERIVATIVE_1ST_UPWIND, U_, csi_, &dcsi_over_dx_);
		grid_.SecondDerivative(rho_times_alpha, csi_, &rho_times_alpha_d2csi_over_dx2_);

		// Equations
		SubEquations_MixtureFraction();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_CounterFlowFlame1D::Print(const double t, const Eigen::VectorXd& y, std::ofstream& fOutput)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = y(count++);
				const double T = y(count++);
				const double m = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				if (is_virtual_chemistry_ == true)
				{
					MW = virtual_chemistry_->MWMix(yy.GetHandle());
					virtual_chemistry_->MoleFractions(MW, yy.GetHandle(), xx.GetHandle());
				}

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()[i];
				fOutput << std::setprecision(9) << std::setw(18) << T;
				fOutput << std::setprecision(9) << std::setw(18) << m;

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum*m/MW;
					}
				}

				// Species
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}

				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		if (type_ == SIMULATION_TYPE_YM)
		{
			// TODO
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{

				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = y(count++);
				const double T = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				if (is_virtual_chemistry_ == true)
				{
					MW = virtual_chemistry_->MWMix(yy.GetHandle());
					virtual_chemistry_->MoleFractions(MW, yy.GetHandle(), xx.GetHandle());
				}

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()[i];
				fOutput << std::setprecision(9) << std::setw(18) << T;
				fOutput << std::setprecision(9) << std::setw(18) << U_(i);

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum * U_(i) / MW;
					}
				}
				
				// Species
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}

				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{

				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				if (is_virtual_chemistry_ == true)
				{
					MW = virtual_chemistry_->MWMix(yy.GetHandle());
					virtual_chemistry_->MoleFractions(MW, yy.GetHandle(), xx.GetHandle());
				}

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()[i];
				fOutput << std::setprecision(9) << std::setw(18) << T_(i);
				fOutput << std::setprecision(9) << std::setw(18) << U_(i);

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum * U_(i) / MW;
					}
				}

				// Species 
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}
				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{

				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = Y_[i](j);
				const double T = y(count++);
				const double m = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				if (is_virtual_chemistry_ == true)
				{
					MW = virtual_chemistry_->MWMix(yy.GetHandle());
					virtual_chemistry_->MoleFractions(MW, yy.GetHandle(), xx.GetHandle());
				}

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()[i];
				fOutput << std::setprecision(9) << std::setw(18) << T;
				fOutput << std::setprecision(9) << std::setw(18) << m;

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum*m/MW;
					}
				}

				// Species
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}
				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			// TODO
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			// TODO
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			// TODO
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			// TODO
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			// TODO
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::Print(const double* y, const double norm_residuals)
	{
		const double max_T = T_.maxCoeff();

		std::cout << std::left << std::setw(14) << std::scientific << std::setprecision(3) << norm_residuals;
		std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << max_T;
		
		std::cout << std::endl;
	}

	void OpenSMOKE_CounterFlowFlame1D::Print(const double t, const double* y)
	{
		if (is_dynamic_boundaries_ == false)
		{
			if (count_video_ == n_steps_video_)
			{
				const double max_T = T_.maxCoeff();

				std::cout << std::left << std::setw(16) << std::scientific << t;
				std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << max_T;
				std::cout << std::endl;

				count_video_ = 0;
			}
		}
		else
		{
			if (iterations_%300 == 0 )
			{
				std::cout << std::left << std::setw(16) << "time[s]";
				std::cout << std::left << std::setw(12) << "Tmax[K]";
				std::cout << std::left << std::setw(12) << "a[1/s]";
				std::cout << std::left << std::setw(12) << "xSt[cm]";
				std::cout << std::left << std::setw(12) << "beta[-]";
				std::cout << std::left << std::setw(12) << "Tfuel[K]";
				std::cout << std::left << std::setw(12) << "Tox[K]";
				std::cout << std::left << std::setw(12) << "vfuel[cm/s]";;
				std::cout << std::left << std::setw(12) << "vox[cm/s]";
				std::cout << std::endl;
			}

			if (count_video_ == dynamic_boundaries_->n_steps_output() )
			{
				const double max_T = T_.maxCoeff();

				std::cout << std::left << std::setw(16) << std::scientific << t;
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << max_T;
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << dynamic_boundaries_->a();
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << dynamic_boundaries_->xst()*100.;
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << dynamic_boundaries_->beta();
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << T_fuel_;
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << T_oxidizer_;
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << v_fuel_ * 100.;
				std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(3) << v_oxidizer_ * 100.;
				std::cout << std::endl;

				dynamic_boundaries_->Print(t, max_T);

				count_video_ = 0;
			}

			UpdateBoundaries(t);
		}

		count_video_++;
		iterations_++;
	}

	void OpenSMOKE_CounterFlowFlame1D::Print(const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "csi[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "v[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "U[kg/m2/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "G[kg/m3/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "H[kg/m3/s2]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "k[W/m/K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "mu[kg/m/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Cp[J/kg/K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Q[W/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qrad[W/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "kPlaGas[1/m]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "kPlaSoot[1/m]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "vTherm[cm/s]", count);

			for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.elements()[j] + "[kmol/m2/s]", count);
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_x", count);
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_w", count);
			fOutput << std::endl;
		}

		if (is_virtual_chemistry_ == false)
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                       Solving Csi                        " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			SetType(SIMULATION_TYPE_CSI);
			const double tEnd = 10.;

			// Solve first the DAE system
			DaeSMOKE::DaeSolver_Parameters dae_parameters;
			int flag = SolveDAE(dae_parameters, tEnd);
		}

		for (int i = 0; i < grid_.Np(); i++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
			OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				yy(j + 1) = Y_[i](j);

			double MW;
			thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

			if (is_virtual_chemistry_ == true)
			{
				MW = virtual_chemistry_->MWMix(yy.GetHandle());
				virtual_chemistry_->MoleFractions(MW, yy.GetHandle(), xx.GetHandle());
			}

			// Calculate thermophoretic velocity [m/s]
			const double vThermophoretic = -0.55*mu_(i) / rho_(i)*dT_over_dx_centered_(i) / T_(i);

			fOutput << std::setprecision(9) << std::setw(20) << csi_(i);
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i]*100.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << (n_geometry_-1.)*U_(i)/rho_(i)*100.;
			fOutput << std::setprecision(9) << std::setw(20) << U_(i);
			fOutput << std::setprecision(9) << std::setw(20) << G_(i);
			fOutput << std::setprecision(9) << std::setw(20) << H_(i);
			fOutput << std::setprecision(9) << std::setw(20) << rho_(i);
			fOutput << std::setprecision(9) << std::setw(20) << lambda_(i);
			fOutput << std::setprecision(9) << std::setw(20) << mu_(i);
			fOutput << std::setprecision(9) << std::setw(20) << cp_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Q_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Q_radiation_(i);
			fOutput << std::setprecision(9) << std::setw(20) << planck_mean_absorption_gas_(i);
			fOutput << std::setprecision(9) << std::setw(20) << planck_mean_absorption_soot_(i);
			fOutput << std::setprecision(9) << std::setw(20) << vThermophoretic*100.;

			// Elements
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				{
					double sum = 0.;
					for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
						sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
					fOutput << std::setprecision(9) << std::setw(20) << sum*U_(i)/MW;	// TODO
				}
			}

			// Sum
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fOutput << std::setprecision(9) << std::setw(20) << xx(j + 1);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fOutput << std::setprecision(9) << std::setw(20) << yy(j + 1);
			}
			fOutput << std::endl;
		}

		fOutput.close();
	}
	
	void OpenSMOKE_CounterFlowFlame1D::PrintSoot(const boost::filesystem::path output_folder)
	{
		// Prepare soot file
		std::ofstream fOutputSoot( (output_folder / "Solution.soot.out").c_str(), std::ios::out);
		fOutputSoot.setf(std::ios::scientific);
		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "csi[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "x[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "v[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "m[kg/m2/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "rho[kg/m3]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "fv(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "x(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "w(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "rho(L)[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "N(L)[#/cm3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "H/C(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/C(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/H(L)[-]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "fv(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "x(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "w(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "rho(S)[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "N(S)[#/cm3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "H/C(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/C(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/H(S)[-]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "depos[kg/m2/s]", count);

			fOutputSoot << std::endl;
		}

		// Prepare soot distribution file
		std::ofstream fOutputSootDistribution((output_folder / "Solution.soot_distribution.out").c_str(), std::ios::out);
		fOutputSootDistribution.setf(std::ios::scientific);
		polimi_soot_analyzer_->WriteDistributionLabel(fOutputSootDistribution);

		// Loop over all the cells
		for (int i = 0; i < grid_.Np(); i++)
		{
			// Gas-phase properties
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << csi_(i);
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 100.;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << T_(i);
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << (n_geometry_ - 1.)*U_(i) / rho_(i)*100.;	
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << U_(i);									
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << rho_(i);

			// Analysis of soot
			polimi_soot_analyzer_->Analysis(T_(i), P_, rho_(i), Y_[i], X_[i]);
			polimi_soot_analyzer_->Distribution();

			// Large sections (soot)
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->fv_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->x_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->omega_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->rho_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->N_large() / 1.e6;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->h_over_c_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_c_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_h_large();

			// Small sections (PAH)
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->fv_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->x_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->omega_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->rho_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->N_small() / 1.e6;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->h_over_c_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_c_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_h_small();

			// Deposition of soot
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << soot_deposition_;
			
			fOutputSoot << std::endl;

			// Distributions
			polimi_soot_analyzer_->WriteDistribution(fOutputSootDistribution, 0., grid_.x()(i), 0., 0., T_(i));
			fOutputSootDistribution << std::endl;
		}

		fOutputSoot.close();
		fOutputSootDistribution.close();
	}

	void OpenSMOKE_CounterFlowFlame1D::PrintHMOM(const boost::filesystem::path output_folder)
	{
		// Prepare soot file
		std::ofstream fOutputHMOM((output_folder / "Solution.hmom.out").c_str(), std::ios::out);
		fOutputHMOM.setf(std::ios::scientific);

		// Labels
		{
			unsigned int count = 1;

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "csi[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "x[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "v[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "m[kg/m2/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "rho[kg/m3]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "fv[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "n[#/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "dp[nm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "dc[nm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "np[nm]", count);

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "M(" + label.str() + ")[mol/m3]";
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Sall(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Snuc(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Sgro(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Soxi(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Scon(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaTot(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDis(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDisSS(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDisSL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDisLL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaCon(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaConSS(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaConSL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaConLL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			fOutputHMOM << std::endl;
		}

		// Indices of relevant species jOH
		const int jOH = thermodynamicsMap_.IndexOfSpecies("OH") - 1;
		const int jH = thermodynamicsMap_.IndexOfSpecies("H") - 1;
		const int jH2O = thermodynamicsMap_.IndexOfSpecies("H2O") - 1;
		const int jH2 = thermodynamicsMap_.IndexOfSpecies("H2") - 1;
		const int jC2H2 = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
		const int jO2 = thermodynamicsMap_.IndexOfSpecies("O2") - 1;
		std::vector<int> jPAH(hmom_->pah_species().size());
		for (unsigned int i = 0; i<hmom_->pah_species().size(); i++)
			jPAH[i] = thermodynamicsMap_.IndexOfSpecies(hmom_->pah_species()[i]) - 1;

		// Loop over all the cells
		for (int i = 0; i < grid_.Np(); i++)
		{
			// Gas-phase properties
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << csi_(i);
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 100.;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << T_(i);
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << (n_geometry_ - 1.)*U_(i) / rho_(i)*100.;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << U_(i);
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << rho_(i);

			// Analysis of soot
			const double ctot = P_ / PhysicalConstants::R_J_kmol / T_[i];	// [kmol/m3]
			const double mass_fraction_OH = Y_[i](jOH);
			const double mass_fraction_H = Y_[i](jH);
			const double conc_OH = ctot*X_[i](jOH);
			const double conc_H = ctot*X_[i](jH);
			const double conc_H2O = ctot*X_[i](jH2O);
			const double conc_H2 = ctot*X_[i](jH2);
			const double conc_C2H2 = ctot*X_[i](jC2H2);
			const double conc_O2 = ctot*X_[i](jO2);
			double conc_PAH = 0.;
			for (unsigned int j = 0; j<hmom_->pah_species().size(); j++)
				conc_PAH += ctot*X_[i](jPAH[j]);

			hmom_->SetNormalizedMoments(hmom_M_[i](0), hmom_M_[i](1), hmom_M_[i](2), hmom_M_[i](3));
			hmom_->SetTemperatureAndPressure(T_[i], P_);
			hmom_->SetMassFractions(mass_fraction_OH, mass_fraction_H);
			hmom_->SetConcentrations("kmol/m3", conc_OH, conc_H, conc_H2O, conc_H2, conc_C2H2, conc_O2, conc_PAH);
			hmom_->SetViscosity(SutherlandViscosity(T_[i]));
			hmom_->CalculateSourceMoments();

			// Soot properties
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootVolumeFraction();
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootParticleNumberDensity();
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootParticleDiameter()*1e9;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootCollisionParticleDiameter()*1e9;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootNumberOfPrimaryParticles();

			// Soot moments
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_M_[i](j);

			// Source terms (overall)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources()(j);

			// Source terms (nucleation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_nucleation()(j);

			// Source terms (surface growth)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_growth()(j);

			// Source terms (oxidation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_oxidation()(j);

			// Source terms (condensation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_condensation()(j);

			// Source terms (coagulation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_overall()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete_ss()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete_sl()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete_ll()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous_ss()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous_sl()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous_ll()(j);

			fOutputHMOM << std::endl;
		}

		fOutputHMOM.close();
	}

	void OpenSMOKE_CounterFlowFlame1D::PrintOnTheFlyPostProcessing()
	{
		// Local residence time
		Eigen::VectorXd tau;
		ResidenceTime(tau);

		// Output files
		std::vector<std::string> additional(1);
		additional[0] = "csi[-]";
		on_the_fly_post_processing_->PrepareOutputFiles(additional);
		for (int i = 0; i < grid_.Np(); i++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble omega(thermodynamicsMap_.NumberOfSpecies());
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				omega[j+1] = Y_[i](j);

			std::vector<double> additional(1);
			additional[0] = csi_(i);

			on_the_fly_post_processing_->WriteOnFile(tau(i), grid_.x()[i], 0., 0., T_(i), P_, omega, additional);
		}
		on_the_fly_post_processing_->CloseOutputFiles();
	}

//	void InitialTemperatureProfile(const double xcen, const double width)
//	{
//
//	}

	void OpenSMOKE_CounterFlowFlame1D::MinimumUnknownsVector(double* v)
	{
		const double minimum_U = -1.e32;				// [kg/m2/s]
		const double minimum_G = -1.e32;				// [kg/m2/s]
		const double minimum_H = -1.e32;				// [kg/m2/s]
		const double minimum_temperature = 100.;		// [K]
		const double zero = 0.;

		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;

				v[count++] = minimum_temperature;
				v[count++] = minimum_U;
				v[count++] = minimum_H;
				v[count++] = minimum_G;
			}
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;

				v[count++] = minimum_U;
				v[count++] = minimum_H;
				v[count++] = minimum_G;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;

				v[count++] = minimum_temperature;
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = minimum_temperature;
				v[count++] = minimum_U;
				v[count++] = minimum_G;
				v[count++] = minimum_H;
			}
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = minimum_U;
				v[count++] = minimum_G;
				v[count++] = minimum_H;
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = zero;
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = zero;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;
				v[count++] = minimum_temperature;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = zero;
			}
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			for (int i = 0; i < grid_.Np(); i++)
				v[i] = zero;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::MaximumUnknownsVector(double* v)
	{
		const double maximum_U = 1e32;				// [kg/m2/s]
		const double maximum_G = 1e32;				// [??]
		const double maximum_H = 1e32;				// [??]
		const double maximum_temperature = 10000.;	// [K]
		const double one = 1.;
		const double big = 1.e16;

		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;

				v[count++] = maximum_temperature;
				v[count++] = maximum_U;
				v[count++] = maximum_G;
				v[count++] = maximum_H;
			}
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;

				v[count++] = maximum_U;
				v[count++] = maximum_G;
				v[count++] = maximum_H;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;

				v[count++] = maximum_temperature;
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = maximum_temperature;
				v[count++] = maximum_U;
				v[count++] = maximum_G;
				v[count++] = maximum_H;
			}
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = maximum_U;
				v[count++] = maximum_G;
				v[count++] = maximum_H;
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = big;
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = big;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;
				v[count++] = maximum_temperature;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = big;
			}
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			for (int i = 0; i < grid_.Np(); i++)
				v[i] = one;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::UnknownsVector(double* v)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);

				v[count++] = T_[i];
				v[count++] = U_[i];
				v[count++] = G_[i];
				v[count++] = H_[i];
			}
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);

				v[count++] = U_[i];
				v[count++] = G_[i];
				v[count++] = H_[i];
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);

				v[count++] = T_[i];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = T_[i];
				v[count++] = U_[i];
				v[count++] = G_[i];
				v[count++] = H_[i];
			}
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = U_[i];
				v[count++] = G_[i];
				v[count++] = H_[i];
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = hmom_M_[i](j);
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = hmom_M_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);
				v[count++] = T_(i);
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = hmom_M_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			for (int i = 0; i < grid_.Np(); i++)
				v[i] = csi_(i);
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::CorrectedUnknownsVector(const double* v)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];

				T_(i) = v[count++];
				U_(i) = v[count++];
				G_(i) = v[count++];
				H_(i) = v[count++];
			}

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];

				U_(i) = v[count++];
				G_(i) = v[count++];
				H_(i) = v[count++];
			}

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];

				T_(i) = v[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				T_(i) = v[count++];
				U_(i) = v[count++];
				G_(i) = v[count++];
				H_(i) = v[count++];
			}

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				U_(i) = v[count++];
				G_(i) = v[count++];
				H_(i) = v[count++];
			}

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = v[count++];

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = v[count++];
			}

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}

		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];
				T_(i) = v[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = v[count++];
			}

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}

		else if (type_ == SIMULATION_TYPE_CSI)
		{
			for (int i = 0; i < grid_.Np(); i++)
				csi_(i) = v[i];

			// Boundary (TODO)
			// m_inlet_ = m_(0);
		}

		// Properties and fluxes
		Properties();
		DiffusionFluxes();
	}

	void OpenSMOKE_CounterFlowFlame1D::SetAlgebraicDifferentialEquations()
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			id_equations_.resize((thermodynamicsMap_.NumberOfSpecies() + 4)*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;

				id_equations_[count++] = true;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_YM)
		{
			id_equations_.resize((thermodynamicsMap_.NumberOfSpecies() + 3)*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				id_equations_[count++] = false;
				id_equations_[count++] = true;
				id_equations_[count++] = false;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			id_equations_.resize((thermodynamicsMap_.NumberOfSpecies() + 1)*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			id_equations_.resize(thermodynamicsMap_.NumberOfSpecies()*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			id_equations_.resize(4*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				id_equations_[count++] = true;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Outlet boundary
			{
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_M)
		{
			id_equations_.resize(3 * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				id_equations_[count++] = false;
				id_equations_[count++] = true;
				id_equations_[count++] = false;
			}

			// Outlet boundary
			{
				id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			id_equations_.resize(hmom_->n_moments() * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				id_equations_[count++] = false;

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;

			// Outlet boundary
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				id_equations_[count++] = false;
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			id_equations_.resize( (hmom_->n_moments()+ thermodynamicsMap_.NumberOfSpecies()) * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			id_equations_.resize((hmom_->n_moments() + 1 + thermodynamicsMap_.NumberOfSpecies()) * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				id_equations_[count++] = true;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_CSI)
		{
			id_equations_.resize(1 * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			id_equations_[count++] = false;

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
				id_equations_[count++] = true;

			// Outlet boundary
			id_equations_[count++] = false;
		}

		unsigned int n_algebraic = std::count_if(id_equations_.begin(), id_equations_.end(), std::bind2nd(std::equal_to<bool>(), false));
		unsigned int n_differential = std::count_if(id_equations_.begin(), id_equations_.end(), std::bind2nd(std::equal_to<bool>(), true));

		algebraic_equations_.resize(n_algebraic);
		differential_equations_.resize(n_differential);

		int count_differential = 0;
		int count_algebraic = 0;
		for (unsigned int i = 0; i < id_equations_.size(); i++)
		{
			if (id_equations_[i] == true)  differential_equations_(count_differential++) = i;
			if (id_equations_[i] == false) algebraic_equations_(count_algebraic++) = i;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::AlgebraicDifferentialVector(double* v)
	{
		int count = 0;
		for (unsigned int i = 0; i < id_equations_.size(); i++)
		{
			if (id_equations_[i] == true)  v[count++] = 1.;
			if (id_equations_[i] == false) v[count++] = 0.;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::CorrectDifferentialEquations(double* upv, double* resv)
	{
		for (int i = 0; i < differential_equations_.size(); i++)
		{
			const int k = differential_equations_[i];
			resv[k] -= upv[k];
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::CorrectAlgebraicEquations(double* yp)
	{
		for (int i = 0; i < algebraic_equations_.size(); i++)
		{
			const int k = algebraic_equations_[i];
			yp[k] = 0.;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::Update(const std::vector<Eigen::VectorXd>& phi)
	{
		MemoryAllocation();

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			for (int j = 0; j < grid_.Np(); j++)
					Y_[j](i) = phi[i](j);
		
		T_ = phi[thermodynamicsMap_.NumberOfSpecies()];
		U_ = phi[thermodynamicsMap_.NumberOfSpecies() + 1];
		G_ = phi[thermodynamicsMap_.NumberOfSpecies() + 2];
		H_ = phi[thermodynamicsMap_.NumberOfSpecies() + 3];

		Properties();
		DiffusionFluxes();
	}

	void OpenSMOKE_CounterFlowFlame1D::SensitivityAnalysis()
	{
		const unsigned int index_temperature = thermodynamicsMap_.NumberOfSpecies();
		const unsigned int index_mass_flow_rate = thermodynamicsMap_.NumberOfSpecies() + 1;

		OpenSMOKE::SensitivityMap_BlockTridiagonal_SteadyState sens(kineticsMap_, NumberOfEquations(), BlockDimensions());
		sens.SetIndexOfTemperature(index_temperature);
		sens.SetIndexOfMassFlowRate(index_mass_flow_rate);
		sens.SetIndicesOfSpecies(indices_of_sensitivity_species_);

		// Scaling factors
		std::vector<Eigen::VectorXd> scaling_Jp(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			scaling_Jp[i].resize(BlockDimensions());	
		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			for (int j = 0; j < grid_.Np(); j++)
				scaling_Jp[j](i) = thermodynamicsMap_.MW(i) / rho_(j);
		for (int j = 0; j < grid_.Np(); j++)
			scaling_Jp[j](index_temperature) = 1. / (rho_(j)*cp_(j));
		for (int j = 0; j < grid_.Np(); j++)
			scaling_Jp[j](index_mass_flow_rate) = 1.;

		// Jacobian matrix
		OpenSMOKE::OpenSMOKEBandMatrixDouble* J;
		J = new OpenSMOKE::OpenSMOKEBandMatrixDouble(NumberOfEquations(), BlockDimensions());
		if (J == NULL)
			OpenSMOKE::FatalErrorMessage("Sensitivity Analysis: Memory allocation for tridiagonal-block matrix");
		J->SetToZero();

		// Calculation of Jacobian matrix
		Jacobian(J);

		// Calculate sensitivity coefficients
		sens.Calculate(T_, P_, X_, *J, scaling_Jp);

		// Destroy the Jacobian matrix
		J->DestroyMat();
		
		// Write sensitivity coefficients
		sens.SaveOnXMLFile(output_folder_.string().c_str());
	}

	// Print XML Files
	void OpenSMOKE_CounterFlowFlame1D::PrintXMLFile(const std::string file_name)
	{
		const unsigned int n_additional = 11;
		
		std::ofstream fXML;
		fXML.open(file_name.c_str(), std::ios::out);
		fXML.setf(std::ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fXML << "<opensmoke version=\"0.1a\">" << std::endl;

		fXML << "<Type> Flame1D </Type>" << std::endl;

		fXML << "<additional>" << std::endl;
		fXML << n_additional << std::endl;
		fXML << "axial-coordinate [cm] 2" << std::endl;
		fXML << "temperature [K] 3" << std::endl;
		fXML << "pressure [Pa] 4" << std::endl;
		fXML << "mol-weight [kg/kmol] 5" << std::endl;
		fXML << "density [kg/m3] 6" << std::endl;
		fXML << "heat-release [W/m3] 7" << std::endl;
		fXML << "axial-velocity [m/s] 8" << std::endl;
		fXML << "mass-flow-rate [kg/m2/s] 9" << std::endl;
		fXML << "G-radial [kg/m3/s] 10" << std::endl;
		fXML << "H-eigenvalue [kg/m3/s2] 11" << std::endl;
		fXML << "csi [-] 12" << std::endl;
		fXML << "</additional>" << std::endl;

		fXML << "<t-p-mw>" << std::endl;
		fXML << "1 2 3" << std::endl;
		fXML << "</t-p-mw>" << std::endl;

		fXML << "<mass-fractions>" << std::endl;
		fXML << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			fXML << thermodynamicsMap_.NamesOfSpecies()[i] << " " << thermodynamicsMap_.MW(i) << " " << n_additional + (i+1) << std::endl;
		fXML << "</mass-fractions>" << std::endl;

		fXML << "<profiles>" << std::endl;
		for (int i = 0; i < grid_.Np(); i++)
		{
			fXML << 1.e2*grid_.x()(i) << " ";
			fXML << T_(i) << " ";
			fXML << P_ << " ";
			fXML << mw_(i) << " ";
			fXML << rho_(i) << " ";
			fXML << Q_(i) << " ";
			fXML << U_(i)/rho_(i) << " ";					// TODO
			fXML << U_(i) << " ";							// TODO
			fXML << G_(i) << " ";							// TODO
			fXML << H_(i) << " ";							// TODO
			fXML << csi_(i) << " ";

			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				fXML << Y_[i](j) << " ";
			fXML << std::endl;
		}
		fXML << "</profiles>" << std::endl;

		fXML << "<profiles-size>" << std::endl;
		fXML << grid_.Np() << " " << thermodynamicsMap_.NumberOfSpecies() + n_additional << std::endl;
		fXML << "</profiles-size>" << std::endl;
		fXML << "</opensmoke>" << std::endl;
		fXML.close();
	}

	void OpenSMOKE_CounterFlowFlame1D::Jacobian(OpenSMOKE::OpenSMOKEBandMatrixDouble* J)
	{
		// Definition of constants
		const int neq = flame_cfdf->NumberOfEquations();
		const int width_ = J->nUpper() + J->nLower() + 1;
		const int ngroups_ = std::min(width_, neq);
		const double ZERO_DER = 1.e-8;
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double DEFAULT_RELATIVE_TOLERANCE = std::sqrt(MachEpsFloat());

		// Memory allocation
		Eigen::VectorXd x(neq);
		Eigen::VectorXd f(neq);
		Eigen::VectorXd f_plus(neq);
		Eigen::VectorXd hJ(neq);

		// Order of magnitude of unknowns
		Eigen::VectorXd x_measure(neq);
		x_measure.setConstant(1.);
		for (int i = 0; i < neq; i++)
		if (x_measure(i) != 0.)
			x_measure(i) = std::fabs(x_measure(i)*DEFAULT_RELATIVE_TOLERANCE);

		// Current solution
		flame_cfdf->UnknownsVector(x.data());
		this->Equations(0., x.data(), f.data());

		// Save the original solution
		Eigen::VectorXd x_plus = x;

		// Loop
		for (int group = 1; group <= ngroups_; group++)
		{
			for (int j = group - 1; j < neq; j += width_)
			{
				const double xh = std::fabs(x(j));
				const double xdh = std::fabs(x_measure(j));
				hJ(j) = ETA2*std::max(xh, xdh);
				hJ(j) = std::max(hJ(j), ZERO_DER);
				hJ(j) = std::min(hJ(j), 0.001 + 0.001*std::fabs(xh));

				x_plus(j) += hJ(j);
			}

			this->Equations(0., x_plus.data(), f_plus.data());

			for (int j = group - 1; j < neq; j += width_)
			{
				x_plus(j) = x(j);

				double* col_j = BAND_COL(J, j);
				int i1 = std::max(0, j - J->nUpper());
				int i2 = std::min(j + J->nLower(), neq - 1);
				for (int i = i1; i <= i2; i++)
					BAND_COL_ELEM(col_j, i, j) = (f_plus(i) - f(i)) / hJ(j);
			}
		}
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveInitialDAE(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double tEnd)
	{
		int flag = -1;

		// Solving the DAE system
		if (dae_parameters.type() != DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
		{
			if (dae_parameters.sparse_linear_algebra() == false)
				flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, dae_parameters, 0., tEnd);
			else
				flag = DaeSMOKE::Solve_Sparse_OpenSMOKEppDae<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, dae_parameters, 0., tEnd);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else // if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
			flag = DaeSMOKE::Solve_TridiagonalBlock_BzzDae<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_CounterFlowFlame1D>(this, dae_object_, dae_parameters, 0., tEnd);
		#endif

		return flag;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveDAE(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double tEnd, const double tStart)
	{
		int flag = -1;

		if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_OPENSMOKEPP)
		{
			if (dae_parameters.sparse_linear_algebra() == false)
				flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, dae_parameters, tStart, tEnd);
			else
				flag = DaeSMOKE::Solve_Sparse_OpenSMOKEppDae<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, dae_parameters, tStart, tEnd);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
			flag = DaeSMOKE::Solve_TridiagonalBlock_BzzDae<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_CounterFlowFlame1D>(this, dae_object_, dae_parameters, tStart, tEnd);
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_IDA)
			flag = DaeSMOKE::Solve_Band_Ida<OpenSMOKE_CounterFlowFlame1D>(this, dae_parameters, tStart, tEnd);
		#endif
		#if OPENSMOKE_USE_DASPK == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_DASPK)
			flag = DaeSMOKE::Solve_Band_Daspk<OpenSMOKE_CounterFlowFlame1D>(this, dae_parameters, tStart, tEnd);
		#endif

		return flag;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveNLS(NlsSMOKE::NonLinearSolver_Parameters& nls_parameters)
	{
		int flag = -1;

		if (nls_parameters.type() == NlsSMOKE::NonLinearSolver_Parameters::NLS_SOLVER_OPENSMOKEPP)
		{
			if (nls_parameters.sparse_linear_algebra() == false)
				flag = NlsSMOKE::Solve_Band_OpenSMOKEppNls<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyNlsSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, nls_parameters);
			else
				flag = NlsSMOKE::Solve_Sparse_OpenSMOKEppNls<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyNlsSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, nls_parameters);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (nls_parameters.type() == NlsSMOKE::NonLinearSolver_Parameters::NLS_SOLVER_BZZNLS)
			flag = NlsSMOKE::Solve_TridiagonalBlock_BzzNls<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyNlsSystem_CounterFlowFlame1D>(this, nls_object_, nls_parameters);
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (nls_parameters.type() == NlsSMOKE::NonLinearSolver_Parameters::NLS_SOLVER_KINSOL)
			flag = NlsSMOKE::Solve_Band_KinSol<OpenSMOKE_CounterFlowFlame1D>(this, nls_parameters);
		#endif

		return flag;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveFalseTransient(NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		if (false_transient_parameters.type() == NlsSMOKE::FalseTransientSolver_Parameters::FALSETRANSIENT_SOLVER_OPENSMOKEPP)
		{
			if (false_transient_parameters.sparse_linear_algebra() == false)
				flag = NlsSMOKE::Solve_Band_OpenSMOKEppFalseTransient<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyFalseTransientSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, false_transient_parameters);
			else
				flag = NlsSMOKE::Solve_Sparse_OpenSMOKEppFalseTransient<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyFalseTransientSystem_OpenSMOKEpp_CounterFlowFlame1D>(this, false_transient_parameters);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (false_transient_parameters.type() == NlsSMOKE::FalseTransientSolver_Parameters::FALSETRANSIENT_SOLVER_BZZNLS)
			flag = NlsSMOKE::Solve_Band_BzzNlsFalseTransient<OpenSMOKE_CounterFlowFlame1D, OpenSMOKE_Flame1D_MyFalseTransientSystem_CounterFlowFlame1D>(this, nls_object_, false_transient_parameters);
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (false_transient_parameters.type() == NlsSMOKE::FalseTransientSolver_Parameters::FALSETRANSIENT_SOLVER_KINSOL)
			flag = NlsSMOKE::Solve_Band_KinSolFalseTransient<OpenSMOKE_CounterFlowFlame1D>(this, false_transient_parameters);
		#endif

		return flag;
	}

	int OpenSMOKE_CounterFlowFlame1D::InitialSolutionCounterFlowDiffusionFlame
		(	DaeSMOKE::DaeSolver_Parameters& dae_parameters,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters	)
	{
		// Solution from scratch:
		// 1. Only momentum
		//    a. NLS: OpenSMOKE++ (default) || BzzNls
		//    b. DAE: OpenSMOKE++ (default) || BzzDae
		// 2. Mass fractions and momentum
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: OpenSMOKE++ (default) || BzzNls
		// 3. Mass fractions, energy, and momentum
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: OpenSMOKE++ (default) || BzzNls

		// Step 1: momentum
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                Initial Solution: Step 1 (M)              " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_M);

			int flag = -1;

			// Solve the NL system
			flag = SolveNLS(nls_parameters);

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeFixedTemperature_);

			// Write solution on ASCII file
			Print((output_folder_ / "Solution.initial.M.out").string().c_str());
		}

		// Step 2: species and momentum
		//if (is_virtual_chemistry_ == false)
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                Initial Solution: Step 2 (Y+M)            " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_YM);

			int flag = -1;

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeFixedTemperature_);

			// Solve the NLS system
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);

			// Write solution on ASCII file
			Print((output_folder_ / "Solution.initial.YM.out").string().c_str());
		}

		// Step 3: species, energy, and temperature
		if (is_userdefined_fixed_temperature_profile_ == false)
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "             Initial Solution: Step 3 (Y+T+M)             " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_YTM);

			int flag = -1;

			// Solve the DAE system
			flag = SolveInitialDAE(dae_parameters, timeComplete_);

			// [DEBUG] If issues were found, try again
			if (flag < 0)
			{
				std::cout << "WARNING: The DAE system was unable to reach a solution..." << std::endl;
				std::cout << "         We try to solve successive DAE systems with increasing final integration times..." << std::endl;

				for (unsigned int k = 1; k <= 5; k++)
				{
					const double tf = timeComplete_ / std::pow(10., 5 + 1 - k);
					std::cout << " * " << k << ". Final time [s]: " << tf << std::endl;
					flag = SolveInitialDAE(dae_parameters, tf);
				}

				if (flag < 0)
				{
					std::cout << "WARNING: The DAE system was unable to reach a solution..." << std::endl;
					std::cout << "         We try to solve successive DAE systems with increasing number of equations..." << std::endl;

					std::cout << " * 1. Fixed mass fractions..." << std::endl;
					SetType(SIMULATION_TYPE_TM);
					flag = SolveInitialDAE(dae_parameters, timeComplete_);

					std::cout << "2. Fixed velocity and temperature..." << std::endl;
					SetType(SIMULATION_TYPE_Y);
					flag = SolveInitialDAE(dae_parameters, timeComplete_);

					std::cout << "3. Complete..." << std::endl;
					SetType(SIMULATION_TYPE_YTM);
					flag = SolveInitialDAE(dae_parameters, timeComplete_);

					if (flag < 0)
					{
						std::cout << "WARNING: The DAE system was unable to reach a solution..." << std::endl;
						std::cout << "         We try to solve successive DAE systems with increasing final integration times..." << std::endl;

						for (unsigned int k = 1; k <= 5; k++)
						{
							const double tf = timeComplete_ / std::pow(10., 5 + 1 - k);
							std::cout << " * " << k << ". Final time [s]: " << tf << std::endl;
							flag = SolveInitialDAE(dae_parameters, tf);
						}
					}

					if (flag < 0)
						OpenSMOKE::FatalErrorMessage("The solver was unable to find a solution on the starting grid. Try to change the setup of you problem.");
				}
			}

			// Then solve the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);

			// Write solution on ASCII file	
			Print((output_folder_ / "Solution.initial.YTM.out").string().c_str());

			// Write soot summary on file
			if (is_polimi_soot_ == true)	
				PrintSoot(output_folder_);

			// On the fly post-processing
			if (is_on_the_fly_post_processing_ == true)
				PrintOnTheFlyPostProcessing();
		}

		return 1;
	}

	int OpenSMOKE_CounterFlowFlame1D::CompleteSolution(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		// Step 2
		//if (is_virtual_chemistry_ == false)
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "               Preparation: Step 2 (T+M)                  " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_TM);

			// Solve first the DAE system
			flag = SolveDAE(dae_parameters, timeFixedComposition_);

			// Solving the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);
		}

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                   Solution (complete)                    " << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		SetType(SIMULATION_TYPE_YTM);

		if (use_dae_solver_ == true)
		{
			// Solve first the DAE system
			flag = SolveDAE(dae_parameters, timeComplete_);

			// [DEBUG] If issues were found, try again
			if (flag < 0)
			{
				std::cout << "WARNING: The DAE system was unable to reach a solution..." << std::endl;
				std::cout << "         We try to solve successive DAE systems with increasing final integration times..." << std::endl;

				for (unsigned int k = 1; k <= 5; k++)
				{
					const double tf = timeComplete_ / std::pow(10., 5 + 1 - k);
					std::cout << " * " << k << ". Final time [s]: " << tf << std::endl;
					flag = SolveDAE(dae_parameters, tf);
				}
			}

			// Then solve the non linear system
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);
		}
		else  // The NLS is solved by definition
		{
			// Then solve the non linear system
			flag = SolveNLS(nls_parameters);

			// In case of failure
			if (flag < 0)
			{
				// Solve the false transient
				SolveFalseTransient(false_transient_parameters);

				// Then solves the non linear system
				flag = SolveNLS(nls_parameters);
			}
		}

		CheckForAdiabaticity();
		AtomicAnalysis();
		CheckForInlet();

		return 1;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveCounterFlowDiffusionFlame(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		// Step 1
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "               Preparation: Step 1 (M)                    " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_M);

			const double tEnd = 1.;

			// Then solve the non linear system
			flag = SolveNLS(nls_parameters);

			// Solve first the DAE system
			flag = SolveDAE(dae_parameters, tEnd);
		}

		// Step 2
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                   Solution (complete)                    " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			
			if (is_userdefined_fixed_temperature_profile_ == false)
				SetType(SIMULATION_TYPE_YTM);
			else
				SetType(SIMULATION_TYPE_YM);

			if (use_dae_solver_ == true)
			{
				// Solve first the DAE system
				flag = SolveDAE(dae_parameters, timeComplete_);

				// [DEBUG] If issues were found, try again
				if (flag < 0)
				{
					std::cout << "WARNING: The DAE system was unable to reach a solution..." << std::endl;
					std::cout << "         We try to solve successive DAE systems with increasing final integration times..." << std::endl;

					for (unsigned int k = 1; k <= 5; k++)
					{
						const double tf = timeComplete_ / std::pow(10., 5 + 1 - k);
						std::cout << " * " << k << ". Final time [s]: " << tf << std::endl;
						flag = SolveDAE(dae_parameters, tf);
					}
				}

				// Then solve the non linear system
				if (use_nls_solver_ == true)
					flag = SolveNLS(nls_parameters);
			}
			else  // The NLS is solved by definition
			{
				// Solve directly the non linear system
				flag = SolveNLS(nls_parameters);

				// In case of failure
				if (flag < 0)
				{
					// Solve the false transient
					SolveFalseTransient(false_transient_parameters);

					// Then solves the non linear system
					flag = SolveNLS(nls_parameters);
				}
			}
		}

		CheckForAdiabaticity();
		AtomicAnalysis();
		CheckForInlet();

		return 1;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveDynamicBoundariesCounterFlowDiffusionFlame(DaeSMOKE::DaeSolver_Parameters& dae_parameters)
	{
		count_video_ = dynamic_boundaries_->n_steps_output();

		dynamic_boundaries_->SetFlame(grid_.L(), P_, n_geometry_);
		dynamic_boundaries_->SetInitialFuelSide(T_fuel_, v_fuel_, Y_fuel_.data());
		dynamic_boundaries_->SetInitialOxidizerSide(T_oxidizer_, v_oxidizer_, Y_oxidizer_.data());
		dynamic_boundaries_->CompleteSetup();
		dynamic_boundaries_->PrepareOutputFiles();
		dynamic_boundaries_->Summary(std::cout);

		// Solve first the DAE system
		for (unsigned int i = 0; i < dynamic_boundaries_->snapshot_list_of_times().size() - 1; i++)
		{
			const double tStart = dynamic_boundaries_->snapshot_list_of_times()[i];
			const double tEnd = dynamic_boundaries_->snapshot_list_of_times()[i+1];

			std::cout << "----------------------------------------------------------"		<< std::endl;
			std::cout << " Solving dynamics from : " << tStart << " to " << tEnd << " s"    << std::endl;
			std::cout << "----------------------------------------------------------"		<< std::endl;

			// Check the temperature profile
			if (is_userdefined_fixed_temperature_profile_ == false)
				SetType(SIMULATION_TYPE_YTM);
			else
				SetType(SIMULATION_TYPE_YM);

			int flag = SolveDAE(dae_parameters, tEnd, tStart);

			// Write current solution
			if (flag >= 0)
			{
				std::stringstream label; label << dynamic_boundaries_->snapshot_list_of_times()[i + 1];
				std::string name_ascii = "Snapshot." + label.str() + ".out";
				Print((output_folder_ / name_ascii).string().c_str());
			}
		}

		return 1;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveHMOMFromExistingSolution(
										OpenSMOKE::HMOM& hmom,
										DaeSMOKE::DaeSolver_Parameters& dae_parameters,
										NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
										NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		// HMOM pointer
		is_hmom_soot_ = true;
		hmom_ = &hmom;

		// -----------------------------------------------------------------------------------
		//				                   Memory allocation
		// -----------------------------------------------------------------------------------

		dhmom_M_over_dx_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dhmom_M_over_dx_[i].resize(hmom_->n_moments());

		dhmom_M_over_dt_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dhmom_M_over_dt_[i].resize(hmom_->n_moments());

		d2hmom_M_over_dx2_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			d2hmom_M_over_dx2_[i].resize(hmom_->n_moments());

		hmom_thermophoresis_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			hmom_thermophoresis_[i].resize(hmom_->n_moments());

		hmom_M_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			hmom_M_[i].resize(hmom_->n_moments());

		// -----------------------------------------------------------------------------------
		//				                   Initial values
		// -----------------------------------------------------------------------------------
		for (int i = 0; i < grid_.Np(); i++)
			hmom_M_[i].setConstant(1e-12);

		// -----------------------------------------------------------------------------------

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                       Solving HMOM                       " << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		
		if (hmom_->PAHConsumption() == false)
		{
			SetType(SIMULATION_TYPE_HMOM);
		}
		else
		{
			SetType(SIMULATION_TYPE_YT_HMOM);
		}

		const double tEnd = 10.;

		// Solve first the DAE system
		int flag = SolveDAE(dae_parameters, tEnd);

		// Print on file
		PrintHMOM(output_folder_);

		// Print main solution on XML file
		PrintXMLFile((output_folder_ / "Output.xml").string().c_str());

		// Write current solution
		Print((output_folder_ / "Solution.final.out").string().c_str());

		// On the fly post-processing
		if (is_on_the_fly_post_processing_ == true)
			PrintOnTheFlyPostProcessing();

		return flag;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_CounterFlowFlame1D::RefineGrid(const unsigned int count)
	{
		OpenSMOKE::Adapter_Grid1D_Status refinement_status;

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                  Grid refinement: " << count               << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 4);
		for (unsigned int i = 0; i < phi.size(); i++)
			phi[i].resize(grid_.Np());

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		for (int j = 0; j < grid_.Np(); j++)
			phi[i](j) = Y_[j](i);
		phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = U_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 2] = G_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 3] = H_;

		std::vector<Eigen::VectorXd> phi_new;
		refinement_status = grid_.Refine(phi, phi_new);

		// If the temperature profile has to be kept fixed, we need to perform
		// interpolation from the user-defined profile, not from the previous profile
		if (is_userdefined_fixed_temperature_profile_ == true && refinement_status != OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
		{
			std::cout << "Interpolating temperature from user-defined profile" << std::endl;
			userdefined_temperature_profile_->Interpolate(grid_.x(), phi_new[thermodynamicsMap_.NumberOfSpecies()]);
		}

		// In case new points have been added
		if (refinement_status != OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
		{
			Update(phi_new);
		}
		
		return refinement_status;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_CounterFlowFlame1D::RefineGrid(const double xA, const double xB, unsigned int count)
	{
		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                  Grid refinement (local): " << count       << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 4);
		for (unsigned int i = 0; i < phi.size(); i++)
			phi[i].resize(grid_.Np());

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		for (int j = 0; j < grid_.Np(); j++)
			phi[i](j) = Y_[j](i);
		phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = U_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = G_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = H_;

		std::vector<Eigen::VectorXd> phi_new;
		grid_.Refine(xA, xB, phi, phi_new);
		Update(phi_new);

		return OpenSMOKE::REGRID_SUCCESS;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_CounterFlowFlame1D::Doubling(unsigned int count)
	{
		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                  Grid doubling: " << count << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 4);
		for (unsigned int i = 0; i < phi.size(); i++)
			phi[i].resize(grid_.Np());

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		for (int j = 0; j < grid_.Np(); j++)
			phi[i](j) = Y_[j](i);
		phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = U_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 2] = G_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 3] = H_;

		std::vector<Eigen::VectorXd> phi_new;
		grid_.Double(phi, phi_new);
		Update(phi_new);

		return OpenSMOKE::REGRID_SUCCESS;
	}

	int OpenSMOKE_CounterFlowFlame1D::SolveFlameFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters
		)
	{
		if (is_dynamic_boundaries_ == false)
		{
			// Initial solution
			InitialSolutionCounterFlowDiffusionFlame(dae_parameters, nls_parameters, false_transient_parameters);

			// Loop
			unsigned int max_refinement_attempts_ = 1000;
			for (unsigned int k = 1; k <= max_refinement_attempts_; k++)
			{
				OpenSMOKE::Adapter_Grid1D_Status refinement_status;
				refinement_status = RefineGrid(k);

				if (refinement_status == OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
				{
					break;
				}
				else
				{
					SolveCounterFlowDiffusionFlame(dae_parameters, nls_parameters, false_transient_parameters);
				}

				// Update statistics
				norm();

				// Write current solution
				Print((output_folder_ / "Solution.current.out").string().c_str());

				// Write soot summary on file
				if (is_polimi_soot_ == true)
					PrintSoot(output_folder_);

				// On the fly post-processing
				if (is_on_the_fly_post_processing_ == true)
					PrintOnTheFlyPostProcessing();

				PrintXMLFile((output_folder_ / "Output.xml").string().c_str());

				if (refinement_status == OpenSMOKE::MAXIMUM_NUMBER_POINTS)
					break;
			}

			// Print final solution
			Print((output_folder_ / "Solution.final.out").string().c_str());

			// Write soot summary on file
			if (is_polimi_soot_ == true)
				PrintSoot(output_folder_);

			// On the fly post-processing
			if (is_on_the_fly_post_processing_ == true)
				PrintOnTheFlyPostProcessing();

			PrintXMLFile((output_folder_ / "Output.xml").string().c_str());

			// Sensitivity Analysis (TODO)
			if (sensitivity_analysis() == true)
			{
				SolveCounterFlowDiffusionFlame(dae_parameters, nls_parameters, false_transient_parameters);
				SensitivityAnalysis();
			}

		}

		if (is_dynamic_boundaries_ == true)
		{
			SolveDynamicBoundariesCounterFlowDiffusionFlame(dae_parameters);
		}
	
		return 1;
	}

	void OpenSMOKE_CounterFlowFlame1D::DiagonalJacobian(const double t, double* y, double* J)
	{
		// Calculated as suggested by Buzzi (private communication)
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;

		double* dy_original = new double[NumberOfEquations()];
		double* dy_plus = new double[NumberOfEquations()];

		double* y_plus = new double[NumberOfEquations()];
		for (int i = 0; i < NumberOfEquations(); i++)
			y_plus[i] = y[i];

		// Fixed values
		Equations(t, y, dy_original);

		// Derivatives with respect to y[kd]
		for (int kd = 0; kd < NumberOfEquations(); kd++)
		{
			double hf = 1.e0;
			double error_weight = 1. / (TOLA + TOLR*std::fabs(y[kd]));
			double hJ = ETA2 * std::fabs(std::max(y[kd], 1. / error_weight));
			double hJf = hf / error_weight;
			hJ = std::max(hJ, hJf);
			hJ = std::max(hJ, ZERO_DER);

			// This is what is done by Buzzi
			double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y[kd]));
			double udy = 1. / dy;
			y_plus[kd] += dy;
			Equations(t, y_plus, dy_plus);

			//for (int j = 1; j <= y.Size(); j++)
			J[kd] = (dy_plus[kd] - dy_original[kd]) * udy;

			y_plus[kd] = y[kd];
		}

		delete dy_original;
		delete dy_plus;
		delete y_plus;
	}

	void OpenSMOKE_CounterFlowFlame1D::DiagonalJacobianForIDA(const double alfa, double* J)
	{
		for (int i = 0; i < differential_equations_.size(); i++)
		{
			const int k = differential_equations_[i];
			J[k] -= alfa*1.;
		}
	}

	// Atomic analysis
	void OpenSMOKE_CounterFlowFlame1D::AtomicAnalysis()
	{
		// TODO
		/*
		// Inlet stream
		std::vector<double> sum_inlet(thermodynamicsMap_.elements().size());
		{
			double MWmix;
			aux_Y.CopyFrom(Y_inlet_.data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, MWmix, aux_Y);

			for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
			{
				sum_inlet[j] = 0.;
				for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
					sum_inlet[j] += thermodynamicsMap_.atomic_composition()(i, j) * aux_X[i + 1];
			}
		}
		
		// Outlet stream
		std::vector<double> sum_final(thermodynamicsMap_.elements().size());
		for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
		{
			sum_final[j] = 0.;
			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
				sum_final[j] += thermodynamicsMap_.atomic_composition()(i, j) * X_[grid_.Ni()](i);
		}

		const double moles_inlet = m_(0) / mw_(0);
		const double moles_final = m_(grid_.Ni()) / mw_(grid_.Ni());

		// Write on the screen
		{

			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout	<< std::setw(16) << std::left << "Atomic balances"
						<< std::setw(15) << std::left << "inlet"
						<< std::setw(15) << std::left << "outlet"
						<< std::setw(15) << std::left << "error(%)"
						<< std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			for (unsigned int j = 0; j<thermodynamicsMap_.elements().size(); j++)
			if (sum_inlet[j] > 0.)
			{
				std::string label = thermodynamicsMap_.elements()[j] + "[kmol/s]";
				std::cout << std::setw(16) << std::left << label
					<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
					<< sum_inlet[j] * moles_inlet
					<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
					<< sum_final[j] * moles_final
					<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
					<< std::fabs((sum_final[j] * moles_final) - (sum_inlet[j] * moles_inlet)) / (sum_inlet[j] * moles_inlet) * 100.
					<< std::endl;
			}

			std::cout << std::endl;
		}
		*/
	}

	void OpenSMOKE_CounterFlowFlame1D::CheckForAdiabaticity()
	{
		// TODO
		/*
		// Velocities and kinetic energy
		const double v_inlet = m_(0) / rho_(0);
		const double v_outlet = m_(grid_.Ni()) / rho_(grid_.Ni());
		const double Ek_inlet = 0.50*m_(0)*v_inlet*v_inlet;
		const double Ek_outlet = 0.50*m_(grid_.Ni())*v_outlet*v_outlet;

		// Enthalpy: inlet 
		double H_inlet;
		{
			// Set thermodynamic map
			thermodynamicsMap_.SetTemperature(T_(0));
			thermodynamicsMap_.SetPressure(P_);

			// Mole fractions and molecular weight
			double MWmix;
			aux_Y.CopyFrom(Y_inlet_.data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, MWmix, aux_Y);

			// Enthalpy [W]
			thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(H_inlet, aux_X);
			H_inlet *= m_(0)/MWmix;
		}

		// Enthalpy: outlet
		double H_outlet;
		{
			// Set thermodynamic map
			thermodynamicsMap_.SetTemperature(T_(grid_.Ni()));
			thermodynamicsMap_.SetPressure(P_);

			// Mole fractions
			aux_X.CopyFrom(X_[grid_.Ni()].data());

			// Enthalpy [W]
			thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(H_outlet, aux_X);
			H_outlet *= m_(grid_.Ni()) / mw_(grid_.Ni());
		}

		// Conductive fluxes [W]
		double Q_inlet  = -lambda_(0)*			(T_(1)-T_(0)) / 
												(grid_.x()(1)-grid_.x()(0));						// implicitly assume area equal to 1 m2
		double Q_outlet = -lambda_(grid_.Ni())*	(T_(grid_.Ni()) - T_(grid_.Ni() - 1)) /
												(grid_.x()(grid_.Ni()) - grid_.x()(grid_.Ni()-1));	// implicitly assume area equal to 1 m2

		// Write on the screen
		{
			std::cout	<< "----------------------------------------------------------" << std::endl;
			std::cout	<< std::setw(16) << std::left << "Energy balances"
						<< std::setw(15) << std::left << "inlet"
						<< std::setw(15) << std::left << "outlet"
						<< std::setw(15) << std::left << "error(%)"
						<< std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			std::cout	<< std::setw(16) << std::left << "Ek[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Ek_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Ek_outlet
						<< std::endl;

			std::cout	<< std::setw(16) << std::left << "Q[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Q_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Q_outlet
						<< std::endl;


			std::cout	<< std::setw(16) << std::left << "H[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_outlet
						<< std::endl;

			std::cout	<< std::setw(16) << std::left << "Etot[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_inlet + Ek_inlet + Q_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_outlet + Ek_outlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< std::fabs((H_inlet + Ek_inlet + Q_inlet) - (H_outlet + Ek_outlet+Q_outlet)) / (H_inlet + Ek_inlet + Q_inlet) * 100.
						<< std::endl;

			std::cout << std::endl;
		}
		*/
	}

	void OpenSMOKE_CounterFlowFlame1D::CheckForInlet()
	{
		// TODO
		/*
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << std::setw(16) << std::left << "Checking inlet"
			<< std::setw(16) << std::left << "nominal"
			<< std::setw(16) << std::left << "current"
			<< std::setw(16) << std::left << "error(%)"
			<< std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		// Check for species
		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		{
			if (Y_inlet_(i) > 0.)
				std::cout	<< std::setw(16) << std::left << thermodynamicsMap_.NamesOfSpecies()[i]
							<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
							<< Y_inlet_(i)
							<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
							<< Y_[0](i)
							<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
							<< std::fabs(Y_inlet_(i) - Y_[0](i)) / Y_inlet_(i) * 100.
							<< std::endl;
		}

		std::cout << std::endl;
		*/
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_CounterFlowFlame1D::Regrid()
	{
		if (grid_.Np() > static_cast<int>(grid_.grid_adapter().regrid_points()))
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                 Regridding...                            " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 4);
			for (unsigned int i = 0; i < phi.size(); i++)
				phi[i].resize(grid_.Np());

			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			for (int j = 0; j < grid_.Np(); j++)
				phi[i](j) = Y_[j](i);
			phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
			phi[thermodynamicsMap_.NumberOfSpecies() + 1] = U_;
			phi[thermodynamicsMap_.NumberOfSpecies() + 2] = G_;
			phi[thermodynamicsMap_.NumberOfSpecies() + 3] = H_;

			std::vector<Eigen::VectorXd> phi_new;
			grid_.Regrid(thermodynamicsMap_.NumberOfSpecies(), phi, phi_new);
			Update(phi_new);
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("Regrid operation cannot be performed because the original number of points is too small.");
		}

		return OpenSMOKE::REGRID_SUCCESS;
	}

	void OpenSMOKE_CounterFlowFlame1D::norm()
	{
		double* y_ = new double[NumberOfEquations()];
		double* f_ = new double[NumberOfEquations()];

		UnknownsVector(y_);
		Equations(0., y_, f_);

		double norm2 = 0.;
		for (int i = 0; i < NumberOfEquations(); i++)
			norm2 += f_[i] * f_[i];
		norm2 = std::sqrt(norm2);

		double yp_mean = 0.;
		for (int i = 0; i < NumberOfEquations(); i++)
			yp_mean += std::fabs(f_[i]);
		yp_mean /= NumberOfEquations();

		std::cout << " n=" << NumberOfEquations() / BlockDimensions() << std::scientific << "  ||f||=" << norm2 << "  |y'|=" << yp_mean << std::endl;
	}

	void OpenSMOKE_CounterFlowFlame1D::InitializeFromBackupFile(const boost::filesystem::path path_file)
	{
		std::vector<double> x_old;
		std::vector<double> T_old;
		std::vector<double> P_old;
		std::vector<double> U_old;
		std::vector<double> G_old;
		std::vector<double> H_old;
		std::vector < std::vector<double> > omega_old;
		ReadFromBackupFile(path_file, thermodynamicsMap_, x_old, T_old, P_old, U_old, G_old, H_old, omega_old);

		// Interpolation (if needed)
		{
		}

		// Grid
		{
			std::cout << " * Building the new mesh from backup data..." << std::endl;

			Eigen::VectorXd x_eigen(x_old.size());
			std::cout << "   Position of fixed point (mm): " << grid_.x_fixed_point()*1e3 << std::endl;
		 	std::cout << "   Index of fixed point:         " << grid_.fixed_point() << "/" << grid_.Np() << std::endl;

			// The grid used is the old one
			for (unsigned int i = 0; i < x_old.size(); i++)
				x_eigen(i) = x_old[i];

			// Update the grid
			grid_.ResetFixedPoint();
			grid_.Update(x_eigen);
		}

		// Setup 
		{
			std::cout << " * Building the first-guess solution from backup data..." << std::endl;

			MemoryAllocation();

			// Pressure
			// The user is free to choose a different pressure (useful to investigate the role of pressure)

			// Temperature
			for (int i = 0; i < grid_.Np(); i++)
				T_(i) = T_old[i];

			// Momentum equation variables
			for (int i = 0; i < grid_.Np(); i++)
			{
				U_(i) = U_old[i];
				G_(i) = G_old[i];
				H_(i) = H_old[i];
			}

			// Composition
			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			{
				for (int j = 0; j < grid_.Np(); j++)
					Y_[j](i) = omega_old[i][j];
			}

			Properties();
			DiffusionFluxes();
		}

		// Messages on the screen
		if (P_ != P_old[0])
		{
			std::cout << "* WARNING: The backup pressure was equal to: " << P_old[0] << " Pa, while the current pressure is: " << P_ << " Pa" << std::endl;
			std::cout << "           Changing the pressure is allowed, but be sure that this is really what you want to do!" << std::endl;
		}

		if (T_fuel_ != T_old[0])
		{
			std::cout << "* WARNING: The backup temperature was equal to: " << T_old[0] << " K, while the current temperature is: " << T_ << " K" << std::endl;
			std::cout << "           Changing the temperature is allowed, but be sure that this is really what you want to do!" << std::endl;
		}

		if (T_oxidizer_ != T_old[grid_.Ni()])
		{
			std::cout << "* WARNING: The backup temperature was equal to: " << T_old[grid_.Ni()] << " K, while the current temperature is: " << T_ << " K" << std::endl;
			std::cout << "           Changing the temperature is allowed, but be sure that this is really what you want to do!" << std::endl;
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		std::vector<unsigned int> rows_single;
		std::vector<unsigned int> cols_single;
		OpenSMOKE::SparsityPatternTridiagonal(NumberOfEquations() / BlockDimensions(), rows_single, cols_single);
		OpenSMOKE::SparsityPatternBlock(NumberOfEquations() / BlockDimensions(), BlockDimensions(), rows_single, cols_single, rows, cols);
	}

	double OpenSMOKE_CounterFlowFlame1D::ReconstructMixtureFraction(const int i)
	{
		std::vector<double> y(thermodynamicsMap_.NumberOfSpecies());
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			y[j] = Y_[i](j);

		std::vector<double> y_fuel;
		std::vector<std::string> names_fuel;
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			if (Y_[0][j] > 0.)
			{
				y_fuel.push_back(Y_[0][j]);
				names_fuel.push_back(thermodynamicsMap_.NamesOfSpecies()[j]);
			}

		std::vector<double> y_oxidizer;
		std::vector<std::string> names_oxidizer;
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			if (Y_[grid_.Ni()][j] > 0.)
			{
				y_oxidizer.push_back(Y_[grid_.Ni()][j]);
				names_oxidizer.push_back(thermodynamicsMap_.NamesOfSpecies()[j]);
			}

		const double csi = thermodynamicsMap_.GetMixtureFractionFromMassFractions(y.data(), names_fuel, y_fuel, names_oxidizer, y_oxidizer);
		return csi;
	}

	void OpenSMOKE_CounterFlowFlame1D::SootDepositionOnTheWall()
	{
		if (is_deposition_wall_ == true)
		{
			if (is_polimi_soot_ == true)
			{
				if (polimi_soot_analyzer_->thermophoretic_effect() == true)
				{
						soot_deposition_ = 0.;	// [kg/m2/s]
						for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
							soot_deposition_ += j_thermophoretic_star_[grid_.Ni()-1](jj);
				}
			}
		}
	}

	void OpenSMOKE_CounterFlowFlame1D::UpdateBoundaries(const double t)
	{
		if (is_dynamic_boundaries_ == true && type_ != Simulation_Type::SIMULATION_TYPE_CSI)
		{
			dynamic_boundaries_->UpdateBoundaryConditions(t,	U_fuel_, U_oxidizer_, T_fuel_, T_oxidizer_,
																rho_fuel_, rho_oxidizer_, 
																Y_fuel_.data(), Y_oxidizer_.data());

			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			{
				boundary_fuel_mass_fractions_(j) = rho_fuel_ * v_fuel_*Y_fuel_(j);
				boundary_oxidizer_mass_fractions_(j) = -rho_oxidizer_ * v_oxidizer_*Y_oxidizer_(j);
			}

			v_fuel_ = U_fuel_ / rho_fuel_ * (n_geometry_ - 1.);
			v_oxidizer_ = -U_oxidizer_/rho_oxidizer_ * (n_geometry_ - 1.);		
		}
	}

	double OpenSMOKE_CounterFlowFlame1D::SutherlandViscosity(const double T)
	{
		return 1.716e-5*std::pow(T / 273.15, 1.5)*(273.15 + 110.4) / (T + 110.4);	// [kg/m/s]
	}
}
