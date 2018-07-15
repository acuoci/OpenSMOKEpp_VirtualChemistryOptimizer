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

#include "Grammar_DynamicBoundaries.h"

namespace OpenSMOKE
{
	DynamicBoundaries::DynamicBoundaries(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
											const boost::filesystem::path path_output) :
	thermodynamicsMap_(thermodynamicsMap)
	{
		path_output_ = path_output;
		is_active_ = false;
		end_time_ = 1e4;
		compact_output_ = false;
		n_steps_output_ = 5;

		type_ = DynamicBoundaries::NONE;
		
		frequency_fuel_ = 0.;
		frequency_ox_ = 0.;
		semi_amplitude_fuel_ = 0.;
		semi_amplitude_ox_ = 0.;
		t_delay_fuel_ = 0.;
		t_delay_ox_ = 0.;

		slope_fuel_velocity_ = 0.;
		slope_ox_velocity_ = 0.;
		slope_fuel_temperature_ = 0.;
		slope_ox_temperature_ = 0.;

		is_virtual_chemistry_ = false;
	}


	void DynamicBoundaries::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_DynamicBoundaries grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Type") == true)
		{
			std::string type;
			dictionary.ReadString("@Type", type);

			if (type == "TEMPERATURE_VELOCITY_SLOPES")
				type_ = DynamicBoundaries::TEMPERATURE_VELOCITY_SLOPES;
			else if (type == "FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE")
				type_ = DynamicBoundaries::FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE;
			else if (type == "FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE")
				type_ = DynamicBoundaries::FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE;
			else if (type == "FIXED_BALANCE_FUEL_VELOCITY_SLOPE")
				type_ = DynamicBoundaries::FIXED_BALANCE_FUEL_VELOCITY_SLOPE;
			else if (type == "FIXED_BALANCE_OX_VELOCITY_SLOPE")
				type_ = DynamicBoundaries::FIXED_BALANCE_OX_VELOCITY_SLOPE;
			else if (type == "SIN_VELOCITY_INPHASE")
				type_ = DynamicBoundaries::SIN_VELOCITY_INPHASE;
			
			else
			{
				const std::string message1 = "Unknown @Type: " + type;
				const std::string message2 = "Available types: TEMPERATURE_VELOCITY_SLOPES | FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE | FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE | FIXED_BALANCE_FUEL_VELOCITY_SLOPE | FIXED_BALANCE_OX_VELOCITY_SLOPE | SIN_VELOCITY_INPHASE";
				FatalErrorMessage(message1 + ". " + message2);
			}
		}

		if (dictionary.CheckOption("@SemiAmplitudeFuel") == true)
		{
			dictionary.ReadDouble("@SemiAmplitudeFuel", semi_amplitude_fuel_);
		}

		if (dictionary.CheckOption("@SemiAmplitudeOxidizer") == true)
		{
			dictionary.ReadDouble("@SemiAmplitudeOxidizer", semi_amplitude_ox_);
		}

		if (dictionary.CheckOption("@FrequencyFuel") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@FrequencyFuel", value, units);
			if (units == "1/s")			frequency_fuel_ = value;
			else if (units == "Hz")		frequency_fuel_ = value;
			else OpenSMOKE::FatalErrorMessage("Unknown frequency units");
		}

		if (dictionary.CheckOption("@FrequencyOxidizer") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@FrequencyOxidizer", value, units);
			if (units == "1/s")			frequency_ox_ = value;
			else if (units == "Hz")		frequency_ox_ = value;
			else OpenSMOKE::FatalErrorMessage("Unknown frequency units");
		}

		if (dictionary.CheckOption("@DelayTimeFuel") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@DelayTimeFuel", value, units);
			if (units == "s")			t_delay_fuel_ = value;
			else if (units == "ms")		t_delay_fuel_ = value/1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		if (dictionary.CheckOption("@DelayTimeOxidizer") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@DelayTimeOxidizer", value, units);
			if (units == "s")			t_delay_ox_ = value;
			else if (units == "ms")		t_delay_ox_ = value / 1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}


		if (dictionary.CheckOption("@SlopeFuelVelocity") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@SlopeFuelVelocity", value, units);
			if (units == "m/s2")			slope_fuel_velocity_ = value;
			else if (units == "cm/s2")		slope_fuel_velocity_ = value / 100.;
			else if (units == "mm/s2")		slope_fuel_velocity_ = value / 1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown slope velocity units");
		}

		if (dictionary.CheckOption("@SlopeFuelTemperature") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@SlopeFuelTemperature", value, units);
			if (units == "K/s")				slope_fuel_temperature_ = value;
			else if (units == "K/min")		slope_fuel_temperature_ = value / 60.;
			else if (units == "C/s")		slope_fuel_temperature_ = value;
			else if (units == "C/min")		slope_fuel_temperature_ = value / 60.;
			else OpenSMOKE::FatalErrorMessage("Unknown slope temperature units");
		}

		if (dictionary.CheckOption("@SlopeOxidizerVelocity") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@SlopeOxidizerVelocity", value, units);
			if (units == "m/s2")		slope_ox_velocity_ = value;
			else if (units == "cm/s2")  slope_ox_velocity_ = value / 100.;
			else if (units == "mm/s2")	slope_ox_velocity_ = value / 1000.;
			else OpenSMOKE::FatalErrorMessage("Unknown slope velocity units");
		}

		if (dictionary.CheckOption("@SlopeOxidizerTemperature") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@SlopeOxidizerTemperature", value, units);
			if (units == "K/s")				slope_ox_temperature_ = value;
			else if (units == "K/min")		slope_ox_temperature_ = value / 60.;
			else if (units == "C/s")		slope_ox_temperature_ = value;
			else if (units == "C/min")		slope_ox_temperature_ = value / 60.;
			else OpenSMOKE::FatalErrorMessage("Unknown slope temperature units");
		}

		if (dictionary.CheckOption("@EndTime") == true)
		{
			double value;
			std::string units;
			dictionary.ReadMeasure("@EndTime", value, units);
			if (units == "s")			end_time_ = value;
			else if (units == "min")	end_time_ = value*60.;
			else if (units == "hr")		end_time_ = value*3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		if (dictionary.CheckOption("@CompactOutput") == true)
		{
			dictionary.ReadBool("@CompactOutput", compact_output_);
		}

		if (dictionary.CheckOption("@StepsOutput") == true)
		{
			dictionary.ReadInt("@StepsOutput", n_steps_output_);
		}

		// Snapshot times
		{
			std::vector<std::string> dummy;
			if (dictionary.CheckOption("@SnapshotTimes") == true)
			{
				dictionary.ReadOption("@SnapshotTimes", dummy);
				const std::string units = dummy[dummy.size() - 1];
				double conversion = 1.;
				if (units == "s")			conversion = 1.;
				else if (units == "min")	conversion = 60.;
				else if (units == "hr")		conversion = 3600.;
				snapshot_list_of_times_.resize(dummy.size() - 1);
				for (unsigned int j = 0; j < dummy.size() - 1; j++)
					snapshot_list_of_times_[j] = boost::lexical_cast<double>(dummy[j])*conversion;
			}
		}
		
		snapshot_list_of_times_.insert(snapshot_list_of_times_.begin(), 0.);
		snapshot_list_of_times_.push_back(end_time_);

		CheckInput();

		is_active_ = true;
	}

	void DynamicBoundaries::CheckInput()
	{
		if (type_ == DynamicBoundaries::TEMPERATURE_VELOCITY_SLOPES)
		{
			if (slope_fuel_velocity_ == 0. && slope_fuel_temperature_ == 0. &&
				slope_ox_velocity_ == 0. && slope_ox_temperature_ == 0.)

				OpenSMOKE::FatalErrorMessage("In case of TEMPERATURE_VELOCITY_SLOPES option, at least 1 slope must be specified.");
		}

		if (type_ == DynamicBoundaries::FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE)
		{
			if (slope_fuel_velocity_ != 0. || slope_ox_velocity_ != 0. || slope_ox_temperature_ != 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE option, \
												only the slope of fuel temperature can be specified.");

			if (slope_fuel_temperature_ == 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE option, \
												the slope of fuel temperature cannot be equal to zero");
		}

		// Autoignition experiments from Prof. Seshadri (University of San Diego, California)
		if (type_ == DynamicBoundaries::FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE)
		{
			if (slope_fuel_velocity_ != 0. || slope_ox_velocity_ != 0. || slope_fuel_temperature_ != 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE option, \
												only the slope of oxidizer temperature can be specified.");

			if (slope_ox_temperature_ == 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE option, \
												the slope of oxidizer temperature cannot be equal to zero");
		}

		if (type_ == DynamicBoundaries::FIXED_BALANCE_FUEL_VELOCITY_SLOPE)
		{
			if (slope_ox_temperature_ != 0. || slope_ox_velocity_ != 0. || slope_fuel_temperature_ != 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_BALANCE_FUEL_VELOCITY_SLOPE option, \
												only the slope of fuel velocity can be specified.");

			if (slope_fuel_velocity_ == 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_BALANCE_FUEL_VELOCITY_SLOPE option, \
												the slope of fuel velocity cannot be equal to zero");
		}

		if (type_ == DynamicBoundaries::FIXED_BALANCE_OX_VELOCITY_SLOPE)
		{
			if (slope_ox_temperature_ != 0. || slope_fuel_velocity_ != 0. || slope_fuel_temperature_ != 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_BALANCE_OX_VELOCITY_SLOPE option, \
												only the slope of oxidizer velocity can be specified.");

			if (slope_ox_velocity_ == 0.)

				OpenSMOKE::FatalErrorMessage("	In case of FIXED_BALANCE_OX_VELOCITY_SLOPE option, \
												the slope of oxidizer velocity cannot be equal to zero");
		}

		if (type_ == DynamicBoundaries::SIN_VELOCITY_INPHASE)
		{
			if (slope_ox_temperature_ != 0. || slope_fuel_velocity_ != 0. || slope_fuel_temperature_ != 0. || slope_ox_velocity_ != 0.)
				OpenSMOKE::FatalErrorMessage("	In case of SIN_VELOCITY_INPHASE, \
												only the frequency and semi-amplitude os velocities can be specified.");

			if (semi_amplitude_fuel_ == 0. && semi_amplitude_ox_ == 0.)
				OpenSMOKE::FatalErrorMessage("	In case of SIN_VELOCITY_INPHASE option, \
												at least one of the two semi-amplitudes must be not equal to 0");

			if (frequency_fuel_ == 0. && frequency_ox_ == 0.)
				OpenSMOKE::FatalErrorMessage("	In case of SIN_VELOCITY_INPHASE option, \
												at least one of the two fequencies must be not equal to 0");
		}
	}

	void DynamicBoundaries::SetVirtualChemistry(OpenSMOKE::VirtualChemistry* virtual_chemistry)
	{
		is_virtual_chemistry_ = true;
		virtual_chemistry_ = virtual_chemistry;
	}

	void DynamicBoundaries::SetFlame(const double length, const double P_Pa, const double n_geometry)
	{
		length_ = length;
		P_Pa_ = P_Pa;
		n_geometry_ = n_geometry;
	}

	void DynamicBoundaries::SetInitialFuelSide(const double T, const double v, const double* omega)
	{
		omega_fuel_0_.resize(thermodynamicsMap_.NumberOfSpecies());
		x_fuel_0_.resize(thermodynamicsMap_.NumberOfSpecies());

		T_fuel_0_ = T;
		v_fuel_0_ = v;
		for (unsigned int i=0;i<thermodynamicsMap_.NumberOfSpecies();i++)
			omega_fuel_0_(i) = omega[i];

		thermodynamicsMap_.MoleFractions_From_MassFractions(x_fuel_0_.data(), mw_fuel_0_, omega_fuel_0_.data());

		if (is_virtual_chemistry_ == true)
		{
			mw_fuel_0_ = virtual_chemistry_->MWMix(omega_fuel_0_.data());
			virtual_chemistry_->MoleFractions(mw_fuel_0_, omega_fuel_0_.data(), x_fuel_0_.data());
		}

		rho_fuel_0_ = P_Pa_ * mw_fuel_0_ / PhysicalConstants::R_J_kmol / T_fuel_0_;
		U_fuel_0_ = rho_fuel_0_ * v_fuel_0_ / (n_geometry_ - 1.);

		T_fuel_ = T_fuel_0_;
		mw_fuel_ = mw_fuel_0_;
		rho_fuel_ = rho_fuel_0_;
		omega_fuel_ = omega_fuel_0_;
		x_fuel_ = x_fuel_0_;
		U_fuel_ = U_fuel_0_;
		v_fuel_ = v_fuel_0_;
	}

	void DynamicBoundaries::SetInitialOxidizerSide(const double T, const double v, const double* omega)
	{
		omega_ox_0_.resize(thermodynamicsMap_.NumberOfSpecies());
		x_ox_0_.resize(thermodynamicsMap_.NumberOfSpecies());

		T_ox_0_ = T;
		v_ox_0_ = v;
		for (unsigned int i = 0; i<thermodynamicsMap_.NumberOfSpecies(); i++)
			omega_ox_0_(i) = omega[i];

		thermodynamicsMap_.MoleFractions_From_MassFractions(x_ox_0_.data(), mw_ox_0_, omega_ox_0_.data());

		if (is_virtual_chemistry_ == true)
		{
			mw_ox_0_ = virtual_chemistry_->MWMix(omega_ox_0_.data());
			virtual_chemistry_->MoleFractions(mw_ox_0_, omega_ox_0_.data(), x_ox_0_.data());
		}

		rho_ox_0_ = P_Pa_ * mw_ox_0_ / PhysicalConstants::R_J_kmol / T_ox_0_;
		U_ox_0_ = -rho_ox_0_ * v_ox_0_ / (n_geometry_ - 1.);

		T_ox_ = T_ox_0_;
		mw_ox_ = mw_ox_0_;
		rho_ox_ = rho_ox_0_;
		omega_ox_ = omega_ox_0_;
		x_ox_ = x_ox_0_;
		U_ox_ = U_ox_0_;
		v_ox_ = v_ox_0_;
	}

	void DynamicBoundaries::CompleteSetup()
	{
		a_0_ = 2.*v_ox_0_ / length_ * (1. + v_fuel_0_ / v_ox_0_ * std::sqrt(rho_fuel_0_/rho_ox_0_));
		a_seshadri_0_ = 2.*v_fuel_0_ / length_ * (1. - v_ox_0_ / v_fuel_0_ * std::sqrt(rho_ox_0_ / rho_fuel_0_));
		xst_0_ = length_ / (1. + rho_ox_0_ / rho_fuel_0_ * std::pow(v_ox_0_ / v_fuel_0_, 2.));
		beta_0_ = v_fuel_0_ / v_ox_0_ * std::sqrt(rho_fuel_0_ / rho_ox_0_);

		a_ = a_0_;
		a_seshadri_ = a_seshadri_0_;
		xst_ = xst_0_;
		beta_ = beta_0_;
	}

	void DynamicBoundaries::PrepareOutputFiles()
	{
		fOutputDynamics_.open((path_output_ / "DynamicSolution.out").c_str(), std::ios::out);
		fOutputDynamics_.setf(std::ios::scientific);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "t[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "a[1/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "a_seshadri[1/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "xst[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "beta[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "Tmax[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "TF[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "TO[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "vF[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "vO[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "rhoF[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputDynamics_, "rhoO[kg/m3]", count);
			fOutputDynamics_ << std::endl;
		}
	}

	void DynamicBoundaries::Print(const double t, const double Tmax)
	{
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << t;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << a_;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << a_seshadri_;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << xst_ * 100.;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << beta_;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << Tmax;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << T_fuel_;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << T_ox_;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << v_fuel_*100.;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << v_ox_*100.;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << rho_fuel_;
		fOutputDynamics_ << std::scientific << std::setprecision(6) << std::setw(20) << rho_ox_;
		fOutputDynamics_ << std::endl;
	}

	void DynamicBoundaries::Summary(std::ostream& out)
	{
		out << "---------------------------------------------------------------" << std::endl;
		out << " Dynamic boundary conditions                                   " << std::endl;
		out << "---------------------------------------------------------------" << std::endl;
		out << " * Fuel side                                                   " << std::endl;
		out << "   temperature (K):            " << T_fuel_0_ << std::endl;
		out << "   velocity (cm/s):            " << v_fuel_0_*100. << std::endl;
		out << "   molecular weight (kg/kmol): " << mw_fuel_0_ << std::endl;
		out << "   density (kg/m3):            " << rho_fuel_0_ << std::endl;
		out << "   mass flow rate (kg/m2/s):   " << U_fuel_0_ << std::endl;
		out << " * Oxidizer side                                               " << std::endl;
		out << "   temperature (K):            " << T_ox_0_ << std::endl;
		out << "   velocity (cm/s):            " << v_ox_0_*100. << std::endl;
		out << "   molecular weight (kg/kmol): " << mw_ox_0_ << std::endl;
		out << "   density (kg/m3):            " << rho_ox_0_ << std::endl;
		out << "   mass flow rate (kg/m2/s):   " << -U_ox_0_ << std::endl;
		out << " * Flame                                                       " << std::endl;
		out << "   strain rate (1/s):          " << a_0_ << std::endl;
		out << "   stagnation plane (cm):      " << xst_0_ * 100. << std::endl;
		out << "   momentum umbalance:         " << beta_0_ << std::endl;
		out << "---------------------------------------------------------------" << std::endl;
		out << std::endl;
	}

	void DynamicBoundaries::UpdateBoundaryConditions(	const double t,
														double& U_fuel, double& U_ox,
														double& T_fuel, double& T_ox,
														double& rho_fuel, double& rho_ox,
														double* omega_fuel, double* omega_ox) 
	{
		// Description: 1. the fuel temperature, the fuel velocity, the oxidizer temperature and the 
		//                 oxidizer velocity are increased/decreased linearly
		//              2. the strain rate changes
		//              3. the momentum balance changes
		if (type_ == DynamicBoundaries::TEMPERATURE_VELOCITY_SLOPES)
		{
			v_fuel_ = v_fuel_0_ + slope_fuel_velocity_ * t;
			T_fuel_ = T_fuel_0_ + slope_fuel_temperature_ * t;
			rho_fuel_ = P_Pa_ * mw_fuel_ / PhysicalConstants::R_J_kmol / T_fuel_;
			U_fuel_ = rho_fuel_ * v_fuel_ / (n_geometry_ - 1.);

			v_ox_ = v_ox_0_ + slope_ox_velocity_ * t;
			T_ox_ = T_ox_0_ + slope_ox_temperature_ * t;
			rho_ox_ = P_Pa_ * mw_ox_ / PhysicalConstants::R_J_kmol / T_ox_;
			U_ox_ = -rho_ox_ * v_ox_ / (n_geometry_ - 1.);
		}
		
		// Description: 1. the fuel temperature is increased/decreased linearly
		//              2. in order to keep the strain rate constant, the fuel velocity is changed accordingly
		//              3. the oxidizer velocity is kept constant
		//              4. the fuel temperature is kept constant
		//              5. the momentum balance and the stagnation plane are kept constant
		if (type_ == DynamicBoundaries::FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE)
		{
			T_fuel_ = T_fuel_0_ + slope_fuel_temperature_ * t;
			rho_fuel_ = P_Pa_ * mw_fuel_ / PhysicalConstants::R_J_kmol / T_fuel_;
			v_fuel_ = beta_0_*v_ox_ * std::sqrt(rho_ox_ / rho_fuel_);
			U_fuel_ = rho_fuel_ * v_fuel_ / (n_geometry_ - 1.);
		}

		// Autoignition experiments from Prof. Seshadri (University of San Diego, California)
		// Description: 1. the oxidizer temperature is increased/decreased linearly
		//              2. in order to keep the strain rate constant, the oxidizer velocity is changed accordingly
		//              3. in order to keep the strain rate constant, the fuel velocity is changed too
		//              4. the fuel temperature is kept constant
		//              5. the momentum balance and the stagnation plane are kept constant
		if (type_ == DynamicBoundaries::FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE)
		{
			T_ox_ = T_ox_0_ + slope_ox_temperature_ * t;
			rho_ox_ = P_Pa_ * mw_ox_ / PhysicalConstants::R_J_kmol / T_ox_;
			v_ox_ = a_0_ * length_ / 2. / (1.+beta_0_);
			U_ox_ = -rho_ox_ * v_ox_ / (n_geometry_ - 1.);

			v_fuel_ = beta_0_*v_ox_*std::sqrt(rho_ox_/rho_fuel_);
			U_fuel_ = rho_fuel_ * v_fuel_ / (n_geometry_ - 1.);
		}

		// Description: 1. the fuel temperature is increased/decreased linearly
		//              2. the momentum balance is kept constant, which means thath the oxidizer velocity is changed
		//              3. the temperatures are kept constant
		//              4. the strain rate changes
		if (type_ == DynamicBoundaries::FIXED_BALANCE_FUEL_VELOCITY_SLOPE)
		{
			v_fuel_ = v_fuel_0_ + slope_fuel_velocity_ * t;
			U_fuel_ = rho_fuel_ * v_fuel_ / (n_geometry_ - 1.);

			v_ox_ = v_fuel_*std::sqrt(rho_fuel_/rho_ox_)/beta_0_;
			U_ox_ = -rho_ox_ * v_ox_ / (n_geometry_ - 1.);
		}

		// Description: 1. the oxidizer temperature is increased/decreased linearly
		//              2. the momentum balance is kept constant, which means thath the fuel velocity is changed
		//              3. the temperatures are kept constant
		//              4. the strain rate changes
		if (type_ == DynamicBoundaries::FIXED_BALANCE_OX_VELOCITY_SLOPE)
		{
			v_ox_ = v_ox_0_ + slope_ox_velocity_ * t;
			U_ox_ = -rho_ox_ * v_ox_ / (n_geometry_ - 1.);

			v_fuel_ = v_ox_ * std::sqrt(rho_ox_ / rho_fuel_) * beta_0_;
			U_fuel_ = rho_fuel_ * v_fuel_ / (n_geometry_ - 1.);
		}

		// Description: 1. The fuel and oxidizer velocities changes according to a sinusoidal function
		//                 v(t) = v0*(1.+A*sin(2*pi*f))
		if (type_ == DynamicBoundaries::SIN_VELOCITY_INPHASE)
		{
			const double tau_fuel_ = 1. / frequency_fuel_;
			const double alpha_fuel = -log(1.e-5) / tau_fuel_;
			const double phi_fuel = (1. - exp(-alpha_fuel * t));
			//double phi_fuel = 1.;

			v_fuel_ = v_fuel_0_*(1.+phi_fuel*semi_amplitude_fuel_*std::sin(2.*PhysicalConstants::pi*frequency_fuel_*(t+t_delay_fuel_)));
			U_fuel_ = rho_fuel_ * v_fuel_ / (n_geometry_ - 1.);

			const double tau_ox_ = 1. / frequency_ox_;
			double alpha_ox = -log(1.e-5) / tau_ox_;
			double phi_ox = (1. - exp(-alpha_ox * t));
			//double phi_ox = 1.;

			v_ox_ = v_ox_0_ * (1. + phi_ox*semi_amplitude_ox_ * std::sin(2.*PhysicalConstants::pi*frequency_ox_*(t+t_delay_ox_)));
			U_ox_ = -rho_ox_ * v_ox_ / (n_geometry_ - 1.);
		}

		// Update flame features
		a_ = 2.*v_ox_ / length_ * (1. + v_fuel_ / v_ox_ * std::sqrt(rho_fuel_ / rho_ox_));
		a_seshadri_ = 2.*v_fuel_ / length_ * (1. - v_ox_ / v_fuel_ * std::sqrt(rho_ox_ / rho_fuel_));
		xst_ = length_ / (1. + rho_ox_ / rho_fuel_ * std::pow(v_ox_ / v_fuel_, 2.));
		beta_ = v_fuel_ / v_ox_ * std::sqrt(rho_fuel_ / rho_ox_);


		// Transfer data
		U_fuel = U_fuel_;
		U_ox = U_ox_;
		T_fuel = T_fuel_;
		T_ox = T_ox_;
		rho_fuel = rho_fuel_;
		rho_ox = rho_ox_;
		for (unsigned int i = 0; i<thermodynamicsMap_.NumberOfSpecies(); i++)
			omega_fuel[i] = omega_fuel_(i);
		for (unsigned int i = 0; i<thermodynamicsMap_.NumberOfSpecies(); i++)
			omega_ox[i] = omega_ox_(i);
	}
}
