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
|   Copyright(C) 2018 Alberto Cuoci and Mattia Bissoli                    |
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
	IgnitionDelayTimes_Analyzer::IgnitionDelayTimes_Analyzer()
	{
		is_active_ = false;
		is_temperature_ = true;
		is_pressure_ = true;
		is_species_slope_ = true;

		x_threshold_ = 1.e-12;
		time_minimum_ = 1.e-12;
		time_minimum_interval_ = 1.e-12;

		temperature_increase_ = 0.;

		Reset();
	}

	void IgnitionDelayTimes_Analyzer::Reset()
	{
		tOld_ = 0.;
		TOld_ = 0.;
		POld_ = 0.;
		std::fill(xOld_.begin(), xOld_.end(), 0.);

		temperature_slope_max_ = 0.;
		pressure_slope_max_ = 0.;
		std::fill(species_max_.begin(), species_max_.end(), 0.);
		std::fill(species_slope_max_.begin(), species_slope_max_.end(), 0.);

		temperature_slope_tau_ = 0.;
		pressure_slope_tau_ = 0.;
		std::fill(species_max_tau_.begin(), species_max_tau_.end(), 0.);
		std::fill(species_slope_max_tau_.begin(), species_slope_max_tau_.end(), 0.);

		T0_ = 0.;
		temperature_increase_tau_ = 0.;
	}

	template<typename Thermodynamics>
	void IgnitionDelayTimes_Analyzer::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMapXML)
	{
		Grammar_IgnitionDelayTimes grammar_idts;
		dictionary.SetGrammar(grammar_idts);

		// Temperature
		if (dictionary.CheckOption("@Temperature") == true)
		{
			dictionary.ReadBool("@Temperature", is_temperature_);
			is_active_ = true;
		}

		// Pressure
		if (dictionary.CheckOption("@Pressure") == true)
		{
			dictionary.ReadBool("@Pressure", is_pressure_);
			is_active_ = true;
		}

		// Max species
		if (dictionary.CheckOption("@Species") == true)
		{
			dictionary.ReadOption("@Species", species_name_);

			for (unsigned int i = 0; i < species_name_.size(); i++)
				index_species_.push_back(thermodynamicsMapXML.IndexOfSpecies(species_name_[i]) - 1);

			species_max_.resize(species_name_.size());
			species_max_tau_.resize(species_name_.size());
			species_slope_max_.resize(species_name_.size());
			species_slope_max_tau_.resize(species_name_.size());
			xOld_.resize(species_name_.size());

			is_active_ = true;
		}

		// Slope species
		if (dictionary.CheckOption("@SpeciesSlope") == true)
		{
			dictionary.ReadBool("@SpeciesSlope", is_species_slope_);
			is_active_ = true;
		}

		// Max species
		if (dictionary.CheckOption("@SpeciesThreshold") == true)
			dictionary.ReadDouble("@SpeciesThreshold", x_threshold_);

		// Minimum time for calculations
		if (dictionary.CheckOption("@MinimumTime") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@MinimumTime", value, units);
			if (units == "s")		  time_minimum_ = value;
			else if (units == "ms")   time_minimum_ = value / 1000.;
			else if (units == "min")  time_minimum_ = value * 60.;
			else if (units == "h")    time_minimum_ = value * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		// Minimum time interval for derivative calculations
		if (dictionary.CheckOption("@MinimumTimeInterval") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@MinimumTimeInterval", value, units);
			if (units == "s")		  time_minimum_interval_ = value;
			else if (units == "ms")   time_minimum_interval_ = value / 1000.;
			else if (units == "min")  time_minimum_interval_ = value * 60.;
			else if (units == "h")    time_minimum_interval_ = value * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		// Temperature
		if (dictionary.CheckOption("@TemperatureIncrease") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@TemperatureIncrease", value, units);
			if (units == "K")		 temperature_increase_ = value;
			else if (units == "C")   temperature_increase_ = value;
			else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
		}
	}

	void IgnitionDelayTimes_Analyzer::Analyze(const double t, const double T, const double P, const double* x)
	{
		// Maximum criteria
		if ( t>time_minimum_ )
		{
			for (unsigned int i = 0; i < index_species_.size(); i++)
			{
				if ( (x[index_species_[i]] > species_max_[i]) && (x[index_species_[i]]>x_threshold_) )
				{
					species_max_[i] = x[index_species_[i]];
					species_max_tau_[i] = t;
				}
			}
		}

		// Slope criteria
		if ( (t>time_minimum_) && (t - tOld_) > time_minimum_interval_  )
		{
			if (is_temperature_ == true)
			{
				const double temperature_slope_current = (T - TOld_) / (t - tOld_);
				if (temperature_slope_current > temperature_slope_max_)
				{
					temperature_slope_max_ = temperature_slope_current;
					temperature_slope_tau_ = t;
				}
			}

			if (is_pressure_ == true)
			{
				const double pressure_slope_current = (P - POld_) / (t - tOld_);
				if (pressure_slope_current > pressure_slope_max_)
				{
					pressure_slope_max_ = pressure_slope_current;
					pressure_slope_tau_ = t;
				}
			}

			if (temperature_increase_ != 0. && temperature_increase_tau_ == 0.)
			{
				if (T0_ == 0.)
					T0_ = T;

				if (T > (T0_+temperature_increase_))
					temperature_increase_tau_ = tOld_ + temperature_increase_/(T-TOld_)*(t-tOld_);
			}

			// Species
			if (is_species_slope_ == true)
			{
				for (unsigned int i = 0; i < index_species_.size(); i++)
				{
					const double species_slope_current = (x[index_species_[i]] - xOld_[i]) / (t - tOld_);

					if ((species_slope_current > species_slope_max_[i]) && (x[index_species_[i]] > x_threshold_))
					{
						species_slope_max_[i] = species_slope_current;
						species_slope_max_tau_[i] = t;
					}
				}
			}
		}

		// Save previous values
		tOld_ = t;
		TOld_ = T;
		POld_ = P;
		for (unsigned int i = 0; i < index_species_.size(); i++)
			xOld_[i] = x[index_species_[i]];
	}

	void IgnitionDelayTimes_Analyzer::PrintOnFile(const boost::filesystem::path file_name)
	{
		std::ofstream fOutput(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		fOutput << "------------------------------------------------------------" << std::endl;
		fOutput << "Criterion                     Tau[s]         Value          " << std::endl;
		fOutput << "------------------------------------------------------------" << std::endl;

		if (is_temperature_ == true)
		{
			fOutput << std::setw(30) << std::left << "T(slope)";
			fOutput << std::setw(15) << std::left << temperature_slope_tau_;
			fOutput << std::setw(15) << std::left << temperature_slope_max_;
			fOutput << std::endl;
		}

		if (is_pressure_ == true)
		{
			fOutput << std::setw(30) << std::left << "P(slope)";
			fOutput << std::setw(15) << std::left << pressure_slope_tau_;
			fOutput << std::setw(15) << std::left << pressure_slope_max_;
			fOutput << std::endl;
		}

		for (unsigned int i = 0; i < index_species_.size(); i++)
		{
			if (is_species_slope_ == true)
			{
				std::string label= species_name_[i] + "(slope)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << species_slope_max_tau_[i];
				fOutput << std::setw(15) << std::left << species_slope_max_[i];
				fOutput << std::endl;
			}

			{
				std::string label = species_name_[i] + "(max)";
				fOutput << std::setw(30) << std::left << label;
				fOutput << std::setw(15) << std::left << species_max_tau_[i];
				fOutput << std::setw(15) << std::left << species_max_[i];
				fOutput << std::endl;
			}
		}

		if (temperature_increase_ != 0.)
		{
			fOutput << std::setw(30) << std::left << "T(increase)";
			fOutput << std::setw(15) << std::left << temperature_increase_tau_;
			fOutput << std::setw(15) << std::left << T0_+temperature_increase_;
			fOutput << std::endl;
		}

		fOutput.close();
	}

	void IgnitionDelayTimes_Analyzer::PrintHeaderLine(unsigned int& counter, std::ostream& fOutput)
	{
		if (is_temperature_ == true)
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_T(slope)[s]", counter);

		if (is_pressure_ == true)
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_P(slope)[s]", counter);

		for (unsigned int i = 0; i < index_species_.size(); i++)
		{
			if (is_species_slope_ == true)
			{
				std::string label = "tau_" + species_name_[i] + "(slope)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}

			{
				std::string label = "tau_" + species_name_[i] + "(max)[s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutput, label, counter);
			}
		}

		if (temperature_increase_ != 0.)
		{
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "tau_T(delta)[s]", counter);
		}
	}

	void IgnitionDelayTimes_Analyzer::Print(std::ostream& fOutput)
	{
		if (is_temperature_ == true)
			fOutput << std::setw(25) << std::left << temperature_slope_tau_;

		if (is_pressure_ == true)
			fOutput << std::setw(25) << std::left << pressure_slope_tau_;

		for (unsigned int i = 0; i < index_species_.size(); i++)
		{
			if (is_species_slope_ == true)
				fOutput << std::setw(30) << std::left << species_slope_max_tau_[i];
			fOutput << std::setw(30) << std::left << species_max_tau_[i];
		}

		if (temperature_increase_ != 0.)
			fOutput << std::setw(25) << std::left << temperature_increase_tau_;
	}
}


