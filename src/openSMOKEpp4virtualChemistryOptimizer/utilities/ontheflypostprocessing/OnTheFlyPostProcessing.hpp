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

#include "Grammar_OnTheFlyPostProcessing.h"

namespace OpenSMOKE
{

	OnTheFlyPostProcessing::OnTheFlyPostProcessing(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
													OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
													const boost::filesystem::path path_output) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	{
		path_output_ = path_output;
		is_active_ = false;
		print_formation_rates_ = false;
		print_reaction_rates_ = false;
		formation_rates_moles_ = true;
	}

	void OnTheFlyPostProcessing::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_OnTheFlyPostProcessing grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@FormationRates") == true)
		{
			std::vector<std::string> formation_rates_species;
			dictionary.ReadOption("@FormationRates", formation_rates_species);

			if (formation_rates_species[0] == "ALL" || formation_rates_species[0] == "all")
			{
				print_formation_rates_ = true;
			}
			else if (formation_rates_species[0] == "NONE" || formation_rates_species[0] == "none")
			{
				print_formation_rates_ = false;
			}
			else
			{
				print_formation_rates_ = true;
				indices_of_formation_rates_species_.resize(formation_rates_species.size());
				for (unsigned int i = 0; i < formation_rates_species.size(); i++)
					indices_of_formation_rates_species_[i] = thermodynamicsMap_.IndexOfSpecies(formation_rates_species[i]);
			}
		}

		if (dictionary.CheckOption("@ReactionRates") == true)
		{
			std::vector<std::string> output_reaction_rates;
			dictionary.ReadOption("@ReactionRates", output_reaction_rates);

			if (output_reaction_rates[0] == "ALL" || output_reaction_rates[0] == "all")
			{
				print_reaction_rates_ = true;
			}
			else if (output_reaction_rates[0] == "NONE" || output_reaction_rates[0] == "none")
			{
				print_reaction_rates_ = false;
			}
			else
			{
				print_reaction_rates_ = true;
				indices_of_reaction_rates_.resize(output_reaction_rates.size());
				for (unsigned int i = 0; i < output_reaction_rates.size(); i++)
					indices_of_reaction_rates_[i] = boost::lexical_cast<unsigned int>(output_reaction_rates[i]);
			}
		}

		if (dictionary.CheckOption("@FormationRatesUnits") == true)
		{
			std::string formation_rates_units;
			dictionary.ReadString("@FormationRatesUnits", formation_rates_units);

			if (formation_rates_units == "mass")
				formation_rates_moles_ = false;
			else if (formation_rates_units == "moles")
				formation_rates_moles_ = true;
			else
				OpenSMOKE::FatalErrorMessage("@FormationRatesUnits: wrong definition. Allowed options are: mass | moles (default)");
		}

		if (print_formation_rates_ == true || print_reaction_rates_ == true)
			is_active_ = true;

		MemoryAllocation();
	}

	void OnTheFlyPostProcessing::MemoryAllocation()
	{
		// Basic variables
		NS_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();

		if (is_active_ == true)
		{
			OpenSMOKE::ChangeDimensions(NS_, &x_, true);
			OpenSMOKE::ChangeDimensions(NS_, &c_, true);
			OpenSMOKE::ChangeDimensions(NS_, &R_, true);
		}
	}

	void OnTheFlyPostProcessing::CloseOutputFiles()
	{
		fFormationRates_.close();
		fReactionRates_.close();
	}

	void OnTheFlyPostProcessing::PrepareOutputFiles(const std::vector<std::string>& additional)
	{
		// Output file
		if (print_formation_rates_ == true)
		{
			if (indices_of_formation_rates_species_.size() != 0)
			{
				widths_of_formation_rates_species_.resize(indices_of_formation_rates_species_.size());
				for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
					widths_of_formation_rates_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[indices_of_formation_rates_species_[i] - 1], NS_);
			}
			else
			{
				widths_of_formation_rates_species_.resize(NS_);
				for (unsigned int i = 0; i < NS_; i++)
					widths_of_formation_rates_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[i], NS_);
			}

			const boost::filesystem::path file_name = path_output_ / "FormationRates.out";
			fFormationRates_.open(file_name.c_str(), std::ios::out);
			fFormationRates_.setf(std::ios::scientific);

			unsigned int counter = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "t[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "x[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "y[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "z[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "P[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, "QR[W/m3]", counter);

			// Additional variables
			for (unsigned int i = 0; i<additional.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fFormationRates_, additional[i], counter);

			// Species
			if (indices_of_formation_rates_species_.size() != 0)
			{
				for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
					OpenSMOKE::PrintTagOnASCIILabel(widths_of_formation_rates_species_[i], fFormationRates_, thermodynamicsMap_.NamesOfSpecies()[indices_of_formation_rates_species_[i] - 1], counter);
			}
			else
			{
				for (unsigned int i = 0; i < NS_; i++)
					OpenSMOKE::PrintTagOnASCIILabel(widths_of_formation_rates_species_[i], fFormationRates_, thermodynamicsMap_.NamesOfSpecies()[i], counter);
			}

			fFormationRates_ << std::endl;
		}

		if (print_reaction_rates_ == true)
		{
			const boost::filesystem::path file_name = path_output_ / "ReactionRates.out";
			fReactionRates_.open(file_name.c_str(), std::ios::out);
			fReactionRates_.setf(std::ios::scientific);

			unsigned int counter = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "t[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "x[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "y[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "z[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "P[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "QR[W/m3]", counter);

			// Additional variables
			for (unsigned int i = 0; i<additional.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, additional[i], counter);

			// Reactions
			if (indices_of_reaction_rates_.size() != 0)
			{
				for (unsigned int i = 0; i < indices_of_reaction_rates_.size(); i++)
					OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "r" + boost::lexical_cast<std::string>(indices_of_reaction_rates_[i]), counter);
			}
			else
			{
				for (unsigned int i = 0; i < NR_; i++)
					OpenSMOKE::PrintTagOnASCIILabel(20, fReactionRates_, "r" + boost::lexical_cast<std::string>(i + 1), counter);
			}

			fReactionRates_ << std::endl;
		}
	}

	void OnTheFlyPostProcessing::PrepareOutputFiles()
	{
		std::vector<std::string> additional(0);
		PrepareOutputFiles(additional);
	}

	void OnTheFlyPostProcessing::WriteOnFile(const double t, const double x, const double y, const double z, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omega, const std::vector<double>& additional)
	{
		double QR = 0;

		// Kinetic data
		{
			thermodynamicsMap_.SetTemperature(T);
			thermodynamicsMap_.SetPressure(P_Pa);

			double MW;
			thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW, omega.GetHandle());
			const double cTot = P_Pa / (PhysicalConstants::R_J_kmol * T);
			Product(cTot, x_, &c_);
			
			kineticsMap_.SetTemperature(T);
			kineticsMap_.SetPressure(P_Pa);
			kineticsMap_.ReactionRates(c_.GetHandle());
			kineticsMap_.FormationRates(R_.GetHandle());
			QR = kineticsMap_.HeatRelease(R_.GetHandle());
		}

		// Print formation rates
		if (print_formation_rates_ == true)
		{
			// Basic variables
			fFormationRates_ << std::setw(20) << std::left << t;
			fFormationRates_ << std::setw(20) << std::left << x;
			fFormationRates_ << std::setw(20) << std::left << y;
			fFormationRates_ << std::setw(20) << std::left << z;
			fFormationRates_ << std::setw(20) << std::left << T;
			fFormationRates_ << std::setw(20) << std::left << P_Pa;

			// Reaction heat [W/m3]
			fFormationRates_ << std::setw(20) << std::left << QR;

			// Additional variables
			for (unsigned int i=0;i<additional.size();i++)
				fFormationRates_ << std::setw(20) << std::left << additional[i];

			// Formation Rates [kmol/m3/s]
			if (formation_rates_moles_ == true)
			{
				if (indices_of_formation_rates_species_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
						fFormationRates_ << std::setw(widths_of_formation_rates_species_[i]) << std::left << R_[indices_of_formation_rates_species_[i]];
				}
				else
				{
					for (unsigned int i = 0; i < NS_; i++)
						fFormationRates_ << std::setw(widths_of_formation_rates_species_[i]) << std::left << R_[i+1];
				}
			}
			else
			{
				if (indices_of_formation_rates_species_.size() != 0)
				{
					for (unsigned int i = 0; i < indices_of_formation_rates_species_.size(); i++)
						fFormationRates_ << std::setw(widths_of_formation_rates_species_[i]) << std::left << thermodynamicsMap_.MW(indices_of_formation_rates_species_[i]-1) * R_[indices_of_formation_rates_species_[i]];
				}
				else
				{
					for (unsigned int i = 0; i < NS_; i++)
						fFormationRates_ << std::setw(widths_of_formation_rates_species_[i]) << std::left << thermodynamicsMap_.MW(i) * R_[i+1];
				}
			}

			fFormationRates_ << std::endl;
		}

		// Print reaction rates
		if (print_reaction_rates_ == true)
		{
			// Reaction rates
			OpenSMOKE::OpenSMOKEVectorDouble r(NR_);
			kineticsMap_.GiveMeReactionRates(r.GetHandle());

			// Basic variables
			fReactionRates_ << std::setw(20) << std::left << t;
			fReactionRates_ << std::setw(20) << std::left << x;
			fReactionRates_ << std::setw(20) << std::left << y;
			fReactionRates_ << std::setw(20) << std::left << z;
			fReactionRates_ << std::setw(20) << std::left << T;
			fReactionRates_ << std::setw(20) << std::left << P_Pa;

			// Reaction heat [W/m3]
			fReactionRates_ << std::setw(20) << std::left << QR;

			// Additional variables
			for (unsigned int i = 0; i<additional.size(); i++)
				fReactionRates_ << std::setw(20) << std::left << additional[i];

			if (indices_of_reaction_rates_.size() != 0)
			{
				for (unsigned int i = 0; i < indices_of_reaction_rates_.size(); i++)
					fReactionRates_ << std::setw(20) << std::left << r[indices_of_reaction_rates_[i]];
			}
			else
			{
				for (unsigned int i = 1; i <= NR_; i++)
					fReactionRates_ << std::setw(20) << std::left << r[i];
			}

			fReactionRates_ << std::endl;
		}
	}

	void OnTheFlyPostProcessing::WriteOnFile(const double t, const double x, const double y, const double z, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omega)
	{
		std::vector<double> additional(0);
		WriteOnFile(t, x, y, z, T, P_Pa, omega, additional);
	}
}

