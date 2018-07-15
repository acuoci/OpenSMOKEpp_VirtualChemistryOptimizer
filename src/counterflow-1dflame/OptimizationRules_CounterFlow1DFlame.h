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

#ifndef OpenSMOKE_OptimizationRules_CounterFlow1DFlame_H
#define OpenSMOKE_OptimizationRules_CounterFlow1DFlame_H

namespace OpenSMOKE
{
	class Grammar_OptimizationRules_CounterFlow1DFlame : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature",
				OpenSMOKE::SINGLE_BOOL,
				"True if temperature has to be optimized (default: true)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Species",
				OpenSMOKE::SINGLE_STRING,
				"Name of the species to be optimized (default: none)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the species profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the temperature profile",
				false));
		}
	};

	class OptimizationRules_CounterFlow1DFlame
	{
	public:

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, OpenSMOKE::OpenSMOKE_DictionaryManager& dictionaries);

		unsigned int species_np() const { return species_np_; }
		double species_x(const unsigned int i) const { return species_x_(i); }
		double species_y(const unsigned int i) const { return species_y_(i); }
		const Eigen::VectorXd& species_x() const { return species_x_; }
		const Eigen::VectorXd& species_y() const { return species_y_; }
		double species_mean_y() const { return species_mean_y_; }
		double species_max_y() const { return species_y_.maxCoeff(); }

		unsigned int t_np() const { return t_np_; }
		double t_x(const unsigned int i) const { return t_x_(i); }
		double t_y(const unsigned int i) const { return t_y_(i); }
		const Eigen::VectorXd& t_x() const { return t_x_; }
		const Eigen::VectorXd& t_y() const { return t_y_; }
		double t_mean_y() const { return t_mean_y_; }
		double t_max_y() const { return t_y_.maxCoeff(); }

		bool is_temperature() const { return is_temperature_; }
		bool is_species() const { return is_species_; }
		std::string species_name() const { return species_name_; }

	private:

		unsigned int species_np_;
		Eigen::VectorXd species_x_;
		Eigen::VectorXd species_y_;
		double species_mean_y_;

		unsigned int t_np_;
		Eigen::VectorXd t_x_;
		Eigen::VectorXd t_y_;
		double t_mean_y_;

		bool is_temperature_;
		bool is_species_;
		std::string species_name_;
	};

	void OptimizationRules_CounterFlow1DFlame::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, OpenSMOKE::OpenSMOKE_DictionaryManager& dictionaries)
	{
		Grammar_OptimizationRules_CounterFlow1DFlame grammar;
		dictionary.SetGrammar(grammar);

		is_species_ = false;
		if (dictionary.CheckOption("@Species") == true)
			dictionary.ReadString("@Species", species_name_);
		if (species_name_ != "none") is_species_ = true;

		if (dictionary.CheckOption("@Temperature") == true)
			dictionary.ReadBool("@Temperature", is_temperature_);

		if (dictionary.CheckOption("@SpeciesProfile") == true)
		{
			std::string name_of_subdictionary;
			dictionary.ReadDictionary("@SpeciesProfile", name_of_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_subdictionary), x, y, x_variable, y_variable);

			species_np_ = x.Size();
			species_x_.resize(species_np_);
			species_y_.resize(species_np_);
			species_mean_y_ = 0.;
			for (unsigned int i = 0; i < species_np_; i++)
			{
				species_x_(i) = x[i + 1];
				species_y_(i) = y[i + 1];
				species_mean_y_ += species_y_(i);
			}
			species_mean_y_ /= static_cast<double>(species_np_);
		}

		if (dictionary.CheckOption("@TemperatureProfile") == true)
		{
			std::string name_of_subdictionary;
			dictionary.ReadDictionary("@TemperatureProfile", name_of_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_subdictionary), x, y, x_variable, y_variable);

			t_np_ = x.Size();
			t_x_.resize(t_np_);
			t_y_.resize(t_np_);
			t_mean_y_ = 0.;
			for (unsigned int i = 0; i < t_np_; i++)
			{
				t_x_(i) = x[i + 1];
				t_y_(i) = y[i + 1];
				t_mean_y_ += t_y_(i);
			}
			t_mean_y_ /= static_cast<double>(t_np_);
		}
	}
}

#endif // OpenSMOKE_OptimizationRules_CounterFlow1DFlame_H