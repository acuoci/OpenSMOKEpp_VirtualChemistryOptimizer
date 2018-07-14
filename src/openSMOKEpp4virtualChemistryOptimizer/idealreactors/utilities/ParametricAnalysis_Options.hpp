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

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"


namespace OpenSMOKE
{
	class Grammar_ParametricAnalysis_Options : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", 
																OpenSMOKE::SINGLE_STRING, 
																"Type of parameter to be analyzed: residence-time || temperature || pressure || moles || masses", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfValues", 
																OpenSMOKE::VECTOR_STRING, 
																"List of values (together with units of measure, if any)",
																true,
																"@MinimumValue",
																"none",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumValue", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Minimum value of parameter", 
																true,
																"@ListOfValues",
																"@MaximumValue",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumValue", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Maximum value of parameter", 
																false,
																"none",
																"@NumberOfPoints",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPoints", 
																OpenSMOKE::SINGLE_INT,
																"Number of points", 
																false,
																"none",
																"@MinimumValue",
																"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfThreads",
															   OpenSMOKE::SINGLE_STRING,
															   "Number of threads (in case OpenMP is enabled)",
															   false) );

		}
	};

	ParametricAnalysis_Options::ParametricAnalysis_Options()
	{
		enabled_        = false;						//!< enabled on/off
		parameter_type_ = PARAMETER_TYPE_TEMPERATURE;	//!< parameter type
		number_of_points_ = 10;							//!< number of points (default 10)
		maximum_value_ = 0.;							//!< minimum value of parameter
		minimum_value_ = 0.;							//!< maximum value of parameter
		list_of_values_.resize(0);						//!< list of values
		number_of_threads_ = 1;							//!< number of threads (is OpenMP is enabled)
	}
	
	void ParametricAnalysis_Options::SetNumberOfPoints(const int number_of_points)
	{
		number_of_points_ = number_of_points;
	}
	
	void ParametricAnalysis_Options::SetMinimumValue(const double minimum_value)
	{
		minimum_value_ = minimum_value;
	}	
	
	void ParametricAnalysis_Options::SetMaximumValue(const double maximum_value)
	{
		maximum_value_ = maximum_value;
	}

	void ParametricAnalysis_Options::SetNumberOfThreads(const int number_of_threads)
	{
		number_of_threads_ = number_of_threads;
	}

	void ParametricAnalysis_Options::SetNumberOfThreads(const std::string number_of_threads)
	{
		if (boost::iequals(number_of_threads, "max"))
		{
			#if defined(_OPENMP)
				SetNumberOfThreads(omp_get_max_threads());
			#endif
		}
		else
		{
			SetNumberOfThreads(boost::lexical_cast<int>(number_of_threads));
		}		
	}

	void ParametricAnalysis_Options::SetListOfValues(const std::vector<double> list_of_values)
	{
		list_of_values_ = list_of_values;
		
		// reorder and remove duplicates
		std::sort(list_of_values_.begin(), list_of_values_.end()); 
		std::vector<double>::iterator last = std::unique(list_of_values_.begin(), list_of_values_.end());
		list_of_values_.erase(last, list_of_values_.end());
		
		minimum_value_ = list_of_values_.front();
		maximum_value_ = list_of_values_.back();
		number_of_points_ = boost::lexical_cast<int>(list_of_values_.size());
	}		

	void ParametricAnalysis_Options::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_ParametricAnalysis_Options grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Type") == true)
		{
			std::string value;
			dictionary.ReadString("@Type", value);
			if (value == "temperature")
				parameter_type_ = PARAMETER_TYPE_TEMPERATURE;
			else if (value == "pressure")
				parameter_type_ = PARAMETER_TYPE_PRESSURE;
			else if (value == "residence-time")
				parameter_type_ = PARAMETER_TYPE_TIME;
			else if (value == "masses")
				parameter_type_ = PARAMETER_TYPE_MASSES;
			else if (value == "moles")
				parameter_type_ = PARAMETER_TYPE_MOLES;
			else
				OpenSMOKE::FatalErrorMessage("Unknown parametric type: " + value);
		}

		if (parameter_type_ == PARAMETER_TYPE_MASSES || parameter_type_ == PARAMETER_TYPE_MOLES)
		{
			if (dictionary.CheckOption("@ListOfValues") == true)
			{
				std::vector<std::string> values;
				dictionary.ReadOption("@ListOfValues", values);

				number_of_points_ = 1;
				for (unsigned int i = 0; i < values.size(); i++)
				if (values[i] == ":")
					number_of_points_++;

				list_of_species_.resize(number_of_points_);
				list_of_compositions_.resize(number_of_points_);

				unsigned int index = 0;
				unsigned int i = 0;
				for (;;)
				{
					if (i >= values.size())
						break;

					if (values[i] == ":")
					{
						index++;
						i++;
						continue;
					}
					else
					{
						list_of_species_[index].push_back(values[i]);
						list_of_compositions_[index].push_back(boost::lexical_cast<double>(values[i + 1]));
						i += 2;
					}
				}

				// Print on the screen the list of values to be analyzed
				std::cout << "Parametric Analysis:" << std::endl;
				for (int i = 0; i < number_of_points_; i++)
				{
					std::cout << i + 1 << "\t";
					for (unsigned int j = 0; j < list_of_species_[i].size(); j++)
						std::cout << list_of_species_[i][j] << " " << list_of_compositions_[i][j] << " ";
					std::cout << std::endl;
				}
			}
			else
				OpenSMOKE::FatalErrorMessage("Parameteric compositions can be specified only using the @ListOfValues options");
		}
		else
		{
			if (dictionary.CheckOption("@ListOfValues") == true)
			{
				std::vector<std::string> values;
				dictionary.ReadOption("@ListOfValues", values);

				std::string units = values.back();
				if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE)
				{
					double conversion_factor = 0.;
					if (values.back() == "K") conversion_factor = 0.;
					else if (values.back() == "C") conversion_factor = 273.15;
					else OpenSMOKE::FatalErrorMessage("Wrong units for temperature");

					list_of_values_.resize(values.size() - 1);
					for (unsigned int i = 0; i < values.size() - 1; i++)
						list_of_values_[i] = (boost::lexical_cast<double>(values[i]) + conversion_factor);
				}
				else if (parameter_type_ == PARAMETER_TYPE_PRESSURE)
				{
					double conversion_factor = 1.;
					if (values.back() == "Pa") conversion_factor = 1.;
					else if (values.back() == "atm") conversion_factor = 101325.;
					else if (values.back() == "bar") conversion_factor = 100000.;
					else OpenSMOKE::FatalErrorMessage("Wrong units for pressure");

					list_of_values_.resize(values.size() - 1);
					for (unsigned int i = 0; i < values.size() - 1; i++)
						list_of_values_[i] = (boost::lexical_cast<double>(values[i])*conversion_factor);
				}
				else if (parameter_type_ == PARAMETER_TYPE_TIME)
				{
					double conversion_factor = 1.;
					if (values.back() == "s") conversion_factor = 1.;
					else if (values.back() == "ms") conversion_factor = 1e-3;
					else OpenSMOKE::FatalErrorMessage("Wrong units for time");

					list_of_values_.resize(values.size() - 1);
					for (unsigned int i = 0; i < values.size() - 1; i++)
						list_of_values_[i] = (boost::lexical_cast<double>(values[i])*conversion_factor);
				}

				SetListOfValues(list_of_values_);
			}
			else
			{
				if (dictionary.CheckOption("@NumberOfPoints") == true)
				{
					dictionary.ReadInt("@NumberOfPoints", number_of_points_);
					if (number_of_points_ <= 1)
						OpenSMOKE::FatalErrorMessage("The number of points must be at least equal to 2");
				}

				if (dictionary.CheckOption("@MinimumValue") == true)
				{
					std::string units;
					dictionary.ReadMeasure("@MinimumValue", minimum_value_, units);
					if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE)
					{
						if (units == "K")			minimum_value_ += 0.;
						else if (units == "C")		minimum_value_ += 273.15;
						else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
					}
					else if (parameter_type_ == PARAMETER_TYPE_PRESSURE)
					{
						if (units == "Pa")			minimum_value_ *= 1.;
						else if (units == "atm")	minimum_value_ *= 101325.;
						else if (units == "bar")	minimum_value_ *= 100000.;
						else OpenSMOKE::FatalErrorMessage("Unknown pressure units");
					}
					else if (parameter_type_ == PARAMETER_TYPE_TIME)
					{
						if (units == "s")			minimum_value_ *= 1.;
						else if (units == "ms")		minimum_value_ *= 1e-3;
						else OpenSMOKE::FatalErrorMessage("Unknown time units");
					}
				}

				if (dictionary.CheckOption("@MaximumValue") == true)
				{
					std::string units;
					dictionary.ReadMeasure("@MaximumValue", maximum_value_, units);
					if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE)
					{
						if (units == "K")			maximum_value_ += 0.;
						else if (units == "C")		maximum_value_ += 273.15;
						else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
					}
					else if (parameter_type_ == PARAMETER_TYPE_PRESSURE)
					{
						if (units == "Pa")			maximum_value_ *= 1.;
						else if (units == "atm")	maximum_value_ *= 101325.;
						else if (units == "bar")	maximum_value_ *= 100000.;
						else OpenSMOKE::FatalErrorMessage("Unknown pressure units");
					}
					else if (parameter_type_ == PARAMETER_TYPE_TIME)
					{
						if (units == "s")			maximum_value_ *= 1.;
						else if (units == "ms")		maximum_value_ *= 1e-3;
						else OpenSMOKE::FatalErrorMessage("Unknown time units");
					}
				}

				// Final Checks
				if (minimum_value_ >= maximum_value_)
					OpenSMOKE::FatalErrorMessage("Parameter minimum value must be strictly smaller than maximum value");

				// Populating the list of parameters
				list_of_values_.resize(number_of_points_);
				const double step = (maximum_value_ - minimum_value_) / double(number_of_points_ - 1);
				list_of_values_[0] = minimum_value_;
				for (unsigned int i = 1; i < list_of_values_.size(); i++)
					list_of_values_[i] = list_of_values_[i - 1] + step;

				// Print on the screen the list of values to be analyzed
				std::cout << "Parametric Analysis:" << std::endl;
				for (unsigned int i = 0; i < list_of_values_.size(); i++)
					std::cout << i + 1 << " " << list_of_values_[i] << std::endl;
			}
		}

		if (dictionary.CheckOption("@NumberOfThreads") == true)
		{
			std::string value;
			dictionary.ReadString("@NumberOfThreads", value);
			SetNumberOfThreads(value);
		}

		enabled_ = true;
	}

}
