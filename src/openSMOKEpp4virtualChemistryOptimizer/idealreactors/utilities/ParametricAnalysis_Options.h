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

#ifndef OpenSMOKE_ParametricAnalysis_Options_H
#define	OpenSMOKE_ParametricAnalysis_Options_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "math/PhysicalConstants.h"

namespace OpenSMOKE
{
	class ParametricAnalysis_Options
	{
	public:

		enum ParametricAnalysisType {	PARAMETER_TYPE_NONE, 
										PARAMETER_TYPE_TEMPERATURE, 
										PARAMETER_TYPE_PRESSURE,
										PARAMETER_TYPE_TIME, 
										PARAMETER_TYPE_MASSES, 
										PARAMETER_TYPE_MOLES };

	public:

		/**
		*@brief Default constructor
		*/
		ParametricAnalysis_Options();

		/**
		*@brief Initializes the object from an external dictionary
		*@param dictionary external dictionary
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		
		/**
		*@brief Sets the number of points to be analyzed
		*/
		void SetNumberOfPoints(const int number_of_points);
		
	    /**
		*@brief Sets the parameter minimum value to be analyzed
		*/
		void SetMinimumValue(const double minimum_value);		
		
	    /**
		*@brief Sets the parameter maximum value to be analyzed
		*/
		void SetMaximumValue(const double maximum_value);
		
	    /**
		*@brief Sets the list of parameter values to be analyzed
		*/
		void SetListOfValues(const std::vector<double> list_of_values);		

		/**
		*@brief Sets the number of threads
		*/
		void SetNumberOfThreads(const int number_of_threads);

		/**
		*@brief Sets the number of threads
		*/
		void SetNumberOfThreads(const std::string number_of_threads);

		/**
		*@brief Returns is parametric analysi was enabled or not
		*/
		bool enabled() const { return enabled_; }

		/**
		*@brief Returns the type of parameter
		*/
		ParametricAnalysisType parameter_type() const { return parameter_type_; }
		
		/**
		*@brief Returns the number of threads to be used if OpenMP is enabled
		*/
		int number_of_threads() const { return number_of_threads_; }

		/**
		*@brief Returns the number of points to be analyzed
		*/
		int number_of_points() const { return number_of_points_; }

		/**
		*@brief Returns the parameter minimum value to be analyzed
		*/
		double minimum_value() const { return minimum_value_; }		
		
		/**
		*@brief Returns the parameter maximum value to be analyzed
		*/
		double maximum_value() const { return maximum_value_; }	
		
	    /**
		*@brief Returns the list of parameter values to be analyzed
		*/
		const std::vector<double>& list_of_values() const { return list_of_values_; }

		/**
		*@brief Returns the list of species
		*/
		const std::vector<std::string>& list_of_species(const unsigned int index) const { return list_of_species_[index]; }
	
		/**
		*@brief Returns the list of compositions
		*/
		const std::vector<double>& list_of_compositions(const unsigned int index) const { return list_of_compositions_[index]; }

	private:

		bool enabled_;							//!< enabled on/off
		ParametricAnalysisType parameter_type_;	//!< parameter type
		int number_of_points_;					//!< number of points (default 10)
		double maximum_value_;					//!< minimum value of parameter
		double minimum_value_;					//!< maximum value of parameter
		std::vector<double>	list_of_values_;	//!< list of values
		int number_of_threads_;					//!< number of threads to be used if OpenMP is enabled
		
		std::vector< std::vector<std::string> >	list_of_species_;		//!< list of species
		std::vector< std::vector<double> >		list_of_compositions_;	//!< list of compositions

	};
}

#include "ParametricAnalysis_Options.hpp"

#endif	/* OpenSMOKE_ParametricAnalysis_Options_H */

