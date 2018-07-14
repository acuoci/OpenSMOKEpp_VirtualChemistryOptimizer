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

#ifndef OpenSMOKE_IgnitionDelayTimes_Analyzer_H
#define	OpenSMOKE_IgnitionDelayTimes_Analyzer_H

#include "Grammar_IgnitionDelayTimes.h"

namespace OpenSMOKE
{

	//!  A class for estimating ignition delay times from temperature, pressure and composition temporal profiles
	/*!
	The purpose of this class is to estimate ignition delay times from temperature, pressure and composition 
	temporal profiles in non-isothermal batch reactors (constant volume, constant pressure or user-defined volume)
	*/

	class IgnitionDelayTimes_Analyzer
	{
	public:

		/**
		*@brief Default constructor
		*/
		IgnitionDelayTimes_Analyzer();

		/**
		*@brief Setup the options from an external dictionary
		*@param dictionary external disctionary containing the options
		*@param thermodynamicsMap map containing the thermodynamic data
		*/
		template<typename Thermodynamics>
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMap);

		/**
		*@brief Returns true if the analyzer is currently active
		*/
		bool is_active() const { return is_active_; }

		/**
		*@brief Extracts the ignition delay time from a provided set of temperature, pressure and composition (mole fractions)
		*@param t the current time (in s)
		*@param T the current temperature (in K)
		*@param P the current pressure (in Pa)
		*@param x the current molar fractions
		*/
		void Analyze(const double t, const double T, const double P, const double* x);

		/**
		*@brief Prints the calculated ignition delay times on a file
		*@param file_name path to the output file
		*/
		void PrintOnFile(const boost::filesystem::path file_name);

		/**
		*@brief Prints the header line for parameteric analysis
		*@param counter index of column
		*@param fOutput the output file
		*/
		void PrintHeaderLine(unsigned int& counter, std::ostream& fOutput);

		/**
		*@brief Prints the ignition delay times on a file for parametric analysis
		*@param fOutputthe output file
		*/
		void Print(std::ostream& fOutput);

		/**
		*@brief Resets the ignition delay times
		*/
		void Reset();

		/**
		*@brief Returns the idt (in s) based on the maximum temperature slope
		*/
		double temperature_slope_tau() const { return temperature_slope_tau_;  }

		/**
		*@brief Returns the idt (in s) based on the temperature increase
		*/
		double temperature_increase_tau() const { return temperature_increase_tau_; }


	private:

		bool is_active_;								//!< true if the analyzer is active
		bool is_pressure_;								//!< true if the user want to estimate the ignition delay time based on the pressure slope
		bool is_temperature_;							//!< true if the user want to estimate the ignition delay time based on the temperature slope
		bool is_species_slope_;							//!< true if the user want to estimate the ignition delay time based on the species slopes

		double x_threshold_;							//!< minimum mole fraction for evaluation of i.d.t. on the basis of max mole fractions of species
		double time_minimum_;							//!< minimum time (in s) for starting the evalution of i.d.t.
		double time_minimum_interval_;					//!< minimum time interval (in s) for the evalution of i.d.t. based on slopes

		std::vector<unsigned int> index_species_;		//!< indices (0-index based) of selected species for evaluation of i.d.t.
		std::vector<std::string> species_name_;			//!< names of selected species for evaluation of i.d.t.

		// Maximum slope criteria
		double temperature_slope_max_;					//!< current max value of temperature slope (in K/s)
		double temperature_slope_tau_;					//!< current value of i.d.t. based on the temperature slope (in s)
		double pressure_slope_max_;						//!< current max value of pressure slope (in Pa/s)
		double pressure_slope_tau_;						//!< current value of i.d.t. based on the pressure slope (in s)
		std::vector<double> species_slope_max_;			//!< current max values of mole fraction slopes (in 1/s)
		std::vector<double> species_slope_max_tau_;		//!< current values of i.d.t. based on the species slope (in s)

		// Maximum criteria
		std::vector<double> species_max_;				//!< current max values of mole fractions
		std::vector<double> species_max_tau_;			//!< current values of i.d.t. based on the species max mole fractions (in s)

		// Old values (for derivatives)
		double tOld_;									//!< previous time (in s)
		double TOld_;									//!< previous temperature (in K)
		double POld_;									//!< previous pressure (in Pa)
		std::vector<double> xOld_;						//!< previous mole fractions (in s)

		// Temperature increase
		double temperature_increase_;					//!< user-defined temperature increase for i.d.t.
		double temperature_increase_tau_;				//!< current value of i.d.t. based on the temperature increase
		double T0_;										//!< inlet/initial temperature (in K)
	};
}

#include "IgnitionDelayTimes_Analyzer.hpp"


#endif	/* OpenSMOKE_IgnitionDelayTimes_Analyzer_H */
