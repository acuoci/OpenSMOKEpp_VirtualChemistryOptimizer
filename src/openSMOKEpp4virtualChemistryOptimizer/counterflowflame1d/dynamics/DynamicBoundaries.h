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

#ifndef OpenSMOKE_DynamicBoundaries_H
#define OpenSMOKE_DynamicBoundaries_H

//#include "counterflowflame1d/OpenSMOKE_CounterFlowFlame1D.h"

namespace OpenSMOKE
{
	//!  A class to force dynamic boundaries in counter-flow diffusion flames according to different laws
	/*!
	A class to force dynamic boundaries in counter-flow diffusion flames according to different laws
	*/

	class DynamicBoundaries
	{

	private:

		enum DynamicBoundaries_Type {	NONE, TEMPERATURE_VELOCITY_SLOPES, 
										FIXED_STRAINRATE_FUEL_TEMPERATURE_SLOPE, FIXED_STRAINRATE_OX_TEMPERATURE_SLOPE,
										FIXED_BALANCE_FUEL_VELOCITY_SLOPE, FIXED_BALANCE_OX_VELOCITY_SLOPE,
										SIN_VELOCITY_INPHASE };
		
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param path_output			path to output folder
		*/
		DynamicBoundaries(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, const boost::filesystem::path path_output);

		/**
		*@brief Setup from a dictionary
		*@param dictionary dictionary name
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Setup (partial)
		*@param length length of the computational domain (in m)
		*@param P_Pa pressure (in Pa)
		*@param n_geometry type of geometry (2=planar, 3=axysymmetric)
		*/
		void SetFlame(const double length, const double P_Pa, const double n_geometry);

		/**
		*@brief Set initial conditions on the fuel side
		*@param T temperature (in K)
		*@param v velocity (in m/s)
		*@param omega mass fractions
		*/
		void SetInitialFuelSide(const double T, const double v, const double* omega);

		/**
		*@brief Set initial conditions on the oxidizer side
		*@param T temperature (in K)
		*@param v velocity (in m/s)
		*@param omega mass fractions
		*/
		void SetInitialOxidizerSide(const double T, const double v, const double* omega);

		/**
		*@brief Set virtual chemistry
		*@param virtual_chemistry virtual chemistry
		*/
		void SetVirtualChemistry(OpenSMOKE::VirtualChemistry* virtual_chemistry);

		/**
		*@brief Completes the setup phase
		*/
		void CompleteSetup();

		/**
		*@brief Prepares the output files
		*/
		void PrepareOutputFiles();

		/**
		*@brief Prints/updates the dynamic output file
		*@param t time (in s)
		*@param Tmax maximum temperature (in K)
		*/
		void Print(const double t, const double Tmax);

		/**
		*@brief Prints a summary on the screen
		*@param out output stream
		*/
		void Summary(std::ostream& out);

		/**
		*@brief Updates the boundary conditions
		*@param t time (in s)
		*@param U_fuel mass flow rate on the fuel side (in kg/m2/s)
		*@param U_ox mass flow rate on the oxidizer side (in kg/m2/s)
		*@param T_fuel temperature on the fuel side (in K)
		*@param T_ox temperature on the oxidizer side (in K)
		*@param rho_fuel density on the fuel side (in kg/m3)
		*@param rho_ox density on the oxidizer side (in kg/m3)
		*@param omega_fuel mass fractions on the fuel side
		*@param omega_ox mass fractions on the oxidizer side
		*/
		void UpdateBoundaryConditions(	const double t, double& U_fuel, double& U_ox, double& T_fuel, double& T_ox,
										double& rho_fuel, double& rho_ox, double* omega_fuel, double* omega_ox);

		/**
		*@brief Returns the strain rate
		*@return the strain rate (in 1/s)
		*/
		double a() const { return a_; }

		/**
		*@brief Returns the position of the stagnation plane (from the fuel nozzle)
		*@return the position of the stagnation plane (in m)
		*/
		double xst() const { return xst_; }

		/**
		*@brief Returns the momentum umbalance factor
		*@return the momentum umbalance factor (1: balance, >1: fuel stronger, <1: oxidizer stronger)
		*/
		double beta() const { return beta_; }

		/**
		*@brief Returns true if compact output is required
		*@return true if compact output is required
		*/
		bool compact_output() const { return compact_output_; }

		/**
		*@brief Returns the number of steps governing the output
		*@return number of steps governing the output
		*/
		int n_steps_output() const { return n_steps_output_; }

		/**
		*@brief Returns the list of times at which flame snapshots are required
		*@return list of times at which flame snapshots are required (in s)
		*/
		const std::vector<double>& snapshot_list_of_times() const { return snapshot_list_of_times_; }

	private:

		/**
		*@brief Check input options
		*/
		void CheckInput();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamics map

		DynamicBoundaries_Type type_;			//!< type of dynamic boundary condition

		boost::filesystem::path path_output_;	//!< path to output folder

		double frequency_fuel_;					//!< frequency of sinusoidal signal (in Hz)
		double semi_amplitude_fuel_;			//!< semi-amplitue (relative) of sinusoidal signal
		double t_delay_fuel_;					//!< delay of sinusoidal signal (in s)

		double frequency_ox_;					//!< frequency of sinusoidal signal (in Hz)
		double semi_amplitude_ox_;				//!< semi-amplitue (relative) of sinusoidal signal
		double t_delay_ox_;						//!< delay of sinusoidal signal (in s)

		double slope_fuel_velocity_;			//!< slope of fuel velocity (in m/s2)
		double slope_ox_velocity_;				//!< slope of oxidizer velocity (in m/s2)
		double slope_fuel_temperature_;			//!< slope of fuel temperature (in K/s)
		double slope_ox_temperature_;			//!< slope of oxidizer temperature (in K/s)

		bool is_active_;						//!< true if the dynamic boundary condition manager is active

		double P_Pa_;							//!< pressure (in Pa)
		double length_;							//!< length of computational domain
		double n_geometry_;						//!< geometry type (2=planar, 3=axysymmetric)

		double a_0_;							//!< initial strain rate (in 1/s)
		double a_;								//!< current strain rate (in 1/s)

		double a_seshadri_0_;					//!< initial strain rate, Seshadri's definition (in 1/s)
		double a_seshadri_;						//!< current strain rate, Seshadri's definition (in 1/s)

		double xst_0_;							//!< initial position of stagnation plane (in m)
		double xst_;							//!< current position of stagnation plane (in m)

		double beta_0_;							//!< initial momentum unbalance coefficient (1=perfect balance)
		double beta_;							//!< current momentum unbalance coefficient (1=perfect balance)

		double T_fuel_0_;						//!< initial fuel temperature (in K)
		double T_ox_0_;							//!< initial oxidizer temperature (in K)
		double v_fuel_0_;						//!< initial fuel velocity (in m/s)
		double v_ox_0_;							//!< initial oxidizer velocity (in m/s)

		double T_fuel_;							//!< current fuel temperature (in K)
		double T_ox_;							//!< current oxidizer temperature (in K)
		double v_fuel_;							//!< current fuel velocity (in m/s)
		double v_ox_;							//!< current oxidizer velocity (in m/s)

		Eigen::VectorXd omega_fuel_0_;			//!< initial fuel mass fractions
		Eigen::VectorXd omega_ox_0_;			//!< initial oxidizer mass fractions
		Eigen::VectorXd x_fuel_0_;				//!< initial fuel mole fractions
		Eigen::VectorXd x_ox_0_;				//!< initial oxidizer mole fractions

		Eigen::VectorXd omega_fuel_;			//!< current fuel mass fractions
		Eigen::VectorXd omega_ox_;				//!< current oxidizer mass fractions
		Eigen::VectorXd x_fuel_;				//!< current fuel mole fractions
		Eigen::VectorXd x_ox_;					//!< current oxidizer mole fractions

		double mw_fuel_0_;						//!< initial fuel molecular weight (in kg/kmol)
		double mw_ox_0_;						//!< initial oxidizer molecular weight (in kg/kmol)
		double mw_fuel_;						//!< current fuel molecular weight (in kg/kmol)
		double mw_ox_;							//!< current oxidizer molecular weight (in kg/kmol)

		double rho_fuel_0_;						//!< initial fuel density (in kg/m3)	
		double rho_ox_0_;						//!< initial oxidizer density (in kg/m3)
		double rho_fuel_;						//!< current fuel density (in kg/m3)
		double rho_ox_;							//!< urrent oxidizer density (in kg/m3)

		double U_fuel_0_;						//!< initial fuel mass flow rate (in kg/m2/s)	
		double U_ox_0_;							//!< initial oxidizer mass flow rate (in kg/m2/s)
		double U_fuel_;							//!< current fuel mass flow rate (in kg/m2/s)
		double U_ox_;							//!< current oxidizer mass flow rate (in kg/m2/s)

		double end_time_;								//!< total time of integration (in s)
		bool compact_output_;							//!< compact output
		int n_steps_output_;							//!< number of steps governing the output
		std::vector<double> snapshot_list_of_times_;	//!< list of output times (in s)

		std::ofstream fOutputDynamics_;					//!< output file

		OpenSMOKE::VirtualChemistry* virtual_chemistry_;
		bool is_virtual_chemistry_;

	};
}

#include "DynamicBoundaries.hpp"

#endif /* OpenSMOKE_DynamicBoundaries_H */

