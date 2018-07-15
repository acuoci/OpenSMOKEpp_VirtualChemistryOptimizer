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

#ifndef OpenSMOKE_CounterFlowFlame1D_H
#define OpenSMOKE_CounterFlowFlame1D_H

#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#endif

// External parameters
#include "rapidxml.hpp"
#include <Eigen/Dense>

// Basic classes
#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "math/OpenSMOKEBandMatrix.h"

// Numerical parameters
#include "math/native-dae-solvers/parameters/DaeSolver_Parameters.h"
#include "math/native-nls-solvers/parameters/NonLinearSolver_Parameters.h"
#include "math/native-nls-solvers/parameters/FalseTransientSolver_Parameters.h"

//Sensitivity analysis
#include "utilities/sensitivityanalysis/SensitivityAnalysis_Options.h"

// Fixed profile
#include "utilities/FixedProfile.h"

// Dynamic boundaries
#include "dynamics/DynamicBoundaries.h"

namespace OpenSMOKE
{
	//!  A class to solve laminar premixed flat (1D) flames
	/*!
	This class provides the tools to solve laminar premixed flat (1D) flames
	*/

	class OpenSMOKE_CounterFlowFlame1D
	{

	public:

		enum Simulation_Type { SIMULATION_TYPE_M, SIMULATION_TYPE_Y, SIMULATION_TYPE_YTM, SIMULATION_TYPE_TM, SIMULATION_TYPE_YT, SIMULATION_TYPE_YM, SIMULATION_TYPE_HMOM, SIMULATION_TYPE_Y_HMOM, SIMULATION_TYPE_YT_HMOM, SIMULATION_TYPE_CSI };
		enum Solver_Type { SOLVER_TYPE_COUNTERFLOW_DIFFUSION, SOLVER_TYPE_BURNER_STABILIZED_STAGNATION };
		
	private:

		enum InitialProfile_Type { INITIAL_PROFILE_TYPE_TRIANGULAR, INITIAL_PROFILE_TYPE_PLATEAU, INITIAL_PROFILE_TYPE_LINEAR};
		enum MassDiffusionCoefficients_Type { MASS_DIFFUSION_COEFFICIENTS_TYPE_MOLECULAR_THEORY_GASES, MASS_DIFFUSION_COEFFICIENTS_TYPE_LEWIS_NUMBERS };

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*@param transportMap		reference to the transport map
		*@param grid				reference to 1D grid
		*/
		OpenSMOKE_CounterFlowFlame1D( OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
										  OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
										  OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
										  OpenSMOKE::Grid1D& grid);

		/**
		*@brief Sets solver type
		*@param solver_type solver type: SOLVER_TYPE_BURNERSTABILIZED | SOLVER_TYPE_FLAMESPEED
		*/
		void SetSolverType(const std::string solver_type);

		/**
		*@brief Sets the fuel side mixture
		*@param Tfuel fuel temperature [K]
		*@param P_Pa pressure [Pa]
		*@param omegaFuel composition (in mass fractions) of fuel mixture
		*@param fuel_velocity velocity [m/s]
		*/
		void SetFuelSide(const double Tfuel, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaFuel, const double fuel_velocity);

		/**
		*@brief Sets the oxidizer side mixture
		*@param Toxidizer oxidizer temperature [K]
		*@param P_Pa pressure [Pa]
		*@param omegaOxidizer composition (in mass fractions) of oxidizer mixture
		*@param oxidizer_velocity velocity [m/s]
		*/
		void SetOxidizerSide(const double Toxidizer, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaOxidizer, const double oxidizer_velocity);

		/**
		*@brief Sets the peak mixture (needed by triangular or plateau profiles)
		*@param Tpeak peak temperature [K]
		*@param omegaPeak composition (in mass fractions) of peak mixture
		*/
		void SetPeakMixture(const double Tpeak, const OpenSMOKE::OpenSMOKEVectorDouble& omegaPeak);

		/**
		*@brief Sets the name of the output folder
		*@param output_folder name of the output folder
		*/
		void SetOutputFolder(const boost::filesystem::path output_folder);

		/**
		*@brief If turned on, after the DAE solution, the solution of the NLS is attempted for accuracy
		*@param use_nls_solver if true, after the DAE solution, the solution of the NLS is attempted for accuracy
		*/
		void SetUseNlsSolver(const bool use_nls_solver);

		/**
		*@brief If turned on, the solution is calculated through successive transient solutions (robust way to get the solution)
		*@param use_dae_solver if true, the transient solutions are calculated, otherwise only non-linear systems are attempted
		*/
		void SetUseDaeSolver(const bool use_dae_solver);

		/**
		*@brief Turns on/off the Soret effect
		*@param soret if true, the Soret effect is accounted for
		*/
		void SetSoret(const bool soret);

		/**
		*@brief Sets the simplified transport properties
		*@param flag if true, simplfiied properties are adopted
		*/
		void SetSimplifiedTransportProperties(const bool flag);

		/**
		*@brief Turns on/off the deposition wall for soot in burner stabilized stagnation flames
		*@param flag if true, the deposition wall is turned on
		*/
		void SetDepositionWall(const bool flag);

		/**
		*@brief Turns on/off the radiative heat transfer
		*@param radiative_heat_transfer if true, the radiative heat transfer is accounted for
		*/
		void SetRadiativeHeatTransfer(const bool radiative_heat_transfer);

		/**
		*@brief Sets the environment temperature
		*@param environment_temperature the environment temperature [K]
		*/
		void SetEnvironmentTemperature(const double environment_temperature);

		/**
		*@brief Sets the soot analyzer fo kinetic mechanisms developed by the CRECK Group at Politecnico di Milano
		*@param polimi_soot_analyzer pointer to the soot analyzer object
		*/
		void SetPolimiSoot(OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_analyzer);

		/**
		*@brief Sets the virtual chemistry approach according to the formulation proposed by B. Fiorina (EC2M)
		*@param virtual_chemistry pointer to the virtual chemistry object
		*/
		void SetVirtualChemistry(OpenSMOKE::VirtualChemistry* virtual_chemistry);

		/**
		*@brief Sets the post processor (reaction and formation rates) "on the fly"
		*@param polimi_soot_analyzer pointer to the post processor (reaction and formation rates) "on the fly"
		*/
		void SetOnTheFlyPostProcessing(OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing);

		/**
		*@brief Sets the dynamic boundaries
		*@param dynamic_boundaries pointer to the dynamic boundary manager
		*/
		void SetDynamicBoundaries(OpenSMOKE::DynamicBoundaries* dynamic_boundaries);
			
		/**
		*@brief Sets a user defined, fixed temperature profile (automatically fixes also the initial profile)
		*@param x axial coordinates [m]
		*@param T temperatures [K]
		*/
		void SetFixedTemperatureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& T);

		/**
		*@brief Sets a user defined, initial temperature profile
		*@param x axial coordinates [m]
		*@param T temperatures [K]
		*/
		void SetInitialTemperatureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& T);

		/**
		*@brief Sets the position of the initial peak of temperature
		*@param x_peak initial position of peak of temperature [m]
		*/
		void SetPeakPosition(const double x_peak);

		/**
		*@brief Sets the width of the initial mixing zone
		*@param width_mixing width of the initial mixing zone [m]
		*/
		void SetMixingZoneWidth(const double width_mixing);

		/**
		*@brief Sets the type of initial profiles of species and temperature
		*@param type type of initial profiles of species and temperature (triangular | linear | plateau)
		*/
		void SetInitialProfileType(const std::string type);

		/**
		*@brief Sets the starting guess value of the eigenvalue
		*@param H the starting guess value [kg/m3/s2]
		*/
		void SetEigenValueStartingGuess(const double H);

		/**
		*@brief Sets the radial gradient on the fuel side
		*@param G the radial gradient on the fuel side [kg/m3/s]
		*/
		void SetRadialGradientFuelSide(const double G);

		/**
		*@brief Sets the radial gradient on the oxidizer side
		*@param G the radial gradient on the oxidizer side [kg/m3/s]
		*/
		void SetRadialGradientOxidizerSide(const double G);

		/**
		*@brief Sets the planar symmetry
		*@param flag if true, the planar symmetry is adopted (otherwise, the cylyndrical symmetry is adopted)
		*/
		void SetPlanarSymmetry(const bool flag);

		/**
		*@brief Sets the optional multiplying factor for all the reaction rates (default: 1)
		*@param value the optional multiplying factor for all the reaction rates
		*/
		void GasReactionRateMultiplier(const double value);

		/**
		*@brief Sets the number of Lewis for each species: Le = alpha/gamma
		*@param lewis_numbers the Lewis numbers for every species
		*/
		void SetLewisNumbers(const std::vector<double> lewis_numbers);

		/**
		*@brief Prepares the solver for burner stabilized calculations
		*/
		void SetupForBurnerStabilized();

		/**
		*@brief Solves the flame problem from scratch
		*@param dae_parameter parameters governing the solution of DAE systems
		*@param nls_parameters parameters governing the solution of NL systems
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int SolveFlameFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameter,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);

		/**
		*@brief Solves the problem from an existing solution
		*@param dae_parameter parameters governing the solution of DAE systems
		*@param nls_parameters parameters governing the solution of NL systems
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int SolveFromExistingSolution(DaeSMOKE::DaeSolver_Parameters& dae_parameter,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);

		/**
		*@brief Calculates the soot formation according to the Hybrid Method of Moments
		*@param hmom class defining the HMOM submodels to be applied
		*@param dae_parameter parameters governing the solution of DAE systems
		*@param nls_parameters parameters governing the solution of NL systems
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int SolveHMOMFromExistingSolution(OpenSMOKE::HMOM& hmom,
			DaeSMOKE::DaeSolver_Parameters& dae_parameters,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);


		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param sensitivity_options options governing the sensitivity analysis
		*/
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);

		/**
		*@brief Updates the current solution/unknowns
		*@param phi new solution/unknowns
		*/
		void Update(const std::vector<Eigen::VectorXd>& phi);

		/**
		*@brief Returns the size of a single block of equations
		*/
		int BlockDimensions() const;

		/**
		*@brief Returns the total number of equations
		*/
		int NumberOfEquations() const;

		/**
		*@brief Returns the grid
		*/
		const OpenSMOKE::Grid1D& grid() const { return grid_; }

		/**
		*@brief Returns the upper band size
		*/
		unsigned int UpperBand() const { return (2*BlockDimensions()-1); }

		/**
		*@brief Returns the lower band size
		*/
		unsigned int LowerBand() const { return (2*BlockDimensions()-1); }

		/**
		*@brief Returns the current solution/unknowns
		*@param v the current solution/unknowns
		*/
		void UnknownsVector(double* v);

		/**
		*@brief Updates the current solution/unknowns
		*@param v new solution/unknowns
		*/
		void CorrectedUnknownsVector(const double* v);

		/**
		*@brief Returns the minimum constraints
		*@param v the minimum constraints for each variable
		*/
		void MinimumUnknownsVector(double* v);

		/**
		*@brief Returns the maximum constraints
		*@param v the maximum constraints for each variable
		*/
		void MaximumUnknownsVector(double* v);

		/**
		*@brief Returns if an equation is differential or algebraic
		*@param v equal to 1 in case of differential equation, equal to 0 in case of algebraic equation
		*/
		void AlgebraicDifferentialVector(double* v);

		/**
		*@brief Returns the equations
		*@param t current time
		*@param y current solution
		*@param dy current residuals (in case of algebratic equations) or derivative with respect to time (in case of differential equations)
		*/
		void Equations(const double t, const double* y, double* dy);

		/**
		*@brief Prints the solution on a file
		*@param t current time
		*@param y current solution
		*@param fOutput file where to write the solution
		*/
		void Print(const double t, const Eigen::VectorXd& y, std::ofstream& fOutput);

		/**
		*@brief This function is needed by the Sundials IDA solver for recognizing the algebraic and differential equations
		*@param yp (TODO)
		*/
		void CorrectAlgebraicEquations(double* yp);

		/**
		*@brief This function is needed by the Sundials IDA solver for recognizing the algebraic and differential equations
		*@param upv (TODO)
		*@param resv (TODO)
		*/
		void CorrectDifferentialEquations(double* upv, double* resv);

		/**
		*@brief Sets simulation type
		*@param solver_type simulation type: SIMULATION_TYPE_Y, SIMULATION_TYPE_YTM, SIMULATION_TYPE_TM, SIMULATION_TYPE_YT
		*/
		void SetType(const Simulation_Type flag) { type_ = flag; SetAlgebraicDifferentialEquations(); }
		
		/**
		*@brief Returns the solver type (SOLVER_TYPE_BURNERSTABILIZED | SOLVER_TYPE_FLAMESPEED)
		*/
		Solver_Type solver_type() const { return solver_type_; }

		/**
		*@brief Returns the simulation type (SIMULATION_TYPE_Y, SIMULATION_TYPE_YTM, SIMULATION_TYPE_TM, SIMULATION_TYPE_YT)
		*/
		Simulation_Type type() const { return type_; }

		/**
		*@brief Returns if the senstivity analysis was enabled or not
		*/
		bool sensitivity_analysis() const { return sensitivity_analysis_; }

		/**
		*@brief Prints the solution on a file
		*@param name_file file where to write the solution
		*/
		void Print(const std::string name_file);

		/**
		*@brief Prints the soot solution on a file
		*@param output_folder name of the output folder where to write the solution
		*/
		void PrintSoot(const boost::filesystem::path output_folder);

		/**
		*@brief Prints the soot (calculated using the HMOM) solution on a file
		*@param output_folder name of the output folder where to write the solution
		*/
		void PrintHMOM(const boost::filesystem::path output_folder);

		/**
		*@brief Prints the post processing data
		*/
		void PrintOnTheFlyPostProcessing();

		/**
		*@brief Functions to be called at the end of every time step of DAE integrators
		*@param t current time
		*@param y current solution
		*/
		void Print(const double t, const double* y);

		/**
		*@brief Print on the screen (function which is called at the end of every time step of DAE integrators)
		*@param y current solution
		*@param norm_residuals norm of residuals
		*/
		void Print(const double* y, const double norm_residuals);

		/**
		*@brief Returns if equations are differential or algebraic
		*@param id_equations_ returns true if the equation is differential, false if it is algebraic
		*/
		const std::vector<bool>&	id_equations() const { return id_equations_;  }

		/**
		*@brief Returns the current mass fractions of species
		*/
		const std::vector<Eigen::VectorXd>&	Y() const { return Y_; }

		/**
		*@brief Returns the current temperatures
		*/
		const Eigen::VectorXd&	T() const { return T_; }

		/**
		*@brief Returns the ???
		*/
		const Eigen::VectorXd&	U() const { return U_; }

		/**
		*@brief Returns the ???
		*/
		const Eigen::VectorXd&	G() const { return G_; }

		/**
		*@brief Returns the ???
		*/
		const Eigen::VectorXd&	H() const { return H_; }

		/**
		*@brief Performs the sensitivity analysis
		*/
		void SensitivityAnalysis();

		/**
		*@brief Prints the solution on a XML file
		*@param file_name file where to write the solution in XML format
		*/
		void PrintXMLFile(const std::string file_name);

		/**
		*@brief Returns the output folder
		*/
		boost::filesystem::path output_folder() const { return output_folder_;  }

		/**
		*@brief Initialize the solver from a backup solution
		*@param name_file name of file corresponding to the backup solution
		*/
		void InitializeFromBackupFile(const boost::filesystem::path name_file);

		/**
		*@brief Calculates only the main diagonal of the Jacobian matrix
		*@param t current time
		*@param y current solution
		*@param J diagonal of the Jacobian matrix to be calculated
		*/
		void DiagonalJacobian(const double t, double* y, double* J);

		/**
		*@brief Corrects the diagonal of the Jacobian matrix to be used by the IDA solver
		*@param alfa parameter from the preconditioner
		*@param J diagonal of the Jacobian matrix to be corrected
		*/
		void DiagonalJacobianForIDA(const double alfa, double* J);

		/**
		*@brief Calculates and prints on the screen the norms of the current residuals
		*/
		void norm();

		/**
		*@brief Returns the sparsity pattern of a tridiagonal block matrix
		*@param rows row indices of non zero elements (1 index based)
		*@param cols column indices of non zero elements (1 index based)
		*/
		void SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

	private:

		/**
		*@brief Allocates the memory
		*/
		void MemoryAllocation();

		/**
		*@brief Calculates the thermodynamic and transport properties together with the formation rates of every species
		*/
		void Properties();

		/**
		*@brief Calculates the corrected diffusion fluxes
		*/
		void DiffusionFluxes();

		/**
		*@brief Calculates the local residence time
		*@param tau the calculated local residence time [s]
		        the residence time is positive if calculated from the fuel side,
				negative if calculated from the oxidizer side
		*/
		void ResidenceTime(Eigen::VectorXd& tau);

		/**
		*@brief Sets the algebraic and differential equations
		*/
		void SetAlgebraicDifferentialEquations();

		/**
		*@brief Returns the equations for momentum only
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_Momentum(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for mass fractions of species and momentum
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MassFractions_Momentum(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for mass fractions of species, temperature, and mass flow rate
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MassFractions_Temperature_Momentum(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for mass fractions of species and temperature
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MassFractions_Temperature(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for mass fractions of species
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MassFractions(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for mass flow rate
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_Temperature_MassFlowRate(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for HMOMs
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_HMOM(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for HMOMs
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MassFractions_HMOM(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for HMOMs
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MassFractions_Temperature_HMOM(const double t, const double* y, double* dy);

		/**
		*@brief Returns the equations for mixture fraction
		*@param t current time
		*@param y current solution
		*@param dy residuals (in case of algebraic equations) or time derivatives (in case of differential equations)
		*/
		void Equations_MixtureFraction(const double t, const double* y, double* dy);

		/**
		*@brief Builds the equations of species mass fractions
		*/
		void SubEquations_MassFractions();

		/**
		*@brief Builds the equations of temperature
		*/
		void SubEquations_Temperature();

		/**
		*@brief Builds the equations of momentum
		*/
		void SubEquations_Momentum();

		/**
		*@brief Builds the equations of HMOMs
		*/
		void SubEquations_HMOM();

		/**
		*@brief Builds the mixture fraction equations
		*/
		void SubEquations_MixtureFraction();

		/**
		*@brief Updates the current solution from the provided vector
		*@param y vector to be used to update the solution
		*/
		void Recover_Unknowns(const double* y);

		/**
		*@brief Returns the residuals (or time derivatives)
		*@param dy the residuals or time-derivatives
		*/
		void Recover_Residuals(double* dy);

		/**
		*@brief Calculates and returns the Jacobian matrix
		*@param J the Jacobian matrix
		*/
		void Jacobian(OpenSMOKE::OpenSMOKEBandMatrixDouble* J);

		/**
		*@brief Calculates the soot deposition rate (in kg/m2/s) along the right wall
		*/
		void SootDepositionOnTheWall();

		/**
		*@brief Updates the boundary conditions in case of dynamic simulations
		*@param the current time (in s)
		*/
		void UpdateBoundaries(const double t);

		/**
		*@brief Calculates the dynamic viscosity according to the (simplified) Sutherland's law
		*@param T temperature (in K)
		*@return the dynamic viscosity (in kg/m/s)
		*/
		double SutherlandViscosity(const double T);

	private:	// solutions

		/**
		*@brief Calculates the initial solution for a counterflow diffusion flame problem
		*@param dae_parameter parameters governing the solution of DAE systems
		*@param nls_parameters parameters governing the solution of NL systems
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int InitialSolutionCounterFlowDiffusionFlame(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);

		/**
		*@brief Calculates the whole solution
		*@param dae_parameter parameters governing the solution of DAE systems
		*@param nls_parameters parameters governing the solution of NL systems
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int CompleteSolution(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);

		/**
		*@brief Calculates the solution for a counterflow diffusion flame problem
		*@param dae_parameter parameters governing the solution of DAE systems
		*@param nls_parameters parameters governing the solution of NL systems
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int SolveCounterFlowDiffusionFlame(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
			NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
			NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);

		/**
		*@brief Calculates the solution for a counterflow diffusion flame problem with dynamic boundaries
		*@param dae_parameter parameters governing the solution of DAE systems
		*/
		int SolveDynamicBoundariesCounterFlowDiffusionFlame(DaeSMOKE::DaeSolver_Parameters& dae_parameters);


	private:

		/**
		*@brief Solves a single DAE system
		*@param dae_parameter parameters governing the solution of DAE systems
		*/
		int SolveInitialDAE(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double tEnd);

		/**
		*@brief Solves a single DAE system
		*@param dae_parameter parameters governing the solution of DAE systems
		*/
		int SolveDAE(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double tEnd, const double tStart = 0.);

		/**
		*@brief Solves a single NL system
		*@param nls_parameters parameters governing the solution of NL systems
		*/
		int SolveNLS(NlsSMOKE::NonLinearSolver_Parameters& nls_parameters);

		/**
		*@brief Solves a single false transient
		*@param false_transient_parameters parameters governing the solution of false-transients
		*/
		int SolveFalseTransient(NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters);

	private:	// additional checks

		/**
		*@brief Performs some tests to check if the grid is sufficiently fine to have a real adiabatic simulation
		*/
		void CheckForAdiabaticity();

		/**
		*@brief Performs some tests to check if the conservation of atomic elements is ensured
		*/
		void AtomicAnalysis();

		/**
		*@brief Quantifies the backdiffusion at the inlet section
		*/
		void CheckForInlet();

		/**
		*@brief Reconstructs the mixture fraction using the Bilger's formula
		*@param i grid point
		*/
		double ReconstructMixtureFraction(const int i);

	private:	// grid refinement

		/**
		*@brief Refines the 1D grid
		*@param count an index tracking the total number of refinements
		*/
		OpenSMOKE::Adapter_Grid1D_Status RefineGrid(const unsigned int count);

		/**
		*@brief Refines the 1D grid locally
		*@param xA defines the interval where the local refinement has to be applied
		*@param xB defines the interval where the local refinement has to be applied
		*@param count an index tracking the total number of refinements
		*/
		OpenSMOKE::Adapter_Grid1D_Status RefineGrid(const double xA, const double xB, const unsigned int count);

		/**
		*@brief Refines the 1D grid by doubling the number of intervals
		*@param count an index tracking the total number of refinements
		*/
		OpenSMOKE::Adapter_Grid1D_Status Doubling(const unsigned int count);

		/**
		*@brief Adapts the current grid (points are only redistributed, but no additions are performed)
		*/
		OpenSMOKE::Adapter_Grid1D_Status Regrid();

	private:

		// References
		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;	//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;		//!< reference to the kinetic map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN&		transportMap_;		//!< reference to the trasport properties map
		OpenSMOKE::Grid1D&								grid_;				//!< reference to the 1D grid

		// Main variables
		double							P_;			//!< pressure [Pa]
		Eigen::VectorXd					T_;			//!< temperature [K]
		Eigen::VectorXd					U_;			//!< scaled momentum [kg/m2/s]
		Eigen::VectorXd					G_;			//!< 
		Eigen::VectorXd					H_;			//!< 
		std::vector<Eigen::VectorXd>	Y_;			//!< mass fractions of species
		std::vector<Eigen::VectorXd>	X_;			//!< mole fractions of species

		// Auxiliary variables
		Eigen::VectorXd					M_;				//!< 
		Eigen::VectorXd					G_over_rho_;	//!< 

		// Mixture properties
		Eigen::VectorXd					rho_;				//!< density [kg/m3]
		Eigen::VectorXd					mw_;				//!< molecular weight [kg/kmol]
		Eigen::VectorXd					cp_;				//!< constant pressure specific heat [J/kg/K]
		std::vector<Eigen::VectorXd>	cp_species_;		//!< constant pressure specific heats of single species [J/kg/K]
		Eigen::VectorXd					lambda_;			//!< thermal conductivity [W/m/K]
		Eigen::VectorXd					mu_;				//!< dynamic viscosity [kg/m/s]
		std::vector<Eigen::VectorXd>	gamma_fick_;		//!< ordinary (Fick) mass diffusion coefficients [m2/s]
		std::vector<Eigen::VectorXd>	gamma_fick_star_;	//!< corrected ordinary (Fick) mass diffusion coefficients [m2/s]
		std::vector<Eigen::VectorXd>	gamma_soret_star_;	//!< corrected Soret mass diffusion coefficients[m2/s]
		Eigen::VectorXd					Q_;					//!< heat released by reactions [W/m3]
		std::vector<Eigen::VectorXd>	omega_;				//!<formation rates of species [kg/m3/s]
		
		// Mass diffusion fluxes
		std::vector<Eigen::VectorXd>	j_star_;					//!< corrected total mass diffusion fluxes [kg/m2/s]
		std::vector<Eigen::VectorXd>	j_fick_star_;				//!< corrected fick mass diffusion fluxes [kg/m2/s]
		std::vector<Eigen::VectorXd>	j_soret_star_;				//!< corrected soret mass diffusion fluxes [kg/m2/s]
		std::vector<Eigen::VectorXd>	j_thermophoretic_star_;		//!< corrected thermophoretic mass diffusion fluxes [kg/m2/s]
		Eigen::VectorXd					jc_star_;					//!< correction flux [kg/m2/s]

		// Spatial derivatives
		Eigen::VectorXd					dT_over_dx_;				//!< spatial derivative of temperature [K/m]
		Eigen::VectorXd					dT_over_dx_centered_;		//!< spatial derivative of temperature (centered) [K/m]
		Eigen::VectorXd					lambda_d2T_over_dx2_;		//!< 2nd order spatial derivative of (lambda*T) 
		std::vector<Eigen::VectorXd>	dX_over_dx_;				//!< spatial derivative of mole fractions
		std::vector<Eigen::VectorXd>	dY_over_dx_;				//!< spatial derivative of mass fractions
		Eigen::VectorXd					dM_over_dx_;				//!< spatial derivative of ???
		Eigen::VectorXd					dU_over_dx_;				//!< spatial derivative of ???
		Eigen::VectorXd					mu_d2G_over_rho_over_dx2_;	//!< 2nd order spatial derivative of ???

		// Time derivatives
		std::vector<Eigen::VectorXd>	dY_over_dt_;				//!< time derivatives of mass fractions
		Eigen::VectorXd					dT_over_dt_;				//!< time derivative of temperature
		Eigen::VectorXd					dU_over_dt_;				//!< time derivative of ???
		Eigen::VectorXd					dG_over_dt_;				//!< time derivative of ???
		Eigen::VectorXd					dH_over_dt_;				//!< time derivative of ???

		// Fuel side
		Eigen::VectorXd					Y_fuel_;					//!< mass fractions of fuel mixture [-]
		double							T_fuel_;					//!< temperature of fuel mixture [K]
		double							U_fuel_;					//!< ???
		double							radial_gradient_fuel_;		//!< radial gradient of fuel mixture [1/s]
		double							rho_fuel_;					//!< ???
		double							v_fuel_;					//!< ???

		// Oxidizer
		Eigen::VectorXd					Y_oxidizer_;				//!< mass fractions of oxidizer mixture [-]
		double							T_oxidizer_;				//!< temperature of oxidizer mixture [K]
		double							U_oxidizer_;				//!< ???
		double							radial_gradient_oxidizer_;	//!< radial gradient of oxidizer mixture [1/s]
		double							rho_oxidizer_;				//!< ???
		double							v_oxidizer_;				//!< ???

		// Algebraic/Differential equations
		Eigen::VectorXi		differential_equations_;	//!< list of differential equations
		Eigen::VectorXi		algebraic_equations_;		//!< list of algebraic equations
		std::vector<bool>	id_equations_;				//!< identifier of differential (true) and algebraic (false) equations

		// Solvers
		Solver_Type solver_type_;	//!< solver type: SOLVER_TYPE_BURNERSTABILIZED | SOLVER_TYPE_FLAMESPEED
		Simulation_Type type_;		//!< simulation type: SIMULATION_TYPE_Y, SIMULATION_TYPE_YTM, SIMULATION_TYPE_TM, SIMULATION_TYPE_YT
		
		// Sensitivity analysis
		bool sensitivity_analysis_;							//!< sensitivity analysis on/off
		Eigen::VectorXi indices_of_sensitivity_species_;	//!< list of species for which the sensitivity coefficients will be written on a file

		// Soot
		bool is_polimi_soot_;									//!< true if the soot analyzer for soot is available
		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_analyzer_;	//!< pointer to the soot analyzer

		// Virtual chemistry
		bool is_virtual_chemistry_;								//!< true if the soot analyzer for virtual chemistry
		OpenSMOKE::VirtualChemistry* virtual_chemistry_;		//!< pointer to the virtual chemistry
		bool is_simplified_transport_properties_;				//!< simplfied transport properties (thermal conductivity and viscosity)


		// OnTheFlyPostProcessing
		bool is_on_the_fly_post_processing_;							//!< true if the post processing (on the fly) is turned on
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing_;	//!< pointer to the post processing (on the fly)

		// Dynamic boundaries
		bool is_dynamic_boundaries_;						//!< true if the dynamic boundaries are turned on
		OpenSMOKE::DynamicBoundaries* dynamic_boundaries_;	//!< pointer to the dynamic boundary object


		// Initial/Fixed temperature profile
		bool is_userdefined_fixed_temperature_profile_;			//!< true if a fixed, user-defined temperature profile has to be used
		bool is_userdefined_initial_temperature_profile_;		//!< true if a initial, user-defined temperature profile has to be used
		FixedProfile* userdefined_temperature_profile_;			//!< userdefined temperature profile
		
		// Radiative heat transfer
		bool radiative_heat_transfer_;					//!< radiative heat transfer on/off
		double environment_temperature_;				//!< temperature of external environment in K
		Eigen::VectorXd	Q_radiation_;					//!< radiative heat [W/m3]
		Eigen::VectorXd	planck_mean_absorption_gas_;	//!< planck mean absorption coefficient gas phase [1/m]
		Eigen::VectorXd	planck_mean_absorption_soot_;	//!< planck mean absorption coefficient soot [1/m]

		// Additional variables
		bool soret_effect_;				//!< soret effect on/off
		bool use_dae_solver_;			//!< true if the transient solutions have to be calculated (more robust approach)
		bool use_nls_solver_;			//!< true if after the DAE solution, the NLS solution is attempted (more accurate approach)
		double timeFixedTemperature_;	//!< time for solving the DAE system in case of fixed temperature (default 1 s)
		double timeFixedComposition_;	//!< time for solving the DAE system in case of fixed composition (default 1 s)
		double timeComplete_;			//!< time for solving the DAE system in case of complete set of eqs (default 1 s)

		// Output
		unsigned int n_steps_video_;				//!< number of steps for updating info on the screen
		unsigned int count_video_;					//!< counter of steps for updating info on the screen
		unsigned int iterations_;					//!< number of iterations
		boost::filesystem::path output_folder_;		//!< name of output folder

		// Auxiliary vectors
		OpenSMOKE::OpenSMOKEVectorDouble	aux_Y;				//!< vector containing the mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble	aux_X;				//!< vector containing the mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble	aux_C;				//!< vector containing the concentration of gaseous species [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble	aux_R;				//!< vector containing the formation rates of gaseous species [kg/m3/s]
		OpenSMOKE::OpenSMOKEVectorDouble    aux_prov;			//!< auxiliary vector

		//HMOM
		OpenSMOKE::HMOM*				hmom_;					//!< HMOM object
		bool is_hmom_soot_;										//!< return true is the HMOM is turned on
		std::vector<Eigen::VectorXd>	dhmom_M_over_dx_;		//!< spatial derivative of normalized HMOMs
		std::vector<Eigen::VectorXd>	dhmom_M_over_dt_;		//!< time derivative of normalized HMOMs
		std::vector<Eigen::VectorXd>	hmom_M_;				//!< normalized HMOMs
		std::vector<Eigen::VectorXd>	d2hmom_M_over_dx2_;		//!< 2nd order spatial derivative of normalized HMOMs
		std::vector<Eigen::VectorXd>	hmom_thermophoresis_;	//!< thermophoretic contribution in HMOM [mol.kg/m6/s] 

		// Mixture fraction
		Eigen::VectorXd rho_times_alpha_d2csi_over_dx2_;	//!< spatial 2nd order derivative of mixture fraction
		Eigen::VectorXd	dcsi_over_dx_;						//!< spatial derivative of mixture fraction
		Eigen::VectorXd	dcsi_over_dt_;						//!< time derivative of mixture fraction
		Eigen::VectorXd	csi_;								//!< mixture fraction

		// To be used only if the BzzMath libraries are available
		#if OPENSMOKE_USE_BZZMATH == 1
		BzzDaeSparseObject dae_object_;
		BzzNonLinearSystemSparseObject nls_object_;
		#endif

		double n_geometry_;
		Eigen::VectorXd boundary_fuel_mass_fractions_;
		Eigen::VectorXd boundary_oxidizer_mass_fractions_;
		
		Eigen::VectorXd Y_peak_;
		double T_peak_;

		double x_peak_;
		double width_mixing_;
		double H_starting_guess_;
		InitialProfile_Type initial_profile_type_;

		OpenSMOKE::derivative_type derivative_type_species_;
		OpenSMOKE::derivative_type derivative_type_temperature_;

		double gas_reaction_rate_multiplier_;

		MassDiffusionCoefficients_Type mass_diffusion_coefficients_type_;
		std::vector<double> lewis_numbers_;

		bool is_deposition_wall_;
		double soot_deposition_;
	};
}

#include "OpenSMOKE_CounterFlowFlame1D.hpp"

#endif /* OpenSMOKE_CounterFlowFlame1D_H */

