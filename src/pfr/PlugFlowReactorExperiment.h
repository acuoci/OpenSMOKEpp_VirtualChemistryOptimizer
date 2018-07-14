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

#ifndef OpenSMOKE_PlugFlowReactorExperiment_H
#define OpenSMOKE_PlugFlowReactorExperiment_H

#include "OptimizationRules_PlugFlowReactor.h"
#include "Grammar_PlugFlowReactorExperiment.h"
#include "idealreactors/plugflow/PlugFlowReactor"

namespace OpenSMOKE
{
	class PlugFlowReactorExperiment
	{
	public:

		void Setup(const std::string input_file_name,
			OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
			OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML,
			OpenSMOKE::VirtualChemistry*			virtual_chemistry);

		void Solve(const bool verbose = false);
		double Solve(const double x);

		double norm2_abs_error() const { return norm2_abs_error_; }
		double norm2_rel_error() const { return norm2_rel_error_; }

		const Eigen::VectorXd& y() const { return y_; }

		const OpenSMOKE::OptimizationRules_PlugFlowReactor* optimization() const { return optimization_; }

	private:

		OpenSMOKE::PlugFlowReactor_Isothermal* plugflow_;

		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML_;

		OpenSMOKE::OptimizationRules_PlugFlowReactor*	optimization_;
		OpenSMOKE::PolimiSoot_Analyzer*					polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing*				on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA*						onTheFlyROPA_;
		OpenSMOKE::PlugFlowReactor_Options*				plugflow_options_;
		OpenSMOKE::ODE_Parameters*						ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options*			sensitivity_options_;
		OpenSMOKE::VirtualChemistry*					virtual_chemistry_;
		OpenSMOKE::IgnitionDelayTimes_Analyzer*			idt;

		double end_value_;

		double norm2_abs_error_;
		double norm2_rel_error_;
		Eigen::VectorXd y_;

	};
}

#include "PlugFlowReactorExperiment.hpp"

#endif // OpenSMOKE_PlugFlowReactorExperiment_H

