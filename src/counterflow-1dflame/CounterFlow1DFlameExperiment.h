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

#ifndef OpenSMOKE_CounterFlow1DFlameExperiment_H
#define OpenSMOKE_CounterFlow1DFlameExperiment_H

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/Utilities.h"
#include "idealreactors/utilities/Grammar_LewisNumbers.h"

// 1D grid
#include "utilities/grids/adaptive/Grid1D.h"
#include "utilities/grids/adaptive/Grammar_Grid1D.h"
#include "utilities/grids/adaptive/Adapter_Grid1D.h"

// Hybrid Method of Moments
#include "utilities/soot/hmom/HMOM.h"

#include "OptimizationRules_CounterFlow1DFlame.h"
#include "Grammar_CounterFlow1DFlameExperiment.h"
#include "OpenSMOKE_CounterFlowFlame1D.h"


namespace OpenSMOKE
{
	class CounterFlow1DFlameExperiment
	{
	public:

		void Setup(const std::string input_file_name,
			OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML,
			OpenSMOKE::KineticsMap_CHEMKIN*				kineticsMapXML,
			OpenSMOKE::TransportPropertiesMap_CHEMKIN*	transportMapXML,
			OpenSMOKE::VirtualChemistry*				virtual_chemistry);

		void Solve(const bool verbose = false);

		double norm2_abs_error() const { return norm2_abs_error_; }
		double norm2_rel_error() const { return norm2_rel_error_; }

		const OpenSMOKE::OptimizationRules_CounterFlow1DFlame* optimization() const { return optimization_; }

	private:

		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*				kineticsMapXML_;
		OpenSMOKE::TransportPropertiesMap_CHEMKIN*	transportMapXML_;

		OpenSMOKE::OptimizationRules_CounterFlow1DFlame*	optimization_;
		DaeSMOKE::DaeSolver_Parameters*					dae_parameters;
		NlsSMOKE::NonLinearSolver_Parameters*			nls_parameters;
		NlsSMOKE::FalseTransientSolver_Parameters*		false_transient_parameters;

		OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
		OpenSMOKE::Grid1D* grid;
		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
		OpenSMOKE::HMOM* hmom;
		OpenSMOKE::DynamicBoundaries* dynamic_boundaries;

		OpenSMOKE::VirtualChemistry*			virtual_chemistry_;

		double end_value_;

		double norm2_abs_error_;
		double norm2_rel_error_;
	};
}

#include "CounterFlow1DFlameExperiment.hpp"

#endif // OpenSMOKE_CounterFlow1DFlameExperiment_H