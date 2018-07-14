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

#ifndef OpenSMOKE_Grammar_XYProfile_H
#define	OpenSMOKE_Grammar_XYProfile_H

#include <string>
#include "boost/filesystem.hpp"
#include "math/OpenSMOKEVector.h"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_XYProfile : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XUnits", 
																OpenSMOKE::SINGLE_STRING, 
																"Units of x variable", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YUnits", 
																OpenSMOKE::SINGLE_STRING, 
																"Units of y variable", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XVariable", 
																OpenSMOKE::SINGLE_STRING, 
																"X variable", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YVariable", 
																OpenSMOKE::SINGLE_STRING, 
																"X variable", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Profile", 
																OpenSMOKE::VECTOR_DOUBLE, 
																"X,Y profile", 
																true) );										
		}
	};

	void GetXYProfileFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, 
									OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y, 
									std::string& x_variable, std::string& y_variable)
	{
		Grammar_XYProfile grammar_xyprofile_status;
		dictionary.SetGrammar(grammar_xyprofile_status);
		
		dictionary.ReadString("@XVariable", x_variable);
		dictionary.ReadString("@YVariable", y_variable);

		std::string x_units;
		dictionary.ReadString("@XUnits", x_units);
		std::string y_units;
		dictionary.ReadString("@YUnits", y_units);

		std::vector<double> xy_profile;
		dictionary.ReadOption("@Profile", xy_profile);

		double conversion_x = 1.;
		if (x_variable == "length")
		{
			if (x_units == "m")				conversion_x = 1.;
			else if (x_units == "cm")		conversion_x = 0.01;
			else if (x_units == "mm")		conversion_x = 0.001;
			else OpenSMOKE::FatalErrorMessage("Unknown length unit");
		}
		else if (x_variable == "time")
		{
			if (x_units == "s")				conversion_x = 1.;
			else if (x_units == "min")		conversion_x = 60.;
			else if (x_units == "ms")		conversion_x = 0.001;
			else OpenSMOKE::FatalErrorMessage("Unknown time unit");
		}
		else OpenSMOKE::FatalErrorMessage("Unknown x variable for xy profile");

		double conversion_y = 1.;
		if (y_variable == "temperature")
		{
			if (y_units == "K")				conversion_y = 1.;
			else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
		}
		else if (y_variable == "velocity")
		{
			if (y_units == "m/s")			conversion_y = 1.;
			else if (y_units == "cm/s")		conversion_y = 0.01;
			else if (y_units == "mm/s")		conversion_y = 0.001;
			else OpenSMOKE::FatalErrorMessage("Unknown velocity unit");
		}
		else if (y_variable == "pressure")
		{
			if (y_units == "Pa")			conversion_y = 1.;
			else if (y_units == "atm")		conversion_y = 101325.;
			else if (y_units == "bar")		conversion_y = 100000.;
			else OpenSMOKE::FatalErrorMessage("Unknown pressure units");
		}
		else if (y_variable == "volume")
		{
			if (y_units == "m3")			conversion_y = 1.;
			else if (y_units == "cm3")		conversion_y = 1e-6;
			else if (y_units == "mm3")		conversion_y = 1e-9;
			else OpenSMOKE::FatalErrorMessage("Unknown volume units");
		}
		else if (y_variable == "specific-mass-flow-rate")
		{
			if (y_units == "kg/m2/s")		conversion_y = 1.;
			else if (y_units == "g/cm2/s")	conversion_y = 10.;
			else OpenSMOKE::FatalErrorMessage("Unknown specific-mass-flow-rate units");
		}
		else if (y_variable == "area-per-unit-of-volume")
		{
			if (y_units == "1/m")			conversion_y = 1.;
			else if (y_units == "1/cm")		conversion_y = 1.e2;
			else if (y_units == "1/mm")		conversion_y = 1.e3;
			else OpenSMOKE::FatalErrorMessage("Unknown area-per-unit-of-volume units");
		}
		else if (y_variable == "dimensionless")
		{
			if (y_units == "dimensionless")	conversion_y = 1.;
			else OpenSMOKE::FatalErrorMessage("Unknown dimensionless units");
		}
		
		else OpenSMOKE::FatalErrorMessage("Unknown y variable for xy profile");

		const int n = boost::lexical_cast<int>(xy_profile.size()) / 2;
		ChangeDimensions(n, &x, true);
		ChangeDimensions(n, &y, true);
		
		unsigned int count = 0;
		for(int i=1;i<=n;i++)
		{
			x[i] = xy_profile[count++] * conversion_x;
			y[i] = xy_profile[count++] * conversion_y;
		}
	}
}

#endif	/* OpenSMOKE_Grammar_XYProfile_H */

