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

#ifndef OpenSMOKE_Grammar_RapidKineticMechanism_H
#define	OpenSMOKE_Grammar_RapidKineticMechanism_H

// Thermodynamics
#include "kernel/thermo/Species.h"
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/thermo/ThermoReader.h"
#include "kernel/thermo/ThermoReaderPolicy_CHEMKIN.h"

// Transport
#include "kernel/transport/TransportPolicy_CHEMKIN.h"
#include "kernel/transport/TransportReader.h"
#include "kernel/transport/TransportReaderPolicy_CHEMKIN.h"

// Kinetics
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"

// Preprocessing
#include "preprocessing/PreProcessorSpecies.h"
#include "preprocessing/PreProcessorKinetics.h"
#include "preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h"

#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/TransportPropertiesMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"

// External libraries
#include <boost/program_options.hpp>
#include "rapidxml.hpp"

#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_RapidKineticMechanism : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Thermodynamics", 
																OpenSMOKE::SINGLE_PATH, 
																"Name of the file containing the thermodynamic data (CHEMKIN format)", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Kinetics", 
																OpenSMOKE::SINGLE_PATH, 
																"Name of the file containing the kinetic mechanism (CHEMKIN format)", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Transport", 
																OpenSMOKE::SINGLE_PATH, 
																"Name of the file containing the transport data (CHEMKIN format)", 
																false) );	

                        
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Output", 
																OpenSMOKE::SINGLE_PATH, 
																"Name of the folder where to write the interpreted kinetic mechanism (XML format)", 
																true) );	
		}
	};

	void RapidKineticMechanismWithoutTransport(	boost::filesystem::path path_output, 
                                                        const boost::filesystem::path thermo_file,
                                                        const boost::filesystem::path kinetics_file,
                                                        const boost::filesystem::path surface_file = "",
                                                        const boost::filesystem::path solid_file = "") 
	{
		if ( !boost::filesystem::exists(thermo_file) )
			OpenSMOKE::FatalErrorMessage("The " + thermo_file.string() + " file does not exist!");

		if ( !boost::filesystem::exists(kinetics_file) )
			OpenSMOKE::FatalErrorMessage("The " + kinetics_file.string() + "file does not exist!");

		if ( !boost::filesystem::exists(path_output) )
			OpenSMOKE::CreateDirectory(path_output);

		typedef OpenSMOKE::Species< OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN > SpeciesCHEMKIN;
		typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN>  > PreProcessorSpecies_CHEMKIN_WithoutTransport;
		typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> > PreProcessorKinetics_CHEMKIN;
		typedef OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > > ThermoReader_CHEMKIN;

		// Log file (Setup)
		boost::filesystem::path file_ascii_log_ = path_output / "log";
		std::ofstream fLog;
		fLog.open(std::string(file_ascii_log_.string()).c_str(), std::ios::out);
		OpenSMOKE::CheckIfFileIsOpen(fLog, file_ascii_log_.string());
		fLog.setf(std::ios::scientific);

		// Readers
		ThermoReader_CHEMKIN* thermoreader;
                
                // Preprocessing the kinetic mechanism
		PreProcessorKinetics_CHEMKIN preprocessor_kinetics(fLog);
		CheckForFatalError( preprocessor_kinetics.ReadFromASCIIFile(kinetics_file.string()) );
                
                // Reading thermodynamic database
		thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
		CheckForFatalError( thermoreader->ReadFromFile(thermo_file.string()) );

		// Preprocessing the thermodynamics
		PreProcessorSpecies_CHEMKIN_WithoutTransport* preprocessor_species_without_transport;
		preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, preprocessor_kinetics, fLog);
		CheckForFatalError( preprocessor_species_without_transport->Setup() );

		// Read kinetics from file
		CheckForFatalError( preprocessor_kinetics.ReadKineticsFromASCIIFile( preprocessor_species_without_transport->AtomicTable()) );

		// Write interpreted kinetics in new style
		{
			std::string author_name = "not provided";
			std::string place_name = "not provided";
			std::string comments = "no comments";

			std::stringstream xml_string;
			xml_string << std::setprecision(8);
			xml_string.setf(std::ios::scientific);

			xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

			xml_string << "<Properties>" << std::endl;
			xml_string << "  <Author>" << author_name <<"</Author>" << std::endl;
			xml_string << "  <Place>" << place_name <<"</Place>" << std::endl;
			xml_string << "  <Date>" << OpenSMOKE::GetCurrentDate() <<"</Date>" << std::endl;
			xml_string << "  <Time>" << OpenSMOKE::GetCurrentTime() <<"</Time>" << std::endl;
			xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" <<"</Comments>" << std::endl;
			xml_string << "</Properties>" << std::endl;

			// Thermodynamics and transport properties
			preprocessor_species_without_transport->WriteXMLFile(xml_string);
				
			// Kinetic mechanism
			preprocessor_kinetics.WriteXMLFile(xml_string);

			xml_string << "</opensmoke>" << std::endl;

			// Write file
			boost::filesystem::path kinetics_xml   = path_output / "kinetics.xml";
			std::ofstream fOutput;
			fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);
			fOutput << xml_string.str();
			fOutput.close();
		}

		// Write the reaction strings
		{
			OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;
			OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;

			// Read kinetic mechanism from file
			{
				rapidxml::xml_document<> doc;
				std::vector<char> xml_string;
				OpenSMOKE::OpenInputFileXML(doc, xml_string, path_output / "kinetics.xml");

				thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc); 
				kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc);
			}

			// Write the reaction names
			{
				std::stringstream xml_string;
				xml_string << std::setprecision(8);
				xml_string.setf(std::ios::scientific);

				xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				xml_string << "<opensmoke version=\"0.1a\">" << std::endl;
				xml_string << "<reaction-number>" << std::endl;
				xml_string << kineticsMapXML->NumberOfReactions() << std::endl;
				xml_string << "</reaction-number>" << std::endl;
				xml_string << "<reaction-names>" << std::endl;
				for (unsigned int k=0;k<kineticsMapXML->NumberOfReactions();k++)
				{
					std::string reaction_string;
					preprocessor_kinetics.reactions()[k].GetReactionString(thermodynamicsMapXML->NamesOfSpecies(), reaction_string);
					boost::erase_all(reaction_string , " ");
					xml_string << reaction_string << std::endl;;
				}
				xml_string << "</reaction-names>" << std::endl;
				xml_string << "</opensmoke>" << std::endl;

				boost::filesystem::path reaction_names_file = path_output / "reaction_names.xml";
				std::ofstream fOutput;
				fOutput.open(std::string(reaction_names_file.string()).c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				fOutput << xml_string.str();
				fOutput.close();
			}

			delete thermodynamicsMapXML; 
			delete kineticsMapXML;
		}
 
                delete thermoreader;
		delete preprocessor_species_without_transport;
	}

        void RapidKineticMechanismWithTransport(boost::filesystem::path path_output,
                                                const boost::filesystem::path transport_file,
                                                const boost::filesystem::path thermo_file,
                                                const boost::filesystem::path kinetics_file,
                                                const boost::filesystem::path surface_file = "",
                                                const boost::filesystem::path solid_file = "")
        {

            if (!boost::filesystem::exists(transport_file))
                OpenSMOKE::FatalErrorMessage("The " + transport_file.string() + " file does not exist!");

            if (!boost::filesystem::exists(thermo_file))
                OpenSMOKE::FatalErrorMessage("The " + thermo_file.string() + " file does not exist!");

            if (!boost::filesystem::exists(kinetics_file))
                OpenSMOKE::FatalErrorMessage("The " + kinetics_file.string() + "file does not exist!");

            if (!boost::filesystem::exists(path_output))
                OpenSMOKE::CreateDirectory(path_output);


            typedef OpenSMOKE::Species< OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN > SpeciesCHEMKIN;
            typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<SpeciesCHEMKIN> > PreProcessorSpecies_CHEMKIN_WithTransport;
            typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> > PreProcessorKinetics_CHEMKIN;
            typedef OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > > ThermoReader_CHEMKIN;
                
            // Log file (Setup)
            boost::filesystem::path file_ascii_log_ = path_output / "log";
            std::ofstream fLog;
            fLog.open(std::string(file_ascii_log_.string()).c_str(), std::ios::out);
            OpenSMOKE::CheckIfFileIsOpen(fLog, file_ascii_log_.string());
            fLog.setf(std::ios::scientific);

            // Readers
            ThermoReader_CHEMKIN* thermoreader;
            OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >* transportreader;

            // Preprocessing the kinetic mechanism
            PreProcessorKinetics_CHEMKIN preprocessor_kinetics(fLog);
            CheckForFatalError( preprocessor_kinetics.ReadFromASCIIFile(kinetics_file.string()) );
            
            // Reading thermodynamic database
            thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
            CheckForFatalError( thermoreader->ReadFromFile(thermo_file.string()) );

            // Reading transport database
            transportreader = new OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >;
            CheckForFatalError( transportreader->ReadFromFile(transport_file.string()) );    

            // Preprocessing the thermodynamic and transport properties
            PreProcessorSpecies_CHEMKIN_WithTransport* preprocessor_species_with_transport;
            preprocessor_species_with_transport = new PreProcessorSpecies_CHEMKIN_WithTransport(*thermoreader, *transportreader, preprocessor_kinetics, fLog);
            CheckForFatalError(preprocessor_species_with_transport->Setup());
            CheckForFatalError(preprocessor_species_with_transport->Fitting());

            // Read kinetics from file
            CheckForFatalError(preprocessor_kinetics.ReadKineticsFromASCIIFile(preprocessor_species_with_transport->AtomicTable()));

            // Write interpreted kinetics in new style
            {
                std::string author_name = "CRECK Modeling";
                std::string place_name = "Milan(Italy)";
                std::string comments = "no comments";

                std::stringstream xml_string;
                xml_string << std::setprecision(8);
                xml_string.setf(std::ios::scientific);

                xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
                xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

                xml_string << "<Properties>" << std::endl;
                xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
                xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
                xml_string << "  <Date>" << OpenSMOKE::GetCurrentDate() << "</Date>" << std::endl;
                xml_string << "  <Time>" << OpenSMOKE::GetCurrentTime() << "</Time>" << std::endl;
                xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "</Comments>" << std::endl;
                xml_string << "</Properties>" << std::endl;

                // Thermodynamics and transport properties
                preprocessor_species_with_transport->WriteXMLFile(xml_string);

                // Kinetic mechanism
                preprocessor_kinetics.WriteXMLFile(xml_string);

                xml_string << "</opensmoke>" << std::endl;

                // Write file
                boost::filesystem::path kinetics_xml = path_output / "kinetics.xml";
                std::ofstream fOutput;
                fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
                fOutput.setf(std::ios::scientific);
                fOutput << xml_string.str();
                fOutput.close();
            }

            // Write the reaction strings
            {
                OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;
                OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;

                // Read kinetic mechanism from file
                {
                    rapidxml::xml_document<> doc;
                    std::vector<char> xml_string;
                    OpenSMOKE::OpenInputFileXML(doc, xml_string, path_output / "kinetics.xml");

                    thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
                    kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc);
                }

                // Write the reaction names
                {
                    std::stringstream xml_string;
                    xml_string << std::setprecision(8);
                    xml_string.setf(std::ios::scientific);

                    xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
                    xml_string << "<opensmoke version=\"0.1a\">" << std::endl;
                    xml_string << "<reaction-number>" << std::endl;
                    xml_string << kineticsMapXML->NumberOfReactions() << std::endl;
                    xml_string << "</reaction-number>" << std::endl;
                    xml_string << "<reaction-names>" << std::endl;
                    for (unsigned int k = 0; k < kineticsMapXML->NumberOfReactions(); k++)
                    {
                        std::string reaction_string;
                        preprocessor_kinetics.reactions()[k].GetReactionString(thermodynamicsMapXML->NamesOfSpecies(), reaction_string);
                        boost::erase_all(reaction_string, " ");
                        xml_string << reaction_string << std::endl;
                    }
                    xml_string << "</reaction-names>" << std::endl;
                    xml_string << "</opensmoke>" << std::endl;

                    boost::filesystem::path reaction_names_file = path_output / "reaction_names.xml";
                    std::ofstream fOutput;
                    fOutput.open(std::string(reaction_names_file.string()).c_str(), std::ios::out);
                    fOutput.setf(std::ios::scientific);
                    fOutput << xml_string.str();
                    fOutput.close();
                }

                        delete thermodynamicsMapXML;
                delete kineticsMapXML;
            }   

            
            delete thermoreader;
            delete transportreader;
            delete preprocessor_species_with_transport;
        }
}


#endif	/* OpenSMOKE_Grammar_RapidKineticMechanism_H */

