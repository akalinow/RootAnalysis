//! Class for extracting data for further ML processing in python.
/*!
  \author Rafal Maselek
  \date May 2018
  
  This class extracts four-momenta of legs and jets from the analysis, together with some global parameters
  and particles' properties. Data is written to a TTree in "Summary" branch of the main analysis file. It
  can be then exported to python and TensorFlow.
*/

#include "MLAnalyzer.h"
#include <utility>
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/tokenizer.hpp"
#include "boost/functional/hash.hpp"

//! Constructor
/*!
*	Constructor that creates an objects and launches parser of config file. To edit the name of the config file
*   edit the variable "cfgFileName". If parser fails, the program execution continues, only warning is thrown.
*/
MLAnalyzer::MLAnalyzer(const std::string & aName, const std::string & aDecayMode = "None") : Analyzer(aName), covMET_(2,2)
{
	std::cout<<"[ML]\tCreation of MLAnalyzer object."<<std::endl;
	decayMode = aDecayMode;
	execution_counter_=0;
	MLTree_=NULL;
	cfgFileName = "ml_Properties.ini";
	try
	{
		parseCfg(cfgFileName);
	}
	catch(std::exception& e)
	{
		std::cerr<<"[ML][WARNING]\t PARTICLE PROPERTY PARSER FAILED! READ THE MESSAGE BELOW:"<<std::endl;
		std::cerr<<e.what()<<std::endl;
	}
}

//! Dummy destructor
MLAnalyzer::~MLAnalyzer()
{
	std::cout<<"[ML]\t MLAnalyzer finished working."<<std::endl;
}

//! Clone maker
/*
* Dummy function requested by multithread framework. Can be developed for some useful purposes.
*/
Analyzer* MLAnalyzer::clone() const 
{
	MLAnalyzer* clone = new MLAnalyzer("MLAnalyzer cloned", decayMode);
    return clone;
};

//! Clears the member field collections
void MLAnalyzer::clear()
{
	legs_p4_.clear();
	jets_p4_.clear();
	params_.clear();
	params_legs_.clear();
	params_jets_.clear();
}

//! Config file parser
/*!
*	Parser for config file containing properties of particles (legs and/or jets) to be extracted.
*	Parses uses boost library, it's the same mechanism as for loading other *.ini files in the analysis.
*	Throws some exceptions, which might be useful for debugging. However, in release version the exceptions
*	are catched (in constructor, this is where the function is called) and a cerr warning is generated, 
*	the execution continues. Fail of the parses should not kill the analysis.
*/
void MLAnalyzer::parseCfg(const std::string & cfgFileName)
{
  std::cout<<"[ML]\tParsing external file with particle properties selections."<<std::endl;
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(cfgFileName, pt);

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(", "); /*!< separator for strings */

  // loop over all possible properties
  for(auto it=PropertyEnumString::enumMap.begin(); it!=PropertyEnumString::enumMap.end(); ++it)
  {
  	std::string property = it->first;
  	boost::optional<std::string> str_opt = pt.get_optional<std::string>(std::string("MLAnalyzer.")+property);
  	if(!str_opt)
  	{
  		std::cerr<<"[ML][WARNING] PARAMETER: "<<property<<" NOT IN PARAM FILE!"<<std::endl;
  		continue;
  	}
  	std::string str = static_cast<std::string>(*str_opt);
  	tokenizer tokens(str, sep);
  	// parse tokens
  	for (auto it: tokens) 
	{
		
		if("l"==it.substr(0,1)) // if token starts with 'l'
		{
			auto iter = params_legs_.find(property);
			if (iter == params_legs_.end())
			{
				std::set<unsigned> s;
			    params_legs_[property]=s; // create new collection for given property
			}

			if(it.size()==1)
				params_legs_[property].insert(1000); // dummy value, will be removed after calling fixSelections()
			else
			{
				unsigned val = std::stoul(it.substr(1));
				params_legs_[property].insert(val);
			}
		}
		else if("j"==it.substr(0,1)) // if token starts with 'j'
		{
			auto iter = params_jets_.find(property);
			if (iter == params_jets_.end())
			{
				std::set<unsigned> s;
			    params_jets_[property]=s;
			}

			if(it.size()==1)
				params_jets_[property].insert(1000); // dummy value, will be removed after calling fixSelections()
			else
			{
				unsigned val = std::stoul(it.substr(1));
				params_jets_[property].insert(val);	// create new collection for given property
			}
		}
		else // unsupported argument
			throw std::invalid_argument(std::string("[ERROR] WRONG SELECTOR FOR PARTICLE PROPERTY! PROPERTY: ")+\
				property+std::string(" SELECTOR: ")+it);
	}
  }
  std::cout<<"[ML]\tParsing successful!"<<std::endl;
}

//! Fixes the selections parsed using parses
/*!
*	Parses accepts an option to extract a given property for all legs/jets. However, when parses is called,
*	the total number of legs/jets per event is unknown. Parses sets a dummy selection number to 1000.
*	This function remove that value and places values for all leg/jets instead.
*/
void MLAnalyzer::fixSelections()
{
	for (auto it=params_legs_.begin(); it!=params_legs_.end(); ++it)
	{
		auto search = (it->second).find(1000); // Check if there is a dummy value.
		if(search != (it->second).end())
		{
			(it->second).erase(search); // Remove dummy value
			for(unsigned ii=1; ii<=legs_no_; ii++)
			{
				(it->second).insert(ii);	// Insert indices of all legs
			}
		}
	}
	for (auto it=params_jets_.begin(); it!=params_jets_.end(); ++it)
	{
		auto search = (it->second).find(1000); // Check if there is a dummy value.
		if(search != (it->second).end())
		{
			(it->second).erase(search);	// Remove dummy value
			for(unsigned ii=1; ii<=jets_no_; ii++)
			{
				(it->second).insert(ii); // Insert indices of all jets
			}
		}
	}
}

//! Adds branches to a TTree in the output file of the analysis
/*! Adds branches to the TTree that will be created in the output file of the analysis.
*	Inherited interface deals with the details.
*	The function will be called twice, one after constructor (default behavior), another one by analyze().
*	It is required, because at the moment of creation MLAnalyzer does not know how many variables are to
*	be extracted to the tree. The if inside function switches different behavior for mentioned situations.
*/
void MLAnalyzer::addBranch(TTree *tree)
{

	if(!MLTree_)
	{
		// Will execute after creation of a MLAnalyzer object
		if(tree)
		{
			// Adding branches for global parameters and saving the address of TTree
			MLTree_ = tree;
			std::cout<<"[ML]\tAdding global parameter branches for ML analysis."<<std::endl;
			tree->Branch("visMass", &visMass_);
			tree->Branch("genMass", &genVisMass_);
			tree->Branch("higgsMassTrans", &higgsMassTrans_);
			tree->Branch("higgsPT", &higgsPT_);
			tree->Branch("covMET00", &covMET_[0][0]);
			tree->Branch("covMET01", &covMET_[0][1]);
			tree->Branch("covMET10", &covMET_[1][0]);
			tree->Branch("covMET11", &covMET_[1][1]);
		}
		else
			std::cerr<<"[ML][WARNING] No tree to write to!"<<std::endl;
	}
	else
	{
		// Will execute once during the first call of MLAnalyzer::analyze()
		std::cout<<"[ML]\tAdding more branches for ML analysis."<<std::endl;
		std::string leg_s ("leg_");
		std::string jet_s ("jet_");
		std::string endings[4]={"_E", "_pX", "_pY", "_pZ"};
		try
		{
			// create branches for four-momenta of legs
			for(unsigned ii=0; ii<legs_p4_.size();ii++)
			{
				std::string name = leg_s+std::to_string(ii/4+1)+endings[ii%4];
				tree->Branch(name.c_str(), &legs_p4_.at(ii));

			}
			// create branches for four-momenta of jets
			for(unsigned ii=0; ii<jets_p4_.size();ii++)
			{
				std::string name = jet_s+std::to_string(ii/4+1)+endings[ii%4];
				tree->Branch(name.c_str(), &jets_p4_.at(ii));
			}

			// create branches for particle properties for legs
			for (auto it=params_legs_.begin(); it!=params_legs_.end(); ++it)
			{
				std::string name = it->first;
				std::set<unsigned> selected_legs = it->second;
				for(auto sel=selected_legs.begin(); sel!=selected_legs.end(); ++sel)
				{
					unsigned no = *sel; /*!< number of selected leg */
					if(no > legs_no_ || no==0) // if the number doesn't match data from MLObjectMessenger
					{
						std::cerr<<"[ML][WARNING]\tRequest to add branch for params for leg with larger number than accessible (or zero). I will ignore the request."<<std::endl;
					}
					else
					{
						std::string full_name = leg_s+std::to_string(no)+std::string("_")+name;
						std::cout<<"[ML]\tAdding "<<full_name<<std::endl;
						tree->Branch(full_name.c_str(), &(params_.at(std::string("legs_")+name).at(no-1)));
					}

				}
			}
			// create branches for particle properties for legs
			for (auto it=params_jets_.begin(); it!=params_jets_.end(); ++it)
			{
				std::string name = it->first;
				std::set<unsigned> selected_jets = it->second;
				for(auto sel=selected_jets.begin(); sel!=selected_jets.end(); ++sel)
				{
					unsigned no = *sel; /*!< number of selected jet */
					if(no > jets_no_ || no==0) // if the number doesn't match data from MLObjectMessenger
					{
						std::cerr<<"[ML][WARNING]\tRequest to add branch for params for jet with larger number than accessible (or zero). I will ignore the request."<<std::endl;
					}
					else
					{
						std::string full_name = jet_s+std::to_string(no)+std::string("_")+name;
						tree->Branch(full_name.c_str(), &(params_.at(std::string("jets_")+name).at(no-1)));
					}

				}
			}
		}
		catch(const std::out_of_range& e)
		{
			std::throw_with_nested(std::runtime_error("[ERROR] OUT OF RANGE IN MLAnalyzer::addBranch!"));
		}
		catch(const std::exception& e)
		{
		     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLAnalyzer::addBranch!"));
		} 

	}
}

//! Populate member field collections with variables that will serve as bufors for saving data to TTree
/*! 
*	This function will check what bufors are needed for saving data to TTree and will populate member field
*	collections with apropriate variables. Deals with four-momenta and particle properties. Global variables
* 	have their own bufor, not in any collection. 
*/
void MLAnalyzer::prepareVariables(const std::vector<const HTTParticle*>* legs, const std::vector<const HTTParticle*>* jets)
{
	try
	{
		std::cout<<"[ML]\tCreating bufors for data extraction."<<std::endl;
		/*** FOUR-MOMENTA ***/
		for(unsigned leg=0; leg<legs->size();leg++)
			for(unsigned coord=0; coord<4; coord++)
			{
				double x = 0.0;
				legs_p4_.push_back(x); // add bufors to legs_p4_

			}
		for(unsigned jet=0; jet<jets->size();jet++)
			for(unsigned coord=0; coord<4; coord++)
			{
				double x = 0.0;
				jets_p4_.push_back(x); // add bufors to jets_p4_
			}
		/*** PARTICLE PROPERTIES ***/
		for (auto it=params_legs_.begin(); it!=params_legs_.end(); ++it)
		{
			std::string name = it->first;
			std::vector<double> v;
			for(unsigned leg=0; leg<legs->size(); leg++)
			{
				double x = 0.0;
				v.push_back(x);
			}
			params_[std::string("legs_")+name] = v; // add bufors to params_
		}
		for (auto it=params_jets_.begin(); it!=params_jets_.end(); ++it)
		{
			std::string name = it->first;
			std::vector<double> v;
			for(unsigned jet=0; jet<jets->size(); jet++)
			{
				double x = 0.0;
				v.push_back(x);
			}
			params_[std::string("jets_")+name] = v; // add bufors to params_
		}
	}
	catch(const std::out_of_range& e)
	{
		std::throw_with_nested(std::runtime_error("[ERROR] OUT OF RANGE IN MLAnalyzer::prepareVariables!"));
	}
	catch(const std::exception& e)
	{
	     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLAnalyzer::prepareVariables!"));
	} 

}

//! Reads particle four-momentum from TLorentzVector and assigns apropriate values to bufors in member field collections.
void MLAnalyzer::extractP4(const TLorentzVector& v, const std::string destination, const unsigned no)
{
	try
	{
		if(destination=="jets")
		{
			jets_p4_.at(4*no) = v.T();
			jets_p4_.at(4*no+1) = v.X();
			jets_p4_.at(4*no+2) = v.Y();
			jets_p4_.at(4*no+3) = v.Z();
		}
		else if(destination=="legs")
		{
			legs_p4_.at(4*no) = v.T();
			legs_p4_.at(4*no+1) = v.X();
			legs_p4_.at(4*no+2) = v.Y();
			legs_p4_.at(4*no+3) = v.Z();
		}
	}
	catch(const std::out_of_range& e)
	{
		std::throw_with_nested(std::runtime_error(std::string("[ERROR] OUT OF RANGE IN MLAnalyzer::extractP4!\nJETS_VECTOR_SIZE: ")\
		 + std::to_string(jets_p4_.size()) +  std::string(", LEGS_VECTOR_SIZE: ") + std::to_string(jets_p4_.size()) + \
		 std::string(", TRYING TO ACCESS: [") + std::to_string(4*no) + std::string(",") + std::to_string(4*no+3)+ std::string("]!")));
	}
	catch(const std::exception& e)
	{
	     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLAnalyzer::extractP4!"));
	} 
}

//! Reads particle properties and assigns apropriate values to bufors in member field collections.
void MLAnalyzer::extractParticleProperties( const std::vector<const HTTParticle*>* legs,  const std::vector<const HTTParticle*>* jets)
{
	try
	{
		std::string legs_s ("legs_");
		std::string jets_s ("jets_"); 
		for (auto it=params_legs_.begin(); it!=params_legs_.end(); ++it)
		{
			std::string name = it->first;
			for(auto iter = (it->second).begin(); iter != (it->second).end(); ++iter)
			{
				unsigned no = *iter;
				// assign property value to variable in params_
				(params_.at(legs_s+name)).at(no-1) = (legs->at(no-1))->getProperty(PropertyEnumString::enumMap.at(name));
			}
		}
	}
	catch(const std::out_of_range& e)
	{
		std::throw_with_nested(std::runtime_error(std::string("[ERROR] OUT OF RANGE IN MLAnalyzer::extractParticleProperties")));
	}
	catch(const std::exception& e)
	{
	     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLAnalyzer::extractParticleProperties!"));
	} 
		
}

//! Calculates and fill member fields (bufors) for quantities specific to HTT analysis -- global parameters.
void MLAnalyzer::globalsHTT(const MLObjectMessenger* mess, const std::vector<const HTTParticle*>* legs, const HTTAnalysis::sysEffects* aSystEffect)
{
	try
	{
		void* p = NULL;
		const HTTParticle* aMET = mess->getObject(static_cast<HTTParticle*>(p), std::string("met"));
		//const float* bs = nullptr;//mess->getObject(static_cast<float*>(p), "beta_score");
		const float* higgs_mass_trans = mess->getObject(static_cast<float*>(p), "higgs_mass_trans");
		covMET_[0][0] = *mess->getObject(static_cast<double*>(p), "covMET00");
		covMET_[0][1] = *mess->getObject(static_cast<double*>(p), "covMET01");
		covMET_[1][0] = *mess->getObject(static_cast<double*>(p), "covMET10");
		covMET_[1][1] = *mess->getObject(static_cast<double*>(p), "covMET11"); 
		
		//TODO: remove when covMET contains real values
		covMET_[0][0] = 1.;
		covMET_[1][1] = 1.;

		//const int* nJets = nullptr;//mess ->getObject(static_cast<int*>(p), "nJets30");

		const HTTParticle* leg1 = legs->at(0);
		const HTTParticle* leg2 = legs->at(1);

		//if(!(leg1 && leg2 && aMET && aSystEffect && bs && higgs_mass && nJets))
		//	throw std::logic_error("[ERROR] NULL POINTERS PRESENT!");
		// Calculation and assignement of global parameters
		const TLorentzVector & aVisSum = leg1->getP4() + leg2->getP4();
		genVisMass_ = *mess->getObject(static_cast<float*>(p),"gen_mass");
		visMass_ = aVisSum.M();
		higgsPT_ =  (aVisSum + aMET->getP4()).Pt();
		higgsMassTrans_ = *higgs_mass_trans;
	}
	catch(const std::out_of_range& e)
	{
		std::throw_with_nested(std::runtime_error("[ERROR] CALCULATING GLOBAL PARAMETERS FOR HTT REQUIRES AT LEAST 2 LEGS!"));
	}
	catch(const std::logic_error& e)
	{
		std::throw_with_nested(std::runtime_error("[ERROR] CANNOT FILL FIELDS IN MLAnalyzer::globalsHTT!"));
	}
	catch(const std::exception& e)
	{
	     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLAnalyzer::globalsHTT!"));
	} 
}

//! Analyzes the event -- extracts necessary data to a TTree
/*!
*	The function executes for each event and takes a pointer to MLObjectMessenger object. If the latter is not NULL,
*	it will proceed to extract data: legs and jets four-momenta, particle properties loaded from external file and
*	a few gloabl parameters for HTT analysis. One can easily get rid of the latter (or change to something else) by
*	commenting out the call of globalsHTT() function.	
*/
void MLAnalyzer::performAnalysis(const EventProxyBase& iEvent, ObjectMessenger *aMessenger)
{
	try
	{
		void* p = NULL; /*!< pointer that will be cast to apropriate type if needed */
		MLObjectMessenger* mess = ((MLObjectMessenger*)aMessenger); 
		const std::vector<const HTTParticle*>* legs = mess->getObjectVector(static_cast<HTTParticle*>(p), std::string("legs"));
		const std::vector<const HTTParticle*>* jets = mess->getObjectVector(static_cast<HTTParticle*>(p), std::string("jets"));
		const auto sE = HTTAnalysis::NOMINAL;
		const HTTAnalysis::sysEffects* aSystEffect = &sE;//mess->getObject(static_cast< HTTAnalysis::sysEffects*>(p), std::string("systEffect"));


		if(mess && legs && jets)
		{
			if(execution_counter_ == 0) // will execute only once at the beginning
			{
				std::cout<<"[ML]\tMLAnalyzer starts working."<<std::endl;
				// Number of leg and jet object must be constant! Objects can be dummy.
				legs_no_ = legs->size();
				jets_no_ = jets->size();
				// Knowing the number of jets and legs fix the selections to match it (remove 1000 dummy values)
				fixSelections();
				// Stuff to do at the first execution, e.g. assigning branches
				prepareVariables(legs, jets); // prepare bufors
				addBranch(MLTree_);	// add branches	
			}
			if(legs->size()>0)
				execution_counter_++; /*!< count the number of calls of this function */

			if(jets->size() != jets_no_ || legs->size() != legs_no_)
				throw std::logic_error("[ERROR] MLAnalyzer REQUIRES CONSTANT NUMBER OF JETS AND LEGS IN EACH CALL!");
			// variable "no" ensures assigning values to correct leg/jet; it has to be set to 0 before for_each
			unsigned no = 0;
			std::for_each(legs->begin(), legs->end(), [&](const HTTParticle* leg)
			{
				extractP4(leg->getP4(), std::string("legs"), no); // copy P4 to bufors
				no++;
			});
			no = 0;
			std::for_each(jets->begin(), jets->end(), [&](const HTTParticle* jet)
			{
				extractP4(jet->getP4(), std::string("jets"), no); // copy P4 to bufors
			});

			// Extract values of properties defined in PropertyEnum.h and selected in external file
			extractParticleProperties(legs, jets);
			// Dealing with global parameters, one can put here another function for other analysis
			globalsHTT(const_cast<const MLObjectMessenger*>(mess), legs, aSystEffect);
			
			// Cleaning the content of ObjectMessenger NO DELETE IS CALLED!!!
			mess->clear();
		}
	}
	catch(const std::exception& e)
	{
		std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLAnalyzer::performAnalysis!"));
	}
}

//! The main function of the class, analyzes the event -- extracts necessary data to a TTree
/*!
*	Executes performAnalysis function.
*/
bool MLAnalyzer::analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger)
{
	if(MLTree_ and aMessenger and std::string("MLObjectMessenger").compare((aMessenger->name()).substr(0,17))==0) // if NULL it will do nothing
	{
		performAnalysis(iEvent, aMessenger);
	}
	return true;
}

//! Dummy function
bool MLAnalyzer::analyze(const EventProxyBase& iEvent)
{
	return analyze(iEvent, NULL);
}





