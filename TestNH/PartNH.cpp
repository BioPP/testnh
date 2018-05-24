//
// File: PartNH.cpp
// Created by: Julien Dutheil
// Created on: Dec Thu 09 11:11 2010
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to describe
the patterns of substitutions along a phylogeny using substitution mapping.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "partnh parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the package manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

class Partition {
  public:
    size_t number;
    size_t size;
};

vector<const Node*> getCandidateNodesForThreshold(map<double, vector<const Node*> >& sortedHeights, double threshold) {
  vector<const Node*> candidates;
  for (map<double, vector<const Node*> >::iterator it = sortedHeights.begin(); it != sortedHeights.end(); ++it) {
    if (it->first <= threshold) {
      for (vector<const Node*>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        for (size_t k = 0; k < (*it2)->getNumberOfSons(); ++k)
          candidates.push_back((*it2)->getSon(k));
    }
  }
  return candidates;
}

vector< vector<int> > getGroups(vector<const Node*>& candidates) {
  map<int, Partition> partitions;
  for (size_t i = 0; i < candidates.size(); ++i) {
    vector<string> ids = TreeTemplateTools::getLeavesNames(*candidates[i]);
    size_t size = ids.size();
    for (size_t j = 0; j < ids.size(); ++j) {
      Partition* part = &partitions[TextTools::toInt(ids[j])];
      if (part->size == 0) { //new partition
        part->number = i;
        part->size = size;
      } else { //Partition already exists.
        //Check if this is a smaller partition:
        if (part->size > size) {
          part->number = i;
          part->size = size;
        }
      }
    }
  }
  //Now reorganize results:
  map<size_t, vector<int> > groups;
  for (map<int, Partition>::iterator it = partitions.begin(); it != partitions.end(); ++it) {
    groups[it->second.number].push_back(it->first);
  }
  //And renumber partitions:
  vector< vector<int> > groups2;
  for (map<size_t, vector<int> >::const_iterator it = groups.begin(); it != groups.end(); it++)
    groups2.push_back(it->second);

  //Return results
  return groups2;
}
 
SubstitutionModelSet* buildModelSetFromPartitions(
    const SubstitutionModel* model,
    const FrequenciesSet* rootFreqs,
    const Tree* tree,
    const vector< vector<int> >& groups,
    const vector<string>& globalParameterNames,
    std::map<int, ParameterList>& initParameters
  ) throw (AlphabetException, Exception)
{
  //Check alphabet:
  if (rootFreqs && model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("SubstitutionModelSetTools::createNonHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  ParameterList globalParameters, branchParameters;
  globalParameters = model->getParameters();
  vector<string> globalParameterPrefs; // vector of the prefixes (when there is a '*' in the declaration)

  vector<string> globalParameterNames2; // vector of the complete names of global parameters

  //First get correct parameter names
  
  for (size_t i = 0; i < globalParameterNames.size(); i++) {
    if (globalParameterNames[i].find("*") != string::npos) {

      for (size_t j = 0; j < globalParameters.size(); j++)
        {
          StringTokenizer stj(globalParameterNames[i], "*", true, false);
          size_t pos1, pos2;
          string parn=globalParameters[j].getName();
          bool flag(true);
          string g=stj.nextToken();
          pos1=parn.find(g);
          if (pos1!=0)
            flag=false;
          pos1+=g.length();
          while (flag && stj.hasMoreToken()){
            g=stj.nextToken();
            pos2=parn.find(g,pos1);
            if (pos2 == string::npos){
              flag=false;
              break;
            }
            pos1=pos2+g.length();
          }
          if (flag &&
              ((g.length()==0) || (pos1==parn.length()) || (parn.rfind(g)==parn.length()-g.length())))
            globalParameterNames2.push_back(parn);
        }
    }
    else if (!globalParameters.hasParameter(globalParameterNames[i]))
      throw Exception("PartNH::buildModelSetFromPartitions. Parameter '" + globalParameterNames[i] + "' is not valid.");
    else
      globalParameterNames2.push_back(globalParameterNames[i]);
  }

  // remove non global parameters
  for (size_t i = globalParameters.size(); i > 0; i--)
    {
      if (find(globalParameterNames2.begin(), globalParameterNames2.end(), globalParameters[i-1].getName()) == globalParameterNames2.end())
        {
          //not a global parameter:
          branchParameters.addParameter(globalParameters[i - 1]);
          globalParameters.deleteParameter(i - 1);
        }
    }


  SubstitutionModelSet* modelSet = rootFreqs ?
    new SubstitutionModelSet(model->getAlphabet(), rootFreqs->clone()) :
    new SubstitutionModelSet(model->getAlphabet());

  //We assign a copy of this model to all nodes in the tree, for each partition, and link all parameters with it.
  for (size_t i = 0; i < groups.size(); ++i) {
    SubstitutionModel* modelC = dynamic_cast<SubstitutionModel*>(model->clone());
    modelC->matchParametersValues(initParameters[groups[i][0]]);
    modelSet->addModel(modelC, groups[i]);
  }

  // Now alias all global parameters on all nodes:
  for (size_t i=0; i < globalParameters.size(); i++)
    {
      string pname=globalParameters[i].getName();

      for (size_t nn = 1; nn < groups.size(); nn++)
        modelSet->aliasParameters(pname+"_1",pname+"_"+TextTools::toString(nn+1));
    }

  return modelSet;
}

ParameterList getParametersToEstimate(const DRTreeLikelihood* drtl, map<string, string>& params) {
  // Should I ignore some parameters?

  if (params.find("optimization.ignore_parameter") != params.end())
    throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");
  
  ParameterList parametersToEstimate = drtl->getParameters();
  vector<string> parNames = parametersToEstimate.getParameterNames();
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", "", true, false);
  StringTokenizer st(paramListDesc, ",");
  while (st.hasMoreToken())
  {
    try
    {
      string param = st.nextToken();
      size_t starpos;
      if (param == "BrLen")
      {
        vector<string> vs = drtl->getBranchLengthsParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if ((starpos = param.find("*")) != string::npos)
      {
        vector<string> vs = ApplicationTools::matchingParameters(param, parNames);
        
        for (vector<string>::iterator it = vs.begin(); it != vs.end(); it++)
        {
          parametersToEstimate.deleteParameter(*it);
          ApplicationTools::displayResult("Parameter ignored", *it);
        }
      }
      else
      {
        parametersToEstimate.deleteParameter(param);
        ApplicationTools::displayResult("Parameter ignored", param);
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
  }

  return parametersToEstimate;
}

void estimateLikelihood(DRTreeLikelihood* drtl, ParameterList& parametersToEstimate, double tolerance, unsigned int nbEvalMax, OutputStream* messageHandler, OutputStream* profiler, bool reparam, bool useClock, unsigned int verbose) {
  // Uses Newton-raphson algorithm with numerical derivatives when required.
  parametersToEstimate.matchParametersValues(drtl->getParameters());
  OptimizationTools::optimizeNumericalParameters2(
    dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(drtl), parametersToEstimate,
    0, tolerance, nbEvalMax, messageHandler, profiler, reparam, useClock, verbose, OptimizationTools::OPTIMIZATION_NEWTON);
}

void outputNHModel(const string& modelPath, double likelihood, const SubstitutionModelSet* modelSet, const DiscreteDistribution* rDist) {
  StlOutputStream out(new ofstream(modelPath.c_str(), ios::out));
  out << "# Log likelihood = ";
  out.setPrecision(20) << likelihood;
  out.endLine();
  out.endLine();
  out << "# Substitution model parameters:";
  out.endLine();
  PhylogeneticsApplicationTools::printParameters(modelSet, out);
  out.endLine();
  PhylogeneticsApplicationTools::printParameters(rDist   , out);
  out.endLine();
}

void outputHModel(const string& modelPath, double likelihood, const SubstitutionModel* model, const DiscreteDistribution* rDist) {
  StlOutputStream out(new ofstream(modelPath.c_str(), ios::out));
  out << "# Log likelihood = ";
  out.setPrecision(20) << likelihood;
  out.endLine();
  out.endLine();
  out << "# Substitution model parameters:";
  out.endLine();
  PhylogeneticsApplicationTools::printParameters(model, out);
  out.endLine();
  PhylogeneticsApplicationTools::printParameters(rDist, out);
  out.endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     PartNH, version 1.3.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 24/11/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  BppApplication partnh(args, argv, "PartNH");
  partnh.startTimer();

  Newick newick;
  string clusterTree = ApplicationTools::getAFilePath("input.cluster_tree.file", partnh.getParams(), true, true);
  ApplicationTools::displayResult("Input cluster tree", clusterTree);
  TreeTemplate<Node>* htree = newick.read(clusterTree);
  
  //We only read NHX tree because we want to be sure to use the correct id:
  //TreeTemplate<Node>* ptree = dynamic_cast<TreeTemplate<Node>*>(PhylogeneticsApplicationTools::getTree(partnh.getParams()));
  string treeIdPath = ApplicationTools::getAFilePath("input.tree.file", partnh.getParams(), true, true);
  ApplicationTools::displayResult("Input tree file", treeIdPath);
  Nhx nhx(true);
  TreeTemplate<Node>* ptree = nhx.read(treeIdPath);

  map<const Node*, double> heights;
  TreeTemplateTools::getHeights(*htree->getRootNode(), heights);
  map<double, vector<const Node*> > sortedHeights;

  string method = ApplicationTools::getStringParameter("partition.method", partnh.getParams(), "threshold");
  for (map<const Node*, double>::iterator it = heights.begin(); it != heights.end(); ++it)
    sortedHeights[max(0., 1. - 2. * it->second)].push_back(it->first);

  //This will contains the final groups of nodes:
  vector< vector<int> > groups;

  if (method == "threshold") {
    double threshold = ApplicationTools::getDoubleParameter("partition.threshold", partnh.getParams(), 0.01);
    ApplicationTools::displayResult("Output partitions for threshold", threshold);
    vector<const Node*> candidates = getCandidateNodesForThreshold(sortedHeights, threshold);
    ApplicationTools::displayResult("Number of nested partitions", candidates.size());
   
    groups = getGroups(candidates);
    ApplicationTools::displayResult("Number of real partitions", groups.size());
    //Display partitions:
    for (size_t i = 0; i < groups.size(); ++i)
      ApplicationTools::displayResult("Partition " + TextTools::toString(i + 1), TextTools::toString(groups[i].size()) + " element(s).");
  
  } else if (method == "auto") {
    //First we need to get the alphabet and data:
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(partnh.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", partnh.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, partnh.getParams());
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, partnh.getParams());
    delete allSites;
    size_t nbSites = sites->getNumberOfSites();

    ApplicationTools::displayResult("Number of sequences", sites->getNumberOfSequences());
    ApplicationTools::displayResult("Number of sites", nbSites);
 
    //Then we need the model to be used, including substitution model, rate distribution and root frequencies set.
    //We also need to specify the parameters that will be shared by all partitions.
    SubstitutionModel* model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, partnh.getParams());
    DiscreteDistribution* rDist = 0;
    if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
    {
      // Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(partnh.getParams());
    }
    //We initialize the procedure by estimating the homogeneous model.
    DRTreeLikelihood* drtl = 0;
    if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
      drtl = new DRHomogeneousTreeLikelihood(*ptree, *sites, model, rDist, true);
    else
      throw Exception("Mixed models not supported so far.");
      //drtl = new DRHomogeneousMixedTreeLikelihood(*ptree, *sites, model, rDist, true);
    drtl->initialize();

    ApplicationTools::displayResult("Log-Likelihood", drtl->getLogLikelihood());

    //Check for saturation:
    double ll = drtl->getValue();
    if (std::isinf(ll))
    {
      ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
      if (codonAlphabet)
      {
        bool f = false;
        size_t s;
        for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
          if (std::isinf(drtl->getLogLikelihoodForASite(i))) {
            const Site& site = sites->getSite(i);
            s = site.size();
            for (size_t j = 0; j < s; j++) {
              if (gCode->isStop(site.getValue(j))) {
                (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                f = true;
              }
            }
          }
        }
        if (f)
          exit(-1);
      }
      bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", partnh.getParams(), false, "", true, 1);
      if (!removeSaturated) {
        ofstream debug ("DEBUG_likelihoods.txt", ios::out);
        for (size_t i = 0; i < sites->getNumberOfSites(); i++)
        {
          debug << "Position " << sites->getSite(i).getPosition() << " = " << drtl->getLogLikelihoodForASite(i) << endl; 
        }
        debug.close();
        ApplicationTools::displayError("!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
        ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
        ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
        exit(1);
      } else {
        for (size_t i = sites->getNumberOfSites(); i > 0; --i) {
          if (std::isinf(drtl->getLogLikelihoodForASite(i - 1))) {
            ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
            sites->deleteSite(i - 1);
          }
        }
        ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());
        drtl->setData(*sites);
        drtl->initialize();
        ll = drtl->getValue();
        if (std::isinf(ll)) {
          throw Exception("Likelihood is still 0 after saturated sites are removed! Looks like a bug...");
        }
        ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-ll, 15));
      }
    }




    //We get this now, so that it does not fail after having done the full optimization!!
    string modelPath = ApplicationTools::getAFilePath("output.model.file", partnh.getParams(), true, false);
    bool outputIntermediateModels = ApplicationTools::getBooleanParameter("output.intermediate.models", partnh.getParams(), false);
    ApplicationTools::displayBooleanResult("Output intermediate models", outputIntermediateModels);
    string logPath = ApplicationTools::getAFilePath("output.log.file", partnh.getParams(), false, false);
    ApplicationTools::displayResult("Output log to", logPath);
    unique_ptr<ofstream> logout;
    if (logPath != "none" && logPath != "None") {
      logout.reset(new ofstream(logPath.c_str(), ios::out));
      *logout << "Threshold\tNbPartitions\tLogL\tDF\tAIC\tAICc\tBIC" << endl;
    }

    //Optimize parameters
    unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", partnh.getParams(), 2);
    
    string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", partnh.getParams(), false, false);
    shared_ptr<OutputStream> messageHandler(0);
    if (mhPath != "none")
    {
      if (mhPath == "std")
        messageHandler= ApplicationTools::message;
      else
        messageHandler=shared_ptr<OutputStream>(new StlOutputStream(new ofstream(mhPath.c_str(), ios::out)));
    }
    
    ApplicationTools::displayResult("Message handler", mhPath + "*");

    string prPath = ApplicationTools::getAFilePath("optimization.profiler", partnh.getParams(), false, false);
    shared_ptr<OutputStream> profiler(0);
    if (prPath != "none")
    {
      if (prPath == "std")
        messageHandler= ApplicationTools::message;
      else
        messageHandler=shared_ptr<OutputStream>(new StlOutputStream(new ofstream(prPath.c_str(), ios::out)));
    }
    
    if (profiler.get()) profiler->setPrecision(20);
    ApplicationTools::displayResult("Profiler", prPath + "*");

    ParameterList parametersToEstimate = getParametersToEstimate(drtl, partnh.getParams()); 

    unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", partnh.getParams(), 1000000);
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

    double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", partnh.getParams(), .000001);
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

    bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", partnh.getParams(), false);
    ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));
    
    string clock = ApplicationTools::getStringParameter("optimization.clock", partnh.getParams(), "None", "", true, false);
    if (clock != "None" && clock != "Global")
      throw Exception("Molecular clock option not recognized, should be one of 'Global' or 'None'.");
    bool useClock = (clock == "Global");
    ApplicationTools::displayResult("Molecular clock", clock);

    estimateLikelihood(drtl, parametersToEstimate, tolerance, nbEvalMax, messageHandler.get(), profiler.get(), reparam, useClock, optVerbose);

    double logL = drtl->getValue();
    double df = static_cast<double>(drtl->getParameters().size());
    double aic  = 2. * (df + logL);
    double aicc = aic + 2 * df * (df + 1) / (static_cast<double>(nbSites) - df - 1);
    double bic  = 2. * logL + df * log(nbSites);
    ApplicationTools::displayResult("* Homogeneous model - LogL", -logL);
    ApplicationTools::displayResult("                    - df", df);
    if (logout.get())
      *logout << 0. << "\t" << 1. << "\t" << -logL << "\t" << df << "\t" << aic << "\t" << aicc << "\t" << bic << endl;
    
    //Get necessary things for building a non-homogeneous model:
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      // Markov-Modulated Markov Model...
      size_t n = static_cast<size_t>(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1. / (double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                     // we should assume a rate distribution for the root also!!!
    }

    bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", partnh.getParams(), false, "", false, false);
    FrequenciesSet* rootFreqs = 0;
    std::map<std::string, std::string> aliasFreqNames;
    if (!stationarity)
    {
      rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, partnh.getParams(), aliasFreqNames, rateFreqs);
      stationarity = !rootFreqs;
    }
    ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);
    if (!stationarity && !ptree->isRooted())
      ApplicationTools::displayWarning("An (most likely) unrooted tree is in use with a non-stationary model! This is probably not what you want...");
   
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous.shared_parameters", partnh.getParams(), ',', "");
    for (size_t i = 0; i < globalParameters.size(); i++)
      ApplicationTools::displayResult("Global parameter", globalParameters[i]);

    string likelihoodComparison = ApplicationTools::getStringParameter("partition.test", partnh.getParams(), "LRT");
    double testThreshold = 0.;
    if (likelihoodComparison == "LRT")
      testThreshold = ApplicationTools::getDoubleParameter("partition.test.threshold", partnh.getParams(), 0.01);
    else if (likelihoodComparison == "AIC") {}
    else if (likelihoodComparison == "AICc") {}
    else if (likelihoodComparison == "BIC") {}
    else
      throw Exception("Unknown likelihood comparison method: " + likelihoodComparison);
    ApplicationTools::displayResult("AIC" , aic);
    ApplicationTools::displayResult("AICc", aicc);
    ApplicationTools::displayResult("BIC" , bic);
    ApplicationTools::displayResult("Likelihood comparison method", likelihoodComparison);
      
    string stopCond = ApplicationTools::getStringParameter("partition.test.stop_condition", partnh.getParams(), "1");
    int stop;
    if (stopCond == "all") {
      stop = static_cast<int>(-log(0)); //Inf
      ApplicationTools::displayResult("Stop condition", string("test all nested models and report global best"));
    } else {
      stop = TextTools::toInt(stopCond);
      if (stop < 1)
        throw Exception("Stop parameter should be at least 1!.");
      ApplicationTools::displayResult("Stop condition", string("test ") + TextTools::toString(stop) + string(" nested models after local best"));
    }

    //In order to save optimization time, we will same parameter estimates after each model was estimated
    //in order to use smart initial values for the next round...
    map<int, ParameterList> currentParameters;
    //At start, all nodes have the same parameter values (homogeneous model):
    vector<int> ids = ptree->getNodesId();
    ids.pop_back();
    for (size_t i = 0; i < ids.size(); ++i) {
      currentParameters[ids[i]] = model->getParameters();
    }

    //Now try more and more complex non-homogeneous models, using the clustering tree set as input.
    vector<const Node*> candidates;
    bool moveForward = true;
    map<double, vector<const Node*> >::iterator it = sortedHeights.begin();
    double currentThreshold = 1.;
    
    SubstitutionModelSet* modelSet = 0;
    SubstitutionModelSet* bestModelSet = 0;
    DiscreteDistribution* bestRDist = 0;
    TreeTemplate<Node>*   bestTree = 0;
    double bestLogL = logL;
    double bestAic  = aic;
    double bestAicc = aicc;
    double bestBic  = bic;
    int previousBest = -1;
    vector< vector<int> > bestGroups;
    size_t modelCount = 0;
    map<int, vector<size_t> > partRecord;

    while (moveForward && it != sortedHeights.end()) {
      modelCount++;
      if (modelCount == 1 && !stationarity) {
        candidates.push_back(htree->getRootNode());
        ApplicationTools::displayBooleanResult("Testing non-stationary model", true);
      } else {
        for (vector<const Node*>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
          for (size_t k = 0; k < (*it2)->getNumberOfSons(); ++k)
            candidates.push_back((*it2)->getSon(k));
        }
        currentThreshold = it->first;
        ApplicationTools::displayResult("Current threshold", currentThreshold);
        it++;
      }
      
      //Get the corresponding partitions:
      vector< vector<int> > newGroups = getGroups(candidates);
      ApplicationTools::displayResult("Number of real partitions", newGroups.size());
      
      //Display and record partitions:
      for (size_t i = 0; i < newGroups.size(); ++i) {
        ApplicationTools::displayResult("Partition " + TextTools::toString(i + 1), TextTools::toString(newGroups[i].size()) + " element(s).");
        for (size_t j = 0; j < newGroups[i].size(); ++j) {
          partRecord[newGroups[i][j]].push_back(i + 1);
        }
      }

      //Now we have to build the corresponding model set:
      SubstitutionModelSet* newModelSet = buildModelSetFromPartitions(model, rootFreqs, ptree, newGroups, globalParameters, currentParameters);
      DiscreteDistribution* newRDist = rDist->clone();
      ParameterList previousParameters = drtl->getBranchLengthsParameters();
      previousParameters.addParameters(drtl->getRateDistributionParameters());
      if (!stationarity && modelCount > 1)
        previousParameters.addParameters(dynamic_cast<DRNonHomogeneousTreeLikelihood*>(drtl)->getSubstitutionModelSet()->getRootFrequenciesParameters());
      delete drtl;
      if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
        drtl = new DRNonHomogeneousTreeLikelihood(*ptree, *sites, newModelSet, newRDist, false);
      else
        throw Exception("Mixed models not supported so far.");
        //drtl = new DRNonHomogeneousMixedTreeLikelihood(*ptree, *sites, modelSet, rDist, true);
      drtl->initialize();
      drtl->matchParametersValues(previousParameters); //This will save some time during optimization!
      
      //Optimize parameters
      
      if (mhPath == "none")
        messageHandler=0;
      else
        if (mhPath == "std") 
          messageHandler = ApplicationTools::message;
        else
          messageHandler.reset(new StlOutputStream(new ofstream((mhPath + TextTools::toString(modelCount)).c_str(), ios::out)));

      if (prPath == "none")
        profiler=0;
      else
        if (prPath == "std")
          profiler = ApplicationTools::message;
        else
          profiler.reset(new StlOutputStream(new ofstream((prPath + TextTools::toString(modelCount)).c_str(), ios::out)));
    
      if (profiler.get()) profiler->setPrecision(20);

      //Reevaluate parameters, as there might be some change when going to a NH model:
      parametersToEstimate = getParametersToEstimate(drtl, partnh.getParams()); 
      
      estimateLikelihood(drtl, parametersToEstimate, tolerance, nbEvalMax, messageHandler.get(), profiler.get(), reparam, useClock, optVerbose);
        
      double newLogL = drtl->getValue();
      double newDf = static_cast<double>(drtl->getParameters().size());
      ApplicationTools::displayResult("* New NH model - LogL", -newLogL);
      ApplicationTools::displayResult("               - df", newDf);
      double newAic  = 2. * (newDf + newLogL);
      double newAicc = newAic + 2 * newDf * (newDf + 1) / (static_cast<double>(nbSites) - newDf - 1);
      double newBic  = 2. * newLogL + newDf * log(nbSites);
      TreeTemplate<Node>* newTree = new TreeTemplate<Node>(drtl->getTree());
      
      //Print model:
      if (outputIntermediateModels) {
        outputNHModel(modelPath + TextTools::toString(modelCount), -newLogL, newModelSet, newRDist);
      }

      double d = 2 * (logL - newLogL);
      double pvalue = 1. - RandomTools::pChisq(d, newDf - df);
      ApplicationTools::displayResult("               - LRT p-value", pvalue);
      ApplicationTools::displayResult("               - AIC", newAic);
      ApplicationTools::displayResult("               - AICc", newAicc);
      ApplicationTools::displayResult("               - BIC", newBic);

      if (logout.get())
        *logout << currentThreshold << "\t" << newGroups.size() << "\t" << -newLogL << "\t" << newDf << "\t" << newAic << "\t" << newAicc << "\t" << newBic << endl;

      //Finally compare new model to the current one:
      bool saveThisModel = false;
      bool test = false;
      if (likelihoodComparison == "LRT") {
        test = (pvalue <= testThreshold);
        if (test) {
          saveThisModel = true;
          moveForward = false;
        }
      } else {
        if (likelihoodComparison == "AIC") {
          test = (newAic < bestAic);
        } else if (likelihoodComparison == "AICc") {
          test = (newAicc < bestAicc);
        } else { //BIC
          test = (newBic < bestBic);
        }
        if (newAic < bestAic)
          bestAic = newAic;
        if (newAicc < bestAicc)
          bestAicc = newAicc;
        if (newBic < bestBic)
          bestBic = newBic;
        if (test) {
          saveThisModel = true;
          previousBest = 0;
        } else {
          previousBest++;
        }
        if (previousBest == stop) {
          //We have to stop here
          moveForward = false;
        }
      }

      if (saveThisModel) {
        if (bestModelSet) {
          delete bestModelSet;
          delete bestRDist;
          delete bestTree;
        }
        bestModelSet = newModelSet->clone();
        bestRDist    = newRDist->clone();
        bestTree     = newTree->clone();
        bestGroups   = newGroups;
        bestLogL     = newLogL;
      }

      //Moving forward:
      if (moveForward) {
        delete ptree;
        if (modelSet)
          delete modelSet;
        delete rDist;
        ptree    = newTree;
        modelSet = newModelSet;
        rDist    = newRDist;
        logL     = newLogL;
        df       = newDf;
        aic      = newAic;
        aicc     = newAicc;
        bic      = newBic;
        groups   = newGroups;
        //Save parameters:
        for (size_t i = 0; i < ids.size(); ++i) {
          currentParameters[ids[i]] = modelSet->getModelForNode(ids[i])->getParameters();
        }
      } else {
        delete newModelSet;
        delete newRDist;
        delete newTree;
      }
    }
    //Write best model to file and output partition tree.
    ApplicationTools::displayResult("Model description output to file", modelPath);
    //We have to distinguish two cases...
    if (bestModelSet) {
      outputNHModel(modelPath, bestLogL, bestModelSet, bestRDist);
    } else {
      outputHModel(modelPath, bestLogL, model, rDist);
    }

    //Write parameter estimates per node:
    string paramPath = ApplicationTools::getAFilePath("output.parameters.file", partnh.getParams(), false, false);
    ApplicationTools::displayResult("Output parameter table to", paramPath);
    if (paramPath != "none") {
      ofstream paramFile(paramPath.c_str(), ios::out);
      vector<string> paramNames = model->getParameters().getParameterNames();
      paramFile << "NodeId";
      for (size_t i = 0; i < paramNames.size(); ++i)
        paramFile << "\t" << paramNames[i];
      paramFile << endl;
      if (bestModelSet) {
        for (unsigned k = 0; k < bestModelSet->getNumberOfModels(); ++k) {
          ParameterList pl = bestModelSet->getSubstitutionModel(k)->getParameters();
          vector<int> idsk = bestModelSet->getNodesWithModel(k);
          for (size_t j = 0; j < idsk.size(); ++j) {
            paramFile << idsk[j];
            for (size_t i = 0; i < paramNames.size(); ++i) {
              paramFile << "\t" << pl.getParameter(paramNames[i]).getValue();
            }
            paramFile << endl;
          }
        }
      } else {
        //All nodes have the same parameters:
        ParameterList pl = model->getParameters();
        for (size_t j = 0; j < ids.size(); ++j) {
          paramFile << ids[j];
          for (size_t i = 0; i < paramNames.size(); ++i) {
            paramFile << "\t" << pl.getParameter(paramNames[i]).getValue();
          }
          paramFile << endl;
        }
      }
      paramFile.close();
    }

    //Print partition record:
    string partRecordPath = ApplicationTools::getAFilePath("output.partitions.record", partnh.getParams(), false, false);
    ApplicationTools::displayResult("Output partitions record to", partRecordPath);
    if (partRecordPath != "none") {
      ofstream partRecordFile(partRecordPath.c_str(), ios::out);
      for (map<int, vector<size_t> >::iterator pit = partRecord.begin(); pit != partRecord.end(); ++pit) {
        partRecordFile << pit->first;
        for (size_t i = 0; i < pit->second.size(); ++i) {
          partRecordFile << "\t" << pit->second[i];
        }
        partRecordFile << endl;
      }
      partRecordFile.close();
    }

    //Cleaning:
    if (bestTree) {
      ptree  = bestTree;
      groups = bestGroups;
    } else {
      groups.clear();
      groups.push_back(ptree->getNodesId());
    }
    delete bestModelSet;
    delete bestRDist;
  
    if (logout.get())
      logout->close();

  } else throw Exception("Unknown option: " + method);

  //Write best tree:
  PhylogeneticsApplicationTools::writeTree(*ptree, partnh.getParams());
 
  //Now write partitions to file:
  string partPath = ApplicationTools::getAFilePath("output.partitions.file", partnh.getParams(), true, false);
  ApplicationTools::displayResult("Partitions output to file", partPath);
  for (size_t i = 0; i < groups.size(); ++i) {
    for (size_t j = 0; j < groups[i].size(); ++j) {
      Node* node = ptree->getNode(groups[i][j]);
      if (node->hasName())
        node->setName(TextTools::toString(i + 1) + "_" + node->getName());
      node->setBranchProperty("partition", BppString(TextTools::toString(i + 1)));
    }
  }
  newick.enableExtendedBootstrapProperty("partition");
  newick.write(*ptree, partPath);

  //Cleaning memory:
  delete htree;
  delete ptree;
  partnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}
