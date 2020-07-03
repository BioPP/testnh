//
// File: StepNH.cpp
// Created by: Julien Dutheil
// Created on: Dec Sat 03 16:41 2005
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to perform
a model selection with non-homogeneous models.

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

// From bpp-seq:
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/NumCalcApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>

using namespace bpp;

/******************************************************************************/

//For debugging:
class MyListener: public OptimizationListener
{
  private:
    Function* f_;

  public:
    MyListener(Function* f): f_(f) {}
    MyListener(const MyListener& ml): f_(ml.f_) {}
    MyListener& operator=(const MyListener& ml)
    {
      f_ = ml.f_;
      return (*this);
    }
    virtual ~MyListener() {}

  public:
    void optimizationInitializationPerformed(const OptimizationEvent & event) {}
    void optimizationStepPerformed(const OptimizationEvent & event) 
    {
      cout << endl << " -ll = " << setprecision(20) << f_->getValue();
    }
    bool listenerModifiesParameters() const { return false; }
};

class SortablePair
{
  public:
    unsigned int i;
    unsigned int j;
    double value;
};

bool cmp(const SortablePair& sp1, const SortablePair& sp2)
{
  return sp1.value > sp2.value;
};

class ParamNameAndIndex
{
  public:
    string name;
    unsigned int index;

  public:
    ParamNameAndIndex(const string & n, unsigned int i): name(n), index(i) {}
};

bool operator<(const ParamNameAndIndex& pnai1, const ParamNameAndIndex& pnai2)
{ 
  return pnai1.name <  pnai2.name;
}
bool operator!=(const ParamNameAndIndex& pnai1, const ParamNameAndIndex& pnai2)
{ 
  return pnai1.name !=  pnai2.name;
}

bool areConnected(const TreeTemplate<Node>& tree, const vector<int>& set1, const vector<int>& set2)
{
  for(unsigned int i = 0; i < set1.size(); i++)
  {
    int id1 = set1[i];
    const Node * node1 = tree.getNode(id1);
    for(unsigned int j = 0; j < set2.size(); j++)
    {
      int id2 = set2[i];
      const Node * node2 = tree.getNode(id2);
      if(node2->getFather()->getId() == id1 || node1->getFather()->getId() == id2) return true;
    }
  }
  return false;
}
    
void stepBackward(NonHomogeneousTreeLikelihood* nhtl,
    double threshold, bool connectedBranchesOnly, bool simultaneous, double precision,
    OutputStream* profiler, OutputStream* messenger, unsigned int verbose)
{
  SubstitutionModelSet* modelSet = nhtl->getSubstitutionModelSet();

  //Create gain matrix:
  ParameterList params = modelSet->getNodeParameters();
  unsigned int nbParams = params.size();
  vector<string> paramModelNames(nbParams);
  vector< vector<int> > nodesId(nbParams);
  string name;
  for (unsigned int i = 0; i < nbParams; i++)
  {
    name = params[i].getName();
    paramModelNames[i] = modelSet->getParameterModelName(name);
    nodesId[i] = modelSet->getNodesWithParameter(name);
  }
  vector<bool> mask(nbParams, true);

  MyListener listener(nhtl);
  //Adjust current model:
  ApplicationTools::displayTask("Optimizing current model");
  OptimizationTools::optimizeNumericalParameters2(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(nhtl),
      nhtl->getParameters(),
      0, precision, 10000, messenger, profiler, false, false, max(0, (int)verbose - 2),
      OptimizationTools::OPTIMIZATION_NEWTON);
  ApplicationTools::displayTaskDone();

  bool mainTest = true;
  while (mainTest)
  {
    double currentL = nhtl->getValue();

    //Fill the matrix:
    vector<SortablePair> preVsp;
    ApplicationTools::displayTask("Computing pairwise amelioration scores", true);
    for (unsigned int i = 1; i < nbParams; i++)
    {
      ApplicationTools::displayGauge(i - 1, nbParams - 2, '=');
      if (!mask[i]) continue;
      string parami = params[i].getName();
      for (unsigned int j = 0; j < i; j++)
      {
        if (!mask[j]) continue;
        string paramj = params[j].getName();
        //Check if parameters i and j are mergeable:
        if (paramModelNames[i] != paramModelNames[j]) continue;
        //Merge only parameters from adjacent nodes?
        if (connectedBranchesOnly && !areConnected(dynamic_cast<const TreeTemplate<Node> &>(nhtl->getTree()), nodesId[i], nodesId[j])) continue;

        //Copy the current model set:
        SubstitutionModelSet* thisModelSet = modelSet->clone();
        //Merge parameters i and j.
        vector<unsigned int> models4j = thisModelSet->getModelsWithParameter(paramj);
        //First remove parameter j:
        thisModelSet->removeParameter(paramj);
        //Then associate parameter i to all models where j was previously found:
        for (size_t k = 0; k < models4j.size(); k++)
        {
          thisModelSet->setParameterToModel(thisModelSet->getParameterIndex(parami), models4j[k]);  
        }

        //Estimate the likelihood with the new model:
        NonHomogeneousTreeLikelihood* thisNhtl = nhtl->clone();
        thisNhtl->setSubstitutionModelSet(thisModelSet);

        BrentOneDimension* optimizer = new BrentOneDimension(thisNhtl);
        optimizer->setProfiler(0);
        optimizer->setMessageHandler(0);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        double x = params[i].getValue();
        if(x < 0.5)
          optimizer->setInitialInterval(x, x+0.01);
        else
          optimizer->setInitialInterval(x, x-0.01);
        optimizer->init(params.subList(i));
        optimizer->getStopCondition()->setTolerance(precision);
        optimizer->setVerbose(0);
        optimizer->optimize();
        delete optimizer;
        double thisL = thisNhtl->getValue();
        double stat = 2*thisL - 2*currentL;
        SortablePair sp;
        sp.i = i;
        sp.j = j;
        sp.value = 1. - RandomTools::pChisq(stat, 1);
        preVsp.push_back(sp);

        delete thisNhtl;
        delete thisModelSet;
      }
    }
    ApplicationTools::displayTaskDone();
    if (preVsp.size() == 0)
    {
      cout << "Stoping here..." << endl;
      mainTest = false; //No more test to perform
      break;
    }

    //Pick the most significant test / significant tests:
    sort(preVsp.begin(), preVsp.end(), cmp);
    double maxpval = preVsp[0].value;

    if (maxpval < threshold)
    {
      if (verbose > 0) cout << "Stoping here..." << endl;
      mainTest = false; //No more test to perform
      break;
    }
    else
    {
      //Check if each pair is not already present by connection of previous pairs:
      vector<SortablePair> vsp;
      for (size_t k = 0; k < preVsp.size(); k++)
      {
        SortablePair sp = preVsp[k];
        unsigned int testPresent = 0;
        for (size_t l = 0; l < vsp.size() && testPresent < 2; l++)
        {
          if (sp.i == vsp[l].i || sp.j == vsp[l].j) testPresent++;
        }
        if (testPresent < 2) vsp.push_back(sp);
      }

      //Compute the number of mergings to perform:
      unsigned int nbTests = 1;
      if (simultaneous)
      {
        while (nbTests < vsp.size() && vsp[nbTests].value > threshold) nbTests++;
      }

      //Merge the parameters:
      bool test = true;
      while (test)
      {
        SubstitutionModelSet *backupSet = modelSet->clone();
        vector<bool> backupMask = mask;
        unsigned int nbp = modelSet->getParameters().size();
        //Get the intersection set:
        vector<ParamNameAndIndex> paramNames;
        for (unsigned int i = 0; i < nbTests; i++)
        {
          ParamNameAndIndex pnai1(params[vsp[i].i].getName(), vsp[i].i);
          ParamNameAndIndex pnai2(params[vsp[i].j].getName(), vsp[i].j);
          paramNames.push_back(pnai1);
          paramNames.push_back(pnai2);
        }
        paramNames = VectorTools::unique<ParamNameAndIndex>(paramNames); //This will also sort the names.
        for (unsigned int i = 1; i < paramNames.size(); i++)
        {
          string parami = paramNames[i].name;
          string paramj = paramNames[i-1].name;
          if (modelSet->getParameterModelName(parami) != modelSet->getParameterModelName(paramj))
            continue;

          if(verbose > 1) cout << "Merging parameter " << parami << " and " << paramj << endl;
          //Then associate parameter i to all models where j was previously found:
          vector<unsigned int> models4j = modelSet->getModelsWithParameter(paramj);
          modelSet->removeParameter(paramj);
          for(unsigned int k = 0; k < models4j.size(); k++)
          {
            modelSet->setParameterToModel(modelSet->getParameterIndex(parami), models4j[k]);  
          }
          mask[paramNames[i-1].index] = false; //This parameter have been merged with i.
        }
        nhtl->setSubstitutionModelSet(modelSet); //Reaffect modified model set.
        
        //Adjust current model:
        ApplicationTools::displayTask("Optimizing current model");
        nhtl->setParameters(nhtl->getParameters());
//        OptimizationTools::optimizeNumericalParameters2(
//            dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(nhtl),
//            nhtl->getParameters(), 0, precision, 10000, messenger, profiler, false,
//            max(0, (int)verbose - 2), OptimizationTools::OPTIMIZATION_NEWTON);
        OptimizationTools::optimizeNumericalParameters(
            dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(nhtl),
            nhtl->getParameters(), 0, 1, precision, 10000, messenger, profiler, false,
            max(0, (int)verbose - 2), OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT);
         ApplicationTools::displayTaskDone();
        unsigned int newNbp = modelSet->getParameters().size();
        double thisL = nhtl->getValue();
        cout << setprecision(20) << currentL << "\t" << thisL << endl;
        double stat = 2*thisL - 2*currentL;
        double pvalue = 1. - RandomTools::pChisq(stat, nbp - newNbp);
        ApplicationTools::displayResult("Global LRT (ddf = " + TextTools::toString(nbp - newNbp) + ")", pvalue);
        if (pvalue < threshold)
        {
          if (verbose > 0) cout << "Moving foreward..." << endl;
          nhtl->setSubstitutionModelSet(backupSet);
          delete modelSet;
          modelSet = backupSet;
          mask = backupMask;
          if (nbTests == 1)
          {
            if(verbose > 0) cout << "Stoping here..." << endl;
            mainTest = false; //No more test to perform
            break;
          }
          nbTests = max((unsigned int)1, nbTests / 2);
        }
        else
        {
          test = false;
          delete backupSet;
        }
      }
    }
  } //while loop.

  //At the end of the procedure, the nhtl object contains the optimal model set.

}

void stepForward(NonHomogeneousTreeLikelihood* nhtl)
{

}

DataTable* getModelTable(const SubstitutionModelSet* modelSet, const vector<int>& ids)
{
  ParameterList pl = modelSet->getNodeParameters();
  DataTable* table = new DataTable(ids.size(),pl.size());
  table->setColumnNames(pl.getParameterNames());
  vector<string> rowNames(ids.size());
  for (size_t i = 0; i < ids.size(); i++)
    rowNames[i] = TextTools::toString(ids[i]);
  table->setRowNames(rowNames);
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    vector<int> thisIds = modelSet->getNodesWithParameter(pl[i].getName());
    for (size_t j = 0; j < thisIds.size(); j++)
    {
      (*table)(TextTools::toString(thisIds[j]), pl[i].getName()) = "X";
    }
  }
  return table;
}

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "stepnh parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the package manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                    Step NH, version 0.1.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  20/03/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 20/03/08 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  BppApplication stepnh(args, argv, "StepNH");
  stepnh.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(stepnh.getParams(), "", false);

  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, stepnh.getParams());
  
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, stepnh.getParams());
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  // Get the initial tree
  Tree* tree = PhylogeneticsApplicationTools::getTree(stepnh.getParams());
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  
  // Try to write the current tree to file. This will be overwritten by the optimized tree,
  // but allow to check file existence before running optimization!
  PhylogeneticsApplicationTools::writeTree(*tree, stepnh.getParams());
  
  RNonHomogeneousTreeLikelihood *tl;

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", stepnh.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  SubstitutionModelSet* modelSet = 0;
  DiscreteDistribution* rDist    = 0;

  if (nhOpt == "one_per_branch")
  {
    SubstitutionModel* model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, stepnh.getParams());
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(stepnh.getParams());
    }
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                   // we should assume a rate distribution for the root also!!!  
    }
    FrequencySet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequencySet(alphabet, sites, stepnh.getParams(), rateFreqs);

    string descGlobal = ApplicationTools::getStringParameter("nonhomogeneous_one_per_branch.shared_parameters", stepnh.getParams(), "", "", true, 1);

    NestedStringTokenizer nst(descGlobal,"[","]",",");
    const deque<string>& descGlobalParameters=nst.getTokens();
    
    map<string, vector<Vint> > globalParameters;
    for (const auto& desc:descGlobalParameters)
    {
      size_t post=desc.rfind("_");
      if (post==std::string::npos || post==desc.size()-1 || desc[post+1]!='[')
        globalParameters[desc]={};
      else
      {
        string key=desc.substr(0,post);
        Vint sint=NumCalcApplicationTools::seqFromString(desc.substr(post+2, desc.size()-post-3));
        if (globalParameters.find(key)==globalParameters.end())
          globalParameters[key]=vector<Vint>(1, sint);
        else
          globalParameters[key].push_back(sint);
      }
    }
    
    for (const auto& globpar:globalParameters)
    {
      ApplicationTools::displayResult("Global parameter", globpar.first);
      if (globpar.second.size()==0)
      {
        string all="All nodes";
        ApplicationTools::displayResult(" shared between nodes", all);
      }
      else
        for (const auto& vint:globpar.second)
          ApplicationTools::displayResult(" shared between nodes", VectorTools::paste(vint,","));
    }
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
    model = 0;
    
    tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, true);
  }
  else if (nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, 0, stepnh.getParams());
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(stepnh.getParams());
    }
    tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, true);
  }
  else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  tl->initialize();
 
    
  double logL = tl->getValue();
  if (std::isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
      if(pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
    }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
  }
  ApplicationTools::displayResult("Initial -log(likelihood)", TextTools::toString(logL, 15));
  if (std::isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    ApplicationTools::displayError("!!! Looking at each site:");
    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
    }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
  }

  double threshold = ApplicationTools::getDoubleParameter("step.lrt_threshold", stepnh.getParams(), 0.05);
  ApplicationTools::displayResult("LRT threshold", threshold);
  bool connected = ApplicationTools::getBooleanParameter("step.connected", stepnh.getParams(), false);
  ApplicationTools::displayResult("Merge only connected branches", connected ? "yes" : "no");
  bool sequential = ApplicationTools::getBooleanParameter("step.sequential", stepnh.getParams(), false);
  ApplicationTools::displayResult("Tests are sequential", sequential ? "yes" : "no");
  string algo = ApplicationTools::getStringParameter("step.algorithm", stepnh.getParams(), "backward");
  unsigned int verbose = ApplicationTools::getParameter<unsigned int>("step.optimization.verbose", stepnh.getParams(), 1);
  double precision = ApplicationTools::getDoubleParameter("step.optimization.tolerance", stepnh.getParams(), 0.000001);
  ApplicationTools::displayResult("Optimization precision", precision);

  string mhPath = ApplicationTools::getAFilePath("step.optimization.message_handler", stepnh.getParams(), false, false);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(mhPath.c_str(), ios::out)));
  if (verbose) ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("step.optimization.profiler", stepnh.getParams(), false, false);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(prPath.c_str(), ios::out)));
  if (profiler) profiler->setPrecision(20);
  if (verbose) ApplicationTools::displayResult("Profiler", prPath);

  if (algo == "backward")
    stepBackward(tl, threshold, connected, !sequential, precision, profiler, messageHandler, verbose);
  else if (algo == "forward")
    stepForward(tl);
  else throw Exception("Unknown algorithm " + algo);
  modelSet = tl->getSubstitutionModelSet();


/*  if(optimizeClock == "global")
  {
    PhylogeneticsApplicationTools::optimizeParameters(dynamic_cast<DiscreteRatesAcrossSitesClockTreeLikelihood *>(tl), stepnh.getParams());
  }
  else
  {
    tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(
        PhylogeneticsApplicationTools::optimizeParameters(tl, stepnh.getParams()));
  }
  
  tree = new TreeTemplate<Node>(*tl->getTree());
  PhylogeneticsApplicationTools::writeTree(* tree, stepnh.getParams()); */
  
  // Write parameters to screen:
  ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
  ParameterList parameters = tl->getSubstitutionModelParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  parameters = tl->getRateDistributionParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  // Write parameters to file:
  string parametersFile = ApplicationTools::getAFilePath("output.estimates", stepnh.getParams(), false, false);
  if (parametersFile != "none")
  {
    auto_ptr<ostream> out_ptr(new ofstream(parametersFile.c_str(), ios::out));
    StlOutputStream out(out_ptr);
    (out << "# Log likelihood = " << tl->getValue()).endLine().endLine();
    (out << "# Substitution model parameters:").endLine().endLine();
    modelSet->matchParametersValues(tl->getParameters());
    PhylogeneticsApplicationTools::printParameters(modelSet, out);
    out.endLine();
    (out << "# Rate distribution parameters:").endLine();
    rDist->matchParametersValues(tl->getParameters());
    PhylogeneticsApplicationTools::printParameters(rDist, out);
  }

  string treeWIdPath = ApplicationTools::getAFilePath("output.tree_wid.path", stepnh.getParams(), false, false);
  if( treeWIdPath != "none")
  {
    TreeTemplate<Node> tt(*tree);
    vector<Node*> nodes = tt.getNodes();
    for (size_t i = 0; i < nodes.size(); i++)
    {
      if (nodes[i]->isLeaf())
        nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
      else
        nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.write(*tree, treeWIdPath);
  }

  string tablePath = ApplicationTools::getAFilePath("output.model_table", stepnh.getParams(), false, false);
  if (tablePath != "none")
  {
    ofstream file(tablePath.c_str(), ios::out);
    DataTable* table = getModelTable(modelSet, tree->getNodesId());
    DataTable::write(*table, file, "\t");
    delete table;
  }

  delete alphabet;
  delete sites;
  delete modelSet;
  delete rDist;
  delete tl;
  delete tree;
  stepnh.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

