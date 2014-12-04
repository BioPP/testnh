//
// File: TestNH.cpp
// Created by: Julien Dutheil
//             Bastien Boussau
// Created on: Mar 18 15:27 2008
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to test the
homogeneity of the substitution process of a given alignment.

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
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "testnh parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the package manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

void simulate(
    NonHomogeneousSequenceSimulator * sim,
    size_t nbSites,
    double threshold,
    size_t nSim,
    size_t observedNbSignif,
    double observedMedian,
    double & nbSignifPValue,
    double & medianPValue,
    const string & distFile)
{
  ApplicationTools::displayTask("Perform simulations", true);

  auto_ptr<ofstream> out;
  if (distFile != "none") {
    out.reset(new ofstream(distFile.c_str(), ios::out));
    *out << "Sim\tNbSignif\tMedian" << endl;
  }

  SiteContainer * sites;
  size_t nbNbTPVal = 0;
  size_t nbMedPVal = 0;
  for (size_t k = 0; k < nSim; k++)
  {
    ApplicationTools::displayGauge(k, nSim-1, '=');

    sites = sim->simulate(nbSites);
    size_t nbSequences = sites->getNumberOfSequences();
    
    size_t nbTest = 0;
    vector<double> bstats;
    BowkerTest * test;
    for (size_t i = 0; i < nbSequences; i++)
    {
      for (size_t j = 0; j < i; j++)
      {
        test = SequenceTools::bowkerTest(sites->getSequence(i), sites->getSequence(j));
        if (test->getPValue() < threshold) nbTest++;
        bstats.push_back(test->getStatistic());
        delete test;
      }
    }
    double median = VectorTools::median(bstats);
    if (out.get())
      *out << k << "\t" << nbTest << "\t" << median << endl;
    if (nbTest >= observedNbSignif) nbNbTPVal++;
    if (median >= observedMedian  ) nbMedPVal++;
   
    delete sites;
  }
  if (out.get())
    out->close();
  nbSignifPValue = static_cast<double>(nbNbTPVal + 1) / static_cast<double>(nSim + 1);
  medianPValue   = static_cast<double>(nbMedPVal + 1) / static_cast<double>(nSim + 1);
  ApplicationTools::displayTaskDone();
}





int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                    Test NH, version 0.1.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  18/03/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 20/03/08 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  BppApplication testnh(args, argv, "TestNH");
  testnh.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(testnh.getParams(), "", false);
  auto_ptr<GeneticCode> gCode;
  CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
  if (codonAlphabet) {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", testnh.getParams(), "Standard", "", true, true);
    ApplicationTools::displayResult("Genetic Code", codeDesc);
      
    gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
  }
  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, testnh.getParams());
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, testnh.getParams());
  delete allSites;

  size_t nbSequences = sites->getNumberOfSequences();
  size_t nbSites = sites->getNumberOfSites();
  ApplicationTools::displayResult("Number of sequences", nbSequences);
  ApplicationTools::displayResult("Number of sites", nbSites);

  // 1) Compute Bowker's statistics
  unsigned nbTest = 0;
  vector<double> bstats;
  double threshold = ApplicationTools::getDoubleParameter("bowker_test.threshold", testnh.getParams(), 0.05);
  ApplicationTools::displayResult("Bowker's test threshold", threshold);
  
  BowkerTest* test;
  ApplicationTools::displayTask("Compute pairwise tests", true);
  for (size_t i = 0; i < nbSequences; i++)
  {
    ApplicationTools::displayGauge(i, nbSequences - 1, '=');
    for(size_t j = 0; j < i; j++)
    {
      test = SequenceTools::bowkerTest(sites->getSequence(i), sites->getSequence(j));
      if(test->getPValue() < threshold) nbTest++;
      bstats.push_back(test->getStatistic());
      delete test;
    }
  }
  double median = VectorTools::median(bstats);
  ApplicationTools::displayTaskDone();
  ApplicationTools::displayResult("Nb. signif. pairwise tests", nbTest);
  ApplicationTools::displayResult("Bowker statistic median", median);

  
  // 2) Run simulations
  NonHomogeneousSequenceSimulator* nhss;
  size_t nbSim = ApplicationTools::getParameter<size_t>("bootstrap.number", testnh.getParams(), 100);
  ApplicationTools::displayResult("Number of simulations", nbSim);
  double nbSignifPValue, medianPValue;
  string distFile;
  
  // 2.1) Get the initial tree and model
  Tree* tree = PhylogeneticsApplicationTools::getTree(testnh.getParams());
  ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
  ApplicationTools::displayBooleanResult("Is rooted", tree->isRooted());
 
  SubstitutionModel   * model    = 0;
  SubstitutionModelSet* modelSet = 0;
  DiscreteDistribution* rDist    = 0;
  DiscreteRatesAcrossSitesTreeLikelihood* tl = 0;

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", testnh.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  if (nhOpt == "no")
  {  
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, testnh.getParams());
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(testnh.getParams());
    }
  
    tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true, true);
  }
  else if (nhOpt == "one_per_branch")
  {
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, testnh.getParams());
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(testnh.getParams());
    }
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      size_t n = static_cast<size_t>(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                   // we should assume a rate distribution for the root also!!!  
    }
    std::map<std::string, std::string> aliasFreqNames;
    FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, testnh.getParams(), aliasFreqNames, rateFreqs);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", testnh.getParams(), ',', "");
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, aliasFreqNames, globalParameters); 
    model = 0;
      
    tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
  }
  else if (nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), 0, testnh.getParams());
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(testnh.getParams());
    }

    tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
  }
  else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  tl->initialize();

  ApplicationTools::displayResult("Likelihood", TextTools::toString(tl->getValue(), 15));
  
  delete tree;
  tree = new TreeTemplate<Node>(tl->getTree());
  delete tl;
   
  // 2.2) Simulate
  distFile = ApplicationTools::getAFilePath("bootstrap.dist_file", testnh.getParams(), false, false);
  ApplicationTools::displayResult("Null distribution in file", distFile);
  if (nhOpt == "no")
    nhss = new NonHomogeneousSequenceSimulator(model, rDist, tree);
  else
    nhss = new NonHomogeneousSequenceSimulator(modelSet, rDist, tree);
  simulate(nhss, nbSites, threshold, nbSim, nbTest, median, nbSignifPValue, medianPValue, distFile);
  ApplicationTools::displayResult("Global Bowker test number signif.", nbSignifPValue);
  ApplicationTools::displayResult("Bowker test median signif.", medianPValue);
  delete nhss;
 

  delete alphabet;
  delete sites;
  if(model)    delete model;
  if(modelSet) delete modelSet;
  delete rDist;
  delete tree;
  testnh.done();
 
  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

