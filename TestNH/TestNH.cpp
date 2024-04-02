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

// From bpp-core:
#include <Bpp/App/NumCalcApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Legacy/Simulation/NonHomogeneousSequenceSimulator.h>

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/NumCalcApplicationTools.h>
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
    NonHomogeneousSequenceSimulator& sim,
    size_t nbSites,
    double threshold,
    size_t nSim,
    size_t observedNbSignif,
    double observedMedian,
    double& nbSignifPValue,
    double& medianPValue,
    const string& distFile)
{
  ApplicationTools::displayTask("Perform simulations", true);

  unique_ptr<ofstream> out;
  if (distFile != "none")
  {
    out.reset(new ofstream(distFile.c_str(), ios::out));
    *out << "Sim\tNbSignif\tMedian" << endl;
  }

  size_t nbNbTPVal = 0;
  size_t nbMedPVal = 0;
  for (size_t k = 0; k < nSim; ++k)
  {
    ApplicationTools::displayGauge(k, nSim - 1, '=');

    auto sites = sim.simulate(nbSites);
    size_t nbSequences = sites->getNumberOfSequences();

    size_t nbTest = 0;
    vector<double> bstats;
    for (size_t i = 0; i < nbSequences; ++i)
    {
      for (size_t j = 0; j < i; ++j)
      {
        auto test = SequenceTools::bowkerTest(sites->sequence(i), sites->sequence(j));
        if (test->getPValue() < threshold)
          nbTest++;
        bstats.push_back(test->getStatistic());
      }
    }
    double median = VectorTools::median(bstats);
    if (out.get())
      *out << k << "\t" << nbTest << "\t" << median << endl;
    if (nbTest >= observedNbSignif)
      nbNbTPVal++;
    if (median >= observedMedian)
      nbMedPVal++;
  }
  if (out)
    out->close();
  nbSignifPValue = static_cast<double>(nbNbTPVal + 1) / static_cast<double>(nSim + 1);
  medianPValue   = static_cast<double>(nbMedPVal + 1) / static_cast<double>(nSim + 1);
  ApplicationTools::displayTaskDone();
}


int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                    Test NH, version 3.0.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  18/03/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 21/02/24 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }

  try
  {
    BppApplication testnh(args, argv, "TestNH");
    testnh.startTimer();

    shared_ptr<const Alphabet> alphabet = SequenceApplicationTools::getAlphabet(testnh.getParams(), "", false);
    shared_ptr<GeneticCode> gCode;
    auto codonAlphabet = std::dynamic_pointer_cast<const CodonAlphabet>(alphabet);
    if (codonAlphabet)
    {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", testnh.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);

      gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
    }
    auto allSites = SequenceApplicationTools::getSiteContainer(alphabet, testnh.getParams());
    shared_ptr<SiteContainerInterface> sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, testnh.getParams());

    size_t nbSequences = sites->getNumberOfSequences();
    size_t nbSites = sites->getNumberOfSites();
    ApplicationTools::displayResult("Number of sequences", nbSequences);
    ApplicationTools::displayResult("Number of sites", nbSites);

    // 1) Compute Bowker's statistics
    unsigned nbTest = 0;
    vector<double> bstats;
    double threshold = ApplicationTools::getDoubleParameter("bowker_test.threshold", testnh.getParams(), 0.05);
    ApplicationTools::displayResult("Bowker's test threshold", threshold);

    ApplicationTools::displayTask("Compute pairwise tests", true);
    for (size_t i = 0; i < nbSequences; i++)
    {
      ApplicationTools::displayGauge(i, nbSequences - 1, '=');
      for (size_t j = 0; j < i; j++)
      {
        auto test = SequenceTools::bowkerTest(sites->sequence(i), sites->sequence(j));
        if (test->getPValue() < threshold)
          nbTest++;
        bstats.push_back(test->getStatistic());
      }
    }
    double median = VectorTools::median(bstats);
    ApplicationTools::displayTaskDone();
    ApplicationTools::displayResult("Nb. signif. pairwise tests", nbTest);
    ApplicationTools::displayResult("Bowker statistic median", median);


    // 2) Run simulations
    unique_ptr<NonHomogeneousSequenceSimulator> nhss;
    size_t nbSim = ApplicationTools::getParameter<size_t>("bootstrap.number", testnh.getParams(), 100);
    ApplicationTools::displayResult("Number of simulations", nbSim);
    double nbSignifPValue, medianPValue;
    string distFile;

    // 2.1) Get the initial tree and model
    shared_ptr<Tree> tree = PhylogeneticsApplicationTools::getTree(testnh.getParams());
    ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
    ApplicationTools::displayBooleanResult("Is rooted", tree->isRooted());

    shared_ptr<SubstitutionModelInterface> model = nullptr;
    shared_ptr<SubstitutionModelSet> modelSet = nullptr;
    shared_ptr<DiscreteDistributionInterface> rDist = nullptr;
    shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl = nullptr;

    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", testnh.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    if (nhOpt == "no")
    {
      map<string, string> unparsed;
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode, sites, testnh.getParams(), unparsed);
      if (model->getNumberOfStates() > model->alphabet().getSize())
      {
        // Markov-modulated Markov model!
        rDist = make_shared<ConstantRateDistribution>();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(testnh.getParams());
      }

      tl = make_shared<RHomogeneousTreeLikelihood>(*tree, *sites, model, rDist, true, true, true);
    }
    else if (nhOpt == "one_per_branch")
    {
      map<string, string> unparsed;
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode, sites, testnh.getParams(), unparsed);
      if (model->getNumberOfStates() > model->alphabet().getSize())
      {
        // Markov-modulated Markov model!
        rDist = make_shared<ConstantRateDistribution>();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(testnh.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        // Markov-Modulated Markov Model...
        size_t n = static_cast<size_t>(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
        // we should assume a rate distribution for the root also!!!
      }
      std::map<std::string, std::string> aliasFreqNames;
      shared_ptr<FrequencySetInterface> rootFreqs = PhylogeneticsApplicationTools::getRootFrequencySet(alphabet, gCode, *sites, testnh.getParams(), aliasFreqNames, rateFreqs);

      string descGlobal = ApplicationTools::getStringParameter("nonhomogeneous_one_per_branch.shared_parameters", testnh.getParams(), "", "", true, 1);

      NestedStringTokenizer nst(descGlobal, "[", "]", ",");
      const deque<string>& descGlobalParameters = nst.getTokens();

      map<string, VVint> globalParameters;
      for (const auto& desc:descGlobalParameters)
      {
        size_t post = desc.rfind("_");
        if (post == std::string::npos || post == desc.size() - 1 || desc[post + 1] != '[')
          globalParameters[desc] = {}
        ;
        else
        {
          string key = desc.substr(0, post);
          Vint sint = NumCalcApplicationTools::seqFromString(desc.substr(post + 2, desc.size() - post - 3));
          if (globalParameters.find(key) == globalParameters.end())
            globalParameters[key] = vector<Vint>(1, sint);
          else
            globalParameters[key].push_back(sint);
        }
      }

      for (const auto& globpar:globalParameters)
      {
        ApplicationTools::displayResult("Global parameter", globpar.first);
        if (globpar.second.size() == 0)
        {
          string all = "All nodes";
          ApplicationTools::displayResult(" shared between nodes", all);
        }
        else
          for (const auto& vint:globpar.second)
          {
            ApplicationTools::displayResult(" shared between nodes", VectorTools::paste(vint, ","));
          }
      }
      modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, *tree, aliasFreqNames, globalParameters);
      model = 0;

      tl = make_shared<RNonHomogeneousTreeLikelihood>(*tree, *sites, modelSet, rDist, true, true);
    }
    else if (nhOpt == "general")
    {
      modelSet = LegacyPhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode, 0, testnh.getParams());
      if (modelSet->getNumberOfStates() > modelSet->alphabet().getSize())
      {
        // Markov-modulated Markov model!
        rDist = make_shared<ConstantRateDistribution>();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(testnh.getParams());
      }

      tl = make_shared<RNonHomogeneousTreeLikelihood>(*tree, *sites, modelSet, rDist, true, true);
    }
    else
      throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    tl->initialize();

    ApplicationTools::displayResult("Likelihood", TextTools::toString(tl->getValue(), 15));

    tree = make_unique<TreeTemplate<Node>>(tl->tree());

    // 2.2) Simulate
    distFile = ApplicationTools::getAFilePath("bootstrap.dist_file", testnh.getParams(), false, false);
    ApplicationTools::displayResult("Null distribution in file", distFile);
    if (nhOpt == "no")
      nhss = make_unique<NonHomogeneousSequenceSimulator>(model, rDist, tree);
    else
      nhss = make_unique<NonHomogeneousSequenceSimulator>(modelSet, rDist, tree);
    simulate(*nhss, nbSites, threshold, nbSim, nbTest, median, nbSignifPValue, medianPValue, distFile);
    ApplicationTools::displayResult("Global Bowker test number signif.", nbSignifPValue);
    ApplicationTools::displayResult("Bowker test median signif.", medianPValue);

    testnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return 0;
}
