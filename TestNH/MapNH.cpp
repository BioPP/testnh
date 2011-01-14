//
// File: MapNH.cpp
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

#include "ChiClustering.h"
#include "MultinomialClustering.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

// From bpp-seq:
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/Mapping.all>
#include <Bpp/Phyl/Simulation.all>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Phyl/Model/JCprot.h>
#include <Bpp/Phyl/Model/CodonNeutralReversibleSubstitutionModel.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Stat/ContingencyTableTest.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "mapnh parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the package manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

//Returns the log probability... but we should use Chisquare instead!
//double dMultinom(const vector<double>& counts, const vector<double>& probs) {
//  unsigned int n = VectorTools::sum(counts);
//  unsigned int l = log(NumTools::fact(n));
//  for (size_t i = 0; i < counts.size(); ++i)
//    l += counts[i] * log(probs[i]) - log(NumTools::fact(counts[i]));
//  return l;
//}

vector< vector<unsigned int> > getCountsPerBranch(
    DRTreeLikelihood& drtl,
    const vector<int>& ids,
    SubstitutionCount& count)
{
  ProbabilisticSubstitutionMapping* mapping = SubstitutionMappingTools::computeSubstitutionVectors(drtl, count, false);
  vector< vector<unsigned int> > counts(ids.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    vector<double> countsf = SubstitutionMappingTools::computeSumForBranch(*mapping, mapping->getNodeIndex(ids[i]));
    counts[i].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j)
      counts[i][j] = static_cast<unsigned int>(floor(countsf[j] + 0.5)); //Round counts
  }
  delete mapping;
  return counts;
}

void buildDnAndDsTrees(const vector< vector<unsigned int> > counts,
                       const vector<int>& ids,
                       Tree* dnTree, 
                       Tree* dsTree)
{
  for (size_t i = 0; i < ids.size(); ++i) {

    if (dnTree->hasFather(i)){
      dnTree->setDistanceToFather(i, counts[i][1]);
      dsTree->setDistanceToFather(i, counts[i][0]);
    }
  }
}


int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     Map NH, version 0.1.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 20/03/08 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  BppApplication mapnh(args, argv, "MapNH");
  mapnh.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(mapnh.getParams(), "", false);
  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, mapnh.getParams());
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, mapnh.getParams());
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  //Get the initial tree
  Tree* tree = PhylogeneticsApplicationTools::getTree(mapnh.getParams());
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  //Convert to NHX if input tree is newick or nexus?
  string treeIdOut = ApplicationTools::getAFilePath("output.tree_with_id.file", mapnh.getParams(), false, false);
  Nhx nhx(true);
  nhx.write(*tree, treeIdOut);

  //Perform the mapping:
  SubstitutionRegister* reg = 0;
  string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "all", "", true, false);
  string regType = "";
  map<string, string> regArgs;
  KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);
  auto_ptr<GeneticCode> geneticCode;
  if (regType == "All")
    reg = new ExhaustiveSubstitutionRegister(alphabet);
  else if (regType == "GC") {
    if (AlphabetTools::isNucleicAlphabet(alphabet)) {
      bool includeZero = ApplicationTools::getBooleanParameter("inc_zeros", regArgs, false);
      reg = new GCSubstitutionRegister(dynamic_cast<NucleicAlphabet*>(alphabet), includeZero);
    } else
      throw Exception("GC categorization is only available for nucleotide alphabet!");
  } else if (regType == "TsTv") {
    if (AlphabetTools::isNucleicAlphabet(alphabet))
      reg = new TsTvSubstitutionRegister(dynamic_cast<NucleicAlphabet*>(alphabet));
    else
      throw Exception("TsTv categorization is only available for nucleotide alphabet!");
  } else if (regType == "DnDs") {
    if (AlphabetTools::isCodonAlphabet(alphabet)) {
      string code = regArgs["code"];
      if (TextTools::isEmpty(code)) {
        code = "Standard";
        ApplicationTools::displayWarning("No genetic code provided, standard code used.");
      }
      geneticCode.reset(SequenceApplicationTools::getGeneticCode(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet(), code));
      reg = new DnDsSubstitutionRegister(geneticCode.get());
    } else
      throw Exception("DnDs categorization is only available for nucleic alphabet!");
  } else
    throw Exception("Unsupported substitution categorization: " + regType);

  //Now perform mapping using a JC model:
  SubstitutionModel* model = 0;
  if (AlphabetTools::isNucleicAlphabet(alphabet)) {
    model = new JCnuc(dynamic_cast<NucleicAlphabet*>(alphabet));
  } else if (AlphabetTools::isProteicAlphabet(alphabet)) {
    model = new JCprot(dynamic_cast<ProteicAlphabet*>(alphabet));
  } else if (AlphabetTools::isCodonAlphabet(alphabet)) {
//    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, mapnh.getParams());
    model = new CodonNeutralReversibleSubstitutionModel(
        dynamic_cast<const CodonAlphabet*>(geneticCode->getSourceAlphabet()),
        new JCnuc(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet()));
  }
  else
    throw Exception("Unsupported alphabet!");

  DiscreteDistribution* rDist = new ConstantDistribution(1., true);
  SubstitutionCount* count = new SimpleSubstitutionCount(reg);

  DRHomogeneousTreeLikelihood drtl(*tree, *sites, model, rDist, false, false);
  drtl.initialize();
  
  //Optimization of parameters:
  //PhylogeneticsApplicationTools::optimizeParameters(&drtl, drtl.getParameters(), mapnh.getParams(), "", true, true);
 
  vector<int> ids = drtl.getTree().getNodesId();
  ids.pop_back(); //remove root id. 
  vector< vector<unsigned int> > counts = getCountsPerBranch(drtl, ids, *count);
    std::cout <<"Size 1: "<<counts.size() << " Size 2: "<<counts[0].size()<<std::endl;
  if (regType == "DnDs") {
    Tree* dnTree = tree->clone();
    Tree* dsTree = tree->clone();
    buildDnAndDsTrees(counts, ids, dnTree, dsTree);
    Newick newick;
    string dnTreeOut = ApplicationTools::getAFilePath("output.dn_tree.file", mapnh.getParams(), false, false);
    string dsTreeOut = ApplicationTools::getAFilePath("output.ds_tree.file", mapnh.getParams(), false, false);
    ApplicationTools::displayResult("Output Dn tree to", dnTreeOut);
    ApplicationTools::displayResult("Output Ds tree to", dsTreeOut);
    newick.write(*dnTree, dnTreeOut);
    newick.write(*dsTree, dsTreeOut);
    delete dnTree;
    delete dsTree;
  }
  //Global homogeneity test:
  bool testGlobal = ApplicationTools::getBooleanParameter("test.global", mapnh.getParams(), true, "", true, false);
  if (testGlobal) {
    vector< vector<unsigned int> > counts2 = counts; 
    //Check if some branches are 0:
    for (size_t i = counts2.size(); i > 0; --i) {
      if (VectorTools::sum(counts2[i - 1]) == 0) {
        ApplicationTools::displayResult("Remove branch with no substitution", ids[i - 1]);
        counts.erase(counts2.begin() + i - 1);
        //ids.erase(ids.begin() + i - 1);
      }
    }
    ApplicationTools::displayResult("Nb. of branches included in test", counts2.size());

    ContingencyTableTest test(counts2, 2000);
    ApplicationTools::displayResult("Global Chi2", test.getStatistic());
    ApplicationTools::displayResult("Global Chi2, p-value", test.getPValue());
  }

  //Branch test!
  bool testBranch = ApplicationTools::getBooleanParameter("test.branch", mapnh.getParams(), false, "", true, false);
  bool testNeighb = ApplicationTools::getBooleanParameter("test.branch.neighbor", mapnh.getParams(), true, "", true, false);
  bool testNegBrL = ApplicationTools::getBooleanParameter("test.branch.negbrlen", mapnh.getParams(), false, "", true, false);
  ApplicationTools::displayBooleanResult("Perform branch clustering", testBranch);
  ApplicationTools::displayBooleanResult("Cluster only neighbor nodes", testNeighb);
  ApplicationTools::displayBooleanResult("Allow len < 0 in clustering", testNegBrL);
  if (testBranch) {
    //ChiClustering htest(counts, ids, true);
    MultinomialClustering htest(counts, ids, drtl.getTree(), testNeighb, testNegBrL, true);
    TreeTemplate<Node>* htree = htest.getTree();
    Newick newick;
    string clusterTreeOut = ApplicationTools::getAFilePath("output.cluster_tree.file", mapnh.getParams(), false, false);
    ApplicationTools::displayResult("Output cluster tree to", clusterTreeOut);
    newick.write(*htree, clusterTreeOut);
    delete htree;
  }

  //Cleaning up:
  delete alphabet;
  delete sites;
  delete tree;
  delete model;
  delete rDist;
  delete count;
  mapnh.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}
