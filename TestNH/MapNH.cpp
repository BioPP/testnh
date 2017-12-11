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

//#include "MultinomialClustering.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/OneProcessSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/SingleProcessSubstitutionMapping.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Text/KeyvalTools.h>

#include "bppTools.h"

using namespace bpp;

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     Map NH, version 2                          *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Modif. 17/12/11      *" << endl;
  cout << "*          L. Guéguen                       Last Modif. 28/11/17 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    bppTools::help("mapnh");
    exit(0);
  }
  
  try
  {
  BppApplication mapnh(args, argv, "MapNH");
  mapnh.startTimer();
  std::map<std::string, std::string> unparsedparams;

  /*********************************/
  /* get Basic objects */
  /*********************************/
  
  unique_ptr<Alphabet> alphabet(bppTools::getAlphabet(mapnh.getParams()));
  unique_ptr<GeneticCode> gCode(bppTools::getGeneticCode(mapnh.getParams(),alphabet.get()));

  map<size_t, AlignedValuesContainer*> mSites=bppTools::getAlignmentsMap(mapnh.getParams(), alphabet.get(), true);

  unique_ptr<SubstitutionProcessCollection> SP(bppTools::getCollection(mapnh.getParams(), alphabet.get(), gCode.get(), mSites, unparsedparams));
                                               
  std::map<size_t, SequenceEvolution*> mProc=bppTools::getProcesses(mapnh.getParams(), *SP, unparsedparams);
  
  unique_ptr<PhyloLikelihoodContainer> plc(bppTools::getPhyloLikelihoods(mapnh.getParams(), mProc, *SP, mSites));

  if (!plc->hasPhyloLikelihood(1))
    throw Exception("Missing first phyloLikelihood.");

  PhyloLikelihood* pl=(*plc)[1];

  /// For now restriction to OneSequencePhyloLikelihoods:

  OneProcessSequencePhyloLikelihood* opspl= dynamic_cast<OneProcessSequencePhyloLikelihood*>(pl);
  SingleProcessPhyloLikelihood* sppl= dynamic_cast<SingleProcessPhyloLikelihood*>(pl);
  
  if (opspl==NULL && sppl==NULL)
    throw Exception("Mapping only for simple phylo.");
  
  // if (model->getName() != "RE08")
  //   SiteContainerTools::changeGapsToUnknownCharacters(*sites);
  
  
  //////////////////////////////////
  // Set register and initialize the parameters for the mapping:
  //////////////

  // Only One alphabet state map
  const StateMap& stmap=opspl?opspl->getSubstitutionProcess().getStateMap():sppl->getSubstitutionProcess().getStateMap();

  string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "All", "", true, false);

  unique_ptr<SubstitutionRegister> reg(PhylogeneticsApplicationTools::getSubstitutionRegister(regTypeDesc, stmap, gCode.get()));

  double thresholdSat = ApplicationTools::getDoubleParameter("count.max", mapnh.getParams(), -1, "", true, 1);
  if (thresholdSat > 0)
    ApplicationTools::displayResult("Saturation threshold used", thresholdSat);

  unique_ptr<PhyloSubstitutionMapping> psm;

  if (opspl)
    psm.reset(new OneProcessSubstitutionMapping(*opspl, *reg, thresholdSat));
  else
    psm.reset(new SingleProcessSubstitutionMapping(*sppl, *reg, thresholdSat));

  //Write categories:
  for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
    ApplicationTools::displayResult("  * Count type " + TextTools::toString(i + 1), reg->getTypeName(i + 1));

  
  // specific parameters to the null models
  string nullProcessParams = ApplicationTools::getStringParameter("nullProcessParams", mapnh.getParams(), "", "", false, 1);
  
  ParameterList nullParams;
  if (nullProcessParams != "")
  {
    string modelName = "";
    map<string, string> npv;
    KeyvalTools::multipleKeyvals(nullProcessParams, npv, ",", false);
    
    map<string, string>::iterator mi(npv.begin());
    while (mi != npv.end())
    {
      nullParams.addParameter(Parameter(mi->first, TextTools::toDouble(mi->second)));
      ApplicationTools::displayResult("null Parameter " + mi->first, mi->second);
      
      mi++;
    }

    psm->computeNormalizations(nullParams);
  }

  
  string outputDesc = ApplicationTools::getStringParameter("output.counts", mapnh.getParams(), "PerType(prefix=)");

  uint siteSize=AlphabetTools::isWordAlphabet(alphabet.get())?dynamic_cast<const CoreWordAlphabet*>(alphabet.get())->getLength():1;

  Vuint ids=psm->getCounts().getAllEdgesIndexes();
    
  ProbabilisticSubstitutionMapping* pCounts;
    
  if (nullProcessParams != "" && outputDesc.find("Branch")!=string::npos)
    pCounts=SubstitutionMappingTools::computeNormalizedCounts(
      &psm->getCounts(), &psm->getNormalizations(), true,  siteSize);
  else
    pCounts=&psm->getCounts();
  
  ////////////////////////////////////////////
  //// OUTPUT
  ////////////////////////////////////////////

    
  string outputType;
  map<string, string> outputArgs;
  KeyvalTools::parseProcedure(outputDesc, outputType, outputArgs);
    
  size_t outputNum=0;
  
  if (outputType.find("Type")!=string::npos)
    outputNum+=1;
  if (outputType.find("Site")!=string::npos)
    outputNum+=2;
  if (outputType.find("Branch")!=string::npos)
    outputNum+=4;
  
  switch(outputNum)
  {
  case 1:
  case 4:
  case 5:
    {
      // Write count trees:
      string treePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_type_", "", true, 1);
      
      Newick newick;
      for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
      {
        string name=reg->getTypeName(i+1);
        if (name=="")
          name=TextTools::toString(i + 1);
        
        unique_ptr<PhyloTree> pt(SubstitutionMappingTools::getTreeForType(*pCounts,i));
        
        string path = treePathPrefix + name + string(".dnd");
        ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
        newick.write(*pt, path);
      }
    }
    break;
  case 6:
    {
      string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_site_per_branch.txt", "", true, 1);
      
      ApplicationTools::displayResult(string("Output counts (branch/site) to file"), perSitenf);
      
      VVdouble counts(SubstitutionMappingTools::getCountsPerSitePerBranch(*pCounts));
      SubstitutionMappingTools::outputPerSitePerBranch(perSitenf, ids, counts);
    }
    break;
  case 2:
  case 3:
    {
      string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_site_per_type.txt", "", true, 1);
      
      ApplicationTools::displayResult(string("Output counts (site/type) to file"), perSitenf);
      
      VVdouble counts;
      
      if (nullProcessParams!="")
        counts=SubstitutionMappingTools::getCountsPerSitePerType(psm->getCounts(), psm->getNormalizations(), true, siteSize);
      else
        counts=SubstitutionMappingTools::getCountsPerSitePerType(psm->getCounts());
      
      SubstitutionMappingTools::outputPerSitePerType(perSitenf, *reg, counts);          
    }
    break;
  case 7:
    {
      string tablePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_site_per_branch_per_type_", "", true, 1);
      
      ApplicationTools::displayResult(string("Output counts (site/branch/type) to files"), tablePathPrefix + "*");
      
      VVVdouble counts3(SubstitutionMappingTools::getCountsPerSitePerBranchPerType(*pCounts));
      SubstitutionMappingTools::outputPerSitePerBranchPerType(tablePathPrefix, ids, *reg, counts3);
    }
    break;
  default:
    throw Exception("Unknown output option: '" + outputType + "'");
  }
  
  //////////////////////////////////////
  /// HOMOGENEITY TESTS
  /////////////////////////////////////
  
  // bool testGlobal = ApplicationTools::getBooleanParameter("test.global", mapnh.getParams(), false, "", true, false);
  // bool testBranch = ApplicationTools::getBooleanParameter("test.branch", mapnh.getParams(), false, "", true, false);

    
    // // Rounded counts
    // vector< vector<size_t> > countsint;
    // for (size_t i = 0; i < counts.size(); i++)
    // {
    //   vector<size_t> countsi2;
    //   for (size_t j = 0; j < counts[i].size(); j++)
    //   {
    //     countsi2.push_back(static_cast<size_t>(floor( counts[i][j] + 0.5)));
    //   }
    //   countsint.push_back(countsi2);
    // }


    // // Global homogeneity test:
    // if (testGlobal)
    // {
    //   vector< vector<size_t> > counts2 = countsint;

    //   // Check if some branches are 0:
    //   for (size_t i = counts2.size(); i > 0; --i)
    //   {
    //     if (VectorTools::sum(counts2[i - 1]) == 0)
    //     {
    //       ApplicationTools::displayResult("Remove branch with no substitution", ids[i - 1]);
    //       counts2.erase(counts2.begin() + static_cast<ptrdiff_t>(i - 1));
    //       // ids.erase(ids.begin() + i - 1);
    //     }
    //   }
    //   ApplicationTools::displayResult("Nb. of branches included in test", counts2.size());


    //   ContingencyTableTest test(counts2, 2000);
    //   ApplicationTools::displayResult("Global Chi2", test.getStatistic());
    //   ApplicationTools::displayResult("Global Chi2, p-value", test.getPValue());
    //   double pvalue = SimpleSubstitutionCountsComparison::multinomialTest(counts2);
    //   ApplicationTools::displayResult("Global heterogeneity test p-value", pvalue);
    // }

    // // Branch test!
    // if (testBranch)
    // {
    //   bool testNeighb = ApplicationTools::getBooleanParameter("test.branch.neighbor", mapnh.getParams(), true, "", true, 1);
    //   bool testNegBrL = ApplicationTools::getBooleanParameter("test.branch.negbrlen", mapnh.getParams(), false, "", true, 2);
    //   ApplicationTools::displayBooleanResult("Perform branch clustering", testBranch);
    //   ApplicationTools::displayBooleanResult("Cluster only neighbor nodes", testNeighb);
    //   ApplicationTools::displayBooleanResult("Allow len < 0 in clustering", testNegBrL);
    //   string autoClustDesc = ApplicationTools::getStringParameter("test.branch.auto_cluster", mapnh.getParams(), "Global(threshold=0)", "", true, 1);
    //   string autoClustName;
    //   map<string, string> autoClustParam;
    //   KeyvalTools::parseProcedure(autoClustDesc, autoClustName, autoClustParam);
    //   ApplicationTools::displayResult("Auto-clustering", autoClustName);
    //   unique_ptr<AutomaticGroupingCondition> autoClust;
    //   if (autoClustName == "None")
    //   {
    //     autoClust.reset(new NoAutomaticGroupingCondition());
    //   }
    //   else if (autoClustName == "Global")
    //   {
    //     size_t threshold = ApplicationTools::getParameter<size_t>("threshold", autoClustParam, 0, "", true, 1);
    //     ApplicationTools::displayResult("Auto-clutering threshold", threshold);
    //     CategorySubstitutionRegister* creg = dynamic_cast<CategorySubstitutionRegister*>(reg.get());
    //     vector<size_t> toIgnore;
    //     if (creg && creg->allowWithin())
    //     {
    //       size_t n = creg->getNumberOfCategories();
    //       for (size_t i = 0; i < n; ++i)
    //       {
    //         toIgnore.push_back(n * (n - 1) + i);
    //       }
    //     }
    //     autoClust.reset(new SumCountsAutomaticGroupingCondition(threshold, toIgnore));
    //   }
    //   else if (autoClustName == "Marginal")
    //   {
    //     size_t threshold = ApplicationTools::getParameter<size_t>("threshold", autoClustParam, 0, "", true, 1);
    //     ApplicationTools::displayResult("Auto-clutering threshold", threshold);
    //     autoClust.reset(new AnyCountAutomaticGroupingCondition(threshold));
    //   }
    //   else
    //   {
    //     throw Exception("Unknown automatic clustering option: " + autoClustName);
    //   }

    //   // ChiClustering htest(counts, ids, true);
    //   MultinomialClustering htest(countsint, ids, drtl->getTree(), *autoClust, testNeighb, testNegBrL, true);
    //   ApplicationTools::displayResult("P-value at root node", *(htest.getPValues().rbegin()));
    //   ApplicationTools::displayResult("Number of tests performed", htest.getPValues().size());
    //   TreeTemplate<Node>* htree = htest.getTree();
    //   Newick newick;
    //   string clusterTreeOut = ApplicationTools::getAFilePath("output.cluster_tree.file", mapnh.getParams(), false, false, "", true, "clusters.dnd", 1);
    //   ApplicationTools::displayResult("Output cluster tree to", clusterTreeOut);
    //   newick.write(*htree, clusterTreeOut);
    //   delete htree;
    // }
  
  // Cleaning up:
  for (auto it : mSites)
    delete it.second;
  
  mapnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }
  
  return 0;
}
