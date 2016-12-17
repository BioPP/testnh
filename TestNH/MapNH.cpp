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

#include "MultinomialClustering.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Io/BppOSubstitutionModelFormat.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Stat/StatTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Stat/ContingencyTableTest.h>

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


void buildCountTree(
  const vector< vector<double> >& counts,
  const vector<int>& ids,
  Tree* cTree,
  size_t type)
{
  for (size_t i = 0; i < ids.size(); ++i)
  {
    if (cTree->hasFather(ids[i]))
    {
      cTree->setDistanceToFather(ids[i], counts[i][type]);
    }
  }
}


SubstitutionRegister* getSubstitutionRegister(const std::string& regTypeDesc, const SubstitutionModel* model)
{
  string regType = "";
  map<string, string> regArgs;
  KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);
  
  SubstitutionRegister* reg = 0;

  if (regType=="Combination")
  {
    VectorOfSubstitionRegisters* vreg= new VectorOfSubstitionRegisters(model);

    size_t i = 0;
    while (++i)
    {
      string regDesc = ApplicationTools::getStringParameter("reg" + TextTools::toString(i), regArgs, "", "", false, 1);
      if (regDesc=="")
        break;
      
      SubstitutionRegister* sreg=getSubstitutionRegister(regDesc, model);

      vreg->addRegister(sreg);
    }
    
    reg=vreg;
  }
  else if (regType == "All")
  {
    reg = new ComprehensiveSubstitutionRegister(model, false);
  }
  else if (regType == "Total")
  {
    reg = new TotalSubstitutionRegister(model);
  }    
  else if (regType == "Selected"){  
    string subsList = ApplicationTools::getStringParameter("substitution.list", regArgs, "All", "", true, false);
    reg = new SelectedSubstitutionRegister(model, subsList);  
  }
  else if (regType == "IntraAA")
  {
    if (AlphabetTools::isCodonAlphabet(model->getAlphabet()))
    {
      reg = new AAInteriorSubstitutionRegister(dynamic_cast<const CodonSubstitutionModel*>(model)); 
    }
    else
      throw Exception("Internal amino-acid categorization is only available for codon alphabet!");
  }
  else if (regType == "InterAA")
  {
    if (AlphabetTools::isCodonAlphabet(model->getAlphabet()))
    {
      reg = new AAExteriorSubstitutionRegister(dynamic_cast<const CodonSubstitutionModel*>(model)); 
    }
    else
      throw Exception("External amino-acid categorization is only available for codon alphabet!");
  }
  else if (regType == "GC")
  {
    if (AlphabetTools::isNucleicAlphabet(model->getAlphabet()))
      reg = new GCSubstitutionRegister(dynamic_cast<const NucleotideSubstitutionModel*>(model), false);
    else if (AlphabetTools::isCodonAlphabet(model->getAlphabet()))
      reg = new GCSynonymousSubstitutionRegister(dynamic_cast<const CodonSubstitutionModel*>(model));
    else
      throw Exception("GC categorization is only available for nucleotide or codon alphabets!");
  }
  else if (regType == "TsTv")
  {
    if (AlphabetTools::isNucleicAlphabet(model->getAlphabet()))
      reg = new TsTvSubstitutionRegister(dynamic_cast<const NucleotideSubstitutionModel*>(model));
    else if (AlphabetTools::isCodonAlphabet(model->getAlphabet()))
      reg = new TsTvSubstitutionRegister(dynamic_cast<const CodonSubstitutionModel*>(model));
    else
      throw Exception("TsTv categorization is only available for nucleotide or codon alphabet!");
  }
  else if (regType == "KrKc")
  {
    if (AlphabetTools::isProteicAlphabet(model->getAlphabet()))
      reg = new KrKcSubstitutionRegister(dynamic_cast<const ProteinSubstitutionModel*>(model));
    else
      throw Exception("KrKc categorization is only available for amino acid alphabet!");
  }
  else if (regType == "DnDs")
  {
    if (AlphabetTools::isCodonAlphabet(model->getAlphabet()))
    {
      reg = new DnDsSubstitutionRegister(dynamic_cast<const CodonSubstitutionModel*>(model), false);
    }
    else
      throw Exception("DnDs categorization is only available for codon alphabet!");
  }
  else
    throw Exception("Unsupported substitution categorization: " + regType);

  CategorySubstitutionRegister* csr=dynamic_cast<CategorySubstitutionRegister*>(reg);
  if (csr)
    csr->setStationarity(ApplicationTools::getBooleanParameter("stationarity", regArgs, true));
    
  return reg;
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     Map NH, version 1.1.1                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Modif. 17/12/11      *" << endl;
  cout << "*          L. Guéguen                       Last Modif. 24/11/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }

  try
  {
    BppApplication mapnh(args, argv, "MapNH");
    mapnh.startTimer();

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(mapnh.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", mapnh.getParams(), "Standard", "", true, 1);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, mapnh.getParams());
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, mapnh.getParams());
    delete allSites;

    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

    // Get the initial tree
    Tree* tree = PhylogeneticsApplicationTools::getTree(mapnh.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
    // Convert to NHX if input tree is newick or nexus?
    string treeIdOut = ApplicationTools::getAFilePath("output.tree_with_id.file", mapnh.getParams(), false, false, "none", 1);
    if (treeIdOut != "none")
    {
      Nhx nhx(true);
      nhx.write(*tree, treeIdOut);
    }

    //
    // Get substitution model and compute likelihood arrays
    //

    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", mapnh.getParams(), "no", "", true, 1);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    DRTreeLikelihood* drtl     = 0;
    SubstitutionModel* model    = 0;
    SubstitutionModelSet* modelSet = 0;
    DiscreteDistribution* rDist    = 0;

    if (nhOpt == "no")
    {
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, mapnh.getParams());
      if (model->getName() != "RE08")
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
      }
      
      drtl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, false, false);
    }
    else if (nhOpt == "one_per_branch")
    {
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, mapnh.getParams());
      if (model->getName() != "RE08")
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        // Markov-Modulated Markov Model...
        size_t n = (size_t)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
        // we should assume a rate distribution for the root also!!!
      }
      std::map<std::string, std::string> aliasFreqNames;
      FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, mapnh.getParams(), aliasFreqNames, rateFreqs);

      vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", mapnh.getParams(), ',', "");
      modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, aliasFreqNames, globalParameters);
      drtl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, false);
    }
    else if (nhOpt == "general")
    {
      modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, mapnh.getParams());
      if (modelSet->getModel(0)->getName() != "RE08")
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
      }
      drtl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, false);
    }
    else
      throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
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
      bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", mapnh.getParams(), false, "", true, 1);
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
    
    //
    // Initialize the parameters for the mapping:
    //

    string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "All", "", true, false);

    SubstitutionRegister* reg = getSubstitutionRegister(regTypeDesc, model ? model : modelSet->getModel(0));
        
    //Write categories:
    for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
      ApplicationTools::displayResult("  * Count type " + TextTools::toString(i + 1), reg->getTypeName(i + 1));

    // specific parameters to the null models
    string nullModelParams = ApplicationTools::getStringParameter("nullModelParams", mapnh.getParams(), "", "", false, 1);

    ParameterList nullParams;
    if (nullModelParams != "")
    {
      string modelName = "";
      map<string, string> npv;
      KeyvalTools::multipleKeyvals(nullModelParams, npv, ",", false);

      map<string, string>::iterator mi(npv.begin());
      while (mi != npv.end())
      {
        nullParams.addParameter(Parameter(mi->first, TextTools::toDouble(mi->second)));
        ApplicationTools::displayResult("null Parameter " + mi->first, mi->second);
        
        mi++;
      }
    }

 
    //
    // Performs mapping
    //

    vector<int> ids = drtl->getTree().getNodesId();
    ids.pop_back(); // remove root id.
    vector< vector<double> > counts;
    double thresholdSat = ApplicationTools::getDoubleParameter("count.max", mapnh.getParams(), -1, "", true, 1);
    if (thresholdSat > 0)
      ApplicationTools::displayResult("Saturation threshold used", thresholdSat);

    if (nullModelParams != "")
    {
      if (model)
      {
        unique_ptr<SubstitutionModel> nullModel(model->clone());
        
        ParameterList pl;
        const ParameterList pl0 = nullModel->getParameters();
        
        for (size_t i = 0; i < nullParams.size(); ++i)
        {
          vector<string> pn = pl0.getMatchingParameterNames(nullParams[i].getName());
          for (size_t j = 0; j < pn.size(); ++j)
          {
            pl.addParameter(Parameter(pn[j], nullParams[i].getValue()));
          }
        }

        nullModel->matchParametersValues(pl);
        
        counts = SubstitutionMappingTools::getNormalizedCountsPerBranch(*drtl, ids, model, nullModel.get(), *reg, true);
      }
      else
      {
        unique_ptr<SubstitutionModelSet> nullModelSet(modelSet->clone());
        ParameterList pl;
        const ParameterList pl0 = nullModelSet->getParameters();

        for (size_t i = 0; i < nullParams.size(); ++i)
        {
          vector<string> pn = pl0.getMatchingParameterNames(nullParams[i].getName());
          for (size_t j = 0; j < pn.size(); ++j)
          {
            pl.addParameter(Parameter(pn[j], nullParams[i].getValue()));
          }
        }

        nullModelSet->matchParametersValues(pl);

        counts = SubstitutionMappingTools::getNormalizedCountsPerBranch(*drtl, ids, modelSet, nullModelSet.get(), *reg, true);
      }
    }
    else
    {
      counts = SubstitutionMappingTools::getRelativeCountsPerBranch(*drtl, ids, model ? model : modelSet->getModel(0), *reg, thresholdSat);
    }
    
    vector<string> outputDesc = ApplicationTools::getVectorParameter<string>("output.counts", mapnh.getParams(), ',', "PerType(prefix=)");
    for (vector<string>::iterator it = outputDesc.begin(); it != outputDesc.end(); ++it) {
      string outputType;
      map<string, string> outputArgs;
      KeyvalTools::parseProcedure(*it, outputType, outputArgs);
      if (outputType == "PerType")
      {
        // Write count trees:
        string treePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_type_", "", true, 1);
        if (treePathPrefix != "none")
        {
          Newick newick;
          for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
          {
            string path = treePathPrefix + TextTools::toString(i + 1) + string(".dnd");
            ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
            Tree* cTree = tree->clone();
            buildCountTree(counts, ids, cTree, i);
            newick.write(*cTree, path);
            delete cTree;
          }
        }
      }
      else if (outputType == "PerBranchPerSite")
      {
        string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_branch_per_site.txt", "", true, 1);
        if (perSitenf != "none")
        {
          ApplicationTools::displayResult(string("Output counts (branch/site) to file"), perSitenf);
          SubstitutionMappingTools::outputTotalCountsPerBranchPerSite(perSitenf, *drtl, ids, model ? model : modelSet->getModel(0), *reg);
        }
      }
      else if (outputType == "PerSitePerType")
      {
        string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_type_per_site.txt", "", true, 1);
        if (perSitenf != "none")
        {
          ApplicationTools::displayResult(string("Output counts (type/site) to file"), perSitenf);
          SubstitutionMappingTools::outputTotalCountsPerTypePerSite(perSitenf, *drtl, ids, model ? model : modelSet->getModel(0), *reg);
        }
      }
      else if (outputType == "PerBranchPerSitePerType")
      {
        string tablePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_branch_per_site_per_type_", "", true, 1);
        if (tablePathPrefix != "none")
        {
          ApplicationTools::displayResult(string("Output counts (branch/site/type) to files"), tablePathPrefix + "*");
          SubstitutionMappingTools::outputIndividualCountsPerBranchPerSite(tablePathPrefix, *drtl, ids, model ? model : modelSet->getModel(0), *reg);
        }
      }
      else
      {
        throw Exception("Unknown output option: '" + outputType + "'");
      }
    }

    // Rounded counts
    vector< vector<size_t> > countsint;
    for (size_t i = 0; i < counts.size(); i++)
    {
      vector<size_t> countsi2;
      for (size_t j = 0; j < counts[i].size(); j++)
      {
        countsi2.push_back(static_cast<size_t>(floor( counts[i][j] + 0.5)));
      }
      countsint.push_back(countsi2);
    }


    // Global homogeneity test:
    bool testGlobal = ApplicationTools::getBooleanParameter("test.global", mapnh.getParams(), false, "", true, false);
    if (testGlobal)
    {
      vector< vector<size_t> > counts2 = countsint;

      // Check if some branches are 0:
      for (size_t i = counts2.size(); i > 0; --i)
      {
        if (VectorTools::sum(counts2[i - 1]) == 0)
        {
          ApplicationTools::displayResult("Remove branch with no substitution", ids[i - 1]);
          counts2.erase(counts2.begin() + static_cast<ptrdiff_t>(i - 1));
          // ids.erase(ids.begin() + i - 1);
        }
      }
      ApplicationTools::displayResult("Nb. of branches included in test", counts2.size());


      ContingencyTableTest test(counts2, 2000);
      ApplicationTools::displayResult("Global Chi2", test.getStatistic());
      ApplicationTools::displayResult("Global Chi2, p-value", test.getPValue());
      double pvalue = SimpleSubstitutionCountsComparison::multinomialTest(counts2);
      ApplicationTools::displayResult("Global heterogeneity test p-value", pvalue);
    }

    // Branch test!
    bool testBranch = ApplicationTools::getBooleanParameter("test.branch", mapnh.getParams(), false, "", true, 1);
    if (testBranch)
    {
      bool testNeighb = ApplicationTools::getBooleanParameter("test.branch.neighbor", mapnh.getParams(), true, "", true, 1);
      bool testNegBrL = ApplicationTools::getBooleanParameter("test.branch.negbrlen", mapnh.getParams(), false, "", true, 2);
      ApplicationTools::displayBooleanResult("Perform branch clustering", testBranch);
      ApplicationTools::displayBooleanResult("Cluster only neighbor nodes", testNeighb);
      ApplicationTools::displayBooleanResult("Allow len < 0 in clustering", testNegBrL);
      string autoClustDesc = ApplicationTools::getStringParameter("test.branch.auto_cluster", mapnh.getParams(), "Global(threshold=0)", "", true, 1);
      string autoClustName;
      map<string, string> autoClustParam;
      KeyvalTools::parseProcedure(autoClustDesc, autoClustName, autoClustParam);
      ApplicationTools::displayResult("Auto-clustering", autoClustName);
      unique_ptr<AutomaticGroupingCondition> autoClust;
      if (autoClustName == "None")
      {
        autoClust.reset(new NoAutomaticGroupingCondition());
      }
      else if (autoClustName == "Global")
      {
        size_t threshold = ApplicationTools::getParameter<size_t>("threshold", autoClustParam, 0, "", true, 1);
        ApplicationTools::displayResult("Auto-clutering threshold", threshold);
        CategorySubstitutionRegister* creg = dynamic_cast<CategorySubstitutionRegister*>(reg);
        vector<size_t> toIgnore;
        if (creg && creg->allowWithin())
        {
          size_t n = creg->getNumberOfCategories();
          for (size_t i = 0; i < n; ++i)
          {
            toIgnore.push_back(n * (n - 1) + i);
          }
        }
        autoClust.reset(new SumCountsAutomaticGroupingCondition(threshold, toIgnore));
      }
      else if (autoClustName == "Marginal")
      {
        size_t threshold = ApplicationTools::getParameter<size_t>("threshold", autoClustParam, 0, "", true, 1);
        ApplicationTools::displayResult("Auto-clutering threshold", threshold);
        autoClust.reset(new AnyCountAutomaticGroupingCondition(threshold));
      }
      else
      {
        throw Exception("Unknown automatic clustering option: " + autoClustName);
      }

      // ChiClustering htest(counts, ids, true);
      MultinomialClustering htest(countsint, ids, drtl->getTree(), *autoClust, testNeighb, testNegBrL, true);
      ApplicationTools::displayResult("P-value at root node", *(htest.getPValues().rbegin()));
      ApplicationTools::displayResult("Number of tests performed", htest.getPValues().size());
      TreeTemplate<Node>* htree = htest.getTree();
      Newick newick;
      string clusterTreeOut = ApplicationTools::getAFilePath("output.cluster_tree.file", mapnh.getParams(), false, false, "", true, "clusters.dnd", 1);
      ApplicationTools::displayResult("Output cluster tree to", clusterTreeOut);
      newick.write(*htree, clusterTreeOut);
      delete htree;
    }

    // Cleaning up:
    delete alphabet;
    delete sites;
    delete tree;
    if (modelSet)
      delete modelSet;
    else
      delete model;
    delete rDist;
    delete reg;
    mapnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return 0;
}
