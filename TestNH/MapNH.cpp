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
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
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
#include <Bpp/Phyl/Io/BppOSubstitutionModelFormat.h>
#include <Bpp/Phyl/Model.all>
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

/*
 * Returns the counts on each branch.
 *
 */

vector< vector<double> > getCountsPerBranch(
                                            DRTreeLikelihood& drtl,
                                            const vector<int>& ids,
                                            SubstitutionModel* model,
                                            const SubstitutionRegister& reg,
                                            double threshold = -1)
{
  auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));

  auto_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(drtl, *count, false));

  vector< vector<double> > counts(ids.size());
  size_t nbSites = mapping->getNumberOfSites();
  size_t nbTypes = mapping->getNumberOfSubstitutionTypes();
  for (size_t k = 0; k < ids.size(); ++k) {
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    size_t nbIgnored = 0;
    bool error = false;
    for (size_t i = 0; !error && i < nbSites; ++i) {
      double s = 0;
      for (size_t t = 0; t < nbTypes; ++t) {
        tmp[t] = (*mapping)(k, i, t);
        error = isnan(tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0) {
        if (s <= threshold)
          countsf += tmp;
        else {
          nbIgnored++;
        }
      } else {
        countsf += tmp;
      }
    }

  ERROR:
    if (error) {
      //We do nothing. This happens for small branches.
      ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (size_t t = 0; t < nbTypes; ++t)
        countsf[t] = 0;
    } else {
      if (nbIgnored > 0) {
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
    }
    
    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j) {
      counts[k][j] = countsf[j]; 
    }
  }
  return counts;
}

/*
 * Returns the normalization factors due to the null model on each
 * branch, for each register
 *
 */

vector< vector<double> > getNormalizationsPerBranch(
                                                    DRTreeLikelihood& drtl,
                                                    const vector<int>& ids,
                                                    const SubstitutionModel* nullModel,
                                                    const SubstitutionRegister& reg)
{
  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = nullModel->getAlphabet()->getSize();
  size_t nbSites = drtl.getNumberOfSites();
  
  // compute the AlphabetIndex for each substitutionType
  vector<UserAlphabetIndex1 > usai(nbTypes, UserAlphabetIndex1(nullModel->getAlphabet()));

  for (size_t i=0; i<nbStates; i++)
    for (size_t j=0; j<nbStates; j++)
      if (i!=j)
        {
          size_t nbt=reg.getType((int)i,(int)j);
          if (nbt!=0)
            usai[nbt-1].setIndex((int)i,usai[nbt-1].getIndex((int)i)+nullModel->Qij((int)i,(int)j)); 
        }
  
  // compute the normalization for each substitutionType
  vector< vector<double> > rewards(ids.size());

  for (size_t k = 0; k < ids.size(); ++k)
    rewards[k].resize(nbTypes);

  for (size_t nbt=0; nbt< nbTypes; nbt++)
    {
      auto_ptr<Reward> reward(new DecompositionReward(nullModel, &usai[nbt]));

      auto_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(drtl, *reward, false));
      
      for (size_t k = 0; k < ids.size(); ++k) {
        double s = 0;
        for (size_t i = 0; i < nbSites; ++i) {
          double tmp = (*mapping)(k, i);
          if (isnan(tmp))
            {
              ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", reward for type " + reg.getTypeName(nbt) + " could not be computed.");
              s=0;
              break;
            }
          s += tmp;
        }
        rewards[k][nbt]=s;
      }
      reward.release();
      mapping.release();
    }
  return rewards;
}

/* Returns counts normalized by a null model */

vector< vector<double> >  getNormalizedCountsPerBranch(
                                                       DRTreeLikelihood& drtl,
                                                       const vector<int>& ids,
                                                       SubstitutionModel* model,
                                                       SubstitutionModel* nullModel,
                                                       const SubstitutionRegister& reg,
                                                       double threshold = -1)
{
  vector< vector<double> > counts=getCountsPerBranch(drtl, ids, model, reg, threshold);
  vector< vector<double> > factors=getNormalizationsPerBranch(drtl, ids, nullModel, reg);

  size_t nbTypes = counts[0].size();
    
  for (size_t k = 0; k < ids.size(); ++k) {
    for (size_t t = 0; t < nbTypes; ++t) {
      counts[k][t] /= factors[k][t];
    }
  }
    
  return counts;
}

/*
 * Mgmt of the non-stationarity to get counts relative to the
 * frequency of the states is case of non-stationarity.
 *
 */

vector< vector<double> >  getRelativeCountsPerBranch(
                                                     DRTreeLikelihood& drtl,
                                                     const vector<int>& ids,
                                                     SubstitutionModel* model,
                                                     const SubstitutionRegister& reg,
                                                     bool stationarity = true,
                                                     double threshold = -1)
{
  vector< vector<double> > counts=getCountsPerBranch(drtl, ids, model, reg, threshold);

  const CategorySubstitutionRegister* creg;  
  if (!stationarity) {
    try {
      creg = &dynamic_cast<const CategorySubstitutionRegister&>(reg);
    } catch (Exception& ex) {
      throw Exception("The stationarity option can only be used with a category substitution register.");
    }
    
    size_t nbTypes = counts[0].size();

    for (size_t k = 0; k < ids.size(); ++k) {
      vector<double> freqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(drtl, ids[k]);
      //Compute frequencies for types:
      vector<double> freqsTypes(creg->getNumberOfCategories());
      for (size_t i = 0; i < freqs.size(); ++i) {
        size_t c = creg->getCategory(static_cast<int>(i));
        freqsTypes[c - 1] += freqs[i];
      }
      
      //We devide the counts by the frequencies and rescale:
      double s = VectorTools::sum(counts[k]);
      for (size_t t = 0; t < nbTypes; ++t) {
        counts[k][t] /= freqsTypes[creg->getCategoryFrom(t + 1) - 1];
      }
       
      double s2 = VectorTools::sum(counts[k]);
      //Scale:
      counts[k] = (counts[k] / s2) * s;
    }
  }
  
  return counts;
}

void outputTotalCountsPerBranchPerSite(
                                       string& filename,
                                       DRTreeLikelihood& drtl,
                                       const vector<int>& ids,
                                       SubstitutionModel* model,
                                       const SubstitutionRegister& reg)
{
  auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));
  auto_ptr<ProbabilisticSubstitutionMapping> smap(SubstitutionMappingTools::computeSubstitutionVectors(drtl, *count, false));

  ofstream file;
  file.open(filename.c_str());

  size_t nbSites = smap->getNumberOfSites();
  size_t nbBr = ids.size();

  vector<int> sdi(nbBr);  //reverse of ids
  for (size_t i = 0; i < nbBr; ++i) {
    for (size_t j = 0; j < nbBr; ++j) {
      if (ids[j] == static_cast<int>(i)) {
        sdi[i] = static_cast<int>(j);
        break;
      }
    }
  }

  file << "sites";
  for (size_t i = 0; i < nbBr; ++i)
    file << "\t" << i ;
  file << endl;
  
  for (size_t k = 0; k < nbSites; ++k) {
    vector<double> countsf = SubstitutionMappingTools::computeTotalSubstitutionVectorForSite(*smap, k);
    file << k;
    for (size_t i = 0; i < nbBr; ++i)
      file << "\t" << countsf[sdi[i]];
    file << endl;
  }
  file.close();
}


void buildCountTree(
                    const vector< vector<double> >& counts,
                    const vector<int>& ids,
                    Tree* cTree, 
                    size_t type)
{
  for (size_t i = 0; i < ids.size(); ++i) {
    if (cTree->hasFather(ids[i])) {
      cTree->setDistanceToFather(ids[i], counts[i][type]);
    }
  }
}


int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     Map NH, version 0.1.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 17/12/11 *" << endl;
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
    if (treeIdOut != "none") {
      Nhx nhx(true);
      nhx.write(*tree, treeIdOut);
    }

    //Initialize the parameters for the mapping:
    SubstitutionRegister* reg = 0;
    string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "All", "", true, false);
    string regType = "";
    map<string, string> regArgs;
    KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);
    auto_ptr<GeneticCode> geneticCode;
    bool stationarity = true;
    if (regType == "All") {
      stationarity = ApplicationTools::getBooleanParameter("stationarity", regArgs, true);
      reg = new ComprehensiveSubstitutionRegister(alphabet, false);
    }
    else if (regType == "GC") {
      stationarity = ApplicationTools::getBooleanParameter("stationarity", regArgs, true);
      if (AlphabetTools::isNucleicAlphabet(alphabet)) 
        reg = new GCSubstitutionRegister(dynamic_cast<NucleicAlphabet*>(alphabet), false);
      else 
        if (AlphabetTools::isCodonAlphabet(alphabet)) {
          string code = regArgs["code"];
          if (TextTools::isEmpty(code)) {
            code = "Standard";
            ApplicationTools::displayWarning("No genetic code provided, standard code used.");
          }
          geneticCode.reset(SequenceApplicationTools::getGeneticCode(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet(), code));
          reg = new GCSynonymousSubstitutionRegister(geneticCode.get());
        }
        else 
          throw Exception("GC categorization is only available for nucleotide or codon alphabets!");
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
        reg = new DnDsSubstitutionRegister(geneticCode.get(), false);
      } else
        throw Exception("DnDs categorization is only available for nucleic alphabet!");
    } else
      throw Exception("Unsupported substitution categorization: " + regType);


    // specific initialization of the null model
    string nullModelDescription = ApplicationTools::getStringParameter("model0", mapnh.getParams(), "");

    SubstitutionModel* nullModel=0;
    if (nullModelDescription!="")
 {
        map<string, string> unparsedParameterValues;

        BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, true);
        nullModel = bIO.read(alphabet, nullModelDescription);
      }


    
    //Now compute likelihood arrays
 
    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", mapnh.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    DRTreeLikelihood     *drtl     = 0;
    SubstitutionModel    *model    = 0;
    SubstitutionModelSet *modelSet = 0;
    DiscreteDistribution *rDist    = 0;

    if (nhOpt == "no")
      { 
        if (ApplicationTools::parameterExists("model", mapnh.getParams())) {
          model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, mapnh.getParams());
          if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
          if(model->getNumberOfStates() > model->getAlphabet()->getSize())
            {
              //Markov-modulated Markov model!
              rDist = new ConstantRateDistribution();
            }
          else
            {
              rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
            }
        } else {
          if (AlphabetTools::isNucleicAlphabet(alphabet)) {
            model = new JCnuc(dynamic_cast<NucleicAlphabet*>(alphabet));
          } else if (AlphabetTools::isProteicAlphabet(alphabet)) {
            model = new JCprot(dynamic_cast<ProteicAlphabet*>(alphabet));
          } else if (AlphabetTools::isCodonAlphabet(alphabet)) {
            model = new CodonRateSubstitutionModel(
                                                   dynamic_cast<const CodonAlphabet*>(geneticCode->getSourceAlphabet()),
                                                   new JCnuc(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet()));
          } else
            throw Exception("Unsupported alphabet!");
          rDist = new ConstantRateDistribution();
        }
        drtl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, false, false);
      }
    else if (nhOpt == "one_per_branch")
      {
        model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, mapnh.getParams());
        if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (model->getNumberOfStates() > model->getAlphabet()->getSize())
          {
            //Markov-modulated Markov model!
            rDist = new ConstantRateDistribution();
          }
        else
          {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
          }
        vector<double> rateFreqs;
        if (model->getNumberOfStates() != alphabet->getSize())
          {
            //Markov-Modulated Markov Model...
            size_t n =(size_t)(model->getNumberOfStates() / alphabet->getSize());
            rateFreqs = vector<double>(n, 1./static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
            // we should assume a rate distribution for the root also!!!  
          }
        FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, mapnh.getParams(), rateFreqs);
        vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", mapnh.getParams(), ',', "");
        modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
        model = 0;
        drtl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, false);
      }
    else if (nhOpt == "general")
      {
        modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, mapnh.getParams());
        if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
          {
            //Markov-modulated Markov model!
            rDist = new ConstantDistribution(1.);
          }
        else
          {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
          }
        drtl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, false);
      }
    else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    drtl->initialize();


    ////
    // Performs the mapping

    vector<int> ids = drtl->getTree().getNodesId();
    ids.pop_back(); //remove root id.
    vector< vector<double> > counts;
    double thresholdSat = ApplicationTools::getDoubleParameter("count.max", mapnh.getParams(), -1);
    if (thresholdSat > 0)
      ApplicationTools::displayResult("Saturation threshold used", thresholdSat);

    if (nullModel)
      counts = getNormalizedCountsPerBranch(*drtl, ids, model ? model : modelSet->getModel(0), nullModel, *reg, thresholdSat);
    else 
      counts = getRelativeCountsPerBranch(*drtl, ids, model ? model : modelSet->getModel(0), *reg, stationarity, thresholdSat);

    string output= ApplicationTools::getStringParameter("output.counts", mapnh.getParams(), "perType");
    if (output=="perType"){
      //Write count trees:
      string treePathPrefix = ApplicationTools::getStringParameter("output.counts.tree.prefix", mapnh.getParams(), "none");
      if (treePathPrefix != "none") {
        Newick newick;
        for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i) {
          string path = treePathPrefix + TextTools::toString(i + 1) + string(".dnd");
          ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
          Tree* cTree = tree->clone();
          buildCountTree(counts, ids, cTree, i);
          newick.write(*cTree, path);
          delete cTree;
        }
      }
    }
    else if (output=="perSite"){
      string perSitenf= ApplicationTools::getStringParameter("output.counts.file",mapnh.getParams(),"none");
      if (perSitenf!="none"){
        outputTotalCountsPerBranchPerSite(perSitenf, *drtl, ids, model? model: modelSet->getModel(0), *reg);
      }
    }

    
    
    // Rounded counts
    vector< vector<size_t> > countsint;
    for (size_t i=0;i< counts.size(); i++){
      vector<size_t> countsi2;
      for (size_t j=0;j< counts[i].size(); j++)
        countsi2.push_back(static_cast<size_t>(floor( counts[i][j]+0.5)));
      countsint.push_back(countsi2);
    }

    
    //Global homogeneity test:
    bool testGlobal = ApplicationTools::getBooleanParameter("test.global", mapnh.getParams(), true, "", true, false);
    if (testGlobal) {

      vector< vector<size_t> > counts2=countsint;
    
      //Check if some branches are 0:
      for (size_t i = counts2.size(); i > 0; --i) {
        if (VectorTools::sum(counts2[i - 1]) == 0) {
          ApplicationTools::displayResult("Remove branch with no substitution", ids[i - 1]);
          counts2.erase(counts2.begin() + i - 1);
          //ids.erase(ids.begin() + i - 1);
        }
      }
      ApplicationTools::displayResult("Nb. of branches included in test", counts2.size());

    
      ContingencyTableTest test(counts2, 2000);
      ApplicationTools::displayResult("Global Chi2", test.getStatistic());
      ApplicationTools::displayResult("Global Chi2, p-value", test.getPValue());
      double pvalue = SimpleSubstitutionCountsComparison::multinomialTest(counts2);
      ApplicationTools::displayResult("Global heterogeneity test p-value", pvalue);
    }

    //Branch test!
    bool testBranch = ApplicationTools::getBooleanParameter("test.branch", mapnh.getParams(), false, "", true, false);
    if (testBranch) {
      bool testNeighb = ApplicationTools::getBooleanParameter("test.branch.neighbor", mapnh.getParams(), true, "", true, false);
      bool testNegBrL = ApplicationTools::getBooleanParameter("test.branch.negbrlen", mapnh.getParams(), false, "", true, false);
      ApplicationTools::displayBooleanResult("Perform branch clustering", testBranch);
      ApplicationTools::displayBooleanResult("Cluster only neighbor nodes", testNeighb);
      ApplicationTools::displayBooleanResult("Allow len < 0 in clustering", testNegBrL);
      string autoClustDesc = ApplicationTools::getStringParameter("test.branch.auto_cluster", mapnh.getParams(), "Global(threshold=0)");
      string autoClustName;
      map<string, string> autoClustParam;
      KeyvalTools::parseProcedure(autoClustDesc, autoClustName, autoClustParam);
      ApplicationTools::displayResult("Auto-clustering", autoClustName);
      auto_ptr<AutomaticGroupingCondition> autoClust;
      if (autoClustName == "None") {
        autoClust.reset(new NoAutomaticGroupingCondition());
      } else if (autoClustName == "Global") {
        size_t threshold = ApplicationTools::getParameter<size_t>("threshold", autoClustParam, 0);
        ApplicationTools::displayResult("Auto-clutering threshold", threshold);
        CategorySubstitutionRegister* creg = dynamic_cast<CategorySubstitutionRegister*>(reg);
        vector<size_t> ignore;
        if (creg && creg->allowWithin()) {
          size_t n = creg->getNumberOfCategories();
          for (size_t i = 0; i < n; ++i) {
            ignore.push_back(n * (n - 1) + i);
          }
        }
        autoClust.reset(new SumCountsAutomaticGroupingCondition(threshold, ignore));
      } else if (autoClustName == "Marginal") {
        size_t threshold = ApplicationTools::getParameter<size_t>("threshold", autoClustParam, 0);
        ApplicationTools::displayResult("Auto-clutering threshold", threshold);
        autoClust.reset(new AnyCountAutomaticGroupingCondition(threshold));
      } else {
        throw Exception("Unknown automatic clustering option: " + autoClustName);
      }

      //ChiClustering htest(counts, ids, true);
      MultinomialClustering htest(countsint, ids, drtl->getTree(), *autoClust, testNeighb, testNegBrL, true);
      ApplicationTools::displayResult("P-value at root node", *(htest.getPValues().rbegin()));
      ApplicationTools::displayResult("Number of tests performed", htest.getPValues().size());
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
    if (model)
      delete model;
    if (modelSet)
      delete modelSet;
    delete rDist;
    delete reg;
    mapnh.done();

  }
  catch (exception& e)
    {
      cout << e.what() << endl;
      exit(-1);
    }

  return (0);
}
