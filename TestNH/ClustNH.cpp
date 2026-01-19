#include "MultinomialClustering.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

// From bpp:
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/OneProcessSequenceSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/SingleProcessSubstitutionMapping.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h>
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Stat/ContingencyTableTest.h>


using namespace bpp;

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                      Clust NH, version 2                       *" << endl;
  cout << "* Authors: J. Dutheil                                            *" << endl;
  cout << "*          B. Boussau                                            *" << endl;
  cout << "*          L. Guéguen                                            *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication clustnh(args, argv, "ClustNH");
    if (args == 1)
    {
      clustnh.help("clustnh");
      exit(0);
    }

    clustnh.startTimer();
    std::map<std::string, std::string> unparsedParams;

    // Read counts from file:
    string mappingPath = ApplicationTools::getAFilePath("input.counts.file", clustnh.getParams(), true, true, "", true,  "mapping_counts_per_branch_per_type", 1);

    // Optionally, read normalizations:
    string normalizationPath = ApplicationTools::getAFilePath("input.norms.file", clustnh.getParams(), false, true, "", true,  "none", 1);
    bool normalize = normalizationPath != "none";

    ApplicationTools::displayResult("Input counts (branch/type) in file", mappingPath);
    ApplicationTools::displayBooleanResult("Normalize counts", normalize);
    shared_ptr<DataTable> normsTable = nullptr;
    if (normalize)
    {
      ApplicationTools::displayResult("Input normalization (branch/type) in file", normalizationPath);
      std::ifstream inputStreamNorm(normalizationPath, ios::in);
      normsTable = DataTable::read(inputStreamNorm, "\t", true, -1);
    }

    std::ifstream inputStreamCounts(mappingPath, ios::in);
    auto countsTable = DataTable::read(inputStreamCounts, "\t", true, -1);

    // Branches are columns, types are rows, need to transpose and convert to numbers:
    VVdouble counts(countsTable->getNumberOfColumns());
    Vint ids(countsTable->getNumberOfColumns());
    for (size_t i = 0; i < countsTable->getNumberOfColumns(); ++i)
    {
      vector<double> v1(countsTable->getNumberOfRows(), 1.);
      vector<double> v2(countsTable->getNumberOfRows(), 1.);
      double delta = 1.;
      if (normalize)
      {
        for (size_t j = 0; j < countsTable->getNumberOfRows(); ++j)
        {
          double c = TextTools::toDouble((*countsTable)(j, i));
          double m = TextTools::toDouble((*normsTable)(j, i));
          v1[j] = c;
          v2[j] = c / m;
        }
        delta = VectorTools::sum(v1) / VectorTools::sum(v2);
      }
      else
      {
        for (size_t j = 0; j < countsTable->getNumberOfRows(); ++j)
        {
          double c = TextTools::toDouble((*countsTable)(j, i));
          v1[j] = c;
        }
      }
      counts[i].resize(countsTable->getNumberOfRows());
      ids[i] = TextTools::toInt(countsTable->getColumnName(i));
      for (size_t j = 0; j < countsTable->getNumberOfRows(); ++j)
      {
        if (normalize)
        {
          counts[i][j] = delta * v2[j];
        }
        else
        {
          counts[i][j] = v1[j];
        }
      }
    }

    //////////////////////////////////////
    /// HOMOGENEITY TESTS
    //////////////////////////////////////

    bool testGlobal = ApplicationTools::getBooleanParameter("test.global", clustnh.getParams(), false, "", true, false);
    bool testBranch = ApplicationTools::getBooleanParameter("test.branch", clustnh.getParams(), false, "", true, false);

    // Rounded counts
    vector<vector<size_t>> countsint;
    for (size_t i = 0; i < counts.size(); ++i)
    {
      vector<size_t> countsi2;
      for (size_t j = 0; j < counts[i].size(); ++j)
      {
        countsi2.push_back(static_cast<size_t>(floor(counts[i][j] + 0.5)));
      }
      countsint.push_back(countsi2);
    }


    // Global homogeneity test:
    if (testGlobal)
    {
      vector< vector<size_t>> counts2 = countsint;

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
    if (testBranch)
    {
      bool testNeighb = ApplicationTools::getBooleanParameter("test.branch.neighbor", clustnh.getParams(), true, "", true, 1);
      bool testNegBrL = ApplicationTools::getBooleanParameter("test.branch.negbrlen", clustnh.getParams(), false, "", true, 2);
      ApplicationTools::displayBooleanResult("Perform branch clustering", testBranch);
      ApplicationTools::displayBooleanResult("Cluster only neighbor nodes", testNeighb);
      ApplicationTools::displayBooleanResult("Allow len < 0 in clustering", testNegBrL);
      string autoClustDesc = ApplicationTools::getStringParameter("test.branch.auto_cluster", clustnh.getParams(), "Global(threshold=0)", "", true, 1);
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
        ApplicationTools::displayResult("Auto-clustering threshold", threshold);
        // which counts to use:
        size_t numberOfCategories = 0;
        numberOfCategories = ApplicationTools::getParameter<size_t>("number_of_categories", clustnh.getParams(), numberOfCategories, "", false, false);
        bool allowWithin = ApplicationTools::getBooleanParameter("allow_within_category_substitutions", clustnh.getParams(), false, "", true, false);
        vector<size_t> toIgnore;
        if (numberOfCategories > 0 && allowWithin)
        {
          size_t n = numberOfCategories;
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
        ApplicationTools::displayResult("Auto-clustering threshold", threshold);
        autoClust.reset(new AnyCountAutomaticGroupingCondition(threshold));
      }
      else
      {
        throw Exception("Unknown automatic clustering option: " + autoClustName);
      }

      // Read the tree:
      shared_ptr<Tree> tree = PhylogeneticsApplicationTools::getTree(clustnh.getParams());
      // Convert to NHX if input tree is newick or nexus?
      string treeIdOut = ApplicationTools::getAFilePath("output.tree_with_id.file", clustnh.getParams(), false, false, "none", 1);
      if (treeIdOut != "none")
      {
        Nhx nhx(true);
        nhx.writeTree(*tree, treeIdOut, true);
      }

      // ChiClustering htest(counts, ids, true);
      MultinomialClustering htest(countsint, ids, *tree, *autoClust, testNeighb, testNegBrL, true);
      ApplicationTools::displayResult("P-value at root node", *(htest.getPValues().rbegin()));
      ApplicationTools::displayResult("Number of tests performed", htest.getPValues().size());
      TreeTemplate<Node>* htree = htest.getTree();
      Newick newick;
      string clusterTreeOut = ApplicationTools::getAFilePath("output.cluster_tree.file", clustnh.getParams(), false, false, "", true, "clusters.dnd", 1);
      ApplicationTools::displayResult("Output cluster tree to", clusterTreeOut);
      newick.writeTree(*htree, clusterTreeOut, true);
      delete htree;
    }

    /////////////////////////////////
    // clean up

    clustnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return 0;
}
