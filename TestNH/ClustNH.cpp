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
#include <Bpp/Phyl/Mapping/PhyloMappings/OneProcessSequenceSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/SingleProcessSubstitutionMapping.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h>
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>

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

    //Read counts from file:
    string perSitenf = ApplicationTools::getStringParameter("input.mapping.file", outputArgs, "mapping_counts_per_site_per_branch", "", true, 1);

    ApplicationTools::displayResult(string("Input counts (branch/site) to file"), perSitenf);


    auto countsTable = DataTable::read(ifstream(perSitenf, ios::in), sep = '\t', header = true, row.names = 1);

    VVdouble counts = 

    //////////////////////////////////////
    /// HOMOGENEITY TESTS
    /////////////////////////////////////


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
