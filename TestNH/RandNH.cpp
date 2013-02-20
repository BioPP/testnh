//
// File: RandNH.cpp
// Created by: Julien Dutheil
// Created on: Jan Tue 25 16:15 2011
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
#include <Bpp/Phyl/Likelihood.all>

// From bpp-seq:
#include <Bpp/Seq/Alphabet.all>
// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "randnh parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the package manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

struct CandidatePart {
  int id;
  size_t depth;

  CandidatePart(int i, size_t d):
    id(i), depth(d) {}
};
   
struct CandidatePartComp {
  bool operator()(const CandidatePart& cp1, const CandidatePart& cp2) const {
    return (cp1.depth > cp2.depth);
  }
};

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     RandNH, version 0.1.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  25/01/11 *" << endl;
  cout << "*          B. Boussau                       Last Modif. 28/04/11 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  BppApplication randnh(args, argv, "RandNH");
  randnh.startTimer();

  //Take an input tree:
  auto_ptr<Tree> tree(PhylogeneticsApplicationTools::getTree(randnh.getParams()));
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

  //Output to file in NHX format:
  string treeIdOut = ApplicationTools::getAFilePath("output.tree.file", randnh.getParams(), false, false);
  Nhx nhx(true);
  nhx.write(*tree, treeIdOut);

  //Now get the number of partitions and assign nodes:
  size_t nbParts = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", randnh.getParams(), 2);
  ApplicationTools::displayResult("Number of partitions", nbParts);
  if (nbParts < 2)
    throw Exception("There should be at least 2 partitions.");

  vector<int> ids = tree->getNodesId();
  ids.pop_back();

  //We create a id vector weighted for node depth:
  vector<double> idWeights;
  map<int, size_t> nodesDepths;
  TreeTools::getDepths(*tree, tree->getRootId(), nodesDepths);
  for (size_t i = 0; i < ids.size(); ++i) {
    idWeights.push_back(static_cast<double>(nodesDepths[ids[i]] + 1));
  }

  map<int, size_t> partitionsIndex;
  
  string modelType = ApplicationTools::getStringParameter("nonhomogeneous.type_of_model", randnh.getParams(), "join");

  if (modelType == "join")
  {
    vector<int> partRootIds(nbParts - 1);
    RandomTools::getSample(ids, idWeights, partRootIds);

    //we have to make sure we get the larger partition first:
    multiset<CandidatePart, CandidatePartComp> depths;
    for (size_t i = 0; i < nbParts - 1; ++i) {
      depths.insert(CandidatePart(partRootIds[i], TreeTools::getDepth(*tree, partRootIds[i])));
    }
    for (size_t i = 0; i < ids.size(); ++i)
      partitionsIndex[ids[i]] = 0;
    size_t p = 1;
    for (set<CandidatePart>::iterator it = depths.begin(); it != depths.end(); ++it) {
      cout << it->id << "\t" << it->depth << endl;
      vector<int> subids = TreeTools::getNodesId(*tree, it->id);
      for (size_t i = 0; i < subids.size(); ++i)
        partitionsIndex[subids[i]] = p;
      p++;
    }
  }
  else if (modelType == "free")
  {
    //For now we assume that all partitions are equally likely.
    for (size_t i = 0; i < ids.size(); ++i)
      partitionsIndex[ids[i]] = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(static_cast<int>(nbParts));
  }
  else throw Exception("Model type should be one of 'free' or 'join'");

  map<size_t, vector<int> > partitions;
  for (map<int, size_t>::iterator it = partitionsIndex.begin(); it != partitionsIndex.end(); ++it) {
    partitions[it->second].push_back(it->first);
  }
    
  //Print partitions:
  string modelOutPath = ApplicationTools::getAFilePath("output.model.file", randnh.getParams(), true, false);
  ofstream modelOut(modelOutPath.c_str(), ios::out);
  size_t index = 1;
  for (map<size_t, vector<int> >::iterator it = partitions.begin(); it != partitions.end(); ++it) {
    modelOut << "model" << index << ".nodes_id = " << it->second[0];
    for (size_t i = 1; i < it->second.size(); ++i)
      modelOut << "," << it->second[i];
    modelOut << endl;
    index++;
  }
  modelOut.close();

  //Clean and exit:
  randnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}
