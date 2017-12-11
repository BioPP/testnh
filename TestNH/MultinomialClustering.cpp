#include "MultinomialClustering.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/PseudoNewtonOptimizer.h>

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

double SimpleSubstitutionCountsComparison::multinomLogL(const vector<size_t>& v, const vector<double>& p)
{
  size_t n = VectorTools::sum(v);
  double l = logFact(n);
  for (size_t i = 0; i < v.size(); ++i) {
    if (p[i] == 0) {
      if (v[i] == 0) {
        //Do nothing, everything ok
      } else {
        //This is impossible
        return log(0.);
      }
    } else {
      l += static_cast<double>(v[i]) * log(p[i]);
      l -= logFact(v[i]);
    }
  }
  return l;
}

double SimpleSubstitutionCountsComparison::multinomialTest(const vector< vector<size_t> >& counts)
{
  size_t nbCounts = counts.size();
  size_t countsSize = counts[0].size();
  vector<size_t> sumCounts(countsSize);
  //Compute likelihood of independent model:
  double logL1 = 0;
  for (size_t i = 0; i < nbCounts; ++i) {
    size_t s = VectorTools::sum(counts[i]);
    if (s == 0)
      throw Exception("SimpleSubstitutionCountsComparison::multinomialTest. Counts with zero observation can't be accounted for.");
    sumCounts += counts[i];
    vector<double> p(countsSize);
    for (size_t j = 0; j < countsSize; ++j) {
      p[j] = static_cast<double>(counts[i][j]) / static_cast<double>(s);
    }
    logL1 += multinomLogL(counts[i], p);
  }
  //Compute likelihood of shared model:
  double logL2 = 0;
  vector<double> p(countsSize);
  size_t s = VectorTools::sum(sumCounts);
  for (size_t j = 0; j < countsSize; ++j) {
    p[j] = static_cast<double>(sumCounts[j]) / static_cast<double>(s);
  }
  for (size_t i = 0; i < nbCounts; ++i) {
    logL2 += multinomLogL(counts[i], p);
  }
  //Now compare models:
  double statistic = -2 * (logL2 - logL1);
  //cout << logL2 << "\t" << logL1 << "\t" << statistic << "\t" << (countsSize - 1) * (nbCounts - 1) << endl;
  double pvalue = 1. - RandomTools::pChisq(statistic, static_cast<double>(countsSize - 1) * static_cast<double>(nbCounts - 1));
  return pvalue;
}

void SimpleSubstitutionCountsComparison::computePValue() {
  size_t s1 = VectorTools::sum(counts1_);
  size_t s2 = VectorTools::sum(counts2_);
  if (s1 == 0 || s2 == 0) {
    statistic_ = 0.;
    pvalue_ = 1.;
    return;
  }
  vector<double> p1(counts1_.size());
  vector<double> p2(counts1_.size());
  vector<double> p12(counts1_.size());
  for (size_t i = 0; i < counts1_.size(); ++i) {
    p1[i]  = static_cast<double>(counts1_[i]) / static_cast<double>(s1);
    p2[i]  = static_cast<double>(counts2_[i]) / static_cast<double>(s2);
    p12[i] = static_cast<double>(counts1_[i] + counts2_[i]) / static_cast<double>(s1 + s2);
  }
  double logL1_1 = multinomLogL(counts1_, p12);
  double logL2_1 = multinomLogL(counts2_, p12);
  double logL1_2 = multinomLogL(counts1_, p1);
  double logL2_2 = multinomLogL(counts2_, p2);
  statistic_ = -2 * (logL1_1 + logL2_1 - logL1_2 - logL2_2);
  pvalue_ = 1. - RandomTools::pChisq(statistic_, static_cast<double>(counts1_.size() - 1));
}

MultinomialClustering::MultinomialClustering(
  const vector< vector<size_t> >& counts,
  const vector<uint>& ids,
  const PhyloTree& tree,
  const AutomaticGroupingCondition& autoGroup,
  bool neighborsOnly, bool negativeBrlen, bool verbose):
  AbstractAgglomerativeDistanceMethod(verbose, true), counts_(counts), neighborsOnly_(neighborsOnly), negativeBrlen_(negativeBrlen), test_(new SimpleSubstitutionCountsComparison()), pvalues_()
{
  size_t n = counts.size();
  matrix_.resize(n);
  MatrixTools::fill(matrix_.asMatrix(), (neighborsOnly_ ? 2. : 1.));
  vector<string> names(n);
  for (size_t i = 0; i < n; ++i)
    names[i] = TextTools::toString(ids[i]);
  for (size_t i = 0; i < n; ++i) {
    if (verbose)
      ApplicationTools::displayGauge(i, n - 1, '=');
    matrix_(i, i) = 0.;
    //Check for small counts:
    if (autoGroup.check(counts_[i])) {
      for (size_t j = 0; j < i; ++j) {
        if (matrix_(i, j) > 0) {
          if (neighborsOnly_) {
            if (tree.getFather(ids[i]) == ids[j])
              matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
            if (tree.getFather(ids[j]) == ids[i])
              matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
            if (tree.getFather(ids[i]) == tree.getRootIndex()
                && tree.getFather(ids[j]) == tree.getRootIndex())
              matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
          } else {
            matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
          }
        }//else distance was already assigned because of autoclustering.
      }
    } else {
      //node i has distance -1 with his father, 1. with brothers and sons and 2. with all others
      for (size_t j = 0; j < n; ++j) {
        if (j != i) {
          uint fatherId_i = tree.getFather(ids[i]);
          uint fatherId_j = tree.getFather(ids[j]);
          if (fatherId_i == ids[j]) {
            matrix_(i, j) = matrix_(j, i) = -1.;
            ApplicationTools::displayResult("Automatically cluster node", names[i] + " with " + names[j]);
          } else if (fatherId_i == fatherId_j || ids[i] == fatherId_j) {
            if (matrix_(i, j) > 0) {
              matrix_(i, j) = matrix_(j, i) = 1.;
            } //otherwise son node was also automatically assigned!
          }
          //else remains 2.
        }
      }
    }
  }
  if (verbose)
    ApplicationTools::displayTaskDone();
  matrix_.setNames(names);
  computeTree();
  if (verbose)
    ApplicationTools::displayTaskDone();
}

ClusterPhyloTree* MultinomialClustering::getTree() const
{
  Node* root = TreeTemplateTools::cloneSubtree<Node>(* dynamic_cast<TreeTemplate<NodeTemplate<ClusterInfos> > *>(tree_)->getRootNode());
  return new TreeTemplate<Node>(root);
}

vector<size_t> MultinomialClustering::getBestPair() throw (Exception)
{
  vector<size_t> bestPair(2);
  double distMin = -std::log(0.);
  for (map<size_t, Node *>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
  {
    size_t id = i->first;
    map<size_t, Node *>::iterator j = i;
    j++;
    for (; j != currentNodes_.end(); j++)
    {
      size_t jd = j->first;
      double dist = matrix_(id, jd);
      if (dist < distMin)
      {
        distMin = dist;
        bestPair[0] = id;
        bestPair[1] = jd;
      }
    }
  }
  // actualize vectors:
  counts_[bestPair[0]] += counts_[bestPair[1]];
  return bestPair;  
}
vector<double> MultinomialClustering::computeBranchLengthsForPair(const vector<size_t>& pair)
{
  vector<double> d(2);
  double dist = max(0., matrix_(pair[0], pair[1])) / 2.; //In case of automatic clustering, -1 corresponds to a distance of 0.
  double d1 = dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[0]])->getInfos().length; 
  double d2 = dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[1]])->getInfos().length; 
  if (negativeBrlen_) {
    d[0] = dist - d1;
    d[1] = dist - d2;
  } else {
    double h = max(dist, max(d1, d2));
    d[0] = h - d1;
    d[1] = h - d2;
  }
  //d[0] = dist; 
  //d[1] = dist; 
  return d;
}

double MultinomialClustering::computeDistancesFromPair(const vector<size_t>& pair, const vector<double>& branchLengths, size_t pos)
{
  //Check if nodes are neighbors:
  if (matrix_(pair[0], pos) > 1. && matrix_(pair[1], pos) > 1.) return 2.;
  //Check if nodes are automatically clustered;
  if (matrix_(pair[0], pos) < 0. || matrix_(pair[1], pos) < 0.) {
    return -1.;
  }
  //Perform a multinomial LRT:
  vector<size_t>* v1 = &counts_[pos];
  vector<size_t>* v2 = &counts_[pair[0]];
  return getDist(*v1, *v2);
}

void MultinomialClustering::finalStep(int idRoot)
{
  NodeTemplate<ClusterInfos>* root = new NodeTemplate<ClusterInfos>(idRoot);
  map<size_t, Node*>::iterator it = currentNodes_.begin();
  size_t i1 = it->first;
  Node* n1        = it->second;
  it++;
  size_t i2 = it->first;
  Node* n2        = it->second;
  double d = matrix_(i1, i2) / 2;
  root->addSon(n1);
  root->addSon(n2);
  double d1 = dynamic_cast<NodeTemplate<ClusterInfos>*>(n1)->getInfos().length; 
  double d2 = dynamic_cast<NodeTemplate<ClusterInfos>*>(n2)->getInfos().length; 
  if (negativeBrlen_) {
    n1->setDistanceToFather(d - d1); 
    n2->setDistanceToFather(d - d2);
  } else {
    double h = max(d, max(d1, d2));
    n1->setDistanceToFather(h - d1); 
    n2->setDistanceToFather(h - d2);
  }
  //n1->setDistanceToFather(d); 
  //n2->setDistanceToFather(d); 
  tree_ = new TreeTemplate< NodeTemplate<ClusterInfos> >(root);
}

Node* MultinomialClustering::getLeafNode(int id, const string& name)
{
  ClusterInfos infos;
  infos.numberOfLeaves = 1;
  infos.length = 0.;
  NodeTemplate<ClusterInfos>* leaf = new NodeTemplate<ClusterInfos>(id, name);
  leaf->setInfos(infos);
  return leaf;
}

Node* MultinomialClustering::getParentNode(int id, Node* son1, Node* son2)
{
  ClusterInfos infos;
  infos.numberOfLeaves = 
    dynamic_cast<NodeTemplate<ClusterInfos> *>(son1)->getInfos().numberOfLeaves
    + dynamic_cast<NodeTemplate<ClusterInfos> *>(son2)->getInfos().numberOfLeaves;
  infos.length = dynamic_cast<NodeTemplate<ClusterInfos> *>(son1)->getInfos().length + son1->getDistanceToFather();
  Node* parent = new NodeTemplate<ClusterInfos>(id);
  dynamic_cast<NodeTemplate<ClusterInfos> *>(parent)->setInfos(infos);
  parent->addSon(son1);
  parent->addSon(son2);
  return parent;
}

double MultinomialClustering::getDist(const vector<size_t>& v1, const vector<size_t>&v2)
{
  test_->setCounts(v1, v2);
  pvalues_.push_back(test_->getPValue());
  return 1. - test_->getPValue();
}

