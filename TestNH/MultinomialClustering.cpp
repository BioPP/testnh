
#include "MultinomialClustering.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/NodeTemplate.h>

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;
	
MultinomialClustering::MultinomialClustering(
    const vector< vector<unsigned int> >& counts,
    const vector<int>& ids,
    const Tree& tree,
    const AutomaticGroupingCondition& autoGroup,
    bool neighborsOnly, bool negativeBrlen, bool verbose):
  AbstractAgglomerativeDistanceMethod(verbose), counts_(counts), neighborsOnly_(neighborsOnly), negativeBrlen_(negativeBrlen)
{
  unsigned int n = counts.size();
  matrix_.resize(n);
  MatrixTools::fill(matrix_, (neighborsOnly_ ? 2. : 1.));
  vector<string> names(n);
  for (unsigned int i = 0; i < n; ++i)
    names[i] = TextTools::toString(ids[i]);
  for (unsigned int i = 0; i < n; ++i) {
    if (verbose)
      ApplicationTools::displayGauge(i, n - 1, '=');
    matrix_(i, i) = 0.;
    //Check for small counts:
    if (autoGroup.check(counts_[i])) {
      for (size_t j = 0; j < i; ++j) {
        if (neighborsOnly_) {
          if (tree.getFatherId(ids[i]) == ids[j])
            matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
          if (tree.getFatherId(ids[j]) == ids[i])
            matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
          if (tree.getFatherId(ids[i]) == tree.getRootId()
           && tree.getFatherId(ids[j]) == tree.getRootId())
            matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
        } else {
          matrix_(i, j) = matrix_(j, i) = getDist(counts_[i], counts_[j]);
        }
      }
    } else {
      //node i has distance 0 with his father, 1. with brothers and sons and 2. with all others
      for (size_t j = 0; j < n; ++j) {
        if (j != i) {
          int fatherId_i = tree.getFatherId(ids[i]);
          int fatherId_j = tree.getFatherId(ids[j]);
          if (fatherId_i == ids[j]) {
            matrix_(i, j) = matrix_(j, i) = 0.;
            ApplicationTools::displayResult("Automatically cluster node", names[i] + " with " + names[j]);
          } else if (fatherId_i == fatherId_j || ids[i] == fatherId_j)
            matrix_(i, j) = matrix_(j, i) = 1.;
          //else remains 2.
        }
      }
    }
  }
  if (verbose)
    ApplicationTools::displayTaskDone();
  matrix_.setNames(names);
  computeTree(true);
  if (verbose)
    ApplicationTools::displayTaskDone();
}

TreeTemplate<Node>* MultinomialClustering::getTree() const
{
	Node* root = TreeTemplateTools::cloneSubtree<Node>(* dynamic_cast<TreeTemplate<NodeTemplate<ClusterInfos> > *>(tree_)->getRootNode());
	return new TreeTemplate<Node>(root);
}

vector<unsigned int> MultinomialClustering::getBestPair() throw (Exception)
{
	vector<unsigned int> bestPair(2);
	double distMin = -std::log(0.);
	for (map<unsigned int, Node *>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
  {
		unsigned int id = i->first;
		map<unsigned int, Node *>::iterator j = i;
		j++;
		for(; j != currentNodes_.end(); j++)
    {
			unsigned int jd = j->first;
			double dist = matrix_(id, jd);
			if(dist < distMin)
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
vector<double> MultinomialClustering::computeBranchLengthsForPair(const vector<unsigned int>& pair)
{
	vector<double> d(2);
	double dist = matrix_(pair[0], pair[1]) / 2.;
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

double MultinomialClustering::computeDistancesFromPair(const vector<unsigned int>& pair, const vector<double>& branchLengths, unsigned int pos)
{
  //Perform a multinomial LRT:
  vector<unsigned int>* v1 = &counts_[pos];
  vector<unsigned int>* v2 = &counts_[pair[0]];
  //Check if nodes are neighbors:
  if (matrix_(pair[0], pos) > 1. && matrix_(pair[1], pos) > 1.) return 2.;
  return getDist(*v1, *v2);
}

void MultinomialClustering::finalStep(int idRoot)
{
	NodeTemplate<ClusterInfos>* root = new NodeTemplate<ClusterInfos>(idRoot);
	map<unsigned int, Node*>::iterator it = currentNodes_.begin();
	unsigned int i1 = it->first;
	Node* n1        = it->second;
	it++;
	unsigned int i2 = it->first;
	Node* n2        = it->second;
	double d = matrix_(i1, i2) / 2;
	root->addSon(n1);
	root->addSon(n2);
	n1->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos> *>(n1)->getInfos().length); 
	n2->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos> *>(n2)->getInfos().length); 
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

double MultinomialClustering::getDist(const vector<unsigned int>& v1, const vector<unsigned int>&v2)
{
  unsigned int s1 = VectorTools::sum(v1);
  unsigned int s2 = VectorTools::sum(v2);
  if (s1 == 0 || s2 == 0) return 0.;
  vector<double> p1(v1.size());
  vector<double> p2(v1.size());
  vector<double> p12(v1.size());
  for (size_t i = 0; i < v1.size(); ++i) {
    p1[i]  = static_cast<double>(v1[i]) / static_cast<double>(s1);
    p2[i]  = static_cast<double>(v2[i]) / static_cast<double>(s2);
    p12[i] = static_cast<double>(v1[i] + v2[i]) / static_cast<double>(s1 + s2);
  }
  double logL1_1 = multinomLogL(v1, p12);
  double logL2_1 = multinomLogL(v2, p12);
  double logL1_2 = multinomLogL(v1, p1);
  double logL2_2 = multinomLogL(v2, p2);
  double stat = -2 * (logL1_1 + logL2_1 - logL1_2 - logL2_2);
	double d = RandomTools::pChisq(stat, v1.size() - 1);
  return d;
}

double MultinomialClustering::multinomLogL(const vector<unsigned int>& v, const vector<double>& p)
{
  unsigned int n = VectorTools::sum(v);
  double l = logFact(n);
  for (size_t i = 0; i < v.size(); ++i) {
    if (p[i] == 0) 
      return 0;
      //throw Exception("Error, expectation has probability 0!");
    l += v[i] * log(static_cast<double>(p[i]));
    l -= logFact(v[i]);
  }
  return l;
}

