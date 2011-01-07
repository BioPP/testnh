
#include "ChiClustering.h"

#include <Bpp/Numeric/Stat/ContingencyTableTest.h>
#include <Bpp/Phyl/NodeTemplate.h>

using namespace bpp;
using namespace std;

ChiClustering::ChiClustering(const vector< vector<unsigned int> >& counts, const vector<int>& ids, bool verbose):
  AbstractAgglomerativeDistanceMethod(verbose), counts_(counts)
{
  unsigned int n = counts.size();
  matrix_.resize(n);
  vector<string> names(n);
  for (unsigned int i = 0; i < n; ++i) {
    ApplicationTools::displayGauge(i, n - 1, '=');
    matrix_(i, i) = 0.;
    names[i] = TextTools::toString(ids[i]);
    vector< vector<unsigned int> > table;
    table.push_back(counts[i]);
    for (size_t j = 0; j < i; ++j) {
      table.push_back(counts[j]);
      double stat = 1.;
      try {
        ContingencyTableTest testbr(table, 0, false);
        //stat = testbr.getStatistic();
        stat = 1. - testbr.getPValue();
      } catch (Exception& ex) {}
      matrix_(i, j) = matrix_(j, i) = stat;
    }
  }
  cout << matrix_.size() << endl;
  matrix_.setNames(names);
  computeTree(true);
}

TreeTemplate<Node>* ChiClustering::getTree() const
{
	Node* root = TreeTemplateTools::cloneSubtree<Node>(*dynamic_cast<TreeTemplate<NodeTemplate<ChiClusterInfos> > *>(tree_)->getRootNode());
	return new TreeTemplate<Node>(root);
}

vector<unsigned int> ChiClustering::getBestPair() throw (Exception)
{
	vector<unsigned int> bestPair(2);
	double distMin = -std::log(0.);
	for (map<unsigned int, Node *>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
  {
		unsigned int id = i->first;
		map<unsigned int, Node *>::iterator j = i;
		j++;
		for (; j != currentNodes_.end(); j++)
    {
			unsigned int jd = j->first;
			double dist = matrix_(id, jd);
			if (dist < distMin)
      {
				distMin = dist;
				bestPair[0] = id;
				bestPair[1] = jd;
			}
		}
	}
	return bestPair;	
}
vector<double> ChiClustering::computeBranchLengthsForPair(const vector<unsigned int>& pair)
{
	vector<double> d(2);
	double dist = matrix_(pair[0], pair[1]) / 2.;
	d[0] = dist - dynamic_cast<NodeTemplate<ChiClusterInfos> *>(currentNodes_[pair[0]])->getInfos().length; 
	d[1] = dist - dynamic_cast<NodeTemplate<ChiClusterInfos> *>(currentNodes_[pair[1]])->getInfos().length; 
	return d;
}

double ChiClustering::computeDistancesFromPair(const vector<unsigned int>& pair, const vector<double>& branchLengths, unsigned int pos)
{
  vector< vector<unsigned int> > table;
  vector<unsigned int> index = dynamic_cast<NodeTemplate<ChiClusterInfos>*>(currentNodes_[pos])->getInfos().index;
  VectorTools::append(index, dynamic_cast<NodeTemplate<ChiClusterInfos>*>(currentNodes_[pair[0]])->getInfos().index);
  VectorTools::append(index, dynamic_cast<NodeTemplate<ChiClusterInfos>*>(currentNodes_[pair[1]])->getInfos().index);
  for (size_t i = 0; i < index.size(); ++i) {
    table.push_back(counts_[index[i]]);
  }
  double stat = 1.;
  try {
    ContingencyTableTest test(table, 0, false);
	  //stat = test.getStatistic();
	  stat = 1. - test.getPValue();
  } catch(Exception& ex) {}
  return stat;
}

void ChiClustering::finalStep(int idRoot)
{
	NodeTemplate<ChiClusterInfos>* root = new NodeTemplate<ChiClusterInfos>(idRoot);
	map<unsigned int, Node*>::iterator it = currentNodes_.begin();
	unsigned int i1 = it->first;
	Node* n1        = it->second;
	it++;
	unsigned int i2 = it->first;
	Node* n2        = it->second;
	double d = matrix_(i1, i2) / 2;
	root->addSon(n1);
	root->addSon(n2);
	n1->setDistanceToFather(d - dynamic_cast<NodeTemplate<ChiClusterInfos>*>(n1)->getInfos().length); 
	n2->setDistanceToFather(d - dynamic_cast<NodeTemplate<ChiClusterInfos>*>(n2)->getInfos().length); 
	//n1->setDistanceToFather(d); 
	//n2->setDistanceToFather(d); 
	tree_ = new TreeTemplate< NodeTemplate<ChiClusterInfos> >(root);
}

Node* ChiClustering::getLeafNode(int id, const string& name)
{
	ChiClusterInfos infos;
	infos.numberOfLeaves = 1;
	infos.length = 0.;
  infos.index.push_back(static_cast<int>(id)); //Leaves are numbered according to their order in the matrix.
	NodeTemplate<ChiClusterInfos>* leaf = new NodeTemplate<ChiClusterInfos>(id, name);
	leaf->setInfos(infos);
	return leaf;
}

Node* ChiClustering::getParentNode(int id, Node* son1, Node* son2)
{
	ChiClusterInfos infos;
	infos.numberOfLeaves = 
		dynamic_cast<NodeTemplate<ChiClusterInfos>*>(son1)->getInfos().numberOfLeaves
	+ dynamic_cast<NodeTemplate<ChiClusterInfos>*>(son2)->getInfos().numberOfLeaves;
	infos.length = dynamic_cast<NodeTemplate<ChiClusterInfos> *>(son1)->getInfos().length + son1->getDistanceToFather();
  VectorTools::append(infos.index, dynamic_cast<NodeTemplate<ChiClusterInfos>*>(son1)->getInfos().index);
  VectorTools::append(infos.index, dynamic_cast<NodeTemplate<ChiClusterInfos>*>(son2)->getInfos().index);
  Node* parent = new NodeTemplate<ChiClusterInfos>(id);
	dynamic_cast<NodeTemplate<ChiClusterInfos>*>(parent)->setInfos(infos);
	parent->addSon(son1);
	parent->addSon(son2);
	return parent;
}


