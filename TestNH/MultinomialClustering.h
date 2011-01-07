
#ifndef _MULTINOMIALCLUSTERING_H_
#define _MULTINOMIALCLUSTERING_H_

#include <Bpp/Phyl/Distance/AbstractAgglomerativeDistanceMethod.h>
#include <Bpp/Phyl/Distance/HierarchicalClustering.h>

#include <vector>

using namespace bpp;
using namespace std;

class MultinomialClustering :
  public AbstractAgglomerativeDistanceMethod
{
  private:
    vector< vector<unsigned int> > counts_;
    bool neighborsOnly_;
    bool negativeBrlen_;
	
  public:
		MultinomialClustering(
        const vector< vector<unsigned int> >& counts,
        const vector<int>& ids,
        const Tree& tree,
        bool neighborsOnly = false,
        bool negativeBrlen = false,
        bool verbose = false);
		
    virtual ~MultinomialClustering() {}

	public:
		TreeTemplate<Node>* getTree() const;

	protected:
		/**
		 * @brief Returns the pair with minimum distance and actualizes the vectors.
		 *
		 * The vector at position bestPair[0] is now the sum of vectors bestPair[0] and bestPair[1].
		 * It is then used for computation of distances.
		 */
		vector<unsigned int> getBestPair() throw (Exception);
		vector<double> computeBranchLengthsForPair(const vector<unsigned int>& pair);
		double computeDistancesFromPair(const vector<unsigned int>& pair, const vector<double>& branchLengths, unsigned int pos);
		void finalStep(int idRoot);	
		virtual Node* getLeafNode(int id, const string& name);
		virtual Node* getParentNode(int id, Node* son1, Node* son2);

    static double multinomLogL(const vector<unsigned int>& v, const vector<double>& p);
    
    static double logFact(unsigned int x) {
      double f = 0;
      for (unsigned int i = 1; i <= x; ++i)
        f += log(static_cast<double>(i));
      return f;
    }
    static double getDist(const vector<unsigned int>& v1, const vector<unsigned int>&v2);
};

#endif //_MULTINOMIALCLUSTERING_H_

