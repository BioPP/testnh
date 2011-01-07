
#include <Bpp/Phyl/Distance/AbstractAgglomerativeDistanceMethod.h>

#include <vector>

using namespace bpp;
using namespace std;

struct ChiClusterInfos
{
  public:
  	unsigned int numberOfLeaves;
	  double length;
    vector<unsigned int> index;

  public:
    ChiClusterInfos():  numberOfLeaves(0), length(0), index() {}
};

class ChiClustering:
  public AbstractAgglomerativeDistanceMethod
{
  private:
    vector< vector<unsigned int> > counts_;
	public:
		ChiClustering(const vector< vector<unsigned int> >& counts, const vector<int>& ids, bool verbose = false);

		virtual ~ChiClustering() {}

	public:
		TreeTemplate<Node>* getTree() const;

	protected:
    std::vector<unsigned int> getBestPair() throw (Exception);
    std::vector<double> computeBranchLengthsForPair(const std::vector<unsigned int>& pair);
		double computeDistancesFromPair(const std::vector<unsigned int>& pair, const std::vector<double> & branchLengths, unsigned int pos);
		void finalStep(int idRoot);	
		virtual Node* getLeafNode(int id, const std::string& name);
		virtual Node* getParentNode(int id, Node* son1, Node* son2);
};

