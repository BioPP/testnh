#ifndef _MULTINOMIALCLUSTERING_H_
#define _MULTINOMIALCLUSTERING_H_

#include <Bpp/Phyl/Distance/AbstractAgglomerativeDistanceMethod.h>
#include <Bpp/Phyl/Distance/HierarchicalClustering.h>
#include <Bpp/Numeric/Stat/StatTest.h>

#include <vector>

using namespace bpp;
using namespace std;

template<class T>
class bpp_ptr
{
  private:
    T* ptr_;

  public:
    bpp_ptr(T* ptr = 0): ptr_(ptr) {}
    
    ~bpp_ptr() { delete ptr_; }

    bpp_ptr(const bpp_ptr& ptr):
      ptr_(ptr->clone()) {}
    
    bpp_ptr& operator=(const bpp_ptr& ptr)
    {
      delete ptr_;
      ptr_ = dynamic_cast<T*>(ptr->clone());
    }

    void reset(T* ptr = 0) {
      if (ptr_) delete ptr_;
      ptr_ = ptr;
    }

    T& operator*() { return *ptr_; } 
    const T& operator*() const { return *ptr_; } 

    T* operator->() { return ptr_; } 
    const T* operator->() const { return ptr_; } 

};

class SubstitutionCountsComparison :
  public StatTest
{
  protected:
    vector<unsigned int> counts1_, counts2_;
    double statistic_, pvalue_;

  public:
    SubstitutionCountsComparison(): counts1_(), counts2_(), statistic_(0), pvalue_(1.) {}

  public:
    void setCounts(const vector<unsigned int>& c1, const vector<unsigned int>& c2) {
      counts1_ = c1;
      counts2_ = c2;
      computePValue();
    }

    double getStatistic() const { return statistic_; }
    double getPValue() const { return pvalue_; }

  protected:
    virtual void computePValue() = 0;
};


class SimpleSubstitutionCountsComparison:
  public SubstitutionCountsComparison
{
  public:
    SimpleSubstitutionCountsComparison* clone() const { return new SimpleSubstitutionCountsComparison(*this); }

  public:
    string getName() const { return "Independent multinomial test."; }

  protected:
    void computePValue();
    static double logFact(unsigned int x) {
      double f = 0;
      for (unsigned int i = 1; i <= x; ++i)
        f += log(static_cast<double>(i));
      return f;
    }

    static double multinomLogL(const vector<unsigned int>& v, const vector<double>& p);
    
};


class AutomaticGroupingCondition
{
  public:
    virtual ~AutomaticGroupingCondition() {}

  public:
    virtual bool check(const vector<unsigned int>& counts) const = 0;
};

class NoAutomaticGroupingCondition:
  public AutomaticGroupingCondition
{
  public:
    ~NoAutomaticGroupingCondition() {}

  public:
    bool check(const vector<unsigned int>& counts) const {
      return true;
    }
};

class SumCountsAutomaticGroupingCondition:
  public AutomaticGroupingCondition
{
  private:
    unsigned int threshold_;
    vector<size_t> ignore_;

  public:
    SumCountsAutomaticGroupingCondition(unsigned int threshold = 0, const vector<size_t>& ignore = vector<size_t>(0)):
      threshold_(threshold),
      ignore_(ignore)
    {}

    ~SumCountsAutomaticGroupingCondition() {}

  public:
    bool check(const vector<unsigned int>& counts) const {
      unsigned int s = 0;
      for (size_t i = 0; i < counts.size(); ++i)
        if (!VectorTools::contains(ignore_, i))
          s += counts[i];
      return s > threshold_;
    }
};

class AnyCountAutomaticGroupingCondition:
  public AutomaticGroupingCondition
{
  private:
    unsigned int threshold_;

  public:
    AnyCountAutomaticGroupingCondition(unsigned int threshold = 0):
      threshold_(threshold)
    {}

    ~AnyCountAutomaticGroupingCondition() {}

  public:
    bool check(const vector<unsigned int>& counts) const {
      for (size_t i = 0; i < counts.size(); ++i) {
        if (counts[i] <= threshold_) return false;
      }
      return true;
    }
};

class MultinomialClustering :
  public AbstractAgglomerativeDistanceMethod
{
  private:
    vector< vector<unsigned int> > counts_;
    bool neighborsOnly_;
    bool negativeBrlen_;
    bpp_ptr<SubstitutionCountsComparison> test_;
	
  public:
		MultinomialClustering(
        const vector< vector<unsigned int> >& counts,
        const vector<int>& ids,
        const Tree& tree,
        const AutomaticGroupingCondition& autoGroup,
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
    double getDist(const vector<unsigned int>& v1, const vector<unsigned int>&v2);
};

#endif //_MULTINOMIALCLUSTERING_H_

