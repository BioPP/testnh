#ifndef _MULTINOMIALCLUSTERING_H_
#define _MULTINOMIALCLUSTERING_H_

#include <Bpp/Phyl/Distance/AbstractAgglomerativeDistanceMethod.h>
#include <Bpp/Phyl/Distance/HierarchicalClustering.h>
#include <Bpp/Numeric/Stat/StatTest.h>

#include <vector>

using namespace bpp;
using namespace std;

class SubstitutionCountsComparison :
  public StatTest
{
protected:
  vector<size_t> counts1_, counts2_;
  double statistic_, pvalue_;

public:
  SubstitutionCountsComparison(): counts1_(), counts2_(), statistic_(0), pvalue_(1.) {}

  SubstitutionCountsComparison* clone() const = 0;

public:
  void setCounts(const vector<size_t>& c1, const vector<size_t>& c2) {
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

public:
  static double logFact(size_t x) {
    double f = 0;
    for (size_t i = 1; i <= x; ++i)
      f += log(static_cast<double>(i));
    return f;
  }

  static double multinomLogL(const vector<size_t>& v, const vector<double>& p);

  static double multinomialTest(const vector< vector<size_t> >& counts);
    
};


class AutomaticGroupingCondition
{
public:
  virtual ~AutomaticGroupingCondition() {}

public:
  virtual bool check(const vector<size_t>& counts) const = 0;
};

class NoAutomaticGroupingCondition:
  public AutomaticGroupingCondition
{
public:
  ~NoAutomaticGroupingCondition() {}

public:
  bool check(const vector<size_t>& counts) const {
    return true;
  }
};

class SumCountsAutomaticGroupingCondition:
  public AutomaticGroupingCondition
{
private:
  size_t threshold_;
  vector<size_t> ignore_;

public:
  SumCountsAutomaticGroupingCondition(size_t threshold = 0, const vector<size_t>& ignoreV = std::vector<size_t>(0)):
    threshold_(threshold),
    ignore_(ignoreV)
  {}

  ~SumCountsAutomaticGroupingCondition() {}

public:
  bool check(const vector<size_t>& counts) const {
    size_t s = 0;
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
  size_t threshold_;

public:
  AnyCountAutomaticGroupingCondition(size_t threshold = 0):
    threshold_(threshold)
  {}

  ~AnyCountAutomaticGroupingCondition() {}

public:
  bool check(const vector<size_t>& counts) const {
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
  vector< vector<size_t> > counts_;
  bool neighborsOnly_;
  bool negativeBrlen_;
  shared_ptr<SubstitutionCountsComparison> test_;
  vector<double> pvalues_;
	
public:
  /**
   * @param counts Total counts for each type, for each branch.
   */
  MultinomialClustering(
    const vector< vector<size_t> >& counts,
    const vector<int>& ids,
    const Tree& tree,
    const AutomaticGroupingCondition& autoGroup,
    bool neighborsOnly = false,
    bool negativeBrlen = false,
    bool verbose = false);
		
  virtual ~MultinomialClustering() {}

  MultinomialClustering* clone() const { return new MultinomialClustering(*this); }

public:
  TreeTemplate<Node>* getTree() const;
  /**
   * @return The list of pvalues, by order of clusters.
   */
  const vector<double> getPValues() const { return pvalues_; }

  std::string getName() const { return "Multinomial clusturing."; }

protected:
  /**
   * @brief Returns the pair with minimum distance and actualizes the vectors.
   *
   * The vector at position bestPair[0] is now the sum of vectors bestPair[0] and bestPair[1].
   * It is then used for computation of distances.
   */
  vector<size_t> getBestPair();
  vector<double> computeBranchLengthsForPair(const vector<size_t>& pair);
  double computeDistancesFromPair(const vector<size_t>& pair, const vector<double>& branchLengths, size_t pos);
  void finalStep(int idRoot);	
  virtual Node* getLeafNode(int id, const string& name);
  virtual Node* getParentNode(int id, Node* son1, Node* son2);
  double getDist(const vector<size_t>& v1, const vector<size_t>&v2);
};

#endif //_MULTINOMIALCLUSTERING_H_

