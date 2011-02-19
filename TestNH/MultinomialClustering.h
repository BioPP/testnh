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


class EquilibriumSubstitutionCountsComparison:
  public SubstitutionCountsComparison
{
  private:
    class LikelihoodFunction:
      public virtual DerivableSecondOrder,
      public AbstractParametrizable
    {
      private:
        size_t s_;
        double mlnL_, mdlnLdalpha_, md2lnLdalpha2_;
        vector<double> mdlnLdr_, md2lnLdr2_;

        // data
        RowMatrix<unsigned int> na_;
        RowMatrix<unsigned int> nb_;

        //data-related variables
        vector<unsigned int> sna_, snb_, nca_, ncb_;
        unsigned int Sna_, Snb_, Sncb_;

        //parameters (analytically calculated)
        vector<double> Pa_, Pb_;					/* Pa[i]=freq of state i in branch a*/
	      RowMatrix<double> g_;	  	/* g[i][j] (j!=i) = proba change to j knowing i changes */
							                    	   		/* (shared by a and b)*/
        //parameters (to be optimized)
        vector<double> r_;	   						/* r[i]=proba of a change in branch a knowing state i */
        double alpha_;						/* alpha*r[i]=proba of a change in branch b knowing state i */

        //helper variables
        RowMatrix<double> pa_, pb_;					/* probabilities of na and nb */

      public:
        LikelihoodFunction(
            const RowMatrix<unsigned int>& na_arg,
            const RowMatrix<unsigned int>& nb_arg,
            double init_alpha,
            const vector<double>& init_r);

        LikelihoodFunction* clone() const { return new LikelihoodFunction(*this); }

        void setParameters(const ParameterList& pl)
          throw (ParameterNotFoundException, ConstraintException, Exception)
        {
          matchParametersValues(pl);
        }

        double getValue() const throw (Exception) { return mlnL_; }

        //Always on
        void enableFirstOrderDerivatives(bool) {}
        bool enableFirstOrderDerivatives() const { return true; }
        void enableSecondOrderDerivatives(bool) {}
        bool enableSecondOrderDerivatives() const { return true; }

        double getFirstOrderDerivative(const std::string& paramName) const throw (Exception) {
          if (paramName == "alpha") return mdlnLdalpha_;
          else return mdlnLdr_[TextTools::toInt(paramName)];
        } 

        double getSecondOrderDerivative(const std::string& paramName) const throw (Exception) {
          if (paramName == "alpha") return md2lnLdalpha2_;
          else return md2lnLdr2_[TextTools::toInt(paramName)];
        }
        
        double getSecondOrderDerivative(const std::string& paramName1, const std::string& paramName2) const throw (Exception) {
          throw Exception("EquilibriumSubstitutionCountsComparison::LikelihoodFunction::getSecondOrderDerivative. Cross derivatives not implemented.");
        }

        void fireParameterChanged(const ParameterList& pl);

      public:
        static double d2LogLikelihoodDri2(
            unsigned int naii,
            unsigned int nbii,
            unsigned int ncai,
            unsigned int ncbi,
            double alpha,
            double ri);

        static double dLogLikelihoodDri(
            unsigned int naii,
            unsigned int nbii,
            unsigned int ncai,
            unsigned int ncbi,
            double alpha,
            double ri);

        static double d2LogLikelihoodDalpha2(
            const RowMatrix<unsigned int>& nb,
            unsigned int Sncb,
            double alpha,
            const vector<double>& r);

        static double dLogLikelihoodDalpha(
            const RowMatrix<unsigned int>& nb,
            unsigned int Sncb,
            double alpha,
            const vector<double>& r);

        static double logLikelihoodEq(
            const RowMatrix<unsigned int>& na,
            const RowMatrix<unsigned int>& nb,
            const vector<double>& Pa,
            const vector<double>& Pb,
            double alpha,
            const vector<double>& r,
            const RowMatrix<double>& g,
            RowMatrix<double>& pa,
            RowMatrix<double>& pb);
    };

  public:
    EquilibriumSubstitutionCountsComparison* clone() const { return new EquilibriumSubstitutionCountsComparison(*this); }

  public:
    string getName() const { return "Equilibrium multinomial test."; }

  protected:
    void computePValue();
    static double logLikelihoodFree(
        const RowMatrix<unsigned int>& na,
        const RowMatrix<unsigned int>& nb,
        const RowMatrix<double>& pa,
        const RowMatrix<double>& pb);

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

  public:
    SumCountsAutomaticGroupingCondition(unsigned int threshold = 0):
      threshold_(threshold)
    {}

    ~SumCountsAutomaticGroupingCondition() {}

  public:
    bool check(const vector<unsigned int>& counts) const {
      return VectorTools::sum(counts) > threshold_;
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
        short clusterType = CLUSTERING_SIMPLE,
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

  public:
    static short CLUSTERING_SIMPLE;
    static short CLUSTERING_EQUILIBRIUM;
};

#endif //_MULTINOMIALCLUSTERING_H_

