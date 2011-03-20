#include "MultinomialClustering.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/AutoParameter.h>
//#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

double SimpleSubstitutionCountsComparison::multinomLogL(const vector<unsigned int>& v, const vector<double>& p)
{
  unsigned int n = VectorTools::sum(v);
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
      l += v[i] * log(static_cast<double>(p[i]));
      l -= logFact(v[i]);
    }
  }
  return l;
}

void SimpleSubstitutionCountsComparison::computePValue() {
  unsigned int s1 = VectorTools::sum(counts1_);
  unsigned int s2 = VectorTools::sum(counts2_);
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
  pvalue_ = 1. - RandomTools::pChisq(statistic_, counts1_.size() - 1);
}

EquilibriumSubstitutionCountsComparison::LikelihoodFunction::LikelihoodFunction(
    const RowMatrix<unsigned int>& na_arg,
    const RowMatrix<unsigned int>& nb_arg,
    double init_alpha,
    const vector<double>& init_r): 
  AbstractParametrizable(""),
  s_(na_arg.getNumberOfRows()),
  mlnL_(0),  mdlnLdalpha_(0), md2lnLdalpha2_(0),
  mdlnLdr_(s_), md2lnLdr2_(s_),
  na_(na_arg), nb_(nb_arg),
  sna_(s_), snb_(s_),
  nca_(s_), ncb_(s_),
  Sna_(0), Snb_(0), Sncb_(0),
  Pa_(s_), Pb_(s_),
  g_(s_, s_),
  r_(s_), alpha_(0),
  pa_(s_, s_), pb_(s_, s_)
{
  if ((na_arg.getNumberOfRows() != nb_arg.getNumberOfRows())
   || (na_arg.getNumberOfRows() != na_arg.getNumberOfColumns())
   || (nb_arg.getNumberOfRows() != nb_arg.getNumberOfColumns()))
    throw Exception("EquilibriumSubstitutionCountsComparison(constructor). Input count vectors should have the same size.");
  
  //define parameters to be optimized
  for (size_t i = 0; i < s_; i++) {
    string paramName = TextTools::toString(i);
    addParameter_(Parameter(paramName, init_r[i], &Parameter::PROP_CONSTRAINT_EX));
  }
  addParameter_(Parameter("alpha", init_alpha, &Parameter::R_PLUS_STAR));
  //calculate analytical parameters
  for (size_t i = 0 ; i < s_; ++i) {
    sna_[i] = snb_[i] = nca_[i] = ncb_[i] = 0;
    for (size_t j = 0; j < s_; ++j) {
      sna_[i] += na_(i, j);
      snb_[i] += nb_(i, j);
      if (j != i) {
        nca_[i] += na_(i, j);
        ncb_[i] += nb_(i, j);
        Sncb_   += nb_(i, j);
      }
      Sna_ += na_(i, j);
      Snb_ += nb_(i, j);
    }
  }

  for (size_t i = 0; i < s_; ++i) {
    Pa_[i] = static_cast<double>(sna_[i]) / static_cast<double>(Sna_);
    Pb_[i] = static_cast<double>(snb_[i]) / static_cast<double>(Snb_);
  }

  for (size_t i = 0; i < s_; ++i) {
    for (size_t j = 0; j < s_; ++j) {
      if (j != i)
        g_(i, j) = static_cast<double>(na_(i, j) + nb_(i, j)) / static_cast<double>(sna_[i] - na_(i, i) + snb_[i] - nb_(i, i));
    }
  }
  fireParameterChanged(getParameters());
}


void EquilibriumSubstitutionCountsComparison::LikelihoodFunction::fireParameterChanged(const ParameterList& pl)
{
  for (size_t i = 0; i < s_; ++i) {
    string paramName = TextTools::toString(i);
    r_[i] = getParameterValue(paramName);
  }
  alpha_ = getParameterValue("alpha");

  mlnL_ = -logLikelihoodEq(na_, nb_, Pa_, Pb_, alpha_, r_, g_, pa_, pb_);
  mdlnLdalpha_ = -dLogLikelihoodDalpha(nb_, Sncb_,  alpha_, r_);
  md2lnLdalpha2_ = -d2LogLikelihoodDalpha2(nb_, Sncb_,  alpha_, r_);
  for (size_t i = 0; i < s_; i++) {
    mdlnLdr_[i]   = -dLogLikelihoodDri(na_(i, i), nb_(i, i), nca_[i], ncb_[i], alpha_, r_[i]);
    md2lnLdr2_[i] = -d2LogLikelihoodDri2(na_(i, i), nb_(i, i), nca_[i], ncb_[i], alpha_, r_[i]);
  }
}

double EquilibriumSubstitutionCountsComparison::logLikelihoodFree(
    const RowMatrix<unsigned int>& na,
    const RowMatrix<unsigned int>& nb,
    const RowMatrix<double>& pa,
    const RowMatrix<double>& pb)
{
  size_t i, j, s = na.getNumberOfRows();
  double lnl = 0;

  for ( i = 0; i < s; ++i) {
    for (j = 0; j < s; ++j) {
      if (pa(i, j) == 0.) {
        if (na(i, j) != 0)
          return log(0.);
      } else {
        lnl += na(i, j) * log(pa(i, j));
      }
      if (pb(i, j) == 0.) {
        if (nb(i, j) != 0)
          return log(0.);
      } else {
        lnl += nb(i, j) * log(pb(i, j));
      }
    }
  }

  //MatrixTools::print(na);
  //MatrixTools::print(nb);
  //MatrixTools::print(pa);
  //MatrixTools::print(pb);
  //cout << lnl << endl;
  return lnl;
}

double EquilibriumSubstitutionCountsComparison::LikelihoodFunction::d2LogLikelihoodDri2(
    unsigned int naii,
    unsigned int nbii,
    unsigned int ncai,
    unsigned int ncbi,
    double alpha,
    double ri)
{
  double d2lnl = 0;

  d2lnl -= static_cast<double>(naii) / ((1 - ri) * (1 - ri));
  d2lnl -= static_cast<double>(nbii) * alpha * alpha/((1 - alpha * ri) * (1 - alpha * ri));
  d2lnl -= static_cast<double>(ncai + ncbi) / (ri * ri);

  return d2lnl;
}

double EquilibriumSubstitutionCountsComparison::LikelihoodFunction::dLogLikelihoodDri(
    unsigned int naii,
    unsigned int nbii,
    unsigned int ncai,
    unsigned int ncbi,
    double alpha,
    double ri)
{
  double dlnl = 0;

  dlnl -= static_cast<double>(naii) / (1 - ri);
  dlnl -= static_cast<double>(nbii) * alpha / (1 - alpha * ri);
  dlnl += static_cast<double>(ncai + ncbi) / ri;

  return dlnl;
}

double EquilibriumSubstitutionCountsComparison::LikelihoodFunction::d2LogLikelihoodDalpha2(
    const RowMatrix<unsigned int>& nb,
    unsigned int Sncb,
    double alpha,
    const vector<double>& r)
{
  double d2lnl = 0;
  double umalphari;

  d2lnl =- static_cast<double>(Sncb) / (alpha*alpha);
  for (unsigned int i = 0; i < nb.getNumberOfRows(); ++i) {
    umalphari = 1 - alpha * r[i];
    d2lnl -= static_cast<double>(nb(i, i)) * r[i] * r[i] / (umalphari * umalphari);
  }

  return d2lnl;
}

double EquilibriumSubstitutionCountsComparison::LikelihoodFunction::dLogLikelihoodDalpha(
    const RowMatrix<unsigned int>& nb,
    unsigned int Sncb,
    double alpha,
    const vector<double>& r)
{
  double dlnl = 0;

  dlnl = static_cast<double>(Sncb) / alpha;
  for (unsigned int i = 0; i < nb.getNumberOfRows(); ++i)
    dlnl -= static_cast<double>(nb(i, i)) * r[i] / (1 - alpha * r[i]);

  return dlnl;
}

double EquilibriumSubstitutionCountsComparison::LikelihoodFunction::logLikelihoodEq(
    const RowMatrix<unsigned int>& na,
    const RowMatrix<unsigned int>& nb,
    const vector<double>& Pa,
    const vector<double>& Pb,
    double alpha,
    const vector<double>& r,
    const RowMatrix<double>& g,
    RowMatrix<double>& pa,
    RowMatrix<double>& pb)
{
  double alphari;

  /* check all alpha*r[i]<1 */
  unsigned int i, j;
  for (i = 0; i < na.getNumberOfRows(); ++i) {
    alphari = alpha * r[i];
    if (alphari >= 1) break;
  }

  /* correct alpha if needed */
  if (i != na.getNumberOfRows()) {
    double maxr = r[0];
    unsigned int maxr_rank = 0;
    printf("correcting alpha\n");
    for (j = 1; j < na.getNumberOfRows(); ++j)
      if (r[j] > maxr) {
        maxr = r[j];
        maxr_rank = j;
      }
    alpha = 1 / maxr - 0.000001;
  }

  /* calculate lnL */
  for (i = 0; i < na.getNumberOfRows(); ++i) {
    alphari = alpha * r[i];
    for (j = 0; j < na.getNumberOfRows(); ++j) {
      if (j == i) {
        pa(i, j) = Pa[i] * (1 - r[i]);
        pb(i, j) = Pb[i] * (1 - alphari);
      }
      else {
        pa(i, j) = Pa[i] * r[i] * g(i, j);
        pb(i, j) = Pb[i] * alphari * g(i, j);
      }
    }
  }

  return logLikelihoodFree(na, nb, pa, pb);
}



void EquilibriumSubstitutionCountsComparison::computePValue()
{
  //Only for GC mapping for now!!!
  size_t s = 2;
  RowMatrix<unsigned int> na(s, s), nb(s, s);
  na(0, 0) = counts1_[3];
  na(0, 1) = counts1_[0];
  na(1, 0) = counts1_[1];
  na(1, 1) = counts1_[2];
  
  nb(0, 0) = counts2_[3];
  nb(0, 1) = counts2_[0];
  nb(1, 0) = counts2_[1];
  nb(1, 1) = counts2_[2];

  unsigned int Sna(0), Snb(0);
  vector<unsigned int> sna(s), snb(s);
  RowMatrix<double> free_pa(s, s), free_pb(s, s);
  double alpha_num(0), alpha_denom(0), init_alpha;
  vector<double> init_r(s);
  RowMatrix<unsigned int> arg1(s, s), arg2(s, s);
  double lnL_free, lnL_eq;

  size_t i, j;

  /* free model */
  for (i = 0; i < s; ++i) {
    for (j = 0; j < s; ++j) {
      Sna += na(i, j);
      Snb += nb(i, j);
    }
  }

  for (i = 0; i < s; ++i) {
    for (j = 0; j < s; ++j) {
      free_pa(i, j) = static_cast<double>(na(i, j)) / static_cast<double>(Sna);
      free_pb(i, j) = static_cast<double>(nb(i, j)) / static_cast<double>(Snb);
    }
  }

  lnL_free = logLikelihoodFree(na, nb, free_pa, free_pb);



  /* equilibrium model */

  //initial conditions
  for (i = 0; i < s; ++i) {
    sna[i] = snb[i] = 0;
    for (j = 0; j < s; ++j) {
      sna[i] += na(i, j);
      snb[i] += nb(i, j);
    }
  }

  for (i = 0; i < s; ++i) {
    alpha_num   += static_cast<double>(sna[i] * (snb[i] - nb(i, i)));
    alpha_denom += static_cast<double>(snb[i] * (sna[i] - na(i, i)));
  }
  init_alpha = alpha_num / alpha_denom;

  for (i = 0; i < s; ++i) {
    init_r[i] = static_cast<double>(sna[i] - na(i, i) + snb[i] - nb(i, i))/static_cast<double>(sna[i] + snb[i] * init_alpha);
    if (init_r[i] == 0.) {
      init_r[i] = 0.0001;
    }
  }

  //  (switch na and nb if needed so that alpha < 1)
  if (init_alpha <= 1.) {
    arg1 = na;
    arg2 = nb;
  } else {
    init_alpha = 1. / init_alpha;
    arg1 = nb;
    arg2 = na;
  }

  if (init_alpha == 0.)
    init_alpha += 0.001;

  // here create optimizer 
  LikelihoodFunction lik(na, nb, init_alpha, init_r);
  PseudoNewtonOptimizer opt(&lik);
  opt.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  opt.setVerbose(0);
  opt.setProfiler(0);
  opt.setMessageHandler(0);
  opt.init(lik.getParameters());
  opt.optimize();

  lnL_eq = -opt.getFunctionValue();

  if (lnL_free < lnL_eq)
    throw Exception("We have a problem here!");

  //cout << lnL_free << "\t" << lnL_eq << endl;
  statistic_ = 2 * (lnL_free - lnL_eq);
  pvalue_ = 1. - RandomTools::pChisq(statistic_, s*s - s - 1);
  //cout << statistic_ << "\t" << pvalue_ << endl;
}


short MultinomialClustering::CLUSTERING_SIMPLE = 0;
short MultinomialClustering::CLUSTERING_EQUILIBRIUM = 1;

MultinomialClustering::MultinomialClustering(
    const vector< vector<unsigned int> >& counts,
    const vector<int>& ids,
    const Tree& tree,
    const AutomaticGroupingCondition& autoGroup,
    short clusterType,
    bool neighborsOnly, bool negativeBrlen, bool verbose):
  AbstractAgglomerativeDistanceMethod(verbose), counts_(counts), neighborsOnly_(neighborsOnly), negativeBrlen_(negativeBrlen), test_()
{
  if (clusterType == CLUSTERING_SIMPLE)
    test_.reset(new SimpleSubstitutionCountsComparison());
  else if (clusterType == CLUSTERING_EQUILIBRIUM)
    test_.reset(new EquilibriumSubstitutionCountsComparison());
  else throw Exception("MultinomialClustering(constructor). Invalid clustering option.");
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

double MultinomialClustering::getDist(const vector<unsigned int>& v1, const vector<unsigned int>&v2)
{
  test_->setCounts(v1, v2);
  return 1. - test_->getPValue();
}

