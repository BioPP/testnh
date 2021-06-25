//
// File: MapNH.cpp
// Created by: Julien Dutheil
// Created on: Dec Thu 09 11:11 2010
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

//#include "MultinomialClustering.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

// From bpp:
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/OneProcessSequenceSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/PhyloMappings/SingleProcessSubstitutionMapping.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h>
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

//// utilitary function for Partition PhyloMappings

VVVdouble assignGoodSites(const VVVVdouble& vcount, const PartitionProcessPhyloLikelihood& ppl)
{
  VVVdouble count(ppl.getNumberOfSites());
  for (auto& c:count)
    c.resize(vcount[0][0].size());
  
  const std::vector<ProcPos>& vpp=ppl.getProcessSiteRelations();

  for (size_t i=0; i<count.size(); i++)
    count[i]=vcount[vpp[i].nProc-1][vpp[i].pos];
  return count;
}

VVdouble assignGoodSites(const VVdouble& vlength, const PartitionProcessPhyloLikelihood& ppl)
{
  VVdouble lengths(ppl.getNumberOfSites());
  for (auto& c:lengths)
    c.resize(vlength[0].size());
  
  const std::vector<ProcPos>& vpp=ppl.getProcessSiteRelations();

  for (size_t i=0; i<lengths.size(); i++)
    lengths[i]=vlength[vpp[i].nProc-1];

  return lengths;
}


//////////////////////////////////////

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     Map NH, version 2                          *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Modif. 17/12/11      *" << endl;
  cout << "*          L. Guéguen                       Last Modif. 28/11/17 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication mapnh(args, argv, "MapNH");
    if (args == 1)
    {
      mapnh.help("mapnh");
      exit(0);
    }
  
    mapnh.startTimer();
    std::map<std::string, std::string> unparsedParams;

    Context context;
    

    /*********************************/
    /* get Basic objects */
    /*********************************/
  
    unique_ptr<Alphabet> alphabet(mapnh.getAlphabet());
    unique_ptr<GeneticCode> gCode(mapnh.getGeneticCode(alphabet.get()));

    map<size_t, AlignedValuesContainer*> mSites = mapnh.getAlignmentsMap(alphabet.get(), true);

    unique_ptr<SubstitutionProcessCollection> SP(mapnh.getCollection(alphabet.get(), gCode.get(), mSites, unparsedParams));
                                               
    std::map<size_t, SequenceEvolution*> mProc = mapnh.getProcesses(*SP, unparsedParams);
  
    auto plc(mapnh.getPhyloLikelihoods(context, mProc, *SP, mSites));

    if (!plc->hasPhyloLikelihood(1))
      throw Exception("Missing first phyloLikelihood.");

    PhyloLikelihood* pl=(*plc)[1];

    mapnh.fixLikelihood(alphabet.get(), gCode.get(), pl);

    double thresholdSat = ApplicationTools::getDoubleParameter("count.max", mapnh.getParams(), -1, "", true, 1);
    if (thresholdSat > 0)
      ApplicationTools::displayResult("Saturation threshold used", thresholdSat);

  
    // Checks all phylolikelihoods fit, and only one alphabet state map

    OneProcessSequencePhyloLikelihood* opspl= dynamic_cast<OneProcessSequencePhyloLikelihood*>(pl);
    SingleProcessPhyloLikelihood* sppl= dynamic_cast<SingleProcessPhyloLikelihood*>(pl);    
    PartitionProcessPhyloLikelihood* sap = dynamic_cast<PartitionProcessPhyloLikelihood*>(pl);
  
    const StateMap* pstmap(0);

    if (opspl==NULL && sppl==NULL && sap==NULL)
      throw Exception("Mapping not possible for this phylo. Ask the developpers");

    if (opspl)
      pstmap=&opspl->getSubstitutionProcess().getStateMap();
    else if (sppl)
      pstmap=&sppl->getSubstitutionProcess().getStateMap();
    else
    {
      const std::vector<size_t>& vpn=sap->getNumbersOfPhyloLikelihoods();
      for (const auto& pn : vpn)
      {
        const AbstractPhyloLikelihood* pap=sap->getAbstractPhyloLikelihood(pn);
        const OneProcessSequencePhyloLikelihood* opspl2= dynamic_cast<const OneProcessSequencePhyloLikelihood*>(pap);
        const SingleProcessPhyloLikelihood* sppl2 = dynamic_cast<const SingleProcessPhyloLikelihood*>(pap);

        if (opspl2==NULL && sppl2==NULL)
          throw Exception("Mapping not possible for a non-single process in container phylo.");

        if (!pstmap)
          pstmap=opspl2?&opspl2->getSubstitutionProcess().getStateMap():&sppl2->getSubstitutionProcess().getStateMap();
        else
        {
          const StateMap* pstmap2=opspl2?&opspl2->getSubstitutionProcess().getStateMap():&sppl2->getSubstitutionProcess().getStateMap();
          if (pstmap2->getAlphabetStates()!=pstmap->getAlphabetStates())
            throw Exception("Discordant alphabet states in container phylo.");
        }
      }
    }
  

    //////////////////////////////////
    // set register and initialize the parameters for the mapping:
    //////////////
  
    string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "All", "", true, false);

    AlphabetIndex2* weights = 0;
    AlphabetIndex2* distances = 0;
  
    unique_ptr<SubstitutionRegister> reg(PhylogeneticsApplicationTools::getSubstitutionRegister(regTypeDesc, *pstmap, gCode.get(), weights, distances));

    shared_ptr<const AlphabetIndex2> sweights(weights);
    shared_ptr<const AlphabetIndex2> sdistances(distances);

    //Write categories:
    for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
      ApplicationTools::displayResult("  * Count type " + TextTools::toString(i + 1), reg->getTypeName(i + 1));

  
    // specific parameters to the null models
    string nullProcessParams = ApplicationTools::getStringParameter("nullProcessParams", mapnh.getParams(), "", "", false, 1);

    // output
  
    string outputDesc = ApplicationTools::getStringParameter("output.counts", mapnh.getParams(), "PerType(prefix=counts_)");

    string outputType;
    map<string, string> outputArgs;
    KeyvalTools::parseProcedure(outputDesc, outputType, outputArgs);
    
    bool perBranchLength(0);
    bool perWord(0);
    bool splitNorm(0);
  
    if (nullProcessParams != "")
    {
      splitNorm = ApplicationTools::getBooleanParameter("splitNorm", outputArgs, false, "", true, 1);

      ApplicationTools::displayResult("Display separate counts and normalizations", splitNorm?"true":"false");
    
      if (!splitNorm)
      {
        perBranchLength=ApplicationTools::getBooleanParameter("perBranchLength", outputArgs, true, "", true, 0);

        ApplicationTools::displayResult("Normalization per branch length", perBranchLength?"true":"false");

        perWord = ApplicationTools::getBooleanParameter("perWordSize", outputArgs, true, "", true, 0);
        
        ApplicationTools::displayResult("Normalization per word size", perWord?"true":"false");
      }
    }
    
    uint siteSize=(perWord && AlphabetTools::isWordAlphabet(alphabet.get()))?dynamic_cast<const CoreWordAlphabet*>(alphabet.get())->getLength():1;
    
    // Stock Phylolikelihoods, needed if perBranchLength 
    vector<const ParametrizablePhyloTree*> vpt;
    
    // Compute phylosubstitutionmapping
    vector<shared_ptr<PhyloSubstitutionMapping> > vpsm;
  
    if (opspl || sppl)
    {
      if (opspl)
      {
        shared_ptr<PhyloSubstitutionMapping> spsm(new OneProcessSequenceSubstitutionMapping(*opspl, *reg, sweights, sdistances));
        spsm->computeCounts(thresholdSat);
        vpsm.push_back(spsm);
        if (perBranchLength)
          vpt.push_back(&opspl->getTree());
      }
      else
      {
        shared_ptr<PhyloSubstitutionMapping> spsm(new SingleProcessSubstitutionMapping(*sppl, *reg, sweights, sdistances));
        spsm->computeCounts(thresholdSat);
        vpsm.push_back(spsm);
        if (perBranchLength)
          vpt.push_back(&sppl->getTree());
      }
    }
    else
    {
      const std::vector<size_t>& vpn=sap->getNumbersOfPhyloLikelihoods();
      size_t comp(0);
      
      for (const auto& pn : vpn)
      {
        OneProcessSequencePhyloLikelihood* opspl2= dynamic_cast<OneProcessSequencePhyloLikelihood*>(sap->getAbstractPhyloLikelihood(pn));
        SingleProcessPhyloLikelihood* sppl2= dynamic_cast<SingleProcessPhyloLikelihood*>(sap->getAbstractPhyloLikelihood(pn));
  
        if (opspl2)
        {
          shared_ptr<PhyloSubstitutionMapping> spsm(new OneProcessSequenceSubstitutionMapping(*opspl2, *reg, sweights, sdistances));
          spsm->computeCounts(thresholdSat,comp<10);
          vpsm.push_back(spsm);
          
          if (perBranchLength)
            vpt.push_back(&opspl2->getTree());
        }
        else
        {
          shared_ptr<PhyloSubstitutionMapping> spsm(new SingleProcessSubstitutionMapping(*sppl2, *reg, sweights, sdistances));
          spsm->computeCounts(thresholdSat,comp<10);
          vpsm.push_back(spsm);
          
          if (perBranchLength)
            vpt.push_back(&sppl2->getTree());
        }

        if (comp>=10)
          ApplicationTools::displayResult("Build Substitution count",TextTools::toString(comp+1));

        comp++;
      }
    }

    ///////////////////////////////////////////////
    //  Compute normalizations if needed

    for (size_t i=0; i<vpsm.size(); i++)
    {
      shared_ptr<PhyloSubstitutionMapping> psm=vpsm[i];

      ParameterList nullParams;
      if (nullProcessParams != "")
      {
        const ParameterList pl0=psm->getParameters();
      
        map<string, string> npv;
        KeyvalTools::multipleKeyvals(nullProcessParams, npv, ",", false);
      
        for (const auto& pv : npv)
        {
          vector<string> pn = pl0.getMatchingParameterNames(pv.first);
          double val=TextTools::toDouble(pv.second);
          for (const auto& n : pn)
          {
            nullParams.addParameter(Parameter(n, val));
            ApplicationTools::displayResult("null Parameter " + n, val);
          }
        }

        psm->computeNormalizations(nullParams,i<10);
        if (i>=10)
          ApplicationTools::displayResult("Build Substitution normalization",TextTools::toString(i+1));
      }
    }
    
  
    ////////////////////////////////////////////
    //// OUTPUT
    ////////////////////////////////////////////

    
    if (outputType.find("Site")==string::npos)
      // OUTPUT PER BRANCH PER TYPE
    {
      // Write count trees:
      string treePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_type_", "", true, 1);
      
      Newick newick;

      PhyloTree pht;               
      if (nullProcessParams!="" && !splitNorm && perBranchLength)
      {
        size_t nbs(0);
        
        for (size_t pr=0;pr<vpt.size();pr++)
        {
          PhyloTree pht2(*vpt[pr]);
          size_t ns=vpsm[pr]->getCounts().getNumberOfSites();
          pht2.scaleTree((double)(ns));
          nbs+=ns;
        // proportional to the number of sites
          
          if (pr==0)
            pht=pht2;
          else
            pht+=pht2;
        }
        if (perBranchLength)
          pht.scaleTree(1./(double)nbs);          
      }
 
      for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
      {
        string name=reg->getTypeName(i+1);
        if (name=="")
          name=TextTools::toString(i + 1);
        
        unique_ptr<PhyloTree> pt, pn, ptt, pnt;
        
        for (auto& psm:vpsm)
        {
          if (pt==0)
            pt.reset(SubstitutionMappingTools::getTreeForType(psm->getCounts(),i));
          else
          {
            ptt.reset(SubstitutionMappingTools::getTreeForType(psm->getCounts(),i));
            (*pt)+=(*ptt);
          }
          
          if (nullProcessParams!="")
          {
            if (pn==0)
              pn.reset(SubstitutionMappingTools::getTreeForType(psm->getNormalizations(),i));
            else
            {
              pnt.reset(SubstitutionMappingTools::getTreeForType(psm->getNormalizations(),i));
              (*pn)+=(*pnt);
            }
          }
        }
        
        if (splitNorm || nullProcessParams=="")
        {
          string path = treePathPrefix + name + string(".dnd");
          ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
          newick.writeTree(*pt, path);
          
          if (splitNorm)
          {
            path = treePathPrefix + name + string("_norm.dnd");
            ApplicationTools::displayResult(string("Output normalizations of type ") + TextTools::toString(i + 1) + string(" to file"), path);
            newick.writeTree(*pn, path);
          }
        }
        else
        {
          (*pt)/=(*pn);
          
          if (perBranchLength)
            (*pt)*=pht;
          pt->scaleTree(1./siteSize);
          
          string path = treePathPrefix + name + string(".dnd");
          ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
          newick.writeTree(*pt, path);
        }
      }
    }
    else
    {
      // ALL PER SITE COUNTS

      // Compute per site per branch per type counts
      VVVVdouble vcounts, vnorm;
      VVdouble vlength;

      bool perBranch=(outputType.find("Branch")!=string::npos);
      bool perType=(outputType.find("Type")!=string::npos);

      vcounts.resize(vpsm.size());
      if (nullProcessParams!="")
      {
        vnorm.resize(vpsm.size());
        if (perBranchLength)
          vlength.resize(vpsm.size());
      }
      
      Vuint ids=vpsm[0]->getCounts().getAllEdgesIndexes();

      for (size_t i=0;i<vpsm.size();i++)
      {
        shared_ptr<PhyloSubstitutionMapping> psm=vpsm[i];

        if (perBranch && i!=0 && ids!=psm->getCounts().getAllEdgesIndexes())
          throw Exception("Branch ids of phylolikelihood " + TextTools::toString(i) + " do not match the ones of the first phylolikelihood.");
          
        vcounts[i]=SubstitutionMappingTools::getCountsPerSitePerBranchPerType(psm->getCounts());
              
        if (nullProcessParams!="")
        {
          vnorm[i]=SubstitutionMappingTools::getCountsPerSitePerBranchPerType(psm->getNormalizations());
          if (perBranchLength)
            vlength[i]=vpt[i]->getBranchLengths();
        }
      }
      
      // build final vectors of counts
      
      VVVdouble counts, norm;
      VVdouble lengths;
      
      if (!sap)
      {
        counts=vcounts[0];
        if (nullProcessParams!="")
        {
          norm=vnorm[0];
          if (perBranchLength)
          {
            lengths.resize(vnorm[0].size());
            for (auto& l:lengths)
              l=vlength[0];
          }
        }
      }
      else
      {
        counts=assignGoodSites(vcounts,*sap);
        if (nullProcessParams!="")
        {
          norm=assignGoodSites(vnorm,*sap);
          if (perBranchLength)
            lengths=assignGoodSites(vlength,*sap);
        }
      }

      // output per site 

      if (!perBranch)
        // PER SITE PER TYPE UNIQUELY
      {
        string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_site_per_type.txt", "", true, 1);
        
        ApplicationTools::displayResult(string("Output counts (site/type) to file"), perSitenf);

        VVdouble counts2(counts.size());

        for (auto& c:counts2)
          c.resize(reg->getNumberOfSubstitutionTypes());

        for (size_t s=0; s<counts.size(); s++)
        {
          VVdouble& counts_s=counts[s];
          Vdouble& counts2_s=counts2[s];
          
          for (const auto& counts_s_br:counts_s)
            counts2_s+=counts_s_br;

          if (nullProcessParams!="" && !splitNorm)
          {
            VVdouble& norm_s=norm[s];
            Vdouble norm2(reg->getNumberOfSubstitutionTypes(),0);
            
            for (const auto& norm_s_br:norm_s)
              norm2+=norm_s_br;

            counts2_s/=(norm2*siteSize);

            if (perBranchLength)
              counts2_s*=VectorTools::sum(lengths[s]);
          }
        }

        SubstitutionMappingTools::outputPerSitePerType(perSitenf, *reg, counts2);
        
        if (nullProcessParams!="" && splitNorm)
        {
          perSitenf += "_norm";
              
          ApplicationTools::displayResult(string("Output normalizations (site/type) to file"), perSitenf);
              
          VVdouble norm2(norm.size());
          for (auto& n:norm2)
            n.resize(reg->getNumberOfSubstitutionTypes());
          
          for (size_t s=0; s<norm.size(); s++)
          {
            VVdouble& norm_s=norm[s];
            Vdouble& norm2_s=norm2[s];
          
            for (const auto& norm_s_br:norm_s)
              norm2_s+=norm_s_br;
          }
          SubstitutionMappingTools::outputPerSitePerType(perSitenf, *reg, norm2);
        }
      }
      else {
        if (!perType)
          // PER SITE PER BRANCH
        {
           string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_site_per_branch.txt", "", true, 1);
        
          ApplicationTools::displayResult(string("Output counts (branch/site) to file"), perSitenf);
          
          VVdouble counts2(counts.size());
          for (auto& c:counts2)
            c.resize(counts[0].size());
          
          for (size_t s=0; s<counts.size(); s++)
          {
            VVdouble& counts_s=counts[s];
            Vdouble& counts2_s=counts2[s];
          
            for (size_t br=0; br<counts_s.size(); br++)
              counts2_s[br]=VectorTools::sum(counts_s[br]);

            if (nullProcessParams!="" && !splitNorm)
            {
              VVdouble& norm_s=norm[s];
              Vdouble norm2(norm_s.size());
              
              for (size_t br=0; br<norm_s.size(); br++)
                norm2[br]=VectorTools::sum(norm_s[br]);

              counts2_s/=(norm2*siteSize);

              if (perBranchLength)
                counts2_s*=lengths[s];
            }
          }

          SubstitutionMappingTools::outputPerSitePerBranch(perSitenf, ids, counts2);

          if (nullProcessParams!="" && splitNorm)
          {
            perSitenf += "_norm";
            
            ApplicationTools::displayResult(string("Output normalizations (branch/site) to file"), perSitenf);
              
            VVdouble norm2(norm.size());
            for (auto& n:norm2)
              n.resize(norm[0].size());
            
            for (size_t s=0; s<norm.size(); s++)
            {
              VVdouble& norm_s=norm[s];
              Vdouble& norm2_s=norm2[s];
          
              for (size_t br=0; br<norm.size(); br++)
                norm2_s[br]=VectorTools::sum(norm_s[br]);
            }
            SubstitutionMappingTools::outputPerSitePerBranch(perSitenf, ids, norm2);
          }
        }
        else
        {
          // OUTPUT PER SITE PER TYPE PER BRANCH
          string tablePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_site_per_branch_per_type_", "", true, 1);
        
          ApplicationTools::displayResult(string("Output counts (site/branch/type) to files"), tablePathPrefix + "*");

          if (nullProcessParams!="" && !splitNorm)
          {
            for (size_t s=0; s<counts.size(); s++)
            {
              VVdouble& counts_s=counts[s];
              VVdouble& norm_s=norm[s];
              Vdouble& lengths_s=lengths[s];
              for (size_t br=0; br<norm_s.size(); br++)
              {
                Vdouble& counts_s_br=counts_s[br];                
                counts_s_br/=(norm_s[br]*siteSize);
                
                if (perBranchLength)
                  counts_s_br*=lengths_s[br];
              }
            }
          }

          SubstitutionMappingTools::outputPerSitePerBranchPerType(tablePathPrefix+"_", ids, *reg, counts);
        
          if (nullProcessParams!="" && splitNorm)
          {
            tablePathPrefix += "_norm_";
          
            ApplicationTools::displayResult(string("Output normalizations (site/branch/type) to files"), tablePathPrefix + "*");
          
            SubstitutionMappingTools::outputPerSitePerBranchPerType(tablePathPrefix, ids, *reg, norm);
          }
        }
      }
    }
    
    /////////////////////////////////
    // clean up

    for (auto it : mSites)
      delete it.second;
    
    mapnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }
  
  return 0;
  
}

