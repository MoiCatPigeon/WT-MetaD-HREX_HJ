/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/CoordinationBase.h"
#include "tools/SwitchingFunction.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/IFile.h"

#include <iostream>

#include <string>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR GHBFIX
/*
Calculate ghbfix interaction energy among GROUPA and GROUPB.

This variable calculates the ghbfix interaction among GROUPA and GROUPB
using a potential defined in ...

This collective variable can be used to analyze or induce hydrogen bond interactions.
Notice that the value of the GHBFIX is returned in plumed units (see \ref UNITS).



\par Examples

\plumedfile
# this is printing the GHBFIX interaction between two groups of atoms
gh: GHBFIX PAIR GROUPA=1,3,4 GROUP=5,6,7 TYPES=types.dat PARAMS=params.dat
PRINT ARG=gh
\endplumedfile

*/
//+ENDPLUMEDOC

class GHBFIX : public CoordinationBase {
    
    double dmax;
    double d0;
    double c;   
    
    std::vector<unsigned> typesTable;
    std::vector<double> etas;
    
    unsigned n;
    bool nlist_flag {false};
    
    double dmax2;
    double A;
    double B;
    double C;
    double D;

public:
  explicit GHBFIX(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  double pairing(double distance,double&dfunc,unsigned i,unsigned j)const override;    
};

PLUMED_REGISTER_ACTION(GHBFIX,"GHBFIX")

void GHBFIX::registerKeywords( Keywords& keys ) {
  CoordinationBase::registerKeywords(keys);
    
  keys.add("compulsory","TYPES","the value of TYPES in the switching function");
  keys.add("compulsory","PARAMS","the value of PARAMS in the switching function");
  keys.add("compulsory","D_MAX","the value of D_MAX in the switching function");
  keys.add("compulsory","D_0","the value of D_0 in the switching function");
  keys.add("compulsory","C","the value of C in the switching function");
  keys.addFlag("ENABLE_NLIST",false,"use functional form adequate for enabling neighbor list");
}

GHBFIX::GHBFIX(const ActionOptions&ao):
  Action(ao),
  CoordinationBase(ao)
{
  std::string types;
  std::string params;
  parse("D_MAX",dmax);
  parse("D_0",d0);
  parse("C",c);
  parse("TYPES",types);
  parse("PARAMS",params);
  parseFlag("ENABLE_NLIST",nlist_flag);
      
  //const calculated once
  dmax2 = dmax-d0;
  
  A = (-c*dmax2*dmax2)/((1-c)*dmax2*dmax2);
  B = (2*dmax2)/((1-c)*dmax2*dmax2);    
  C = -1/((1-c)*dmax2*dmax2);
  D = 1/(c*dmax2*dmax2);
    
  //setup typesTable 
  IFile typesfile;
  typesfile.link(*this);
  typesfile.open(types);
  int itype;
  while(typesfile.scanField("itype",itype).scanField()) {
      plumed_assert(itype>=0)<<"itype ="<<itype<<", should be non-negative";
      typesTable.push_back(itype);
  }
      
  n = (int)*max_element(std::begin(typesTable), std::end(typesTable));
  n+=1;
      
  //scalingParameters  
  etas.resize(n*n,0.0);
  IFile etafile;
  etafile.open(params);
  int it,jt;
  double eta;
  while(etafile.scanField("itype",it).scanField("jtype",jt).scanField("eta",eta).scanField()){
      plumed_assert(it>=0)<<"itype ="<<it<<", should be non-negative";
      plumed_assert(jt>=0)<<"jtype ="<<jt<<", should be non-negative";
      plumed_assert(it<n)<<"itype ="<<it<<", should be smaller than "<<n;
      plumed_assert(jt<n)<<"jtype ="<<jt<<", should be smaller than "<<n;
      etas[n*it+jt]=eta;
  }
      
}
    
    
double GHBFIX::pairing(double distance2,double&dfunc,unsigned i,unsigned j)const {
        
    const auto i1=getAbsoluteIndex(i).index();
    plumed_assert(i1<typesTable.size())<<"your types table only covers "<<typesTable.size()<<" atoms, but you are trying to access atom number "<<(i1+1);
    const auto t1=typesTable[i1];
    
    const auto i2=getAbsoluteIndex(j).index();
    plumed_assert(i2<typesTable.size())<<"your types table only covers "<<typesTable.size()<<" atoms, but you are trying to access atom number "<<(i2+1);
    const auto t2=typesTable[i2];
    
    const double scale=etas[n*t1+t2]; 

    double distance=std::sqrt(distance2);
    double result;
    
    const double rdist = (distance-d0);
    
    if(distance>dmax) {
      if (nlist_flag==false) result=4.184*scale;
      else result=1*scale;
      dfunc=0.0;
      return result;
    }
    
    if(rdist<=0.) {
     result=0.0;
     dfunc=0.0;
     } else {
    result=0.0;
    dfunc=0.0; 

    if (rdist > c*dmax2) {
    result+=1-(A + B*rdist + C*rdist*rdist);
    dfunc-=B+2*C*rdist;
    } else if (rdist > 0.0) {
    result+=1-D*(rdist*rdist);
    dfunc-=2*D*rdist; 
    }
    
    if (nlist_flag==false) {    
        result-=1.;
        result*=-4.184;
        dfunc*=-4.184;
    }
        
    dfunc/=distance;
        
    result*=scale;
    dfunc*=scale;
        
  }
  return result;
}

}

}
