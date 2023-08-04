#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!


typedef float scalar;
typedef string word;

/****************************************************** Functions ***********************************************/
vector<string> lookupDictStringLine(word dictName, word lineName)
{
  ifstream dict(dictName);
  if(!dict)
  {
    cout << dictName << "could not be opened!" << endl;
  }

  string line;
  string stringLine;
  while(getline(dict, line))
  {
    if (line.find(lineName) != std::string::npos)
    {
      stringLine = line;
    }
  }
  dict.close();


  vector<string> stringList;
  string word = "";
  for (auto x : stringLine)
  {
    if ( x != ' ' && x != ';' && x != '\t')
    {
      word = word + x;
    }
    else
    {
      if(word != "")
      {
        stringList.push_back(word);
        word = "";
      }
    }
  }

  return stringList;
}


scalar lookupDictScalar(word dictName, word wordName)
{
  return stof(lookupDictStringLine(dictName, wordName)[1]);
}


int lookupDictInt(word dictName, word wordName)
{
  return stoi(lookupDictStringLine(dictName, wordName)[1]);
}


word lookupDictWord(word dictName, word wordName)
{
  return lookupDictStringLine(dictName, wordName)[1];
}


class material
{
  private:
    vector<string> materialProp;

  public:
    //- Constructor
    material (word dictName, word lineName)
    :
    materialProp(lookupDictStringLine(dictName, lineName))
    {}

    //- Destructor
    virtual ~material(){}

    //- Member Functions
      //- name
      virtual word name() const {return materialProp[0];}

      //- Density [kg/m^3]
      virtual scalar rho() const {return stof(materialProp[1]);}

      //- SpecificHeatCapacity [J/(Kg.K)]
      virtual scalar Cp() const {return stof(materialProp[2]);}

      //- ThermalConductivity [W/(m.K)]
      virtual scalar kappa() const {return stof(materialProp[3]);}

      //- DynamicViscosity [Pa.s]
      virtual scalar mu() const {return stof(materialProp[4]);}

      //- ThermalExpansionCoefficient [1/K]
      virtual scalar beta() const {return stof(materialProp[5]);}
};




/****************************************************** Properties ***********************************************/
word thermophysicalPropertiesDict_ = "thermophysicalProperties.c";


//- baseFluid
word baseFluidName_ = "water";
material baseFluid_(thermophysicalPropertiesDict_, baseFluidName_);
material* baseFluidPtr_ = &baseFluid_;
scalar rhobf_ = baseFluidPtr_->rho();
scalar Cpbf_ = baseFluidPtr_->Cp();
scalar kappabf_ = baseFluidPtr_->kappa();
scalar mubf_ = baseFluidPtr_->mu();
scalar betabf_ = baseFluidPtr_->beta();
scalar nubf_ = mubf_ / rhobf_;
scalar alphabf_ = kappabf_ / (rhobf_ * Cpbf_);
scalar Prbf_ = nubf_ / alphabf_;


//- nanoParticle1
word NP1Name_ = "SiO2";
material nanoParticle1_(thermophysicalPropertiesDict_, NP1Name_);
material* nanoParticle1Ptr_ = &nanoParticle1_;
scalar rhonp1_ = nanoParticle1Ptr_->rho();
scalar Cpnp1_ = nanoParticle1Ptr_->Cp();
scalar kappanp1_ = nanoParticle1Ptr_->kappa();
scalar munp1_ = nanoParticle1Ptr_->mu();
scalar betanp1_ = nanoParticle1Ptr_->beta();
scalar lambdanp1_ = lookupDictScalar(thermophysicalPropertiesDict_, "lambdanp1");
scalar phinp1_ = lookupDictScalar(thermophysicalPropertiesDict_, "phinp1");


//- nanoParticle2
word NP2Name_ = "MWCNT";
material nanoParticle2_(thermophysicalPropertiesDict_, NP2Name_);
material* nanoParticle2Ptr_ = &nanoParticle2_;
scalar rhonp2_ = nanoParticle2Ptr_->rho();
scalar Cpnp2_ = nanoParticle2Ptr_->Cp();
scalar kappanp2_ = nanoParticle2Ptr_->kappa();
scalar munp2_ = nanoParticle2Ptr_->mu();
scalar betanp2_ = nanoParticle2Ptr_->beta();
scalar lambdanp2_ = lookupDictScalar(thermophysicalPropertiesDict_, "lambdanp2");
scalar phinp2_ = lookupDictScalar(thermophysicalPropertiesDict_, "phinp2");


//- hybrid nanoFluid
scalar phinp_ = phinp1_ + phinp2_;
scalar rho_ = (1 - phinp_) * rhobf_ + phinp1_ * rhonp1_ + phinp2_ * rhonp2_;
scalar Cp_ = ((1 - phinp_) * rhobf_ * Cpbf_ + phinp1_ * rhonp1_ * Cpnp1_ + phinp2_ * rhonp2_ * Cpnp2_)/rho_;
scalar kappa_ = ( ( kappanp1_ + ( lambdanp1_ - 1.0 ) * kappabf_ - ( lambdanp1_ - 1.0 ) * ( kappabf_ - kappanp1_ ) * phinp1_ ) / ( kappanp1_ + ( lambdanp1_ - 1.0 ) * kappabf_ + ( kappabf_ - kappanp1_ ) * phinp1_ ) ) * ( ( kappanp2_ + ( lambdanp2_ - 1.0 ) * kappabf_ - ( lambdanp2_ - 1.0 ) * ( kappabf_ - kappanp2_ ) * phinp2_ ) / ( kappanp2_ + ( lambdanp2_ - 1.0 ) * kappabf_ + ( kappabf_ - kappanp2_ ) * phinp2_ ) ) * kappabf_;
scalar mu_ = mubf_ / pow( 1.0 - phinp_, 2.5);
scalar beta_ = ((1 - phinp_) * rhobf_ * betabf_ + phinp1_ * rhonp1_ * betanp1_ + phinp2_ * rhonp2_ * betanp2_)/rho_;
scalar nu_ = mu_ / rho_;
scalar alpha_ = kappa_ / (rho_ * Cp_);
scalar Pr_ = nu_ / alpha_;

//- 
//scalar CRa_ = (beta_ / betabf_) * (nubf_ / nu_) * (alphabf_ / alpha_);



/****************************************************** auxiliaries parameters ***********************************************/
scalar Tcold_ = lookupDictScalar(thermophysicalPropertiesDict_, "Tcold");
scalar Thot_ = lookupDictScalar(thermophysicalPropertiesDict_, "Thot");
//scalar length_ = lookupDictScalar(thermophysicalPropertiesDict_, "length");
scalar g_ = lookupDictScalar(thermophysicalPropertiesDict_, "g");
scalar Rabf_ = lookupDictScalar(thermophysicalPropertiesDict_, "Ra");
//scalar Rabf_ = (g_ * betabf_ * (Thot_ - Tcold_) * pow(length_, 3.0)) / (nubf_ * alphabf_);
scalar length_ = pow(((Rabf_ * nubf_ * alphabf_) / (betabf_ * (Thot_ - Tcold_) * g_)), (scalar) 1/3);

int pointNum_ = lookupDictInt(thermophysicalPropertiesDict_, "pointNum");
scalar endTime_ = lookupDictScalar(thermophysicalPropertiesDict_, "endTime");
scalar residual_ = lookupDictScalar(thermophysicalPropertiesDict_, "residual");





