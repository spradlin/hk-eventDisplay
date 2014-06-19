#ifndef THKGAMMA_H
#define THKGAMMA_H
#include <TEveLine.h>
#include <TEveStraightLineSet.h>
/*
Classes to allow modification of usual TEveLine behaviour
*/
/// When a THKGamma is selected so are its parents AND its daughters
class THKGamma: 
public TEveLine
{
protected:
    float theta,phi,momentum,x,y,z;
public:
    THKGamma(const char* n = "THKGamma")
    : TEveLine( n ) 
    { };        
    /// When this item is selected, select also its parent.
    void FillImpliedSelectedSet(Set_t& impSelSet);
    /// Standard creation method for a straight line of a given length/theta/phi and starting
    /// at a particular position with a certain colour.
    void Create(TString title,Color_t color,double pos[3],double length,double Momentum,double theta,double phi);
    void Describe();
    
private:
    ClassDef(THKGamma, 0)  
    
};
/// When a THKLine is selected, so are its daughters
class THKLine: 
public THKGamma
{
    
public:
    THKLine(const char* n = "THKLine")
    : THKGamma( n ) 
    { }; 
    
    /// When this item is selected, select all its daughters.
    void FillImpliedSelectedSet(Set_t& impSelSet);
    
private:
    ClassDef(THKLine, 0)  
    
};
/// When a THKCerenkov is selected so are its parents
class THKCerenkov: 
public THKGamma
{
    
public:
    THKCerenkov(const char* n = "THKCerenkov")
    : THKGamma( n ) 
    { };        
    /// When this item is selected, select also its parent.
    void FillImpliedSelectedSet(Set_t& impSelSet);
    void Describe();
  
private:
    ClassDef(THKCerenkov, 0)  
    
};
class THKMCTrack: public TEveLine
{
public:
      THKMCTrack(const char* n = "THKMCTrack")
    : TEveLine( n ) 
    { }; 
    float start[3],stop[3],mass,momentum,energy,theta,phi;
    void SetValues(float Start[3],float Stop[3],float mass,float momentum,
    	float energy,float theta,float phi);
    void Describe();
     ClassDef(THKMCTrack, 0)  
};
class THKCerenkov2D: public TEveStraightLineSet
{
public:
      THKCerenkov2D(const char* n = "THKCerenkov2D")
    : TEveStraightLineSet( n ) 
    { }; 
    float start[3],momentum,theta,phi;
    void SetValues(double Start[3],double momentum,double theta,double phi);
    void Describe();
     ClassDef(THKCerenkov2D, 0)  
};
#endif
