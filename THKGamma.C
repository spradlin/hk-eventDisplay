#include "THKGamma.h"
#include <iostream>
using std::cout;
using std::endl;
ClassImp(THKLine)
void THKLine::FillImpliedSelectedSet(Set_t& impSelSet)
{
	TEveElement::FillImpliedSelectedSet(impSelSet);
	
	for (List_ci i=fChildren.begin(); i!=fChildren.end(); ++i)
	{
		impSelSet.insert(*i);
	}
}
ClassImp(THKCerenkov)
void THKCerenkov::FillImpliedSelectedSet(Set_t& impSelSet)
{
	TEveElement::FillImpliedSelectedSet(impSelSet);
	
	for (List_ci i=fParents.begin(); i!=fParents.end(); ++i)
	{
		impSelSet.insert(*i);
	}
}
void THKCerenkov::Describe()
{
	cout<<"Cerenkov ring from "<<GetTitle()<<endl;
	List_ci i=fParents.begin();
	TEveElement* obj=(*i);
	THKGamma *g = dynamic_cast<THKGamma*> (obj);
	if(g){
		g->Describe();
	}
	else
	{
		cout<<Form("User clicked on: \"%s\"", obj->GetElementName())<<endl;
		//cout<<" of class "<<obj->ClassElementName()<<endl;
	}
}
ClassImp(THKGamma)
/// Override selection, select a gamma and you select my whole family!
void THKGamma::FillImpliedSelectedSet(Set_t& impSelSet)
{
	TEveElement::FillImpliedSelectedSet(impSelSet);
	
	for (List_ci i=fParents.begin(); i!=fParents.end(); ++i)
	{
		impSelSet.insert(*i);
	}
	
	for (List_ci i=fChildren.begin(); i!=fChildren.end(); ++i)
	{
		impSelSet.insert(*i);
	}
}
void THKGamma::Create(TString title,Color_t color,double pos[3],double length,double Momentum,double Theta,double Phi){
	theta=Theta;
	phi=Phi;
	x=pos[0];
	y=pos[1];
	z=pos[2];
	momentum=Momentum;
	SetElementTitle(title);
	SetMainColor(color);
	SetNextPoint(pos[0],pos[1],pos[2]);
	double xend = pos[0]  + length*cos(phi)*sin(theta);
	double yend = pos[1]  + length*sin(phi)*sin(theta);
	double zend = pos[2]  + length*cos(theta);
	SetNextPoint(xend,yend,zend);
}
void THKGamma::Describe()
{
	cout<<GetTitle()<<": Momentum: "<<momentum<<" Theta: "<<theta<<" Phi: "<<phi<<" Vertex ("<<x<<","<<y<<","<<z<<")"<<endl;
}
ClassImp(THKMCTrack)
void THKMCTrack::SetValues(float Start[3],float Stop[3],float Mass,float Momentum,
    	float Energy,float Theta,float Phi)
{
	for(int i=0;i<3;i++)
	{
		start[i]=Start[i];
		stop[i]=Stop[i];
	}
	mass=Mass;
	momentum=Momentum;
	energy=Energy;
	theta=Theta;
	phi=Phi;
}
void THKMCTrack::Describe()
{
	cout<<" MC truth track: "<<GetTitle()<<endl;
	cout<<" Start: (";
	for(int i=0;i<3;i++){cout<<start[i];if(i<2)cout<<",";}
	cout<<")";
	cout<<" Stop: (";
	for(int i=0;i<3;i++){cout<<stop[i];if(i<2)cout<<",";}
	cout<<")";
	cout<<", Mass: "<<mass<<", Momentum: "<<momentum<<", Energy: "<<energy<<", Theta: "<<theta<<", Phi: "<<phi<<endl;	
}
ClassImp(THKCerenkov2D)
void THKCerenkov2D::SetValues(double Start[3],double Momentum,double Theta,double Phi)
{
	for(int i=0;i<3;i++)
	{
		start[i]=Start[i];
	}
	momentum=Momentum;
	theta=Theta;
	phi=Phi;
}
void THKCerenkov2D::Describe()
{
	cout<<GetTitle()<<": Momentum: "<<momentum<<" Theta: "<<theta<<" Phi: "<<phi<<" Vertex ("<<start[0]<<","<<start[1]<<","<<start[2]<<")"<<endl;
}

