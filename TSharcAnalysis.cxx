#include "TSharcAnalysis.h"

#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>

ClassImp(TSharcAnalysis)

// we don't want this to run if you don't need it
TSRIM* TSharcAnalysis::p_in_targ									= 0;
TSRIM* TSharcAnalysis::p_in_si  									= 0;
TSRIM* TSharcAnalysis::d_in_targ									= 0;
TSRIM* TSharcAnalysis::d_in_si  									= 0;
TSRIM* TSharcAnalysis::c_in_targ									= 0;
TSRIM* TSharcAnalysis::c_in_si  									= 0;
TSRIM* TSharcAnalysis::a_in_targ									= 0;
TSRIM* TSharcAnalysis::a_in_si  									= 0;
TSRIM* TSharcAnalysis::beam_in_targ  							= 0;

// segmentation and pitch sizes (warning, QQQ back pitch is angular [deg])
int TSharcAnalysis::number_of_detectors   = 16;	
int TSharcAnalysis::frontstrips[16]       = {16,16,16,16,	24,24,24,24,	24,24,24,24,	16,16,16,16};
int TSharcAnalysis::backstrips[16]        = {24,24,24,24,	48,48,48,48,	48,48,48,48,	24,24,24,24};		
double TSharcAnalysis::frontpitches[16]   = {2.0,2.0,2.0,2.0,	3.0,3.0,3.0,3.0,	3.0,3.0,3.0,3.0,	2.0,2.0,2.0,2.0};
double TSharcAnalysis::backpitches[16]    = {3.4,3.4,3.4,3.4,	1.0,1.0,1.0,1.0,	1.0,1.0,1.0,1.0,	3.4,3.4,3.4,3.4}; 	
// various thicknesses in um
double TSharcAnalysis::DELdeadlayers[16]  = {0.7,0.7,0.7,0.7,	0.1,0.1,0.1,0.1,	0.1,0.1,0.1,0.1,	0.7,0.7,0.7,0.7};
double TSharcAnalysis::PADdeadlayers[16]  = {0.0,0.0,0.0,0.0,	1.0,1.0,1.0,1.0,	0.0,0.0,0.0,0.0,	0.0,0.0,0.0,0.0}; 
double TSharcAnalysis::DELthicknesses[16] = {998,998,998,1001,	141,142,133,143,	999,1001,1001,1002,	390,390,383,385}; 
double TSharcAnalysis::PADthicknesses[16] = {0.0,0.0,0.0,0.0,	1534,1535,1535,1539,	0.0,0.0,0.0,0.0,	0.0,0.0,0.0,0.0};
// target properties (um)
TVector3 TSharcAnalysis::position_offset					; 
TVector3 TSharcAnalysis::position_error						; 
double TSharcAnalysis::targetthickness            = 0.0; // THIS IS TRUE FOR 94Sr (2013) AND 95,96Sr (2014)   
double TSharcAnalysis::targetradius               = 2000; 
std::string TSharcAnalysis::targmat								= "";
// strips and acceptances
TReaction *TSharcAnalysis::reaction 							= 0;
TList *TSharcAnalysis::coveragelist 							= 0;
Int_t TSharcAnalysis::accres								 			= 0;
std::string TSharcAnalysis::badstripsfile					= "";
Bool_t TSharcAnalysis::resetbadstrips							= true;
UInt_t TSharcAnalysis::nbadstrips									= 0;
double TSharcAnalysis::Omega_min								 	= 0.0000;
double TSharcAnalysis::Omega_max 									= 0.0000;



TSharcAnalysis::TSharcAnalysis(double Xoff, double Yoff, double Zoff)	{	
  InitializeSRIMInputs();	
  SetTarget(Xoff,Yoff,Zoff);
  SetOmegaMaxMin(true);
}

TSharcAnalysis::~TSharcAnalysis()	{	}


void TSharcAnalysis::SetTarget(double x, double y, double z, double thickness, const char *material, double dx, double dy, double dz){
	
	TVector3 vec(x,y,z);
	SetTargetPosition(vec);
	
	TVector3 errvec(dx,dy,dz);
	SetTargetPositionError(errvec); // assign default ±1mm error to target position
	
	targetthickness = thickness;
	targmat.assign(material);
	
}


double TSharcAnalysis::GetDetectorThickness(double dist, int det, double theta, double phi)		{

	double phi_90 =  fmod(std::fabs(phi),PI/2); // get angle in range 0 to 90 degress
	double phi_45 = phi_90;
	if(phi_90>PI/4)
		phi_45 = PI/2-phi_90; // if it's greater than 45 degrees, set it to be 90 degrees minus ang_capped_90. 

	if(det>=5 && det<=12 ) // BOX
		return dist/(TMath::Sin(theta)*TMath::Cos(phi_45));
	else // QQQ
		return std::fabs(dist/(TMath::Cos(theta)));
}


double TSharcAnalysis::GetTargetThickness(double theta, double phi, double fractional_depth)		{ 
   double thickness;
   
   if(targetthickness==0){
   	printf("\n{TSharcAnalysis} Target thickness hasn't been set.\n\n");
   	return 0;
	}

   if(theta<PI/2) 
 		thickness = std::fabs((1-fractional_depth)*targetthickness/TMath::Cos(theta));
   else if(theta>PI/2) 
 		thickness = std::fabs(fractional_depth*targetthickness/TMath::Cos(theta));	

  if(thickness>targetradius || theta==PI/2)
    return targetradius; 
	else 
		return thickness;
}


void TSharcAnalysis::InitializeSRIMInputs()	{

	if(!targmat.size())
		printf("\n\tError :  Target Material Must Be Set!\n");

	p_in_targ  = new TSRIM(Form("p_in_%s.txt",targmat.c_str()),100e3,0,false);
	p_in_si    = new TSRIM("p_in_si.txt",100e3,0,false); 

	d_in_targ  = new TSRIM(Form("d_in_%s.txt",targmat.c_str()),100e3,0,false);
	d_in_si    = new TSRIM("d_in_si.txt",100e3,0,false);

	a_in_targ  = new TSRIM(Form("a_in_%s.txt",targmat.c_str()),200e3,0,false);
	a_in_si    = new TSRIM("a_in_si.txt",200e3,0,false);

	c_in_targ  = new TSRIM(Form("c_in_%s.txt",targmat.c_str()),300e3,0,false);
	c_in_si    = new TSRIM("c_in_si.txt",300e3,0,false);

}


TSRIM *TSharcAnalysis::GetSRIM(char ion, std::string material){

  if(p_in_targ==0)
    TSharcAnalysis::InitializeSRIMInputs();

  TSRIM *srim = 0;
  switch(ion){
    case 'p':
      if    (material.find("targ")==0)
        srim = p_in_targ;
      else if(material.find("si")==0)
        srim = p_in_si;
      else printf("{TSharcAnalysis} Error: Material '%s' not recognised.\n",material.c_str());
      break;
    case 'd':
      if    (material.find("targ")==0)
        srim = d_in_targ;
      else if(material.find("si")==0)
        srim = d_in_si;
      else printf("{TSharcAnalysis} Error: Material '%s' not recognised.\n",material.c_str());
      break;
    case 'c':
      if    (material.find("targ")==0)
        srim = c_in_targ;
      else if(material.find("si")==0)
        srim = c_in_si;
      else printf("{TSharcAnalysis} Error: Material '%s' not recognised.\n",material.c_str());
      break;
    case 'a':
      if    (material.find("targ")==0)
        srim = a_in_targ;
      else if(material.find("si")==0)
        srim = a_in_si;
      else printf("{TSharcAnalysis} Error: Material '%s' not recognised.\n",material.c_str());
      break;
    default:
      printf("{TSharcAnalysis} Error: Ion '%c' not recognised.\n",ion);
  }

  return srim;
}


double TSharcAnalysis::GetReconstructedEnergy(TVector3 position, int det, double edel, double epad, char ion){

  TSRIM *srim_targ,*srim_si;
  srim_targ= GetSRIM(ion,"targ");
  srim_si  = GetSRIM(ion,"si");
  if(srim_targ==0 || srim_si==0)
    return 0;

  double theta, phi;
  theta = position.Theta();
  phi   = position.Phi();

  double dist_target,dist_deadlayer,dist_pdeadlayer,dist_pad;
  double Ekin,Eafter_target,Eafter_deadlayer,Eafter_del,Eafter_pdeadlayer,Eafter_pad;

  dist_target     = GetTargetThickness(theta,phi); 
  dist_deadlayer  = GetDelDeadLayer(det,theta,phi); 
  dist_pdeadlayer = GetPadDeadLayer(det,theta,phi); 
  dist_pad			  = GetPadThickness(det,theta,phi);

  Eafter_pad        = 0.0; // this is not true when protons punch through DEL+PAD, although it'll be harmless because we can see it distinctly
  Eafter_pdeadlayer = Eafter_pad + epad; // use data
  if(Eafter_pdeadlayer>srim_si->GetEmin())
        Eafter_del  = srim_si->GetEnergy(Eafter_pdeadlayer,-dist_pdeadlayer);
  else  Eafter_del  = 0.0;
  Eafter_deadlayer  = Eafter_del + edel; // use data
  Eafter_target     = srim_si->GetEnergy(Eafter_deadlayer,-dist_deadlayer);
  Ekin              = srim_targ->GetEnergy(Eafter_target,-dist_target);
  
  return Ekin;
}


std::vector<double> TSharcAnalysis::GetMeasuredEnergy(TVector3 position, int det, double ekin, char ion, Option_t *opt, double edel){
  TSRIM *srim_targ,*srim_si;
  std::vector<double> Emeas;
  srim_targ= GetSRIM(ion,"targ");
  srim_si  = GetSRIM(ion,"si");
  if(srim_targ==0 || srim_si==0)
    return Emeas;

	// use options flag to remove target or deadlayer energy loss 
  bool no_target = false, no_deadlayer = false;
  if(strcmp(opt,"no_target")==0)
    no_target = true;
  else if(strcmp(opt,"no_deadlayer")==0)
    no_deadlayer = true;

  double theta, phi;
  theta = position.Theta();
  phi   = position.Phi();

  double dist_target,dist_deadlayer,dist_del,dist_pdeadlayer,dist_pad;
  double Eafter_target,Eafter_deadlayer,Eafter_del,Eafter_pdeadlayer,Eafter_pad;

  // adjust thicknesses of insensitive regions so that their energy loss effects can be switched on and off
  if(no_target)
       dist_target = 0.0;
  else dist_target = GetTargetThickness(theta,phi); 
  
  if(no_deadlayer){
      dist_deadlayer    = 0.0; 
      dist_pdeadlayer   = 0.0;  
  } else {
      dist_deadlayer    = GetDelDeadLayer(det,theta,phi); 
      dist_pdeadlayer   = GetPadDeadLayer(det,theta,phi);  
  }      

  dist_del			  = GetDelThickness(det,theta,phi);
  dist_pad			  = GetPadThickness(det,theta,phi);

  Eafter_target   = srim_targ->GetEnergy(ekin,dist_target);

  if(Eafter_target>srim_targ->GetEmin())
       Eafter_deadlayer = srim_si->GetEnergy(Eafter_target,dist_deadlayer);  
  else Eafter_deadlayer = 0.0;

  if(Eafter_deadlayer>srim_si->GetEmin()){
  	if(edel>0.0)
  		Eafter_del = Eafter_deadlayer - edel;
		else
			Eafter_del = srim_si->GetEnergy(Eafter_deadlayer,dist_del);			
	} else Eafter_del = 0.0;

  if(Eafter_del>srim_si->GetEmin())
       Eafter_pdeadlayer = srim_si->GetEnergy(Eafter_del,dist_pdeadlayer);  
  else Eafter_pdeadlayer = 0.0;

  if(Eafter_pdeadlayer>srim_si->GetEmin())
       Eafter_pad     = srim_si->GetEnergy(Eafter_pdeadlayer,dist_pad);  
  else Eafter_pad     = 0.0;

//  printf(DRED"theta = %.02f\n",theta*180/3.1415);
//  printf(DGREEN"dist_target = %.02f\tdist_deadlayer = %.02f\tdist_dE = %.02f\tdist_pdeadlayer = %.02f\tdist_pad = %.02f\n",dist_target,dist_deadlayer,dist_dE,dist_pdeadlayer,dist_pad);
//  printf(DBLUE"Eafter_target = %.02f\tEafter_deadlayer = %.02f\tEafter_dE = %.02f\tEafter_pdeadlayer = %.02f\tEafter_pad = %.02f\n\n"RESET_COLOR"",Eafter_target,Eafter_deadlayer,Eafter_dE,Eafter_pdeadlayer,Eafter_pad);
  Emeas.push_back( Eafter_deadlayer  - Eafter_del );  // we measure energy lost in dE sensitive region
  Emeas.push_back( Eafter_pdeadlayer - Eafter_pad ); // ... and also in pad sensitive region
  return Emeas;
}


double TSharcAnalysis::GetBeamEnergyInTarget(const char *beam, double ebeam, double frac_depth_end){
  beam_in_targ  = new TSRIM(Form("%s_in_%s.txt",beam,targmat.c_str()),-1,0,false);

  return beam_in_targ->GetEnergy(ebeam,GetTargetThickness(0.0,0.0,frac_depth_end));
}


double TSharcAnalysis::GetResEnergyAfterTarget(const char *res, double eres, double frac_depth_beg){
  TSRIM *res_in_targ  = new TSRIM(Form("%s_in_%s.txt",res,targmat.c_str()),-1,0,false);

  return res_in_targ->GetEnergy(eres,GetTargetThickness(0.0,0.0,1-frac_depth_beg));
}  


TVector3 TSharcAnalysis::GetEdgePoints(int edgenum, int det, int fs, int bs) { // edgenum = 1,2,3,4

	TVector3 pos_lim[4];
	TVector3 position = GetPosition(det,fs,bs);
	double frontpitch = frontpitches[det-1];
	double backpitch  = backpitches[det-1];
	
	if(det>=5 && det <=12) {// BOX
		if(det%2==1){ // 5 7 9 11
			pos_lim[0].SetXYZ(	position.X()-frontpitch/2.0,	position.Y(),	position.Z()-backpitch/2.0	);
			pos_lim[1].SetXYZ(	position.X()-frontpitch/2.0,	position.Y(),	position.Z()+backpitch/2.0	);
			pos_lim[2].SetXYZ(	position.X()+frontpitch/2.0,	position.Y(),	position.Z()-backpitch/2.0	);
			pos_lim[3].SetXYZ(	position.X()+frontpitch/2.0,	position.Y(),	position.Z()+backpitch/2.0	);			
		}
		else if(det%2==0){ // 6 8 10 12
			pos_lim[0].SetXYZ(	position.X(),	position.Y()-frontpitch/2.0,	position.Z()-backpitch/2.0	);
			pos_lim[1].SetXYZ(	position.X(),	position.Y()-frontpitch/2.0,	position.Z()+backpitch/2.0	);
			pos_lim[2].SetXYZ(	position.X(),	position.Y()+frontpitch/2.0,	position.Z()-backpitch/2.0	);
			pos_lim[3].SetXYZ(	position.X(),	position.Y()+frontpitch/2.0,	position.Z()+backpitch/2.0	);	
		}
	}
	else if(det<5 || det>12) {// QQQ
		backpitch*=D2R; // now convert to radians
		double rho = GetPosition(det,fs,bs).Perp();		
		double phi = GetPosition(det,fs,bs).Phi();
				
		pos_lim[0].SetXYZ(	(rho-frontpitch/2)*TMath::Cos(phi-backpitch/2),	(rho-frontpitch/2)*TMath::Sin(phi-backpitch/2),	position.Z()	);
		pos_lim[1].SetXYZ(	(rho-frontpitch/2)*TMath::Cos(phi+backpitch/2),	(rho-frontpitch/2)*TMath::Sin(phi+backpitch/2),	position.Z()	);
		pos_lim[2].SetXYZ(	(rho+frontpitch/2)*TMath::Cos(phi-backpitch/2),	(rho+frontpitch/2)*TMath::Sin(phi-backpitch/2),	position.Z()	);
		pos_lim[3].SetXYZ(	(rho+frontpitch/2)*TMath::Cos(phi+backpitch/2),	(rho+frontpitch/2)*TMath::Sin(phi+backpitch/2),	position.Z()	);			
	}
	return pos_lim[edgenum];
}


double TSharcAnalysis::GetSolidAngle(int det, int fs, int bs){

  double Omega = 0; 
  
  if(det<5 || det>12){
    std::vector<double> rho;
    double rholim[2], z;
    
    for(int i=0;i<4;i++)
      rho.push_back(GetEdgePoints(i,det,fs,bs).Perp());
    
     
    z = GetPosition(det,fs,bs).Z();
    rholim[0] = TMath::MinElement(rho.size(),&rho[0]);
    rholim[1] = TMath::MaxElement(rho.size(),&rho[0]);
    rho.clear();
    
    Omega = fabs(z)*backpitches[det-1]*(TMath::Power(rholim[0]*rholim[0]+z*z,-0.5)-TMath::Power(rholim[1]*rholim[1]+z*z,-0.5));
    
  } else {        
    std::vector<double> x,y,z;
    double xlim[2],ylim[2],zlim[2],perp[2],t[2][2],fixed;        
    
    for(int i=0;i<4;i++){
      x.push_back(GetEdgePoints(i,det,fs,bs).X());  
      y.push_back(GetEdgePoints(i,det,fs,bs).Y());            
      z.push_back(GetEdgePoints(i,det,fs,bs).Z());       
    }
    
    xlim[0] = TMath::MinElement(x.size(),&x[0]);
    xlim[1] = TMath::MaxElement(x.size(),&x[0]);
    ylim[0] = TMath::MinElement(y.size(),&y[0]);
    ylim[1] = TMath::MaxElement(y.size(),&y[0]);          
    zlim[0] = TMath::MinElement(z.size(),&z[0]);
    zlim[1] = TMath::MaxElement(z.size(),&z[0]);       
    
    if(xlim[1]-xlim[0] < 0.001){
     perp[0] = ylim[0];
     perp[1] = ylim[1];
     fixed = xlim[0];
    } else {
     perp[0] = xlim[0];
     perp[1] = xlim[1];
     fixed = ylim[0];
    }
    for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
        t[i][j]=TMath::ATan(perp[i]*zlim[j]/fabs(fixed)*TMath::Power(perp[i]*perp[i]+zlim[j]*zlim[j]+fixed*fixed,-0.5));

  // printf("\nx[] = [%.4f, %.4f]\ty[] = [%.4f, %.4f]\tz[] = [%.4f, %.4f]\nperp[] = [%.4f, %.4f]\nt[] = [%.4f, %.4f; %.4f, %.4f]\n",xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1],perp[0],perp[1],t[0][0],t[0][1],t[1][0],t[1][1]);  
              
    Omega=t[0][0]+t[1][1]-t[1][0]-t[0][1];              
  }
  
  return Omega;  
}


TH2F *TSharcAnalysis::GetSolidAngleMatrix(int det, Bool_t badstrips){

	TH2F *h = new TH2F(Form("SolidAngles_Det%02i",det),Form("Solid Angles for Det%02i; Back Strip; Front Strip",det),48,0,48,24,0,24);

	for(int bs=0; bs<48; bs++){
		if(badstrips && BadStrip(det,-1,bs))
			continue;	
		for(int fs=0; fs<24; fs++){
			if(badstrips && BadStrip(det,bs,-1))
				continue;
			h->SetBinContent(bs,fs,GetSolidAngle(det,fs,bs));
		}
	}
			
	return h;			
}


double TSharcAnalysis::GetSumSolidAngle(int det_min, int fs_min, int bs_min, int det_max, int fs_max, int bs_max){
  
  double Omega = 0; 
  int F_min, F_max, B_min, B_max;
    
  if(det_min<1)
    det_min = 1;
  if(det_max<0)
    det_max = 16;
  else if(det_max<det_min)
    det_max = det_min;
    
  for(int d=det_min; d<=det_max; d++){

    if(fs_min<0)
      F_min = 0;
    if(fs_max<0 || fs_max>=frontstrips[d-1])
      F_max = frontstrips[d-1]-1;
    else if(fs_max<fs_min)
      F_max = F_min;

    for(int f=F_min; f<=F_max; f++){
  
      if(bs_min<0)
        B_min = 1;
      if(bs_max<0 || bs_max>=backstrips[d-1])
        B_max = backstrips[d-1]-1;
      else if(bs_max<bs_min)
        B_max = B_min;
        
      for(int b=B_min; b<=B_max; b++)
          Omega += GetSolidAngle(d,f,b);
   
    }
  }
  
  return Omega;
}


// generates acceptance curve by simulating a set of reactions which are distributed randomly across the target(± position error)
// and land randomly within the pixel area
TH1D *TSharcAnalysis::SetAcceptance(int nmax, TReaction *r, const char *stripsfilename){

	coveragelist = new TList();

	int mmax = 10;
	if(!position_error.Mag()){
		nmax = 1;
		accres = 0;
		printf("\n\t Warning :  Acceptance has no uncertainty when target error is not set.\n\n");
	} else{
		accres = nmax;
	}
	
	if(!r) // reset reaction
		reaction = NULL;
		
	TVector3 pos, tmp_pos;
	SetBadStrips(stripsfilename); 
	double omega, thetalab;	

	TH2F *hcovlab2,*hcovcm2;
	TH1D *hcovlab, *hcovcm, *hcorlab, *hcorcm;
		
	GetRid("SharcCoverageLabVsRun");	
	hcovlab2 = new TH2F("SharcCoverageLabVsRun","SharcCoverageLabMatrix",180,0,180,nmax,0,nmax);
	hcovlab2->SetTitle(Form("%s; Theta Lab [deg]; Simulation Number",hcovlab2->GetTitle()));	
	coveragelist->Add(hcovlab2);
	
	GetRid("SharcCoverageLab");	
	hcovlab = new TH1D("SharcCoverageLab","SharcCoverageLab",180,0,180);
	hcovlab->SetTitle(Form("%s; Theta Lab [deg]; Coverage [sr]",hcovlab->GetTitle()));	
	coveragelist->Add(hcovlab);
	
	GetRid("SharcCorrectionLab");	
	hcorlab = new TH1D("SharcCorrectionLab","SharcCorrectionLab",180,0,180);
	hcorlab->SetTitle(Form("%s; Theta Lab [deg]; Correction factor ",hcorlab->GetTitle()));	
	coveragelist->Add(hcorlab);
	
	std::string rname="";	
	if(r){
		reaction = r;
  	rname.assign(r->GetNameFull());	
	
		hcovcm2 = new TH2F("SharcCoverageCmVsRun","SharcCoverageCmMatrix",180,0,180,nmax,0,nmax);
		hcovcm2->SetTitle(Form("%s; Theta Cm [deg]; Simulation Number",hcovcm2->GetTitle()));	
		coveragelist->Add(hcovcm2);

		hcovcm = new TH1D("SharcCoverageCm",Form("SharcCoverageCm : %s",rname.c_str()),180,0,180);
		hcovcm->SetTitle(Form("%s; Theta Cm [deg]; Coverage [sr]",hcovcm->GetTitle()));	
		coveragelist->Add(hcovcm);

		hcorcm = new TH1D("SharcCorrectionCm",Form("SharcCorrectionCm : %s",rname.c_str()),180,0,180);
		hcorcm->SetTitle(Form("%s; Theta Cm [deg]; Correction factor ",hcorcm->GetTitle()));	
		coveragelist->Add(hcorcm);		
	}
	
	TVector3 targ_pos = GetTargetPosition();

	double x,y,z;		
	printf("\nMaking acceptance curve by randomizing hits..\n");
	for(int n=0; n<nmax; n++){
		printf("\tProgress = %4i/%4i\r",n,nmax);

		tmp_pos.Clear();
		// randomise origin coordinate across uncertainty range
		x = targ_pos.X()+position_error.X()*gRandom->Uniform();
		y = targ_pos.Y()+position_error.Y()*gRandom->Uniform();
		z = targ_pos.Z()+position_error.Z()*gRandom->Uniform();
		// randomly generate a new target position within specified position uncertainty		
		tmp_pos.SetXYZ(x,y,z);
		SetTargetPosition(tmp_pos);

		for(int det=5; det<=16; det++){
			for(int bs=0; bs<GetBackStrips(det); bs++){

				if(BadStrip(det,-1,bs))
					continue; 	
			
				for(int fs=0; fs<GetFrontStrips(det); fs++){
			
					if(BadStrip(det,fs,-1))
						continue;
								
			 		 omega = GetSolidAngle(det,fs,bs)/mmax;
			 		 // now randomize across the pixel to smooth out theta curve
				 	 for(int m=0; m<mmax; m++){
					 	thetalab = RandomizeThetaLab(det,fs,bs);
					 	hcovlab2->Fill(thetalab,n,omega);	
					 	if(r) hcovcm2->Fill(r->ConvertThetaLabToCm(thetalab*D2R,2)*R2D,n,omega);						 		 
					}
				}
			}		
		}
	}
	printf("\tProgress = Complete !!  \n\n");
	
	// restore target position	
	SetTargetPosition(targ_pos);
			
  TF1 *fcovlab, *fcovcm; 
	fcovlab = new TF1("MaxThetaCoverageLab",ThetaCoverage,0,180,2);
  fcovlab->SetParameters(2*PI*D2R,0);

	double theta_val, coverage, variance, cov_err, correction, rel_err;
	TH1D *htmpy;	
	
	// use statistical spread to determine error bars
	for(int i=1;i<=hcovlab->GetNbinsX();i++){

		htmpy = hcovlab2->ProjectionY("",i,i);
		coverage = htmpy->Integral()/((double)nmax); // coverage is average across runs
		if(coverage==0)
			continue;

		theta_val  =  hcovlab->GetBinCenter(i);			
		correction = 	fcovlab->Eval(theta_val)/coverage;
		
		variance =0.0;	
		for(int j=1; j<=htmpy->GetNbinsX(); j++)
			variance  += pow(coverage-htmpy->GetBinContent(j),2.0);	
		variance /=(double)nmax;	
		cov_err = pow(variance,0.5);
		
		rel_err = cov_err/coverage;
		if(correction>10.0 || rel_err>0.9)
			continue;
			
		hcovlab->SetBinContent(i,coverage);
		hcovlab->SetBinError(i,cov_err);
		
		hcorlab->SetBinContent(i,correction);
		hcorlab->SetBinError(i,rel_err * correction);	// same relative error				
	}
	
	if(r){
    fcovcm = new TF1(Form("MaxThetaCoverageCm_%s",rname.c_str()),ThetaCoverage,0,180,2);
	  fcovcm->SetParameters(2*PI*D2R,1);
		coveragelist->Add(fcovcm);
			
		for(int i=1;i<=hcovcm->GetNbinsX();i++){

			htmpy = hcovcm2->ProjectionY("",i,i);
			
			coverage = htmpy->Integral()/((double)nmax); // coverage is average across runs
			if(coverage==0)
				continue;

			theta_val  =  hcovcm->GetBinCenter(i);			
			correction = 	fcovcm->Eval(theta_val)/coverage;
		
			variance =0.0;	
			for(int j=1; j<=htmpy->GetNbinsX(); j++)
				variance  += pow(coverage-htmpy->GetBinContent(j),2.0);	
			variance /=(double)nmax;	
			cov_err = pow(variance,0.5);
		
			rel_err = cov_err/coverage;
			if(correction>5.0 || rel_err>0.9)
				continue;
			
			hcovcm->SetBinContent(i,coverage);
			hcovcm->SetBinError(i,cov_err);
		
			hcorcm->SetBinContent(i,correction);
			hcorcm->SetBinError(i,rel_err * correction);	// same relative error				
		}
	}
		
	TH1D *h = hcorlab;
		
	int cmint = 0;
	if(r){
		h = hcorcm;
		cmint = 1;
	}

	gStyle->SetOptStat(0);
	GetRid("SharcAcceptanceCurves");	
	TCanvas *c = new TCanvas("SharcAcceptanceCurves","SharcAcceptanceCurves");
	TLegend *leglab = new TLegend(0.7,0.7,0.9,0.9);
	leglab->AddEntry(fcovlab,"Full phi coverage","lp");
	leglab->AddEntry(hcovlab,"Sharc coverage","lp");	
		
	c->Divide(1+cmint,2);
	c->cd(1);
	hcovlab->GetYaxis()->SetRangeUser(0,fcovlab->GetMaximum());
	hcovlab->Draw();
	fcovlab->Draw("same");
	leglab->Draw();
	c->cd(2+cmint);
	gPad->SetGridy(1);
	hcorlab->Draw();		
	
	if(r){
		TLegend *legcm = new TLegend(0.7,0.7,0.9,0.9);
		legcm->AddEntry(fcovcm,"Full phi coverage","lp");
		legcm->AddEntry(hcovcm,"Sharc coverage","lp");		

		c->cd(2);
		hcovcm->GetYaxis()->SetRangeUser(0,fcovcm->GetMaximum());		
		hcovcm->Draw();
		fcovcm->Draw("same");				
		legcm->Draw();
		c->cd(4);
		gPad->SetGridy(1);		
		hcorcm->Draw();
	}
		
	coveragelist->Add(c);
			
	/*	
	TCanvas *c = new TCanvas;
	
	c->Divide(1,2);
	c->cd(1);
	hcovlab->Draw();
	c->cd(2);
	hcovlab2->Draw("colz");
		*/
		
	return h;
}


TH1D *TSharcAnalysis::RebinAcceptance(TH1D *h, int binsz){

	TH1D *hreb = new TH1D(Form("%s_%iDegBins",h->GetName(),binsz),"",180/(double)binsz,0,180);
	hreb->SetTitle(Form("%s;%s;%s",h->GetTitle(),h->GetXaxis()->GetTitle(),h->GetYaxis()->GetTitle()));
	
	Double_t val, err, sumval, sumerr;
	Int_t nval = 0, n = 0;
	
	printf("\n\tRebinning acceptance curve [%i deg. bins]\n",binsz);
	
	for(int i=1; i<=180; i++){
	
		// only average over non-zero bins
		val = h->GetBinContent(i);		
		err = h->GetBinError(i);		
	//  printf("\n\t theta = %i\t val = %.2e +/- %.2e",i,val,err);
		if(val){// && err){
			sumval += val;
			sumerr += err;
			nval++; 
		}
				
		if(i%binsz==0){
			n++;
			if(nval){
		  	printf("\n\t Bin %i. Theta = %.1f\t   content = %.3f +/- %.3f  [npts = %i]",n,(double)i-binsz*0.5,sumval/nval,sumerr/nval,nval);
				hreb->SetBinContent(n,sumval/nval);
				hreb->SetBinError(n,sumerr/nval);				
			}
			nval = 0;
			sumval = 0;
			sumerr = 0;		
		}
	}
	printf("\n\n");
	
	return hreb;			
}


TH1D *TSharcAnalysis::MakeSin(Double_t binsz){

	TH1D *hsin = new TH1D("sin","sin",180,0,180);
	double integral=0.0;
  for(int i=0; i<180; i++){
		hsin->SetBinContent(i,TMath::Sin(D2R*hsin->GetBinCenter(i)));
		hsin->SetBinError(i,0.0);
    integral+=hsin->GetBinContent(i)*hsin->GetBinWidth(i)*D2R;
  }
  hsin->Scale(2./integral);	// integral must be 2.
	
	hsin->Rebin(binsz);				// rebin theta
	hsin->Scale(1/binsz);	
	
	return hsin;	
}


TH1D *TSharcAnalysis::SimulateAngDist(const char *sname, double &cntres, double &reserr, double t1, double t2, double counts, double counterr){

	static TH1D *hacc;
	if(!hacc){
		const char *stripsfile = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";
		TReaction *r = new TReaction("sr95","d","p","sr96",510.9,0,true);
		SetTarget(0.0,0.0,0.0,4.5,"cd2",0.5,0.5,0.5);
	
		TList *acclist = TSharcAnalysis::GetAcceptanceList(r,stripsfile,30);	  
		hacc = (TH1D*)acclist->FindObject("SharcCorrectionCm");	
		Double_t er, ct = hacc->IntegralAndError(1,180,er);
		
		printf("\n\n\nAcceptance correction integral = %.3e +/- %.3e [rel err = %.2f %%]\n\n",ct,er,er/ct*100.0);
	}
	GetRid("sin");
	TH1D *hsin = MakeSin(1.0);
	
	std::string str;
	str.assign(sname);
	std::ifstream infile;
	if(str.find(".txt")!=std::string::npos)
		infile.open(str.c_str());
	else
		infile.open(Form("%s/Fresco%s.txt",SHCDIR,sname));
  TH1D *hsigma = new TH1D(sname,Form("%s; Theta Cm [deg.]; Differential Cross Section [mb/sr]",sname),180,0,180);
  hsigma->SetLineColor(kRed);
  hsigma->Sumw2();
  Double_t thetacm, dsig;
  while(infile.good()){
    infile >> thetacm >> dsig;
   // printf("\n\tthetacm = %.2f\tdsig = %.2e",thetacm,dsig);
    hsigma->SetBinContent((int)thetacm,dsig);
    hsigma->SetBinError((int)thetacm,0.0);
  }
  infile.close();
  
	GetRid("SimulatedCounts");  
  TCanvas *canv = new TCanvas("SimulatedCounts","Counts, Simulated using DWBA",1200,600);
	canv->Divide(2,1);
	canv->cd(1);
	hsigma->Draw();
		
	Double_t sigtot = hsigma->Integral();  
  if(sigtot==0){
  	printf("\n\t Warning :  No FRESCO data was found for %s.\n\n",str.c_str());
  	return 0;
	}

	GetRid("SineCorrectedSigma");
  TH1D *hsincor = (TH1D*)hsigma->Clone("SineCorrectedSigma");
  hsincor->SetTitle("Sine Corrected Cross Section; Theta Cm [deg.]; Cross Section / Sine");
  hsincor->SetLineColor(kGreen);  
  hsincor->Multiply(hsin);  // sine correction  
  hsincor->Draw("same");

	GetRid("SineAccCorrectedSigma");    
  TH1D *hacccor = (TH1D*)hsincor->Clone("SineAccCorrectedSigma");
  hacccor->SetTitle("Acceptance & Sine Corrected Cross Section; Theta Cm [deg.]; Cross Section / Sine * Acceptance ");  
  hacccor->SetMarkerColor(kBlue);
  hacccor->SetLineColor(kBlue);
  hacccor->Divide(hacc);    // reduce counts at each angle to reflect dead strips 
  hacccor->Draw("same e1");
  
  TLegend *leg = new TLegend(0.5,0.6,1.0,0.9);
  leg->AddEntry(hsigma,"FRESCO Cross Seciton","lp");
  leg->AddEntry(hsincor,"FRESCO * Sin","lp");  
  leg->AddEntry(hacccor,"FRESCO * Sin * Acc","lp");  
  leg->SetLineColor(kBlack);  
  leg->Draw("same");
  gPad->SetLogz();
  
  TH1D *hcounts = (TH1D*)hacccor->Clone(Form("SimulatedCounts_%s",sname));
  hcounts->SetTitle("Simulated Counts; Theta Cm [deg.]; Counts ");  
  hcounts->SetLineColor(kBlue); 
  hcounts->SetMarkerColor(kBlue);


  canv->cd(2);   
  ////////////////////  SF DETERMINED USING INTEGRAL OVER ALL ANGLES   //////////////////// 
  // only one angular dependent uncertainty, comes form acceptance
  // we need to use IntegralAndError(t1,t2,err) to get the error as it depends on t1,t2
  // then we apply global scaling for normalization and SF. 
  // The errors on norm and SF these are not independent, so are added to err in quadrature
  // err_tot = sqrt( pow(err/tot,2.0) + pow(normerr/norm,2.0) + pow(sferr/sf,2.0) )
  
	GetRid("Normalization");
	// scalar hists have no error : it is applied to integral error at the end		
//	Double_t norm = 7.81e-4, normerr = 0.06e-4;
	Double_t norm = 8.8e-4, normerr = 0.1e-4;
 	TH1D *hnorm = MakeScaleHist(1,norm,0.0,"Normalization"); 
  hcounts->Divide(hnorm);   // convert mb/sr into counts  
  
	hcounts->Draw("e1");

	// predicted counts for SF=1
	Double_t cntmaxerr, cntmax = hcounts->IntegralAndError(1,180,cntmaxerr); 
	// now calculate norm error	
	cntmaxerr = cntmax*sqrt(pow(cntmaxerr/cntmax,2.0)+pow(normerr/norm,2.0)); 
	
	TH1D *hsfcounts = (TH1D*)hcounts->Clone(Form("FixedCounts_%s",sname));
	hsfcounts->SetLineColor(kBlack);
  hsfcounts->SetMarkerColor(kBlack);
	
	GetRid("SpecFactor");		
	// determine SF required to reproduce input total counts
	Double_t sf = counts/cntmax; 
	// now calculate SF error	
	Double_t sferr = sf*sqrt(pow(cntmaxerr/cntmax,2)+pow(counterr/counts,2)); // error
		// scalar hists have no error : it is applied to integral error at the end	
 	TH1D *hsf = MakeScaleHist(1,sf,0.0,"SpecFactor"); 
 	hsfcounts->Multiply(hsf);
 	
	hsfcounts->Draw("same e1");
	// predicted counts for SF!=1	
	Double_t cntset = hsfcounts->Integral(1,180); 
	// add experimental relative error in quadrature	
	Double_t cntseterr = cntset*sqrt(pow(cntmaxerr/cntmax,2.0)+pow(counterr/counts,2)); 
	
 // predicted counts for SF!=1 over theta range, t1-t2
	cntres = hsfcounts->IntegralAndError((int)t1,(int)t2,reserr); 	
 // now apply all errors in quadrature
//	printf("\n\n reserr/cntres = %.3f    cntseterr/cntset = %.3f",reserr/cntres,cntseterr/cntset);
	reserr = cntres*sqrt(pow(reserr/cntres,2.0)+pow(cntseterr/cntset,2.0));
//	printf("  ---> reserr/cntres = %3f\n",reserr/cntres);

	// PHEW !
	////////////////////////////////////////////////////////////////////////////////////////	
	
	
	TLegend *leg2 = new TLegend(0.4,0.6,0.95,0.85);
	leg2->AddEntry(hcounts,Form("Total counts [SF=1.0] = %.2e +/- %.2e",cntmax,cntmaxerr),"lp");
	leg2->AddEntry(hsfcounts,Form("Counts [SF=%.3f] = %.2e +/- %.2e",sf,cntset,cntseterr),"lp");

	leg2->Draw();
  	
  printf("\n\n Simulation Results for %s, Measured counts = %.2e +/- %.2e :-",sname,counts,counterr);
  printf("\n\t -> Total Cross Section             = %.2e mb/sr",sigtot);
  printf("\n\t -> Total Proton Counts [SF=1]      = %.2e +/- %.2e",cntmax,cntmaxerr); 
  printf("\n\t -> Spectroscopic Factor            = %.3f +/- %.3f",sf,sferr);  
  printf("\n\t -> Total Proton Counts [SF=%.2f]   = %.2e +/- %.2e",sf,cntset,cntseterr); 
  printf("\n\t -> Counts In Rng [%.1f-%.1f]       = %.2e +/- %.2e\n\n",t1,t2,cntres,reserr);  
  	
  return hsfcounts;
}


TList *TSharcAnalysis::GetAcceptanceList(TReaction *r, const char *stripsfile, Int_t resolution){

	Bool_t same_reac=false, same_file=false, same_res=false;
	if(r && reaction && strcmp(r->GetName(),reaction->GetName())==0){
		same_reac = true;
	} else if(!r && !reaction){
		same_reac = true;	
	}
	
	if(strcmp(stripsfile,badstripsfile.c_str())==0){
		same_file = true;
	}
	
	if(!position_error.Mag() || resolution == accres){
		same_res = true;
	}
	
	if(coveragelist && same_reac && same_file && same_res){
		printf("\n\t Using existing acceptance results.\n\n");
		return coveragelist;
	}

	SetAcceptance(resolution,r,stripsfile);
	
	return coveragelist;	
}


void TSharcAnalysis::SetBadStrips(const char *stripsfilename){
	badstripsfile.assign(stripsfilename);
	resetbadstrips = true;
	coveragelist = new TList();
	BadStrip();
}


int TSharcAnalysis::BadStrip(int det, int fs, int bs){
 	
 	static unsigned int badfrontstrip[16][24],badbackstrip[16][48];
	static TH2F *h[12];
 	
 	if(resetbadstrips){
 	
 		nbadstrips=0;

		for(int i=0; i<15; i++){
			for(int j=0; j<24; j++){
				badfrontstrip[i][j]   = 0;
				badbackstrip[i][j]    = 0;		
				badbackstrip[i][24+j] = 0;		
			}
		}		
		
		ifstream infile;
    infile.open(badstripsfile.c_str());
    printf("\nTSharcAnalysis :  Bad Strips File '%s' will be used.. ",badstripsfile.c_str());
    if(infile.is_open()){
    	printf(" Successfully opened!"); 	
			int d,f,b;	
			
			while(infile.good()){
				infile >> d >> f >> b;
	
				nbadstrips++;
				if(d>0 && d<=16){
					if(f>-1) badfrontstrip[d-1][f] = 1;
					else if(b>-1) badbackstrip[d-1][b] = 1;
				}else nbadstrips--;
			}
	 
   	}	else { 
			printf("File openening FAILED"); 
		}
		
	 	coveragelist = new TList();
 		int maxstrips = 0;		
		for(int dd=5; dd<=16; dd++){
			maxstrips+=GetFrontStrips(dd)+GetBackStrips(dd);
			GetRid(Form("SharcHitPattern_Det%02i",dd));
			h[dd-5] = new TH2F(Form("SharcHitPattern_Det%02i",dd),Form("SharcHitPattern_Det%02i; Back Strip; Front Strip",dd),48,0,48,24,0,24);
			coveragelist->Add(h[dd-5]);
						
			for(int ff=0; ff<GetFrontStrips(dd); ff++){
				if(badfrontstrip[dd-1][ff])
					continue;
				
				for(int bb=0; bb<GetBackStrips(dd); bb++){
					if(badbackstrip[dd-1][bb])
						continue;					
				
					h[dd-5]->Fill(bb,ff,1); 
				}
			}
		}			
		infile.close();	
 		resetbadstrips = false;	
	
 		printf("\n\t-> Made bad strips list. %i / %i strips [ %4.2f %% ] will be excluded.\n\n",nbadstrips,maxstrips,((Double_t)nbadstrips/(Double_t)maxstrips)*100.0);
	}
	
	if(det<5 || (fs<0 && bs<0)){
		return -1; // invalid argument
  } else if(fs>-1 && fs<GetFrontStrips(det)){
  	h[det-5]->Fill(bs,fs);
  	return badfrontstrip[det-1][fs];
	} else if(bs>-1 && bs<GetBackStrips(det)){
		h[det-5]->Fill(bs,fs);
    return badbackstrip[det-1][bs];
  } else {			 
  	return -1;
  }	
}


double TSharcAnalysis::RandomizeThetaLab(int det, int fs, int bs){
   // radnomizes hits across pixel centers. Assumes continuum limit, so grassiness will still be present when there are few pixels in a given theta range.
	 // assumes flat weighting, and thus solid angle, between theta limits in pixel... Only accurate for small pixels
   double theta[4], theta_max, theta_min;
   theta_max=0.0;
   theta_min=180.0;
   
   for(int i=0;i<4;i++){ 
    theta[i]= TSharcAnalysis::GetEdgePoints(i,det,fs,bs).Theta();  
    if(theta[i]>theta_max)
       theta_max=theta[i];
    if(theta[i]<theta_min)
       theta_min=theta[i];
   }
  
  return gRandom->Uniform(theta_min,theta_max)*R2D;
}


double TSharcAnalysis::ThetaCoverage(Double_t *x, Double_t *pars){
	
	Double_t val = pars[0]*TMath::Sin(D2R*x[0]);
  if(pars[1]==1){ // CM frame = 1
     static TGraph *g;
     if(!g){
        g = reaction->OmegaVsTheta(0,180,2,false);    
        coveragelist->Add(g);
     }
     val*=g->Eval(x[0]);
  }
  return val;
}


void TSharcAnalysis::SetOmegaMaxMin(bool printout){

  double omega_temp = 0.0, omega_max = 0.0000, omega_min = 1000.0;
  for(int d=1;d<=16;d++){
		for(int fs=0; fs<frontstrips[d-1]; fs++){	
			for(int bs=0; bs<backstrips[d-1]; bs++){
			
				omega_temp = GetSolidAngle(d,fs,bs);
				if(omega_temp > omega_max)
					omega_max = omega_temp;
				else if(omega_temp < omega_min)
					omega_min = omega_temp;
			}
		}
  }

  Omega_max = omega_max;
  Omega_min = omega_min;

  if(printout)  printf("{TSharcAalysis}:\t Omega_max = %f\tOmega_min = %f\n",Omega_max,Omega_min);

  return;
}


TH1D *TSharcAnalysis::MakeScaleHist(Double_t binsz, Double_t scale, Double_t err, const char *name){

  Int_t nbins = (Int_t)180.0/binsz;
  TH1D *hscale = new TH1D(name,name,nbins,0,180);
  hscale->SetTitle(Form("Scale Histogram for %s; Theta [deg]; Scale",name));
  for(int i=1; i<=nbins; i++){
    hscale->SetBinContent(i,scale);
    hscale->SetBinError(i,err);
  }
  return hscale;
}    


void TSharcAnalysis::Print(Option_t *opt) {

	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
	printf("\n\n\t____TSharcAnalysis____");
	
	if(strcmp(opt,"all")==0){
		printf("\n\n\tSHARC configuration table :-");
		printf("\n\n DET : FS / pitch : BS / pitch : DEL[um] : DEL[dead][um] : PAD[um] : PAD[dead][um] \n");
		for(int det=1; det<=16; det++){
		
			printf("\n%4i %4i%5.1f[mm]%4i%5.1f%s%9.2f%10.2f",det,frontstrips[det-1],frontpitches[det-1],
				backstrips[det-1],backpitches[det-1],det<4||det>12?"[rad]":"[mm] ",DELthicknesses[det-1],DELdeadlayers[det-1]);
			if(det>4 && det<=8)
				printf("%15.2f%10.2f",PADthicknesses[det-1],PADdeadlayers[det-1]);
		}
	}
	printf("\n\n\n\t Target  :-\n\n Material = ' %s '\n Thickness =  %.2f [um]\n Position =  [%.2f,%.2f,%.2f] [mm]",
			targmat.c_str(),targetthickness,position_offset.X(),position_offset.Y(),position_offset.Z());
	
	if(position_error.Mag())
		printf(" +/- [%.2f,%.2f,%.2f] [mm]",position_error.X(),position_error.Y(),position_error.Z());

	if(badstripsfile.length())
		printf("\n\n\t BadStrips   :-\n\n File Name = ' %s '\n Number Of Bad Strips = %i",badstripsfile.c_str(),nbadstrips);
	if(reaction){
		printf("\n\n\t Reaction   :-\n\n  Name = ' %s '",reaction->GetNameFull());
		if(strcmp(opt,"all")==0)
			reaction->Print("all");
	}
	printf("\n\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
		
}


void TSharcAnalysis::Clear(Option_t *opt) {

// we don't want this to run if you don't need it
	p_in_targ = 0;
	p_in_si   = 0;
	d_in_targ = 0;
	d_in_si   = 0;
	c_in_targ = 0;
	c_in_si   = 0;
	a_in_targ = 0;
	a_in_si   = 0;
	beam_in_targ  = 0;

	position_offset.SetXYZ(0,0,0);
	position_error.SetXYZ(0,0,0);
	
	ResetDetectors();

// target properties (um)
	targetthickness    = 0.0;    
	targetradius       = 2000; 
	targmat					   = "";
	badstripsfile			 = "";
	resetbadstrips     = true;

	Omega_min = 0.0000;
	Omega_max = 0.0000;

}


void TSharcAnalysis::ResetDetectors(){

// segmentation and pitch sizes (warning, QQQ back pitch is angular)
	number_of_detectors= 16;	
	for(int i=0; i<4; i++){
	
		frontstrips[i]  = 16;		 frontstrips[4+i]  = 24;  frontstrips[8+i]	= 24;  frontstrips[12+i]= 16;
		backstrips[i]		= 24;	 	 backstrips[4+i]	 = 48;  backstrips[8+i]	  = 48;	 backstrips[12+i]= 24; 
		frontpitches[i]	=	2.0;	 frontpitches[4+i] = 3.0; frontpitches[8+i] =	3.0; frontpitches[12+i]=	2.0;
		backpitches[i]	=	PI/48; backpitches[4+i]  = 1.0; backpitches[8+i]  =	1.0; backpitches[12+i]=	PI/48;

// various thicknesses in um
		DELdeadlayers[i] = 0.7; DELdeadlayers[4+i] = 0.1; DELdeadlayers[8+i] = 0.1; DELdeadlayers[12+i] = 0.7; 
		PADdeadlayers[i] = 0.0; PADdeadlayers[4+i] = 1.0; PADdeadlayers[8+i] = 0.0; PADdeadlayers[12+i] = 0.0; 
		
		PADthicknesses[i] = 0.0; PADthicknesses[8+i] = 0.0; PADthicknesses[12+i] = 0.0;
	}
	
	DELthicknesses[0]  = 998.0; DELthicknesses[1] = 998.0;  DELthicknesses[2] = 998.0;   DELthicknesses[3] = 1001.0;
	DELthicknesses[4]  = 141.0; DELthicknesses[5] = 142.0;  DELthicknesses[6] = 133.0;   DELthicknesses[7] = 143.0;
	DELthicknesses[8]  = 999.0; DELthicknesses[9] = 1001.0; DELthicknesses[10] = 1001.0; DELthicknesses[11] = 1002.0;
	DELthicknesses[12] = 390.0; DELthicknesses[13] = 390.0; DELthicknesses[14] = 383.0;  DELthicknesses[15] = 385.0;
	
	PADthicknesses[4]  = 1534.0; PADthicknesses[5] = 1535.0;  DELthicknesses[6] = 1535.0;   DELthicknesses[7] = 1539.0;
}


void TSharcAnalysis::GetRid(const char *name, Bool_t delete_all){
	TObject *obj;

	if(coveragelist){
		obj = coveragelist->FindObject(name);
		if(obj) coveragelist->Remove(obj);
	}
	
	if(delete_all){
		obj = gROOT->FindObjectAny(name);
		if(!obj) return;
		if(strcmp(obj->ClassName(),"TCanvas")==0){
			TCanvas *canv = (TCanvas*)obj;
			canv->Clear();
			canv->Closed();
	//		canv = NULL;
			gSystem->ProcessEvents();
			delete canv;
		}
		else obj->Delete();
	}
	
}


