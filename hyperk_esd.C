// Based on alice_esd.C from ROOT tutorials
// Modified for hyperk by Alex Finch


WCSimRootGeom * wcsimrootgeom;
WCSimRootEvent * wcsimrootEvent;
WCSimRootTrigger * wcsimrootTrigger;
TFile* f;
TTree * wcsimGeoT;
TTree * wcsimT;
//TH1F* fHistogram;

void       make_gui();
void       load_event();
void       createGeometry();
void       hyperk_load_event();


Int_t event_id       = 0; // Current event id.

TEveTrackList *gTrackList = 0;

TEveGeoTopNode* geomRoot;


/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/

//______________________________________________________________________________
void hyperk_esd()
{
	// Main function, initializes the application.
	
	
	const TString weh("alice_esd()");
	
	const char *filetypes[] = {
		"ROOT files",    "*.root",    	
		"All files",     "*",
		0,               0
	};
	
	TString CurrentDirectory=gSystem->pwd();
	TString originalDirectory=CurrentDirectory;
	TGFileInfo fi;
	fi.fFileTypes = filetypes;
	fi.fIniDir    = StrDup(CurrentDirectory);
	new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, &fi);
	if (!fi.fFilename) {
		cout<<" No file chosen "<<endl;
		return;
	}
	f = new TFile(fi.fFilename);
	gSystem->cd(originalDirectory);
	
	
    	wcsimGeoT=(TTree *) f->Get("wcsimGeoT");
    	wcsimrootgeom = new WCSimRootGeom (); 
    	wcsimGeoT -> SetBranchAddress ("wcsimrootgeom" ,&wcsimrootgeom );
    	
    	wcsimGeoT->GetEntry(0) ;
    	
    	TEveManager::Create(kTRUE,"V");
    	TEveManager::Create();
    	gROOT->LoadMacro("hyperk_esd_html_summary.C");
    	fgHtmlSummary = new HtmlSummary("HyperK Event Display Summary Table");
    	slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    	fgHtml = new TGHtml(0, 100, 100);
    	TEveWindowFrame *wf = slot->MakeFrame(fgHtml);
    	fgHtml->MapSubwindows();
    	wf->SetElementName("Summary");
    	
    	
    	gSystem->Load("libGeom");
    	new TGeoManager("HyperK", "HyperK Detector");
    	createGeometry();
    	gGeoManager->SetTopVisible();
    	TGeoNode* node = gGeoManager->GetTopNode();
    	geomRoot = new TEveGeoTopNode(gGeoManager, node);
    	geomRoot->SetVisLevel(4);
    	geomRoot->GetNode()->GetVolume()->SetVisibility(kFALSE);
    	gEve->AddGlobalElement(geomRoot);   	

    	
    	gEve->GetBrowser()->GetTabRight()->SetTab(1);
    	
    	make_gui();
    	
    	wcsimT = (TTree*) f->Get("wcsimT");
    	wcsimrootEvent = new WCSimRootEvent ();
    	wcsimT -> SetBranchAddress ("wcsimrootevent" ,&wcsimrootEvent);   	
    	wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
    	wcsimT->GetEvent(0);
    	cout<<" there are "<<wcsimT->GetEntries()<<" events "<<endl;
    	
    	//TCanvas*  Hists = new TCanvas("Hists","Histograms",-1,300,400,200);
    	
    	
    	load_event();
    	
    	
    	gEve->GetDefaultGLViewer()->UpdateScene();
    	gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.
    	
}

//______________________________________________________________________________

//______________________________________________________________________________
void load_event()
{
	// Load event specified in global event_id.
	// The contents of previous event are removed.
	
	printf("Loading event %d.\n", event_id);
	
	gEve->GetViewers()->DeleteAnnotations();
	TEveEventManager* CurrentEvent =gEve->GetCurrentEvent();
	if(CurrentEvent != 0)CurrentEvent->DestroyElements();	
	
	wcsimT->GetEvent(event_id);	
	//	cout<<" this event has "<<wcsimrootEvent->GetNumberOfEvents()<<" triggers "<<endl;
	wcsimrootTrigger = wcsimrootEvent->GetTrigger(0); 
	//    	cout<<"trigger 0: "<< vertex "<<wcsimrootTrigger->GetVtx(0)<<" "<<wcsimrootTrigger->GetVtx(1)<<" "<<wcsimrootTrigger->GetVtx(2)<<endl;
	hyperk_load_event();
	TEveElement* top = gEve->GetCurrentEvent();
	//	update_html_summary();
	gEve->Redraw3D(kFALSE, kTRUE);
}
/******************************************************************************/
// GUI
/******************************************************************************/
//______________________________________________________________________________
// 
// EvNavHandler class is needed to connect GUI signals.

class EvNavHandler
{
public:
	void Fwd()
	{
		if (event_id < wcsimT->GetEntries()) {
			++event_id;
			load_event();
		} else {
			printf("Already at last event.\n");
		}
	}
	void Bck()
	{
		if (event_id > 0) {
			--event_id;
			load_event();
		} else {
			printf("Already at first event.\n");
		}
	}
};
//______________________________________________________________________________
void make_gui()
{
	// Create minimal GUI for event navigation.
	
	TEveBrowser* browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kLeft);
	
	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
	frmMain->SetWindowName("XX GUI");
	frmMain->SetCleanup(kDeepCleanup);
	
	TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
	{
		
		TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
		TGPictureButton* b = 0;
		EvNavHandler    *fh = new EvNavHandler;
		
		b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
		hf->AddFrame(b);
		b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");
		
		b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
		hf->AddFrame(b);
		b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
	}
	frmMain->AddFrame(hf);
	
	frmMain->MapSubwindows();
	frmMain->Resize();
	frmMain->MapWindow();
	
	browser->StopEmbedding();
	browser->SetTabTitle("Event Control", 0);
}
//______________________________________________________________________________
void hyperk_load_event()
{
	//fHistogram->Reset();
	/*
	Loop over all cherenkov hits and add them to Eve event
	*/
	gStyle->SetPalette(1, 0);
	
	TEveRGBAPalette* pal = new TEveRGBAPalette(0, 10);
	
	
	TEveBoxSet*  CherenkovHits= new TEveBoxSet("CherenkovHits");
	CherenkovHits->SetPalette(pal);
	//    CherenkovHits->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
	CherenkovHits->Reset(TEveBoxSet::kBT_Cone, kFALSE, 64);
	
	float PMTRadius = wcsimrootgeom->GetWCPMTRadius() ;
	int max=wcsimrootTrigger->GetNcherenkovdigihits();
	for (int i = 0; i<max; i++){
		WCSimRootCherenkovDigiHit *cDigiHit = (WCSimRootCherenkovDigiHit*) wcsimrootTrigger->GetCherenkovDigiHits()->At(i);
		int tubeId=cDigiHit->GetTubeId();
		float Q=cDigiHit->GetQ(); 
		WCSimRootPMT pmt = wcsimrootgeom -> GetPMT ( tubeId );
		double pmtX = pmt.GetPosition (0);
		double pmtY = pmt.GetPosition (1);
		double pmtZ = pmt.GetPosition (2);
		int location=pmt.GetCylLoc();
		if(location==1)
		{
			float lengthXY=sqrt(pmtX*pmtX+pmtY*pmtY);
			float xd=pmtX/lengthXY;
			float yd=pmtY/lengthXY;
			CherenkovHits->AddCone(TEveVector(pmtX,pmtY,pmtZ),TEveVector(xd,yd,0.0) ,PMTRadius);
		}
		else
		{
			CherenkovHits->AddCone(TEveVector(pmtX,pmtY,pmtZ),TEveVector(0.0,0.0,1.0) ,PMTRadius ); 
		}
		CherenkovHits->DigitValue(Q);
	}
	CherenkovHits->RefitPlex();
	
	TEveTrans& t = CherenkovHits->RefMainTrans();
	t.SetPos(0.0,0.0,0.0);
	gEve->AddElement(CherenkovHits);
	
	/*
	Loop over truth tracks and add them to Eve
	*/
	
	// Get the number of tracks
	int ntrack = wcsimrootTrigger->GetNtrack();
	printf("ntracks=%d\n",ntrack);
	
	int npar = wcsimrootTrigger->GetNpar();
	
	
	int i;
	// Loop through elements in the TClonesArray of WCSimTracks
	for (i=0; i<ntrack; i++)
	{
		TObject *element = (wcsimrootTrigger->GetTracks())->At(i);
		
		WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
		int pdgCode=wcsimroottrack->GetIpnu();
		TString Name;
		if(pdgCode==11)Name="electron";
		if(pdgCode==-11)Name="positron";
		if(pdgCode==12)Name="electron neutrino";
		if(pdgCode==-12)Name="electron anti-neutrino";
		if(pdgCode==13)Name="muon";
		if(pdgCode==-13)Name="anti-muon";
		if(pdgCode==14)Name="muon neutrino";
		if(pdgCode==-14)Name="muon anti-neutrino";
		if(pdgCode==2212)Name="proton";
		if(pdgCode==-2212)Name="anti-proton";
		TEveLine* track=new TEveLine(Name);
		track->SetElementTitle(Name);
		if(abs(pdgCode)==11)track->SetMainColor(kYellow);
		if(abs(pdgCode)==12)track->SetMainColor(kBlue);
		if(abs(pdgCode)==13)track->SetMainColor(kMagenta);
		if(abs(pdgCode)==14)track->SetMainColor(kGreen);
		if(abs(pdgCode)==2212)track->SetMainColor(kRed);
		track->SetNextPoint(wcsimroottrack->GetStart(0),wcsimroottrack->GetStart(1),wcsimroottrack->GetStart(2));
		track->SetNextPoint(wcsimroottrack->GetStop(0),wcsimroottrack->GetStop(1),wcsimroottrack->GetStop(2));
		gEve->AddElement(track);
		
		
		
	}  // End of loop over tracks
	//Hists->Update();
	update_html_summary(wcsimrootTrigger);
	
	
}
void createGeometry()
{
    	float maxX=0;
    	float maxY=0;
    	float maxZ=0;
	for(int tubeId=0;tubeId<wcsimrootgeom->GetWCNumPMT();tubeId++)
    	{
    		WCSimRootPMT pmt = wcsimrootgeom -> GetPMT ( tubeId );
    		double pmtX = pmt.GetPosition (0);
    		double pmtY = pmt.GetPosition (1);
    		double pmtZ = pmt.GetPosition (2);
    		if(pmtX>maxX)maxX=pmtX;
    		if(pmtY>maxY)maxY=pmtY;
    		if(pmtZ>maxZ)maxZ=pmtZ;   		
    	}
	
	// MATERIALS, MIXTURES AND TRACKING MEDIA
	// Material: world
	Double_t a       = 0.000000;
	Double_t z       = 0.000000;
	Double_t density = 0.000000;
	Double_t radl    = 1000000000000000019884624838656.000000;
	Double_t absl    = 1000000000000000019884624838656.000000;
	TGeoMaterial* pMat1 = new TGeoMaterial("world", a,z,density,radl,absl);
	pMat1->SetIndex(0);
	pMat1->SetTransparency(60);
	// Medium: medium0
	Int_t numed   = 0;  // medium number
        Double_t par[20];
	par[0]  = 0.000000; // isvol
	par[1]  = 0.000000; // ifield
	par[2]  = 0.000000; // fieldm
	par[3]  = 0.000000; // tmaxfd
	par[4]  = 0.000000; // stemax
	par[5]  = 0.000000; // deemax
	par[6]  = 0.000000; // epsil
	par[7]  = 0.000000; // stmin
	pMed1 = new TGeoMedium("medium0", numed,pMat1, par);
	
	
	dx = 2*maxX;
	dy = 2*maxY;
	dz = 2*maxZ;
	TGeoShape *pworldbox_1 = new TGeoBBox("worldbox", dx,dy,dz);
	// Volume: volume0
	WorldVolume = new TGeoVolume("volume0",pworldbox_1, pMed1);
	WorldVolume->SetVisLeaves(kTRUE);
	
	// SET TOP VOLUME OF GEOMETRY
	gGeoManager->SetMaxVisNodes(wcsimrootgeom->GetWCNumPMT()+1000);
	gGeoManager->SetTopVolume(WorldVolume);
	
	// SHAPES, VOLUMES AND GEOMETRICAL HIERARCHY
	// Shape: phototube type: TGeoSphere
	float PMTRadius = wcsimrootgeom->GetWCPMTRadius() ;
	cout<<" Phototube radius is "<<PMTRadius<<endl;
	Double_t rmin   = 0.900000*PMTRadius;
	Double_t rmax   = PMTRadius;
	Double_t theta1 = 0.000000;
	Double_t theta2 = 90.000000;
	Double_t phi1   = 0.000000;
	Double_t phi2   = 360.000000;
	TGeoShape *PhotoTubeShape = new TGeoSphere("phototube",rmin,rmax,theta1, theta2,phi1,phi2);
	// Volume: phototube
	PhotoTubeVolume = new TGeoVolume("phototube",PhotoTubeShape, pMed1);
	PhotoTubeVolume->SetVisLeaves(kTRUE);
	PhotoTubeVolume->SetLineColor(kYellow-5);
	TGeoRotation rota("rot",10,20,30);
	TGeoTranslation trans;
	TGeoCombiTrans *c1 = new TGeoCombiTrans(trans,rota);
	
	/*
	Read in location of ALL phototubes
	*/
	//bool first=true;
	Double_t theta, phi;
	
    	for(int tubeId=0;tubeId<wcsimrootgeom->GetWCNumPMT();tubeId++)
    	{
    		WCSimRootPMT pmt = wcsimrootgeom -> GetPMT ( tubeId );
    		double pmtX = pmt.GetPosition (0);
    		double pmtY = pmt.GetPosition (1);
    		double pmtZ = pmt.GetPosition (2);
    		
    		int location=pmt.GetCylLoc();
    		theta=0.0;
    		phi=0.0;
    		float rad2deg=57.3;
    		if(location==1)
    		{
    			float lengthXY=sqrt(pmtX*pmtX+pmtY*pmtY);
    			float xd=pmtX/lengthXY;
    			float yd=pmtY/lengthXY;
    			TVector3 D(xd,yd,0.0);
    			theta=D.Theta()*rad2deg;
    			phi=D.Phi()*rad2deg;
    		}
    		else
    		{
    			theta=0.0;
    			phi=0.0;
    			if(location==1)theta=180.0;
    		}
    		if(location==0 || location ==1 || location ==2){
    			TGeoRotation TubeRotation("rotation",phi-90.0,theta,0.0);//D.Phi(),D.Theta,0.0);
    			TGeoTranslation PhototubePositionMatrix("shift",pmtX,pmtY,pmtZ);
    			TGeoCombiTrans *ShiftAndTwist = new TGeoCombiTrans(PhototubePositionMatrix,TubeRotation);
    			WorldVolume->AddNode(PhotoTubeVolume, tubeId, ShiftAndTwist);
    		}
    		
    	}
	// CLOSE GEOMETRY
	gGeoManager->CloseGeometry();
}
