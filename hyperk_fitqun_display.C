// Based on alice_esd.C from ROOT tutorials
// Modified for hyperk by Alex Finch
TFile* f;
TTree * fiTQun;

void load_fitqun_event();


/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/

//______________________________________________________________________________
void hyperk_fitqun_display()
{
	// Main function, initializes the application.
	
	gROOT->LoadMacro("selector.C");
	
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
	
	
    	fiTQun=(TTree *) f->Get("fiTQun");
    	TEveManager::Create(kTRUE,"V");
    	selector S(fiTQun);
    	S.Init(fiTQun);
    	S.Process(0);
    	//fiTQun->GetEvent(0);
    //	load_fitqun_event();
    	gEve->GetDefaultGLViewer()->UpdateScene();
    	gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.
}
void load_fitqun_event()
{
	cout<<" hi "<<endl;
}


