// Mainframe macro generated from application: /opt/ROOT/6.18-04/bin/root.exe
// By ROOT version 6.18/04 on 2022-04-20 23:52:30

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGTextEntry
#include "TGTextEntry.h"
#endif
#ifndef ROOT_TRootCanvas
#include "TRootCanvas.h"
#endif
#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif

#include "Riostream.h"

void ()
{

   // main frame
   TGMainFrame *fRootCanvas49 = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);

   // menu bar
   TGMenuBar *fMenuBar59 = new TGMenuBar(fRootCanvas49,700,22,kHorizontalFrame);

   // "File" menu
   TGPopupMenu *fPopupMenu51 = new TGPopupMenu(gClient->GetDefaultRoot(),110,144,kOwnBackground);
   fPopupMenu51->AddEntry("&New Canvas",0);
   fPopupMenu51->AddEntry("&Open...",1);
   fPopupMenu51->AddEntry("&Close Canvas",13);
   fPopupMenu51->AddSeparator();

   // cascaded menu "Save"
   TGPopupMenu *fPopupMenu50 = new TGPopupMenu(gClient->GetDefaultRoot(),186,168,kOwnBackground);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&ps",5);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&eps",6);
   fPopupMenu50->AddEntry("c_total_occupancy_U.p&df",7);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&tex",11);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&gif",8);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&jpg",9);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&png",10);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&C",4);
   fPopupMenu50->AddEntry("c_total_occupancy_U.&root",3);
   fPopupMenu51->AddPopup("&Save",fPopupMenu50);
   fPopupMenu51->AddEntry("Save &As...",2);
   fPopupMenu51->AddSeparator();
   fPopupMenu51->AddEntry("&Print...",12);
   fPopupMenu51->AddSeparator();
   fPopupMenu51->AddEntry("&Quit ROOT",14);
   fMenuBar59->AddPopup("&File",fPopupMenu51, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Edit" menu
   TGPopupMenu *fPopupMenu53 = new TGPopupMenu(gClient->GetDefaultRoot(),74,144,kOwnBackground);
   fPopupMenu53->AddEntry("&Style...",15);
   fPopupMenu53->AddSeparator();
   fPopupMenu53->AddEntry("Cu&t",16);
   fPopupMenu53->DisableEntry(16);
   fPopupMenu53->AddEntry("&Copy",17);
   fPopupMenu53->DisableEntry(17);
   fPopupMenu53->AddEntry("&Paste",18);
   fPopupMenu53->DisableEntry(18);
   fPopupMenu53->AddSeparator();

   // cascaded menu "Clear"
   TGPopupMenu *fPopupMenu52 = new TGPopupMenu(gClient->GetDefaultRoot(),74,42,kOwnBackground);
   fPopupMenu52->AddEntry("&Pad",19);
   fPopupMenu52->AddEntry("&Canvas",20);
   fPopupMenu53->AddPopup("C&lear",fPopupMenu52);
   fPopupMenu53->AddSeparator();
   fPopupMenu53->AddEntry("&Undo",21);
   fPopupMenu53->DisableEntry(21);
   fPopupMenu53->AddEntry("&Redo",22);
   fPopupMenu53->DisableEntry(22);
   fMenuBar59->AddPopup("&Edit",fPopupMenu53, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "View" menu
   TGPopupMenu *fPopupMenu55 = new TGPopupMenu(gClient->GetDefaultRoot(),128,180,kOwnBackground);
   fPopupMenu55->AddEntry("&Editor",23);
   fPopupMenu55->AddEntry("&Toolbar",24);
   fPopupMenu55->AddEntry("Event &Statusbar",25);
   fPopupMenu55->AddEntry("T&oolTip Info",26);
   fPopupMenu55->AddSeparator();
   fPopupMenu55->AddEntry("&Colors",27);
   fPopupMenu55->AddEntry("&Fonts",28);
   fPopupMenu55->DisableEntry(28);
   fPopupMenu55->AddEntry("&Markers",29);
   fPopupMenu55->AddSeparator();
   fPopupMenu55->AddEntry("&Iconify",30);
   fPopupMenu55->AddSeparator();

   // cascaded menu "View With"
   TGPopupMenu *fPopupMenu54 = new TGPopupMenu(gClient->GetDefaultRoot(),78,42,kOwnBackground);
   fPopupMenu54->AddEntry("&X3D",31);
   fPopupMenu54->AddEntry("&OpenGL",32);
   fPopupMenu55->AddPopup("&View With",fPopupMenu54);
   fMenuBar59->AddPopup("&View",fPopupMenu55, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Options" menu
   TGPopupMenu *fPopupMenu56 = new TGPopupMenu(gClient->GetDefaultRoot(),151,216,kOwnBackground);
   fPopupMenu56->AddEntry("&Auto Resize Canvas",33);
   fPopupMenu56->CheckEntry(33);
   fPopupMenu56->AddEntry("&Resize Canvas",34);
   fPopupMenu56->AddEntry("&Move Opaque",35);
   fPopupMenu56->CheckEntry(35);
   fPopupMenu56->AddEntry("Resize &Opaque",36);
   fPopupMenu56->CheckEntry(36);
   fPopupMenu56->AddSeparator();
   fPopupMenu56->AddEntry("&Interrupt",37);
   fPopupMenu56->AddEntry("R&efresh",38);
   fPopupMenu56->AddSeparator();
   fPopupMenu56->AddEntry("&Pad Auto Exec",39);
   fPopupMenu56->AddSeparator();
   fPopupMenu56->AddEntry("&Statistics",40);
   fPopupMenu56->CheckEntry(40);
   fPopupMenu56->AddEntry("Histogram &Title",41);
   fPopupMenu56->CheckEntry(41);
   fPopupMenu56->AddEntry("&Fit Parameters",42);
   fPopupMenu56->AddEntry("Can Edit &Histograms",43);
   fMenuBar59->AddPopup("&Options",fPopupMenu56, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Tools" menu
   TGPopupMenu *fPopupMenu57 = new TGPopupMenu(gClient->GetDefaultRoot(),123,114,kOwnBackground);
   fPopupMenu57->AddEntry("&Inspect ROOT",44);
   fPopupMenu57->AddEntry("&Class Tree",45);
   fPopupMenu57->AddEntry("&Fit Panel",46);
   fPopupMenu57->AddEntry("&Start Browser",47);
   fPopupMenu57->AddEntry("&Gui Builder",48);
   fPopupMenu57->AddEntry("&Event Recorder",49);
   fMenuBar59->AddPopup("&Tools",fPopupMenu57, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Help" menu
   TGPopupMenu *fPopupMenu58 = new TGPopupMenu(gClient->GetDefaultRoot(),126,158,kOwnBackground);
   fPopupMenu58->AddLabel("Basic Help On...");
   fPopupMenu58->DefaultEntry(-1);
   fPopupMenu58->AddSeparator();
   fPopupMenu58->AddEntry("&Canvas",51);
   fPopupMenu58->AddEntry("&Menus",52);
   fPopupMenu58->AddEntry("&Graphics Editor",53);
   fPopupMenu58->AddEntry("&Browser",54);
   fPopupMenu58->AddEntry("&Objects",55);
   fPopupMenu58->AddEntry("&PostScript",56);
   fPopupMenu58->AddSeparator();
   fPopupMenu58->AddEntry("&About ROOT...",50);
   fMenuBar59->AddPopup("&Help",fPopupMenu58, new TGLayoutHints(kLHintsRight | kLHintsTop));
   fRootCanvas49->AddFrame(fMenuBar59, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,0,0,1,1));
   TGHorizontal3DLine *fHorizontal3DLine67 = new TGHorizontal3DLine(fRootCanvas49,700,2);
   fRootCanvas49->AddFrame(fHorizontal3DLine67, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // dockable frame
   TGDockableFrame *fDockableFrame68 = new TGDockableFrame(fRootCanvas49);

   // next lines belong to the dockable frame widget
   fDockableFrame68->EnableUndock(kTRUE);
   fDockableFrame68->EnableHide(kFALSE);
   fDockableFrame68->SetWindowName("ToolBar: ");
   fDockableFrame68->DockContainer();

   fRootCanvas49->AddFrame(fDockableFrame68, new TGLayoutHints(kLHintsExpandX));
   TGHorizontal3DLine *fHorizontal3DLine73 = new TGHorizontal3DLine(fRootCanvas49,700,2);
   fRootCanvas49->AddFrame(fHorizontal3DLine73, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // composite frame
   TGCompositeFrame *fCompositeFrame74 = new TGCompositeFrame(fRootCanvas49,700,476,kHorizontalFrame);

   // composite frame
   TGCompositeFrame *fCompositeFrame75 = new TGCompositeFrame(fCompositeFrame74,175,476,kFixedWidth);
   fCompositeFrame75->SetLayoutManager(new TGVerticalLayout(fCompositeFrame75));

   fCompositeFrame74->AddFrame(fCompositeFrame75, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   // canvas widget
   TGCanvas *fCanvas76 = new TGCanvas(fCompositeFrame74,700,476,kSunkenFrame);

   // canvas viewport
   TGViewPort *fViewPort77 = fCanvas76->GetViewPort();

   // canvas container
   Int_t canvasID = gVirtualX->InitWindow((ULong_t)fCanvas76->GetId());
   Window_t winC = gVirtualX->GetWindowID(canvasID);
   TGCompositeFrame *fCompositeFrame86 = new TGCompositeFrame(gClient,winC,fViewPort77);
   fViewPort77->AddFrame(fCompositeFrame86);
   fCompositeFrame86->SetLayoutManager(new TGVerticalLayout(fCompositeFrame86));
   fCompositeFrame86->MapSubwindows();
   fCanvas76->SetContainer(fCompositeFrame86);
   fCanvas76->MapSubwindows();
   fCompositeFrame74->AddFrame(fCanvas76, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY));

   fRootCanvas49->AddFrame(fCompositeFrame74, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   // status bar
   TGStatusBar *fStatusBar89 = new TGStatusBar(fRootCanvas49,10,18);
   Int_t partsusBar89[] = {33,10,10,47};
   fStatusBar89->SetParts(partsusBar89,4);
   fRootCanvas49->AddFrame(fStatusBar89, new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandX,2,2,1,1));

   fRootCanvas49->SetIconPixmap("macro_s.xpm");
   fRootCanvas49->SetClassHints("ROOT","Canvas");
   fRootCanvas49->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fRootCanvas49->MapSubwindows();
   fHorizontal3DLine67->UnmapWindow();
   fDockableFrame68->UnmapWindow();
   fHorizontal3DLine73->UnmapWindow();
   fCompositeFrame75->UnmapWindow();
   fStatusBar89->UnmapWindow();

   fRootCanvas49->Resize(fRootCanvas49->GetDefaultSize());
   fRootCanvas49->MapWindow();
   fRootCanvas49->Resize(700,500);
}  
