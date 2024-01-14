#include <PalmOS.h>
#include "resids.h"
#include "u_struct.h"

#define UNUSED  __attribute__((__unused__))

#define md_Select 0
#define md_Draw 1

int orgx,orgy,lastx,lasty;
Boolean drawstarted;
TMolecule* mol;
TAtom* curatom=NULL;

RectangleType selrect;
RectangleType WorkingRec;


TAtom* newatom;
TBond* bond;
int i;
MemHandle fldHandle;
FieldType* fldALabel;
ListType *list;
char *ps;
char curdbname[32];

int curmode=md_Draw;

TAtom* FindNearbyAtom(TMolecule* mol, int x, int y)
{
 int i;
 TAtom* atom;

  for (i=0;i<mol->atoms->Count;i++)
    {
     atom=Get_Atom(mol->atoms,i);
     if ((abs(atom->x-x)<5) && (abs(atom->y-y)<5)) return atom;
    }
 return NULL;
}


void Draw_Atom(TAtom *atom)
{
  if (StrCompare(atom->ALabel,"C"))
    WinDrawChars(atom->ALabel,StrLen(atom->ALabel),atom->x-FntLineWidth(atom->ALabel,StrLen(atom->ALabel))/2,atom->y-FntLineHeight()/2);
}


int sign(int n)
{
 if (n>=0) return 1; else return -1;
}

void GetShifts(int dx,int dy,int *shiftx,int *shifty) // Get X,Y shifts of 1st point along the specified line
{
 float dydx;
 if (dx==0)
 {
  *shiftx=0;
  *shifty=sign(dy);
  return;
 }

 dydx=abs(dy/dx);

 if (dydx<0.4) { *shiftx=sign(dx); *shifty=0; } // tan(22.5 deg)
 else if (dydx<2.4) { *shiftx=sign(dx); *shifty=sign(dy); } // tan(67.5 deg)
 else { *shiftx=0; *shifty=sign(dx); }
}


void Draw_Molecule(TMolecule* mol)
{
 int i, shiftx,shifty;
 TBond *bond;

 for (i=0;i<mol->bonds->Count;i++)
 {
  bond=Get_Bond(mol->bonds,i);

  if (bond->btype==bt_double)
  {
   GetShifts(-(bond->ToAtom->y-bond->FromAtom->y),bond->ToAtom->x-bond->FromAtom->x,&shiftx,&shifty);
   WinDrawLine(bond->FromAtom->x-shiftx,bond->FromAtom->y-shifty,bond->ToAtom->x-shiftx,bond->ToAtom->y-shifty);
   WinDrawLine(bond->FromAtom->x+shiftx,bond->FromAtom->y+shifty,bond->ToAtom->x+shiftx,bond->ToAtom->y+shifty);
  } 
  else if (bond->btype==bt_triple) 
  {
   GetShifts(-(bond->ToAtom->y-bond->FromAtom->y),bond->ToAtom->x-bond->FromAtom->x,&shiftx,&shifty);
   shiftx*=2;
   shifty*=2;
   WinDrawLine(bond->FromAtom->x-shiftx,bond->FromAtom->y-shifty,bond->ToAtom->x-shiftx,bond->ToAtom->y-shifty);
   WinDrawLine(bond->FromAtom->x,bond->FromAtom->y,bond->ToAtom->x,bond->ToAtom->y);
   WinDrawLine(bond->FromAtom->x+shiftx,bond->FromAtom->y+shifty,bond->ToAtom->x+shiftx,bond->ToAtom->y+shifty);
  } 
  else WinDrawLine(bond->FromAtom->x,bond->FromAtom->y,bond->ToAtom->x,bond->ToAtom->y);
 }

 for (i=0;i<mol->atoms->Count;i++)
  Draw_Atom(Get_Atom(mol->atoms,i));
 
}

void HideSelection()
{
 if (selrect.topLeft.x==-1) return;

 WinInvertRectangleFrame(simpleFrame,&selrect);
}

void DrawSelectedAtom(TAtom* atom)
{
 if (!atom) { selrect.topLeft.x=-1; return;}

 Draw_Atom(atom);

if (StrCompare(atom->ALabel,"C")==0)
{
 selrect.topLeft.x=atom->x-2;
 selrect.topLeft.y=atom->y-2;
 selrect.extent.y=5;
 selrect.extent.x=5;
} else
{
 selrect.topLeft.x=atom->x-FntLineHeight()/2;
 selrect.topLeft.y=atom->y-FntLineHeight()/2;
 selrect.extent.y=FntLineHeight();
 selrect.extent.x=FntLineHeight();
}


 WinInvertRectangleFrame(simpleFrame,&selrect);
 
}

void SetCurAtom(TAtom* atom)
{
 HideSelection();
 curatom=atom; 
 DrawSelectedAtom(curatom);
}

void Redraw()
{
 HideSelection();
 WinEraseRectangle(&WorkingRec,0);
 Draw_Molecule(mol);
 DrawSelectedAtom(curatom);
}


Boolean GetNewDBName(char *dbname)
{
 for (i=0;i<100;i++) 
 { 
  StrPrintF(dbname,"Molecule%d",i);
  if (DmFindDatabase (0, dbname)==0) return true; 
 }
 return false;
}

void EraseAll()
{
 LocalID dbid;
 if ((dbid=DmFindDatabase (0, curdbname))!=0) DmDeleteDatabase (0, dbid);
 Clear_Molecule(mol);
 curatom=NULL;
 Redraw();
}


Boolean SaveToDB(char *dbname)
{
 int err;
 DmOpenRef dbRef;
 LocalID dbid;

 if (mol->atoms->Count==0) return false;

// If exists, delete
 
 if ((dbid=DmFindDatabase (0, dbname))!=0) DmDeleteDatabase (0, dbid);

// Create

 err = DmCreateDatabase (0, dbname, myCreator, myDBType, false);
 if (err) return false;


 dbRef=DmOpenDatabase (0, DmFindDatabase (0, dbname), dmModeReadWrite);


// dbRef = DmOpenDatabaseByTypeCreator(myCreator, myDBType, dmModeReadWrite);

 if (! dbRef) return false;


 dmWrite_Molecule(dbRef, mol);

 DmCloseDatabase (dbRef);
 return true;
}

void LoadFromDB(char *dbname)
{
 LocalID dbid;
 DmOpenRef dbRef;
 if (!(dbid=DmFindDatabase (0, dbname))) return;

 dbRef=DmOpenDatabase (0, dbid, dmModeReadOnly);

 dmRead_Molecule(dbRef,mol);
 DmCloseDatabase (dbRef);

 StrCopy(curdbname,dbname);
 curatom=NULL;
 Redraw();
}

	
Boolean MainFormHandler(EventType *e)
{
  FormType *form = FrmGetActiveForm();

//  FldHandleEvent(FrmGetObjectPtr(form,FrmGetObjectIndex(form, fld_Label)), e));


  switch (e->eType) {

  case frmOpenEvent:
//    if (FtrGet('gdbS', 0, &feature) == 0 && feature == 0x12BEEF34)
//      CtlSetValue(FrmGetObjectPtr(form, FrmGetObjectIndex(form, EnableB)), 1);

    selrect.topLeft.x=-1;
    WinGetBounds (FrmGetWindowHandle(FrmGetActiveForm()),&WorkingRec);

    WorkingRec.topLeft.x+=11;
    WorkingRec.extent.y-=10;

    mol=New_Molecule();

    fldALabel=FrmGetObjectPtr(form,FrmGetObjectIndex(form, fld_Label));
    list=FrmGetObjectPtr(form,FrmGetObjectIndex(form, List_ID));

   fldHandle=MemHandleNew(4);
   StrCopy(MemHandleLock(fldHandle),"O");
   MemHandleUnlock(fldHandle);

    FldSetTextHandle(fldALabel,fldHandle);
    
    GetNewDBName(curdbname);

    FrmDrawForm(form);

    return false;

  case appStopEvent: 
         SaveToDB(curdbname);
         Free_Molecule(mol);
       break;

  case menuEvent: 
        switch (e->data.menu.itemID) {
         case Menu_New:
               break;
         case Menu_Eraseall:
               EraseAll();
               break;
         case Menu_Open: 
                         {
			  int n, err;	
// Find all molecule databases

 DmSearchStateType state;
 UInt16 card;
 LocalID currentDB = 0;
 TList* names;
 MemPtr dbname;
 int i;

 SaveToDB(curdbname);

 names=New_List();
 err = DmGetNextDatabaseByTypeCreator(true, &state, myDBType, myCreator, false, &card, &currentDB);
 while (!err && currentDB) {
 dbname=MemPtrNew(32);
 List_Add(names,dbname);
 err= DmDatabaseInfo (card, currentDB, (char*)dbname, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

// UInt16 *attributesP,
//UInt16 *versionP, UInt32 *crDateP,
//UInt32 *modDateP, UInt32 *bckUpDateP,
//UInt32 *modNumP, LocalID *appInfoIDP,
//LocalID *sortInfoIDP, UInt32 *typeP,
//UInt32 *creatorP)

 
 err = DmGetNextDatabaseByTypeCreator(false, &state, myDBType, myCreator, false, &card, &currentDB);
 }




  			  LstSetListChoices(list,(char**)names->Items,names->Count);
	                  n=LstPopupList (list);

			  if (n!=-1) LoadFromDB((char*)names->Items[n]);

// Free name list
	                  for (i=0;i<names->Count;i++) MemPtrFree(names->Items[i]);
                          Free_List(names);
			 }
                         break;

         case Menu_Send: {
                          ExgSocketType sk;
                          MemPtr data;
                          Err er;

                          data=GetMol(mol);

                          MemSet(&sk,sizeof(sk),0);

                          sk.count=1;
                          sk.length=StrLen(data);
                          sk.type="text/plain";
                          sk.name="?_send;_local:a.mol";
                          if (ExgPut(&sk)==errNone)
                          {
                           ExgSend(&sk,data,StrLen(data),&er);
                           ExgDisconnect(&sk,er);
                          }

                          MemPtrFree(data);
                         }
                        break;
         case Menu_About: FrmAlert(ID_About); break;
	                
        }
        break;

  case ctlSelectEvent:
    switch (e->data.ctlSelect.controlID) {
     case btn_Select: curmode=md_Select; break;  // Button

     case btn_Draw: curmode=md_Draw; break;  // Button

     case btn_Delete: if (!curatom) break;
       
      
//      select an atom connected to this one if any

     newatom=curatom;
     if (curatom->bonds->Count>0) curatom=Bond_NextAtom(Get_Bond(curatom->bonds,0),curatom);
     if (curatom==newatom) curatom=NULL;

     Molecule_DeleteAtom(mol,newatom);

     Redraw();

     break;  // Button

     case btn_Label: 
      if (!curatom) break;
      fldHandle=FldGetTextHandle(fldALabel);
      ps=MemHandleLock(fldHandle);
      if (*ps=='\0') StrCopy(curatom->ALabel,"C"); 
      else StrNCopy(curatom->ALabel,ps,4);      
      MemHandleUnlock(fldHandle);
      SetCurAtom(curatom);
     break;
    }
    break;

  case keyDownEvent:

//   FldHandleEvent(FrmGetObjectPtr(form,FrmGetObjectIndex(form, fld_Label)), e);

/* struct _KeyDownEventType {
WChar chr;
UInt16 keyCode;
UInt16 modifiers;
};
chr The character code.
keyCode Unused.
modifiers 0, or one or more of the following values:
shiftKeyMask Graffiti is in case-shift mode.
capsLockMask Graffiti is in cap-shift mode.
numLockMask Graffiti is in numeric-shift mode.
commandKeyMask The Graffiti glyph was the menu
command glyph or a virtual key code.
optionKeyMask Not implemented. Reserved.
controlKeyMask Not implemented. Reserved.
autoRepeatKeyMask Event was generated due to auto-repeat.
doubleTapKeyMask Not implemented. Reserved.
poweredOnKeyMask The key press caused the system to be
powered on.
appEvtHookKeyMask System use only.
libEvtHookKeyMask System use only. */
  break;

  case penDownEvent:
         if (!RctPtInRectangle(e->screenX,e->screenY,&WorkingRec)) break;
         SetCurAtom(FindNearbyAtom(mol,e->screenX,e->screenY));
         orgx=e->screenX; orgy=e->screenY;
         drawstarted=false;
  break;

  case penMoveEvent:
         if (!RctPtInRectangle(e->screenX,e->screenY,&WorkingRec)) break;

   if (curmode==md_Select)
   {
   if (curatom) {
    curatom->x=e->screenX;
    curatom->y=e->screenY;
    Redraw();
    }
    break;    
   }

 // only start drawing if pen moved at least 3 pixels away

   if (drawstarted) {} else
   if ((e->screenX-orgx)*(e->screenX-orgx)+(e->screenY-orgy)*(e->screenY-orgy)>9)
   {
    drawstarted=true;
    lastx=-1;
    if (!curatom)
    {     
     curatom=New_Atom();
     List_Add(mol->atoms,curatom);
     curatom->x=orgx;
     curatom->y=orgy;     
     SetCurAtom(curatom);
    }
   } else break; // not yet drawing

   if (lastx!=-1) WinInvertLine(curatom->x,curatom->y,lastx,lasty); // hide former line

    
    if ((newatom=FindNearbyAtom(mol,e->screenX,e->screenY))) { lastx=newatom->x; lasty=newatom->y;}
    else { lastx=e->screenX; lasty=e->screenY;}
   
    WinInvertLine(curatom->x,curatom->y,lastx,lasty) ;
 
  break;

  case penUpEvent: 
    
   if (!RctPtInRectangle(e->screenX,e->screenY,&WorkingRec)) break;

   if (curmode!=md_Draw) break;

   if (!drawstarted) break;
   
   newatom=FindNearbyAtom(mol,e->screenX,e->screenY);
   if (newatom==curatom) break;
  
   if (!newatom)
   {
     newatom=New_Atom();
     List_Add(mol->atoms,newatom);
     newatom->x=e->screenX;
     newatom->y=e->screenY;
   }

 // check if atoms already have a bond. If so, change its order


   for (i=0;i<newatom->bonds->Count;i++)
   {
    bond=Get_Bond(newatom->bonds,i);

    if (Bond_NextAtom(bond,newatom)==curatom)
    {
     switch (bond->btype) {
      case bt_single: bond->btype=bt_double; break;
      case bt_double: bond->btype=bt_triple; break;
     default: bond->btype=bt_single;
     }
     Redraw();
     return false;
    }

   }


     bond=New_Bond();
     List_Add(mol->bonds,bond);
     bond->FromAtom=curatom;
     bond->ToAtom=newatom;
     List_Add(curatom->bonds,bond);
     List_Add(newatom->bonds,bond);
     curatom=newatom;
     SetCurAtom(curatom);

  break;
  
  default:
    break;
  }

  return false;
}

Boolean MiscHandler(EventType *e)
{
  switch (e->eType) {

  case frmOpenEvent:
    return false;

  default:
    return false;
  }
}

UInt32 PilotMain(UInt16 cmd, void *cmdPBP UNUSED, UInt16 launchFlags UNUSED)
{
  EventType e;
  UInt16 err;

  if (cmd == sysAppLaunchCmdNormalLaunch) {

    FrmGotoForm(MainForm);

    do {
      EvtGetEvent(&e, evtWaitForever);

      if (SysHandleEvent(&e)) continue;
      if (MenuHandleEvent(NULL, &e, &err)) continue;
      if (e.eType == frmLoadEvent) {
	FormType *form = FrmInitForm(e.data.frmLoad.formID);
	FrmSetActiveForm(form);
	FrmSetEventHandler(form, MainFormHandler);
	continue;
      }
      if (FrmDispatchEvent(&e)) continue;
      if (MiscHandler(&e)) continue;
    } while (e.eType != appStopEvent);

    FrmCloseAllForms();
  }

  return 0;
}

// void FrmPopupForm (UInt16 formId)
// void FrmReturnToForm (UInt16 formId) 0 - last form

/*
	if (event->eType == frmLoadEvent)
	{
		// Load the form resource specified in  event then activate the form.
		formId = event->data.frmLoad.formID;
		frm = FrmInitForm(formId);
		FrmSetActiveForm(frm);

		// Set the event handler for the form.  The handler of the currently 
		// active form is called by FrmDispatchEvent each time it receives  
		// an event.
		switch (formId)
		{
		case OrderForm:
			FrmSetEventHandler(frm, OrderHandleEvent);
			break;
			
		case CustomersForm:
			FrmSetEventHandler(frm, CustomersHandleEvent);
			break;
			
		}
		handled = true;
	}
*/

/*
UInt16 volRefNum;
UInt32 volIterator = vfsIteratorStart;
while (volIterator != vfsIteratorStop) {
err = VFSVolumeEnumerate(&volRefNum, &volIterator);
if (err == errNone) {
// Do something with the volRefNum
} else {
// handle error... possibly by
// breaking out of the loop
}
}*/

/*Err VFSFileOpen (UInt16 volRefNum,
const Char *pathNameP, UInt16 openMode,
FileRef *fileRefP)*/


// ErrNonFatalDisplayIf (condition, msg)

// ErrTry {
// Do something which may fail. Call ErrThrow to signal
// failure and force a jump to the following ErrCatch
// block.
//}
// ErrCatch(inErr) {
// Recover or cleanup after a failure in the above ErrTry
// block. "inErr" is an exception code identifying the
// reason for the failure.
// Call ErrThrow if you want to jump out to the next
// ErrCatch block.
// The code in this block doen’t execute if the above
// ErrTry block completes without a call to ErrThrow.
//} ErrEndCatch

//ErrDisplay (msg)
//Call this macro to display an error message, source code
//and line number. This macro is compiled into the code only
//compiler define ERROR_CHECK_LEVEL is set to 1 or 2
//(ERROR_CHECK_PARTIAL or ERROR_CHECK_FULL).

//ErrAlert (err)
