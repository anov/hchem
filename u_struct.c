#include <PalmOS.h>

#include "u_struct.h"

void ErrNumFmt(char *fmt, int num)
{
 char err[32];
 StrPrintF(err,fmt,num);
 ErrDisplay(err);
}

TList* New_List()
{
 MemHandle mh;
 TList* list;
 list=(TList*)MemPtrNew(sizeof(TList));
 mh=MemHandleNew(16*sizeof(MemPtr));
 list->Count=0;
 list->Capacity=16; 
 list->mh=mh;
 list->Items=MemHandleLock(mh);
 return list;
}

void Free_List(TList* list)
{
 MemHandleUnlock(list->mh);
 MemHandleFree(list->mh);
 MemPtrFree((MemPtr)list);
}


void List_Add(TList* list,void* item)
{
 if (list->Count<list->Capacity) {}
 else
 {
  list->Capacity+=16;
  MemHandleUnlock(list->mh);
  ErrFatalDisplayIf(MemHandleResize(list->mh,list->Capacity*sizeof(MemPtr))!=errNone,"Couldn't expand the list");  
  list->Items=MemHandleLock(list->mh);
 }
 list->Items[list->Count++]=item;
}

MemPtr Get_Item(TList* list,int n)
{
 if ((n<0) || (n>=list->Count)) ErrNumFmt("Invalid index %d",n);
 return list->Items[n];
}

void List_Clear(TList* list)
{
 list->Count=0;
}


// List_get()


void List_Delete(TList* list, int n)
{
 MemMove(list->Items+n,list->Items+n+1,(list->Count-n)*sizeof(MemPtr));
 list->Count--;
}

int  List_IndexOf(TList* list, void *item)
{
 int i;
 for (i=0;i<list->Count;i++) {
  if (list->Items[i]==item) return i;
 }
 return -1;
}


void Stack_Push(TStack* stack, int p)
{
 List_Add(stack,(void*)(Int32)p);
}

int Stack_Pop(TStack *stack)
{
 int r;
 r=(Int32)stack->Items[stack->Count-1];
 List_Delete(stack,stack->Count-1);
 return r;
}

TAtom* Bond_NextAtom(TBond* bond, TAtom* from)
{
 if (from==bond->FromAtom) return bond->ToAtom; else return bond->FromAtom;
}

TBond* New_Bond()
{
 TBond* bond;
 bond=(TBond*)MemPtrNew(sizeof(TBond));
 bond->btype=bt_single;
 return bond;
}

TAtom* New_Atom()
{
 TAtom* atom;
 atom=(TAtom*)MemPtrNew(sizeof(TAtom));
 atom->bonds=New_List();
 StrCopy(atom->ALabel,"C");
 return atom;
}

void Free_Atom(TAtom* atom)
{
 Free_List(atom->bonds);
 MemPtrFree(atom);
}


void Clear_Molecule(TMolecule* mol)
{
 int i;
 for (i=0;i<mol->bonds->Count;i++) Free_Bond(Get_Bond(mol->bonds,i));
 List_Clear(mol->bonds);
 for (i=0;i<mol->atoms->Count;i++) Free_Atom(Get_Atom(mol->atoms,i));
 List_Clear(mol->atoms);
}


void Molecule_Aromatize(TMolecule *mol) // marks all bonds between sp2 atoms aromatic
{
 int i;
 TBond* bond;
 for (i=0;i<mol->bonds->Count;i++)
 {
  bond=Get_Bond(mol->bonds,i);
  if ((bond->FromAtom->hybridization==h_sp2) && (bond->ToAtom->hybridization==h_sp2)) bond->btype=bt_aromatic;
 }
}


TMolecule* New_Molecule()
{
 TMolecule* mol;
 mol=(TMolecule*)MemPtrNew(sizeof(TMolecule));
 mol->bonds=New_List();
 mol->atoms=New_List();
 mol->atomsplaced=true;
 return mol;
}

void Free_Molecule(TMolecule* mol)
{
 Clear_Molecule(mol);
 Free_List(mol->bonds);
 Free_List(mol->atoms);
 MemPtrFree((MemPtr)mol); 
}


void Molecule_DeleteAtom(TMolecule * mol,TAtom* atom)
{
 int i,n;
 TBond* bond;


 n=List_IndexOf(mol->atoms,atom);
 if (n==-1) return;
 List_Delete(mol->atoms,n);

 for (i=mol->bonds->Count-1;i>=0;i--)
 {
  bond=Get_Bond(mol->bonds,i);
  if ((bond->FromAtom==atom) || (bond->ToAtom==atom))
  {
   n=List_IndexOf(bond->FromAtom->bonds,bond); if (n!=-1) List_Delete(bond->FromAtom->bonds,n);
   n=List_IndexOf(bond->ToAtom->bonds,bond); if (n!=-1) List_Delete(bond->ToAtom->bonds,n);
   List_Delete(mol->bonds,i);
   Free_Bond(bond);
  }
 }

 Free_Atom(atom);
}

MemPtr GetMol(TMolecule* mol)
{
 int size,i;
 MemPtr data;
 TAtom* atom;
 TBond* bond;
 size=StrLen("MOL")+2+StrLen("CDRAW")+2+0+2+40+2+mol->atoms->Count*(69+2)+mol->bonds->Count*(21+2)+StrLen("M  END")+1;
 data=MemPtrNew(size);
 if (!data) return NULL;
 StrCopy(data,"MOL\r\nCDRAW\r\n\r\n");
 StrPrintF(data+StrLen(data),"%3d%3d%3d%3d%3d%3d%3d%3d%3d  0999 V2000\r\n",mol->atoms->Count,mol->bonds->Count,0,0,0,0,0,0,0);

 for (i=0;i<mol->atoms->Count;i++)
 {
  atom=Get_Atom(mol->atoms,i);
  atom->num=i+1;
  StrPrintF(data+StrLen(data),"%5d.0   %5d.0       0.0    %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\r\n",atom->x,atom->y,atom->ALabel,0,0,0,0,0,0,0,0,0,0,0,0);
 }

 for (i=0;i<mol->bonds->Count;i++)
 {
  int bo;
  bond=Get_Bond(mol->bonds,i);
  switch (bond->btype)
  {
   case bt_single: bo=1; break;
   case bt_double: bo=2; break;
   case bt_triple: bo=3; break;
   case bt_aromatic: bo=4; break;
   default: bo=8;
  }

  StrPrintF(data+StrLen(data),"%3d%3d%3d%3d      \r\n",bond->FromAtom->num,bond->ToAtom->num,bo,0);
 }

 StrCat(data,"M  END");

 return data;
}


void dmWrite_Molecule(DmOpenRef db,TMolecule* mol)
{
 short d; 
 int i, base;
 UInt16 nrec;
 MemHandle hrec;
 MemPtr prec;
 TAtom *atom;
 TBond *bond;

ErrTry{

// Record 1 - version, natoms, nbonds

 nrec = dmMaxRecordIndex;
 hrec = DmNewRecord(db, &nrec, sizeof(short)*3);
// if ( !hrec ) return;

 prec=MemHandleLock(hrec);

 d=1;
 DmWrite(prec,0,&d,sizeof(d));

 d=mol->atoms->Count;
 DmWrite(prec,sizeof(short),&d,sizeof(d));

 d=mol->bonds->Count;
 DmWrite(prec,sizeof(short)*2,&d,sizeof(d));

 MemHandleUnlock(hrec); 
 DmReleaseRecord(db, nrec, 1);

// number atoms and bonds

for (i=0;i<mol->atoms->Count;i++)
{
 Get_Atom(mol->atoms,i)->num=i;
}

for (i=0;i<mol->bonds->Count;i++)
{
 Get_Bond(mol->bonds,i)->num=i;
}

// record 2 atoms: label, x,y

 nrec = dmMaxRecordIndex;
 hrec = DmNewRecord(db, &nrec, (sizeof(TALabel)+sizeof(short)*2)*mol->atoms->Count);

 if ( !hrec ) return;

 prec=MemHandleLock(hrec);

for (i=0;i<mol->atoms->Count;i++)
{
 base=i*(sizeof(atom->ALabel)+sizeof(short)*2);

 atom=Get_Atom(mol->atoms,i);
 DmStrCopy(prec,base,atom->ALabel);
 d=atom->x;
 DmWrite(prec,base+sizeof(TALabel),&d,sizeof(short));
 d=atom->y;
 DmWrite(prec,base+sizeof(TALabel)+sizeof(short),&d,sizeof(short));
}
 MemHandleUnlock(hrec); 
 DmReleaseRecord(db, nrec, 1);


//  record 3, bonds. Each - type, from, to

 nrec = dmMaxRecordIndex;
 hrec = DmNewRecord(db, &nrec, sizeof(short)*3*mol->bonds->Count);
 if ( !hrec ) return;
 prec=MemHandleLock(hrec);

for (i=0;i<mol->bonds->Count;i++)
{
 base=i*(sizeof(short)*3);
 bond=Get_Bond(mol->bonds,i);
 d=bond->btype;
 DmWrite(prec,base,&d,sizeof(short));
 d=bond->FromAtom->num;
 DmWrite(prec,base+sizeof(short),&d,sizeof(short));
 d=bond->ToAtom->num;
 DmWrite(prec,base+sizeof(short)*2,&d,sizeof(short));
}

 MemHandleUnlock(hrec); 
 DmReleaseRecord(db, nrec, 1);

}

ErrCatch(inErr) {
 ErrAlert(inErr);
} ErrEndCatch

}

void dmRead_Molecule(DmOpenRef db,TMolecule* mol)
{
 int i, nbonds,natoms;
 MemHandle hrec;
 short* prec;
 TAtom *atom;
 TBond *bond;

 Clear_Molecule(mol);

ErrTry{

// record 1, molecule header


 hrec=DmGetRecord (db, 0); prec=MemHandleLock(hrec);

 natoms=*(prec+1);
 nbonds=*(prec+2);

 MemHandleUnlock(hrec);  DmReleaseRecord(db, 0, 0);

 if (natoms==0) return;

// record 2, atoms


 hrec=DmGetRecord (db, 1); prec=MemHandleLock(hrec);


for (i=0;i<natoms;i++)
{
 atom=New_Atom();
 List_Add(mol->atoms,atom);
 StrNCopy(atom->ALabel,(char*)prec,10);
 atom->x=*(prec+5);
 atom->y=*(prec+6);
 prec+=7;
}

 MemHandleUnlock(hrec); DmReleaseRecord(db, 1, 0);

// record 3, bonds


 hrec=DmGetRecord (db, 2); prec=MemHandleLock(hrec);


for (i=0;i<nbonds;i++)
{
 bond=New_Bond();
 List_Add(mol->bonds,bond);

 bond->btype=*prec;
 
 bond->FromAtom=Get_Atom(mol->atoms, *(prec+1));

 bond->ToAtom=Get_Atom(mol->atoms, *(prec+2));

 List_Add(bond->FromAtom->bonds,bond);

 List_Add(bond->ToAtom->bonds,bond);
 prec+=3; 

}

 MemHandleUnlock(hrec); DmReleaseRecord(db, 2, 0);
}

ErrCatch(inErr) {
  ErrDisplay ("Error cauth");
 ErrAlert(inErr);
} ErrEndCatch

}
