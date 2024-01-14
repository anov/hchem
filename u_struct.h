#include <PalmOS.h>


#define bt_single 1
#define bt_double 2
#define bt_triple 4
#define bt_aromatic 8
#define bt_nobond 0

#define h_unknown 0
#define h_sp3 1
#define h_sp2 2
#define h_sp 3

typedef struct
{
 int Count;
 int Capacity;
 MemHandle mh; 
 MemPtr* Items;
} TList, TStack;

typedef char TALabel[10];

typedef struct {
     int tag;
     TALabel ALabel;
     TList * bonds;
     int hybridization;
     int num;
     int x;
     int y;
} TAtom;

typedef struct {
     int btype;
     TAtom* FromAtom;
     TAtom* ToAtom;
     int num;
     int MatchingBond;
} TBond;

typedef struct {
     TList* atoms;
     TList* bonds;
     Boolean atomsplaced;
} TMolecule;


#define Free_Bond(bond) MemPtrFree(bond)

TList* New_List();
void List_Add(TList* list,void* item);
void List_Clear(TList* list);
void List_Delete(TList* list, int n);
int  List_IndexOf(TList* list, void *item);

void Free_List(TList* list);
TAtom* Bond_NextAtom(TBond* bond,TAtom* atom);
void Molecule_Aromatize(TMolecule *mol); // marks all bonds between sp2 atoms aromatic
TMolecule* LoadSmiles(char* str);
char* GetSmiles(TMolecule* mol);
TMolecule* LoadMol(char* fname);
TMolecule LoadMolString(char* str);
MemPtr GetMol(TMolecule* mol);
int FindSubStructure(TMolecule* mol, TMolecule* inmol);
TAtom* New_Atom();
TBond* New_Bond();
void Free_Atom(TAtom* atom);
void Clear_Molecule(TMolecule* mol);
void Free_Molecule(TMolecule* mol);
TMolecule* New_Molecule();
void Molecule_DeleteAtom(TMolecule * mol,TAtom* atom);

void dmWrite_Molecule(DmOpenRef db,TMolecule* mol);
void dmRead_Molecule(DmOpenRef db,TMolecule* mol);

MemPtr Get_Item(TList* list,int n);

#define Get_Bond(bonds,n) ((TBond*)Get_Item(bonds,n))
#define Get_Atom(atoms,n) ((TAtom*)Get_Item(atoms,n))
#define Get_Molecule(molecules,n) ((TMolecule*)Get_Item(molecules,n))

