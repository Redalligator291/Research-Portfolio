# AB42 Drug Binding Project  

## Objective  
This project aims to propose potential binding sites on AB42 for nine unique compounds tested in clinical trials. These compounds are identified as candidates for modifying or disrupting AB42 nucleation and aggregation, contributing to Alzheimer's disease therapeutics and diagnostics.  

## Reference Paper  
The ligands utilized in this project are referenced from the following paper:  
Kim, Hye Yun, and YoungSoo Kim.  
*"Chemical-Driven Amyloid Clearance for Therapeutics and Diagnostics of Alzheimer’s Disease."*  
*Accounts of Chemical Research (2024): 8997.*  

Please read this paper to understand the hypotheses proposed for each drug target.  

---

## Compounds and Hypotheses  

| Compound          | Hypothesized Interaction                                                                                     | SMILES Code                                                                                  |
|--------------------|-------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| **Oxytetracycline** | Interacts with Aβ dimers                                                                                    | `C[C@@]1([C@H]2[C@@H]([C@H]3[C@@H](C(=O)C(=C([C@]3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O)O` |
| **Sunitinib**      | Interacts with Aβ dimers                                                                                    | `CCN(CC)CCNC(=O)C1=C(NC(=C1C)/C=C\2/C3=C(C=CC(=C3)F)NC2=O)C`                                |
| **Rhizolutin**     | Binds to AB42 fibrils (plaque, aggregates)                                                                  | `CC[C@H]1C[C@H]2[C@@H](C=C(C=C[C@@H]3C[C@@H](CC(=O)O[C@H]3C=C2C)O)C)C(=O)O1`                |
| **Borrelidin**     | Binds to AB42 fibrils (plaque, aggregates)                                                                  | `C[C@H]1C[C@H](C[C@@H]([C@H](/C(=C\C=C\C[C@H](OC(=O)C[C@@H]([C@H](C1)C)O)[C@@H]2CCC[C@H]2C(=O)O)/C#N)O)C)C` |
| **YB-09**          | Binds to AB42 fibrils and dissolves toxic AB42 oligomers (dimers or trimers)                                | `COC1=CC=C(C=C1)C1=COC2=CC(OCC(O)=O)=C(OC)C=C12`                                              |
| **EPPS**           | Dissolves toxic AB42 oligomers (dimers or trimers)                                                         | `C1CN(CCN1CCCS(=O)(=O)O)CCO.O` |
| **Quinacrine**     | Dissolves toxic AB42 oligomers (dimers or trimers)                                                         | `CCN(CC)CCCC(C)NC1=C2C=C(C=CC2=NC3=C1C=CC(=C3)Cl)OC`                                        |
| **YIAD-0205**      | Interacts with beta-sheets and specific binding sites                                                      | `C1=CC=C(C=C1)C2=C(N3C=CC=C3C4=NC5=CC(=C(C=C5N24)Br)Br)CC6=CC=C(C=C6)Cl`                  |
| **YIAD-0121**      | Targets specific hydrophobic sites within AB42, including 14-20, 28-35, and 36-42                          | `FC1=CC=C(C=C1)C(=O)C1C(N=CC2=CC=CN12)C1=CC=CO1`                                           |

---

Since EPPS is particularly notable for its ability to cross the blood-brain barrier, we suggest starting with this drug.  
The literature hypothesizes that EPPS binds with oligomer species, specifically dimers or trimers.  
We recommend focusing on the six conformers of the trimer to evaluate how EPPS binds.  

---

## AB42 Species  

We have collected AB42 species in monomeric, trimeric, and pentameric forms. The files are located in the `AB42` folder and follow this naming convention:  

- `specie_conformernumber.pdb`  
  - Example: `Monomer_1.pdb`  

---
## Updates on the Docking 

Currently the SwissDock doesn't support Sunitinib and Borrelidin docking. Switch to DiffDock.
