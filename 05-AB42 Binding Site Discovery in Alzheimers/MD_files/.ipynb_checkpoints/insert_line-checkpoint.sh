#!/bin/bash

# Define the text to be inserted
insert_text='
; Include Position restraint file
#ifdef POSRES_A
#include "posre_Protein_chain_A.itp"
#endif

#ifdef POSRES_B
#include "posre_Protein_chain_B.itp"
#endif

#ifdef POSRES_C
#include "posre_Protein_chain_C.itp"
#endif

#ifdef POSRES_D
#include "posre_Protein_chain_D.itp"
#endif

#ifdef POSRES_E
#include "posre_Protein_chain_E.itp"
#endif

#ifdef POSRES_P
#include "posre_Protein_chain_P.itp"
#endif
'

# Use sed to insert the text after the specified line in the file
sed -i '/#include "topol_Protein_chain_P.itp"/ r /dev/stdin' topol.top <<<"$insert_text"

