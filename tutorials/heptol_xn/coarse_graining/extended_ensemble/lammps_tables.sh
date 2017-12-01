#!/bin/bash


#table_a0.xvg  table_b0.xvg  table_b1.xvg  table_b2.xvg  table_CB_CB.xvg  table_CB_CM.xvg  table_CB_CT.xvg  table_CF_CB.xvg  table_CF_CF.xvg  table_CF_CM.xvg  table_CF_CT.xvg  table_CM_CM.xvg  table_CT_CM.xvg  table_CT_CT.xvg

#ERROR: Too few arguments.  Accepted modes are:
#'translate_table.py [gromacs_table] [table_type] [lammps_table] [interaction_name]'
#[table_type] = {bond, nb, angle, dih}


translate_table.py table_a0.xvg angle lammps_angle_CT-CM-CT.table LJLJ

translate_table.py table_b0.xvg bond lammps_bond_CT-CM.table LJLJ

translate_table.py table_b1.xvg bond lammps_bond_CF-CB.table LJLJ

translate_table.py table_b2.xvg bond lammps_bond_CB-CB.table LJLJ

translate_table.py table_CB_CB.xvg nb lammps_nb_CBCB.table LJLJ

translate_table.py table_CB_CM.xvg nb lammps_nb_CBCM.table LJLJ

translate_table.py table_CB_CT.xvg nb lammps_nb_CBCT.table LJLJ

translate_table.py table_CF_CB.xvg nb lammps_nb_CFCB.table LJLJ

translate_table.py table_CF_CF.xvg nb lammps_nb_CFCF.table LJLJ

translate_table.py table_CF_CM.xvg nb lammps_nb_CFCM.table LJLJ

translate_table.py table_CF_CT.xvg nb lammps_nb_CFCT.table LJLJ

translate_table.py table_CM_CM.xvg nb lammps_nb_CMCM.table LJLJ

translate_table.py table_CT_CM.xvg nb lammps_nb_CTCM.table LJLJ

translate_table.py table_CT_CT.xvg nb lammps_nb_CTCT.table LJLJ

