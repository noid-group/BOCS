#!/bin/bash

EXE="tables"

mkdir table

INPUT="f_forces.Pair_Interaction.total.CFCF.dat"
OUTPUT="table/table_CF_CF.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.Pair_Interaction.total.CFCB.dat"
OUTPUT="table/table_CF_CB.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.Pair_Interaction.total.CBCB.dat"
OUTPUT="table/table_CB_CB.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.BondStretch.total.Bond_CF-CB.dat"
OUTPUT="table/table_b1.xvg"
TYPE="bond"
BASIS="linear"
RMIN="0.0"
RMAX="0.800"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMIN $RMAX


INPUT="f_forces.BondStretch.total.Bond_CB-CB.dat"
OUTPUT="table/table_b2.xvg"
TYPE="bond"
BASIS="linear"
RMIN="0.0"
RMAX="0.800"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMIN $RMAX

INPUT="f_forces.Pair_Interaction.total.CTCT.dat"
OUTPUT="table/table_CT_CT.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.Pair_Interaction.total.CMCM.dat"
OUTPUT="table/table_CM_CM.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.Pair_Interaction.total.CTCM.dat"
OUTPUT="table/table_CT_CM.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.BondStretch.total.Bond_CT-CM.dat"
OUTPUT="table/table_b0.xvg"
TYPE="bond"
BASIS="linear"
RMIN="0.0"
RMAX="0.800"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMIN $RMAX


INPUT="f_forces.Angle.total.Angle_CT-CM-CT.dat"
OUTPUT="table/table_a0.xvg"
TYPE="angle"
BASIS="linear"
$EXE $INPUT $OUTPUT $TYPE $BASIS


INPUT="f_forces.Pair_Interaction.total.CTCF.dat"
OUTPUT="table/table_CF_CT.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="f_forces.Pair_Interaction.total.CMCF.dat"
OUTPUT="table/table_CF_CM.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="f_forces.Pair_Interaction.total.CTCB.dat"
OUTPUT="table/table_CB_CT.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


INPUT="f_forces.Pair_Interaction.total.CMCB.dat"
OUTPUT="table/table_CB_CM.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="6.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

