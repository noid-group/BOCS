#GO GO GO
#Check if args passed
if [ "$#" -lt 1 ] || [ ${1} == "-h" ]
then
        echo "PASS THE ITERATION WE'RE ON"
        exit 0
fi

iteration=${1}
backward=$(($iteration - 1))

mkdir iter"$iteration"
cd iter"$iteration"

cp ../cgwater.lmp .  # Template - LAMMPS Simulation file 
cp ../cgwater.data . # Sim input (template based)
cp ../config.ini .  # Template, I always add a 0 to the filename for the first iteration 
cp ../iter"$backward"/vol.dat . # Previous vol data
cp ../iter"$backward"/press.dat . # Previous press data
cp ../lammps_nb_SOLSOL.dat . # Sim input (template based)
cp ../lammps.sh . # Template - Cluster Submission Script

# Change previous psi
# I always add a 0 manually for launching the first round 
sed -i "s:NUMBER:iter"$backward"/psi.dat:g" config.ini

#Generate new Psi's: (Must have pmatching loaded)
~/group_storage/shared_data/Software/BOCS_Developers/pressure-matching/build/bin/p_match config.ini

#Read in Results
PSI_1=$(awk 'FNR == 4 {print}' psi.dat)
PSI_2=$(awk 'FNR == 5 {print}' psi.dat)

#Convert PSIs to correct units (Assuming units started in bar nm^3) 
NEWPSI_1=$(echo | awk -v awkpsi1="$PSI_1" '{print awkpsi1 * 986.9232667}')
NEWPSI_2=$(echo | awk -v awkpsi2="$PSI_2" '{print awkpsi2 * 986.9232667}')

echo $NEWPSI_1

#Place the converted PSIs in the lammps file

sed -i "s:UNO:"$NEWPSI_1":g" cgwater.lmp
sed -i "s:DOS:"$NEWPSI_2":g" cgwater.lmp

#Label submission script
sed -i "s:NAME:iter"$iteration":g" lammps.sh

sbatch lammps.sh

