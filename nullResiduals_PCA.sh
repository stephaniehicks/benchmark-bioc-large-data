#$ -pe local 2
#$ -l mem_free=5G,h_vmem=7G
#$ -cwd
#$ -o log/
#$ -e log/
#$ -m e
#$ -M shicks19@jhu.edu
module load conda_R/4.0.x

run_id="stephanie_cluster"
data_name="tenxbraindata"
data_type="inmem"
fam="poisson"
type="deviance"
num_workers=2

Rscript --slave nullResiduals_PCA.R --args $run_id $data_name $data_type $fam $type $num_workers