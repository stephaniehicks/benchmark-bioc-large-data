# -pe local 1
# -l mem_free=10G,h_vmem=12G
#$ -l mem_free=10G,h_vmem=12G
#$ -cwd
#$ -o log/
#$ -e log/
#$ -m e
#$ -M shicks19@jhu.edu
module load conda_R/4.0.x

run_id="stephanie_cluster"
data_name="tenxbraindata"
data_type="ondisk"
fam="poisson"
model_type="pearson"
num_workers=1
mode="mem"

Rscript --slave nullResiduals_PCA.R --args $run_id $data_name $data_type $fam $model_type $num_workers $mode