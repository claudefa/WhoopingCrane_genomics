# All modern, gen 13

java -cp /home/qvw641/bin/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder wc_two-epoch_realfold_4.blueprint

dir=/projects/mjolnir1/people/qvw641/WhoopingCrane/StairwayPlot/
jobName=$dir/out/stairwayplotFolded_4.sh
echo '#!/bin/bash' > $jobName
echo "bash wc_two-epoch_realfold_4.blueprint.sh" >> $jobName
sbatch -c 1 --mem-per-cpu 100G --time 23:00:00 -o ${dir}/out/stairwayplotFolded_13.log --job-name 13StairwayPlotF -- $jobName

# Only wild, gen 13

java -cp /home/qvw641/bin/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder wc_two-epoch_wild_u145.blueprint


dir=/projects/mjolnir1/people/qvw641/WhoopingCrane/StairwayPlot/
jobName=$dir/out/stairwayplot_u145.sh
echo '#!/bin/bash' > $jobName
echo "bash wc_two-epoch_wild_u145.blueprint.sh" >> $jobName
sbatch -c 1 --mem-per-cpu 100G --time 23:00:00 -o ${dir}/out/stairwayplot_u145.log --job-name 13Stairway145 -- $jobName


