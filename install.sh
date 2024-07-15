conda create env create -f function_annotation.yml
pip install goatools
pip install pandas
conda install -n function_annotation r-ggplot2
conda install -n function_annotation dplyr
path=$(conda info --base)
cp config.yml ${path}/envs/function_annotation/bin/.
