#packages 
my_packages <- c('dplyr', 'tidyr', 'marked', 'ggplot2', 'R2ucare')
new_packages <- my_packages[!(my_packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages)
