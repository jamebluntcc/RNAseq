original_table=$1
plot_table=$2

awk '$0 ~ /[0-9]+:[0-9]+/' $original_table | sed -re 's/ \((\w+):(\w+)\)/\t\1\t\2/g' | cut -f1,4,5,8,9 > $plot_table
