name_base="PredictedDataReal"
run_num="05641"
space=" "
addstring=""
for i in {0..19}
	do
		filename=$name_base$run_num"_"$i".root"
		addstring=$addstring$space$filename
#		echo $filename
	done
hadd $name_base$run_num".root" $addstring

