name_base="PredictedDataReal"
run_num="05641_"
space=" "
addstring=""
for i in {0..19}
	do
		filename=$name_base$run_num$i".root"
		$addstring=$addstring$space$filename
		echo $filename
	done
hadd $name_base".root" $addstring

