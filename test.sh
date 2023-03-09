for i in {1..2}
do
	slim -d REP="$i" "stickles - Copy.slim"
	C:/ProgramData/Anaconda3/envs/SLiM_3.7/python.exe "localancestry_proportions - F.py" "E:/git/SLiM_stickles/output/sticklebacks_{$i}.trees" 
sleep 60
done

