for i in {1..1000}; do
	slim -d REP="$i" "stickles - 4P.slim"
	python "localancestry_proportions.py" ./output/sticklebacks_ancestryproportions_$i.trees & test
done

