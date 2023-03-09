for i in {1..20}; do
	slim -d REP="$i" "stickles - Copy.slim"
	python "localancestry_proportions - F.py" ./output/sticklebacks_$i.trees & test
done

