for f in /rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios/*.yml; do
    echo $f
    qR.rb -c 4 -t "01:00:00" -r R/sim-cmp-distances-cw.R --args  scenario_file=$f
done
