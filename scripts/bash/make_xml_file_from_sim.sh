

Rscript scripts/R/make_xml_file_from_sim_cells.R outputs/sim_CC_cells.csv outputs/CC_final_sampled_cells.csv outputs/CC_sim_tmp.xml

sed -n '/<!-- The patterns from 1 to end/,$p' outputs/sim_template.xml | awk 'NR!=1{print $0}' >> outputs/CC_sim_tmp.xml

sed -e 's/CC_sim_test/CC_sim/g' outputs/CC_sim_tmp.xml > outputs/CC_sim.xml
