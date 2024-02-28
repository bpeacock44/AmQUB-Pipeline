```sh Now run the rest of this code without making changes. You can ignore the warning about "rank" being depreciated.
echo -e 'beta_diversity:metrics\'${B_MET} > bd_parameters.txt # creates parameter file
beta_diversity_through_plots.py -p bd_parameters.txt -i nc.${reg_asv_table} -m ${mapping_file} -o bdivs.nc.${reg_asv_table} # runs analysis
cd bdivs.nc.${reg_asv_table}
condense_multiple_PCOA_plots
cd ..
```
