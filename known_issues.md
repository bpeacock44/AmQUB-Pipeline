# Known Issues

Sometimes there will be an incomplete read while retrieving taxonomy during the taxonomy assignment step in part 4. You will see an error like this one:

```
#Traceback (most recent call last):
  #File "/sw/miniconda3/lib/python3.10/http/client.py", line 566, in _get_chunk_left
    #chunk_left = self._read_next_chunk_size()
  #File "/sw/miniconda3/lib/python3.10/http/client.py", line 533, in _read_next_chunk_size
    #return int(line, 16)
#ValueError: invalid literal for int() with base 16: b''
#
#During handling of the above exception, another exception occurred:
#
#Traceback (most recent call last):
  #File "/sw/miniconda3/lib/python3.10/http/client.py", line 583, in _read_chunked
    #chunk_left = self._get_chunk_left()
  #File "/sw/miniconda3/lib/python3.10/http/client.py", line 568, in _get_chunk_left
    #raise IncompleteRead(b'')
#http.client.IncompleteRead: IncompleteRead(0 bytes read)
#
#During handling of the above exception, another exception occurred:
#
#Traceback (most recent call last):
  #File "/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/blast_assign_taxonomy.py", line 1160, in <module>
    #main(__doc__)
  #File "/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/blast_assign_taxonomy.py", line 1134, in main
    #taxonomy_dict = get_taxonomy(opts, asv_list)
  #File "/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/blast_assign_taxonomy.py", line 603, in get_taxonomy
    #download_eposted_taxonIDs_to_XML(opts)
  #File "/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/blast_assign_taxonomy.py", line 546, in download_eposted_taxonIDs_to_XML
    #data = fetch_handle.read()
  #File "/sw/miniconda3/lib/python3.10/http/client.py", line 460, in read
    #return self._read_chunked(amt)
  #File "/sw/miniconda3/lib/python3.10/http/client.py", line 598, in _read_chunked
    #raise IncompleteRead(b''.join(value))
#http.client.IncompleteRead: IncompleteRead(1478171 bytes read)
#Error on line blast_assign_taxonomy.py -i "${output_dir}/asvs/r
```

I find that I can run the script again once or twice and this is resolved, so if you get this error please run the part 4 script again with the -s option (to skip the blast).