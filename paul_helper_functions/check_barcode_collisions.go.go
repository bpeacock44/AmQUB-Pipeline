package main

import (
	"fmt"
	"os"
	"strings"
)

func main() {
	// get options from command line
	opt := make(map[string]string)
	opt["i"] = ""
	opt["m"] = ""
	opt["o"] = ""
	opt["M"] = "2"
	opt["t"] = ""
	opt["C"] = ""
	opt["f"] = ""
	opt["h"] = ""

	args := os.Args[1:]
	for i := 0; i < len(args); i++ {
		if strings.HasPrefix(args[i], "-") {
			opt[args[i][1:]] = args[i+1]
			i++
		}
	}

	if opt["i"] == "" || opt["m"] == "" {
		usage()
	}

	if opt["h"] != "" {
		usage()
	}

	opt["M"] = "2"

	// get dirs of input files
	mapdir := ""
	if strings.Contains(opt["m"], "/") {
		mapdir = opt["m"][:strings.LastIndex(opt["m"], "/")+1]
	} else {
		mapdir = "./"
	}

	fasqdir := ""
	fastqfile := ""
	if strings.Contains(opt["i"], "/") {
		fasqdir = opt["i"][:strings.LastIndex(opt["i"], "/")+1]
		fastqfile = opt["i"][strings.LastIndex(opt["i"], "/")+1:]
	} else {
		fasqdir = "./"
		fastqfile = opt["i"]
	}

	// get barcodes from mapping file
	mapBCs := getMappingFile(opt["m"])
	bclen := len(mapBCs[0])
	mapBCcount := len(mapBCs)

	// create an output filename, saving in the same directory as the fastq file
	if opt["o"] == "" {
		opt["o"] = fmt.Sprintf("%suFQBC_%s_BC%d_M%s.txt", fasqdir, fastqfile, mapBCcount, opt["M"])
	} else {
		if _, err := os.Stat(opt["o"]); err == nil && opt["f"] == "" {
			fmt.Printf("** Warning ** Output file [%s] already exists!\n", opt["o"])
			os.Exit(1)
		}
	}

	// make a barcode->SampleID lookup
	bcIDlu := make(map[string]string)
	for i := 0; i < len(mapBCs); i++ {
		bcIDlu[mapBCs[i]] = mapBCs[i]
	}

	fmt.Printf("Mapping file [%s] contains %d barcodes\n", opt["m"], mapBCcount)
	fmt.Printf("Map barcode lengths=%d\n", bclen)

	// determine input fastq type
	fqtype := ""
	fq, err := os.Open(opt["i"])
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer fq.Close()

	buf := make([]byte, 1024)
	n, err := fq.Read(buf)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	fqtype = "sequence"
	if strings.Contains(string(buf[:n]), ":") {
		fqtype = "barcode"
	}

	fmt.Printf("Fastq type: %s\n", fqtype)

	// get .fastq barcodes
	fmt.Printf("Reading fastq file [%s]\n", opt["i"])

	Fq := make(map[string]int)
	if fqtype == "sequence" {
		for {
			_, err := fq.Read(buf)
			if err != nil {
				break
			}

			lines := strings.Split(string(buf), "\n")
			for i := 1; i < len(lines); i += 4 {
				if lines[i] == "" {
					continue
				}
				barcode := lines[i][strings.LastIndex(lines[i], ":")+1:]
				Fq[barcode]++
			}
		}
	} else if fqtype == "barcode" {
		for {
			_, err := fq.Read(buf)
			if err != nil {
				break
			}

			lines := strings.Split(string(buf), "\n")
			for i := 1; i < len(lines); i += 4 {
				if lines[i] == "" {
					continue
				}
				barcode := lines[i]
				Fq[barcode]++
			}
		}
	}

	// get a list of all barcodes found
	FqBCs := make([]string, 0, len(Fq))
	for k := range Fq {
		FqBCs = append(FqBCs, k)
	}

	// determine which FqBCs have mismatches (1 .. $$opt{M}) to any mapping file barcode (BCs)
	fmt.Printf("Finding 0-%s mismatches to mapfile barcodes\n", opt["M"])
	C := findMismatches(mapBCs, FqBCs, opt["M"])

	// save usable fastq barcodes to file
	sings, err := os.Create(opt["o"])
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer sings.Close()

	fmt.Printf("Checking for usable and unusable barcodes and saving results: ")
	sings.WriteString(fmt.Sprintf("# %s [Version: %s]\n", os.Args[0], "VERSION"))
	sings.WriteString(fmt.Sprintf("# %d mapfile barcodes\n", mapBCcount))
	sings.WriteString("# Mapbc\tCount\t[SampleID]/\n")
	sings.WriteString("# FQbc\tCount\tmmPos\n")

	for _, bc := range mapBCs {
		// initialize any mapfile bcs not found in the .fastq file
		if _, ok := Fq[bc]; !ok {
			Fq[bc] = 0
		}

		bccount := Fq[bc]
		bcmmcount := make(map[int]int)

		sings.WriteString(fmt.Sprintf("\nbc\t%s\t%d\t[%s]\n", bc, Fq[bc], bcIDlu[bc]))

		// check each mismatch (1bpmm, 2bpmm, ...) up to whatever $$opt{M} is
		for d := 1; d <= len(opt["M"]); d++ {
			pmmm := fmt.Sprintf("pm%dm", d)
			fqbcs := sortSub(C[pmmm]["BC"][bc]["FQ"], Fq)

			de := d - 1
			dees := make([]int, de)
			for i := 0; i < de; i++ {
				dees[i] = i + 1
			}

			for _, fqbc := range fqbcs {
				// skip if it's a mapfile barcode (because those are handled by $bc)
				if _, ok := C["PMC"][fqbc]; ok {
					continue
				}

				mapbcs := make([]string, 0, len(C[pmmm]["FQ"][fqbc]["BC"]))
				for k := range C[pmmm]["FQ"][fqbc]["BC"] {
					mapbcs = append(mapbcs, k)
				}

				alreadyprinted := false
				for _, dee := range dees {
					if _, ok := C["multibc"][fqbc][dee]; ok {
						alreadyprinted = true
						break
					}
				}

				if alreadyprinted {
					continue
				}

				if len(mapbcs) > 1 {
					sings.WriteString(fmt.Sprintf("*%d\t%s\t%d\t%s", d, fqbc, Fq[fqbc], strings.Join(C[pmmm]["BC"][bc]["FQ"][fqbc], " ")))
					sings.WriteString("\t[ ")
					for _, mapbc := range mapbcs {
						sings.WriteString(fmt.Sprintf("%s:%s ", bcIDlu[mapbc], mapbc))
					}
					sings.WriteString("]\n")
				}

				if len(mapbcs) == 1 {
					sings.WriteString(fmt.Sprintf("m%d\t%s\t%d\t%s\n", d, fqbc, Fq[fqbc], strings.Join(C[pmmm]["BC"][bc]["FQ"][fqbc], " ")))
					bcmmcount[d] += Fq[fqbc]
				}
			}
		}

		bcmmtots := make([]int, 0, len(bcmmcount))
		for _, v := range bcmmcount {
			bcmmtots = append(bcmmtots, v)
		}

		colsNeeded := len(opt["M"]) - len(bcmmtots)
		if colsNeeded > 0 {
			bcmmtots = append(bcmmtots, make([]int, colsNeeded)...)
		}

		sings.WriteString(fmt.Sprintf("tot\t%s\t%d\t%s\t%s\n", bc, bccount, strings.Join(bcmmtots, "\t"), bcIDlu[bc]))
	}

	fmt.Printf("[%s]\n", opt["o"])
}

func usage() {
	message := `
Determine usable and unusable mismatched barcodes in a fastq file (compared to barcodes in a mapping file)
Usage: $0 [options]
   [REQUIRED]
   -i <input file>   # The input fastq filepath (barcodes OR reads with barcodes in headers. See NOTE below)
   -m <map file>     # Barcodes in mapping file will be denoted in the output
   [OPTIONAL]
   -o <output file>  # Output filepath (default=uFQBC_<input_fastq_filename>_BC<barcode_length>_M<mismatches>.txt
   -M <integer>      # Mismatches to check (default=2)
   -t <integer>      # Terminate after this many barcodes have been loaded (for quick checks)
   -C                # Print mismatch Collisions (fastq barcodes with >= 1 bp mismatch to multiple mapping file barcodes)
   -f                # Force overwrite of output file
   -h                # This help
NOTE: Read fastq ID lines must end with a ':' and a barcode (eg., \":CTCGACTACTGA\")
IMPORTANT: The map file should contain the barcodes for ALL SAMPLES present in the fastq file or the results may contain undetected errors.
`
	fmt.Println(message)
	os.Exit(0)
}

func getMappingFile(file string) []string {
	f, err := os.Open(file)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer f.Close()

	headers := make([]string, 0)
	var line string
	fmt.Fscanln(f, &line)
	headers = strings.Split(line, "\t")

	mapping := make([]string, 0)
	for {
		if _, err := fmt.Fscanln(f, &line); err != nil {
			break
		}
		mapping = append(mapping, strings.Split(line, "\t")[0])
	}

	return mapping
}

func findMismatches(BCs []string, FqBCs []string, M string) map[string]map[string]map[string][]int {
	C := make(map[string]map[string]map[string][]int)

	for _, bc := range BCs {
		C["PMC"][bc] = make(map[string][]int)

		for _, fqbc := range FqBCs {
			diffpos := findDiffs(bc, fqbc)

			if len(diffpos) == 0 {
				C["PMC"][fqbc] = append(C["PMC"][fqbc], 0)
			} else if len(diffpos) >= 1 && len(diffpos) <= len(M) {
				d := len(diffpos)
				pmmm := fmt.Sprintf("pm%dm", d)
				C[pmmm]["FQ"][fqbc]["BC"][bc] = diffpos
				C[pmmm]["BC"][bc]["FQ"][fqbc] = diffpos
				C["multibc"][fqbc][d]++
			}
		}
	}

	return C
}

func sortSub(unsorted map[string][]int, Fq map[string]int) []string {
	bcs := make([]string, 0, len(unsorted))
	for k := range unsorted {
		bcs = append(bcs, k)
	}

	idxs2sort := make([]int, len(bcs))
	for i := 0; i < len(bcs); i++ {
		idxs2sort[i] = i
	}

	sorted := make([]string, len(bcs))
	copy(sorted, bcs)

	for i := 0; i < len(bcs); i++ {
		for j := i + 1; j < len(bcs); j++ {
			for _, idx := range idxs2sort {
				cmp := unsorted[sorted[i]][idx] - unsorted[sorted[j]][idx]
				if cmp != 0 {
					return cmp
				}
			}
		}
	}

	return sorted
}

func findDiffs(a string, b string) []int {
	diffpos := make([]int, 0)

	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			diffpos = append(diffpos, i)
		}
	}

	return diffpos
}
