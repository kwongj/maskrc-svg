# maskrc-svg
Masks recombinant regions in an alignment based on ClonalFrameML or Gubbins output  
Option to draw SVG of recombinant regions

## Authors
Jason Kwong (@kwongjc)  ::  [kwongj](https://github.com/kwongj)  
Torsten Seemann (@torstenseemann)  ::  [tseemann](https://github.com/tseemann)  

## Dependencies
* [Python 3.x](https://www.python.org/downloads/)
* [BioPython](http://biopython.org/wiki/Main_Page)
* [ete3](http://etetoolkit.org/)
* [bcbio-gff](https://github.com/chapmanb/bcbb/tree/master/gff)
* [svgwrite](https://pypi.python.org/pypi/svgwrite/)

## Usage
`$ maskrc-svg.py -h`  
```
usage: 
  maskrc-svg.py --aln FASTA --out OUTFILE [--gubbins] <PREFIX>

Mask recombination from ClonalFrameML/Gubbins output and draw SVG of recombinant regions

positional arguments:
  PREFIX               prefix used for CFML/Gubbins input files (required)

optional arguments:
  -h, --help           show this help message and exit
  --gubbins            parse as Gubbins instead of ClonalFrameML
  --aln FASTA          multiFASTA alignment used as input for CFML (required)
  --out OUTFILE        output file for masked alignment (default="maskrc.aln")
  --symbol CHAR        symbol to use for masking (default="?")
  --regions FILE       output recombinant regions to file
  --svg FILE           draw SVG output of recombinant regions and save as specified file
  --svgsize WIDExHIGH  specify width and height of SVG in pixels (default="800x600")
  --svgorder FILE      specify file containing list of taxa (1 per line) in desired order
  --svgcolour COLOUR   specify colour of recombination regions in HEX format (default=black)
  --consensus          add consensus row of recombination hotspots
  --version            show program's version number and exit
```

**Requires:**
* Output from ClonalFrameML or Gubbins
* MultiFASTA genome alignment used as input

**Options:**
* Specify output file using `--out OUTFILE`
* Specify symbol to use in the alignment for masking `--symbol ?`
* Save tab-separated file of recombinant region coordinates `--regions FILE`
* Draw SVG of recombinant regions `--svg FILE`
* Specify size of SVG in pixels eg. 800x600 `--svgsize WIDTHxHEIGHT`
* Specify desired order of taxa in SVG `--svgorder FILE`
* Specify colour to show extant recombination (ancestral recombination is shown in grey) `--svgcolour COLOUR`
* Add consensus row of recombination hotspots to SVG `--consensus`

## Bugs
Please submit via the GitHub issues page: [https://github.com/kwongj/maskrc-svg/issues](https://github.com/kwongj/maskrc-svg/issues)  

## Software Licence
GPLv3: [https://github.com/kwongj/maskrc-svg/blob/master/LICENSE](https://github.com/kwongj/maskrc-svg/blob/master/LICENSE)

## Other links
* [Gubbins](https://github.com/sanger-pathogens/gubbins)
* [ClonalFrameML](https://github.com/xavierdidelot/clonalframeml)
