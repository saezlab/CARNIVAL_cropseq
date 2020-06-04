This code runs CARNIVAL on CRISPR-CAS9 perturbation data: [Crop-Seq](http://www.medical-epigenomics.org/papers/datlinger2017/).

## How to use: 
Change settings and all of the paths in <b>settings.yml</b> <br/> in the corresponding profile in the file (e.g. local) 
Then run in console:

### Local run
```
Rscript run_pipeline_cropseq.R -f settings_yml_path -l  
```

### Help with description of all options
```
Rscript run_pipeline_cropseq.R -h
```

## Contacts: 
olga.ivanova@bioquant.uni-heidelberg.de

## License:

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt]() or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).


## Tools used:
[Seurat](https://satijalab.org/seurat/) <br/>
[VIPER](https://bioconductor.org/packages/release/bioc/html/viper.html) 
[DoRothEA](https://saezlab.github.io/DoRothEA/) <br/>
[CARNIVAL](https://github.com/saezlab/CARNIVAL) <br/>

## Citations: 
[Carnival: Liu, Trairatphisan, Gjerga et al.](https://www.biorxiv.org/content/10.1101/541888v1) <br/>

[Crop-Seq paper](https://www.nature.com/articles/nmeth.4177)
Paul Datlinger, André F Rendeiro*, Christian Schmidl*, Thomas Krausgruber, Peter Traxler, Johanna Klughammer, Linda C Schuster, Amelie Kuchler, Donat Alpar, Christoph Bock (2017). Pooled CRISPR screening with single-cell transcriptome readout. Nature Methods. DOI: 10.1038/nmeth.4177. 
 
