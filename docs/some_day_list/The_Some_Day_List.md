# The "some day" List

Stuff that I'd like to get to at some point that isn't important enough to prioritize *now*.

## Development
- Cram convenience/plotting functions into a R package for ease of maintenance/sharing.
- Find a way to summarize lots of GO terms into a human interpretable output that is both succinct and doesn't suck.
	- Fine tune an LLM and use RAG with the term descriptions?
	- Preferably from within R without relying on external setup to run models.
- Write up CRISPRball for JOSS or such & submit.
- Make IBET more performant.
	- Avoid complete re-plotting when adding labels.
	- Add option to add labels automatically on highlighted gene(sets).
- Make the nf-core Cut & Run pipeline script for fragment length histograms more memory efficient by changing from using `np.array` for each iteration to a list and concatenating the list objects into an array at the end.
	- This sucker takes like 500+ GB of memory for large datasets unnecessarily.

## Data Processing
- Typical `params.yaml` files for nf-core pipelines pointing to our reference data rather than having that crammed into configs.
- Template notebooks for ATAC-seq.
- Template notebooks for ChIP-seq/CnR.

## Reference Data 
- Get genesets from: https://www.science.org/doi/10.1126/science.add7046
- Make combined reference genomes for nf-core pipelines, e.g. mm10 + hg38. GTFs and FASTAs, with clear species indicators.
	- At least for RNA-seq, this takes advantage of the fractional counts based on probability (EM) from salmon for ambiguous reads.
## Organization
- List all interactive reports/apps on wiki.