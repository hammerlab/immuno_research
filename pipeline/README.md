## Immunotherapy Pipeline


### Todo

- VCF input support (VCF -> annotations -> transcript IDs -> insert mutation on trancript -> translation)
- MAF input support
- Check peptides against all self-peptides - save index instead of indexing every run
- Add proteasomal cleavage model (NetChop or SMM method (smm matrices are available but only seem appropriate for 6mers))
- Investigate TAP models


### Usage

```sh

python cancer_pipeline.py --input <fasta_file> --output <output_file> 

```

### Adding a new scorer

```python

from pipeline import PipelineElement

class NewScorer(PipelineElement):
  def __init__(name):

  def apply(self, data):
    return transform(data)

```
Then in `cancer_pipeline.py`

```python
pipeline.add_scorer(NewScorer(name="new scoring mech"))

```


