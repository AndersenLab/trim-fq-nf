# trimmomatic-nf

[Documentation](http://andersenlab.org/dry-guide/pipeline-trimming/)

## Debug locally

You can test the pipeline locally with the following command. Note that you must
use dsl2 with version `20.01.0` of nextflow.

```bash
NXF_VER=20.01.0 nextflow run main.nf -resume -profile debug
```