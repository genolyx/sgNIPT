# WES path parity with carrier-screening

Single Gene NIPT shares the same biological WES step (FASTQ → aligned BAM → germline-style VCF on targets) as the **carrier-screening** pipeline (sibling repo, e.g. `carrier_screening/`). The **execution model** differs: carrier-screening processes often bootstrap **micromamba** inside each task (for nested Docker / portable installs); sgNIPT runs inside **one conda-based image** and calls the same tools from `PATH`.

## Behaviour alignment (implemented in sgNIPT)

| Step | carrier-screening | sgNIPT |
|------|-------------------|--------|
| Align | `bwa mem \| samtools sort -m 768M` | Same (`params.samtools_sort_mem`, default `768M`) |
| Mark duplicates | `gatk MarkDuplicates` | `picard MarkDuplicates` (same Picard engine); `use_picard_markdup=true` by default |
| BAM stats | `samtools stats` | `SAMTOOLS_STATS` (`samtools stats`) |
| Flagstat | (downstream) | `SAMTOOLS_FLAGSTAT` |
| Coverage | mosdepth + samtools bedcov / depth | `MOSDEPTH` + panel-specific `bam_qc.py` |
| **Germline VCF** | `HaplotypeCaller` + `VariantFiltration` | `variant_caller=gatk`: same HaplotypeCaller + VariantFiltration filters |

## Choosing variant caller

- **`bcftools`** (default): `mpileup` + `call` — lighter, good for cfDNA NIPT pipelines.
- **`gatk`**: carrier-screening **CALL_VARIANTS** equivalent (HaplotypeCaller + VariantFiltration on `target_bed`).
- **`mutect2`**: legacy Mutect2 path (not typical for germline WES; use only if you need it).

Example:

```bash
nextflow run main.nf ... --variant_caller gatk
```

## Shared Nextflow modules (options)

carrier-screening’s `.nf` files are **not** drop-in for sgNIPT because they embed micromamba bootstrap and host-specific paths (`/home/sam/.cache/...`). To share logic across repos:

1. **Extract a third repo** (e.g. `genomics-nf-modules`) with processes that only use `params`, `task`, and tools on `PATH` (no per-task conda install). Both pipelines `include` from a path or a submodule.
2. **Git submodule**: add the shared repo under `sgNIPT/submodules/nf-wes` and `include { ... } from "${projectDir}/submodules/nf-wes/align.nf"`.
3. **Symlink** (single machine): `ln -s /home/ken/carrier_screening/bin/modules/align.nf /home/ken/sgNIPT/bin/modules/align_carrier.nf` — fragile; prefer submodule or published package.
4. **Monorepo**: one repository with `pipelines/sgNIPT`, `pipelines/carrier-screening`, `nf-modules/wes/` shared by both.

**Nextflow `include` does not resolve packages from PyPI**; sharing is done via **git submodules**, **relative paths**, or **copying** maintained modules.

## Maintenance decision

**sgNIPT and carrier-screening stay separate repositories.** Behaviour is kept aligned where it matters (documented in this file); there is no shared Nextflow submodule or monorepo merge planned. Revisit only if duplicate maintenance becomes a real burden.
