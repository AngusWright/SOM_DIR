# SOM_DIR

<p align="center">
  <strong>Self-organising-map direct calibration for photometric-redshift studies.</strong>
</p>

<p align="center">
  Train a self-organising map (SOM), match photometric and reference catalogues, and produce catalogue-level DIR weights for cosmic-shear redshift calibration.
</p>

<p align="center">
  <img src="Demo/Outputs/TrainingCat_zpaint.png" alt="Example SOM maps of reference and photometric redshift" width="800">
</p>

## What it does

`SOM_DIR.R` is a command-line R workflow for direct photometric-redshift calibration. It learns a SOM from supplied photometric features, maps a reference sample and a training sample onto the same cells, and derives `SOMweight` values that reweight the catalogues by their relative cell occupancy.

The workflow can also:

- create diagnostic maps and redshift plots;
- save trained SOMs for reuse;
- tune the number of hierarchical clusters;
- apply quality-control rejection; and
- process corresponding sets of training and reference catalogues.

## Quick start

### 1. Get the code and dependencies

```bash
git clone https://github.com/AngusWright/SOM_DIR.git
cd SOM_DIR
bash INSTALL.sh
```

You need a working installation of R and `Rscript` on your `PATH`. The install script obtains the required R packages, including `kohonen`, `data.table`, `FITSio`, and `helpRfuncs`.

### 2. Run the bundled example

The repository includes small example catalogues and their expected products. Recreate them with:

```bash
Rscript R/SOM_DIR.R \
  -r Demo/Inputs/ReferenceCat.csv \
  -t Demo/Inputs/TrainingCat.csv \
  --zr.label z_phot \
  --zt.label z_spec \
  -k MAG_u MAG_g MAG_r MAG_i MAG_g-MAG_r \
  -ct MAG_r \
  -cr MAG_r \
  -o Demo/Outputs/ \
  --som.dim 54 54 \
  -pp \
  --som.cores 100 \
  --som.iter 1000 \
  --force \
  --noqc
```

The generated weighted catalogues and diagnostics are written to `Demo/Outputs/`. The `--force` option is necessary here because that directory already contains example outputs.

## Inputs

Provide a reference catalogue with `-r` and a training catalogue with `-t`. Each can accept one or more corresponding catalogue paths.

| Input | Required content |
| --- | --- |
| Reference catalogue | A photometric-redshift column, supplied with `--zr.label`, plus every feature named after `-k`. |
| Training catalogue | A spectroscopic/training redshift column, supplied with `--zt.label`, plus every feature named after `-k`. |
| Both catalogues | Compatible feature columns or expressions for the SOM. Use `-cr` and `-ct` to select the reference and training count/weight variables when needed. |

The included example uses `z_phot` in the reference catalogue, `z_spec` in the training catalogue, four magnitude columns, and a colour expression (`MAG_g-MAG_r`).

## Outputs

For each input pair, SOM_DIR writes augmented catalogues in the selected output directory. Their names are based on the training catalogue name:

| Output | Contents |
| --- | --- |
| `*_DIRsom_*.csv` | Training catalogue with `GroupFactor` and `SOMweight`. |
| `*_refr_DIRsom_*.csv` | Reference catalogue with `GroupFactor` and `SOMweight`. |
| `*_SOMdata.Rdata` | Saved trained SOM. |
| `*_refr_SOMdata.Rdata` | Saved reference-sample SOM when requested. |
| `*.png` / `*.pdf` | Optional SOM, weight, occupancy, and redshift diagnostics. |

The output catalogue format follows the input extension; CSV, FITS, and RData outputs are supported.

## Common options

Run `Rscript R/SOM_DIR.R -h` for the full command reference and `Rscript R/SOM_DIR.R -hd` for default values.

| Option | Purpose |
| --- | --- |
| `-r`, `-t` | Reference and training catalogue paths. |
| `-k` | Feature columns or expressions used to train the SOM. |
| `--zr.label`, `--zt.label` | Reference photometric-redshift and training redshift labels. |
| `-cr`, `-ct` | Reference and training count/weight variables. |
| `-o`, `-of` | Output directory and output filename(s). |
| `--som.dim`, `--som.iter`, `--som.cores` | SOM geometry, iterations, and parallel workers. |
| `--old.som` | Reuse a previously saved SOM. |
| `--optimise` | Optimise the number of hierarchical clusters. |
| `--noqc` | Skip quality-control rejection. |
| `-p`, `-pp`, `-np` | Produce some plots, extensive plots, or no plots. |
| `--only.som` | Train and save the SOM without producing calibration weights (does not require reference catalogue). |

## Repository layout

```text
R/SOM_DIR.R       Main command-line workflow
Demo/Inputs/      Example reference and training catalogues
Demo/Outputs/     Example weighted catalogues, saved SOMs, and diagnostics
INSTALL.sh        Dependency installer
DEMO              Standalone copy of the demonstration command
```

## License

Distributed under the [MIT License](LICENSE).
