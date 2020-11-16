## Filter neoepitopes based on threshold from Wells et al. (2020)

### Clone the repository:
```
git clone https://github.com/tanyaphung/neoantigens_prioritization.git
```

### Prequisites
- Python version 3.6.11
- Pandas version 1.1.3

### Directory structure
- Each sample has its own directory, for example: `patient_3466` (under the directory `example_data`)
    - Sub-directories inside the sample directory are the directories of the HLAs, for example: `HLA-A01-01`. The directory name of the HLA must match with the hla file for use as argument to the script.
        - Sub-directories inside the HLA directory are called `9_mers` or `10_mers` or `11_mers`

### Usage:

```
python neoepitope_prediction.py --hla_types_fn example_data/patient_3466/hla.txt --quant_fn example_data/patient_3466/3466_quant.sf --sample_id patient_3466 --mers 9 --data_dir example_data/
```
- This script outputs a file called `filtered_neoepitopes.tsv` inside the sample directory

**If there's no TPM file available:**
- If you don't have RNAseq data and therefore do not have the TPM file, you don't have to provide any input to `quant_fn`. Instead, this script will filter based on MHC binding affinity and MHC stability only.
```
python neoepitope_prediction.py --hla_types_fn example_data/patient_3466/hla.txt --sample_id patient_3466 --mers 9 --data_dir example_data/
```
- This script outputs a file called `filtered_neoepitopes_no_tpm.tsv` inside the sample directory

### Required inputs
- `--hla_types`: this is a file listing each HLA type per line. The name of the HLA type must match with the subdirectory inside the sample directory
- `--sample_id`: sample id must match the sample directory
- `--mers`: Enter a number (eg: 9) or a list of numbers 9,10,11
- `--data_dir`: Enter the parent directory that hosts all the data. For example, because the directory `patient_3466` is under the directory `example_data`, give the path to the `example_data` directory.

### Optional inputs
- `--quant_fn`: give the path to the output from salmon
- `binding_threshold`: Default is 34 nM (Wells et al. 2020). Input a different number if you want to change the default threshold
- `stability_threshold`: Default is 1.4h (Wells et al. 2020). Input a different number if you want to change the default threshold
- `tpm_threshold`: Default is 33TPM (Wells et al. 2020). Input a different number if you want to change the default threshold
