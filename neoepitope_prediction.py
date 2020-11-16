import argparse
import os
import pandas as pd
pd.options.mode.chained_assignment = None

class ParseExpression:
    def __init__(self, quant_fn):
        self.quant_fn = quant_fn
    def parse_expression(self):
        """

        :return: a dictionary where key is the gene/transcript name and value is TPM
        """
        tpm_dict = {}
        with open(self.quant_fn, "r") as f:
            for line in f:
                if not line.startswith("Name"): #TODO: check header
                    cols = line.rstrip("\n").split("\t")
                    gene_name = cols[0].split(".")[0][4:]
                    tpm_dict[gene_name] = float(cols[4])
        return tpm_dict

class FilterPeptide:
    def __init__(self, binding_fn, stability_fn, hla, sample_id):
        self.binding_fn = binding_fn
        self.stability_fn = stability_fn
        self.hla = hla
        self.sample_id = sample_id
    def clean_gene_name(self, gene_name):
        """
        Clean gene name
        :param gene_name:
        :return:
        """
        return gene_name.split("_")[0][1:]
    def binding_stability_filtering(self, binding_thresold, stability_threshold):
        """
        Filtering based on MHC binding affinity and MHC stability by using pandas
        :return:
        """
        binding_df = pd.read_csv(self.binding_fn, skiprows=1, sep="\t")
        # print(len(binding_df.index))
        binding_df_rmdups = binding_df.drop_duplicates(subset=["Peptide", "ID"])
        # print(len(binding_df_rmdups.index))


        stability_df = pd.read_csv(self.stability_fn, skiprows=1, sep="\t")
        # print(len(stability_df.index))
        stability_df_rmdups = stability_df.drop_duplicates(subset=["Peptide", "ID"])
        # print(len(stability_df_rmdups.index))

        combined = pd.merge(binding_df_rmdups, stability_df_rmdups, how="inner", on=["Peptide", "ID"])
        # print(combined.head())
        # print(len(combined.index))

        filtered = combined[(combined["nM"] < binding_thresold) & (combined["Thalf(h)"] > stability_threshold)]
        filtered_clean = filtered[["Peptide", "ID", "nM", "Thalf(h)"]]

        return filtered_clean

    def tpm_filtering(self, tpm_dict, filtered_df, tpm_threshold):
        if not tpm_dict:
            return None
        filtered_df["gene"] =  filtered_df.apply(lambda x: self.clean_gene_name(x["ID"]), axis=1)

        filtered_df["tpm"] = filtered_df.apply(lambda x: tpm_dict.get(x["gene"]), axis=1)

        filter_df_tpm = filtered_df[(filtered_df["tpm"] > tpm_threshold)]

        filter_df_tpm["HLA"] = self.hla

        return filter_df_tpm

def main(args):

    # Check if there's the TPM file available
    if args.quant_fn:
        parsing_expression = ParseExpression(args.quant_fn)
        tpm_dict = parsing_expression.parse_expression()

    # Proces the HLAs into a list
    with open(args.hla_types_fn, "r") as f:
        hlas = [line.rstrip("\n") for line in f]

    # Get the list for mers
    mers = args.mers.split(",")

    for hla in hlas:
        for mer in mers:
            netmhc_fn = os.path.join(args.data_dir, args.sample_id, hla, mer + "_mers", "netmhc.xsl")
            netmhcstab_fn = os.path.join(args. data_dir, args.sample_id, hla, mer + "_mers", "netmhcstab.xsl")
            if not os.path.exists (netmhc_fn) or not os.path.exists(netmhcstab_fn): #Check if there's a file for netmhc and netmhcstab
                print("The netmhc file and/or netmhc stability file could not be found. Please check the readme for the directory structure. Existing")
                exit()
            else:
                filtering_peptide = FilterPeptide(binding_fn=netmhc_fn, stability_fn=netmhcstab_fn, hla=hla,
                                                  sample_id=args.sample_id)
                if args.quant_fn: #If there's the TPM file available
                    filtered_df = filtering_peptide.binding_stability_filtering(args.binding_threshold, args.stability_threshold)
                    results = filtering_peptide.tpm_filtering(tpm_dict=tpm_dict, filtered_df=filtered_df, tpm_threshold=args.tpm_threshold)
                    results.to_csv(os.path.join(args.data_dir, args.sample_id, hla + "_" + mer + "_mers_filtered_neoepitopes.tsv"), sep="\t", index=False)
                else:
                    results = filtering_peptide.binding_stability_filtering(args.binding_threshold, args.stability_threshold)
                    results.to_csv(os.path.join(args.data_dir, args.sample_id, hla + "_" + mer + "_mers_filtered_neoepitopes_no_tpm.tsv"), sep="\t", index=False)

def parse_args():
    parser = argparse.ArgumentParser(description="Filter neoepitopes based on threshold defined by Wells et al. 2020")
    parser.add_argument("--hla_types_fn", required=True, help="Path to a file listing HLA types (each HLA per line). See readme for directory structure.")
    parser.add_argument("--quant_fn", required=False, help="Path to salmon output")
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--data_dir", required=True, help="Path to the directory where the data is. This is the parent directory to the sample directory")
    parser.add_argument("--mers", type=str, required=True, help="Enter a number (eg: 9) or a list of numbers (9,10,11 or \"9, 10, 11\") ")
    parser.add_argument("--binding_threshold", required=False, type=float, default=float(34), help="Enter a number for binding threshold. The default is 34 nM")
    parser.add_argument("--stability_threshold", required=False, type=float, default=float(1.4),
                        help="Enter a number for stability threshold. The default is 1.4 hr")
    parser.add_argument("--tpm_threshold", required=False, type=float, default=float(33),
                        help="Enter a number for binding threshold. The default is 33 TPM")

    return parser.parse_args()

main(parse_args())