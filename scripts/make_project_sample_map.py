#!/usr/bin/env python3

import subprocess, argparse, os, sys, re
import numpy as np
import pandas as pd
from pathlib import Path

CONFIG = {
	"RUN_PLATE_FORMATTERS": {
		"IDT for Illumina": lambda row, run: ".".join([run, "IDT" + row["primer plate name"].split("XTR_Plate ")[1]]),
		"Nextera XT": lambda row, run: ".".join([run, row["XT barcode plate"].replace("set ", "UDI")])
	},
	"RUN_FORMATTER": lambda x: re.sub(pattern="[^.a-zA-Z0-9]", repl=".", string=x),
	"CTRL_WELLS": ["E12", "F12", "G12", "H12"]
}

def main(args):
	# main columns: ["gDNA plate ID", "WELL", "SAMPLE ID"]
	sample_to_gdna = get_gdna(args.gdna)
	blanks = pd.isna(sample_to_gdna["SAMPLE ID"])

	if True in blanks.value_counts().index:
		print("Removing " + str(blanks.value_counts()[True]) + " wells with no sample. If this is too many, check that the 'SAMPLE ID' column is filled.\n")
	sample_to_gdna = sample_to_gdna.loc[~blanks, ]
	if args.debug:
		print("gDNA-to-barcode plates:")
		print(sample_to_gdna)
	
	if(args.manifest):
		sample_to_gdna = apply_manifest(sample_to_gdna, args.manifest)
		if args.debug: 
			print("After applying manifest:")
			print(sample_to_gdna)
	
	gdna_to_barcode_plate = get_pooling(args.pooling)
	if args.debug:
		print("Pooling detail:")
		print(gdna_to_barcode_plate)
	project_map = gdna_to_barcode_plate.merge(sample_to_gdna, how="inner", on=["gDNA plate ID"])
	if sample_to_gdna.shape[0] > project_map.shape[0]:
		plate_diff = set(sample_to_gdna["gDNA plate ID"]).difference(
			set(gdna_to_barcode_plate["gDNA plate ID"])
			)
		print("gDNA plates not in pooling details: " + ",".join(plate_diff))
		quit("Error: Lost rows when joining gDNA maps to pooling detail. gDNA plate IDs did not match.")
	if args.debug: 
		print("Joined map:")
		print(project_map)
	
	project_map.loc[:, "RUN.PLATEPOSITION"] = project_map.apply(
		lambda row: ".".join([row["RUN.PLATE"], row["WELL"]]), 
		axis=1
		)
	project_map = project_map.rename(columns={"SAMPLE ID": "sampleID"})
	project_map["sampleID"] = project_map["sampleID"].str.replace("[^a-zA-Z0-9_\.]", "_", regex=True)
	project_map["sampleID"] = make_unique(project_map["sampleID"])


	print("Final counts:")
	counts = (project_map["WELL"]
       	.isin(CONFIG["CTRL_WELLS"])
		.value_counts())
	counts.index = counts.index.map({False: "Client samples", True: "MSL Controls"})
	for index in counts.index:
		print(str(index) + "\t" + str(counts.loc[index]))
	print()

	project_map = project_map[['RUN.PLATEPOSITION', "sampleID"]]
	
	project_map.to_csv(args.output, sep="\t", index=False)



def apply_manifest(sample_to_gdna, manifest_file):
	sample_to_gdna = sample_to_gdna.copy()
	manifest = get_manifest(manifest_file)["sample_name"]

	# temporarily detach wells where controls are expected
	ctrls_bool = sample_to_gdna["WELL"].isin(CONFIG["CTRL_WELLS"])
	ctrls = sample_to_gdna.loc[ctrls_bool, :]
	sample_to_gdna = sample_to_gdna.loc[~ctrls_bool, :]

	not_in_gDNA = ~manifest.index.isin(sample_to_gdna["SAMPLE ID"].values)
	replacements = sample_to_gdna.loc[:, "SAMPLE ID"].map(manifest) # NaN if not in manifest
	to_repl = ~replacements.isna()
	sample_to_gdna.loc[to_repl, "SAMPLE ID"] = replacements[to_repl]
	if(any(not_in_gDNA)):
		print("Error: Manifest samples missing from gDNA maps:")
		print(manifest.loc[list(not_in_gDNA)])
		quit()
	if(not all(to_repl)):
		print("Warning: gDNA samples not in manifest:")
		sample_to_gdna.loc[~to_repl, "SAMPLE ID"].apply(print)
		print("(Controls are expected in " + ",".join(CONFIG["CTRL_WELLS"]) + ")")
		print("This is okay if these samples were added by MD Genomics.\n")

	sample_to_gdna = pd.concat([sample_to_gdna, ctrls]) # add back the controls
	return sample_to_gdna

def get_pooling(indir):
	# load config from json??
	

	pooling_dfs = []
	# process files in indir
	for pooling_file in indir.iterdir():
		if pooling_file.suffix == ".xlsx":
			# get gDNA plate IDs, IDT plate numbers, and run name
			
			# "124 Ill_2step_XTR_16S PCR MSL_NS_124 Pooling detail 1-11-2023.xlsx"
			# "/local/projects-t3/MSL/pipelines/test/pooling_detail/118 Ill_2step XT_16S PCR MSL_MS_118 Pooling detail 10-20-2022.xlsx"
			pooling_df = pd.read_excel(pooling_file, skiprows=1, header=None, dtype=str)
			# could also get run name from the "Illumina 16S PCR pool ID" column (column 11)
			run = pooling_df.iloc[0, 4]
			run = CONFIG["RUN_FORMATTER"](run)

			# remove unneeded rows and one empty column
			pooling_df = pooling_df.iloc[3:, 1:]
			# get header
			pooling_df.columns = pooling_df.iloc[0, :]
			pooling_df = pooling_df.iloc[1:, :]
			
			pooling_df = (pooling_df.iloc[np.arange(int(pooling_df.shape[0] / 2)) * 2, :] # remove odd rows
				.loc[-pd.isna(pooling_df["gDNA plate ID"]), :]) # remove rows with empty field

			# get col immediately left of this column
			maybe_idt_plates = pooling_df.iloc[:, pooling_df.columns.get_loc("Sample Description") - 1]
			
			# does every element start with "XTR_Plate"?
			if all([str(x).startswith("XTR_Plate") for x in maybe_idt_plates]):		
				# process as IDT UDI plates	
				pooling_df["RUN.PLATE"] = pooling_df.apply(CONFIG["RUN_PLATE_FORMATTERS"]["IDT for Illumina"], axis=1, run=run) 
				pooling_df = pooling_df[["RUN.PLATE", "gDNA plate ID"]]
				
			elif (pooling_df.columns[1] == "XT barcode plate" 
				and all([re.search("^set [A-D]$", plate) for plate in pooling_df["XT barcode plate"]])):
					# process as Nextera XT UDI plates
					pooling_df["RUN.PLATE"] = pooling_df.apply(CONFIG["RUN_PLATE_FORMATTERS"]["Nextera XT"], axis=1, run=run) 
					pooling_df = pooling_df[["RUN.PLATE", "gDNA plate ID"]]
			else:
				msg = "Error: Unable to detect UDI plate configuration from pooling: " + pooling_file
				msg = msg + """
					Acceptable configurations:
						IDT for Illumina: The column left of "Sample Description" has values that match the pattern "^XTR_Plate.*"
						Nextera XT: The third column is called "XT barcode plate" and has values that match the pattern "^set [A-D]$"
				"""
				quit(msg)

		pooling_dfs.append(pooling_df)

	return pd.concat(pooling_dfs)

def make_unique(x, sep=".", wrap_in_brackets=False):
	if not all([type(val) is str for val in x]):
		print(x)
		quit("Input to `make_unique` must be a string iterable. Got " )
	x = pd.Series(x, dtype=str)
	
	def fun(a):
		if len(a) > 1:
			suffixes = range(len(a))
			if wrap_in_brackets:
				suffixes = ["(" + str(s) + ")" for s in suffixes]
			return [val + sep + str(suffix) for val, suffix in zip(a, suffixes)]
		else:
			return a
	deduplicated = x.groupby(x).transform(fun)
	
	# probably terribly unoptimized
	values_still_duplicated = deduplicated.loc[deduplicated.duplicated()]
	if len(pd.Series(values_still_duplicated).dropna()) > 0:
		quit("""
make_unique failed to make iterable unique.\n
This is because appending '<dup_number>' to duplicate values led to
creation of term(s) that were in the original dataset: \n[
""" + ", ".join(values_still_duplicated) + """
]\n\nPlease try again with a different argument for either `wrap_in_brackets` or `sep`
			"""
			)
		
	return deduplicated.array



def get_gdna(indir):
	gdna_dfs = []

	for gdna_file in indir.iterdir():
		if gdna_file.suffix == ".xlsx":
			gdna_df = pd.read_excel(gdna_file, header=0)
			if(not all([x in gdna_df.columns for x in ["WELL", "SAMPLE ID"]])):
				raise ValueError("gDNA map does not contain header columns \"WELL\" and \"SAMPLE ID\": " + gdna_file.name)
			
			gdna_plate = str(gdna_file.name)
			gdna_plate_parts = gdna_plate.split("_")
			gdna_plate = gdna_plate_parts[len(gdna_plate_parts)-1].split(".")[0]
			gdna_df["gDNA plate ID"] = gdna_plate
			gdna_dfs.append(gdna_df)
			
			if args.debug:
				print(str(gdna_file) + ":")
				print(gdna_df)
	gdna_df = pd.concat(gdna_dfs)


	return gdna_df

def get_manifest(filepath):
	manifest_df = pd.read_excel(filepath, skiprows=14, header=0, index_col=1)

	return manifest_df

if __name__=="__main__":
	ap = argparse.ArgumentParser(
		epilog='''Defaults: --gdna gDNA_maps/ --pooling pooling/ --output ./project_map.txt'''
	)
	#-d $models -i $fasta -o 
	ap.add_argument("--pooling", "-p", type=Path, required=False, default=Path("./pooling"))
	ap.add_argument("--gdna", "-g", type=Path, required=False, default=Path("./gDNA_maps"), help="Directory containing gDNA plate maps named \"*_<id>.xlsx\"")
	ap.add_argument("--manifest", "-m", type=Path, required=False)
	ap.add_argument("--control", type=str, required=False)
	ap.add_argument("--debug", action="store_true")

	ap.add_argument("--output", "-o", type=Path, required=False, default=Path("./project_map.txt"))

	args = ap.parse_args()
	main(args)
