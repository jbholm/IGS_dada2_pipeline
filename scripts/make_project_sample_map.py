#!/usr/bin/env python3

import subprocess, argparse, os, sys, re
import numpy as np
import pandas as pd
from pathlib import Path

CONFIG = {
	"RUN_PLATE_FORMATTERS": {
		"IDT for Illumina": lambda row, run: ".".join([run, "IDT" + row["primer plate name"].split("XTR_Plate ")[1]]),
		"Nextera XT": lambda row, run: ".".join([run, "UDI" + row["XT barcode plate"].split("set ")[1]])
	},
	"RUN_FORMATTER": lambda x: re.sub(pattern="[^.a-zA-Z0-9]", repl=".", string=x)
}

def main(args):
	gdna_to_barcode_plate = get_pooling(args.pooling)

	# main columns: ["gDNA plate ID", "WELL", "SAMPLE ID"]
	sample_to_gdna = get_gdna(args.gdna)
	sample_to_gdna = sample_to_gdna.loc[-pd.isna(sample_to_gdna["SAMPLE ID"]), ]

	project_map = gdna_to_barcode_plate.merge(sample_to_gdna, how="inner")
	join_fun = lambda row: ".".join([row["RUN.PLATE"], row["WELL"]])
	project_map["RUN.PLATEPOSITION"] = project_map.apply(join_fun, axis=1)
	project_map = project_map.rename(columns={"SAMPLE ID": "sampleID"})
	project_map = project_map[['RUN.PLATEPOSITION', "sampleID"]]

	if(args.manifest):
		project_map = project_map.merge(get_manifest(args.manifest))

	project_map.to_csv(args.output, sep="\t", index=False)

def get_pooling(indir):
	# load config from json??
	

	pooling_dfs = []
	# process files in indir
	for pooling_file in indir.iterdir():
		# get gDNA plate IDs, IDT plate numbers, and run name
		
		# "124 Ill_2step_XTR_16S PCR MSL_NS_124 Pooling detail 1-11-2023.xlsx"
		# "/local/projects-t3/MSL/pipelines/test/pooling_detail/118 Ill_2step XT_16S PCR MSL_MS_118 Pooling detail 10-20-2022.xlsx"
		pooling_df = pd.read_excel(pooling_file, skiprows=1, header=None, dtype=str)
		run = pooling_df.iloc[0, 4]
		run = CONFIG["RUN_FORMATTER"](run)

		# remove unneeded rows and one empty column
		pooling_df = pooling_df.iloc[3:, 1:]

		pooling_df.columns = pooling_df.iloc[0, :]
		pooling_df = pooling_df.iloc[1:, :]
		
		# remove odd rows that occur because of merged rows
		# any way to do this with just [] not iloc?
		pooling_df = pooling_df.iloc[np.arange(int(pooling_df.shape[0] / 2)) * 2, :]

		# remove rows we aren't interested in
		pooling_df = pooling_df.loc[-pd.isna(pooling_df["gDNA plate ID"]), :]

		# get col immediately left of this column
		maybe_idt_plates = pooling_df.iloc[:, pooling_df.columns.get_loc("Sample Description") - 1]
		
		# does every element start with "XTR_Plate"?
		if all([str(x).startswith("XTR_Plate") for x in maybe_idt_plates]):		
			# process as IDT UDI plates	
			pooling_df["RUN.PLATE"] = pooling_df.apply(CONFIG["RUN_PLATE_FORMATTERS"]["IDT for Illumina"], axis=1, run=run) 
			pooling_df = pooling_df[["RUN.PLATE", "gDNA plate ID"]]
			
		elif pooling_df.columns[1] == "XT barcode plate" :
			if all([re.search("^set [A-D]$", plate) for plate in pooling_df["XT barcode plate"]]):
				# process as XT UDI plates
				pooling_df["RUN.PLATE"] = pooling_df.apply(CONFIG["RUN_PLATE_FORMATTERS"]["Nextera XT"], axis=1, run=run) 
				pooling_df = pooling_df[["RUN.PLATE", "gDNA plate ID"]]

		pooling_dfs.append(pooling_df)

	return pd.concat(pooling_dfs)

def get_gdna(indir):
	gdna_dfs = []

	for gdna_file in indir.iterdir():
		# get gDNA plate IDs, IDT plate numbers, and run name
		
		# "124 Ill_2step_XTR_16S PCR MSL_NS_124 Pooling detail 1-11-2023.xlsx"
		# "/local/projects-t3/MSL/pipelines/test/pooling_detail/118 Ill_2step XT_16S PCR MSL_MS_118 Pooling detail 10-20-2022.xlsx"
		gdna_df = pd.read_excel(gdna_file, header=0)
		gdna_plate = str(gdna_file.name)
		gdna_plate_parts = gdna_plate.split("_")
		gdna_plate = gdna_plate_parts[len(gdna_plate_parts)-1].split(".")[0]
		gdna_df["gDNA plate ID"] = gdna_plate
		gdna_dfs.append(gdna_df)
	
	return pd.concat(gdna_dfs)

def get_manifest(filepath):
	return pd.DataFrame()

if __name__=="__main__":
	ap = argparse.ArgumentParser()
	#-d $models -i $fasta -o 
	ap.add_argument("--pooling", "-p", type=Path, required=False, default=Path("./pooling"))
	ap.add_argument("--gdna", "-g", type=Path, required=False, default=Path("./gDNA_maps"))
	ap.add_argument("--manifest", "-m", type=Path, required=False)

	ap.add_argument("--output", "-o", type=Path, required=False, default=Path("./project_map.txt"))

	args = ap.parse_args()
	main(args)
