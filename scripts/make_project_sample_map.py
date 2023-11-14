#!/usr/bin/env python3

import subprocess, argparse, os, sys, re
import numpy as np
import pandas as pd
from pathlib import Path
from enum import Enum

CONFIG = {
	"RUN_PLATE_FORMATTERS": {
		"IDT for Illumina": lambda row, run: ".".join([run, "IDT" + row["primer plate name"].split("XTR_Plate ")[1]]),
		"Nextera XT": lambda row, run: ".".join([run, row["XT barcode plate"].replace("set ", "UDI")])
	},
	"RUN_FORMATTER": lambda x: re.sub(pattern="[^.a-zA-Z0-9]", repl=".", string=x),
	"CTRL_WELLS": {
		"E12": "Zymo Positive Library Control", 
		"F12": "Vaginal Positive Library Control", 
		"G12": "Negative Library Control", 
		"H12": "Negative Library Control"
		},
	"gDNA_COLUMNS": ["WELL_ID", "SAMPLE_NAME"]
}

Manifest_type = Enum('Manifest_type', ['TUBES', 'PLATES'])

def main(args):
	if args.verbose:
		pd.set_option('display.max_columns', None)
		pd.set_option('display.max_rows', None)
		pd.set_option('display.width', 2000)
		pd.set_option('display.max_colwidth', None)
	
	# main columns: ["PLATE ID, "WELL", "SAMPLE ID"]
	if args.platform == "illumina":
		sample_to_plate_wells = get_gdna(args.gdna)
		run_plate_name = "PLATE ID"

		if args.verbose:
			print("gDNA-to-barcode plates:")
			print(sample_to_plate_wells)
			print("")

		# Pooling will have "RUN.PLATE", "plate ID", "run"
		gdna_to_run_barcode_plate = get_pooling(args.pooling, args.platform)
		if args.verbose:
			print("Pooling detail:")
			print(gdna_to_run_barcode_plate)
			print("")
		project_map = gdna_to_run_barcode_plate.merge(sample_to_plate_wells, how="inner", left_on="gDNA plate ID", right_on="PLATE ID")
		if sample_to_plate_wells.shape[0] > project_map.shape[0]:
			print(sample_to_plate_wells["PLATE ID"].dtype)
			plate_diff = set(sample_to_plate_wells["PLATE ID"]).difference(
				set(gdna_to_run_barcode_plate["gDNA plate ID"])
				)
			print([type(x) for x in plate_diff])
			print(plate_diff)
			print("gDNA plates not in pooling details: " + ",".join(plate_diff))
			quit("Error: Lost rows when joining gDNA maps to pooling detail. gDNA plate IDs did not match.")
		if args.verbose: 
			print("Joined map:")
			print(project_map)
			print("")

	else:
		# pacbio demux map already links the run and plate positions to the tube ID
		# main columns: "RUN ID", "well", "sample ID" (sample ID is ignored if client manifest provided and in plate format)
		# in doing this, the run name will be picked up from the filename, and the plate position columns concatenated
		if args.pooling and args.pooling.exists():
			print("PacBio demux maps in use. Ignoring pooling file.\n")

		project_map = get_pacbio_demux_maps(args.pb_demux_maps)
		project_map["RUN.PLATE"] = project_map["RUN"] # artifact of conforming PacBio to Illumina pipeline

	# if args.assigned_labels:
	# 	project_map = apply_assigned_labels(
	# 		project_map, 
	# 		args.assigned_labels
	# 		)

	if args.manifest:
		validate_manifest(args)
		project_map = apply_manifest(
			project_map, 
			manifest_file = args.manifest, 
			manifest_type = args.manifest_format,
			tubeid_col = args.manifest_tubeid_column, 
			sampleid_col = args.manifest_samplename_column,
			plateid_col = args.manifest_plateid_column,
			plate_row_column = args.manifest_platerow_column,
			plate_remap = args.manifest_plate_remap
			)
		if args.verbose: 
			print("After applying manifest:")
			print(project_map)
			print("")

	project_map.loc[:, "RUN.PLATEPOSITION"] = project_map.apply(
		lambda row: ".".join([row["RUN.PLATE"], row["WELL"]]), 
		axis=1
		)

	project_map = project_map.rename(columns={"SAMPLE ID": "sampleID"})

	# Make samples unique within each run, disambiguating by gDNA plate name
	fun = lambda df: disambiguate(df, ambig_col="sampleID", disambig_col="RUN.PLATE")
	project_map = project_map.groupby("RUN", as_index=False).apply(fun)
	
	# At this point, if any samples still have the same name, disambiguate by the run name
	# for PacBio, there will never be any ambiguous samples at this stage
	project_map = disambiguate(project_map, ambig_col="sampleID", disambig_col="RUN")

	project_map["sampleID"] = project_map["sampleID"].str.replace("[^a-zA-Z0-9_\.]", "_", regex=True)

	print("Final counts:")
	counts = (project_map["WELL"]
       	.isin(CONFIG["CTRL_WELLS"].keys())
		.value_counts())
	counts.index = counts.index.map({False: "Client samples", True: "MSL Controls"})
	for index in counts.index:
		print(str(index) + "\t" + str(counts.loc[index]))
	print()

	project_map = project_map[['RUN.PLATEPOSITION', "sampleID"]]
	
	project_map.to_csv(args.output, sep="\t", index=False)

	return


def get_gdna(indir):
	gdna_dfs = []

	for gdna_file in indir.iterdir():
		if gdna_file.suffix == ".xlsx":
			gdna_df = (pd.read_excel(gdna_file, header=0, dtype=str)
				.rename(str.upper, axis='columns'))
			if not all([x in gdna_df.columns for x in CONFIG["gDNA_COLUMNS"]]):
				raise ValueError("gDNA map does not contain header columns [" + ", ".join(CONFIG["gDNA_COLUMNS"]) + "]: " + gdna_file.name)			
			
			print(str(gdna_file) + ":")
			if args.verbose:
				print(gdna_df)
				print("")

			gdna_df["PLATE ID"] = str(gdna_file.name).split("-")[0]
			gdna_df = gdna_df.rename(
				columns={
					CONFIG["gDNA_COLUMNS"][0]: "WELL", 
					CONFIG["gDNA_COLUMNS"][1]: "SAMPLE ID"
					}
					)

			# log stats for this plate
			blanks = pd.isna(gdna_df["SAMPLE ID"])
			if True in blanks.value_counts().index:
				print("Removing " + str(blanks.value_counts()[True]) + " wells with no sample. If this is too many, check that the 'SAMPLE ID' column is filled.")
			gdna_df = gdna_df.loc[~blanks, ]

			waters = gdna_df["SAMPLE ID"] == "water"
			print("Removing " + str(waters.sum()) + " wells containing water.")
			gdna_df = gdna_df.loc[~waters, ]

			counts = (gdna_df["WELL"]
				.isin(CONFIG["CTRL_WELLS"].keys())
				.value_counts())
			counts.index = counts.index.map({False: "Client samples", True: "MSL Controls"})
			for index in counts.index:
				print(str(index) + "\t" + str(counts.loc[index]))
			if "MSL Controls" not in counts.index:
				print("Adding MSL controls...")
				ctrl_df = pd.DataFrame.from_dict(
					CONFIG["CTRL_WELLS"], 
					orient='index', 
					columns=['SAMPLE ID']
					)
				ctrl_df["WELL"] = ctrl_df.index
				ctrl_df["PLATE ID"] = gdna_df["PLATE ID"][0]
				gdna_df = pd.concat([gdna_df, ctrl_df])
			print()


			# needed because gDNA maps have MSL control names in the SAMPLE ID column
			# (and sometimes they have client sample names too)
			gdna_df["SAMPLE ID"] = gdna_df["SAMPLE ID"].apply(str).str.replace("[^a-zA-Z0-9_\.]", "_", regex=True)
			gdna_df["SAMPLE ID"] = make_unique(gdna_df["SAMPLE ID"])
			
			gdna_dfs.append(gdna_df)

	gdna_df = pd.concat(gdna_dfs)

	return gdna_df

def get_pacbio_demux_maps(indir):
	# main columns: "run ID", "well", "sample ID"
	demux_maps = []

	for map_file in indir.iterdir():
		if map_file.suffix == ".txt":
			try:
				df = pd.read_table(map_file, header=0, dtype=str)
			except UnicodeDecodeError as e:
				print(str(e))
				print("Please try running dos2unix on the pacbio demux maps.")
				quit()
		elif map_file.suffix == ".csv":
			df = pd.read_csv(map_file, header=0, dtype=str)
		else:
			pass
		
		df = df.rename(str.upper, axis='columns')
		if(not all([x in df.columns for x in ["ROW", "WELL", "SAMPLEID"]])):
			raise ValueError("Demux map does not contain header columns \"ROW\", \"WELL\", and \"sampleID\": " + map_file.name)
		df = df.rename(columns={"SAMPLEID": "SAMPLE ID"})
		
		run = str(map_file.name).split("-")[0]
		df["RUN"] = run

		df["WELL"] = df.apply(lambda row: row["ROW"] + row["WELL"], axis=1)

		if args.verbose:
			print(str(map_file) + ":")
			print(df)
			print("")

		counts = (df["WELL"]
			.isin(CONFIG["CTRL_WELLS"].keys())
			.value_counts())
		counts.index = counts.index.map({False: "Client samples", True: "MSL Controls"})
		for index in counts.index:
			print(str(index) + "\t" + str(counts.loc[index]))
		if "MSL Controls" not in counts.index:
			print("Adding MSL controls...")
			ctrl_df = pd.DataFrame.from_dict(
				CONFIG["CTRL_WELLS"], 
				orient='index', 
				columns=['SAMPLE ID']
				)
			ctrl_df["WELL"] = ctrl_df.index
			ctrl_df["RUN"] = run
			df = pd.concat([df, ctrl_df])
		print()
		
		demux_maps.append(df)

	demux_df = pd.concat(demux_maps)

	return demux_df

# def try_override_ctrls(gDNA_map, pacbio_demux_map):
# 	# demux_maps have run name but gDNA_maps have plate id. what to do? Wait until after merging with pooling?


# def apply_assigned_labels(sample_map, label_file):
# 	label_df = pd.read_excel(label_file, header=0, index_col=None, dtype=str)
# 	label_df = label_df.set_index(label_df.columns[1])
# 	if args.verbose: 
# 		print("Labels assigned to client tubes:")
# 		print(label_df)
# 		print("")
# 	replacements = sample_map["SAMPLE ID"].map(label_df.iloc[:, 0])
# 	if not all(replacements.isna()):
# 		to_repl = ~replacements.isna()
# 		sample_map.loc[to_repl, "SAMPLE ID"] = replacements[to_repl]
# 	else:
# 		print("Failed to assign client sample names from the second column of \"assigned labels\" spreadsheet. Trying inverse column order...")
# 		label_df = pd.read_excel(label_file, header=0, index_col=None, dtype=str)
# 		label_df = label_df.set_index(label_df.columns[0])
# 		if args.verbose: 
# 			print("Labels assigned to client tubes:")
# 			print(label_df)
# 			print("")
# 		replacements = sample_map["SAMPLE ID"].map(label_df.iloc[:, 0])
# 		to_repl = ~replacements.isna()
# 		sample_map.loc[to_repl, "SAMPLE ID"] = replacements[to_repl]
# 	if args.verbose: 
# 		print("With client sample names:")
# 		print(sample_map)
# 		print("")
# 	return(sample_map)

def validate_manifest(args):
	manifest = pd.read_excel(args.manifest, skiprows=14, header=0, index_col=None, dtype=str)
	# validate command line options (this isn't helpful if defaults are being used though)
	if args.manifest_format == Manifest_type.TUBES:
		for argname in ["manifest_samplename_column", "manifest_tubeid_column"]:
			if getattr(args, argname) not in manifest.columns:
				raise KeyError("\"%s\" (which may be given by --%s) not in %s" \
		   			% (getattr(args, argname), argname, str(args.manifest)))
	elif args.manifest_format == Manifest_type.PLATES:
		for argname in ["manifest_samplename_column", "manifest_plateid_column", "manifest_platerow_column"]:
			if getattr(args, argname) not in manifest.columns:
				raise KeyError("\"%s\" (which may be given by --%s) not in %s" \
		   			% (getattr(args, argname), argname, str(args.manifest)))

def apply_manifest(sample_to_gdna, manifest_file, manifest_type, tubeid_col, sampleid_col, plateid_col, plate_row_column, plate_remap):

	sample_to_gdna = sample_to_gdna.copy() # idk why
	run_plate_name = "PLATE ID" if "PLATE ID" in sample_to_gdna.columns else "RUN"
	
	manifest = pd.read_excel(manifest_file, skiprows=14, header=0, index_col=None, dtype=str)
	# remove empty cols
	manifest = manifest.loc[:, [any(-pd.isna(manifest[col])) for col in manifest.columns]]

	# temporarily detach wells where controls are expected
	ctrls_bool = sample_to_gdna["WELL"].isin(CONFIG["CTRL_WELLS"].keys())
	ctrls = sample_to_gdna.loc[ctrls_bool, :]
	sample_to_gdna = sample_to_gdna.loc[~ctrls_bool, :]

	if manifest_type == Manifest_type.TUBES:
		manifest = manifest.set_index(tubeid_col)

		left_on = "SAMPLE ID"
	elif manifest_type == Manifest_type.PLATES:
		plate_remap_series = pd.Series(dtype=str)
		if plate_remap:
			for pair in plate_remap:
				msl_id, client_id = pair.split("=", 1)
				plate_remap_series.loc[client_id] = msl_id
		
		replacements = manifest[plateid_col].map(plate_remap_series)
		to_repl = ~replacements.isna()
		manifest.loc[to_repl, plateid_col] = replacements[to_repl]
		row_col = manifest.columns.get_loc(plate_row_column)
		manifest["WELL"] = manifest.apply(lambda row: row.iloc[row_col] + row.iloc[row_col + 1], axis=1)
		manifest = manifest.set_index([plateid_col, "WELL"])

		left_on = [run_plate_name, "WELL"]
	
	if args.verbose:
		print(str(manifest_file) + ":")
		print(manifest)
		print("")
	
	not_in_gDNA = manifest.index.difference(pd.Index(sample_to_gdna[left_on]))
	if any(not_in_gDNA):
		print("Error: Manifest samples missing from gDNA maps:")
		print(manifest.loc[not_in_gDNA])
		quit()

	# merge(how="left"...) yields NaN in manifest columns if not in gDNA plates
	merged = sample_to_gdna.merge(manifest, how="left", left_on=left_on, right_index=True, suffixes=(None, "_y"))
	if sampleid_col + "_y" in merged.columns:
		sampleid_col = sampleid_col + "_y"
	
	
	# replace the names of samples found in the manifest; leave the rest untouched
	to_repl = ~merged[sampleid_col].isna()
	merged.loc[to_repl, "SAMPLE ID"] = merged.loc[to_repl, sampleid_col]
	
	merged = merged[sample_to_gdna.columns]

	if not all(to_repl):
		print("Warning: gDNA samples not in manifest:")
		merged.loc[~to_repl, "SAMPLE ID"].apply(print)
		print("(Controls are expected in " + ",".join(CONFIG["CTRL_WELLS"].keys()) + ")")
		print("This is okay if these samples were added by MD Genomics.\n")

	merged = pd.concat([merged, ctrls]) # add back the controls
	return merged

def get_pooling(indir, platform):
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
			
			pooling_df["RUN"] = run

			if args.verbose:
				print(str(pooling_file) + ":")
				print(pooling_df)
				print("")

			if platform == "illumina":
				# get col immediately left of "Sample Description" column
				maybe_idt_plates = pooling_df.iloc[:, pooling_df.columns.get_loc("Sample Description") - 1]
				
				# does every element start with "XTR_Plate"?
				if all([str(x).startswith("XTR_Plate") for x in maybe_idt_plates]):		
					# process as IDT UDI plates	
					pooling_df["RUN.PLATE"] = pooling_df.apply(CONFIG["RUN_PLATE_FORMATTERS"]["IDT for Illumina"], axis=1, run=run) 
					
				elif (pooling_df.columns[1] == "XT barcode plate" 
					and all([re.search("^set [A-D]$", plate) for plate in pooling_df["XT barcode plate"]])):
						# process as Nextera XT UDI plates
						pooling_df["RUN.PLATE"] = pooling_df.apply(CONFIG["RUN_PLATE_FORMATTERS"]["Nextera XT"], axis=1, run=run) 
						
				else:
					msg = "Error: Unable to detect UDI plate configuration from pooling: " + pooling_file
					msg = msg + """
						Acceptable configurations:
							IDT for Illumina: The column left of "Sample Description" has values that match the pattern "^XTR_Plate.*"
							Nextera XT: The third column is called "XT barcode plate" and has values that match the pattern "^set [A-D]$"
					"""
					quit(msg)
			elif platform == "pacbio":
				pooling_df["RUN.PLATE"] = pooling_df["RUN"]

			pooling_df = pooling_df.loc[:, ["RUN.PLATE", "gDNA plate ID", "RUN"]]

			pooling_dfs.append(pooling_df)

	return pd.concat(pooling_dfs)

def make_names(x):
	# A syntactically valid name consists of letters, numbers and the dot or underline 
	# characters and starts with a letter or the dot not followed by a number.
	x = re.sub(r'[a-zA-Z0-9\._]', r'\.', x)

	# I'm going to use R's make.names() strategy of prepending one or two .'s because
	# this is the easiest to undo, on the client's side
	x = re.sub(r'^(.[0-9])', r'\.\1', x)
	x = re.sub(r'^(.[0-9])', r'\.\1', x)
	return x

def make_unique(x, sep=".", wrap_in_brackets=False):
	if not (all([type(val) is str for val in x]) or type(x) is pd.Series):
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

def disambiguate(df, ambig_col, disambig_col, sep = "."):
	def fun(a):
		if a.shape[0] > 1:
			suffixes = a[disambig_col]
			a[ambig_col] = [val + sep + str(suffix) for val, suffix in zip(a[ambig_col], suffixes)]
			return a
		else:
			return a
	
	return df.groupby(ambig_col, as_index=False).apply(fun)

if __name__=="__main__":
	ap = argparse.ArgumentParser(
		epilog='''Defaults: --gdna gDNA_maps/ --pooling pooling/ --output ./project_map.txt'''
	)
	#-d $models -i $fasta -o 
	ap.add_argument("--pooling", "-p", type=Path, required=False, default=Path("./pooling"), help="Directory containing pooling details. Default: ./pooling/")
	ap.add_argument("--gdna", "-g", type=Path, required=False, default=Path("./gDNA_maps"), help="Directory containing gDNA plate maps named \"*_<id>.xlsx\". Default: ./gDNA_maps")
	ap.add_argument("--pb_demux_maps", type=Path, required=False, default=Path("./demux_maps/"), help="Directory containing demux maps named \"<run>_*.csv\". Default: ./demux_maps/")
	ap.add_argument("--manifest", "-m", type=Path, required=False)
	ap.add_argument("--manifest_samplename_column", type=str, required=False, default="sample_name", help="The name of the manifest column giving the client's sample names. Default: sample_name")
	ap.add_argument("--manifest_format", required=False, choices=Manifest_type.__members__.keys(), help="A tube-formatted manifest must have a tubeid column. A plate-formatted manifest has adjacent row and column position columns.")
	ap.add_argument("--manifest_tubeid_column", required=False, default="tube_label", type=str, help="For a tube-formatted manifest, the name of the column that gives the tube ID. Default: tube_label")
	ap.add_argument("--manifest_plateid_column", type=str, required=False, default="plate id", help="The name of the column giving client plate IDs. Default: \"plate id\"")
	ap.add_argument("--manifest_plate_remap", required=False, action="append", help="Use this to map MSL gDNA plate or runs to client manifest plate IDs. Format: --manifest_plate_remap=\"<MSL_ID>=<client_id>\". Can be repeated for multiple plates.")
	ap.add_argument("--manifest_platerow_column", required=False, default="row position", type=str, help="For a plate-formatted manifest, the name of the column that gives the well row positions. The well column positions will be taken from the immediate next column. Default: \"row position\"")
	# ap.add_argument("--assigned_labels", required=False, type=Path, help="The \"ASSIGNED LABELS\" spreadsheet, if MD Genomics assigned tube barcodes.")
	ap.add_argument("--control", type=str, required=False)
	ap.add_argument("--verbose", action="store_true")
	ap.add_argument("--platform", choices=["pacbio", "illumina"], required=True)

	ap.add_argument("--output", "-o", type=Path, required=False, default=Path("./project_map.txt"))

	args = ap.parse_args()
	if args.manifest:
		if not args.manifest_format:
			if re.search("plates.xlsx$", str.lower(args.manifest.name)):
				args.manifest_format = Manifest_type.PLATES
				print("Inferred manifest type was PLATES.")
			elif re.search("tubes.xlsx$", str.lower(args.manifest.name)):
				args.manifest_format = Manifest_type.TUBES
				print("Inferred manifest type was TUBES.")
			else:
				quit("Need --manifest_format")
		else:
			args.manifest_format = Manifest_type[args.manifest_format]
	elif (args.manifest_samplename_column or args.manifest_format or 
      		args.manifest_tubeid_column or args.manifest_plateid_column or
		      args.manifest_plate_remap or args.manifest_platerow_column
		):
			print("--manifest not given, so ignoring all other --manifest* options.")
			
	main(args)
