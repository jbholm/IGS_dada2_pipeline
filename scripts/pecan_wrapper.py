import os, subprocess, argparse, csv, json
from pathlib import Path

def run_commands(cmds, verbose=True):
	if verbose:
		print("Running external command line application(s). This may print "
			  "messages to stdout and/or stderr.")
		print("The command(s) being run are below. These commands cannot "
			  "be manually re-run as they will depend on temporary files that "
			  "no longer exist.")
	for cmd in cmds:
		if verbose:
			print("\nCommand:", end=' ')
			print(" ".join(cmd), end='\n\n')
		subprocess.run(cmd, check=True)

pipeline_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
with (Path(pipeline_dir) / Path("config.json")).open("r") as fh:
    CONFIG = json.load(fh)

_SUBS = {
	"Paenibacillus_polymyxa_Paenibacillus_peoriae_Paenibacillus_borealis_Paenibacillus_durus_Paenibacillus_barcinonensis_Paenibacillus_chibensis_Paenibacillus_glacialis_Paenibacillus_lautus_Paenibacillus_amylolyticus_Paenibacillus_graminis_Paenibacillus_agarexedens_Paenibacillus_lactis_Paenibacillus_nanensis_Paenibacillus_cineris_Paenibacillus_xylanexedens_Paenibacillus_cookii": "Paenibacillus_Cluster_1",
	"Bacillus_subtilis_Bacillus_vallismortis_Bacillus_mojavensis_Bacillus_amyloliquefaciens_Bacillus_tequilensis_Bacillus_licheniformis_Bacillus_sonorensis_Bacillus_velezensis": "Bacillus_Cluster_1",
	"Pseudomonas_reactans_Pseudomonas_fluorescens_Pseudomonas_salomonii_Pseudomonas_moraviensis_Pseudomonas_veronii_Pseudomonas_arsenicoxydans_Pseudomonas_koreensis_Pseudomonas_mandelii_Pseudomonas_marginalis_Pseudomonas_gessardii_Pseudomonas_jessenii": "Pseudomonas_Cluster_1"
}
_VAGINAL_SUBS = {
	"Lactobacillus_crispatus_Lactobacillus_helveticus": "Lactobacillus_crispatus",
	"Lactobacillus_gasseri_Lactobacillus_johnsonii": "Lactobacillus_gasseri",
	"Lactobacillus_gasseri_johnsonii": "Lactobacillus_gasseri",
	"Proteobacteria_bacterium": "Candidatus_Saccharibacteria_genomosp_TM7_H1",
	"Lactobacillus_kitasatonis": "Lactobacillus_crispatus",
	"Lactobacillus_acidophilus": "Lactobacillus_crispatus",
	"Lactobacillus_jensenii_1": "Lactobacillus_jensenii",
	"Lactobacillus_jensenii_3": "Lactobacillus_jensenii",
	"f_Beggiatoaceae": "Candidatus_Saccharibacteria_genomosp_TM7_H1",
	"Lactobacillus_kitasatonis_Lactobacillus_gallinarum_Lactobacillus_crispatus": "Lactobacillus_crispatus",
	"Atopobium_sp_Atopobium_parvulum_Atopobium_rimae_Atopobium_vaginae": "Atopobium_vaginae",
	"Lactobacillus_acidophilus_Lactobacillus_sp_Lactobacillus_taiwanensis_Lactobacillus_johnsonii": "Lactobacillus_johnsonii",
	"Lactobacillus_crispatus_helveticus": "Lactobacillus_crispatus"
}

def postprocess(taxonomy_file, vaginal):
	assignments = []
	
	subs = _SUBS
	if vaginal:
		subs = subs | _VAGINAL_SUBS
	
	out_filename = Path(taxonomy_file).stem + ".clean.txt"
	
	with open(out_filename, "w+") as file_out:
		with open(taxonomy_file, 'r') as file_in:
			writer = csv.writer(file_out, delimiter='\t')
			for in_line in file_in:
				row_values = in_line.strip().split("\t")
				print(row_values[0])
				print(row_values[1])
				if row_values[1] in subs.keys():
					row_values[1] = subs[row_values[1]]
				writer.writerow(row_values)
	
	return

def main(models, in_file, out_file, vaginal):
	cmd = [CONFIG['speciateIT_bin'],
			'-d', str(models),
			'-i', str(in_file),
			'-o', str(out_file)
	]
	try:
		run_commands([cmd])
		postprocess(Path(out_file) / Path("MC_order7_results.txt"), args.vaginal)
	except subprocess.CalledProcessError as e:
		raise Exception("An error was encountered while running PECAN"
						" (return code %d), please inspect stdout"
						" and stderr to learn more." % e.returncode)
	return

if __name__=="__main__":
	ap = argparse.ArgumentParser()
	#-d $models -i $fasta -o 
	ap.add_argument("-d", type=Path, default=CONFIG['pecan_models'])
	ap.add_argument("-i", type=Path)
	ap.add_argument("-o", type=Path)

	ap.add_argument("--vaginal", action="store_true")

	args = ap.parse_args()
	main(args.d, args.i, args.o, args.vaginal)
