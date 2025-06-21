# Parse ESS file to generate doppelganger driver files
import argparse
import glob
import os

ASSESSMENT_LEVELS = {
	"pass" : 0,
	"questionable" : 1,
	"questionable_critical" : 2,
	"fail" : 3,
	"dnu" : 4,
}

def read_ESS(ESS, assessment_header_keys, threshold, doppelganger_n=2):
	sID_2_libs = {}
	with open(ESS, 'r') as f:
		header_fields = f.readline().lower().strip().split('\t')
		header_fields = [header.rstrip('-') for header in header_fields]
		lib_id_idx = header_fields.index('library_id')
		batch_idx = header_fields.index('library_batch')
		exp_idx = header_fields.index('experiment')
		for idx, field in enumerate(header_fields):
			if len([s for s in assessment_header_keys if s in field]) > 0:
				assessment_idx = idx
				break
		for line in f:
			fields = line.split('\t')
			library_id = fields[lib_id_idx]
			experiment = fields[exp_idx].lower()
			sID = library_id.split('.')[0]
			batch_id = fields[batch_idx]
			if not(fields[assessment_idx]) or fields[assessment_idx].lower() == "control":
				assessment = ASSESSMENT_LEVELS["pass"]
			else:
				for key in ASSESSMENT_LEVELS:
					if key in fields[assessment_idx].lower():
						assessment = ASSESSMENT_LEVELS[key]
			if assessment < ASSESSMENT_LEVELS[threshold] and experiment != "raw" and sID[0] == "S":
				if sID in sID_2_libs:
					if len(sID_2_libs[sID]) >= doppelganger_n:
						raise ValueError(f'more than allowed number of doppelgangers for {sID}!')
					for lib in sID_2_libs[sID]:
						if lib[0] == library_id:
							raise ValueError(f'Duplicate library ID for {sID}!')
						if lib[1] == batch_id and doppelganger_n < 3:
							raise ValueError(f'Duplicate batch ID for {sID}!')
					sID_2_libs[sID].append((library_id, batch_id))
				else:
					sID_2_libs[sID] = [(library_id, batch_id)]
	return sID_2_libs

def get_singleton_libs(doppelganger_dict):
	return [doppelganger_dict[sID][0] for sID in doppelganger_dict if len(doppelganger_dict[sID]) == 1]

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	# parser.add_argument("--custom_fail_list", default=[], nargs='+', help="list custom values indicating library failure. Default is to fail when assessment contains 'fail' or 'dnu'")
	parser.add_argument('-t', "--threshold", default="fail", choices=['dnu', 'fail', 'questionable_critical', 'questionable'], help="threshold at which to consider a library failed")
	parser.add_argument("--assessment_header_keys", default=['call', 'assess'], nargs="+", help="keywords to use to dermine assessment column in ESS")
	parser.add_argument('-o', "--output_name", required=True, help="name for the output files, typically the name of the sequencing run")
	parser.add_argument('-m', "--max_doppelgangers", default=2, help="maximum number of doppelgangers per plate")
	parser.add_argument("ESS", help="ESS with assessment calls")

	args = parser.parse_args()

	doppelganger_dict = read_ESS(args.ESS, args.assessment_header_keys, args.threshold, int(args.max_doppelgangers))
	with open(f"{args.output_name}.singletons", 'w') as f:
		for lib in get_singleton_libs(doppelganger_dict):
			f.write('\t'.join(lib) + '\n')
	with open(f"{args.output_name}.doppelgangers", 'w') as dg:
		with open(f"{args.output_name}.bam_lookup", 'w') as bams:
			for sID in doppelganger_dict:
				libs = [lib[0] for lib in doppelganger_dict[sID]]
				dg.write('\t'.join([f'{sID}_{args.output_name}', sID.replace('S', "I")] + libs) + '\n')
				for lib in libs:
					nuc_bams = glob.glob(f'/n/groups/reich/matt/pipeline/released_libraries/{lib}/{lib}*wist*hg19*bam')
					mt_bams = glob.glob(f'/n/groups/reich/matt/pipeline/released_libraries/{lib}/{lib}*wist*rsrs*bam')
					lastest_nuc_bam = max(nuc_bams, key=os.path.getctime)
					lastest_mt_bam = max(mt_bams, key=os.path.getctime)
					bams.write('\t'.join([lib, str(lastest_nuc_bam), str(lastest_mt_bam)]) + '\n')
