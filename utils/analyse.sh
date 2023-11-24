#!/usr/bin/env python
'''
USAGE: ./analyse.sh <output file to analyse> <k - number of results (k of top-k)> <n - number of top results wanted (5, 10 etc)> <relative path to utils folder>
Call from the folder where the output file to be analysed is stored, the answer files would be stored in the same folder
'''

import sys
import subprocess as sp

if(len(sys.argv)!=5):
	print("USAGE: ./analyse.sh <output file to analyse> <k - number of results (k of top-k)> <n - number of top results wanted (5, 10 etc)> <relative path to utils folder>")
	exit(0)

ifile = sys.argv[1]
k = sys.argv[2]
n = sys.argv[3]
utils_path = sys.argv[4]

# find n best of top-k
# calling sort_output.sh
# taking top n results

tvar = int(k)+1
sort_output = sp.Popen(["./"+utils_path+"utils/sort_output.sh", ifile, str(tvar)], stdout=sp.PIPE, stderr=sp.STDOUT)
head_so = sp.Popen(["head", "-n", str(n)], stdin=sort_output.stdout, stdout=sp.PIPE, stderr=sp.STDOUT)
# stdout, stderr = head_so.communicate()
# print (stdout)

outfile = open(ifile+"_analysis.txt", 'w')

# dividing output in n lines

topn_sublist = list()
for l in head_so.stdout.readlines():
	topn_sublist.append(l)
#print(topn_sublist)

# splitting each line at the end of brackets and sorting based on the v_id in the bracket 
# cleaning brackets away

subg_dict = dict()
subg_dict_flag = 0	# signify that the dict has been initialised once (to check erroneous addition of more query vertex ids)
for item in topn_sublist:
#	print(item)
#	print(item.lstrip().replace(",  ==>", ""))
	line_split = item.lstrip().replace("  ==>\t", "").replace(", ", "\t").split('\t')
#	print(line_split)
# last element in line split holds the chi square value of the subgraph, 1st element holds its rank
	ele_count = len(line_split)
#	print(ele_count)
	rank = int(line_split[0])
	chisq = float(line_split[ele_count-1])
#	print(rank, chisq)
# re-initialise subg_dict without changing the keys
#	print(subg_dict)
	for key in subg_dict.keys():
		subg_dict[key]=None
#	print(subg_dict)
	for vids in line_split[1:ele_count-1]:
#		print(vids)
		ip_vid = vids.split()[0]
		q_vid = vids.split()[1].strip("()")
#		print(ip_vid, q_vid)
		if(subg_dict_flag==1):
			if(q_vid not in subg_dict):
				print("ERROR! New query vertex:" + q_vid)
				exit(0)
		subg_dict[q_vid] = ip_vid
#	print (subg_dict)
	output_header = "query_vid," + str(sorted(subg_dict)).strip("[]").replace(" ", " ,").replace("'", "") + "\n"
#	print(output_header)
# for each v_id of output subgraph run find_about_v.sh
	vid_line = "vid,"
	degree_line = "degree,"
	neigh_count_line = "#neighbours,"
	chisq_line = "chisq,"
	for key in sorted(subg_dict):
		if(subg_dict[key]!=None):
			vid_line = vid_line + subg_dict[key] + ", ,"
			vid_data = sp.Popen(["./"+utils_path+"utils/find_about_v.sh", ifile, subg_dict[key]], stdout=sp.PIPE, stderr=sp.STDOUT)
#			stdout, stderr = vid_data.communicate()
#			print(stdout)
			retrieve_data = list()
# interpretting the output
			for l in vid_data.stdout.readlines():
				if(":" in  l):
#					print(l.split(":")[1].strip())
					retrieve_data.append(l.split(":")[1].strip())
				else:
					retrieve_data.append(l.strip())
#			print(retrieve_data)
			degree_line = degree_line + retrieve_data[0] + ", ,"
			neigh_count_line = neigh_count_line = neigh_count_line + retrieve_data[2] + ", ,"
			chisq_line = chisq_line + retrieve_data[1].split()[0].strip() + ", ,"
		else:	# if value of subg_dict[key] is set to None
			vid_line = vid_line + " , ,"
			degree_line = degree_line + " , ,"
			neigh_count_line = neigh_count_line + " , ,"
			chisq_line = chisq_line + " , ,"
# appending chi-sq value of the subgraph and its rank in top-k (sorted by number of vertices in the answer followed by sort on chi-sq)
	vid_line = vid_line + str(chisq) + "," + line_split[0]
# original rank of the subgraph (in outupt sorted only on chi-sq)
	orig_rank_topk = sp.Popen(["tail", "-n", str(tvar), ifile], stdout=sp.PIPE, stderr=sp.STDOUT)
#	stdout, stderr = orig_rank_topk.communicate()
#	print(stdout)	
#	print(line_split[1].split()[0]+" (")
	orig_rank = sp.Popen(["grep", "-n", line_split[1].split()[0]+" ("], stdin=orig_rank_topk.stdout, stdout=sp.PIPE, stderr=sp.STDOUT)
	stdout, stderr = orig_rank.communicate()
#	print(stdout.split(":")[0])
#	print("stderr:\t"+str(stderr))
	vid_line = vid_line + "," + stdout.split(":")[0]
# storing the output in a file
	if(subg_dict_flag==0):
#		print(output_header)
		outfile.write(output_header)
#	print(vid_line + "\n" + degree_line + "\n" + neigh_count_line + "\n" + chisq_line)
	outfile.write(vid_line + "\n" + degree_line + "\n" + neigh_count_line + "\n" + chisq_line + "\n")
	outfile.flush()
	subg_dict_flag = 1
# store each set of subgraph vertices in a file to be searched for subgraph in actual graph
	rankfile = open(ifile+"_rank"+line_split[0], 'w')
	for  key in sorted(subg_dict):
		if(subg_dict[key]!=None):
			rankfile.write(subg_dict[key] + "\n")
	rankfile.close()

# find subgraph for each such file
	sp.call(["./"+utils_path+"plugins/grep_qg_edges_withp.sh", ifile+"_rank"+line_split[0], utils_path+"ppi_dataset/protein_edges_perturbed.txt"])
	sp.call(["mv", utils_path+"ppi_dataset/tmp_"+ifile+"_rank"+line_split[0], ifile+"_rank"+line_split[0]+"_edges"])
	print("rank" + line_split[0] + " done")
