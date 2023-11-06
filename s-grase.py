import igraph as ig
import pandas as pd
import sys
import os

def check_fromGTF_input(commandLineArg, argNum):
	"""
	Checks command line arguments of fromGTF.event.txt files (args 5 - 8). If the files open properly, then they are
	returned to main for further processing. These files will be read and parsed in order to map the coordinates of
	called events to the igraph object that is imported from the graphml file. rMATS and DEXSeq will both be mapped to
	this igraph object. Other command line arguments (2 - 4) are checked manually in main.

	:param commandLineArg: argv[] values passed in to check args 5 - 8
	:param argNum: The corresponding argument number of argv[] to give error codes.
	:return: fromGTF.event.txt file and the type of event (i.e. A3SS, A5SS, SE, RI)
	"""
	if "fromGTF" not in commandLineArg:
		print("\nError, command line argument " + argNum + " must be a fromGTF.event.txt file")
		exit(0)
	if "fromGTF.A3SS.txt" in commandLineArg:
		try:
			fromGTF_A3SS = open(commandLineArg)
			eventType = "fromGTF_A3SS"
			return fromGTF_A3SS, eventType
		except IOError:
			print("Oops! That fromGTF.A3SS.txt file does not exist. Try again...\n")
			exit(0)
	elif "fromGTF.A5SS.txt" in commandLineArg:
		try:
			fromGTF_A5SS = open(commandLineArg)
			eventType = "fromGTF_A5SS"
			return fromGTF_A5SS, eventType
		except IOError:
			print("Oops! That fromGTF.A5SS.txt file does not exist. Try again...\n")
			exit(0)
	elif "fromGTF.SE.txt" in commandLineArg:
		try:
			fromGTF_SE = open(commandLineArg)
			eventType = "fromGTF_SE"
			return fromGTF_SE, eventType
		except IOError:
			print("Oops! That fromGTF.SE.txt file does not exist. Try again...\n")
			exit(0)
	elif "fromGTF.RI.txt" in commandLineArg:
		try:
			fromGTF_RI = open(commandLineArg)
			eventType = "fromGTF_RI"
			return fromGTF_RI, eventType
		except IOError:
			print("Oops! That fromGTF.RI.txt file does not exist. Try again...\n")
			exit(0)
	else:
		print("Error, command line argument " + str(argNum) + " must be in the format 'fromGTF.event.txt'\n Event options: A3SS, A5SS, SE, RI")
		exit(0)


def map_DEXSeq_from_gff(g, gff):
	"""
	Takes a gff DEXSeq output file and reads it. The function will take the coordinates of DEXSeq exon fragments
	in order to create edges on the igraph object that map to those fragments. The fragments are labelled with
	a "dexseq_fragment" attribute with the value of the corresponding exonic part number from the gff file (i.e. E001).

	:param g: igraph object that has been imported from the graphml object read into this program
	:param gff: DEXSeq gff file that will be used to create fragment edges on the igraph object
	:return: igraph object after the DEXSeq edges have been added
	"""
	leftCoords = []
	rightCoords = []
	dex_frag = []

	g.es["dexseq_fragment"] = ''
	for x in gff:
		if x.split()[2] == "aggregate_gene":
			g["strand"] = x.split()[6]
			g["gene"] = x.split()[-1].strip('\"')
		if x.split()[2] == "exonic_part":
			leftCoords.append(x.split()[3])
			rightCoords.append(x.split()[4])
			dex_frag.append(x.split()[-1].strip('\"'))

	if g["strand"] == '-':
		for x in range(len(rightCoords)):
			rightCoords[x] = str(int(rightCoords[x]) + 1)
			g.add_edges([(rightCoords[x], leftCoords[x])])
			g.es[-1]["dexseq_fragment"] = dex_frag[x]

	if g["strand"] == '+':
		for x in range(len(rightCoords)):
			rightCoords[x] = str(int(rightCoords[x]) + 1)
			g.add_edges([(leftCoords[x], rightCoords[x])])
			g.es[-1]["dexseq_fragment"] = dex_frag[x]

	return g


def map_rMATS_event_overhang(g, fromGTF, eventType):
	"""
	Takes a fromGTF.event.txt rMATS output file and reads it. This function will take the coordinates of rMATS events in
	order to create edges on the igraph object that map those events with corresponding DEXSeq fragments. The goal is to
	map DEXSeq fragments that should be differentially expressed when alternative splicing (AS) events occur. This
	function specifically works with AS events that produce an overhang (A3SS and A5SS). Sometimes, an overhang may span
	over multiple DEXSeq fragments. This function accounts for that by iterating over the indices of nodes between the
	start and end coordinates of the rMATS event. Bifurcating paths of A3SS or A5SS events will be labelled "rmats
	short" or "rmats long", matching the fromGTF.event.txt file. The difference between the long exon and short exon
	will represent the overhang, and DEXSeq fragments spanning that overhang will be labelled with an edge attribute
	matching the rMATS event. In addition, this function will take the original fromGTF.event.txt input file and append
	a new column "DexseqFragment" that lists which DEXSeq exonic part number maps specifically to each rMATS event.
	This modified file will be output to the output directory specified in the command line arguments.

	:param g: igraph object that has been imported from the graphml object read into this program
	:param fromGTF: rMATS fromGTF.event.txt file that will be used to label DEXSeq edges on the igraph object with
					corresponding rMATS events. This will then be converted to a dataframe, and a column will be
					appended that will map rMATS event ID to DEXSeq fragment(s)
	:param eventType: Tracks rMATS event type (A3SS or A5SS) to label edges on the igrpah object appropriately
	:return: igraph object after the rMATS labels have been added to the DEXSeq edges appropriately.
	"""
	df = pd.read_table(fromGTF, dtype=str)
	fromGTF.seek(0)

	dx_ID = {}
	ID = []
	longES = []
	longEE = []
	shortES = []
	shortEE = []

	for x in fromGTF:
		if x.split()[0] == "ID":
			continue
		dx_ID[x.split()[0]] = []
		ID.append(x.split()[0])
		longES.append(x.split()[5])
		longEE.append(x.split()[6])
		shortES.append(x.split()[7])
		shortEE.append(x.split()[8])

	for x in range(len(longES)):
		longES[x] = str(int(longES[x]) + 1)
		longEE[x] = str(int(longEE[x]) + 1)
		shortES[x] = str(int(shortES[x]) + 1)
		shortEE[x] = str(int(shortEE[x]) + 1)
		if eventType == "A3SS":
			g.es.find(_within=(g.vs.find(longES[x]).index, g.vs.find(longEE[x]).index))["rmats"] = "rmats long"
			g.es.find(_within=(g.vs.find(shortES[x]).index, g.vs.find(shortEE[x]).index))["rmats"] = "rmats short"
			if longES[x] == shortES[x]:
				for i in range(g.vs.find(longEE[x]).index, g.vs.find(shortEE[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.find(_within=(i, i+1))[eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))

			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(longES[x]).index, g.vs.find(shortES[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.find(_within=(i, i+1))[eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
		if eventType == "A5SS":
			g.es.find(_within=(g.vs.find(longEE[x]).index, g.vs.find(longES[x]).index))["rmats"] = "rmats long"
			g.es.find(_within=(g.vs.find(shortEE[x]).index, g.vs.find(shortES[x]).index))["rmats"] = "rmats short"
			if longES[x] == shortES[x]:
				for i in range(g.vs.find(shortEE[x]).index, g.vs.find(longEE[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.find(_within=(i, i+1))["A5SS"] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(shortES[x]).index, g.vs.find(longES[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.find(_within=(i, i+1))["A5SS"] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	df['DexseqFragment'] = df['ID'].map(dx_ID)
	df.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)

	return g



def map_rMATS_event_full_fragment(g, fromGTF, eventType):
	"""
	Takes a fromGTF.event.txt rMATS output file and reads it. This function will take the coordinates of rMATS events in
	order to create edges on the igraph object that map those events with corresponding DEXSeq fragments. The goal is to
	map DEXSeq fragments that should be differentially expressed when alternative splicing (AS) events occur. This
	function specifically works with AS events that span a full exon (RI and SE). Sometimes, an exon may span over
	multiple DEXSeq fragments. This function accounts for that by iterating over the indices of nodes between the start
	and end coordinates of the rMATS event. DEXSeq fragments spanning the SE/RI exon will be labelled with an edge
	attribute matching the rMATS event. In addition, this function will take the original fromGTF.event.txt input file
	and append a new column "DexseqFragment" that lists which DEXSeq exonic part number maps specifically to each rMATS
	event. This modified file will be output to the output directory specified in the command line arguments.

	:param g: igraph object that has been imported from the graphml object read into this program
	:param fromGTF: rMATS fromGTF.event.txt file that will be used to label DEXSeq edges on the igraph object with
					corresponding rMATS events. This will then be converted to a dataframe, and a column will be
					appended that will map rMATS event ID to DEXSeq fragment(s)
	:param eventType: Tracks rMATS event type (SE or RI) to label edges on the igrpah object appropriately
	:return: igraph object after the rMATS labels have been added to the DEXSeq edges appropriately.
	"""
	df = pd.read_table(fromGTF, dtype=str)
	fromGTF.seek(0)

	dx_ID = {}
	ID = []
	exonStart = []
	exonEnd = []

	for x in fromGTF:
		if x.split()[0] == "ID":
			continue
		dx_ID[x.split()[0]] = []
		ID.append(x.split()[0])
		exonStart.append(x.split()[5])
		exonEnd.append(x.split()[6])

	for x in range(len(exonStart)):
		exonStart[x] = str(int(exonStart[x]) + 1)
		exonEnd[x] = str(int(exonEnd[x]) + 1)
		if g["strand"] == '+':
			for i in range(g.vs.find(exonStart[x]).index, g.vs.find(exonEnd[x]).index):
				for k in range(len(g.es.select(_within=(i, i+1)))):
					if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
						g.es.select(_within=(i, i+1))[k][eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
		if g["strand"] == '-':
			for i in range(g.vs.find(exonEnd[x]).index, g.vs.find(exonStart[x]).index):
				for k in range(len(g.es.select(_within=(i, i+1)))):
					if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
						g.es.select(_within=(i, i+1))[k][eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	df['DexseqFragment'] = df['ID'].map(dx_ID)
	df.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)

	return g









fromGTF_A3SS = fromGTF_A5SS = fromGTF_SE = fromGTF_RI = ""
args = len(sys.argv)

if args > 8 or args < 5:
	print("\nError, command line arguments should be:\n /path/to/output_directory   graphmlFile   gffFile  fromGTF.event.txtFile(s)")
	exit(0)

if not os.path.exists(sys.argv[1]):
	print("\nError, output_directory does not exist")
	exit(0)

if not sys.argv[2].endswith(".graphml"):
	print("\nError, first command line argument must be a *.graphml file")
	exit(0)
try:
	g =	ig.Graph.Read_GraphML(sys.argv[2])
except IOError:
	print("Oops! The graphml file does not exist. Try again...\n")
	exit(0)


if not sys.argv[3].endswith(".gff"):
	print("\nError, second command line argument must be a *.gff file")
	exit(0)
try:
	gff = open(sys.argv[3])
except IOError:
	print("Oops! The gff file does not exist. Try again...\n")
	exit(0)

fromGTF, eventType = check_fromGTF_input(sys.argv[4], 4)
locals()[eventType] = fromGTF

if args > 5:
	fromGTF, eventType = check_fromGTF_input(sys.argv[5], 5)
	locals()[eventType] = fromGTF

if args > 6:
	fromGTF, eventType = check_fromGTF_input(sys.argv[6], 6)
	locals()[eventType] = fromGTF

if args > 7:
	fromGTF, eventType = check_fromGTF_input(sys.argv[7], 7)
	locals()[eventType] = fromGTF



g = map_DEXSeq_from_gff(g, gff)
gff.close()



g.es["rmats"] = g.es["event"] = g.es["A3SS"] = g.es["A5SS"] = g.es["SE"] = g.es["RI"] = ""

if fromGTF_A3SS:
	g = map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS")
	fromGTF_A3SS.close()

if fromGTF_A5SS:
	g = map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS")
	fromGTF_A5SS.close()

if fromGTF_SE:
	g = map_rMATS_event_full_fragment(g, fromGTF_SE, "SE")
	fromGTF_SE.close()

if fromGTF_RI:
	g = map_rMATS_event_full_fragment(g, fromGTF_RI, "RI")
	fromGTF_RI.close()



edge_labels = []

for x in range(len(g.es)):
	if g.es[x]["A3SS"] == True:
		A3SS = "A3SS"
	else:
		A3SS = ""

	if g.es[x]["A5SS"] == True:
		A5SS = "A5SS"
	else:
		A5SS = ""

	if g.es[x]["SE"] == True:
		SE = "SE"
	else:
		SE = ""

	if g.es[x]["RI"] == True:
		RI = "RI"
	else:
		RI = ""

	edge_labels.append(g.es["dexseq_fragment"][x] + " " + A3SS + " " + A5SS + " " + SE + " " + RI)

color_dict = {"ex": "orange", "in": "grey", "NA": "grey", None: "dark green"}
curved_dict = {"ex": -0.3, "in": False, "NA": False, None: 0}
width_dict = {"ex": 5, "in": 2, "NA": 2, None: 5}
order_dict = {}
for name in g.vs['name']:
	order_dict[name] = name

order_dict["L"] = "100000000000"
order_dict["R"] = "0"

for name in order_dict:
	order_dict[name] = int(order_dict[name])

visual_style = {}

visual_style["edge_curved"] = [curved_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]]
visual_style["edge_color"] = [color_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]]
visual_style["edge_width"] = [width_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]]
visual_style["order"] = [order_dict[order] for order in g.vs["name"]]
visual_style["vertex_label"] = g.vs["id"]
visual_style["edge_arrow_size"] = 0.001
visual_style["edge_label"] = edge_labels
visual_style["vertex_shape"] = "hidden"
visual_style["vertex_label_size"] = 25
visual_style["edge_label_size"] = 25
visual_style["bbox"] = (3500, 1000)
visual_style["margin"] = 100


layout = g.layout_sugiyama()
layout.rotate(270)

ig.plot(g, sys.argv[1] + "/graph_" + g["gene"] + ".png", layout=layout, **visual_style)

#g.write_graphml("updated_" + g["gene"] + ".graphml")




######################################################################## code that has been converted to function calls

########################## gff mapping
#lines = gff.readlines()
#
#leftCoords = []
#rightCoords = []
#dex_frag = []
#
#g.es["dexseq_fragment"] = ''
#
#for x in lines:
#	if x.split()[2] == "aggregate_gene":
#		g["strand"] = x.split()[6]
#		g["gene"] = x.split()[-1].strip('\"')
#	if x.split()[2] == "exonic_part":
#		leftCoords.append(x.split()[3])
#		rightCoords.append(x.split()[4])
#		dex_frag.append(x.split()[-1].strip('\"'))
#
#if g["strand"] == '-':
#	for x in range(len(rightCoords)):
#		rightCoords[x] = str(int(rightCoords[x]) + 1)
#		g.add_edges([(rightCoords[x], leftCoords[x])])
#		g.es[-1]["dexseq_fragment"] = dex_frag[x]
#
#if g["strand"] == '+':
#	for x in range(len(rightCoords)):
#		rightCoords[x] = str(int(rightCoords[x]) + 1)
#		g.add_edges([(leftCoords[x], rightCoords[x])])
#		g.es[-1]["dexseq_fragment"] = dex_frag[x]


########################## A3SS mapping
#lines = fromGTF_A3SS.readlines()
#
#dx_ID_A3SS = {}
#ID_A3SS = []
#longES_A3SS = []
#longEE_A3SS = []
#shortES_A3SS = []
#shortEE_A3SS = []
#
#for x in lines:
#	if x.split()[0] == "ID":
#		continue
#	dx_ID_A3SS[x.split()[0]] = []
#	ID_A3SS.append(x.split()[0])
#	longES_A3SS.append(x.split()[5])
#	longEE_A3SS.append(x.split()[6])
#	shortES_A3SS.append(x.split()[7])
#	shortEE_A3SS.append(x.split()[8])
#
#for x in range(len(longES_A3SS)):
#	longES_A3SS[x] = str(int(longES_A3SS[x]) + 1)
#	longEE_A3SS[x] = str(int(longEE_A3SS[x]) + 1)
#	shortES_A3SS[x] = str(int(shortES_A3SS[x]) + 1)
#	shortEE_A3SS[x] = str(int(shortEE_A3SS[x]) + 1)
#	g.es.find(_within=(g.vs.find(longES_A3SS[x]).index, g.vs.find(longEE_A3SS[x]).index))["rmats"] = "rmats long"
#	g.es.find(_within=(g.vs.find(shortES_A3SS[x]).index, g.vs.find(shortEE_A3SS[x]).index))["rmats"] = "rmats short"
#	if longES_A3SS[x] == shortES_A3SS[x]:
#		for i in range(g.vs.find(longEE_A3SS[x]).index, g.vs.find(shortEE_A3SS[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.find(_within=(i, i+1))["A3SS"] = True
#					dx_ID_A3SS[ID_A3SS[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#
#	if longEE_A3SS[x] == shortEE_A3SS[x]:
#		for i in range(g.vs.find(longES_A3SS[x]).index, g.vs.find(shortES_A3SS[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.find(_within=(i, i+1))["A3SS"] = True
#					dx_ID_A3SS[ID_A3SS[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#
#for x in dx_ID_A3SS:
#	dx_ID_A3SS[x] = ','.join(dx_ID_A3SS[x])
#
#df_A3SS['DexseqFragment'] = df_A3SS['ID'].map(dx_ID_A3SS)
#df_A3SS.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + ".A3SS.txt", sep = '\t', index = False)


########################## A5SS mapping
#lines = fromGTF_A5SS.readlines()
#
#dx_ID_A5SS = {}
#ID_A5SS = []
#longES_A5SS = []
#longEE_A5SS = []
#shortES_A5SS = []
#shortEE_A5SS = []
#
#for x in lines:
#	if x.split()[0] == "ID":
#		continue
#	dx_ID_A5SS[x.split()[0]] = []
#	ID_A5SS.append(x.split()[0])
#	longES_A5SS.append(x.split()[5])
#	longEE_A5SS.append(x.split()[6])
#	shortES_A5SS.append(x.split()[7])
#	shortEE_A5SS.append(x.split()[8])
#
#for x in range(len(longES_A5SS)):
#	longES_A5SS[x] = str(int(longES_A5SS[x]) + 1)
#	longEE_A5SS[x] = str(int(longEE_A5SS[x]) + 1)
#	shortES_A5SS[x] = str(int(shortES_A5SS[x]) + 1)
#	shortEE_A5SS[x] = str(int(shortEE_A5SS[x]) + 1)
#	g.es.find(_within=(g.vs.find(longEE_A5SS[x]).index, g.vs.find(longES_A5SS[x]).index))["rmats"] = "rmats long"
#	g.es.find(_within=(g.vs.find(shortEE_A5SS[x]).index, g.vs.find(shortES_A5SS[x]).index))["rmats"] = "rmats short"
#	if longES_A5SS[x] == shortES_A5SS[x]:
#		for i in range(g.vs.find(shortEE_A5SS[x]).index, g.vs.find(longEE_A5SS[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.find(_within=(i, i+1))["A5SS"] = True
#					dx_ID_A5SS[ID_A5SS[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#	if longEE_A5SS[x] == shortEE_A5SS[x]:
#		for i in range(g.vs.find(shortES_A5SS[x]).index, g.vs.find(longES_A5SS[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.find(_within=(i, i+1))["A5SS"] = True
#					dx_ID_A5SS[ID_A5SS[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#
#for x in dx_ID_A5SS:
#	dx_ID_A5SS[x] = ','.join(dx_ID_A5SS[x])
#
#df_A5SS['DexseqFragment'] = df_A5SS['ID'].map(dx_ID_A5SS)
#df_A5SS.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + ".A5SS.txt", sep = '\t', index = False)


########################## SE mapping
#lines = fromGTF_SE.readlines()
#
#dx_ID_SE = {}
#ID_SE = []
#exonStart_SE = []
#exonEnd_SE = []
#
#for x in lines:
#	if x.split()[0] == "ID":
#		continue
#	dx_ID_SE[x.split()[0]] = []
#	ID_SE.append(x.split()[0])
#	exonStart_SE.append(x.split()[5])
#	exonEnd_SE.append(x.split()[6])
#
#for x in range(len(exonStart_SE)):
#	exonStart_SE[x] = str(int(exonStart_SE[x]) + 1)
#	exonEnd_SE[x] = str(int(exonEnd_SE[x]) + 1)
#	if g["strand"] == '+':
#		for i in range(g.vs.find(exonStart_SE[x]).index, g.vs.find(exonEnd_SE[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.select(_within=(i, i+1))[k]["SE"] = True
#					dx_ID_SE[ID_SE[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#	if g["strand"] == '-':
#		for i in range(g.vs.find(exonEnd_SE[x]).index, g.vs.find(exonStart_SE[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.select(_within=(i, i+1))[k]["SE"] = True
#					dx_ID_SE[ID_SE[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#
#for x in dx_ID_SE:
#	dx_ID_SE[x] = ','.join(dx_ID_SE[x])
#
#df_SE['DexseqFragment'] = df_SE['ID'].map(dx_ID_SE)
#df_SE.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + ".SE.txt", sep = '\t', index = False)


########################## RI mapping
#lines = fromGTF_RI.readlines()
#
#dx_ID_RI = {}
#ID_RI = []
#exonStart_RI = []
#exonEnd_RI = []
#
#for x in lines:
#	if x.split()[0] == "ID":
#		continue
#	dx_ID_RI[x.split()[0]] = []
#	ID_RI.append(x.split()[0])
#	exonStart_RI.append(x.split()[5])
#	exonEnd_RI.append(x.split()[6])
#
#for x in range(len(exonStart_RI)):
#	exonStart_RI[x] = str(int(exonStart_RI[x]) + 1)
#	exonEnd_RI[x] = str(int(exonEnd_RI[x]) + 1)
#	if g["strand"] == '+':
#		for i in range(g.vs.find(exonStart_RI[x]).index, g.vs.find(exonEnd_RI[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.select(_within=(i, i+1))[k]["RI"] = True
#					dx_ID_RI[ID_RI[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#	if g["strand"] == '-':
#		for i in range(g.vs.find(exonEnd_RI[x]).index, g.vs.find(exonStart_RI[x]).index):
#			for k in range(len(g.es.select(_within=(i, i+1)))):
#				if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
#					g.es.select(_within=(i, i+1))[k]["RI"] = True
#					dx_ID_RI[ID_RI[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
#
#for x in dx_ID_RI:
#	dx_ID_RI[x] = ','.join(dx_ID_RI[x])
#
#df_RI['DexseqFragment'] = df_RI['ID'].map(dx_ID_RI)
#df_RI.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + ".RI.txt", sep = '\t', index = False)


########################## argument checking for fromGTF files
#if "fromGTF" not in sys.argv[4]:
#	print("\nError, third command line argument must be a fromGTF.event.txt file")
#	exit(0)
#if "fromGTF.A3SS.txt" in sys.argv[4]:
#	try:
#		fromGTF_A3SS = open(sys.argv[4])
#		#df_A3SS = pd.read_table(sys.argv[4], dtype = str)
#	except IOError:
#		print("Oops! That fromGTF.A3SS.txt file does not exist. Try again...\n")
#		exit(0)
#elif "fromGTF.A5SS.txt" in sys.argv[4]:
#	try:
#		fromGTF_A5SS = open(sys.argv[4])
#		df_A5SS = pd.read_table(sys.argv[4], dtype = str)
#	except IOError:
#		print("Oops! That fromGTF.A5SS.txt file does not exist. Try again...\n")
#		exit(0)
#elif "fromGTF.SE.txt" in sys.argv[4]:
#	try:
#		fromGTF_SE = open(sys.argv[4])
#		df_SE = pd.read_table(sys.argv[4], dtype = str)
#	except IOError:
#		print("Oops! That fromGTF.SE.txt file does not exist. Try again...\n")
#		exit(0)
#elif "fromGTF.RI.txt" in sys.argv[4]:
#	try:
#		fromGTF_RI = open(sys.argv[4])
#		df_RI = pd.read_table(sys.argv[4], dtype = str)
#	except IOError:
#		print("Oops! That fromGTF.RI.txt file does not exist. Try again...\n")
#		exit(0)
#else:
#	print("Error, file must be in the format 'fromGTF.event.txt'\n Event options: A3SS, A5SS, SE, RI")
#	exit(0)
#
#
#if args > 5:
#	if "fromGTF" not in sys.argv[5]:
#		print("\nError, third command line argument must be a fromGTF.event.txt file")
#		exit(0)
#	if "fromGTF.A3SS.txt" in sys.argv[5]:
#		try:
#			fromGTF_A3SS = open(sys.argv[5])
#			df_A3SS = pd.read_table(sys.argv[5], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.A3SS.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.A5SS.txt" in sys.argv[5]:
#		try:
#			fromGTF_A5SS = open(sys.argv[5])
#			df_A5SS = pd.read_table(sys.argv[5], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.A5SS.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.SE.txt" in sys.argv[5]:
#		try:
#			fromGTF_SE = open(sys.argv[5])
#			df_SE = pd.read_table(sys.argv[5], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.SE.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.RI.txt" in sys.argv[5]:
#		try:
#			fromGTF_RI = open(sys.argv[5])
#			df_RI = pd.read_table(sys.argv[5], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.RI.txt file does not exist. Try again...\n")
#			exit(0)
#	else:
#		print("Error, file must be in the format 'fromGTF.event.txt'\n Event options: A3SS, A5SS, SE, RI")
#		exit(0)
#
#
#if args > 6:
#	if "fromGTF" not in sys.argv[6]:
#		print("\nError, third command line argument must be a fromGTF.event.txt file")
#		exit(0)
#	if "fromGTF.A3SS.txt" in sys.argv[6]:
#		try:
#			fromGTF_A3SS = open(sys.argv[6])
#			df_A3SS = pd.read_table(sys.argv[6], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.A3SS.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.A5SS.txt" in sys.argv[6]:
#		try:
#			fromGTF_A5SS = open(sys.argv[6])
#			df_A5SS = pd.read_table(sys.argv[6], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.A5SS.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.SE.txt" in sys.argv[6]:
#		try:
#			fromGTF_SE = open(sys.argv[6])
#			df_SE = pd.read_table(sys.argv[6], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.SE.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.RI.txt" in sys.argv[6]:
#		try:
#			fromGTF_RI = open(sys.argv[6])
#			df_RI = pd.read_table(sys.argv[6], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.RI.txt file does not exist. Try again...\n")
#			exit(0)
#	else:
#		print("Error, file must be in the format 'fromGTF.event.txt'\n Event options: A3SS, A5SS, SE, RI")
#		exit(0)
#
#
#if args > 7:
#	if "fromGTF" not in sys.argv[7]:
#		print("\nError, third command line argument must be a fromGTF.event.txt file")
#		exit(0)
#	if "fromGTF.A3SS.txt" in sys.argv[7]:
#		try:
#			fromGTF_A3SS = open(sys.argv[7])
#			df_A3SS = pd.read_table(sys.argv[7], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.A3SS.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.A5SS.txt" in sys.argv[7]:
#		try:
#			fromGTF_A5SS = open(sys.argv[7])
#			df_A5SS = pd.read_table(sys.argv[7], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.A5SS.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.SE.txt" in sys.argv[7]:
#		try:
#			fromGTF_SE = open(sys.argv[7])
#			df_SE = pd.read_table(sys.argv[7], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.SE.txt file does not exist. Try again...\n")
#			exit(0)
#	elif "fromGTF.RI.txt" in sys.argv[7]:
#		try:
#			fromGTF_RI = open(sys.argv[7])
#			df_RI = pd.read_table(sys.argv[7], dtype = str)
#		except IOError:
#			print("Oops! That fromGTF.RI.txt file does not exist. Try again...\n")
#			exit(0)
#	else:
#		print("Error, file must be in the format 'fromGTF.event.txt'\n Event options: A3SS, A5SS, SE, RI")
#		exit(0)
