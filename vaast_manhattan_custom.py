#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, random

# Copyright (c) Brett J. Kennedy, Gordon Lemmon, W. Scott Watkins 2017 MIT License

__authors__="Brett J. Kennedy and Gordon Lemmon and Scott Watkins"
__date__="January 24 2015"
__lastModified__="Scott Watkins for PCGC figure printing, Dec 1, 2016"

#Example used to create a small phevor plot for the PCGC data:
#vaast_manhattan_custom.py --phevor2 02504.phevor2 02504.phevor2.png  ~/bin/Phevor2/ontos/genes.gff3 --title 02504 --mm_width 50 --mm_height 50 --red_genes MYH6


plt.rcParams['figure.autolayout'] = True

def rotate(nist):
    """Spins a list"""
    return nist[1:] + nist[:1]

def parse_info(info):
    """creates a dictionary of the info field in a GFF3"""
    infoD = {}
    info = info.split(';')
    for e in info:
            e = e.split('=')
            try:    infoD[e[0]] = e[1]
            except: continue
    return infoD

def parse_scores(input,phevor2):
    """Parses either .simple or .phevor output and returns dictonary\
    of sores by gene."""
    scores = {}
    highs = 0
    if phevor2 == False:
        with open(input) as t:
            for line in csv.reader(t,delimiter = "\t"):
                if line[0] == "RANK": continue
                sco = -np.log10(float(line[2]))
                scores[line[1].strip()] = sco
                if sco>highs:
                    highs = sco
    
    if phevor2 == True:
        with open(input) as t:
            for line in csv.reader(t,delimiter = "\t"):
                if "#" in line[0]:  continue
                sco = float(line[2])
                scores[line[1].strip()] = sco
                if sco>highs:
                    highs = sco
    return scores,highs

def populate_scores_coords(gff3,scores):
    """Parses the gff3 and creates a dictionary of genomic\
     coordinates 'coord' and 'scores'"""
    genes = {}
    with open(gff3) as q:
        for line in csv.reader(q,delimiter = "\t"):
            if "#" in line[0]:  continue
            if line[2] !=  "gene": continue
            info = parse_info(line[8])
            name = info['Name']
            chrom = line[0]
            start = int(line[3])
            stop = int(line[4])
            if line[0] not in genes.keys():
                genes[line[0]] = {}
            genes[line[0]][name] = {}
            if name in scores.keys():
                genes[line[0]][name]['score'] = scores[name]
            else:
                genes[line[0]][name]['score'] = 0
            coord = random.randint(start,stop)
            genes[line[0]][name]['coord'] = coord
    return genes

def last_list(genes,chroms):
    last = [0]
    for i in range(len(chroms)):
    	alcor = []
    	for j in genes[chroms[i]].keys():
    		alcor.append(genes[chroms[i]][j]['coord'])
    	if len(last) > 1:
    		last.append(max(alcor)+last[i])
    	else:
    		last.append(max(alcor))
    return last

def parse_args():
		parser = argparse.ArgumentParser(description = "Creates a manhattan plot from VAAST simple  or Phevor output", epilog = "Author: Brett J. Kennedy / January 24 2015")
		parser.add_argument("input", help = "simple vaast or pVAAST output")
		parser.add_argument("output", help = "save file to this handle",type = str)
		parser.add_argument("gff3", help = "GFF file of genomic coordinates of genes")
		parser.add_argument("--red_genes",help = "comma seperated genes of interest to label in red",default = None)
		parser.add_argument("--blue_genes",help = "comma seperated genes of interest to label in blue",default = None)
		parser.add_argument("--phevor2",help = "use this arguemnt if the input is phevor2 rather than VAAST/pVAAST",action = 'store_true',default = False)
		parser.add_argument("--xlab",help = "x-axis label",type = str,default = None)
		parser.add_argument("--ylab",help = "y-axis label",type = str,default = None)
		parser.add_argument("--title",help = "title to use above the plot",type = str,default = None)
		width = parser.add_mutually_exclusive_group()
		width.add_argument("--inch_width", help = "width in inches",type = float,default = 2)
		width.add_argument("--mm_width", help="width in mm",type = float,default = None)
		height = parser.add_mutually_exclusive_group()
		height.add_argument("--inch_height", help = "height in inches",type = float,default = 2)
		height.add_argument("--mm_height", help = "height in mm",type = float,default = None)
		parser.add_argument("--point_size", help = "size of points",type = float, default = 100)
		parser.add_argument("--red_label", help = "red gene label1",type = str, default = None)
		parser.add_argument("--red_label2", help = "red gene label2",type = str, default = None)	
		parser.add_argument("--blue_label", help = "blue gene label",type = str, default = None)
		parser.add_argument("--rpos", help = "align red gene label, left|right",type = str, default = None)
		parser.add_argument("--bpos", help = "align blue gene label, left|right",type = str, default = None)
		
		return parser.parse_args()

# genes is a gene name or list of names
def plot_stuff(
		input,
		output,
		gff3,
		phevor2 = False,
		blue_genes = [],
		red_genes = [],
		xlab = None,
		width = 15,
		height = 10,
		point_size = 100
):
	if not isinstance(blue_genes, list):
		blue_genes = [blue_genes]
	if not isinstance(red_genes, list):
		red_genes = [red_genes]
	chroms = list(np.arange(1,23))+['X','Y']
	chroms = ["chr"+str(i) for i in chroms]
		
	scores,highS = parse_scores(input,phevor2)
	genes = populate_scores_coords(gff3,scores)
	last = last_list(genes,chroms)

	corr = []
	scor = []
	t = 0
	blueD = {}
	redD = {}
	for k in chroms:
		for j in genes[k].keys():
			if j in blue_genes:
				blueD[j] = [genes[k][j]['coord']+last[t], genes[k][j]['score']]
			if j in red_genes:
				redD[j] = [genes[k][j]['coord']+last[t], genes[k][j]['score']]
			corr.append(genes[k][j]['coord']+last[t])
			scor.append(genes[k][j]['score'])
		t+= 1

	centering = [(last[i]+last[i+1])/2 for i in range(len(last)) if i < len(last)-1]
	colors = ['b','g','r','y']

	#plt.figure(figsize = (width,height))
	#ax = plt.subplot()

        fig, ax = plt.subplots(1, figsize=(width, height))

#set INSIDE tick parameters
	plt.tick_params(top='off', bottom='off')



# Set y-high limit to highS+1 for auto scaling
	ax.set_ylim([-0.05,10])
	ax.tick_params(axis='y', labelsize=6)

	ax.set_xlim([min(corr),max(corr)])

#	ax.set_xticks(centering)

# white out the chr on x axis with nada
# Set nada back to chroms for display of all chromes	
	nada = ''
	ax.set_xticklabels(nada,rotation = 45)



	if xlab:
            xlab  =  xlab
	else:
            xlab = 'Chromosome'
	if phevor2 == True:
		ylab = args.ylab if args.ylab else "Phevor Score"
		title = args.title if args.title else "Genes Re-Ranked by Phevor"
	else:
		ylab = args.ylab if args.ylab else "$-log_{10}$ pVAAST p-value"
		title = args.title if args.title else "Genes Scored by pVAAST"
        if args.title == "NONE":
                title = ""
        if args.xlab == "NONE":
                xlab = ""
        if args.ylab == "NONE":
                ylab = ""
        plt.xlabel(xlab,fontsize = 8, fontweight='bold')
	plt.ylabel(ylab,fontsize = 8, fontweight='bold')
	plt.title(title,fontsize = 12, fontweight='bold')
	t = 0
	for k in chroms:
		zero_coor = []
		zero_scor = []
		corr = []
		scor = []
		for j in genes[k].keys():
			if genes[k][j]['score'] > 0:
				corr.append(genes[k][j]['coord']+last[t])
				scor.append(genes[k][j]['score'])
			else:
				zero_coor.append(genes[k][j]['coord']+last[t])
				zero_scor.append(genes[k][j]['score'])
		plt.scatter(zero_coor,zero_scor,c = colors[0],marker = '.',s = 200,edgecolors = 'none')
		plt.scatter(corr,scor,c = colors[0],marker = '.',s = point_size,edgecolors = 'black',)
		colors = rotate(colors)
		t+=1
	for g in blue_genes:
		if args.blue_label != None:
			plt.text(blueD[g][0],blueD[g][1]+1.2,g,fontsize = 8,color = 'blue',weight = 'bold', horizontalalignment = bpos )
                	plt.text(blueD[g][0],blueD[g][1]+0.5,blue_label,fontsize = 6,color = 'black',weight = 'normal',horizontalalignment = bpos, style = 'italic')

		else:
			plt.text(blueD[g][0]+0.5,blueD[g][1]+0.4,g,fontsize = 6,color = 'blue',horizontalalignment = bpos)

	for g in red_genes:
		if args.red_label != None and args.red_label2 != None:
			plt.text(redD[g][0],redD[g][1]+1.5,g,fontsize = 8,color = 'red',weight = 'bold', horizontalalignment = rpos )
                	plt.text(redD[g][0],redD[g][1]+1.0,red_label,fontsize = 5,color = 'black',weight = 'normal',horizontalalignment = rpos, style = 'italic')
                	plt.text(redD[g][0],redD[g][1]+0.5,red_label2,fontsize = 5,color = 'black',weight = 'normal',horizontalalignment = rpos, style = 'italic')

		elif args.red_label != None:
			plt.text(redD[g][0],redD[g][1]+1.2,g,fontsize = 8,color = 'red',weight = 'bold', horizontalalignment = rpos )
                	plt.text(redD[g][0],redD[g][1]+0.5,red_label,fontsize = 6,color = 'black',weight = 'normal',horizontalalignment = rpos, style = 'italic')

		else:
			plt.text(redD[g][0],redD[g][1]+0.5,g,fontsize = 8,color = 'red',weight = 'bold', horizontalalignment = rpos )



	#plot a reference line
	plt.axhline(y=2.3, linewidth=0.5, linestyle='--', color = 'r')		

	plt.subplots_adjust()


	plt.savefig(output,dpi = 1200)


	#plt.show()

if __name__ == "__main__":
	args = parse_args()

	blue_genes = []
	if args.blue_genes != None:
		blue_genes = args.blue_genes.split(',')

	red_genes = []
	if args.red_genes != None:
		red_genes = args.red_genes.split(',')

	width = args.inch_width
	height = args.inch_height

	red_label = args.red_label
	red_label2 = args.red_label2

	rpos = 'center'
	if args.rpos != None:
		rpos = args.rpos

	blue_label = args.blue_label

	bpos = 'center'
	if args.bpos != None:
		bpos = args.bpos

	if args.mm_width:
		width = args.mm_width / 25.4
	if args.mm_height:
		height = args.mm_height / 25.4
	plot_stuff(
			args.input,
			args.output,
			args.gff3,
			args.phevor2,
			blue_genes,
			red_genes,
			args.xlab,
			width,
			height,
			args.point_size
			
	)
