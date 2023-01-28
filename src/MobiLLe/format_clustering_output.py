
import sys
#===================================================================================
#Global variables
Input_file = ""
Output_file= ""


#Usage
usage = "python format_clustering_output.py -i Input_file -o Output_file \n"

#===================================================================================
#Read parameters
def readParameters(args):
	global Input_file
	global Output_file

	for i in range(1,len(args)):
		if (args[i] == "-i"):
			Input_file = args[i+1]
		elif (args[i] == "-o"):
			Output_file = args[i+1]
		elif (args[i] == "-h"):
			print (usage)
#===================================================================================
### Check parameters
def checkParameters():
	if (Input_file == ""):
		print ("ERROR::Parameter -i Input_file is required\n")
		sys.exit(1);
	elif (Output_file == ""):
		print ("ERROR::Parameter -o Output_file is required\n")
		sys.exit(1);

#===================================================================================
def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

#===================================================================================
def dico_format(lines):
	Dico_first={}
	for l in range(0,len(lines)):
		split=lines[l].split("\t")
		if split[0] in Dico_first.keys():
			Dico_first[split[0] ].append(split[1])	
		else:
			Dico_first[split[0] ]=[split[1]]
	return Dico_first
#===================================================================================

#===================================================================================				
def write_file(dico):
	f=open(Output_file,"w")
	for key in dico.keys():
		name=str(key)+"\t"
		f.write(name)
		for seq in range(len(dico[key])-1):
			sequence=str(dico[key][seq])+" "
			f.write(str(sequence))
		sequence=str(dico[key][-1])
		f.write(str(sequence))
		f.write("\n")
	f.close()
	return 0

#===================================================================================
#			    		Main
#===================================================================================
readParameters(sys.argv)
checkParameters()
lines=read_output_file(Input_file)
Dico_first= dico_format(lines)
write_file(Dico_first)
