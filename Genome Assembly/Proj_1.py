#/usr/bin/python
import sys
# encode_hash takes a string and returns an integer representing 
# a sequence
def encode_hash(sequence):
	n_index = 0;
	sequence_num = 0;
	for letter in sequence:
		if(letter == 'A'):
			n = 0
		elif(letter == 'C'):
			n = 1
		elif(letter == 'G'):
			n = 2
		elif(letter == 'T'):
			n = 3
		sequence_num = sequence_num + n*(4**n_index)
		n_index+=1
	return sequence_num

def load_reads(filename):
	print "Loading Reads.."
	hash_table = {}
	read_table = {}
	alpha_read_table = {}
	degree_table = {}
	read_count = 0
	edge_string = ""
	with open(filename) as reads_file:
		for line in reads_file:
			if(line[0] != '>' and line[0] != '\n'):
				# Append Read to the string "Data"
				cur_read = line.replace('\n','')
				# split read into a prefix and a suffix
				prefix, suffix = cur_read[:len(cur_read)/2], cur_read[len(cur_read)/2:]
				# convert prefix to hash code
				prefix_hash = encode_hash(prefix)
				suffix_hash = encode_hash(suffix)
				# add encoded prefix and suffix to read_table hash table
				read_table[read_count] = (prefix_hash,suffix_hash)
				alpha_read_table[read_count] = (prefix,suffix) 
				# convert read number to string
				edge_string = str(read_count);
				if (hash_table.has_key(prefix_hash) and hash_table[prefix_hash].has_key(suffix_hash)):
					hash_table[prefix_hash][suffix_hash] = hash_table[prefix_hash][suffix_hash] + "," + edge_string
					##print "Same prefix/suffix"
					##print hash_table[prefix_hash]
				elif (hash_table.has_key(prefix_hash)):
					hash_table[prefix_hash].update({suffix_hash: edge_string})
					##print "Same prefix"
					##print hash_table[prefix_hash].get(suffix_hash, default = None)
				else:
					hash_table[prefix_hash] = {suffix_hash: edge_string}
					##print "unique"
					##print hash_table[prefix_hash].get(suffix_hash, default = None)
				#Populate degree_table storing indegrees and outdegrees
				if(degree_table.has_key(prefix_hash)==0):
					degree_table[prefix_hash] = [1,0]
				elif(degree_table.has_key(prefix_hash) == 1):
					outdegree = degree_table[prefix_hash][0]
					outdegree = outdegree+1
					indegree = degree_table[prefix_hash][1]
					degree_table[prefix_hash]=[outdegree, indegree]	
				if(degree_table.has_key(suffix_hash)==1):
					outdegree = degree_table[suffix_hash][0]
					indegree = degree_table[suffix_hash][1]
					indegree = indegree + 1
					degree_table[suffix_hash]=[outdegree, indegree]
				elif(degree_table.has_key(suffix_hash)==0):
					degree_table[suffix_hash] = [0,1]
				read_count+=1
	print "Done."
	return [hash_table, read_table, degree_table, alpha_read_table]

class Path_Assembler(object):

	def __init__(self, filename):
		self.GraphHash, self.ReadHash, self.DegreeHash, self.AlphaReadHash = load_reads(filename)
		self.visited_nodes = []
		self.paths = []
		self.cycles = set()
		self.currentpath = []
		self.contigs = []

	def print_reads(self):
		for key in self.ReadHash.keys():
			#print "line number: " + str(key)
			print self.ReadHash[key]
	
	# finds the first instance where outdegree > indegree
	def start_node(self):
		#degree_hash = self.DegreeHash
		for key in self.DegreeHash.keys():
			if self.DegreeHash[key][0] > self.DegreeHash[key][1]:
				return key
		for key in self.DegreeHash.keys():
			if self.DegreeHash[key][0] == self.DegreeHash[key][1] and self.DegreeHash[key][0] != 0:
				return key
		return -1
		#found_nodes = (key for key, val in degree_hash.items() if val[0] > val[1] in val)

	def print_paths(self):
		print self.paths

	def start_assembly(self):
		print "Assembling Paths..."
		path_bool = 1
		key = self.start_node()
		while(key != -1 and path_bool != -1):
			path_bool = self.find_paths(key)
			key = self.start_node()
		print "Done."

	#chooses a value not equal to '', avoids choice of case where prefix == suffix if same is 1 does not if 0
	def suffix_picker(self, prefix, same):
		if(same == 0):
			for suffix in self.GraphHash[prefix].keys():
				if(self.GraphHash[prefix][suffix] != '' and prefix != suffix):
					return self.GraphHash[prefix][suffix]
				elif(self.GraphHash[prefix][suffix] != '' and prefix == suffix):
					# Add case where prefix == suffix
					rindex = self.GraphHash[prefix][suffix]
					indexlist = rindex.split(',')
					new_edge_list = indexlist[0:len(indexlist)-1]
					rindex = indexlist[len(indexlist)-1]
					new_suffix_value = ','.join(new_edge_list)
					new_suffix_value = str(new_suffix_value)
					rindexint = int(rindex)
					# store suffix value
					start_suffix = self.ReadHash[rindexint][1]
					self.GraphHash[prefix][suffix] = new_suffix_value
					self.visited_nodes.append(self.ReadHash[rindexint][0])
					self.currentpath.append(self.AlphaReadHash[rindexint][0])
					self.visited_nodes.append(self.ReadHash[rindexint][1])
					self.currentpath.append(self.AlphaReadHash[rindexint][1])
					Prefix_out = self.DegreeHash[prefix][0]-1
					Prefix_in = self.DegreeHash[prefix][1]
					self.DegreeHash[suffix] = [Prefix_out, Prefix_in]
					Suffix_out = self.DegreeHash[suffix][0]
					Suffix_in = self.DegreeHash[suffix][1]-1
					self.DegreeHash[suffix] = [Suffix_out, Suffix_in]
		elif(same ==1):
			for suffix in self.GraphHash[prefix].keys():
				if(self.GraphHash[prefix][suffix] != ''):
					return self.GraphHash[prefix][suffix]



	def find_paths(self, start):
		if(len(self.GraphHash[start].values()) > 0 and self.GraphHash[start].values()[len(self.GraphHash[start].values())-1] != '' and self.DegreeHash[start][0] > 0):				
			rindex = self.suffix_picker(start,1)
			# if the edgebag contains multiple edges, split it take off the last edge
			# and reassign the new edgebag to the suffix
			if(len(rindex) > 1):
				indexlist = rindex.split(',')
				new_edge_list = indexlist[1:len(indexlist)]
				rindex = indexlist[len(indexlist)-1]
				new_suffix_value = ','.join(new_edge_list)
				new_suffix_value = str(new_suffix_value)
				rindexint = int(rindex)
				# store suffix value
				start_suffix = self.ReadHash[rindexint][1]
				#if the prefix and the suffix are the same, try other edges first
				if(start == start_suffix):
					rindex = self.suffix_picker(start,1)
					indexlist = rindex.split(',')
					new_edge_list = indexlist[0:len(indexlist)-1]
					rindex = indexlist[len(indexlist)-1]
					new_suffix_value = ','.join(new_edge_list)
					new_suffix_value = str(new_suffix_value)
					rindexint = int(rindex)
					# store suffix value
					start_suffix = self.ReadHash[rindexint][1]
					self.GraphHash[start][start_suffix] = new_suffix_value
				else:
					self.GraphHash[start][start_suffix] = new_suffix_value
			# if the suffix only has one value store it in the start_suffix variable and delete the suffix
			elif(len(rindex) == 1):	
				rindexint = int(rindex)
				# store suffix value
				start_suffix = self.ReadHash[rindexint][1]
				self.GraphHash[start].pop(start_suffix, None)
			else:
				print "ERROR"

			#add prefix to edge set
			self.visited_nodes.append(self.ReadHash[rindexint][0])
			# adds prefix to current path traveled
			self.currentpath.append(self.AlphaReadHash[rindexint][0])
			# update Prefix and suffix in/out degree
			Prefix_out = self.DegreeHash[start][0]-1
			Prefix_in = self.DegreeHash[start][1]
			self.DegreeHash[start] = [Prefix_out, Prefix_in]
			Suffix_out = self.DegreeHash[start_suffix][0]
			Suffix_in = self.DegreeHash[start_suffix][1]-1
			self.DegreeHash[start_suffix] = [Suffix_out, Suffix_in]


			# if no edges are availale to be used from the first test, combine edges with same prefix 
			# and suffix into one list for storage in the maths member
			if(start==start_suffix):
				self.visited_nodes.append(self.ReadHash[rindexint][1])
				self.currentpath.append(self.AlphaReadHash[rindexint][1])
				addlist =  self.currentpath[:]
				self.paths.append(addlist) 
				del self.currentpath[:]
				del self.visited_nodes[:]
				return 1
			elif start_suffix in self.visited_nodes or start==start_suffix:
				addlist =  self.currentpath[:]
				self.paths.append(addlist) 
				del self.currentpath[:]
				del self.visited_nodes[:]
				return 1
			elif(self.DegreeHash[start_suffix][0] <= 0):
				self.visited_nodes.append(self.ReadHash[rindexint][1])
				self.currentpath.append(self.AlphaReadHash[rindexint][1])
				addlist =  self.currentpath[:]
				self.paths.append(addlist) 
				del self.currentpath[:]
				del self.visited_nodes[:]
				return 1
			self.find_paths(start_suffix)
		else:
			return -1

	def return_paths(self):
		return self.paths
	
	def remove_dups(self):
		return 0

	def print_Graph(self):
		print "Alpha       ",
		print self.AlphaReadHash
		print "Read Table: ",
		print self.ReadHash
		print "Hash: ",
		print self.GraphHash
		print "Degree Table",
		print self.DegreeHash

	def print_key(self, start):
		print "Key: "+ str(start), 
		print self.GraphHash[start]
		print self.GraphHash[start].values()

	def remove_dups(self):
		unique =  list()
		for contig in self.paths:
			if contig not in unique:
				unique.append(contig)
		self.paths = unique[:]
	def order_paths(self):
		maxcontig = list()
		for contig in self.paths:
			size = len(contig)
	def combine_paths(self):
		for contig in self.paths:
			self.contigs.append(''.join(contig))
		self.contigs = sorted(self.contigs, key=len)
		self.contigs.reverse()
		
	#checks each list item to see if string is contained within the list
	def remove_contained_contig(self):
		s = ''
		remove_list = []
		remove = 0
		for i in xrange(len(self.contigs)):
			if any(self.contigs[i] in s for s in self.contigs[i+1: len(self.contigs)]):
				remove = 1
			if any(self.contigs[i] in s for s in self.contigs[0:i]):
				remove = 1
			if remove == 1:
				remove_list.append(self.contigs[i])

			remove = 0
		for element in self.contigs:
			if (element in remove_list):
				self.contigs.remove(element)

	def print_contigs(self, filename):
		for contig in self.contigs:
			print contig
		with open("output.fasta", "w+") as file:
			for contig in self.contigs:
				file.write("{}\n".format(contig))

	def return_contig_len(self):
		return len(self.contigs)

filename = sys.argv[-1]

G = Path_Assembler(filename)
G.start_assembly()
G.remove_dups()
G.combine_paths()
G.print_contigs(filename)










