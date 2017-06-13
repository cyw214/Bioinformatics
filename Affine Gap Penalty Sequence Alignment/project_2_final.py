#/usr/bin/python
import math
import os
import sys

# Gap start penalty
S = -10
# Gap extention penalty
E = -.5
# Match score
Match = 1
# Mistmatch score
Mismatch = -4

def readFasta(filename):
  try:
    afile = file(filename)
  except IOError, TypeError:                     
    print "Cannot Find File: " + filename
    sys.exit()
    return
  seqTitle = []
  sequences = {}
    
  for line in afile:
    if line.startswith('>'):
      name = line[1:].rstrip('\n')
      seqTitle.append(name)
      sequences[name] = ''
    else:
      sequences[name] += line.rstrip('\n').rstrip('*')
  return seqTitle, sequences

def changeSize(string, maxval):
  dif = maxval - len(string)
  if(len(string) < maxval):
    for i in range(dif):
      string = string + ' '
    return string
  else:
    return string


class Align:
  def __init__(self, filename):
    self.seqNames, self.sequences = readFasta(filename)
    self.contig1 = []
    self.contig2 = []
    self.name1 = ""
    self.name2 = ""
    self.rows = 1
    self.columns = 1
    self.M = [[0 for x in range(self.columns)] for x in range(self.rows)]
    self.X = [[0 for x in range(self.columns)] for x in range(self.rows)]
    self.Y = [[0 for x in range(self.columns)] for x in range(self.rows)]

  def printSeqNames(self):
    print "Option #\t" + "Sequence Description"
    print "--------\t" + "--------------------"
    print "(" + str(0) + ")\t\t"  + "EXIT"
    for i in xrange(len(self.seqNames)):
      print "(" + str(i+1) + ")\t\t"  + ">" + self.seqNames[i]

  def getSeqNamesLength(self):
    return len(self.seqNames)

  def SelectSequences(self, seq1, seq2):
    if(len(self.sequences.get(self.seqNames[seq1]))  >= len(self.sequences.get(self.seqNames[seq2]))):
      self.contig1 = list(self.sequences.get(self.seqNames[seq1]))
      self.name1 = self.seqNames[seq1]
      connum1 = seq1+1
      self.contig2 = list(self.sequences.get(self.seqNames[seq2]))
      self.name2 = self.seqNames[seq2]
      connum2 = seq2+1
    elif(len(self.sequences.get(self.seqNames[seq1])) < len(self.sequences.get(self.seqNames[seq2]))):
      self.contig1 = list(self.sequences.get(self.seqNames[seq2]))
      self.name1 = self.seqNames[seq2]
      connum1 = seq2+1
      self.contig2 = list(self.sequences.get(self.seqNames[seq1]))
      self.name2 = self.seqNames[seq1]
      connum2 = seq1+1
    else:
      print "Error: Invalid Sequence"
      connum1 = 0
      connum2 = 0     
    print "Selected: "
    print str(connum1) + ")" + self.seqNames[connum1-1]
    print str(connum2) + ")" + self.seqNames[connum2 -1]
    self.contig1.reverse()
    self.contig1.append('#')
    self.contig1.reverse()
    self.contig2.reverse()
    self.contig2.append('#')
    self.contig2.reverse()
    self.columns = len(self.contig1)
    self.rows = len(self.contig2)

  def Xscore(self, i, j):
    return max(S + E + self.M[i][j-1], E + self.X[i][j-1], S + E + self.Y[i][j-1])

  def Yscore(self, i, j):
    return max(S + E + self.M[i-1][j], S + E + self.X[i-1][j], E + self.Y[i-1][j])

  def Match(self, i, j):
    if(self.contig2[i] == self.contig1[j]):
      return 1
    else:
      return -4

  def Mscore(self, i, j):
    return self.Match(i, j) + max(self.M[i-1][j-1], self.X[i-1][j-1], self.Y[i-1][j-1])

  def InitMatricies(self):
    print "Scoring..."
    self.M = [[0 for x in range(self.columns)] for x in range(self.rows)]
    self.X = [[0 for x in range(self.columns)] for x in range(self.rows)]
    self.Y = [[0 for x in range(self.columns)] for x in range(self.rows)]
    self.X[0][1] = S + E
    self.Y[1][0] = S + E
    for i in xrange(self.columns):
      if(i > 0):
        self.M[0][i] = float('-infinity');
        self.Y[0][i] = float('-infinity');
      if(i > 1):
        self.X[0][i] = self.X[0][i-1] + E
      if(i > 0 and i < self.rows):
        self.M[i][0] = float('-infinity');
        self.X[i][0] = float('-infinity');
      if(i > 1 and i < self.rows):
        self.Y[i][0] = self.Y[i-1][0] + E

  def ScoreMatricies(self):
    for i in xrange(1, self.rows):
      for j in xrange(1, self.columns):
        self.X[i][j] = self.Xscore(i, j);
        self.Y[i][j] = self.Yscore(i, j);
        self.M[i][j] = self.Mscore(i, j);

  def backM(self, i, j):
    m = self.M[i-1][j-1] + self.Match(i,j)
    x = self.X[i-1][j-1]
    y = self.Y[i-1][j-1]
    if(max(m,x,y) == m):
      return ['m',i-1,j-1]
    elif(max(m,x,y) == x):
      return ['x',i-1,j-1]
    elif(max(m,x,y) == y):
      return ['y',i-1,j-1]

  def backX(self, i, j):
    m = self.M[i][j-1] + S + E
    x = self.X[i][j-1] + E
    y = self.Y[i][j-1] + S + E
    if(max(m,x,y) == m):
      return ['m',i,j-1]
    elif(max(m,x,y) == x):
      return ['x',i,j-1]
    elif(max(m,x,y) == y):
      return ['y',i,j-1]

  def backY(self, i, j):
    m = self.M[i-1][j] + S + E
    x = self.X[i-1][j] + S + E
    y = self.Y[i-1][j] + S
    if(max(m,x,y) == m):
      return ['m',i-1,j]
    elif(max(m,x,y) == x):
      return ['x',i-1,j]
    elif(max(m,x,y) == y):
      return ['y',i-1,j]

  def BackTrack(self):
    print "Back Tracking..."
    i = self.rows-1
    j = self.columns-1
    graph = ''
    count = 0
    if(max(self.M[i-1][j-1] + self.Match(i-1, j-1), self.X[i-1][j-1], self.Y[i-1][j-1]) == self.M[i-1][j-1] + self.Match(i-1, j-1)):
      graph = 'm'
    elif(max(self.M[i-1][j-1] + self.Match(i-1, j-1), self.X[i][j], self.Y[i][j]) == self.X[i][j]):
      graph = 'x'
    elif(max(self.M[i-1][j-1] + self.Match(i-1, j-1), self.Y[i][j]) == self.Y[i][j]):
      graph = 'y'
    while(i+j > 0):
      if(graph == 'm'):
        [graph, i, j] = self.backM(i, j)
      elif(graph == 'x'):
        self.contig2.insert(i+1, "_")
        [graph, i, j] = self.backX(i, j)
      elif(graph == 'y'):
        self.contig1.insert(j+1, "_")
        [graph, i, j] = self.backY(i, j)
      else:
        print "Failed"
        break
      count = count +1

  def outputAlignment(self):
    maxsize = max(len(self.name1),len(self.name2))
    name1 = changeSize(self.name1 , maxsize)
    name2 = changeSize(self.name2 , maxsize)
    matches = changeSize(" ", maxsize)
    if(len(name1) > 20):
      name1 = name1[0:20]
    if(len(name2) > 20):
      name2 = name2[0:20]
    if(len(matches) > 20):
      matches = matches[0:20]
    alignment1 = ''.join(self.contig1[1:len(self.contig1)])
    alignment2 = ''.join(self.contig2[1:len(self.contig2)])
    index1 = len(alignment1)
    index2 = len(alignment2)
    k = 0
    i = 0
    j = 0
    l = 0
    printflag1 = 1
    printflag2 = 1
    printflag3 = 1
    alignstring1 = ""
    alignstring2 = ""
    outputstring = ""
    print "(Append Output To Existing .fasta File or Create New Output File)"
    outname = raw_input("Enter Name For .fasta Output File: ")
    outname = outname + ".fasta"
    file = open(outname,"a+")   # Trying to create a new file or open one
    file.close()
    
    while(i < index1):
      while(i < index1):
        if(printflag1 == 1):
          print ">" + name1 +  " ",
          outputstring = outputstring + ">" + name1 + ": "
        printflag1 = 0      
        print alignment1[i],
        outputstring = outputstring + alignment1[i]
        i = i + 1
        if(i%60 == 0):
          printflag1 = 1  
          break
      print
      outputstring = outputstring + "\n"

      while(j < index2):
        if(printflag2 == 1):
          print ">" + name2 + " ",
          outputstring = outputstring + ">" + name2 + ": "
        printflag2 = 0  
        print alignment2[j],
        outputstring = outputstring + alignment2[j]
        j = j + 1
        if(j%60 == 0):
          printflag2 = 1  
          break
      print
      outputstring = outputstring + "\n"

      while(l < index2):
        if(printflag3 == 1):
          print " " + matches + " ",
          outputstring = outputstring + " " + matches + "  "
        printflag3 = 0  
        if(alignment1[l] == alignment2[l]):
            outputstring = outputstring + "*"
            print "*",
        else:
          outputstring = outputstring + " "
          print " ", 
        l = l + 1
        if(l%60 == 0):
          printflag3 = 1  
          print
          print
          break
      outputstring = outputstring + "\n\n"
    with open(outname, "a") as output:
      output.write(outputstring)

filename = sys.argv[-1]
A = Align(filename)
nameslen = A.getSeqNamesLength()
os.system('clear')  
while(1):
  # select first sequence to be aligned
  while(1):
    try:
      A.printSeqNames()
      print "Select First Sequence To Be Aligned"
      selection1 = raw_input("Option #: ")
      if(int(selection1) > nameslen or int(selection1) < 0):
        os.system('clear')  
        raise ValueError("INVALID. ENTER A NUMERIC OPTION BETWEEN 0 AND " + str(nameslen))
      if(int(selection1) <= nameslen):
        os.system('clear')
        break
    except ValueError:
      os.system('clear')  
      print "INVALID. ENTER A NUMERIC OPTION BETWEEN 0 AND " + str(nameslen)
  if(int(selection1) == 0):
    break
  os.system('clear')
  # select second sequence to be aligned
  while(1):
    try:
      A.printSeqNames()
      print "Select Second Sequence To Be Aligned"
      selection2 = raw_input("Option #: ")
      if(int(selection2) > nameslen or int(selection2) < 0):
        os.system('clear')  
        raise ValueError("INVALID. ENTER A NUMERIC OPTION BETWEEN 0 AND " + str(nameslen))
      if(int(selection2) <= nameslen):
        os.system('clear')
        break
    except ValueError:
      os.system('clear')  
      print "INVALID. ENTER A NUMERIC OPTION BETWEEN 0 AND " + str(nameslen)
  if(int(selection2) == 0):
    break

  A.SelectSequences(int(selection1)-1, int(selection2)-1)
  A.InitMatricies()
  A.ScoreMatricies()
  A.BackTrack()
  print "Alignment:"
  A.outputAlignment()
  print
  raw_input("Press Enter To Continue")
  os.system('clear')



















