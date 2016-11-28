import Bio
import sys
from Bio import SwissProt as sp
from  Bio import SeqIO
from urllib.parse import urlparse
import urllib
import re
from Bio.SubsMat.MatrixInfo import blosum62 as bs
import sys
import os
#Methods for web application
import tornado.ioloop
import tornado.web
import json
import tornado.httpserver



# Use secondary structure prediction and blosum 62 sum across a column

#Yields a list of topological annotations for each item in the text file 
class Residue:
    def __init__(self,annotation,aa):
        self.annotation = annotation
        self.aa = aa
#A data structure to represent the a sequence of Residues/indels
class Sequence:
    def __init__(self,identifier):
        self.sequence = [] #Linked list may have been better data structure, allows for insertion and deletion in constant time
        self.id = identifier
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self, item):
        return self.sequence[item]
    def __repr__(self):
        output = ""
        for item in self.sequence:
            output += item.aa
        return output


LINK = "http://www.uniprot.org/uniprot/{0}.txt"

def msaprocess(): #MSA must be provided in fasta format
    processed = SeqIO.parse("msa2.txt","fasta")
    output = [] #A list of sequence objects all annotated
    for record in processed:
        identifier = record.name
        #Parse the web record for a specific protein's SP annotation
        textrecord = urllib.request.urlopen(LINK.format(record.id))
        topology = process(textrecord)
        gapped_sequence = str(record.seq)
        output.append(generate_rep(gapped_sequence,topology,identifier))
    return output


def process(files): #the file is a text file of a swissprot protein
    parsed = sp.parse(files)
    record = next(parsed)
    sequence = record.sequence
    length = len(sequence)
    a = record.features
    finished = False
    output = []
    for item in a:
        if finished:
            break
        else:
            if(not(item[0] == "CHAIN") and not(item[0] == "DOMAIN")): #this requires more refinement to allow for a range of domains
                output.append([item[0],item[1],item[2]])
                if item[2] == length:
                    finished = True
    return output
def generate_rep(sequence, topology,identifier):
    #Initialize the classes to store sequence data and annotations
    overall = Sequence(identifier)

    #Keep track of all gaps and append positions to list
    gaps = []
    for i in range(len(sequence)):
        if (sequence[i] == "-"):
            gaps.append(i)

    #Remove gaps from sequence
    sequence = re.sub(r'-',"",sequence)

    #Initialize each protein without annotations
    for item in sequence:
        x = Residue("",item)
        overall.sequence.append(x)

    #Annotate those amino acids with Swissprot data
    for item in topology:
        start = item[1]
        end = item[2]
        notation = item[0]
        for i in range(start - 1, end):
            overall.sequence[i].annotation = notation


    #Reinsert gap characters into representation
    gap = Residue("","-")
    for index in gaps:
        overall.sequence.insert(index,gap)


    return overall

def scoreAcids(residue1, residue2):
    try:
        return bs[(residue1,residue2)]
    except:
        try:
            return bs[(residue2,residue1)]
        except:
            return 0 #Score of a gap character is set to 0


""" Indicates the position of the query in an MSA """
def indexer(lst, query):
    for i in range(len(lst)):
        if lst[i].id == query:
            return i


""" Takes a list of domains at a site and returns a boolean whether the
query is similar enough to the target at that point """
def similarity(hashmap, percentage):
    if (bool(hashmap)):
        mostcommon = max(hashmap,key = lambda x: hashmap[x])
        if (hashmap[mostcommon] / sum(hashmap.values()) >= percentage):
            return True
    return False

def mostcommon(hashmap):
    return max(hashmap,key = lambda x: hashmap[x])



def MDA(query, accuracy = 0.50): #Query is the target identifier
    #This is where scoring functions will be integrated with the matrix generated
    #from msaprocess()
    topology_matrix = msaprocess()

    #Identifies the number of changes that were made 
    changes = 0

    #identifies which chain you want to transfer the consensus annotation to
    target = indexer(topology_matrix,query)

    #Iterates through every residue in a sequence
    for sequence in range(len(topology_matrix[0])):

        #Collects the different domains in the column
        domains = {}

        #Collects the values of the different amino acids in the matrix
        acids = {}

        #Iterates through the number of columns
        for column in range(len(topology_matrix)):
            x = topology_matrix[column][sequence]
            if x.annotation in domains:
                domains[x.annotation] += 1
            else:
                domains[x.annotation] = 1
            if (x.aa,x.annotation) in acids:
                acids[(x.aa,x.annotation)] += 1
            else:
                acids[(x.aa,x.annotation)] = 1
        if similarity(domains,accuracy) and mostcommon(acids) != ('-', ''):
            targetDomain = mostcommon(domains)
            if blosumsop(acids,topology_matrix[target][sequence].aa,targetDomain) > 0:
                topology_matrix[target][sequence].annotation = targetDomain
                changes += 1
    print(str(changes) + " " + "Changes were Made")
    return topology_matrix[target]

""" Returns the blosum average score from comparing target annotation to query"""
def blosumsop(acids,targetAA,targetDomain):
    targetAcids = [x[0] for x in acids if x[1] == targetDomain]
    blosumscore = 0
    for residue in targetAcids:
        blosumscore += scoreAcids(targetAA,residue)
    return blosumscore



def printMDA(query):
    target = MDA(query)
    target.sequence = [x for x in target.sequence if x.aa != '-']
    overall = [0]
    count = 0
    for i in range(len(target)):
        if i == 0:
            start = [str(target[i].annotation), "start: " + str(i)]
            overall[count] = start
        elif target[i].annotation != target[i-1].annotation:
            overall[count].append("end: " + str(i-1))
            count += 1
            overall.append(0)
            start = [str(target[i].annotation), "start: " + str(i)]
            overall[count] = start

    overall[count].append("end: " + str(len(target)))
    #Prints the annotation if it's not empty or a dashed region 
    return [x for x in overall if x[0] != '']


"""The Framework for web routing (ignores the upload of the MSA file) """ 

class MainHandler(tornado.web.RequestHandler):
    def get(self):
        self.render("static/index.html")
    def post(self):
        query = self.get_argument("query").strip()
        result = printMDA(query)
        self.write("<h1 text-align='center'>The Toplogy for " + query + "</h1>")
        self.write('<table id="gradient-style" class="Experiments" align="center">')
        self.write('<tr><td>' + "Domain</td>" + "<td>Start</td>" + "<td>End</td></tr>")
        for item in result:
            self.write('<tr>')
            self.write('<td>'+item[0]+'</td><td>'+item[1][7:]+'</td><td>'+item[2][5:]+'</td>')
            self.write('</tr>')
        self.write('</table>')
        self.write("""<style> 

#gradient-style
{
    font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif,serif;
    font-size: 12px;
    margin: 45px;
    width: 480px;
    text-align: left;
    border-collapse: collapse;
}
#gradient-style th
{
    font-size: 13px;
    font-weight: normal;
    padding: 8px;
    background: #b9c9fe url('static/img/table-images/gradhead.png') repeat-x;
    border-top: 2px solid #d3ddff;
    border-bottom: 1px solid #fff;
    color: #039;
}
#gradient-style td
{
    padding: 8px;
    border-bottom: 1px solid #fff;
    color: #669;
    border-top: 1px solid #fff;
    background: #e8edff url('static/img/table-images/gradback.png') repeat-x;
}
#gradient-style tfoot tr td
{
    background: #e8edff;
    font-size: 12px;
    color: #99c;
}
#gradient-style tbody tr:hover td
{
    background: #d0dafd url('static/img/table-images/gradhover.png') repeat-x;
    color: #339;
}
</style>'""")


def make_app():
    return tornado.web.Application([
        (r"/",MainHandler)
    ])


if __name__ == "__main__":
    app = make_app()
    http_server = tornado.httpserver.HTTPServer(app)
    port = int(os.environ.get("PORT", 5000))
    http_server.listen(port) # hosts on localhost:5000
    print("Running Brainspell at http://localhost:5000...")
    tornado.ioloop.IOLoop.current().start()

""" End Web Framework """ 

"""
Conclusions: 
This project may have had further room for improvement by using a windowing technique to identify outliers in sets of data.
If there is a non-matching residue in a secondary structure, it must adopt the secondary structure of its neighbors. This approach
Can be conducted across columns and inserted as an additional pass to identify the true MDA of a protein (Pseudocode below). As it 
stands, a misfit in a domain will be shown in PrintMDA as a two identical domains each half the length of the true domain 
one after the other. The foundations for making this program into a web application is also provided. This would handle the 
difficulties of the extensive secondary structure prediction times, using asynchronous request handling. 


This assumes that secondary structure can be parsed as a list and the residue matrix is already annotated
def secondary structure(residue_matrix, secondary_structure_list):

    #Keep track of all gaps and append positions to list
    gaps = []
    for i in range(len(sequence)):
        if (sequence[i] == "-"):
            gaps.append(i)

    #Remove gaps from sequence
    sequence = re.sub(r'-',"",sequence)

    #Annotate those amino acids with secondary structure data
    for item in secondary_structure_list: #A list of type [secondary structure, start, end]
        start = item[1]
        end = item[2]
        notation = item[0]
        amino_acids = {} #Used to generate the most common amino acid 
        for i in range(start,end):
            overall.sequence[i].annotation.append(notation) #Adding a second annotation to each acid
            aa = overall.sequence[i].aa
            if aa in amino_acids:
                amino_acids[aa] += 1 
            else:
                amino_acids[aa] = 1
        maximum = mostcommon(acids) #This should represent a class of amino acids 

        #Ensure that the entire body the secondary structure is labeled the same 
        for i in range(start,end):
            overall.sequence[i].annotation = maximum #This assigns the value to a member of a class of amino acids that represents that secondary structure 


    #Reinsert gap characters into representation
    gap = Residue("","-")
    for index in gaps:
        overall.sequence.insert(index,gap)


    #Now iterate over a window k and ensure outliers are well matched 
    i = 0
    while (i +k) < len(residue_matrix[0]):



"""





