#!/usr/bin/python
import os 
import re
from os import listdir
from os.path import isfile, join


def geticon(classname):
    if classname.startswith("BC"): return ":fas:`border-top-left;fa-fw` "
    if classname.startswith("IC"): return ":fas:`circle-right;fa-fw` "
    if classname.startswith("IO"): return ":fas:`print;fa-fw` "
    if classname.startswith("Integrator"): return ":fas:`gear;fa-fw` "
    if classname.startswith("Model"): return ":fas:`panorama;fa-fw` "
    if classname.startswith("Numeric"): return ":fas:`calculator;fa-fw` "
    if classname.startswith("Operator"): return ":far:`map;fa-fw` "
    if classname.startswith("Set"): return ":fas:`braille;fa-fw` "
    if classname.startswith("Solver"): return ":fas:`diamond-turn-right;fa-fw` "
    if classname.startswith("Util"): return ":fas:`sliders;fa-fw` "
    else: return ""

def getdocumentation(filename):
    sourcefile = open(filename+".H")
    ret = ""
    for line in sourcefile.readlines():
        if line.startswith(r"///"): # special provision for legacy doxygen comments
            ret += line.split(r"///")[1]
        elif line.startswith(r"// "):
            ret += line.split(r"// ")[1]
        elif line.startswith(r"//"):
            ret += line.split(r"//")[1]
        else:
            return ret
    return ret


def getParseDef(line):
    return "static void Parse" in line


def getParmParseDef(line):
    info = re.findall("ParmParse pp\(\"(.*)\"\)", line.split(r'//')[0])
    if not len(info): return None
    #if len(info) != 1: return None
    return info[0]
    
def getParmParseInfo(line):
    # Try "query" or "queryarr"
    info = re.findall(r'pp.(queryarr\b|query\b)\(\"(.*)\",(.*)\);', line.split(r'//')[0])
    if info: return info[0]

    # Try "queryclass"
    info = re.findall(r'pp.(queryclass\b)\(\"(.*)\",', line.split(r'//')[0])
    if info:
        return [info[0][0],info[0][1],None]

    return None


def getDocs(lines,i):
    ret = None
    for j in range(0,len(lines)-i):
        if j>0 and r'//' not in lines[i-j]: break
        if j>0 and getParseDef(lines[i-j]): break
        if j>0 and getParmParseDef(lines[i-j]): break
        if j>0 and getParmParseInfo(lines[i-j]): break
        #if j>0 and not lines[i-j].split(r'/')[0].isspace(): break
        docs = re.findall(r'.?\/\/(.*)', lines[i-j])
        if len(docs):
            if not ret: ret = ""
            if docs[0].endswith('/'): ret = docs[0] + ret
            else:                     ret = docs[0] + "\n" + ret
        elif j==0: continue
        else: break
    return ret


def extract(basefilename):
    rets = dict()
    for filename in [basefilename+".H",basefilename+".cpp"]:
        if not os.path.isfile(filename): continue
        sourcefile = open(filename)
        lines = sourcefile.readlines()
        prefix = None
        for i in range(len(lines)):
            try:
                if ("ParmParse pp" in lines[i].split(r'//')[0]):
                    prefix = getParmParseDef(lines[i])
                    if prefix in rets.keys(): continue
                    # Initialize the record for this ParmParser
                    rets[prefix] = dict()
                    rets[prefix]["items"] = []
                    docs = getDocs(lines,i)
                    if docs: rets[prefix]["docs"] = docs
                if (getParseDef(lines[i])):
                    prefix = "[prefix]"
                    rets[prefix] = dict()
                    rets[prefix]["items"] = []
                    docs = getDocs(lines,i)
                    if docs: rets[prefix]["docs"] = docs
                if ('pp.query' in lines[i]):
                    ret = dict()
                    docs = getDocs(lines,i)
                    if docs: ret["docs"] = docs
                    info = getParmParseInfo(lines[i])
                    if not info: continue
                    ret["query"], ret["string"], ret["variable"] = info
                    if prefix: 
                        ret["string"] = prefix + "." + ret["string"]
                    if prefix not in rets.keys(): 
                        rets[prefix] = dict()
                        rets[prefix]["items"] = []
                    rets[prefix]['items'].append(ret)
            except Exception as e:
                print("ERROR: reading file ", filename, " at line ", i)
                print("ERROR: ")
                print("ERROR:      ",lines[i])
                print("ERROR: Tried to read prefix, got prefix=",prefix)
                raise
    return rets

    
#    if classname.startswith("BC"): return ":icon:`border_outer` "
#    if classname.startswith("IC"): return ":icon:`start` "
#    if classname.startswith("Integrator"): return ":icon:`settings` "
#    if classname.startswith("Model"): return ":icon:`vrpano` "
#    if classname.startswith("Numeric"): return ":icon:`full_stacked_bar_chart` "
#    if classname.startswith("Util"): return ":icon:`settings` "
#    if classname.startswith("Solver"): return ":icon:`directions` "
#    if classname.startswith("IO"): return ":icon:`print` "
#    if classname.startswith("Operator"): return ":icon:`rebase_edit` "    

docfile    = open("Inputs.rst","w")
docfile.write(r"""
.. _inputs: 

=============================
:fas:`cube;fa-fw` Inputs
=============================


""")

headerchar = ["=","*","-","~","."]
written_headers = []

num_tot = 0
num_doc = 0

for dirname, subdirlist, filelist in sorted(os.walk("../../src/")):
    hdrname = dirname.replace("../../src/","").replace("/","::")
    depth = len(hdrname.split("::")) 

    srcfileset = set()
    for f in filelist:
        if f.endswith(".cpp"): srcfileset.add(f.replace(".cpp",""))
        if f.endswith(".H"): srcfileset.add(f.replace(".H",""))
    srcfilelist = list(srcfileset)
    
    #
    # This function makes sure pure abstract classes get
    # listed first.
    #
    def alphabetize_with_abstract_first(key):
        if key == hdrname.split("::")[-1]:
            return "0"
        return(key[0])
    for f in sorted(srcfilelist,key=alphabetize_with_abstract_first):
        
        if True: 
            try:
                inputs = extract(dirname+"/"+f)
            except Exception as e:
                print("ERROR: problem reading",dirname)
                raise
            documentation = getdocumentation(dirname+"/"+f)
            if not len(inputs) and not documentation:
                continue

            classname = dirname.replace("../../src/","").replace("/","::") + "::" + f.replace(".H","").replace(".cpp","")

            for i in range(len(classname.split('::'))):
                subhdr = '::'.join(classname.split('::')[:i])
                if subhdr not in written_headers:
                    if '::' not in subhdr and subhdr != "":
                        docfile.write("--------------------\n\n\n")
                        docfile.write(geticon(subhdr) + subhdr+"\n")
                        docfile.write("".ljust(len(geticon(subhdr)+subhdr),headerchar[i-1]))
                    else:
                        docfile.write("\n" + subhdr+"\n")
                        docfile.write("".ljust(len(subhdr),headerchar[i-1]))
                    docfile.write("\n\n\n")
                    written_headers.append(subhdr)

            if classname.split("::")[-1] != classname.split("::")[-2]:
                docfile.write(classname + "\n")
                lev = len(classname.split('::'))-1
                docfile.write("".ljust(len(classname),headerchar[lev])+"\n\n")
            
            if documentation:
                docfile.write(documentation)
            
            if not len(inputs): continue
            if len(inputs) == 1 and not list(inputs)[0]: continue


            docfile.write("\n\n")
            docfile.write(".. rst-class:: api-inputs-table\n\n")
            docfile.write(".. flat-table:: \n")
            docfile.write("    :widths: 20 10 70\n")
            docfile.write("    :header-rows: 1\n\n")
            docfile.write("    * - Parameter\n")
            docfile.write("      - Type\n")
            docfile.write("      - Description\n")

            for prefix in inputs:
                #if len(inputs[prefix]['items']) == 0: continue

                if ('docs' in inputs[prefix].keys()):
                    docfile.write("    * - {}  \n".format(inputs[prefix]['docs'].replace("\n", "\n        ")))

                for item in inputs[prefix]['items']:
                    num_tot += 1
                    docfile.write("    * - :code:`{}`\n".format(item['string']))
                    docfile.write("      - {}\n".format(item['query']))
                    if 'docs' in item.keys():
                        num_doc += 1
                        docfile.write("      - {}\n".format(item['docs'].replace('\n','\n        ')))
                    else:
                        docfile.write("      - \n")
                docfile.write("\n")
            docfile.write("\n")


print("\n{} of {} inputs documented\n".format(num_doc,num_tot))


docfile.close()


