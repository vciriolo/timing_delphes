from sys import argv

inputfile = "../Cards/JetStudies_Phase_II_140PileUp_conf4.tcl"

if len(argv) > 1:
  inputfile = argv[1]

templatefile = inputfile + "_TEMPLATE"
formulafile = inputfile + "_FORMULAE"
  
writingeffformula = False
writingresformula = False
writingcalorimeter = False
effformdict = {}
resformdict = {}
etadict = {}
currentmodule = ""
calodict = {}

outfile = open(templatefile,'w')

for line in open(inputfile,'r').readlines():
  writeline = True
  if line.split("#")[0].count("module"):
    currentmodule = line.split()[2]
    effkeynum = 0
    reskeynum = 0
    etakeynum = 0
  if line.count("module Calorimeter"):
    key = "CALO_%s" % currentmodule
    outfile.write(key+"\n")
    writingcalorimeter = True
    calodict[key] = ""
  if writingcalorimeter:
    calodict[key] += line
    brackets += line.split("#")[0].count("{") - line.split("#")[0].count("}")
    if brackets == 0:
      writingcalorimeter = False
    writeline = False  
    continue
  if line.count("set TrackerEta") or line.count("set TauEta") or line.count("set PartonEta"):
    key = "ETA_%s_%i"%(currentmodule,etakeynum)
    outfile.write(key+"\n")
    etadict[key] = line
    writeline = False
    etakeynum += 1
  if line.split("#")[0].count("EfficiencyFormula {"):
    brackets = 0
    writingeffformula = True
    key= "EFF_%s_%i"%(currentmodule,effkeynum)
    effformdict[key] = ""
    outfile.write(key+"\n")
    effkeynum += 1
  if writingeffformula:
    effformdict[key] += line
    brackets += line.split("#")[0].count("{") - line.split("#")[0].count("}")
#    print brackets,line
    if brackets == 0:
      writingeffformula = False
    writeline = False
  if line.split("#")[0].count("ResolutionFormula {"):
    brackets = 0
    writingresformula = True
    key = "RES_%s_%i"%(currentmodule,reskeynum)
    resformdict[key] = ""
    outfile.write(key+"\n")
    reskeynum =+ 1
  if writingresformula:
    resformdict[key] += line
    brackets += line.split("#")[0].count("{") - line.split("#")[0].count("}")
    if brackets == 0:
      writingresformula = False
    writeline = False
  if writeline:
    outfile.write(line)  

outfile.close()

ffile =open(formulafile,'w')
ffile.write("etadict = %s\n" % etadict)
ffile.write("calodict = %s\n" % calodict)
ffile.write("effformdict = %s\n" % effformdict)
ffile.write("resformdict = %s\n" % resformdict)
ffile.close()

