from sys import argv

templatefile = "../Cards/JetStudies_Phase_II_140PileUp_conf4.tcl_TEMPLATE"
formulafile = "../Cards/JetStudies_Phase_II_140PileUp_conf4.tcl_FORMULAE"
newfile = "test.tcl"

if len(argv) > 1 and not len(argv) == 4:
  raise Exception,"Wrong arguments"

if len(argv) > 3:
  templatefile = argv[1]
  formulafile = argv[2]
  newfile = argv[3]
  
writingeffformula = False
writingresformula = False
#effformdict = {}
#resformdict = {}
#etadict = {}
currentmodule = ""
keynum = -1

for line in open(formulafile,'r'):
  exec line.strip()

template = open(templatefile,'r').readlines()
outfile = open(newfile,'w')

effkeynum = 0
reskeynum = 0
etakeynum = 0

for line in template:
  writeline = True
  if line.split("#")[0].count("module"):
    currentmodule = line.split()[2]
    effkeynum = 0
    reskeynum = 0
    etakeynum = 0
#  print keynum,line.strip()
  etakey = "ETA_%s_%i"%(currentmodule,etakeynum)
  reskey = "RES_%s_%i"%(currentmodule,reskeynum)
  effkey = "EFF_%s_%i"%(currentmodule,effkeynum)
  if line.count("CALO_"):
    calokey = line.rstrip()
    outfile.write(calodict[calokey])
    writeline = False
  if line.count(etakey):
    outfile.write(etadict[etakey])
    writeline = False
    etakeynum += 1
  if line.count(reskey):
    outfile.write(resformdict[reskey])
    writeline = False
    reskeynum += 1
  if line.count(effkey):
    outfile.write(effformdict[effkey])
    writeline = False
    effkeynum += 1
  if writeline:
    outfile.write(line)  

outfile.close()
