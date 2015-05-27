# -*- coding: utf-8 -*-
# PyKnot: A PyMOL plugin for discovery and analysis of knots in protein structures
# Script/plugin by Rhonald Lua (rhonald.lua@gmail.com)
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

#Possible locations of the PyMOL plugins directory:
#On Windows, this is usually:
#   C:\Program Files\DeLano Scientific\PyMOL\modules\pmg_tk\startup
#On Macintosh, with the X11/Hybrid version, the location is probably:
#   PyMOLX11Hybrid.app/pymol/modules/pmg_tk/startup
#On a linux machine, it might be under:
#   /var/lib/python-support/python2.4/pmg_tk/startup/

from pymol import cmd #This launches pymol
from pymol.cgo import *
from pymol.vfont import plain
from Tkinter import *
import Tkinter,Pmw,math,random
from numpy import linalg, zeros

def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Launch PyKnot',
                             label='PyKnot...',
                             command = lambda s=self: PyKnotTools(s))

class KAControlGroup:
    MINSHAPE=0.52 #threshold used in reducing backbones and finding the knotted core in the statistical/heuristic approach
    BACKBONE_RADIUS=1
    CROSSING_RADIUS=1
    WARNING_MISSINGRESIDUES="(WARNING:ARTIFICIAL BONDS ARE PRESENT)"
    knottable=[("UNKNOTTED (0_1)",0,1.0,0,0),
               ("TREFOIL (3_1)",3,3.0,1,1),
               ("FIGURE EIGHT (4_1)",4,5.0,-1,0),
               ("PENTAFOIL (5_1)",5,5.0,3,5),
               ("THREE-TWIST KNOT (5_2)",5,7.0,2,3),
               ("STEVEDOREs KNOT (6_1)",6,9.0,-2,1),
               ("KNOTTED (6_2)",6,11.0,-1,1),
               ("KNOTTED (6_3)",6,13.0,1,0),
               ("SEPTAFOIL (7_1)",7,7.0,6,14),
               ("KNOTTED (7_2)",7,11.0,3,6),
               ("KNOTTED (7_3)",7,13.0,5,11),
               ("KNOTTED (7_4)",7,15.0,4,8),
               ("KNOTTED (7_5)",7,17.0,4,8),
               ("KNOTTED (7_6)",7,19.0,1,2),
               ("KNOTTED (7_7)",7,21.0,-1,1)]
        
    def __init__(self,
                 page,
                 groupname='Knot Analyzer',
                 defaultstructurename='(pdb)',
                 defaultclosureoption=0,
                 defaultchain=''):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group
        #self.groupscrolled=group
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.structureframe=Tkinter.Frame(group.interior())
        #Field for entering name of structure or model
        self.structure = Pmw.EntryField(self.structureframe,
                                        labelpos='w',
                                        label_text='structure to use: ',
                                        value=defaultstructurename,
                                        )

        self.chain = Pmw.EntryField(self.structureframe,
                                    labelpos='w',
                                    label_text='chain: ',
                                    value=defaultchain,
                                    )
        self.structure.grid(column=0,row=0)
        self.chain.grid(column=1,row=0)

        self.analyze_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.analyze_buttonbox.add('Analyze knot',command=self.analyzeKnot)

        self.test_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.test_buttonbox.add('Run tests',command=self.runTests)

        self.advancedgroup = Pmw.Group(group.interior(),tag_text='Fine-tune your analysis')

        self.residuenumbersframe=Tkinter.Frame(self.advancedgroup.interior())
        #Fields for entering starting and ending residue number in chain
        self.startres = Pmw.EntryField(self.residuenumbersframe,
                                        labelpos='w',
                                        label_text='starting residue number: ',
                                        #value='',
                                        )

        self.endres = Pmw.EntryField(self.residuenumbersframe,
                                    labelpos='w',
                                    label_text='ending residue number: ',
                                    #value='',
                                    )
        self.startres.grid(column=0,row=0)
        self.endres.grid(column=1,row=0)

        self.closure_options_tuple = ('Outward connection of terminals',
                                      'Direct connection of terminals',
                                      'Custom extension from terminals (enter vectors below)')
        self.show_closure_options = Pmw.OptionMenu(self.advancedgroup.interior(),
                                              labelpos = 'w',
                                              label_text = 'Knot closure method',
                                              items = self.closure_options_tuple,
                                              initialitem = self.closure_options_tuple[defaultclosureoption],
                                              #command=self.analyzeKnot
                                              )
        ########## Bioinformatics Reviewer 4 request
        self.customvectorsgroup = Pmw.Group(self.advancedgroup.interior(),tag_text='Custom extension from terminals')

        self.customvectorsframe=Tkinter.Frame(self.customvectorsgroup.interior())
        #Fields for entering starting and ending residue number in chain
        self.Nx = Pmw.EntryField(self.customvectorsframe,
                                 labelpos='w',
                                 label_text='N-terminal vector: x',
                                 validate={'validator':'real'},
                                 value='1'
                                 )
        self.Ny = Pmw.EntryField(self.customvectorsframe,
                                 labelpos='w',
                                 label_text='y',
                                 validate={'validator':'real'},
                                 value='0',
                                 )
        self.Nz = Pmw.EntryField(self.customvectorsframe,
                                 labelpos='w',
                                 label_text='z',
                                 validate={'validator':'real'},
                                 value='0',
                                 )
        self.Cx = Pmw.EntryField(self.customvectorsframe,
                                 labelpos='w',
                                 label_text='C-terminal vector: x',
                                 validate={'validator':'real'},
                                 value='-1',
                                 )
        self.Cy = Pmw.EntryField(self.customvectorsframe,
                                 labelpos='w',
                                 label_text='y',
                                 validate={'validator':'real'},
                                 value='0',
                                 )
        self.Cz = Pmw.EntryField(self.customvectorsframe,
                                 labelpos='w',
                                 label_text='z',
                                 validate={'validator':'real'},
                                 value='0',
                                 )

        self.Nx.grid(column=0,row=0)
        self.Ny.grid(column=1,row=0)
        self.Nz.grid(column=2,row=0)
        self.Cx.grid(column=0,row=1)
        self.Cy.grid(column=1,row=1)
        self.Cz.grid(column=2,row=1)
        self.customvectorsframe.pack(fill='x',padx=4,pady=1) # vertical
        #########
        self.backboneatom_options_tuple = ('CA (Protein)',
                                            'P (DNA/RNA)')
        self.show_backboneatom_options = Pmw.OptionMenu(self.advancedgroup.interior(),
                                              labelpos = 'w',
                                              label_text = 'Use as backbone atom',
                                              items = self.backboneatom_options_tuple,
                                              initialitem = self.backboneatom_options_tuple[0],
                                              #command=self.updateInterface #Need two arguments for this method
                                              )
        
        self.knotlabel_var=IntVar()
        self.knotlabel_var.set(1)
        self.knotlabel_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Show label indicating the knot type suggested by the knot invariants",
                                                     variable = self.knotlabel_var)

        self.gausscodeframe=Tkinter.Frame(self.advancedgroup.interior())
        self.gausscode_var=IntVar()
        self.gausscode_var.set(1)
        self.gausscode_checkbutton = Checkbutton(self.gausscodeframe,
                                                     text = "Print the equivalent gauss code of the knot",
                                                     variable = self.gausscode_var)
        self.gausscode_options_tuple = ('After reduction by Reidemeister moves',
                                            'Before reduction by Reidemeister moves')
        self.show_gausscode_options = Pmw.OptionMenu(self.gausscodeframe,
                                              labelpos = 'w',
                                              label_text = '',
                                              items = self.gausscode_options_tuple,
                                              initialitem = self.gausscode_options_tuple[0]
                                              )
        self.gausscode_checkbutton.grid(column=0,row=0)
        self.show_gausscode_options.grid(column=1,row=0)

        self.reidemeister_var=IntVar()
        self.reidemeister_var.set(1)
        self.reidemeister_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Show progress of Reidemeister moves in reducing crossings in the knot projection",
                                                     variable = self.reidemeister_var)

        self.altloc_var=IntVar()
        self.altloc_var.set(0)
        self.altloc_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Use first atom variant when alternate location indicators exist",
                                                     variable = self.altloc_var)

        self.ss_var=IntVar()
        self.ss_var.set(0)
        self.ss_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Use PyMOL secondary structure assignments to reduce crossings in the knot projection",
                                                     variable = self.ss_var)

        self.reducebackboneframe=Tkinter.Frame(self.advancedgroup.interior())
        self.reducebackbone_var=IntVar()
        self.reducebackbone_var.set(0)
        self.reducebackbone_checkbutton = Checkbutton(self.reducebackboneframe,
                                                     text = "Reduce backbone",
                                                     variable = self.reducebackbone_var)
        self.reducebackbone_options_tuple = ('Tight corners first',
                                            'Random')
        self.show_reducebackbone_options = Pmw.OptionMenu(self.reducebackboneframe,
                                              labelpos = 'w',
                                              label_text = '',
                                              items = self.reducebackbone_options_tuple,
                                              initialitem = self.reducebackbone_options_tuple[0]
                                              )
        self.reducebackbone_checkbutton.grid(column=0,row=0)
        self.show_reducebackbone_options.grid(column=1,row=0)

        self.crossings_var=IntVar()
        self.crossings_var.set(0)
        self.crossings_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Mark Reidemeister-reduced knot crossings in the x-y plane projection of the structure",
                                                     variable = self.crossings_var)

        self.printbackbone_var=IntVar() #Also motivated by Bioinformatics Reviewer 4
        self.printbackbone_var.set(0)
        self.printbackbone_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Print the coordinates of the backbone",
                                                     variable = self.printbackbone_var)

        self.findknotcoreframe=Tkinter.Frame(self.advancedgroup.interior())
        self.findknotcore_var=IntVar() #Also motivated by Bioinformatics Reviewer 4
        self.findknotcore_var.set(0)
        self.findknotcore_checkbutton = Checkbutton(self.findknotcoreframe,
                                                     text = "Estimate the location of the knotted core.",
                                                     variable = self.findknotcore_var)
        self.findknotcoresamples = Pmw.EntryField(self.findknotcoreframe,
                                 labelpos='w',
                                 label_text='Number of samples:',
                                 validate={'validator':'integer','min':1,'max':100},
                                 value='10',
                                 )
        self.findknotcore_checkbutton.grid(column=0,row=0)
        self.findknotcoresamples.grid(column=1,row=0)

        self.residuenumbersframe.pack(fill='x',padx=4,pady=1)
        self.show_closure_options.pack(fill='x',padx=4,pady=1)
        self.customvectorsgroup.pack(fill='x',padx=4,pady=1)
        self.show_backboneatom_options.pack(fill='x',padx=4,pady=1)
        self.knotlabel_checkbutton.pack(fill='x',padx=4,pady=1)
        #self.gausscode_checkbutton.pack(fill='x',padx=4,pady=1)
        self.gausscodeframe.pack(fill='x',padx=4,pady=1)
        self.crossings_checkbutton.pack(fill='x',padx=4,pady=1)
        self.reidemeister_checkbutton.pack(fill='x',padx=4,pady=1)
        self.reducebackboneframe.pack(fill='x',padx=4,pady=1)
        self.altloc_checkbutton.pack(fill='x',padx=4,pady=1)
        self.printbackbone_checkbutton.pack(fill='x',padx=4,pady=1)
        self.findknotcoreframe.pack(fill='x',padx=4,pady=1)
        #self.ss_checkbutton.pack(fill='x',padx=4,pady=1)
        
        for entry in (self.structureframe,
                      self.analyze_buttonbox,
                      #self.test_buttonbox,
                      self.advancedgroup):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        #Tool tips (help) for user-interface components and boxes
        self.balloon.bind(self.structure,"Enter a name from the list of objects (PDB structures) in the PyMOL Viewer")
        #self.balloon.bind(self.chain,"Enter a chain identifier (A,B,etc.). Leave blank to use all chains in the structure (not recommended).")
        self.balloon.bind(self.chain,"Enter a chain identifier (A,B,etc.)")
        self.balloon.bind(self.analyze_buttonbox,"Click here to perform knot analysis on your structure. A possible knot type will be reported")
        self.balloon.bind(self.startres,"The backbone begins at this residue number (optional)")
        self.balloon.bind(self.endres,"The backbone ends at this residue number (optional)")
        self.balloon.bind(self.show_closure_options,"Choose how the terminals are connected for the backbone to form a closed loop")
        self.balloon.bind(self.show_backboneatom_options,"Choose an atom type to define the backbone")
        self.balloon.bind(self.knotlabel_checkbutton,"Check to indicate the knot type as a label at the C-terminal end of the backbone")
        self.balloon.bind(self.gausscode_checkbutton,"Check to print a gauss code for the knot (not a unique representation of the knot)")
        self.balloon.bind(self.crossings_checkbutton,"Check to show the backbone crossings in the x-y projection as cylinders")
        self.balloon.bind(self.reidemeister_checkbutton,"Check to print the number of knot crossings remaining after each application of Reidemeister moves")
        self.balloon.bind(self.reducebackbone_checkbutton,"Check to reduce the number of backbone atoms by removing bends and kinks")
        self.balloon.bind(self.altloc_checkbutton,"Check to use the first atom variant when alternate location indicators exist in the PDB")
        self.balloon.bind(self.printbackbone_checkbutton,"Print the x,y,z coordinates of the backbone, followed by the coordinates of the closing segments")
        self.balloon.bind(self.Nx,"Enter the vector components of the direction of the N-terminal extension to the surface of the bounding sphere")
        self.balloon.bind(self.Cx,"Enter the vector components of the direction of the C-terminal extension to the surface of the bounding sphere")
        self.balloon.bind(self.findknotcore_checkbutton,"Estimate the start and end of the knot's location in the sequence by sampling random reduced backbones")
        self.balloon.bind(self.findknotcoresamples,"Enter the number of random reduced backbones to sample (from 1 to 100)")
##The proteins with trefoil knots (31) are: 1lugA, 1v2xA, 1o6dA,
##1mxiA, 1ualA, 1vhyA, 1t0hB, 1js1X, 1k3rB, 1x7oA, 1p7lA, 1gz0E,
##1vhkD, 1gkuB, 1xi4C.
##The proteins with figure-eight knots (41) are: 1qmgA, 1u2zC,  
##1m72B.
##The protein with the knot 52: 1xd3A. (01/18/12 Delta=21 for chain A. Why? j-i!=N-1 in getCrossings fixes this. This was
## actually caused by a subtle error in Alexander matrix computation)

##Additional: 2efv (Trefoil, 3_1), 3bjx (Stevedore, 6_1), 1ns5 (trefoil)

##KNOT Alexander, Vassiliev, Vassiliev, CHIRAL?
##|Delta(−1)K| v2(K) |v3(K)|
##0_1 1 0 0 NO
##3_1 3 1 1 YES
##4_1 5 -1 0 NO
##5_1 5 3 5 YES
##5_2 7 2 3 YES
##6_1 9 -2 1
##6_2 11 -1 1
##6_3 13 1 0

#1xd3A and other PDB have lines with alternate location indicators like
#ATOM   1911  CA ASER A 228      17.315  41.317  34.251  0.50  8.25           C  
#ATOM   1912  CA BSER A 228      17.315  41.236  34.289  0.50 10.88           C

    #This is for development only
    def runTests(self):
        numtests=0
        numverified=0

        for gausscode, invariants in [("b-1,a-2,b-3,a-1,b-2,a-3",(3.0,1,1)),
            ("a+1,b+2,a+3,b+1,a+2,b+3",(3.0,1,1)),
            ("b+1,a+1,a+2,a-3,a+4,b+5,a+6,a-7,a+8,b+2,b-3,b+4,a+5,b+6,b-7,b+8",(3.0,1,1)),
            ("a-1,b-2,a-3,b-4,a-5,b-6,a-7,b-8,a-9,b-1,a-2,b-3,a-4,b-5,a-6,b-7,a-8,b-9", (9,10,30)),
            ("a+1,b-2,a-3,b-4,a-5,a-6,b-7,b-8,a-2,a-9,b-10,b-3,a-4,b-5,a-8,b+1,b-9,a-10,b-6,a-7", (5,7,18)), #PERKO PAIR
            ("b-1,a+2,b+3,a+4,b+5,b+6,a+7,a+8,b+2,b+9,a+10,a+3,b+4,a+5,b+8,a-1,a+9,b+10,a+6,b+7", (5,7,18))]:
            numtests+=1
            l=gausscode.split(',')
            ablist=[]
            signedlabellist=[]
            for c in l:
                ablist.append(c[0])
                signedlabellist.append(int(c[1:]))
                
            if True:
                print "Verifying knot invariants for %s..." % gausscode
                if len(ablist)>5:
                    Delta=self.computeAlexander(ablist,signedlabellist)
                    (v2,arf)=self.computeVassiliev2(ablist,signedlabellist)
                    v3=self.computeVassiliev3(ablist,signedlabellist)
                else:
                    Delta=1.0
                    v2=0
                    v3=0
                if math.fabs(Delta-invariants[0])<0.000001 and \
                   v2==invariants[1] and math.fabs(v3)==invariants[2]:
                    print "...Verified"
                    numverified+=1
                else:
                    print "...Invariants do not match!!!!!!!!!!!!!!!!!!!!", (Delta, v2, "|%d|" % v3), "!=", invariants

        for pdbc, invariants in [('2phyA',(1.0,0,0)),
                                 ('2efvA',(3.0,1,1)),
                                 ('1lugA',(3.0,1,1)),
                                 ('1yveI',(5.0,-1,0)),
                                 ('1qmgA',(5.0,-1,0)),
                                 ('1xd3A',(7.0,2,3)),
                                 ('3bjxA',(9.0,-2,1))]:
            
            numtests+=1
            #Use ET Server to get the structures
            cmd.load('http://mammoth.bcm.tmc.edu/ETserver2/pdbeasytrace/%s.etvx' % pdbc, pdbc)
            self.structurename=pdbc
            self.chainindicator=''

            structurename=pdbc
            try:
                model=cmd.get_model(structurename)
            except Exception:
                print "Structure \'%s\' does not exist!" % structurename
                return

            #TODO: Restrict to amino acids
            if self.getBackbone(model,'CA'):
                self.closeBackbone()
                try:
                    (ablist,signedlabellist)=self.getCrossings()
                except Exception, e:
                    print "Error (Exception) encountered in processing knot crossings!"
                    print e
                
                cmd.hide("lines",structurename)
                cmd.show("ribbon",structurename)
                #cmd.zoom("all")

                print "Verifying knot invariants for %s..." % structurename
                if len(ablist)>5:
                    Delta=self.computeAlexander(ablist,signedlabellist)
                    (v2,arf)=self.computeVassiliev2(ablist,signedlabellist)
                    v3=self.computeVassiliev3(ablist,signedlabellist)
                else:
                    Delta=1.0
                    v2=0
                    v3=0
                if math.fabs(Delta-invariants[0])<0.000001 and \
                   v2==invariants[1] and math.fabs(v3)==invariants[2]:
                    print "...Verified"
                    numverified+=1
                else:
                    print "...Invariants do not match!!!!!!!!!!!!!!!!!!!!", (Delta, v2, v3), "!=", invariants
            else:
                Delta=1.0
                v2=0
                v3=0
                if math.fabs(Delta-invariants[0])<0.000001 and \
                   v2==invariants[1] and math.fabs(v3)==invariants[2]:
                    print "...Verified"
                    numverified+=1
                else:
                    print "...Invariants do not match!!!!!!!!!!!!!!!!!!!!", (Delta, v2, v3), "!=", invariants
                print 'Unexpected error!!!!!!!!!!!!!!!!!!! The backbone is too short and therefore unknotted!'

        if numtests==numverified:
            print "PASSED ALL (%d) TESTS" % numtests
        else:
            print "%d DID NOT PASS THE TESTS!!!!!!!!!!" % (numtests-numverified)

    def analyzeKnot(self):
        #Get user-provided PyMOL structure selection.
        #Get structurename and check if it exists in PyMOL.
        structurename=self.structure.getvalue().strip()
        supper=structurename.upper()
        for n in cmd.get_names('all'):
            if supper==n.upper():
                break
        else:
            print 'structure name must be in PyMOL viewer list:'
            print cmd.get_names('all')
            return

        if len(structurename)<1:
            print 'Provide structure name!'
            return
        try:
            model=cmd.get_model(structurename)
        except Exception:
            print "Structure \'%s\' does not exist!" % structurename
            return
        
        chainindicator=self.chain.getvalue().strip() #TODO add chain to structurename?
        startresnum=self.startres.getvalue().strip()
        endresnum=self.endres.getvalue().strip()

        self.structurename=structurename
        if len(chainindicator)>0:
            self.chainindicator=chainindicator[0]
        else:
            self.chainindicator=''

        if len(self.chainindicator)>0:
            if self.chainindicator!='?':
                structurename='%s and chain %s' % (self.structurename,self.chainindicator)
            else:
                structurename=self.structurename
        else:
            #Chain is blank, check if there are multiple chains in the structure
            structurename=self.structurename
##            chaindict={}
##            try:
##                model=cmd.get_model(structurename)
##            except Exception:
##                print "Structure \'%s\' does not exist!" % structurename
##                return
##            for atom in model.atom:
##                chaindict[atom.chain]=1
##            if len(chaindict)>1:
            chainlist=cmd.get_chains(structurename)
            if len(chainlist)>1:
                print "Please select one from the following chains of %s:" % structurename, ','.join(chainlist)
                return

        if len(startresnum)>0 or len(endresnum)>0:
            structurename+=' and resi %s-%s' % (startresnum,endresnum)
        try:
            model=cmd.get_model(structurename)
        except Exception:
            print "Structure \'%s\' does not exist!" % structurename
            return

        print "BEGIN KNOT ANALYSIS OF %s" % structurename

        #TODO: Restrict to amino acids
        backboneatom_option=self.show_backboneatom_options.getvalue()
        backboneatomname='CA'
        if backboneatom_option==self.backboneatom_options_tuple[0]:
            backboneatomname='CA'
        elif backboneatom_option==self.backboneatom_options_tuple[1]:
            backboneatomname='P'
        if self.getBackbone(model,backboneatomname):
            if self.printbackbone_var.get()==1:
                #Print backbone coordinates (Requested by Bioinformatics Reviewer 4)
                print "Printing backbone (%s-atom) coordinates:" % backboneatomname
                print "   x       y       z"   
                i=0
                for x,y,z in self.backbone:
                    #(x,y,z)=self.backbone[i]
                    if i==0:
                        print "%7.3f %7.3f %7.3f (N-terminal)" % (x,y,z)
                    elif i==len(self.backbone)-1:
                        print "%7.3f %7.3f %7.3f (C-terminal)" % (x,y,z)
                    else:
                        print "%7.3f %7.3f %7.3f" % (x,y,z)
                    i+=1
                #print 'N-terminal C-alpha coordinates:', self.backbone[0]
                #print 'C-terminal C-alpha coordinates:', self.backbone[-1]
            self.closeBackbone()
            try:
                (ablist,signedlabellist)=self.getCrossings()
            except Exception, e:
                print "Error (Exception) encountered in processing knot crossings!"
                print e
            
            cmd.hide("lines",structurename)
            cmd.show("ribbon",structurename)
            cmd.zoom("all")

            if self.gausscode_var.get()==1:
                if len(self.gausscode)>0:
                    print '***** Gauss code for %s *****' % structurename
                    print ','.join(self.gausscode)

            if len(ablist)>5:
                if False: #self.gausscode_var.get()==1:
                    gausscode=[None]*len(ablist)
                    for i in xrange(len(ablist)):
                        slabel=signedlabellist[i]
                        if slabel>0:
                            sign='+'
                        else:
                            slabel=-slabel
                            sign='-'
                        gausscode[i]='%s%s%d' % (ablist[i],sign,slabel)
                    print '***** Gauss code for %s *****' % structurename
                    print ','.join(gausscode)

                print '***** Alexander invariant for %s *****' % structurename
                Delta=None
                try:
                    Delta=self.computeAlexander(ablist,signedlabellist)
                    print '|Delta(-1)|=', Delta
                except Exception, e:
                    print "Error (Exception) encountered in computing the Alexander invariant!"
                    print e

                print '***** Vassiliev invariants for %s *****' % structurename
                try:
                    (v2,arf)=self.computeVassiliev2(ablist,signedlabellist)
                    print 'degree 2 (v2) and arf:', v2, arf
                    v3=self.computeVassiliev3(ablist,signedlabellist)
                    print 'degree 3 (v3):', v3
                except Exception, e:
                    print "Error (Exception) encountered in computing vassiliev invariants!"
                    print e

                #Print terminals of knotted core
                if (not (Delta==1 and v2==0 and v3==0)) and self.findknotcore_var.get()==1:
                    try:
                        numsamples=int(self.findknotcoresamples.getvalue())
                    except Exception:
                        numsamples=10
                    print "Setting the number of samples of the random reduced backbones to %d" % numsamples
                    self.findKnotCore(numsamples,backboneatomname)
                
                foundmatch=False
                #TODO How to handle negative resi? It is interpreted as a range of residue numbers by PyMOL
                #use lastresnum instead
                if self.knotlabel_var.get()==1 and self.lastresnum is not None:
                    if (Delta is not None) and (v2 is not None) and (v3 is not None):
                        for textk, crossingk, Dk, v2k, v3k in self.knottable:
                            if math.fabs(Delta-Dk)<0.000001 and v2==v2k and math.fabs(v3)==v3k:
                                if len(ablist)/2>crossingk: #TODO This is too strict.
                                    print "The knot type is possibly %s" % textk
                                    print "The upper bound on the number of crossings is %d" % (len(ablist)/2)
                                    if self.break_resi:
                                        print self.WARNING_MISSINGRESIDUES
                                    textk+="?" #Indicate uncertainty of knot type
                                else:
                                    print "The knot type is %s with %d crossings" % (textk,len(ablist)/2)
                                    if self.break_resi:
                                        print self.WARNING_MISSINGRESIDUES
                                if self.break_resi:
                                    textk+=self.WARNING_MISSINGRESIDUES
                                #It is ok if self.chainindicator is blank
                                cmd.label("%s and chain %s and resi %s and name %s" % (self.structurename,
                                                                             self.chainindicator,
                                                                             self.lastresnum,
                                                                               backboneatomname),
                                        "'%s'" % textk) #Pay attention to the quotes within quotes
                                foundmatch=True
                                break
                        if not foundmatch:
                            cmd.label("%s and chain %s and resi %s and name %s" % (self.structurename,
                                                                         self.chainindicator,
                                                                         self.lastresnum,
                                                                           backboneatomname),
                                    "'KNOTTED'") #Pay attention to the quotes within quotes
                            print "The backbone is knotted"
                            print "The upper bound on the number of crossings is %d" % len(ablist)/2
                            if self.break_resi:
                                print self.WARNING_MISSINGRESIDUES
            else:
                print "The backbone of %s is unknotted!" % structurename
                if self.break_resi:
                    print self.WARNING_MISSINGRESIDUES
                if self.knotlabel_var.get()==1 and self.lastresnum is not None:
                    textk=self.knottable[0][0]
                    if self.break_resi:
                        textk+=self.WARNING_MISSINGRESIDUES
                    cmd.label("%s and chain %s and resi %s and name %s" % (self.structurename,
                                                                         self.chainindicator,
                                                                         self.lastresnum,
                                                                           backboneatomname),
                              "'%s'" % textk) #Pay attention to the quotes within quotes
        else:
            self.break_resi=False #TODO No need to check break_resi in this block
            print 'The backbone of %s is too short and therefore unknotted!' % structurename
            if self.break_resi:
                print self.WARNING_MISSINGRESIDUES
            if self.knotlabel_var.get()==1 and self.lastresnum is not None:
                textk=self.knottable[0][0]
                if self.break_resi:
                    textk+=self.WARNING_MISSINGRESIDUES
                cmd.label("%s and chain %s and resi %s and name %s" % (self.structurename,
                                                                     self.chainindicator,
                                                                     self.lastresnum,
                                                                       backboneatomname),
                          "'%s'" % textk) #Pay attention to the quotes within quotes

        print "END KNOT ANALYSIS OF %s" % structurename

    #The backbone must be knotted to call this method
    #Warning: this modifies self.backbone
    def findKnotCore(self,numsamples,backboneatomname):
        i=0
        startcore=[]
        endcore=[]
        while i<numsamples:
            self.backbone=self.original_backbone[:]
            self.reduceBackbone_random(self.MINSHAPE)
            if len(self.backbone)>3:
                startcore.append(self.backbone[1]) #Collect approximate (candidate) starting residue for the core
                endcore.append(self.backbone[-2])  #Collect approximate ending residue for the core
            else:
                print "Warning: Unexpected reduced backbone sample in findKnotCore!"
            i+=1
        if len(startcore)<1:
            print "Warning: Unexpected result in findKnotCore!"
            return
        #Find the CM, then find the coordinate closest to the CM
        #This is PyKnot's estimate of the start (N-term) of the knotted core
        (xcm,ycm,zcm)=(0,0,0)
        for nv in startcore:
            xcm+=nv[0]
            ycm+=nv[1]
            zcm+=nv[2]
        xcm*=1.0/len(startcore)
        ycm*=1.0/len(startcore)
        zcm*=1.0/len(startcore)
        tmp=startcore[0]
        mind2=(tmp[0]-xcm)*(tmp[0]-xcm)+\
             (tmp[1]-ycm)*(tmp[1]-ycm)+\
             (tmp[2]-zcm)*(tmp[2]-zcm)
        startcore_best=tmp
        for nv in startcore:
            d2=(nv[0]-xcm)*(nv[0]-xcm)+(nv[1]-ycm)*(nv[1]-ycm)+(nv[2]-zcm)*(nv[2]-zcm)
            if d2<mind2:
                startcore_best=nv
                mind2=d2
        #Find the CM, then find the coordinate closest to the CM
        #This is PyKnot's estimate of the end (C-term) of the knotted core
        (xcm,ycm,zcm)=(0,0,0)
        for nv in endcore:
            xcm+=nv[0]
            ycm+=nv[1]
            zcm+=nv[2]
        xcm*=1.0/len(endcore)
        ycm*=1.0/len(endcore)
        zcm*=1.0/len(endcore)
        tmp=endcore[0]
        mind2=(tmp[0]-xcm)*(tmp[0]-xcm)+\
             (tmp[1]-ycm)*(tmp[1]-ycm)+\
             (tmp[2]-zcm)*(tmp[2]-zcm)
        endcore_best=tmp
        for nv in endcore:
            d2=(nv[0]-xcm)*(nv[0]-xcm)+(nv[1]-ycm)*(nv[1]-ycm)+(nv[2]-zcm)*(nv[2]-zcm)
            if d2<mind2:
                endcore_best=nv
                mind2=d2

        #Find the residue numbers corresponding to the estimated start and end of the core
        try:
            model=cmd.get_model(self.structurename)
        except Exception:
            print "Error in findKnotCore: Structure \'%s\' does not exist!" % self.structurename
            return

        tmp=startcore_best
        mind2=1000000
        for atom in model.atom:
            if atom.name==backboneatomname:
                c=atom.coord
                d2=(tmp[0]-c[0])*(tmp[0]-c[0])+\
                     (tmp[1]-c[1])*(tmp[1]-c[1])+\
                     (tmp[2]-c[2])*(tmp[2]-c[2])
                if d2<mind2:
                    startres=atom.resi
                    mind2=d2

        tmp=endcore_best
        mind2=1000000
        for atom in model.atom:
            if atom.name==backboneatomname:
                c=atom.coord
                d2=(tmp[0]-c[0])*(tmp[0]-c[0])+\
                     (tmp[1]-c[1])*(tmp[1]-c[1])+\
                     (tmp[2]-c[2])*(tmp[2]-c[2])
                if d2<mind2:
                    mind2=d2
                    endres=atom.resi

        if len(self.chainindicator)>0:
            structurename='%s and chain %s' % (self.structurename,self.chainindicator)
        else:
            structurename=self.structurename
        structurename_chain=self.structurename+self.chainindicator
        cmd.select("KNOTCORE_%s" % structurename_chain, "%s and resi %s-%s" % (structurename,startres,endres))

        print "Estimated starting residue number of the knotted core:", startres
        print "Estimated ending residue number of the knotted core:", endres
        
        #TODO cone or cylinder. Need coordinates of adjacent backbone atom then
        #obj=[SPHERE, startcore_best[0], startcore_best[1], startcore_best[2], self.BACKBONE_RADIUS]
##        obj=[SPHERE, startcore_best[0], startcore_best[1], startcore_best[2], 2*self.BACKBONE_RADIUS, 0.0, 1.0, 0.0]
##        structurename_chain=self.structurename+self.chainindicator
##        cmd.delete("STARTKNOT_%s" % structurename_chain)
##        cmd.load_cgo(obj,"STARTKNOT_%s" % structurename_chain)
##
##        obj=[SPHERE, endcore_best[0], endcore_best[1], endcore_best[2], 2*self.BACKBONE_RADIUS, 1.0, 0.0, 0.0]
##        structurename_chain=self.structurename+self.chainindicator
##        cmd.delete("ENDKNOT_%s" % structurename_chain)
##        cmd.load_cgo(obj,"ENDKNOT_%s" % structurename_chain)

    #Get coordinates of backbone
    def getBackbone(self,model,backboneatomname):
        use_firstaltloc=self.altloc_var.get()
        use_reducebackbone=self.reducebackbone_var.get()
        numCA=0
        firstresnum=None
        lastresnum=None
        for atom in model.atom:
            if atom.name==backboneatomname: #Use alpha-carbon for protein. For DNA, use 'P'.
                #if atom.__dict__.has_key("alt") and \
                if use_firstaltloc==1 and \
                   atom.alt!='' and atom.alt!='A': #Accept the first alternate location for the atom and ignore others
                    continue
                if firstresnum is None:
                    firstresnum=atom.resi
                lastresnum=atom.resi
                numCA+=1
        self.firstresnum=firstresnum
        self.lastresnum=lastresnum

        if numCA>4:
            self.backbone=[None]*numCA
            #self.ss=[None]*numCA
            i=0
            ssnum=0
            #prevss=None
            previntresi=None
            break_resi=False
            breakstop=[]
            for atom in model.atom:
                if atom.name==backboneatomname:
                    #if atom.__dict__.has_key("alt") and \
                    if use_firstaltloc==1 and \
                       atom.alt!='' and atom.alt!='A': #Accept the first alternate location for the atom and ignore others
                        continue
                    self.backbone[i]=atom.coord
                    #TODO Add option to change projection planes (x-y, y-z, z-x) by permutation, or use rotate axis, angle, selection.
                    #self.backbone[i]=(atom.coord[1],atom.coord[2],atom.coord[0])
##                    if atom.ss!=prevss: #No need to call cmd.dss to make sure. DNA (P backbone) has ss==''.
##                        ssnum+=1
##                    self.ss[i]=(atom.ss,ssnum) #Store secondary structure to help with knot projection simplification
##                    prevss=atom.ss
                    try:
                        intresi=int(atom.resi)
                    except ValueError:
                        intresi=int(atom.resi[:-1])
                    if previntresi is not None:
                        if intresi-previntresi>1:
                            print 'Warning: Skipped residue numbers (between %d and %d)! Artificial bond will be added to the backbone.' % (previntresi,intresi)
                            break_resi=True
                            breakstop.append(i)
                    previntresi=intresi
                    i+=1

            #Save the unreduced backbone, to be used for estimating the knotted core
            self.original_backbone=self.backbone[:]
            
            if break_resi:
                self.break_resi=True
                #Draw virtual bonds between residues that have non-consecutive residue numbers
                obj=[]
                for i in breakstop:
                    x1=self.backbone[i-1][0]
                    y1=self.backbone[i-1][1]
                    z1=self.backbone[i-1][2]
                    x2=self.backbone[i][0]
                    y2=self.backbone[i][1]
                    z2=self.backbone[i][2]
                    obj.extend([CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                               0.0, 1.0, 1.0, 0.0, 1.0, 1.0])
                structurename_chain=self.structurename+self.chainindicator
                cmd.delete("VIRTUALBONDS_%s" % structurename_chain)
                cmd.load_cgo(obj,"VIRTUALBONDS_%s" % structurename_chain)

            else:
                self.break_resi=False

            if use_reducebackbone==1:
                #TODO Test more, and change self.ss too
                print "Original backbone length: %d atoms" % len(self.backbone)
                try:
                    reducebackbone_option=self.show_reducebackbone_options.getvalue()
                    if reducebackbone_option==self.reducebackbone_options_tuple[0]:
                        self.reduceBackbone_sort(self.MINSHAPE) #Operates on the backbone that does NOT include the closure segments
                    elif reducebackbone_option==self.reducebackbone_options_tuple[1]:
                        self.reduceBackbone_random(self.MINSHAPE)
                except Exception, e:
                    print e
                print "Reduced backbone length: %d atoms" % len(self.backbone)
                #Draw reduced backbone
                i=0
                obj=[]
                while i+1<len(self.backbone):
                    x1=self.backbone[i][0]
                    y1=self.backbone[i][1]
                    z1=self.backbone[i][2]
                    x2=self.backbone[i+1][0]
                    y2=self.backbone[i+1][1]
                    z2=self.backbone[i+1][2]
                    #obj.extend([CYLINDER, x1, y1, z1, x2, y2, z2, 1,
                    #           0.0, 1.0, 0.0, 0.0, 1.0, 0.0])
                    obj.extend([SPHERE, x1, y1, z1, self.BACKBONE_RADIUS,
                                CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                               0.0, 1.0, 0.0, 0.0, 1.0, 0.0])
                    i+=1
                structurename_chain=self.structurename+self.chainindicator
                cmd.delete("BACKBONE_%s" % structurename_chain)
                cmd.load_cgo(obj,"BACKBONE_%s" % structurename_chain)
                if len(self.backbone)<5:
                    self.backbone=None
                    return False
            
            return True
        else:
            if numCA==0: #Warning: atom named P exists in protein 2efv
                print "No %s atoms found!" % backboneatomname
            self.backbone=None
            return False

    #Get the knot closure options from the UI
    def closeBackbone(self):
        closure_option=self.show_closure_options.getvalue()
        if closure_option==self.closure_options_tuple[0]:
            self.centerClosure()
        elif closure_option==self.closure_options_tuple[1]:
            self.directClosure()
        elif closure_option==self.closure_options_tuple[2]:
            self.customClosure()
            
    def directClosure(self):
        structurename_chain=self.structurename+self.chainindicator
        cmd.delete('DIRECTCLOSE_%s' % structurename_chain)
        #cmd.delete('CENTERCLOSE*%s' % structurename_chain)
        
        (x1,y1,z1)=self.backbone[0]
        (x2,y2,z2)=self.backbone[-1]
        radius=self.BACKBONE_RADIUS
        obj = [CYLINDER, x1, y1, z1, x2, y2, z2, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        cmd.load_cgo(obj, 'DIRECTCLOSE_%s' % structurename_chain)

    def centerClosure(self):
        structurename_chain=self.structurename+self.chainindicator
        cmd.delete('OUTWARDCLOSE_%s' % structurename_chain)
        #cmd.delete('DIRECTCLOSE_%s' % structurename_chain)

        #Get geometric center of backbone
        (xcm,ycm,zcm)=(0,0,0)
        for x,y,z in self.backbone:
            xcm+=x
            ycm+=y
            zcm+=z
        xcm*=1.0/len(self.backbone)
        ycm*=1.0/len(self.backbone)
        zcm*=1.0/len(self.backbone)

        #Get "radius" of protein by finding the backbone atom with largest distance form the CM
        maxr2=-1
        for x,y,z in self.backbone:
            r2=(x-xcm)*(x-xcm)+(y-ycm)*(y-ycm)+(z-zcm)*(z-zcm)
            if r2>maxr2:
                maxr2=r2
        maxr=1.5*math.sqrt(maxr2) #At least sqrt(2) times the most distant atom...

        radius=self.BACKBONE_RADIUS
        obj=[]

        #Create new nodes (atoms) to close the knot
        #(x3,y3,z3) connected to first backbone atom
        #(x4,y4,z4) connected to last backbone atom
        #(x5,y5,z5) connected to the previous two nodes
        (x1,y1,z1)=self.backbone[0]
        mag=math.sqrt((x1-xcm)*(x1-xcm)+(y1-ycm)*(y1-ycm)+(z1-zcm)*(z1-zcm))
        (x3,y3,z3)=(xcm+maxr*(x1-xcm)/mag,
                    ycm+maxr*(y1-ycm)/mag,
                    zcm+maxr*(z1-zcm)/mag)
        obj.extend([CYLINDER, x1, y1, z1, x3, y3, z3, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER1_%s' % structurename_chain)
        
        (x1,y1,z1)=self.backbone[-1]
        mag=math.sqrt((x1-xcm)*(x1-xcm)+(y1-ycm)*(y1-ycm)+(z1-zcm)*(z1-zcm))
        (x4,y4,z4)=(xcm+maxr*(x1-xcm)/mag,
                    ycm+maxr*(y1-ycm)/mag,
                    zcm+maxr*(z1-zcm)/mag)
        obj.extend([CYLINDER, x1, y1, z1, x4, y4, z4, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER2_%s' % structurename_chain)

        #TODO Handle exceptional cases
        #Find end-point (from CM) of angular bisector of previous two cylinders
        (x1,y1,z1)=(0.5*(x3+x4),0.5*(y3+y4),0.5*(z3+z4))
        mag=math.sqrt((x1-xcm)*(x1-xcm)+(y1-ycm)*(y1-ycm)+(z1-zcm)*(z1-zcm))
        (x5,y5,z5)=(xcm+maxr*(x1-xcm)/mag,
                    ycm+maxr*(y1-ycm)/mag,
                    zcm+maxr*(z1-zcm)/mag)

        obj.extend([CYLINDER, x5, y5, z5, x3, y3, z3, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER3_%s' % structurename_chain)

        obj.extend([CYLINDER, x5, y5, z5, x4, y4, z4, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER4_%s' % structurename_chain)

        cmd.load_cgo(obj, 'OUTWARDCLOSE_%s' % structurename_chain)
        
        #Add the coordinates of the new segments to the backbone
        self.backbone.append((x4,y4,z4))
        self.backbone.append((x5,y5,z5))
        self.backbone.append((x3,y3,z3))

        if self.printbackbone_var.get()==1:
            print "%7.3f %7.3f %7.3f (OUTWARDCLOSE_%s)" % (x4,y4,z4,structurename_chain)
            print "%7.3f %7.3f %7.3f (OUTWARDCLOSE_%s)" % (x5,y5,z5,structurename_chain)
            print "%7.3f %7.3f %7.3f (OUTWARDCLOSE_%s)" % (x3,y3,z3,structurename_chain)
        #Add bogus secondary structure entries for the new backbone nodes
##        self.ss.append(('X',-1))
##        self.ss.append(('X',-1))
##        self.ss.append(('X',-1))

    #Suggested by Bioinformatics Reviewer number 4.
    def customClosure(self):
        structurename_chain=self.structurename+self.chainindicator
        cmd.delete('CUSTOMCLOSE_%s' % structurename_chain)
        #cmd.delete('DIRECTCLOSE_%s' % structurename_chain)

        #Get geometric center of backbone
        (xcm,ycm,zcm)=(0,0,0)
        for x,y,z in self.backbone:
            xcm+=x
            ycm+=y
            zcm+=z
        xcm*=1.0/len(self.backbone)
        ycm*=1.0/len(self.backbone)
        zcm*=1.0/len(self.backbone)

        #Get "radius" of protein by finding the backbone atom with largest distance from the CM
        maxr2=-1
        for x,y,z in self.backbone:
            r2=(x-xcm)*(x-xcm)+(y-ycm)*(y-ycm)+(z-zcm)*(z-zcm)
            if r2>maxr2:
                maxr2=r2
        maxr=1.5*math.sqrt(maxr2) #At least sqrt(2) times the most distant atom...

        radius=self.BACKBONE_RADIUS
        obj=[]

        #Get direction of extension from N and C terminals provided by the user
        (nx,ny,nz)=(float(self.Nx.getvalue().strip()),
                    float(self.Ny.getvalue().strip()),
                    float(self.Nz.getvalue().strip()))
        (cx,cy,cz)=(float(self.Cx.getvalue().strip()),
                    float(self.Cy.getvalue().strip()),
                    float(self.Cz.getvalue().strip()))

        #Create new nodes (atoms) to close the knot
        #(x3,y3,z3) connected to first backbone atom
        #(x4,y4,z4) connected to last backbone atom
        #(x5,y5,z5) connected to the previous two nodes
        (x1,y1,z1)=self.backbone[0]

        #Let n be the unit vector extension from the terminal (N or C)
        #Let t be the magnitude of this extension that intersects the sphere of radius r=maxr about the CM
        #Let e be the vector from the CM to the terminal
        #Then to find t, solve t^2+2t*dotproduct(n,e)+dotproduct(e,e)=r^2
        mag=math.sqrt(nx*nx+ny*ny+nz*nz)
        nx/=mag
        ny/=mag
        nz/=mag
        (ex,ey,ez)=(x1-xcm,y1-ycm,z1-zcm)
        ndote=(nx*ex+ny*ey+nz*ez)
        e2=ex*ex+ey*ey+ez*ez
        r2=2.25*maxr2 #maxr^2
        t=math.sqrt(ndote*ndote+r2-e2)-ndote

        (x3,y3,z3)=(x1+t*nx,
                    y1+t*ny,
                    z1+t*nz)
        obj.extend([CYLINDER, x1, y1, z1, x3, y3, z3, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER1_%s' % structurename_chain)
        
        (x1,y1,z1)=self.backbone[-1]

        mag=math.sqrt(cx*cx+cy*cy+cz*cz)
        cx/=mag
        cy/=mag
        cz/=mag
        (ex,ey,ez)=(x1-xcm,y1-ycm,z1-zcm)
        ndote=(cx*ex+cy*ey+cz*ez)
        e2=ex*ex+ey*ey+ez*ez
        r2=2.25*maxr2 #maxr^2
        t=math.sqrt(ndote*ndote+r2-e2)-ndote

        (x4,y4,z4)=(x1+t*cx,
                    y1+t*cy,
                    z1+t*cz)
        obj.extend([CYLINDER, x1, y1, z1, x4, y4, z4, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER2_%s' % structurename_chain)

        #TODO Handle exceptional cases
        #Find end-point (from CM) of angular bisector of previous two cylinders
        (x1,y1,z1)=(0.5*(x3+x4),0.5*(y3+y4),0.5*(z3+z4))
        mag=math.sqrt((x1-xcm)*(x1-xcm)+(y1-ycm)*(y1-ycm)+(z1-zcm)*(z1-zcm))
        (x5,y5,z5)=(xcm+maxr*(x1-xcm)/mag,
                    ycm+maxr*(y1-ycm)/mag,
                    zcm+maxr*(z1-zcm)/mag)

        obj.extend([CYLINDER, x5, y5, z5, x3, y3, z3, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER3_%s' % structurename_chain)

        obj.extend([CYLINDER, x5, y5, z5, x4, y4, z4, radius,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #cmd.load_cgo(obj, 'CENTER4_%s' % structurename_chain)

        cmd.load_cgo(obj, 'CUSTOMCLOSE_%s' % structurename_chain)
        
        #Add the coordinates of the new segments to the backbone
        self.backbone.append((x4,y4,z4))
        self.backbone.append((x5,y5,z5))
        self.backbone.append((x3,y3,z3))

        if self.printbackbone_var.get()==1:
            print "%7.3f %7.3f %7.3f (CUSTOMCLOSE_%s)" % (x4,y4,z4,structurename_chain)
            print "%7.3f %7.3f %7.3f (CUSTOMCLOSE_%s)" % (x5,y5,z5,structurename_chain)
            print "%7.3f %7.3f %7.3f (CUSTOMCLOSE_%s)" % (x3,y3,z3,structurename_chain)
        
    class SegmentCrossings:
        def __init__(self,r1,r2):
            x1=r1[0]
            x2=r2[0]
            y1=r1[1]
            y2=r2[1]
            if x1!=x2:
                self.xy=True
                #Sort according to x
                #self.start=x1
                #self.end=x2
                if x1<x2:
                    self.dir=True
                else:
                    self.dir=False
            else:
                self.xy=False
                #Sort according to y
                #self.start=y1
                #self.end=y2
                if y1<y2:
                    self.dir=True
                else:
                    self.dir=False
            self.crossings=[]

        def addCrossing(self,x,y,ab,sign,ncrossinglabel):
            i=0
            while i<len(self.crossings):
                if self.xy:
                    xi=self.crossings[i][0]
                    if self.dir:
                        if x<xi:
                            break
                    else:
                        if xi<x:
                            break
                else:
                    yi=self.crossings[i][1]
                    if self.dir:
                        if y<yi:
                            break
                    else:
                        if yi<y:
                            break
                i+=1
##            if ncrossinglabel==12 or ncrossinglabel==8:
##                print "crossings.insert",i,(x,y,ab,sign,ncrossinglabel)
            #Insert, sort crossings accordingly (according to distance from the start of the segment)
            self.crossings.insert(i,(x,y,ab,sign,ncrossinglabel))
        
    #Brute force O(N^2)
    #Build gauss code and simplify knot projection using Reidemeister moves 1 and 2.
    def getCrossings(self):
        use_ss=self.ss_var.get()
        show_reidcrossings=self.reidemeister_var.get()
        show_crossings=self.crossings_var.get()
        #savecrossingsxyz=self.reducedknot_var.get()
        #cmd.delete('crossing*')

        crossingsxyz={}
        crossings={}
        N=len(self.backbone)
        ncrossinglabel=1
        for i in xrange(N):
##            ssi=self.ss[i][0]
##            ssnumi=self.ss[i][1]
            for j in xrange(N):
                #if use_ss==1 and (ssi=='H' or ssi=='S') and ssnumi==self.ss[j][1]:
##                if use_ss==1 and ssnumi==self.ss[j][1]:
##                    #Reidemeister moves already very effected at drastically reducing the number of crossings.
##                    continue #Residues on the same stretch of secondary structure do not intersect
                if i<j-1 and j-i!=N-1: #Check intersection of unconnected(not adjacent) segments
                    #Without j-i!=N-1, the knot invariants in different projections of some proteins
                    #e.g. 1xd3A and 1yveI, are inconsistent (probably not true anymore)
                    jnext=j+1
                    if jnext==N:
                        jnext=0
                    crossingdata=self.findIntersection(self.backbone[i],
                                                       self.backbone[i+1],
                                                       self.backbone[j],
                                                       self.backbone[jnext])
                    if crossingdata is not None:
                        (x,y,z1,z2,sign)=crossingdata
##                        if sign==+1:
##                            obj = [CYLINDER, x, y, z1, x, y, z2, 1,
##                                   1.0, 1.0, 0.0, 1.0, 1.0, 0.0]
##                        else:
##                            obj = [CYLINDER, x, y, z1, x, y, z2, 1,
##                                   0.0, 1.0, 0.0, 0.0, 1.0, 0.0]
                        
                        if z1>z2:
                            #obj.extend([SPHERE, x, y, z1, 1])
                            #cmd.load_cgo(obj,'crossing%d_%d' % (i,j))
                            ab=True #Is above (over crossing)
                            #TODO make this optional?
                            crossingsxyz[ncrossinglabel]=(x,y,z1,z2)
                        else:
                            #obj.extend([SPHERE, x, y, z2, 1])
                            #cmd.load_cgo(obj,'crossing%d_%d' % (j,i))
                            ab=False #Is below (under crossing)
                            #TODO make this optional?
                            crossingsxyz[ncrossinglabel]=(x,y,z2,z1)

                        try:
                            crossings[i].addCrossing(x,y,ab,sign,ncrossinglabel)
                        except KeyError:
                            crossings[i]=KAControlGroup.SegmentCrossings(self.backbone[i],self.backbone[i+1])
                            crossings[i].addCrossing(x,y,ab,sign,ncrossinglabel)

                        try:
                            crossings[j].addCrossing(x,y,not ab,sign,ncrossinglabel)
                        except KeyError:
                            crossings[j]=KAControlGroup.SegmentCrossings(self.backbone[j],self.backbone[jnext])
                            crossings[j].addCrossing(x,y,not ab,sign,ncrossinglabel)

                        ncrossinglabel+=1

        if show_reidcrossings==1:
            print "Building Gauss code while applying Reidemeister move 1:"
        #Automatically perform Reidemeister move 1
        ablist=[]
        signedlabellist=[]
        prevlabel=None
        for i in xrange(N):
            try:
                clist=crossings[i].crossings
                for x,y,ab,sign,label in clist:
                    if ab:
                        abtext='a'
                    else:
                        abtext='b'
                    if sign==+1:
                        signtext='+'
                    else:
                        signtext='-'

##                    if True:
##                        obj=[SPHERE, x, y, 0, 1]
##                        cmd.load_cgo(obj,'%s%s%d_%d' % (abtext, signtext, label, i))

                    signedlabel=sign*label
                    if signedlabel!=prevlabel:
                        ablist.append(abtext)
                        signedlabellist.append(signedlabel)
                        prevlabel=signedlabel
                        #DEBUG
                        #print abtext, signtext, label
                        #
                    else:
                        #Complete Reidemeister move 1 by removing the last crossing
                        ablist.pop()
                        signedlabellist.pop()
                        #DEBUG
                        #print 'Ignored via Reidemeister move 1'
                        if len(signedlabellist)>0:
                            prevlabel=signedlabellist[-1]
                        else:
                            prevlabel=None

                    if False:
                        axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]
                        obj=[]
                        pos = [x,y,0.0]
                        wire_text(obj,plain,pos,'%s%s%s' % (abtext,signtext,label),axes)
                        cmd.load_cgo(obj,'label%d' % i)
                    
            except KeyError:
                pass
        #TODO Check Reidemeister move 1 on first and last crossing
        if show_reidcrossings==1:
            print "crossing number=", len(ablist)/2

        if self.gausscode_var.get() and self.show_gausscode_options.getvalue()==self.gausscode_options_tuple[1]:
            #Build gauss code, with only most of the "nugatory" crossings removed
            self.gausscode=[None]*len(ablist)
            for i in xrange(len(ablist)):
                slabel=signedlabellist[i]
                if slabel>0:
                    sign='+'
                else:
                    slabel=-slabel
                    sign='-'
                self.gausscode[i]='%s%s%d' % (ablist[i],sign,slabel)

        if show_reidcrossings==1:
            print "Alternate application of a series of Reidemeister moves:"
        prev=len(ablist)+1
        while len(ablist)<prev:
            prev=len(ablist)
            #if show_reidcrossings==1:
            #    print "Reidemeister 2"
            self.Reidemeister2(ablist,signedlabellist)
            if show_reidcrossings==1:
                print "crossing number=", len(ablist)/2
                #print "Macro move with Reidemeister 1"
            self.macromove1(ablist,signedlabellist)
            if show_reidcrossings==1:
                print "crossing number=", len(ablist)/2
##            self.macromove2(ablist,signedlabellist)
##            if show_reidcrossings==1:
##                print "crossing number=", len(ablist)/2

            self.Reidemeister3_1(ablist,signedlabellist)
            if show_reidcrossings==1:
                print "crossing number=", len(ablist)/2
                #print [(ab,s) for ab, s in zip(ablist,signedlabellist)]

            self.Reidemeister3_2(ablist,signedlabellist)
            if show_reidcrossings==1:
                print "crossing number=", len(ablist)/2

            self.Reidemeister3_3(ablist,signedlabellist)
            if show_reidcrossings==1:
                print "crossing number=", len(ablist)/2

        if False: #if savecrossingsxyz==1:
            #Draw simplified knot based on reduced crossings
            structurename_chain=self.structurename+self.chainindicator
            #Get center of mass of crossings
            (xcm,ycm,zcm)=(0,0,0)
            for x,y,z1,z2 in crossingsxyz.values():
                xcm+=x+x
                ycm+=y+y
                zcm+=z1+z2
            xcm*=0.5/len(crossingsxyz)
            ycm*=0.5/len(crossingsxyz)
            zcm*=0.5/len(crossingsxyz)

            i=0
            obj=[]
            while i<len(ablist):
                if i==len(ablist)-1:
                    inext=0
                else:
                    inext=i+1
                label1=signedlabellist[i]
                if label1<0:
                    label1=-label1
                (x1,y1,z1a,z1b)=crossingsxyz[label1]
                if ablist[i]=='a':
                    z1=z1a
                else:
                    z1=z1b
                label2=signedlabellist[inext]
                if label2<0:
                    label2=-label2
                (x2,y2,z2a,z2b)=crossingsxyz[label2]
                if ablist[inext]=='a':
                    z2=z2a
                else:
                    z2=z2b
                #obj.extend([CYLINDER, x1, y1, z1, x2, y2, z2, 1,
                #       1.0, 0.0, 0.0, 1.0, 0.0, 0.0])

                #Test 2
##                dx=x2-x1
##                dy=y2-y1
##                dz=z2-z1
##                dx*=0.1
##                dy*=0.1
##                dz*=0.1
##                obj=[CYLINDER, x1-dx, y1-dy, z1-dz, x2+dx, y2+dy, z2+dz, 1,
##                       1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
                #Test 1
##                inc=2*i
                x=(x1+x2)-xcm
                y=(y1+y2)-ycm
                z=(z1+z2)-zcm
                obj=[CYLINDER, x1, y1, z1, x, y, z, 1,
                       1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                      CYLINDER, x2, y2, z2, x, y, z, 1,
                       1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
                i+=1
                cmd.load_cgo(obj,"KNOT_%s_%d" % (structurename_chain,i))

        if show_crossings==1:
            #TODO/DEBUG Simplify crossings by taking two crossings separated by an even number of crossings
            #then checking if there is a segment of the chain in between
##            success=False
##            for i1 in xrange(len(ablist)):
##                for i2 in xrange(len(ablist)):
##                    if i1+2<i2 and (i2-i1-1)%2==0:
##                        label1=signedlabellist[i1]
##                        if label1<0:
##                            label1=-label1
##                        label2=signedlabellist[i2]
##                        if label2<0:
##                            label2=-label2
##                        x1,y1,z11,z12=crossingsxyz[label1]
##                        x2,y2,z21,z22=crossingsxyz[label2]
##                        r1=(x1,y1)
##                        r2=(x2,y2)
##                        empty=True
##                        for j in xrange(N):
##                            jnext=j+1
##                            if jnext==N:
##                                jnext=0
##                            if self.isIntersecting(self.backbone[j],self.backbone[jnext],r1,r2):
##                                empty=False
##                                break
##                        if empty: #Remove the crossings between label1 and label2
##                            abpure=True
##                            for j in range(i1+1,i2):
##                                if ablist[i1+1]!=ablist[j]: #Crossings should be all 'a' or all 'b'
##                                    abpure=False
##                                    break
##                            if abpure:
##                                success=Tue
##                                ablist[i1+1:i2]=[]
##                                tmplist=signedlabellist[i1+1:i2] #Make copy first
##                                signedlabellist[i1+1:i2]=[] #Remove the crossings
##                                for s in tmplist: #Find and remove the partner crossings
##                                    index=signedlabellist.index(s)
##                                    signedlabellist.pop(index)
##                                    ablist.pop(index)
##                                break
##                if success:
##                    print "Got one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
##                    break

            
            structurename_chain=self.structurename+self.chainindicator
            obj=[]
            i=0
            while i<len(ablist):
                if ablist[i]=="a":
                    label=signedlabellist[i]
                    if label<0:
                        label=-label
                    x,y,z1,z2=crossingsxyz[label]
                    obj.extend([CYLINDER, x, y, z1, x, y, z2, self.CROSSING_RADIUS,
                               1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
                i+=1
            cmd.delete("CROSSINGS_%s" % structurename_chain)
            cmd.load_cgo(obj,"CROSSINGS_%s" % structurename_chain)


        if self.gausscode_var.get() and self.show_gausscode_options.getvalue()==self.gausscode_options_tuple[0]:
            #Build gauss code, after Reidemeister moves have been performed
            self.gausscode=[None]*len(ablist)
            for i in xrange(len(ablist)):
                slabel=signedlabellist[i]
                if slabel>0:
                    sign='+'
                else:
                    slabel=-slabel
                    sign='-'
                self.gausscode[i]='%s%s%d' % (ablist[i],sign,slabel)

        return ablist, signedlabellist

############################### Reidemeister moves
    
    #Under construction
    def macromove2(self,ablist,signedlabellist):
        i=0
        while i+1<len(ablist):
            label1=signedlabellist[i]
            label2=signedlabellist[i+1]
            if label1!=label2:
                j=0
                matchandremove=False
                nextpair=False
                while j+1<i:
                    if label1==signedlabellist[j] or label2==signedlabellist[j]:
                        j+=1
                        firstab=ablist[j]
                        j1=j
                        j+=1
                        while j<i:
                            if firstab!=ablist[j]:
                                if label1==signedlabellist[j] or label2==signedlabellist[j]:
                                    ablist[j1:j]=[]
                                    tmplist=signedlabellist[j1:j]
                                    signedlabellist[j1:j]=[]
                                    for s in tmplist:
                                        index=signedlabellist.index(s) #TODO Use .remove?
                                        signedlabellist.pop(index)
                                        ablist.pop(index)
                                    matchandremove=True
                                break
                            elif label1==signedlabellist[j] or label2==signedlabellist[j]:
                                ablist[j1:j]=[]
                                tmplist=signedlabellist[j1:j]
                                signedlabellist[j1:j]=[]
                                for s in tmplist:
                                    index=signedlabellist.index(s) #TODO Use .remove?
                                    signedlabellist.pop(index)
                                    ablist.pop(index)
                                matchandremove=True
                            j+=1
                        break
                    j+=1
                if matchandremove:
                   break
                #TODO Find loop between i and i+1
            i+=1

    #Slide a loop that overlaps a portion of the knot
    #
    #        __
    #       /  \
    #    XXXXXXXXXX           XXXXXXXXXX
    #   XXXXXXXXXXXX ---->   XXXXXXXXXXXX
    #    XXXXXXXXXX           XXXXXXXXXX
    #      \    /
    #       \  /                  /\
    #        \/                   \/
    #        /\                   /\
    #
    def macromove1(self,ablist,signedlabellist):
        i=0
        while i+1<len(ablist):
            matchandremove=False
            startcrossing=signedlabellist[i]
            if True:
                j=i+1
                firstab=ablist[j]
                while j<len(ablist):
                    if firstab!=ablist[j]:
                        if startcrossing==signedlabellist[j]:
                            #DEBUG
                            #print "Removing", startcrossing, signedlabellist[j]
                            ablist[i:j]=[]
                            tmplist=signedlabellist[i:j] #Make copy first
                            signedlabellist[i:j]=[] #Remove the crossings
                            for s in tmplist: #Find and remove the partner crossings
                                index=signedlabellist.index(s)
                                signedlabellist.pop(index)
                                ablist.pop(index)
                            matchandremove=True
                            #if i>0:
                            #    i=i-1
                            i=0
                        break
                    elif startcrossing==signedlabellist[j]:
                        #DEBUG
                        #print "Removing", startcrossing, signedlabellist[j]
                        ablist[i:j]=[]
                        tmplist=signedlabellist[i:j] #Make copy first
                        signedlabellist[i:j]=[] #Remove the crossings
                        for s in tmplist: #Find and remove the partner crossings
                            index=signedlabellist.index(s)
                            signedlabellist.pop(index)
                            ablist.pop(index)
                        matchandremove=True
                        #if i>0:
                        #    i=i-1
                        i=0
                        break
                    j+=1
            if matchandremove:
                pass
            else:
                i+=1

        #If applicable, apply Reidemeister move 1 on first and last crossing
        if len(signedlabellist)>1 and signedlabellist[0]==signedlabellist[-1]:
            signedlabellist.pop()
            signedlabellist.pop(0)
            ablist.pop()
            ablist.pop(0)

    #Not used anymore because macromove1 is more general
    def Reidemeister1(self,ablist,signedlabellist):
        i=0
        while i+1<len(ablist):
            matchandremove=False
            if signedlabellist[i]==signedlabellist[i+1]:
                #DEBUG
                #print "Removing", signedlabellist[i], signedlabellist[i+1]
                ablist[i:i+2]=[] #or use pop twice
                signedlabellist[i:i+2]=[]
                matchandremove=True
            if matchandremove:
                if i>0:
                    i=i-1
            else:
                i+=1

    #    /\
    # ---------- ---->  ---------- Net reduction of two crossings
    #  /    \             _____
    # /      \           /     \
    def Reidemeister2(self,ablist,signedlabellist):
        i=0
        while i+1<len(ablist):
            #Find adjacent crossings having both 'a' or both 'b'
            matchandremove=False
            if ablist[i]==ablist[i+1]:
                j=i+2
                while j+1<len(ablist):
                    if (signedlabellist[i]==signedlabellist[j] and \
                       signedlabellist[i+1]==signedlabellist[j+1]) or \
                        (signedlabellist[i]==signedlabellist[j+1] and \
                       signedlabellist[i+1]==signedlabellist[j]):
                        #DEBUG
                        #print "Removing", signedlabellist[i], signedlabellist[i+1]
##                        ablist.pop(j+1)
##                        ablist.pop(j)
##                        ablist.pop(i+1)
##                        ablist.pop(i)
##                        signedlabellist.pop(j+1)
##                        signedlabellist.pop(j)
##                        signedlabellist.pop(i+1)
##                        signedlabellist.pop(i)
                        ablist[j:j+2]=[]
                        ablist[i:i+2]=[]
                        signedlabellist[j:j+2]=[]
                        signedlabellist[i:i+2]=[]
                        matchandremove=True
                        break
                    j+=1
            if matchandremove:
                pass
            else:
                i+=1

    #Find four consecutive undercrossings (or overcrossings) "surrounding" a fifth crossing
    #
    #     __|_                      |
    #    \  | \                     |
    #  -----|-----     ---->   -----|----- Net reduction of four crossings
    #      \|   \                   |
    #       |    \                  |
    #       |\    \                 | _____
    #       | \    \                | \    \
    #
    def Reidemeister3_1(self,ablist,signedlabellist):
        i=0
        while i+5<len(ablist):
            matchandremove=False
            if ablist[i]==ablist[i+1] and \
               ablist[i]==ablist[i+2] and \
               ablist[i]==ablist[i+3]:
                j=0
                centrallabel=None
                while j<len(ablist):
                    if j<i or i+3<j:
                        jnext=j+1
                        if jnext==len(ablist):
                            jnext=0
                        jprev=j-1
                        if jprev<0:
                            jprev=len(ablist)-1
                        if (signedlabellist[i]==signedlabellist[jnext] and \
                           signedlabellist[i+2]==signedlabellist[jprev]) or \
                           (signedlabellist[i+2]==signedlabellist[jnext] and \
                           signedlabellist[i]==signedlabellist[jprev]):
                            if centrallabel is None:
                                centrallabel=signedlabellist[j]
                            elif centrallabel==signedlabellist[j]:
                                matchandremove=True
                                break
                            else:
                                break
                        if (signedlabellist[i+1]==signedlabellist[jnext] and \
                           signedlabellist[i+3]==signedlabellist[jprev]) or \
                           (signedlabellist[i+3]==signedlabellist[jnext] and \
                           signedlabellist[i+1]==signedlabellist[jprev]):
                            if centrallabel is None:
                                centrallabel=signedlabellist[j]
                            elif centrallabel==signedlabellist[j]:
                                matchandremove=True
                                break
                            else:
                                break
                    j+=1
            if matchandremove:
                ablist[i:i+4]=[]
                tmplist=signedlabellist[i:i+4]
                signedlabellist[i:i+4]=[]
                for s in tmplist: #Find and remove the partner crossings
                    index=signedlabellist.index(s)
                    signedlabellist.pop(index)
                    ablist.pop(index)
                #if i>0:
                #    i=i-1
            else:
                i+=1

    #Find three consecutive undercrossings (or overcrossings) "surrounding" a fourth crossing
    #WARNING: This method swaps the order of two crossings, so when the crossings are marked with
    #cylinders in PyMOL, the knot diagram may not make sense.\
    #
    #      _|_                      |
    #     / | \                     |
    #  ------------    ---->   ------------ Net reduction of two crossings
    #   /   |   \               ____|____
    #  /    |    \             /    |    \
    #
    def Reidemeister3_2(self,ablist,signedlabellist):
        i=0
        while i+4<len(ablist):
            matchandremove=False
            if ablist[i]==ablist[i+1] and \
               ablist[i]==ablist[i+2]:
                j=0
                centrallabel=None
                while j<len(ablist):
                    if j<i or i+2<j:
                        jnext=j+1
                        if jnext==len(ablist):
                            jnext=0
                        jprev=j-1
                        if jprev<0:
                            jprev=len(ablist)-1
                        if (signedlabellist[i]==signedlabellist[jnext] and \
                           signedlabellist[i+2]==signedlabellist[jprev]) or \
                           (signedlabellist[i+2]==signedlabellist[jnext] and \
                           signedlabellist[i]==signedlabellist[jprev]):
                            centrallabel=signedlabellist[j]
                            break
                    j+=1
                if centrallabel is not None:
                    j=0
                    while j<len(ablist):
                        if j<i or i+2<j:
                            jnext=j+1
                            if jnext==len(ablist):
                                jnext=0
                            jprev=j-1
                            if jprev<0:
                                jprev=len(ablist)-1
                            if (signedlabellist[i+1]==signedlabellist[jnext] and \
                               centrallabel==signedlabellist[j]):
                                matchandremove=True
                                j1=j
                                j2=jnext
                                break
                            if (signedlabellist[i+1]==signedlabellist[jprev] and \
                               centrallabel==signedlabellist[j]):
                                matchandremove=True
                                j1=j
                                j2=jprev
                                break
                        j+=1
            if matchandremove:
                #Swap order of crossing j1 and crossing j2
                #Delete crossing i+2 and i (and their partners)
                tmp=ablist[j1]
                ablist[j1]=ablist[j2]
                ablist[j2]=tmp
                tmp=signedlabellist[j1]
                signedlabellist[j1]=signedlabellist[j2]
                signedlabellist[j2]=tmp
                ablist.pop(i+2)
                ablist.pop(i)
                tmplist=(signedlabellist[i],signedlabellist[i+2])
                signedlabellist.pop(i+2)
                signedlabellist.pop(i)
                for s in tmplist: #Find and remove the partner crossings
                    index=signedlabellist.index(s)
                    signedlabellist.pop(index)
                    ablist.pop(index)
                #if i>0:
                #    i=i-1
            else:
                i+=1

    #WARNING: This method swaps the order of two crossings, so when the crossings are marked with
    #cylinders in PyMOL, the knot diagram may not make sense.
    #Combine Reidemeister move 1 and 3:
    #      __                             _
    #     /  \                           / \
    # ---/--------    ---->    -------------\--  Net reduction of one crossing
    #    \    /                        /    /
    #     \  /                        /    /
    #      \                         /    /
    #      /\                       /    /
    def Reidemeister3_3(self,ablist,signedlabellist):
        i=0
        while i+3<len(ablist):
            matchandremove=False
            if signedlabellist[i]==signedlabellist[i+3] and \
               ablist[i]==ablist[i+1] and \
               ablist[i+2]==ablist[i+3] and \
               signedlabellist[i+1]!=signedlabellist[i+2]:
                j=0
                while j<len(ablist):
                    if j<i or i+3<j:
                        jnext=j+1
                        if jnext==len(ablist):
                            jnext=0
                        if (signedlabellist[i+1]==signedlabellist[j] and \
                           signedlabellist[i+2]==signedlabellist[jnext]) or \
                           (signedlabellist[i+2]==signedlabellist[j] and \
                           signedlabellist[i+1]==signedlabellist[jnext]):
                            j1=j
                            j2=jnext
                            matchandremove=True
                            break
                    j+=1
            if matchandremove:
                #Swap order of crossing j and crossing jnext
                #Delete crossing i and i+3
                tmp=ablist[j1]
                ablist[j1]=ablist[j2]
                ablist[j2]=tmp
                tmp=signedlabellist[j1]
                signedlabellist[j1]=signedlabellist[j2]
                signedlabellist[j2]=tmp
                ablist.pop(i+3)
                ablist.pop(i)
                signedlabellist.pop(i+3)
                signedlabellist.pop(i)
                #if i>0:
                #    i=i-1
            else:
                i+=1

################################ Computational geometry
                
    #Iteration one: brute force
    #1. Find intersection of x-y projection between segment r1-r2 and segment r3-r4. r1 is a vector (x,y,z).
    #2. If they intersect, find the z-coordinate of each segment (used to find which segment goes over(under) the other).
    #3. Return sign of the crossing.
    def findIntersection(self,r1,r2,r3,r4):
        #Solve for x and y in
        #y=m1*(x-x1)+y1, m1=(y2-y1)/(x2-x1)
        #y=m2*(x-x3)+y3, m2=(y4-y3)/(x4-x3)
        dx1=r2[0]-r1[0]
        dy1=r2[1]-r1[1]
        dx2=r4[0]-r3[0]
        dy2=r4[1]-r3[1]

        #Handle exceptional cases
        #TODO Check that the case when the segment endpoints may coincide is handled properly
        if dx1==0 and dx2==0:
            if r1[0]!=r3[0]:
                return None
            print "Warning: Encountered two overlapping vertical segments! Perturbing one segment..."
            r1[0]+=0.01
            #r3[0]+=0.01
            dx1=r2[0]-r1[0]
            dx2=r4[0]-r3[0] #TODO r3 not modified anymore, so why bother?
        if dy1==0 and dy2==0:
            if r1[1]!=r3[1]:
                return None
            print "Warning: Encountered two overlapping horizontal segments! Perturbing one segment..."
            r1[1]+=0.01
            #r3[1]+=0.01
            dy1=r2[1]-r1[1]
            dy2=r4[1]-r3[1] #TODO r3 not modified anymore, so why bother?
        if dx1==0 and dy2==0:
            x=r1[0]
            y=r3[1]
        elif dy1==0 and dx2==0:
            x=r3[0]
            y=r1[1]
        elif dx1==0:
            x=r1[0]
            y=(x-r3[0])*dy2/dx2+r3[1]
        elif dx2==0:
            x=r3[0]
            y=(x-r1[0])*dy1/dx1+r1[1]
        elif dy1==0:
            y=r1[1]
            x=(y-r3[1])*dx2/dy2+r3[0]
        elif dy2==0:
            y=r3[1]
            x=(y-r1[1])*dx1/dy1+r1[0]
        else:
            m1=dy1/dx1
            m2=dy2/dx2
            #What if the slopes are the same? (parallel segments)
            if m1==m2:
                print "Warning: Encountered two parallel segments!  Perturbing one segment..."
                r1[0]+=0.01
                dx1=r2[0]-r1[0]
                m1=dy1/dx1
            x=(m1*r1[0]-m2*r3[0]+r3[1]-r1[1])/(m1-m2) #TODO reverse roles of x and y for better numerical precision?
            y=m1*(x-r1[0])+r1[1]
            #y=m2*(x-r3[0])+r3[1] #TODO take the average of this with the other one
        #Check that (x,y) lies on the segments for proper intersection
        if self.isBetween(x,r1[0],r2[0]) and self.isBetween(y,r1[1],r2[1]):
            if self.isBetween(x,r3[0],r4[0]) and self.isBetween(y,r3[1],r4[1]):
                if dx1!=0:
                    z1=(r2[2]-r1[2])*(x-r1[0])/dx1+r1[2]
                else:
                    z1=(r2[2]-r1[2])*(y-r1[1])/dy1+r1[2]
                if dx2!=0:
                    z2=(r4[2]-r3[2])*(x-r3[0])/dx2+r3[2]
                else:
                    z2=(r4[2]-r3[2])*(y-r3[1])/dy2+r3[2]
                #Use cross product to get the sign of the crossing
                if dx1*dy2-dx2*dy1>0:
                    if z1>z2:
                        sign=-1
                    else:
                        sign=+1
                else: #TODO check if cross product is zero
                    if z1>z2:
                        sign=+1
                    else:
                        sign=-1
                return (x,y,z1,z2,sign)
        return None

    #2D segments
    #TODO Handle case when m1==m2, if this method gets put back to use.
    def isIntersecting(self,r1,r2,r3,r4):
        #Solve for x and y in
        #y=m1*(x-x1)+y1, m1=(y2-y1)/(x2-x1)
        #y=m2*(x-x3)+y3, m2=(y4-y3)/(x4-x3)
        dx1=r2[0]-r1[0]
        dy1=r2[1]-r1[1]
        dx2=r4[0]-r3[0]
        dy2=r4[1]-r3[1]

        #Handle exceptional cases
        #TODO Check that the case when the segment endpoints may coincide is handled properly
        if dx1==0 and dx2==0:
            if r1[0]!=r3[0]:
                return None
            print "Warning: Encountered two overlapping vertical segments! Perturbing both segments..."
            r1[0]+=0.01
            #r3[0]+=0.01
            dx1=r2[0]-r1[0]
            dx2=r4[0]-r3[0]
        if dy1==0 and dy2==0:
            if r1[1]!=r3[1]:
                return None
            print "Warning: Encountered two overlapping horizontal segments! Perturbing one segment..."
            r1[1]+=0.01
            #r3[1]+=0.01
            dy1=r2[1]-r1[1]
            dy2=r4[1]-r3[1]
        if dx1==0 and dy2==0:
            x=r1[0]
            y=r3[1]
        elif dy1==0 and dx2==0:
            x=r3[0]
            y=r1[1]
        elif dx1==0:
            x=r1[0]
            y=(x-r3[0])*dy2/dx2+r3[1]
        elif dx2==0:
            x=r3[0]
            y=(x-r1[0])*dy1/dx1+r1[1]
        elif dy1==0:
            y=r1[1]
            x=(y-r3[1])*dx2/dy2+r3[0]
        elif dy2==0:
            y=r3[1]
            x=(y-r1[1])*dx1/dy1+r1[0]
        else:
            m1=dy1/dx1
            m2=dy2/dx2
            x=(m1*r1[0]-m2*r3[0]+r3[1]-r1[1])/(m1-m2) #TODO reverse roles of x and y for better numerical precision?
            y=m1*(x-r1[0])+r1[1]
            #y=m2*(x-r3[0])+r3[1] #TODO take the average of this with the other one
        #Check that (x,y) lies on the segments for proper intersection
        if self.isBetween(x,r1[0],r2[0]) and self.isBetween(y,r1[1],r2[1]):
            if self.isBetween(x,r3[0],r4[0]) and self.isBetween(y,r3[1],r4[1]):
                return True
        else:
            return False

    def isBetween(self,x,x1,x2):
        if x1<x2:
            if x1<=x and x<=x2:
                return True
        else: #if x2<x1:
            if x2<=x and x<=x1:
                return True
        return False

    #en.wikipedia.org/wiki/Line-plane_intersection
    #t1-t2-t3 form vertices of the triangle
    #r1-r2 form a line segment
    def findIntersectionTriangle(self,t1,t2,t3,r1,r2):
        x0=r1[0]-t1[0]
        y0=r1[1]-t1[1]
        z0=r1[2]-t1[2]

        x1=r1[0]-r2[0]
        y1=r1[1]-r2[1]
        z1=r1[2]-r2[2]

        x2=t2[0]-t1[0]
        y2=t2[1]-t1[1]
        z2=t2[2]-t1[2]

        x3=t3[0]-t1[0]
        y3=t3[1]-t1[1]
        z3=t3[2]-t1[2]

        det=x1*y2*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1-x1*y3*z2-x2*y1*z3
        if det==0.0:
            print "Warning: Found a backbone segment parallel to a triplet's plane!"
            return False

        #Solve for intersection in terms of parametric equations for a line and a plane
        #s is the parameter along the line/segment
        #u and v are parameters of the plane
        s=(x0*y2*z3+x2*y3*z0+x3*y0*z2-x3*y2*z0-x0*y3*z2-x2*y0*z3)/det
        u=(x1*y0*z3+x0*y3*z1+x3*y1*z0-x3*y0*z1-x1*y3*z0-x0*y1*z3)/det
        v=(x1*y2*z0+x2*y0*z1+x0*y1*z2-x0*y2*z1-x1*y0*z2-x2*y1*z0)/det

        #if t1==r1 or t1==r2 or t3==r1 or t3==r2:

        #if 0<=s and s<=1 and 0<=u and u<=1 and 0<=v and v<=1 and (u+v)<=1:
        #   return True
        #margins (upper and lower bounds)
        #U=1.2
        #L=-0.2
        U=1.0
        L=0.0
        if L<s and s<U and L<u and u<U and L<v and v<U and (u+v)<U:
            return True
        return False

    #I don't like the reduced backbone when I use getArea2 to sort below (maybe adjust parameters?)
    def getArea2(self,t1,t2,t3):
        x2=t2[0]-t1[0]
        y2=t2[1]-t1[1]
        z2=t2[2]-t1[2]

        x3=t3[0]-t1[0]
        y3=t3[1]-t1[1]
        z3=t3[2]-t1[2]

        x=y2*z3-z2*y3
        y=z2*x3-x2*z3
        z=x2*y3-y2*x3

        #returns the square of the area of the parallelogram
        return x*x+y*y+z*z

    #TODO Define shape function as in FE?
    def getShape(self,tm,t,tp):
        x1=tp[0]-tm[0]
        y1=tp[1]-tm[1]
        z1=tp[2]-tm[2]

        x2=tm[0]-t[0]
        y2=tm[1]-t[1]
        z2=tm[2]-t[2]

        x3=tp[0]-t[0]
        y3=tp[1]-t[1]
        z3=tp[2]-t[2]

        #The value is 1 if t is a right angle (sum of squares of two sides, divided by the square of the edge tm-tp)
        return ((x2*x2+y2*y2+z2*z2)+(x3*x3+y3*y3+z3*z3))/(x1*x1+y1*y1+z1*z1)


    #TODO Include the closure segments
    #Leave the terminals fixed
    def reduceBackbone(self,minshape):
        N=len(self.backbone)
        #areas=[0]*N
        shapes=[0]*N
        i=1
        while i+1<N:
            iprev=i-1
            inext=i+1
            #areas[i]=self.getArea2(self.backbone[iprev],self.backbone[i],self.backbone[inext])
            shapes[i]=self.getShape(self.backbone[iprev],self.backbone[i],self.backbone[inext])
            i+=1

        prevlength=None
        while prevlength!=len(self.backbone):
            prevlength=len(self.backbone)
            i=1
            #matchandremove=False
            while i+1<len(self.backbone):
                iprev=i-1
                inext=i+1
                if shapes[i]>minshape:
                    #Check if it is safe to remove the triangle
                    j=0
                    while j+1<len(self.backbone):
                        jnext=j+1
                        if jnext<iprev or inext<j: #Don't check intersection with segments adjacent to the triangle
                            if self.findIntersectionTriangle(self.backbone[iprev],
                                                             self.backbone[i],
                                                             self.backbone[inext],
                                                             self.backbone[j],
                                                             self.backbone[jnext]):
                                #foundintersection=True
                                break
                        j+=1
                    if j+1==len(self.backbone):
                        #Remove triangle
                        self.backbone.pop(i)
                        shapes.pop(i)
                        #Compute shape of affected triangles
                        if iprev>0:
                            shapes[iprev]=self.getShape(self.backbone[iprev-1],self.backbone[iprev],self.backbone[iprev+1])
                        if inext<len(self.backbone):
                            shapes[i]=self.getShape(self.backbone[iprev],self.backbone[i],self.backbone[inext])
                        break
                i+=1

    #Sort triangles according to area or shape
    #TODO Include the closure segments
    #Leave the terminals fixed
    def reduceBackbone_sort(self,minshape):
        N=len(self.backbone)
        #areas=[0]*N
        shapes=[0]*N #Warning: terminals included in the list. The triangles at the terminals are not inspected.
        i=1
        while i+1<N:
            iprev=i-1
            inext=i+1
            #areas[i]=self.getArea2(self.backbone[iprev],self.backbone[i],self.backbone[inext])
            shapes[i]=self.getShape(self.backbone[iprev],self.backbone[i],self.backbone[inext])
            i+=1
        tmplist=zip(range(0,N),shapes)
        i_shapesorted=sorted(tmplist,key=lambda x: x[1], reverse=True)
        #i_shapesorted=sorted(tmplist,key=lambda x: x[1], reverse=False) #Ascending
        #print i_shapesorted
        
        prevlength=None
        while prevlength!=len(self.backbone):
            prevlength=len(self.backbone)
            #i=1
            #matchandremove=False
            #while i+1<len(self.backbone):
            for i,s in i_shapesorted:
                iprev=i-1
                inext=i+1
                if shapes[i]>minshape:
                    #Check if it is safe to remove the triangle
                    j=0
                    while j+1<len(self.backbone):
                        jnext=j+1
                        if jnext<iprev or inext<j: #Don't check intersection with segments adjacent to the triangle
                            if self.findIntersectionTriangle(self.backbone[iprev],
                                                             self.backbone[i],
                                                             self.backbone[inext],
                                                             self.backbone[j],
                                                             self.backbone[jnext]):
                                #foundintersection=True
                                break
                        j+=1
                    if j+1==len(self.backbone):
                        #Remove triangle
                        self.backbone.pop(i)
                        shapes.pop(i)

                        tmplist=zip(range(0,len(shapes)),shapes)
                        i_shapesorted=sorted(tmplist,key=lambda x: x[1], reverse=True) #TODO Make sorting calls more efficient?
                        #i_shapesorted=sorted(tmplist,key=lambda x: x[1], reverse=False) #TODO Make sorting calls more efficient?

                        #Compute shape of affected triangles
                        if iprev>0:
                            shapes[iprev]=self.getShape(self.backbone[iprev-1],self.backbone[iprev],self.backbone[iprev+1])
                        if inext<len(self.backbone):
                            shapes[i]=self.getShape(self.backbone[iprev],self.backbone[i],self.backbone[inext])
                        break
                #i+=1

    #Random shuffling was able to help reduce the number of crossings in the 3bjxA (Stevedore's) knot to six (6)!
    #TODO Include the closure segments
    #Leave the terminals fixed
    def reduceBackbone_random(self,minshape):
        N=len(self.backbone)
        #areas=[0]*N
        shapes=[0]*N #Warning: terminals included in the list. The triangles at the terminals are not inspected.
        i=1
        while i+1<N:
            iprev=i-1
            inext=i+1
            #areas[i]=self.getArea2(self.backbone[iprev],self.backbone[i],self.backbone[inext])
            shapes[i]=self.getShape(self.backbone[iprev],self.backbone[i],self.backbone[inext])
            i+=1
        ishapes_list=zip(range(0,N),shapes)
        #Try random order
        random.shuffle(ishapes_list)
        
        prevlength=None
        while prevlength!=len(self.backbone):
            prevlength=len(self.backbone)
            i=1
            #matchandremove=False
            #while i+1<len(self.backbone):
            for i,s in ishapes_list:
                iprev=i-1
                inext=i+1
                if shapes[i]>minshape:
                    #Check if it is safe to remove the triangle
                    j=0
                    while j+1<len(self.backbone):
                        jnext=j+1
                        if jnext<iprev or inext<j: #Don't check intersection with segments adjacent to the triangle
                            if self.findIntersectionTriangle(self.backbone[iprev],
                                                             self.backbone[i],
                                                             self.backbone[inext],
                                                             self.backbone[j],
                                                             self.backbone[jnext]):
                                #foundintersection=True
                                break
                        j+=1
                    if j+1==len(self.backbone):
                        #Remove triangle
                        self.backbone.pop(i)
                        shapes.pop(i)

                        ishapes_list=zip(range(0,len(shapes)),shapes)
                        #Try random order
                        random.shuffle(ishapes_list)

                        #Compute shape of affected triangles
                        if iprev>0:
                            shapes[iprev]=self.getShape(self.backbone[iprev-1],self.backbone[iprev],self.backbone[iprev+1])
                        if inext<len(self.backbone):
                            shapes[i]=self.getShape(self.backbone[iprev],self.backbone[i],self.backbone[inext])
                        break
                #i+=1

            
################################# Compute knot invariants

    #Computation based on a paper by Vologodskii, et al
    def computeAlexander(self,ablist,num):
	numunderpasses=len(ablist)/2; #numoverpasses too
	numcrossings=len(ablist)
	underpassnumlist=[0]*numcrossings
	gennumlist=[0]*numcrossings
	#Start at the first underpass ('b')
        istart=0
        while istart<numcrossings:
            if ablist[istart]=='b':
                break
            istart+=1
        if istart==numcrossings:
            print "No underpass found. Protein is unknotted or there is an error!"
            return None
        #Assign underpass and generator numbers
	gennum=1 #This is important
	underpassnum=1
	for i in xrange(numcrossings):
            j=(istart+i)%numcrossings
            if ablist[j]=='b':
                underpassnumlist[j]=underpassnum
                underpassnum+=1
                if gennum>=numunderpasses:
                    #print 'gennum>=numunderpasses'
                    gennum=1
                else:
                    gennum+=1
            else:
                gennumlist[j]=gennum
        #print "generator numbers", gennumlist
        #print "underpass numbers", underpassnumlist

	for i in xrange(numcrossings):
            if ablist[i]=='b':
                for j in xrange(numcrossings):
                    if i!=j and num[i]==num[j]:
                        gennumlist[i]=gennumlist[j] #Assign generator number of the overpass to the underpass
                        #print "underpass and generator numbers", underpassnumlist[i], gennumlist[i]
                        break
        
        #Construct Alexander matrix
        am=zeros((numunderpasses,numunderpasses))
        tvar=-1.0
	for ii in xrange(numcrossings):
            jj=(istart+ii)%numcrossings
            if ablist[jj]=='b':
                #Assign elements based on Type I or Type II underpass
                k=underpassnumlist[jj]
                i=gennumlist[jj] #overpassing generator number
                if i==k or i==k+1: #Rule 1
                    am[k-1][k-1]=-1.0
                    if k<numunderpasses:
                        am[k-1][k]=1.0
##                    m=1
##                    while m<=numunderpasses:
##                        if m==k:
##                            am[k-1][m-1]=-1.0
##                        elif m==k+1:
##                            am[k-1][m-1]=1.0
##                        else:
##                            am[k-1][m-1]=0
##                        m+=1
                else: #Rule 2
                    if num[jj]>0: #num[jj]<0 works also. Choice of which sign is Type I or Type II is a convention
                        am[k-1][k-1]=-tvar
                        if k<numunderpasses:
                            am[k-1][k]=1.0
                        am[k-1][i-1]=tvar-1.0
##                        m=1
##                        while m<=numunderpasses:
##                            if m==k:
##                                am[k-1][m-1]=-tvar
##                            elif m==k+1:
##                                am[k-1][m-1]=1.0
##                            elif m==i:
##                                am[k-1][m-1]=tvar-1.0
##                            else:
##                                am[k-1][m-1]=0
##                            m+=1
                    else:
                        am[k-1][k-1]=1.0
                        if k<numunderpasses:
                            am[k-1][k]=-tvar
                        am[k-1][i-1]=tvar-1.0
##                        m=1
##                        while m<=numunderpasses:
##                            if m==k:
##                                am[k-1][m-1]=1.0
##                            elif m==k+1:
##                                am[k-1][m-1]=-tvar
##                            elif m==i:
##                                am[k-1][m-1]=tvar-1.0
##                            else:
##                                am[k-1][m-1]=0
##                            m+=1

        #Compute determinant of n-1 minor. Any n-1 minor should give the same answer, but that is not what I see. TODO fixed already?
        minor=am[:numunderpasses-1,:numunderpasses-1]
        #minor=am[1:numunderpasses,1:numunderpasses] #This one is good too.
        #print am
        #print minor
        det=linalg.det(minor)
        #eig=linalg.eig(minor)
        #print eig
        #(sign, logdet)=linalg.slogdet(am[:numunderpasses-1,:numunderpasses-1])
        #print "slogdet:", sign, logdet
        #det=numpy.exp(logdet)
        if det<0:
            det=-det
        return det

##b-1,a-2,b-3,a-1,b-2,a-3
##
##b-1,a-2,b+3,a+4,b-2,a-1,b+4,a+3
##
##a-1,b-2,a-3,b-4,a-5,b-1,a-2,b-3,a-4,b-5
##
##a-1,b-2,a-3,b-4,a-5,b-1,a-2,b-5,a-4,b-3
##
##a+1,b+2,a-3,b-4,a-5,b-6,a+2,b+1,a-6,b-5,a-4,b-3
##
##b+1,a+2,b-3,a-4,b-5,a-6,b+2,a+1,b-6,a-3,b-4,a-5
##
##a-1,b-2,a+3,b+4,a+5,b+3,a-6,b-1,a-2,b-6,a+4,b+5
##
##a-1,b-2,a-3,b-4,a-5,b-6,a-7,b-1,a-2,b-3,a-4,b-5,a-6,b-7
##
##a+1,b+2,a-3,b-4,a-5,b-6,a-7,b-8,a+2,b+1,a-8,b-7,a-6,b-5,a-4,b-3
##
##a+1,b-2,b+3,a-4,b-5,b-6,a-7,b+1,b-8,a+3,a-6,b-7,a-2,a-8,b-4,a-5
##
##a-1,b-2,a-3,b-4,a-5,b-6,a-7,b-8,a-9,b-1,a-2,b-3,a-4,b-5,a-6,b-7,a-8,b-9
##
##a+1,b-2,a-3,b-4,a-5,a-6,b-7,b-8,a-2,a-9,b-10,b-3,a-4,b-5,a-8,b+1,b-9,a-10,b-6,a-7

    #Thd computation of the vassiliev invariants is a translation of the gauss diagrams formulas
    # in the paper by Polyak and Viro.
    #The input is a gauss code, e.g. b-1,a-2,b+3,a+4,b-2,a-1,b+4,a+3
    #ab contains the sequence of 'a','b' (above/below or over/under)
    #num contains the crossing labels (numbers) with sign
    def computeVassiliev2(self,ab,num):
	# Vassiliev Invariant of Degree 2
	# pattern: b1 B b2 a1 a2
	# place base point B at end of Gauss code
	(sum,reps)=(0,0)
	#(ci,ci2,ci3,ci4,se,pc,pc2,pc3,pc4)=(0,0,0,0,0,0,0,0,0)
	se=len(ab)

	# translated from C++ to PERL to Python

        #for($ci=0;$ci<$se;$ci++)
        #{
	ci=0
	while ci<se:
            pc=ci
            reject=False
            if ab[pc]=='b':	#b1
                #for($ci2=$ci+1;$ci2<$se;$ci2++)	#scan code up to end to check for appearance of a1
                #{
                ci2=ci+1
                while ci2<se:
                    pc2=ci2
                    if ab[pc2]=='a' and num[pc2]==num[pc]:
                        reject=True
                        break
                    ci2+=1

                if not reject:
                    #for($ci2=0;$ci2<$ci;$ci2++)	#search for a b2
                    #{
                    ci2=0
                    while ci2<ci:
                        pc2=ci2;
                        if ab[pc2]=='b':	#b2
                            #for($ci3=$ci2+1;$ci3<$ci;$ci3++)	#search for an a1
                            #{
                            ci3=ci2+1
                            while ci3<ci:
                                pc3=ci3
                                if ab[pc3]=='a':
                                    if num[pc3]==num[pc2]:
                                        break	# reject
                                    elif num[pc3]==num[pc]:	#a1
                                        #for($ci4=$ci3+1;$ci4<$ci;$ci4++)	#search for an a2
                                        #{
                                        ci4=ci3+1
                                        while ci4<ci:
                                            pc4=ci4
                                            if ab[pc4]=='a' and num[pc4]==num[pc2]:	#a2
                                                if num[pc]*num[pc2]>0:
                                                    sum+=1
                                                else:
                                                    sum-=1
                                                reps+=1
                                            ci4+=1
                                        #}
                                ci3+=1
                            #}
                        elif ab[pc2]=='a' and num[pc2]==num[pc]:
                            break	# reject
                        ci2+=1
                    #}
                #}
            #}
            ci+=1
        #}

	#$v2=$sum;
	#$arf=$reps%2;
	return (sum,reps%2)

    # Vassiliev Invariant of Degree 3, translated from C++
    # pattern1: b1 b2 a1 b3 a2 a3
    # pattern2: b1 a2 b3 a1 b2 a3
    # v3: 1/2(pattern1) + pattern2
    def computeVassiliev3(self,ab,num):
	(sum1,sum2)=(0,0)
	#my ($ci,$ci2,$ci3,$ci4,$ci5,$ci6,$sb,$se);
	#my ($pc,$pc2,$pc3,$pc4,$pc5,$pc6);
	sb=0
	se=len(ab)

        # pattern 1

        #for($ci=$sb;$ci<$se;$ci++)	#search for a b1
        #{
        ci=sb
        while ci<se:
            pc=ci
            if ab[pc]=='b':	#b1
                #for($ci2=$ci+1;$ci2!=$ci;$ci2++)	#search for a b2
                #{
                ci2=ci+1
                while ci2!=ci:
                    if ci2==se:
                        if ci==sb:
                            break
                        ci2=sb
                    pc2=ci2
                    if ab[pc2]=='b':	#b2
                        #for($ci3=$ci2+1;$ci3!=$ci;$ci3++)	#search for a1
                        #{
                        ci3=ci2+1
                        while ci3!=ci:
                            if ci3==se:
                                if ci==sb:
                                    break
                                ci3=sb
                            pc3=ci3
                            if ab[pc3]=='a':	#a1
                                if num[pc3]==num[pc]:
                                    #for($ci4=$ci3+1;$ci4!=$ci;$ci4++)	#search for a b3
                                    #{
                                    ci4=ci3+1
                                    while ci4!=ci:
                                        if ci4==se:
                                            if ci==sb:
                                                break
                                            ci4=sb
                                        pc4=ci4
                                        if ab[pc4]=='b':	#b3
                                            #for($ci5=$ci4+1;$ci5!=$ci;$ci5++)	#search for a2
                                            #{
                                            ci5=ci4+1
                                            while ci5!=ci:
                                                if ci5==se:
                                                    if ci==sb:
                                                        break
                                                    ci5=sb
                                                pc5=ci5
                                                if ab[pc5]=='a':	#a2
                                                    if num[pc5]==num[pc2]:
                                                        #for($ci6=$ci5+1;$ci6!=$ci;$ci6++)	#search for a3
                                                        #{
                                                        ci6=ci5+1
                                                        while ci6!=ci:
                                                            if ci6==se:
                                                                if ci==sb:
                                                                    break
                                                                ci6=sb
                                                            pc6=ci6
                                                            if ab[pc6]=='a' and num[pc6]==num[pc4]:	#a3
                                                                #$sum1+=(($num[$pc]*$num[$pc2]*$num[$pc4])>0?1:-1);
                                                                if num[pc]*num[pc2]*num[pc4]>0:
                                                                    sum1+=1
                                                                else:
                                                                    sum1-=1
                                                            ci6+=1
                                                        #}
                                                    elif num[pc5]==num[pc4]:
                                                        break	#reject
                                                #}
                                                ci5+=1
                                            #}
                                        #}
                                        ci4+=1
                                    #}
                                #}
                                elif num[pc3]==num[pc2]:
                                    break	#reject
                            #}
                            ci3+=1
                        #}
                    #}
                    elif ab[pc2]=='a' and num[pc2]==num[pc]:
                        break	#reject
                    ci2+=1
                #}
            #}
            ci+=1
        #}

        # pattern 2

        #for($ci=0;$ci<$se;$ci++)	#search for a b1
        #{
        ci=0
        while ci<se:
            pc=ci
            if ab[pc]=='b':	#b1
                #for($ci2=$ci+1;$ci2!=$ci;$ci2++)	#search for an a2
                #{
                ci2=ci+1
                while ci2!=ci:
                    if ci2==se:
                        if ci==sb:
                            break
                        ci2=sb
                    pc2=ci2
                    if ab[pc2]=='a':	#a2
                        if num[pc2]==num[pc]:
                            break	#reject
                        #for($ci3=$ci2+1;$ci3!=$ci;$ci3++)	#search for a b3
                        #{
                        ci3=ci2+1
                        while ci3!=ci:
                            if ci3==se:
                                if ci==sb:
                                    break
                                ci3=sb
                            pc3=ci3
                            if ab[pc3]=='b':	#b3
                                if num[pc3]==num[pc2]:
                                    break	#reject
                                #for($ci4=$ci3+1;$ci4!=$ci;$ci4++)	#search for a1
                                #{
                                ci4=ci3+1
                                while ci4!=ci:
                                    if ci4==se:
                                        if ci==sb:
                                            break
                                        ci4=sb
                                    pc4=ci4
                                    if ab[pc4]=='a':	#a1
                                        if num[pc4]==num[pc3]:
                                            break	#reject
                                        elif num[pc4]==num[pc]:
                                            #for($ci5=$ci4+1;$ci5!=$ci;$ci5++)	#search for b2
                                            #{
                                            ci5=ci4+1
                                            while ci5!=ci:
                                                if ci5==se:
                                                    if ci==sb:
                                                        break
                                                    ci5=sb
                                                pc5=ci5
                                                if ab[pc5]=='b' and num[pc5]==num[pc2]:	#b2
                                                    #for($ci6=$ci5+1;$ci6!=$ci;$ci6++)	#search for a3
                                                    #{
                                                    ci6=ci5+1
                                                    while ci6!=ci:
                                                        if ci6==se:
                                                            if ci==sb:
                                                                break
                                                            ci6=sb
                                                        pc6=ci6
                                                        if ab[pc6]=='a' and num[pc6]==num[pc3]:	#a3
                                                            #$sum2+=(($num[$pc]*$num[$pc2]*$num[$pc3])>0?1:-1);
                                                            if num[pc]*num[pc2]*num[pc3]>0:
                                                                sum2+=1
                                                            else:
                                                                sum2-=1
                                                        ci6+=1
                                                    #}
                                                #}
                                                ci5+=1
                                            #}
                                        #}
                                    #}
                                    ci4+=1
                                #}
                            #}
                            ci3+=1
                        #}
                    #}
                    ci2+=1
                #}
            #}
            ci+=1
        #}

	v3=sum1/2+sum2/3	# the 1/3 is due to the three-fold symmetry of the chord diagram for pattern 2
	#print "Signed v3", v3
	#if v3<0:
        #    v3=-v3
        return v3 #Sign of v3 flip when the mirror image of the knot is taken

#############################################################################

class LinkControlGroup:
    def __init__(self,
                 page,
                 groupname='Link Analyzer',
                 defaultstructurename1='1ihf',
                 defaultclosureoption1=0,
                 defaultchain1='C',
                 defaultstructurename2='1ihf',
                 defaultclosureoption2=0,
                 defaultchain2='E'):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group
        #self.groupscrolled=group
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.structureframe=Tkinter.Frame(group.interior())
        #Field for entering name of structure or model
        self.structure1 = Pmw.EntryField(self.structureframe,
                                        labelpos='w',
                                        label_text='1st structure to use: ',
                                        value=defaultstructurename1,
                                        )

        self.chain1 = Pmw.EntryField(self.structureframe,
                                    labelpos='w',
                                    label_text='chain: ',
                                    value=defaultchain1,
                                    )
        self.structure2 = Pmw.EntryField(self.structureframe,
                                        labelpos='w',
                                        label_text='2nd structure to use: ',
                                        value=defaultstructurename2,
                                        )

        self.chain2 = Pmw.EntryField(self.structureframe,
                                    labelpos='w',
                                    label_text='chain: ',
                                    value=defaultchain2,
                                    )

        self.structure1.grid(column=0,row=0)
        self.chain1.grid(column=1,row=0)
        self.structure2.grid(column=0,row=1)
        self.chain2.grid(column=1,row=1)

        self.analyze_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.analyze_buttonbox.add('Analyze link',command=self.analyzeLink)

        self.advancedgroup = Pmw.Group(group.interior(),tag_text='Fine-tune your analysis')

        self.backboneatom_options_tuple = ('CA (Protein)',
                                            'P (DNA/RNA)')
        self.show_backboneatom_options1 = Pmw.OptionMenu(self.advancedgroup.interior(),
                                              labelpos = 'w',
                                              label_text = '1st structure backbone atom',
                                              items = self.backboneatom_options_tuple,
                                              initialitem = self.backboneatom_options_tuple[1],
                                              )
        self.show_backboneatom_options2 = Pmw.OptionMenu(self.advancedgroup.interior(),
                                              labelpos = 'w',
                                              label_text = '2nd structure backbone atom',
                                              items = self.backboneatom_options_tuple,
                                              initialitem = self.backboneatom_options_tuple[1],
                                              )

        #self.show_backboneatom_options1.pack(fill='x',padx=4,pady=1)
        #self.show_backboneatom_options2.pack(fill='x',padx=4,pady=1)
        self.show_backboneatom_options1.grid(column=0,row=0)
        self.show_backboneatom_options2.grid(column=0,row=1)

        self.altloc_var=IntVar()
        self.altloc_var.set(0)
        self.altloc_checkbutton = Checkbutton(self.advancedgroup.interior(),
                                                     text = "Use first atom variant when alternate location indicators exist",
                                                     variable = self.altloc_var)
        #self.altloc_checkbutton.pack(fill='x',padx=4,pady=1) #TODO Not compatible with .grid!??????
        self.altloc_checkbutton.grid(column=0,row=2)
        
        for entry in (self.structureframe,
                      self.analyze_buttonbox,
                      self.advancedgroup):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        #TODO balloons (tooltips)

    def analyzeLink(self):
        #Use directclosure by default
        m1=self.getModel(self.structure1.getvalue().strip(),
                                               self.chain1.getvalue().strip())
        m2=self.getModel(self.structure2.getvalue().strip(),
                                               self.chain2.getvalue().strip())

        if m1 and m2:
            composite_structurename1,model1=m1
            composite_structurename2,model2=m2
        else:
            return

        print "BEGIN LINK ANALYSIS OF %s AND %s" % (composite_structurename1,
                                                    composite_structurename2)

        #TODO: Restrict to amino acids
        backboneatom_option1=self.show_backboneatom_options1.getvalue()
        backboneatomname1='CA'
        if backboneatom_option1==self.backboneatom_options_tuple[0]:
            backboneatomname1='CA'
        elif backboneatom_option1==self.backboneatom_options_tuple[1]:
            backboneatomname1='P'
        backboneatom_option2=self.show_backboneatom_options2.getvalue()
        backboneatomname2='CA'
        if backboneatom_option2==self.backboneatom_options_tuple[0]:
            backboneatomname2='CA'
        elif backboneatom_option2==self.backboneatom_options_tuple[1]:
            backboneatomname2='P'

        backbone1=self.getBackbone(model1,backboneatomname1)
        backbone2=self.getBackbone(model2,backboneatomname2)
        if backbone1 and backbone2:
            #Calculate linking number
            N1=len(backbone1)
            N2=len(backbone2)
            if N1>1 and N2>1:
                Lk=0
                for i in range(N1):
                    for j in range(N2):
                        inext=i+1
                        if inext==N1:
                            inext=0
                        jnext=j+1
                        if jnext==N2:
                            jnext=0
                        crossingdata=self.findIntersection(backbone1[i],
                                                           backbone1[inext],
                                                           backbone2[j],
                                                           backbone2[jnext])
                        if crossingdata is not None:
                            (x,y,z1,z2,sign)=crossingdata
                            Lk+=sign
                Lk*=0.5
                print "Linking number is ", Lk


        print "END LINK ANALYSIS OF %s AND %s" % (composite_structurename1,
                                                  composite_structurename2)

    def getModel(self,structurename,chainindicator,startresnum='',endresnum=''):
        #Get user-provided PyMOL structure selection.
        #Get structurename and check if it exists in PyMOL.
        #structurename=self.structure1.getvalue().strip()
        supper=structurename.upper()
        for n in cmd.get_names('all'):
            if supper==n.upper():
                break
        else:
            print 'structure name must be in PyMOL viewer list:'
            print cmd.get_names('all')
            return

        if len(structurename)<1:
            print 'Provide structure name!'
            return
        try:
            model=cmd.get_model(structurename)
        except Exception:
            print "Structure \'%s\' does not exist!" % structurename
            return
        
        #chainindicator=self.chain.getvalue().strip() #TODO add chain to structurename?
        #startresnum=self.startres.getvalue().strip()
        #endresnum=self.endres.getvalue().strip()

        composite_structurename=structurename
        if len(chainindicator)>0:
            char_chainindicator=chainindicator[0]
        else:
            char_chainindicator=''

        if len(char_chainindicator)>0:
            if char_chainindicator!='?':
                composite_structurename='%s and chain %s' % (structurename,char_chainindicator)
            else:
                composite_structurename=structurename
        else:
            #Chain is blank, check if there are multiple chains in the structure
            #structurename=self.structurename
##            chaindict={}
##            try:
##                model=cmd.get_model(structurename)
##            except Exception:
##                print "Structure \'%s\' does not exist!" % structurename
##                return
##            for atom in model.atom:
##                chaindict[atom.chain]=1
##            if len(chaindict)>1:
            chainlist=cmd.get_chains(structurename)
            if len(chainlist)>1:
                print "Please select one from the following chains of %s:" % structurename, ','.join(chainlist)
                return

        if len(startresnum)>0 or len(endresnum)>0:
            composite_structurename+=' and resi %s-%s' % (startresnum,endresnum)
        try:
            model=cmd.get_model(composite_structurename)
        except Exception:
            print "Structure \'%s\' does not exist!" % composite_structurename
            return
        
        return (composite_structurename,model)

    #Get coordinates of backbone
    def getBackbone(self,model,backboneatomname):
        use_firstaltloc=self.altloc_var.get()
        numCA=0
        firstresnum=None
        lastresnum=None
        for atom in model.atom:
            if atom.name==backboneatomname: #Use alpha-carbon for protein. For DNA, use 'P'.
                #if atom.__dict__.has_key("alt") and \
                if use_firstaltloc==1 and \
                   atom.alt!='' and atom.alt!='A': #Accept the first alternate location for the atom and ignore others
                    continue
                if firstresnum is None:
                    firstresnum=atom.resi
                lastresnum=atom.resi
                numCA+=1
        #self.firstresnum=firstresnum
        #self.lastresnum=lastresnum

        if numCA>0:
            backbone=[None]*numCA
            #self.ss=[None]*numCA
            i=0
            ssnum=0
            #prevss=None
            previntresi=None
            break_resi=False
            breakstop=[]
            for atom in model.atom:
                if atom.name==backboneatomname:
                    #if atom.__dict__.has_key("alt") and \
                    if use_firstaltloc==1 and \
                       atom.alt!='' and atom.alt!='A': #Accept the first alternate location for the atom and ignore others
                        continue
                    backbone[i]=atom.coord
                    #TODO Add option to change projection planes (x-y, y-z, z-x) by permutation, or use rotate axis, angle, selection.
                    #self.backbone[i]=(atom.coord[1],atom.coord[2],atom.coord[0])
##                    if atom.ss!=prevss: #No need to call cmd.dss to make sure. DNA (P backbone) has ss==''.
##                        ssnum+=1
##                    self.ss[i]=(atom.ss,ssnum) #Store secondary structure to help with knot projection simplification
##                    prevss=atom.ss
                    try:
                        intresi=int(atom.resi)
                    except ValueError:
                        intresi=int(atom.resi[:-1])
                    if previntresi is not None:
                        if intresi-previntresi>1:
                            print 'Warning: Skipped residue numbers (between %d and %d)! Artificial bond will be added to the backbone.' % (previntresi,intresi)
                            break_resi=True
                            breakstop.append(i)
                    previntresi=intresi
                    i+=1
            return backbone
        else:
            if numCA==0: #Warning: atom named P exists in protein 2efv
                print "No %s atoms found!" % backboneatomname
            #backbone=None
            return None

    #Borrow from KAControlGroup
    #findIntersection=KAControlGroup.findIntersection
    #isBetween=KAControlGroup.isBetween

    #Iteration one: brute force
    #1. Find intersection of x-y projection between segment r1-r2 and segment r3-r4. r1 is a vector (x,y,z).
    #2. If they intersect, find the z-coordinate of each segment (used to find which segment goes over(under) the other).
    #3. Return sign of the crossing.
    def findIntersection(self,r1,r2,r3,r4):
        #Solve for x and y in
        #y=m1*(x-x1)+y1, m1=(y2-y1)/(x2-x1)
        #y=m2*(x-x3)+y3, m2=(y4-y3)/(x4-x3)
        dx1=r2[0]-r1[0]
        dy1=r2[1]-r1[1]
        dx2=r4[0]-r3[0]
        dy2=r4[1]-r3[1]

        #Handle exceptional cases
        #TODO Check that the case when the segment endpoints may coincide is handled properly
        if dx1==0 and dx2==0:
            if r1[0]!=r3[0]:
                return None
            print "Warning: Encountered two overlapping vertical segments! Perturbing one segment..."
            r1[0]+=0.01
            #r3[0]+=0.01
            dx1=r2[0]-r1[0]
            dx2=r4[0]-r3[0] #TODO r3 not modified anymore, so why bother?
        if dy1==0 and dy2==0:
            if r1[1]!=r3[1]:
                return None
            print "Warning: Encountered two overlapping horizontal segments! Perturbing one segment..."
            r1[1]+=0.01
            #r3[1]+=0.01
            dy1=r2[1]-r1[1]
            dy2=r4[1]-r3[1] #TODO r3 not modified anymore, so why bother?
        if dx1==0 and dy2==0:
            x=r1[0]
            y=r3[1]
        elif dy1==0 and dx2==0:
            x=r3[0]
            y=r1[1]
        elif dx1==0:
            x=r1[0]
            y=(x-r3[0])*dy2/dx2+r3[1]
        elif dx2==0:
            x=r3[0]
            y=(x-r1[0])*dy1/dx1+r1[1]
        elif dy1==0:
            y=r1[1]
            x=(y-r3[1])*dx2/dy2+r3[0]
        elif dy2==0:
            y=r3[1]
            x=(y-r1[1])*dx1/dy1+r1[0]
        else:
            m1=dy1/dx1
            m2=dy2/dx2
            #What if the slopes are the same? (parallel segments)
            if m1==m2:
                print "Warning: Encountered two parallel segments!  Perturbing one segment..."
                r1[0]+=0.01
                dx1=r2[0]-r1[0]
                m1=dy1/dx1
            x=(m1*r1[0]-m2*r3[0]+r3[1]-r1[1])/(m1-m2) #TODO reverse roles of x and y for better numerical precision?
            y=m1*(x-r1[0])+r1[1]
            #y=m2*(x-r3[0])+r3[1] #TODO take the average of this with the other one
        #Check that (x,y) lies on the segments for proper intersection
        if self.isBetween(x,r1[0],r2[0]) and self.isBetween(y,r1[1],r2[1]):
            if self.isBetween(x,r3[0],r4[0]) and self.isBetween(y,r3[1],r4[1]):
                if dx1!=0:
                    z1=(r2[2]-r1[2])*(x-r1[0])/dx1+r1[2]
                else:
                    z1=(r2[2]-r1[2])*(y-r1[1])/dy1+r1[2]
                if dx2!=0:
                    z2=(r4[2]-r3[2])*(x-r3[0])/dx2+r3[2]
                else:
                    z2=(r4[2]-r3[2])*(y-r3[1])/dy2+r3[2]
                #Use cross product to get the sign of the crossing
                if dx1*dy2-dx2*dy1>0:
                    if z1>z2:
                        sign=-1
                    else:
                        sign=+1
                else: #TODO check if cross product is zero
                    if z1>z2:
                        sign=+1
                    else:
                        sign=-1
                return (x,y,z1,z2,sign)
        return None

    def isBetween(self,x,x1,x2):
        if x1<x2:
            if x1<=x and x<=x2:
                return True
        else: #if x2<x1:
            if x2<=x and x<=x1:
                return True
        return False

#############################################################################

class KnotCreatorControlGroup:
    def __init__(self,
                 page,
                 groupname='Knot Creator'):
        group = Pmw.Group(page,tag_text=groupname)
##        group=Pmw.ScrolledFrame(page,
##                                labelpos='nw',
##                                label_text=groupname)
        self.groupname=groupname
        self.group=group
        #self.groupscrolled=group
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.notebook = Pmw.NoteBook(group.interior())
        self.notebook.pack(fill='both',expand=1,padx=4,pady=1)

        page=self.notebook.add('Torus')
        kc1=TorusKnotControlGroup(page)
        self.balloon.bind(page,"Create a backbone on the surface of a torus. The backbone atoms are created\n with PyMOL's pseudoatom command, available in version 1.3.")

        page=self.notebook.add('Twist')
        kc2=TwistKnotControlGroup(page)
        self.balloon.bind(page,"Create a backbone forming a twist knot. The backbone atoms are created\n with PyMOL's pseudoatom command, available in version 1.3.")

        page=self.notebook.add('Hilbert 2D')
        kc3=Hilbert2DControlGroup(page)
        self.balloon.bind(page,"Create a backbone forming a 2D Hilbert curve. The backbone atoms are created\n with PyMOL's pseudoatom command, available in version 1.3.")

        page=self.notebook.add('Hilbert 3D')
        kc4=Hilbert3DControlGroup(page)
        self.balloon.bind(page,"Create a backbone forming a 3D Hilbert curve. The backbone atoms are created\n with PyMOL's pseudoatom command, available in version 1.3.")

        page=self.notebook.add('Peano 2D')
        kc5=Peano2DControlGroup(page)
        self.balloon.bind(page,"Create a backbone forming a 2D Peano curve. The backbone atoms are created\n with PyMOL's pseudoatom command, available in version 1.3.")



class TorusKnotControlGroup:
    BACKBONE_RADIUS=1
    def __init__(self,
                 page,
                 groupname='Create a torus knot T(M,N)'):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group

        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.M = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of turns through the tube (M):',
                                 validate={'validator':'integer','min':1,'max':100},
                                 value='2',
                                 )
        
        self.N = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of turns around the axis of the tube (N):',
                                 validate={'validator':'integer','min':1,'max':100},
                                 value='3',
                                 )

        self.tuberadius = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Tube radius:',
                                 validate={'validator':'real','min':1,'max':100},
                                 value='2',
                                 )

        self.holeradius = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Hole radius:',
                                 validate={'validator':'real','min':1,'max':100},
                                 value='4',
                                 )

        self.numsegments = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of segments:',
                                 validate={'validator':'integer','min':1,'max':1000},
                                 value='33',
                                 )

        self.create_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.create_buttonbox.add('Create knot',command=self.createTorusKnot)
        
##        self.M.grid(column=0,row=0)
##        self.N.grid(column=0,row=1)
##        self.tuberadius.grid(column=0,row=2)
##        self.holeradius.grid(column=0,row=3)

        for entry in (self.M,
                      self.N,
                      self.tuberadius,
                      self.holeradius,
                      self.numsegments,
                      self.create_buttonbox):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        #TODO balloons (tooltips)

    def createTorusKnot(self):
        #TODO When numsegments in even, the knot analysis is wrong sometimes!?
        N=int(self.N.getvalue())
        M=int(self.M.getvalue())
        R1=float(self.tuberadius.getvalue())
        R2=float(self.holeradius.getvalue())
        numsegments=int(self.numsegments.getvalue())
        knotname="TORUS_%d_%d" % (M,N)

        gcd=self.GCD(M,N)
        if gcd>1:
            print "Warning! GCD(%d,%d)=%d" % (M,N,gcd)

        backbone=[None]*numsegments
        pN=2*math.pi*N
        pM=2*math.pi*M
        for i in range(numsegments):
            t=i*1.0/numsegments
            t1=pN*t
            t2=pM*t
            r=R2+R1*(1+math.cos(t1))
            backbone[i]=(r*math.cos(t2),r*math.sin(t2),R1*math.sin(t1))

        self.drawBackbone(backbone,knotname+"_TUBES")

        #Create PyMOL object that can be accessed by PyMOL's get_model
        #pseudoatom is not available in version 0.99rc!
        knotname_atoms=knotname+"_ATOMS"
        cmd.delete(knotname_atoms)
        for i in range(numsegments):
            cmd.pseudoatom(knotname_atoms,pos=tuple(backbone[i]),
                           resi="%d" % (i+1),name="CA",hetatm=0,vdw=1,segi="",b=1.0,q=1.0,elem="C")

        print "Knot Creator has created a backbone with %d atoms" % len(backbone)

    def GCD(self,M,N):
        if M>N:
            b=M
            r=N
        else:
            b=N
            r=M
        while r>0:
            s=r
            r=b%r
            b=s
        return b
    
    def drawBackbone(self,backbone,knotname):
        i=0
        obj=[]
        while i+1<len(backbone):
            x1=backbone[i][0]
            y1=backbone[i][1]
            z1=backbone[i][2]
            x2=backbone[i+1][0]
            y2=backbone[i+1][1]
            z2=backbone[i+1][2]
            obj.extend([SPHERE, x1, y1, z1, self.BACKBONE_RADIUS,
                        CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                       1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
            i+=1
        cmd.delete(knotname)
        cmd.load_cgo(obj,knotname)

class TwistKnotControlGroup:
    BACKBONE_RADIUS=1
    def __init__(self,
                 page,
                 groupname='Create a twist knot'):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group

        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.N = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of half twists:',
                                 validate={'validator':'integer','min':1,'max':100},
                                 value='1',
                                 )
        
        self.tuberadius = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Tube radius:',
                                 validate={'validator':'real','min':1,'max':100},
                                 value='4',
                                 )

        self.numsegments = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of segments:',
                                 validate={'validator':'integer','min':1,'max':1000},
                                 value='10',
                                 )

        self.create_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.create_buttonbox.add('Create knot',command=self.createTwistKnot)
        
##        self.M.grid(column=0,row=0)
##        self.N.grid(column=0,row=1)
##        self.tuberadius.grid(column=0,row=2)
##        self.holeradius.grid(column=0,row=3)

        for entry in (self.N,
                      self.tuberadius,
                      self.numsegments,
                      self.create_buttonbox):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        #TODO balloons (tooltips)

    def createTwistKnot(self):
        N=int(self.N.getvalue())
        R1=float(self.tuberadius.getvalue())
        numsegments=int(self.numsegments.getvalue())
        knotname="TWIST_%d" % (N)

        #1. Create the first backbone of the central helix (direction to be reversed)
        backbone1=[None]*numsegments
        pN=math.pi*N
        tubelength=2*N*R1
        for i in range(numsegments):
            t=(i+1)*1.0/numsegments
            backbone1[i]=(R1*math.cos(pN*t),tubelength*t,R1*math.sin(pN*t))

        #self.drawBackbone(backbone1,knotname+"_TUBES1")

        #2. Create the second backbone of the central helix
        backbone2=[None]*numsegments
        pN=math.pi*N
        tubelength=2*N*R1
        for i in range(numsegments):
            t=i*1.0/numsegments
            backbone2[i]=(-R1*math.cos(pN*t),tubelength*t,-R1*math.sin(pN*t))

        #self.drawBackbone(backbone2,knotname+"_TUBES2")

        #3. Create the small connecting arc
        n=numsegments/N
        backbone3=[None]*n
        w=math.pi
        for i in range(n):
            t=i*1.0/n
            backbone3[i]=(R1*math.cos(w*t),-R1*math.sin(w*t),0)

        #self.drawBackbone(backbone3,knotname+"_TUBES3")

        #4. Create the short segment of the first closing arc (direction to be reversed if N%2==0)
        backbone4=[None]*(n/3)
        for i in range(n/3):
            t=(i+1)*3.0/n
            backbone4[i]=(R1*(1+t),tubelength,0)

        #self.drawBackbone(backbone4,knotname+"_TUBES4")

        #5. Create the short segment of the second closing arc (direction to be reversed if N%2==1)
        backbone5=[None]*(n/3)
        for i in range(n/3):
            t=(i+1)*3.0/n
            backbone5[i]=(-R1*(1+t),tubelength,0)

        #self.drawBackbone(backbone5,knotname+"_TUBES5")

        #6. Create the first closing arc (direction to be reversed if N%2==0)
        backbone6=[None]*numsegments
        w=math.pi
        for i in range(numsegments):
            t=(i+1)*1.0/(numsegments+1)
            backbone6[i]=(R1+R1*math.cos(w*t),tubelength*(1-t),2*N*R1*math.sin(w*t))

        #self.drawBackbone(backbone6,knotname+"_TUBES6")

        #7. Create the second closing arc (direction to be reversed if N%2==1)
        backbone7=[None]*numsegments
        w=math.pi
        for i in range(numsegments):
            t=(i+1)*1.0/(numsegments+1)
            backbone7[i]=(-R1-R1*math.cos(w*t),tubelength*(1-t),-2*N*R1*math.sin(w*t))

        #self.drawBackbone(backbone7,knotname+"_TUBES7")
        
        backbone=[None]*len(backbone1)
        #Read backbone1 in reverse
        i=0
        for b in reversed(backbone1):
            backbone[i]=b
            i+=1
        backbone.extend(backbone3)
        backbone.extend(backbone2)
        if N%2==1:
            backbone.extend(backbone4)
            backbone.extend(backbone6)
            for b in reversed(backbone7):
                backbone.append(b)
            for b in reversed(backbone5):
                backbone.append(b)
        else:
            backbone.extend(backbone5)
            backbone.extend(backbone7)
            for b in reversed(backbone6):
                backbone.append(b)
            for b in reversed(backbone4):
                backbone.append(b)

        self.drawBackbone(backbone,knotname+"_TUBES")

        #Create PyMOL object that can be accessed by PyMOL's get_model
        #pseudoatom is not available in version 0.99rc!
        knotname_atoms=knotname+"_ATOMS"
        cmd.delete(knotname_atoms)
        for i in range(len(backbone)):
            cmd.pseudoatom(knotname_atoms,pos=tuple(backbone[i]),
                           resi="%d" % (i+1),name="CA",hetatm=0,vdw=1,segi="",b=1.0,q=1.0,elem="C")

        print "Knot Creator has created a backbone with %d atoms" % len(backbone)

    def drawBackbone(self,backbone,knotname):
        i=0
        obj=[]
        while i+1<len(backbone):
            x1=backbone[i][0]
            y1=backbone[i][1]
            z1=backbone[i][2]
            x2=backbone[i+1][0]
            y2=backbone[i+1][1]
            z2=backbone[i+1][2]
            obj.extend([SPHERE, x1, y1, z1, self.BACKBONE_RADIUS,
                        CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                       1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
            i+=1
        cmd.delete(knotname)
        cmd.load_cgo(obj,knotname)

class Hilbert2DControlGroup:
    BACKBONE_RADIUS=1
    def __init__(self,
                 page,
                 groupname='Create a 2D Hilbert curve'):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group

        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.N = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Order (Number of iterations):',
                                 validate={'validator':'integer','min':1,'max':10},
                                 value='2',
                                 )
        
        self.stepsize = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Step size:',
                                 validate={'validator':'real','min':1,'max':100},
                                 value='2',
                                 )

        self.numatoms = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of atoms per segment:',
                                 validate={'validator':'integer','min':1,'max':10},
                                 value='1',
                                 )

        self.create_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.create_buttonbox.add('Create curve',command=self.createHilbert2D)

        for entry in (self.N,
                      self.stepsize,
                      self.numatoms,
                      self.create_buttonbox):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        self.rotateCWmap={'N':'E',
                          'E':'S',
                          'S':'W',
                          'W':'N'}
        self.rotateCCWmap={'E':'N',
                          'S':'E',
                          'W':'S',
                          'N':'W'}
        self.reversemap={'N':'S',
                          'S':'N',
                          'E':'W',
                          'W':'E'}

    def createHilbert2D(self):
        N=int(self.N.getvalue())
        stepsize=float(self.stepsize.getvalue())
        numatoms=int(self.numatoms.getvalue())
        knotname="HILBERT2D_%d" % (N)

        self.step={'N':(0,stepsize),
                   'S':(0,-stepsize),
                   'E':(stepsize,0),
                   'W':(-stepsize,0)}

        sequence=['N','E','S']
        while N>1:
            tmp=self.rotateCW(sequence)
            s1=self.reverse(tmp)
            s2=sequence
            s3=sequence
            tmp=self.rotateCCW(sequence)
            s4=self.reverse(tmp)
            sequence=s1+['N']+s2+['E']+s3+['S']+s4
            N-=1
            
        backbone=[(0,0,0)]
        for s in sequence:
            current=backbone[-1]
            backbone.append((current[0]+self.step[s][0],
                             current[1]+self.step[s][1],
                             current[2]))

        self.drawBackbone(backbone,knotname+'_TUBES')

        #Create PyMOL object that can be accessed by PyMOL's get_model
        #pseudoatom is not available in version 0.99rc!
        knotname_atoms=knotname+"_ATOMS"
        cmd.delete(knotname_atoms)
        atomcount=0
        for i in range(len(backbone)-1):
            #Interpolate the coordinates of atoms to place between "joints"
            for j in range(numatoms): #Number of atoms per segment or step
                x=backbone[i][0]+self.step[sequence[i]][0]*j*1.0/numatoms
                y=backbone[i][1]+self.step[sequence[i]][1]*j*1.0/numatoms
                z=backbone[i][2]
                #print x,y,z
                atomcount+=1
                cmd.pseudoatom(knotname_atoms,pos=(x,y,z),
                               resi="%d" % (atomcount),name="CA",hetatm=0,vdw=1,segi="",b=1.0,q=1.0,elem="C")

    def drawBackbone(self,backbone,knotname):
        i=0
        obj=[]
        while i+1<len(backbone):
            x1=backbone[i][0]
            y1=backbone[i][1]
            z1=backbone[i][2]
            x2=backbone[i+1][0]
            y2=backbone[i+1][1]
            z2=backbone[i+1][2]
            obj.extend([SPHERE, x1, y1, z1, self.BACKBONE_RADIUS,
                        CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                       1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
            i+=1
        cmd.delete(knotname)
        cmd.load_cgo(obj,knotname)

    def rotateCW(self,sequence):
        newsequence=[]
        for s in sequence:
            newsequence.append(self.rotateCWmap[s])
        return newsequence

    def rotateCCW(self,sequence):
        newsequence=[]
        for s in sequence:
            newsequence.append(self.rotateCCWmap[s])
        return newsequence

    def reverse(self,sequence):
        newsequence=[]
        sequence.reverse()
        for s in sequence:
            newsequence.append(self.reversemap[s])
        return newsequence

class Hilbert3DControlGroup:
    BACKBONE_RADIUS=1
    def __init__(self,
                 page,
                 groupname='Create a 3D Hilbert curve'):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group

        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.N = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Order (Number of iterations):',
                                 validate={'validator':'integer','min':1,'max':10},
                                 value='2',
                                 )
        
        self.stepsize = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Step size:',
                                 validate={'validator':'real','min':1,'max':100},
                                 value='2',
                                 )

        self.numatoms = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of atoms per segment:',
                                 validate={'validator':'integer','min':1,'max':10},
                                 value='1',
                                 )

        self.create_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.create_buttonbox.add('Create curve',command=self.createHilbert3D)

        for entry in (self.N,
                      self.stepsize,
                      self.numatoms,
                      self.create_buttonbox):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        self.transformmap1={'N':'F',
                         'S':'B',
                         'F':'S',
                         'B':'N',
                         'E':'E',
                         'W':'W'}
        self.transformmap2={'N':'F',
                         'E':'N',
                         'S':'B',
                         'F':'E',
                         'W':'S',
                         'B':'W'}
        self.transformmap3={'N':'S',
                         'E':'W',
                         'S':'N',
                         'F':'F',
                         'W':'E',
                         'B':'B'}
        self.transformmap4={'N':'N',
                         'E':'B',
                         'S':'S',
                         'F':'E',
                         'W':'F',
                         'B':'W'}
        self.transformmap5={'N':'B',
                         'E':'E',
                         'S':'F',
                         'F':'N',
                         'W':'W',
                         'B':'S'}
        self.reversemap={'N':'S',
                         'S':'N',
                         'E':'W',
                         'W':'E',
                         'F':'B',
                         'B':'F'}

    def createHilbert3D(self):
        N=int(self.N.getvalue())
        stepsize=float(self.stepsize.getvalue())
        numatoms=int(self.numatoms.getvalue())
        knotname="HILBERT3D_%d" % (N)

        self.step={'N':(0,stepsize,0),
                   'S':(0,-stepsize,0),
                   'E':(stepsize,0,0),
                   'W':(-stepsize,0,0),
                   'F':(0,0,stepsize),
                   'B':(0,0,-stepsize)}

        sequence=['N','E','S','F','N','W','S']
        while N>1:
            tmp=self.transform(sequence,1)
            s1=self.reverse(tmp)
            s2=self.transform(sequence,2)
            s3=s2
            s4=self.transform(sequence,3)
            s5=s4
            tmp=self.transform(sequence,4)
            s6=self.reverse(tmp)
            s7=s6
            tmp=self.transform(sequence,5)
            s8=self.reverse(tmp)
            sequence=s1+['N']+s2+['E']+s3+['S']+\
                      s4+['F']+s5+['N']+s6+['W']+\
                      s7+['S']+s8
            N-=1
            
        backbone=[(0,0,0)]
        for s in sequence:
            current=backbone[-1]
            backbone.append((current[0]+self.step[s][0],
                             current[1]+self.step[s][1],
                             current[2]+self.step[s][2]))

        self.drawBackbone(backbone,knotname+'_TUBES')

        #Create PyMOL object that can be accessed by PyMOL's get_model
        #pseudoatom is not available in version 0.99rc!
        knotname_atoms=knotname+"_ATOMS"
        cmd.delete(knotname_atoms)
        atomcount=0
        for i in range(len(backbone)-1):
            #Interpolate the coordinates of atoms to place between "joints"
            for j in range(numatoms): #Number of atoms per segment or step
                x=backbone[i][0]+self.step[sequence[i]][0]*j*1.0/numatoms
                y=backbone[i][1]+self.step[sequence[i]][1]*j*1.0/numatoms
                z=backbone[i][2]+self.step[sequence[i]][2]*j*1.0/numatoms
                #print x,y,z
                atomcount+=1
                cmd.pseudoatom(knotname_atoms,pos=(x,y,z),
                               resi="%d" % (atomcount),name="CA",hetatm=0,vdw=1,segi="",b=1.0,q=1.0,elem="C")

    def drawBackbone(self,backbone,knotname):
        i=0
        obj=[]
        while i+1<len(backbone):
            x1=backbone[i][0]
            y1=backbone[i][1]
            z1=backbone[i][2]
            x2=backbone[i+1][0]
            y2=backbone[i+1][1]
            z2=backbone[i+1][2]
            obj.extend([SPHERE, x1, y1, z1, self.BACKBONE_RADIUS,
                        CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                       1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
            i+=1
        cmd.delete(knotname)
        cmd.load_cgo(obj,knotname)

    def transform(self,sequence,option):
        newsequence=[]
        if option==1:
            for s in sequence:
                newsequence.append(self.transformmap1[s])
        elif option==2:
            for s in sequence:
                newsequence.append(self.transformmap2[s])
        elif option==3:
            for s in sequence:
                newsequence.append(self.transformmap3[s])
        elif option==4:
            for s in sequence:
                newsequence.append(self.transformmap4[s])
        elif option==5:
            for s in sequence:
                newsequence.append(self.transformmap5[s])
        return newsequence

    def reverse(self,sequence):
        newsequence=[]
        sequence.reverse()
        for s in sequence:
            newsequence.append(self.reversemap[s])
        return newsequence
        

class Peano2DControlGroup:
    BACKBONE_RADIUS=1
    def __init__(self,
                 page,
                 groupname='Create a 2D Peano curve'):
        #group = Pmw.Group(page,tag_text=groupname)
        group=Pmw.ScrolledFrame(page,
                                labelpos='nw',
                                label_text=groupname)
        self.groupname=groupname
        self.group=group

        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.balloon=Pmw.Balloon(group.interior())

        self.N = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Order (Number of iterations):',
                                 validate={'validator':'integer','min':1,'max':10},
                                 value='2',
                                 )
        
        self.stepsize = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Step size:',
                                 validate={'validator':'real','min':1,'max':100},
                                 value='2',
                                 )

        self.numatoms = Pmw.EntryField(group.interior(),
                                 labelpos='w',
                                 label_text='Number of atoms per segment:',
                                 validate={'validator':'integer','min':1,'max':10},
                                 value='1',
                                 )

        self.create_buttonbox = Pmw.ButtonBox(group.interior(), padx=0)
        self.create_buttonbox.add('Create curve',command=self.createPeano2D)

        for entry in (self.N,
                      self.stepsize,
                      self.numatoms,
                      self.create_buttonbox):
            entry.pack(fill='x',padx=4,pady=1) # vertical

        self.rotateCWmap={'N':'E',
                          'E':'S',
                          'S':'W',
                          'W':'N'}
        self.rotateCCWmap={'E':'N',
                          'S':'E',
                          'W':'S',
                          'N':'W'}
        self.reversemap={'N':'S',
                          'S':'N',
                          'E':'W',
                          'W':'E'}

    def createPeano2D(self):
        N=int(self.N.getvalue())
        stepsize=float(self.stepsize.getvalue())
        numatoms=int(self.numatoms.getvalue())
        knotname="PEANO2D_%d" % (N)

        self.step={'N':(0,stepsize),
                   'S':(0,-stepsize),
                   'E':(stepsize,0),
                   'W':(-stepsize,0)}

            
        sequence=self.Peano(N)
        sequence=sequence.replace('U','N')
        sequence=sequence.replace('D','S')
        sequence=sequence.replace('R','E')
        sequence=sequence.replace('L','W')
        #print sequence

        backbone=[(0,0,0)]
        for s in sequence:
            current=backbone[-1]
            backbone.append((current[0]+self.step[s][0],
                             current[1]+self.step[s][1],
                             current[2]))

        self.drawBackbone(backbone,knotname+'_TUBES')

        #Create PyMOL object that can be accessed by PyMOL's get_model
        #pseudoatom is not available in version 0.99rc!
        knotname_atoms=knotname+"_ATOMS"
        cmd.delete(knotname_atoms)
        atomcount=0
        for i in range(len(backbone)-1):
            #Interpolate the coordinates of atoms to place between "joints"
            for j in range(numatoms): #Number of atoms per segment or step
                x=backbone[i][0]+self.step[sequence[i]][0]*j*1.0/numatoms
                y=backbone[i][1]+self.step[sequence[i]][1]*j*1.0/numatoms
                z=backbone[i][2]
                #print x,y,z
                atomcount+=1
                cmd.pseudoatom(knotname_atoms,pos=(x,y,z),
                               resi="%d" % (atomcount),name="CA",hetatm=0,vdw=1,segi="",b=1.0,q=1.0,elem="C")

    def drawBackbone(self,backbone,knotname):
        i=0
        obj=[]
        while i+1<len(backbone):
            x1=backbone[i][0]
            y1=backbone[i][1]
            z1=backbone[i][2]
            x2=backbone[i+1][0]
            y2=backbone[i+1][1]
            z2=backbone[i+1][2]
            obj.extend([SPHERE, x1, y1, z1, self.BACKBONE_RADIUS,
                        CYLINDER, x1, y1, z1, x2, y2, z2, self.BACKBONE_RADIUS,
                       1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
            i+=1
        cmd.delete(knotname)
        cmd.load_cgo(obj,knotname)

    #RH,RV and Peano by Dana
    def RH(self,s):
    	s = s.replace('R','l')
    	s = s.replace('L','r')
    	return s.upper()

    def RV(self,s):
    	s = s.replace('U','d')
    	s = s.replace('D','u')
    	return s.upper()

    def Peano(self,n):
    	target = 'DRURD'
    	for i in range(n-1):
            add = ''
            add+='D'+self.RH(target)+'D'
            add+=target+'R'
            add+=self.RV(target)+'U'
            add+=self.RH(self.RV(target))+'U'
            add+=self.RV(target)+'R'
            add+=target+'D'
            add+=self.RH(target)+'D'
            add+=target
            target+=add
    	return target


class PyKnotTools:
    
    def __init__(self,app):
        parent = app.root
        self.parent = parent

        # Create the main dialog.
        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('Exit PyKnot',),
                                 title = 'PyKnot PyMOL Knot Analysis Tool 1.2',
                                 command = self.close)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        page = self.notebook.add('Knot Analyzer')
        self.KA=KAControlGroup(page,
                               groupname='Analyze a structure',
                               #defaultstructurename='2efv',
                               defaultstructurename='',
                               #defaultchain='A',
                               defaultchain='')

        page = self.notebook.add('Link Analyzer')
        self.LA=LinkControlGroup(page,
                               groupname='Analyze linked structures',
                               #defaultstructurename='2efv',
                               #defaultstructurename='',
                               #defaultchain='A',
                               #defaultchain=''
                               )

        page = self.notebook.add('Knot Creator')
        self.KC=KnotCreatorControlGroup(page,
                               groupname='Create a knot',
                               #defaultstructurename='2efv',
                               #defaultstructurename='',
                               #defaultchain='A',
                               #defaultchain=''
                               )

        page = self.notebook.add('Credits')
        group = Pmw.Group(page,tag_text='Credits')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        creditstext='''
        PyKnot PyMOL Knot Analysis Tool
        Rhonald Lua
        Baylor College of Medicine
        Houston, Texas

        Thank you to Alexander Grosberg, for introducing me to the world of knots.
        Thank you to W. L. DeLano for PyMOL, and to the developers of
        apbs_tools.py and remote_pdb_load.py for providing these examples
        on how to write a PyMOL plugin, and to others who made their
        python scripts available on the Web.
        
        '''
        credits = Tkinter.Label(group.interior(),
                                text = creditstext,
                                background = 'black',
                                foreground = 'white',
                                #pady = 20,
                                )
        credits.pack(expand = 1, fill = 'both', padx = 4, pady = 4)

        #Sometimes, we get a name that ends with an '_', e.g. '1tem_'.
        #Underscores are not allowed in tab names for a notebook, so we have to get rid of it.
        #if cmd._et_tools_selectpage.replace('_',' ') in self.notebook.pagenames():
        #    self.notebook.selectpage(cmd._et_tools_selectpage.replace('_',' '))

        self.notebook.pack(fill='both',expand=1,padx=10,pady=10)
        self.notebook.setnaturalsize()
        self.dialog.show()


    def close(self, result):
        if 1:
        #if result=='Exit PyKnot tools':
            if __name__ == '__main__':
                #
                # dies with traceback, but who cares
                #
                self.parent.destroy()
            else:
                #self.dialog.deactivate(result)
                self.dialog.withdraw() #This just hides the window, does not destroy and free resources?

# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = Tk()
    Pmw.initialise(app.root)
    app.root.title('Exit button for __main__')
    
    widget = PyKnotTools(app)
    exitButton = Button(app.root, text = 'Exit', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
